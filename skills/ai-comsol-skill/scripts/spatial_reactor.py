"""Native 2D/3D laminar flow, dilute organic transport and secondary current scenario."""
from pathlib import Path
import argparse,json,math,os,subprocess,sys,csv,hashlib
from datetime import datetime,timezone
from reaction_core import InputError,owned_path,F
from comsol_native import jstr,jarray,decode
from numerical_checks import verify_readback
from spatial_contract import validate_spatial,applicability,anode_current_density,conservation_checks

SPECIES=("EG","GLY","FMA","OX")
from workspace import bind_workspace, comsol_install
def load_config(path, project_root=None):
 c=json.loads(Path(path).read_text(encoding="utf-8-sig"))
 c=bind_workspace(c, project_root)
 return validate_spatial(c)

def generate(c,out):
 d=c["dimension"];bd=d-1;depth="W" if d==2 else "1";fd=c.get("flow_model")=="fully-developed-poiseuille"
 p={"L":f'{c["geometry"]["length_m"]}[m]',"H":f'{c["geometry"]["gap_m"]}[m]',"W":f'{c["geometry"]["width_m"]}[m]',"Q":f'{c["operation"]["flow_m3_s"]}[m^3/s]',"cin":f'{c["operation"]["feed_EG_mol_m3"]}[mol/m^3]',"Tiso":f'{c["operation"]["temperature_K"]}[K]',"Ea":f'{c["operation"]["anode_metal_potential_V_RHE"]}[V]',"Ec":f'{c["operation"]["cathode_metal_potential_V_RHE"]}[V]',"rho_liq":f'{c["properties"]["density_kg_m3"]}[kg/m^3]',"mu_liq":f'{c["properties"]["viscosity_Pa_s"]}[Pa*s]',"kappa":f'{c["properties"]["conductivity_S_m"]}[S/m]',"Fnum":f"{F}[C/mol]","ja_ref":f'{c["kinetics"]["anode_anchor_current_A_m2"]}[A/m^2]',"Ea_ref":f'{c["kinetics"]["anode_anchor_potential_V_RHE"]}[V]',"ba":f'{c["kinetics"]["anode_tafel_V_dec"]}[V]',"mEG":str(c["kinetics"]["EG_reaction_order"]),"jc_ref":f'{c["kinetics"]["cathode_anchor_current_A_m2"]}[A/m^2]',"ec_ref":f'{c["kinetics"]["cathode_anchor_overpotential_V"]}[V]',"bc":f'{c["kinetics"]["cathode_tafel_V_dec"]}[V]'}
 for n in SPECIES:p["D"+n]=f'{c["properties"][n+"_diffusivity_m2_s"]}[m^2/s]'
 p["cEGref"]=f'{c["kinetics"]["anode_EG_reference_mol_m3"]}[mol/m^3]'
 for n in SPECIES[1:]:p["fe"+n]=str(c["faradaic_efficiency"][n])
 ratio=c["geometry"]["gap_m"]/c["geometry"]["width_m"]
 p["f_rect"]=str(1-192*ratio/math.pi**5*sum(math.tanh(n*math.pi/(2*ratio))/n**5 for n in range(1,80,2)) if d==3 else 1)
 lines=["import com.comsol.model.*;","import com.comsol.model.util.*;","import java.util.Arrays;","public class SpatialReactor {","public static Model run() throws Exception {",'Model m=ModelUtil.create("EGSpatial");',f'm.label({jstr(c["case_id"]+"-"+str(d)+"D-literature-proxy")});',f'm.modelPath({jstr(out)});']
 def line(x):lines.append(x)
 for key,val in p.items():line(f'm.param().set({jstr(key)},{jstr(val)});')
 line('m.component().create("comp1",true);')
 line(f'm.component("comp1").geom().create("geom1",{d});')
 line('m.component("comp1").geom("geom1").lengthUnit("m");')
 shape="Rectangle" if d==2 else "Block"
 line(f'm.component("comp1").geom("geom1").create("g1","{shape}");')
 line(f'm.component("comp1").geom("geom1").feature("g1").set("size",{jarray(["L","H"] if d==2 else ["L","H","W"])});')
 line('m.component("comp1").geom("geom1").run();')
 line('System.out.println("GEOMETRY_BOUNDS "+Arrays.toString(m.component("comp1").geom("geom1").getBoundingBox()));')
 def box(tag,dim,mins,maxs):
  line(f'm.component("comp1").selection().create("{tag}","Box");')
  line(f'm.component("comp1").selection("{tag}").set("entitydim",{dim});')
  line(f'm.component("comp1").selection("{tag}").set("condition","inside");')
  for axis,lo,hi in zip(("x","y","z"),mins,maxs):
   def bound(expr):
    dimensions={"L":c["geometry"]["length_m"],"H":c["geometry"]["gap_m"],"W":c["geometry"]["width_m"]}
    if expr[0] in dimensions:return dimensions[expr[0]]+float(expr[1:] or 0)
    return float(expr)
   lo=lo.replace("1e-10","1e-7");hi=hi.replace("1e-10","1e-7")
   line(f'm.component("comp1").selection("{tag}").set("{axis}min",{bound(lo)});')
   line(f'm.component("comp1").selection("{tag}").set("{axis}max",{bound(hi)});')
 def union(tag,dim,inputs):
  line(f'm.component("comp1").selection().create("{tag}","Union");')
  line(f'm.component("comp1").selection("{tag}").set("entitydim",{dim});')
  line(f'm.component("comp1").selection("{tag}").set("input",{jarray(inputs)});')
 eps="1e-10"
 box("inlet",bd,["-1e-10","-1e-10","-1e-10"],["1e-10","H+1e-10","W+1e-10"])
 box("outlet",bd,["L-1e-10","-1e-10","-1e-10"],["L+1e-10","H+1e-10","W+1e-10"])
 box("anode",bd,["-1e-10","-1e-10","-1e-10"],["L+1e-10","1e-10","W+1e-10"])
 box("cathode",bd,["-1e-10","H-1e-10","-1e-10"],["L+1e-10","H+1e-10","W+1e-10"])
 if d==3:box("span0",2,["-1e-10","-1e-10","-1e-10"],["L+1e-10","H+1e-10","1e-10"])
 if d==3:box("span1",2,["-1e-10","-1e-10","W-1e-10"],["L+1e-10","H+1e-10","W+1e-10"])
 union("walls",bd,["anode","cathode"]+(["span0","span1"] if d==3 else []))
 if d==2:
  union("gap_edges",1,["inlet","outlet"])
  union("axial_edges",1,["anode","cathode"])
 for tag in ("inlet","outlet","anode","cathode"):
  line(f'System.out.println("SELECTION {tag} "+Arrays.toString(m.component("comp1").selection("{tag}").entities({bd})));')
  line(f'if(m.component("comp1").selection("{tag}").entities({bd}).length!=1)throw new RuntimeException("Wrong entity count: {tag}");')
 line('m.component("comp1").variable().create("v1");')
 variables={"j_an":"ja_ref*(cEG/cEGref)^mEG*10^((Ea-phil-Ea_ref)/ba)","j_ca":"-jc_ref*10^((phil-Ec-ec_ref)/bc)","rEG":"j_an/Fnum*(feGLY/4+feFMA/6+feOX/8)","rGLY":"j_an*feGLY/(4*Fnum)","rFMA":"2*j_an*feFMA/(6*Fnum)","rOX":"j_an*feOX/(8*Fnum)"}
 variables.update({"uF":f'wfd*Q/(intin(wfd)*{depth})',"vF":"0[m/s]","wF":"0[m/s]","pF":f'12*mu_liq*Q/(W*H^3*f_rect)*(H*W)/(intin(wfd)*{depth})*(L-x)' } if fd else {"uF":"u","vF":"v","wF":"w" if d==3 else "0[m/s]","pF":"p"})
 variables["speedF"]="sqrt(uF^2+vF^2+wF^2)"
 for key,val in variables.items():line(f'm.component("comp1").variable("v1").set({jstr(key)},{jstr(val)});')
 line('m.component("comp1").physics().create("spf","LaminarFlow","geom1");')
 line('m.component("comp1").physics("spf").prop("ShapeProperty").set("order_fluid",2);')
 for key,val in (("rho_mat","userdef"),("rho","rho_liq"),("mu_mat","userdef"),("mu","mu_liq")):line(f'm.component("comp1").physics("spf").feature("fp1").set("{key}","{val}");')
 initial_u="6*Q/(W*H)*(y/H)*(1-y/H)" if d==2 else "36*Q/(W*H)*(y/H)*(1-y/H)*(z/W)*(1-z/W)"
 line('m.component("comp1").physics("spf").feature("init1").set("u_init",'+jarray([initial_u,"0","0"])+');')
 line('m.component("comp1").physics("spf").feature("init1").set("p_init","12*mu_liq*Q/(W*H^3)*(L-x)");')
 line(f'm.component("comp1").physics("spf").create("inl1","Inlet",{bd});')
 line('m.component("comp1").physics("spf").feature("inl1").selection().named("inlet");')
 line('m.component("comp1").physics("spf").feature("inl1").set("BoundaryCondition","LaminarInflow");')
 line('m.component("comp1").physics("spf").feature("inl1").set("Uav","Q/(H*W)");')
 line(f'm.component("comp1").physics("spf").create("out1","Outlet",{bd});')
 line('m.component("comp1").physics("spf").feature("out1").selection().named("outlet");')
 line('m.component("comp1").physics("spf").feature("out1").set("p0","0[Pa]");')
 if fd:
  line('m.component("comp1").physics().create("fd","CoefficientFormPDE","geom1");')
  line('m.component("comp1").physics("fd").field("dimensionless").field("wfd");')
  line('m.component("comp1").physics("fd").field("dimensionless").component(new String[]{"wfd"});')
  tensor=["0","0","0","1"] if d==2 else ["0","0","0","0","1","0","0","0","1"]
  line('m.component("comp1").physics("fd").feature("cfeq1").set("c",'+jarray(tensor)+');')
  line('m.component("comp1").physics("fd").feature("cfeq1").set("a","0");')
  line('m.component("comp1").physics("fd").feature("cfeq1").set("f","12/(H^2*f_rect)");')
  line(f'm.component("comp1").physics("fd").create("wall","DirichletBoundary",{bd});')
  line('m.component("comp1").physics("fd").feature("wall").selection().named("walls");')
  line('m.component("comp1").physics("fd").feature("wall").set("r","0");')
 line('m.component("comp1").physics().create("tds","DilutedSpecies","geom1");')
 line('m.component("comp1").physics("tds").prop("ShapeProperty").set("order_concentration",2);')
 line('m.component("comp1").physics("tds").prop("AdvancedSettings").set("ConvectiveTerm","noncons");')
 if d==3:line('m.component("comp1").physics("tds").prop("MassConsistentStabilization").set("massCrosswindDiffusion",false);')
 line('m.component("comp1").physics("tds").field("concentration").component('+jarray(["c"+n for n in SPECIES])+');')
 line('m.component("comp1").physics("tds").feature("cdm1").set("u_src","userdef");')
 line('m.component("comp1").physics("tds").feature("cdm1").set("u",'+jarray(["uF","vF","wF"])+');')
 for n in SPECIES:line(f'm.component("comp1").physics("tds").feature("cdm1").set("D_c{n}",'+jarray(["D"+n,"0","0","0","D"+n,"0","0","0","D"+n])+');')
 line('m.component("comp1").physics("tds").feature("init1").set("initc",'+jarray(["cin","0","0","0"])+');')
 for tag,kind,sel in (("conc1","Concentration","inlet"),("out1","Outflow","outlet"),("flux1","FluxBoundary","anode")):
  line(f'm.component("comp1").physics("tds").create("{tag}","{kind}",{bd});')
  line(f'm.component("comp1").physics("tds").feature("{tag}").selection().named("{sel}");')
 line('m.component("comp1").physics("tds").feature("conc1").set("species",new String[]{"1","1","1","1"});')
 line('m.component("comp1").physics("tds").feature("conc1").set("c0",'+jarray(["cin","0","0","0"])+');')
 line('m.component("comp1").physics("tds").feature("flux1").set("species",new String[]{"1","1","1","1"});')
 line('m.component("comp1").physics("tds").feature("flux1").set("J0",'+jarray(["-rEG","rGLY","rFMA","rOX"])+');')
 line('m.component("comp1").physics().create("cd","SecondaryCurrentDistribution","geom1");')
 line('m.component("comp1").physics("cd").prop("ShapeProperty").set("order_electricpotentialionicphase",2);')
 line('m.component("comp1").physics("cd").feature("ice1").set("sigmal_mat","userdef");')
 line('m.component("comp1").physics("cd").feature("ice1").set("sigmal","kappa");')
 line('m.component("comp1").physics("cd").feature("init1").set("phil","0.04[V]*(1-y/H)-0.02[V]");')
 for tag,sel,potential,expr in (("an","anode","Ea","j_an"),("ca","cathode","Ec","j_ca")):
  line(f'm.component("comp1").physics("cd").create("{tag}","ElectrodeSurface",{bd});')
  line(f'm.component("comp1").physics("cd").feature("{tag}").selection().named("{sel}");')
  line(f'm.component("comp1").physics("cd").feature("{tag}").set("phisext0","{potential}");')
  line(f'm.component("comp1").physics("cd").feature("{tag}").set("T_src","userdef");')
  line(f'm.component("comp1").physics("cd").feature("{tag}").set("T","Tiso");')
  line(f'm.component("comp1").physics("cd").feature("{tag}").feature("er1").set("Eeq_mat","userdef");')
  line(f'm.component("comp1").physics("cd").feature("{tag}").feature("er1").set("Eeq","0[V]");')
  line(f'm.component("comp1").physics("cd").feature("{tag}").feature("er1").set("ilocmat_mat","userdef");')
  line(f'm.component("comp1").physics("cd").feature("{tag}").feature("er1").set("ilocmat","{expr}");')
 if d==2:
  for phy in ("tds","cd"):line(f'm.component("comp1").physics("{phy}").prop("dz").set("dz","W");')
 for tag,sel in (("intin","inlet"),("intout","outlet"),("inta","anode"),("intc","cathode")):
  line(f'm.component("comp1").cpl().create("{tag}","Integration");')
  line(f'm.component("comp1").cpl("{tag}").selection().named("{sel}");')
 for tag,sel,kind in (("avgin","inlet","Average"),("avgout","outlet","Average"),("mina","anode","Minimum"),("maxa","anode","Maximum")):
  line(f'm.component("comp1").cpl().create("{tag}","{kind}");')
  line(f'm.component("comp1").cpl("{tag}").selection().named("{sel}");')
 line('m.component("comp1").mesh().create("mesh1");')
 if d==3:
  box("xedge",1,["-1e-10","-1e-10","-1e-10"],["L+1e-10","1e-10","1e-10"])
  box("xedge2",1,["-1e-10","H-1e-10","-1e-10"],["L+1e-10","H+1e-10","1e-10"])
  box("yedge",1,["-1e-10","-1e-10","-1e-10"],["1e-10","H+1e-10","1e-10"])
  box("yedge2",1,["L-1e-10","-1e-10","-1e-10"],["L+1e-10","H+1e-10","1e-10"])
  union("axial_edges",1,["xedge","xedge2"])
  union("gap_edges",1,["yedge","yedge2"])
 line('m.component("comp1").mesh("mesh1").create("map1","Map");')
 if d==3:line('m.component("comp1").mesh("mesh1").feature("map1").selection().named("span0");')
 for tag,sel,n,ratio in (("disx","axial_edges",c["mesh"]["axial_elements"],None),("disy","gap_edges",c["mesh"]["gap_elements"],c["mesh"]["gap_ratio"])):
  line(f'm.component("comp1").mesh("mesh1").feature("map1").create("{tag}","Distribution");')
  line(f'm.component("comp1").mesh("mesh1").feature("map1").feature("{tag}").selection().named("{sel}");')
  if ratio:
   line(f'm.component("comp1").mesh("mesh1").feature("map1").feature("{tag}").set("type","predefined");')
   line(f'm.component("comp1").mesh("mesh1").feature("map1").feature("{tag}").set("elemcount",{int(n)});')
   line(f'm.component("comp1").mesh("mesh1").feature("map1").feature("{tag}").set("elemratio",{float(ratio)});')
   line(f'm.component("comp1").mesh("mesh1").feature("map1").feature("{tag}").set("symmetric",true);')
  else:line(f'm.component("comp1").mesh("mesh1").feature("map1").feature("{tag}").set("numelem",{int(n)});')
 if d==3:
  line('m.component("comp1").mesh("mesh1").create("swe1","Sweep");')
  line('m.component("comp1").mesh("mesh1").feature("swe1").create("disz","Distribution");')
  line('m.component("comp1").mesh("mesh1").feature("swe1").feature("disz").set("type","predefined");')
  line(f'm.component("comp1").mesh("mesh1").feature("swe1").feature("disz").set("elemcount",{int(c["mesh"]["span_elements"])});')
  line('m.component("comp1").mesh("mesh1").feature("swe1").feature("disz").set("elemratio",10.0);')
  line('m.component("comp1").mesh("mesh1").feature("swe1").feature("disz").set("symmetric",true);')
 line('m.component("comp1").mesh("mesh1").run();')
 line('System.out.println("MESH_BUILT");')
 line('m.study().create("flow");m.study("flow").create("stat","Stationary");')
 activation=["spf","off","fd","on","tds","off","cd","off"] if fd else ["spf","on","tds","off","cd","off"]
 line('m.study("flow").feature("stat").set("activate",'+jarray(activation)+');')
 line('m.study("flow").createAutoSequences("all");')
 line('m.sol("sol1").feature("s1").set("stol","1e-8");')
 line('m.sol("sol1").feature("s1").feature("fc1").set("maxiter",50);')
 line('m.sol("sol1").runAll();System.out.println("FLOW_SOLVED");')
 line('m.result().numerical().create("flowcheck","EvalGlobal");')
 line('m.result().numerical("flowcheck").set("data","dset1");')
 vn_flow="uF*nx+vF*ny"+("+wF*nz" if d==3 else "")
 line('m.result().numerical("flowcheck").set("expr",'+jarray([f'-intin({vn_flow})*{depth}',f'intout({vn_flow})*{depth}',"avgin(pF)-avgout(pF)"])+');')
 line('m.result().table().create("flowtbl","Table");m.result().numerical("flowcheck").set("table","flowtbl");m.result().numerical("flowcheck").setResult();')
 line(f'm.result().table("flowtbl").save({jstr(out/"flow_metrics.csv")});')
 line(f'm.save({jstr(out/"flow_solved.mph")});')
 line('m.study().create("std1");m.study("std1").create("stat","Stationary");')
 activation=["spf","off","fd","off","tds","on","cd","on"] if fd else ["spf","off","tds","on","cd","on"]
 line('m.study("std1").feature("stat").set("activate",'+jarray(activation)+');')
 line('m.study("std1").feature("stat").set("usesol","on");')
 line('m.study("std1").feature("stat").set("notsolmethod","sol");')
 line('m.study("std1").feature("stat").set("notstudy","flow");')
 line('m.study("std1").feature("stat").set("notsol","sol1");')
 line('m.study("std1").createAutoSequences("all");')
 line(f'm.sol("sol2").feature("s1").set("stol",{jstr(str(c["solver"]["relative_tolerance"]))});')
 line('if(!Arrays.asList(m.sol("sol2").feature("s1").feature().tags()).contains("fc1")){m.sol("sol2").feature("s1").create("fc1","FullyCoupled");if(Arrays.asList(m.sol("sol2").feature("s1").feature().tags()).contains("se1"))m.sol("sol2").feature("s1").feature("se1").active(false);}')
 line('m.sol("sol2").feature("s1").create("dSpatial","Direct");m.sol("sol2").feature("s1").feature("dSpatial").set("linsolver","pardiso");m.sol("sol2").feature("s1").feature("fc1").set("linsolver","dSpatial");')
 line('m.sol("sol2").feature("s1").feature("fc1").set("maxiter",80);')
 for tag,scale in (("comp1_cEG","1000"),("comp1_cGLY","100"),("comp1_cFMA","30"),("comp1_cOX","10"),("comp1_phil","0.05")):
  line(f'm.sol("sol2").feature("v1").feature("{tag}").set("scalemethod","manual");')
  line(f'm.sol("sol2").feature("v1").feature("{tag}").set("scaleval","{scale}");')
 line('m.sol("sol2").runAll();System.out.println("SPATIAL_SOLVED");')
 vn="uF*nx+vF*ny"+("+wF*nz" if d==3 else "")
 grad=lambda n:"c"+n+"x*nx+c"+n+"y*ny"+("+c"+n+"z*nz" if d==3 else "")
 metrics={"Qin_m3_s":f"-intin({vn})*{depth}","Qout_m3_s":f"intout({vn})*{depth}","I_anode_A":f"inta(j_an)*{depth}","I_cathode_A":f"intc(j_ca)*{depth}","EG_surface_min_mol_m3":"mina(cEG)","EG_surface_max_mol_m3":"maxa(cEG)","j_anode_min_A_m2":"mina(j_an)","j_anode_max_A_m2":"maxa(j_an)","pressure_drop_Pa":"avgin(pF)-avgout(pF)","EG_consumption_mol_s":f"inta(rEG)*{depth}","GLY_generation_mol_s":f"inta(rGLY)*{depth}","FMA_generation_mol_s":f"inta(rFMA)*{depth}","OX_generation_mol_s":f"inta(rOX)*{depth}"}
 for n in SPECIES:
  flux=f"({vn})*c{n}-D{n}*({grad(n)})"
  metrics[n+"_net_boundary_out_mol_s"]=f"(intin({flux})+intout({flux}))*{depth}"
 def table(tag,path):
  line(f'm.result().numerical().create("{tag}","EvalGlobal");')
  line(f'm.result().numerical("{tag}").set("data","dset2");')
  line(f'm.result().numerical("{tag}").set("expr",{jarray(metrics.values())});')
  line(f'm.result().numerical("{tag}").set("descr",{jarray(metrics.keys())});')
  line(f'm.result().table().create("{tag}tbl","Table");m.result().numerical("{tag}").set("table","{tag}tbl");m.result().numerical("{tag}").setResult();m.result().table("{tag}tbl").save({jstr(path)});')
 plot_start=len(lines)
 for tag,expr,label in (("pgEG","cEG","EG concentration mol/m3"),("pgU","speedF","Velocity m/s"),("pgV","phil","Electrolyte potential V")):
  line(f'm.result().create("{tag}","PlotGroup{d}D");m.result("{tag}").label("{label}");')
  line(f'm.result("{tag}").set("data","dset2");')
  feature="Surface" if d==2 else "Slice"
  line(f'm.result("{tag}").create("view1","{feature}");m.result("{tag}").feature("view1").set("expr","{expr}");')
 plot_lines=lines[plot_start:];del lines[plot_start:]
 table("check",out/"native_metrics.csv")
 line('m.result().export().create("fields","Data");')
 line('m.result().export("fields").set("data","dset2");')
 line('m.result().export("fields").set("expr",'+jarray(["cEG","cGLY","cFMA","cOX","speedF","pF","phil","j_an"])+');')
 line(f'm.result().export("fields").set("filename",{jstr(out/"native_fields.txt")});')
 line('m.result().export("fields").run();')
 line(f'm.save({jstr(out/"reactor_spatial.mph")});')
 line('try {')
 lines.extend(plot_lines)
 line(f'm.save({jstr(out/"reactor_spatial.mph")});')
 line('}catch(Exception plotError){System.out.println("NATIVE_PLOT_WARNING "+plotError.getMessage());}')
 line('ModelUtil.remove("EGSpatial");')
 line(f'm=ModelUtil.load("Reopened",{jstr(out/"reactor_spatial.mph")});')
 table("recheck",out/"native_reopened_metrics.csv")
 line(f'm.result().export("fields").set("filename",{jstr(out/"native_reopened_fields.txt")});m.result().export("fields").run();')
 line('System.out.println("SPATIAL_SAVED_REOPENED");return m;}public static void main(String[] a)throws Exception{run();}}')
 return "\n".join(lines),list(metrics)

def numeric_table(path,n):
 rows=[]
 for row in csv.reader(decode(path).splitlines()):
  if not row or row[0].lstrip().startswith("%"):continue
  try:vals=[float(x) for x in row]
  except ValueError:raise InputError("invalid native metric row")
  if len(vals)!=n:raise InputError("native metric width mismatch")
  rows.append(vals)
 if len(rows)!=1 or not all(math.isfinite(v) for v in rows[0]):raise InputError("native metric width/count/nonfinite")
 return rows[0]

def field_rows(path,dimension):
 values=[]
 for row in decode(path).splitlines():
  if not row.strip() or row.lstrip().startswith("%"):continue
  try:nums=[float(x) for x in row.split()]
  except ValueError:raise InputError("invalid native field row")
  if len(nums)!=dimension+8 or not all(math.isfinite(x) for x in nums):raise InputError("field width or finite-value check failed")
  values.append(nums)
 if not values:raise InputError("empty native field export")
 return values

def field_checks(path,dimension,feed,config=None):
 values=field_rows(path,dimension)
 result={"exported_points":len(values),"species_min_mol_m3":{n:min(row[dimension+i] for row in values) for i,n in enumerate(SPECIES)},"species_max_mol_m3":{n:max(row[dimension+i] for row in values) for i,n in enumerate(SPECIES)}}
 if result["species_min_mol_m3"]["EG"]<=0 or result["species_max_mol_m3"]["EG"]>feed*1.005:raise InputError("EG maximum principle or positivity violated")
 if any(v < -feed*1e-5 for v in result["species_min_mol_m3"].values()):raise InputError("negative species field beyond numerical tolerance")
 if config:
  error=0.0
  for row in values:
   expected=anode_current_density(config,row[dimension],row[dimension+6])
   error=max(error,abs(row[dimension+7]-expected)/max(abs(expected),1e-12))
  if error>1e-7:raise InputError("exported current disagrees with fixed-reference kinetic law")
  result["anode_kinetic_field_relative_residual"]=error
  fraction=max((row[dimension+1]+row[dimension+2]+2*row[dimension+3])/config["operation"]["support_electrolyte_mol_m3"] for row in values)
  result["maximum_product_charge_equivalent_fraction_of_support"]=fraction
  result["support_electrolyte_approximation_warning"]=fraction>.1
  result["support_fraction_screen_note"]="0.1 is an engineering review indicator, not a validated accuracy threshold; no spatial ionic balance is solved"
 return result

def run(c,out,install,timeout=1200):
 validate_spatial(c)
 root=Path(c["project_root"]);out=owned_path(root,out,"07_handoff")
 if out.exists():raise InputError("choose a new result directory")
 work=owned_path(root,root/"AGENTCODING/04_work/comsol_spatial_runtime"/out.name)
 qa=owned_path(root,root/"AGENTCODING/05_qa/comsol_spatial_runtime"/out.name)
 if work.exists() or qa.exists():raise InputError("run ID already exists")
 for path in (out,work,qa,work/"temp"):path.mkdir(parents=True)
 env=os.environ.copy()
 for key in ("TEMP","TMP","TMPDIR"):env[key]=str(work/"temp")
 java,labels=generate(c,out);source=work/"SpatialReactor.java";source.write_text(java,encoding="utf-8")
 (out/"config_used.json").write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding="utf-8")
 (qa/"metric-labels.json").write_text(json.dumps(labels,indent=2),encoding="utf-8")
 install=Path(install)
 args=[str(install/"java/win64/jre/bin/javac.exe"),"-J-Duser.language=en","-J-Dfile.encoding=UTF-8","-proc:none","-encoding","UTF-8","-cp",str(install/"plugins/*"),"-d",str(work),str(source)]
 with (qa/"compile.log").open("wb") as f:p=subprocess.run(args,cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=90)
 if p.returncode:raise InputError("compile failed; "+str(qa))
 args=[str(install/"bin/win64/comsolbatch.exe"),"-np","2","-prefsdir",str(work/"prefs"),"-tmpdir",str(work/"temp"),"-configuration",str(work/"configuration"),"-data",str(work/"workspace"),"-recoverydir",str(work/"recovery"),"-inputfile",str(work/"SpatialReactor.class"),"-nosave","-batchlog",str(qa/"batch.log")]
 (qa/"command.json").write_text(json.dumps(args,ensure_ascii=False,indent=2),encoding="utf-8")
 with (qa/"launcher.log").open("wb") as f:
  proc=subprocess.Popen(args,cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT)
  try:proc.wait(timeout=timeout)
  except subprocess.TimeoutExpired:
   if os.name=="nt":subprocess.run(["taskkill","/PID",str(proc.pid),"/T","/F"],stdout=f,stderr=subprocess.STDOUT,timeout=30)
   else:proc.kill()
   proc.wait(timeout=30)
   raise InputError("native timeout; owned process tree terminated; see "+str(qa))
  p=proc
 needed=["native_metrics.csv","native_reopened_metrics.csv","reactor_spatial.mph","native_fields.txt","native_reopened_fields.txt"]
 if p.returncode or any(not (out/n).is_file() or not (out/n).stat().st_size for n in needed):raise InputError("native solve/save/export/reopen incomplete; "+str(qa))
 values=numeric_table(out/"native_metrics.csv",len(labels));reread=numeric_table(out/"native_reopened_metrics.csv",len(labels))
 readback=verify_readback([values],[reread])
 field_readback=verify_readback(field_rows(out/"native_fields.txt",c["dimension"]),field_rows(out/"native_reopened_fields.txt",c["dimension"]))
 v=dict(zip(labels,values));ia=v["I_anode_A"]
 fields=field_checks(out/"native_fields.txt",c["dimension"],c["operation"]["feed_EG_mol_m3"],c)
 checks={**conservation_checks(c,v),**readback,**{"field_"+k:v for k,v in field_readback.items()}}
 card={"schema_version":1,"case_id":c["case_id"],"dimension":c["dimension"],"status":"native-spatial-solve-and-readback-passed","flow_model":c["flow_model"],"model_scope":"fully-developed Stokes-Poiseuille velocity PDE with integral flow normalization + passive dilute organic transport + secondary current distribution, membrane-free planar channel" if c["flow_model"]=="fully-developed-poiseuille" else "Navier-Stokes + passive dilute organic transport + secondary current distribution","engineering_acceptance":"not-assessed","project_experimental_calibration":"not-performed","metrics":v,"checks":checks,"field_checks":fields,"assumptions":c["assumptions"],"provenance":c["provenance"],"mesh":c["mesh"],"native_model_sha256":hashlib.sha256((out/"reactor_spatial.mph").read_bytes()).hexdigest(),"source_java_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),"evidence_directory":str(qa),"created_utc":datetime.now(timezone.utc).isoformat()}
 card.update({"artifact_type":"spatial-reactor-result","interface_version":"spatial-2","applicability":applicability(c),"data_status":c["data_status"],"config_sha256":hashlib.sha256((out/"config_used.json").read_bytes()).hexdigest(),"native_fields_sha256":hashlib.sha256((out/"native_fields.txt").read_bytes()).hexdigest(),"numerical_scheme":{"species":"quadratic SUPG","crosswind_species":c["dimension"]==2},"native_plot_status":"warning" if "NATIVE_PLOT_WARNING" in decode(qa/"launcher.log") else "saved"})
 (out/"result_card.json").write_text(json.dumps(card,ensure_ascii=False,indent=2),encoding="utf-8")
 return card

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("command",choices=("validate","native"));p.add_argument("--config",required=True);p.add_argument("--out");p.add_argument("--project-root",default=None);p.add_argument("--dimension",type=int,choices=(2,3));p.add_argument("--mesh-factor",type=float,default=1);p.add_argument("--flow-ml-min",type=float);p.add_argument("--comsol-root",default=None)
 a=p.parse_args()
 target=None;is_new=False
 try:
  c=load_config(a.config,a.project_root)
  if a.dimension:c["dimension"]=a.dimension
  if not math.isfinite(a.mesh_factor) or a.mesh_factor<=0:raise InputError("positive finite mesh-factor required")
  for k in ("axial_elements","gap_elements","span_elements"):c["mesh"][k]=max(4,round(c["mesh"][k]*a.mesh_factor))
  if a.flow_ml_min is not None:
   if not math.isfinite(a.flow_ml_min) or a.flow_ml_min<=0:raise InputError("positive flow required")
   c["operation"]["flow_m3_s"]=a.flow_ml_min*1e-6/60
  validate_spatial(c)
  if a.command=="validate":print(json.dumps({"status":"input-checked","dimension":c["dimension"],"applicability":applicability(c),"experimental_calibration":False}));return
  if not a.out:raise InputError("--out required")
  target=owned_path(c["project_root"],a.out,"07_handoff");is_new=not target.exists()
  card=run(c,target,comsol_install(a.comsol_root))
  print(json.dumps({"status":card["status"],"dimension":card["dimension"],"metrics":card["metrics"],"checks":card["checks"]},ensure_ascii=False))
 except (InputError,ValueError,KeyError,subprocess.TimeoutExpired) as ex:
  if is_new and target is not None and target.is_dir() and not (target/"result_card.json").exists():
   (target/"run-status.json").write_text(json.dumps({"status":"rejected","reason":str(ex),"engineering_use":False},ensure_ascii=False,indent=2),encoding="utf-8")
  print(json.dumps({"status":"rejected","reason":str(ex)},ensure_ascii=False),file=sys.stderr);raise SystemExit(2)
if __name__=="__main__":main()
