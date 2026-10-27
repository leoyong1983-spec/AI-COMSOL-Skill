"""Independent analytic velocity checks and controlled mesh comparison for spatial cases."""
from pathlib import Path
import argparse,hashlib,json,math,sys
from reaction_core import InputError,owned_path,finite_tree
from spatial_reactor import field_checks,field_rows
from spatial_contract import validate_spatial,conservation_checks
from numerical_checks import verify_readback

def read(path):return json.loads(Path(path).read_text(encoding="utf-8-sig"))
def write(path,obj,project_root=None):
 path=owned_path(project_root or Path.cwd(),path,"07_handoff")
 if path.exists():raise InputError("review output exists")
 path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")

def physics_key(c):
 return {k:c[k] for k in ("project_id","dimension","flow_model","geometry","operation","properties","kinetics","faradaic_efficiency","solver")}
def accepted_case(folder):
 folder=Path(folder);card=read(folder/"result_card.json");c=read(folder/"config_used.json")
 if card.get("status")!="native-spatial-solve-and-readback-passed" or card.get("interface_version")!="spatial-2":raise InputError("case must pass the current native solve/readback contract; rerun historical cases")
 validate_spatial(c);finite_tree(card,"result card")
 for name,key in (("config_used.json","config_sha256"),("reactor_spatial.mph","native_model_sha256"),("native_fields.txt","native_fields_sha256")):
  path=folder/name
  if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=card.get(key):raise InputError("case artifact hash mismatch: "+name)
 conservation_checks(c,card["metrics"])
 # Verify actual exported tables as well as the summary card, including every metric.
 from spatial_reactor import numeric_table
 values=numeric_table(folder/"native_metrics.csv",len(card["metrics"]))
 verify_readback([list(card["metrics"].values())],[values])
 verify_readback([values],[numeric_table(folder/"native_reopened_metrics.csv",len(values))])
 verify_readback(field_rows(folder/"native_fields.txt",c["dimension"]),field_rows(folder/"native_reopened_fields.txt",c["dimension"]))
 field_checks(folder/"native_fields.txt",c["dimension"],c["operation"]["feed_EG_mol_m3"],c)
 return card,c

def compare_data(a,b,ca,cb,limit=.02):
 if not math.isfinite(limit) or limit<=0:raise InputError("positive finite comparison limit required")
 finite_tree(a,"coarse card");finite_tree(b,"fine card")
 if physics_key(ca)!=physics_key(cb):raise InputError("mesh comparison requires identical physical inputs and dimension")
 for card in (a,b):
  if card.get("status")!="native-spatial-solve-and-readback-passed":raise InputError("rejected or incomplete case cannot pass mesh comparison")
 active=("axial_elements","gap_elements")+(("span_elements",) if ca["dimension"]==3 else ())
 if ca["mesh"]["gap_ratio"]!=cb["mesh"]["gap_ratio"]:raise InputError("mesh grading must remain fixed")
 if any(cb["mesh"][k]<ca["mesh"][k] for k in active) or not any(cb["mesh"][k]>ca["mesh"][k] for k in active):raise InputError("fine mesh must refine active dimensions without coarsening another dimension")
 names=("I_anode_A","EG_surface_min_mol_m3","j_anode_max_A_m2","GLY_generation_mol_s","pressure_drop_Pa")
 if any(not isinstance(card["metrics"].get(k),(int,float)) or card["metrics"][k]<=0 for card in (a,b) for k in names):raise InputError("missing/nonpositive comparison metric")
 differences={k:abs(a["metrics"][k]-b["metrics"][k])/max(abs(b["metrics"][k]),1e-15) for k in names}
 return {"status":"mesh-comparison-within-declared-limits" if max(differences.values())<=limit else "further-refinement-required","relative_limit":limit,"differences":differences,"coarse_mesh":ca["mesh"],"fine_mesh":cb["mesh"],"scope":"numerical convergence only; not experimental validation"}

def compare(coarse,fine,limit=.02):
 a,ca=accepted_case(coarse);b,cb=accepted_case(fine)
 return {**compare_data(a,b,ca,cb,limit),"coarse":str(Path(coarse).resolve()),"fine":str(Path(fine).resolve())}
def analytic_velocity(c,y,z=0):
 g=c["geometry"];q=c["operation"]["flow_m3_s"];h=g["gap_m"];w=g["width_m"]
 if c["dimension"]==2:return 6*q/(h*w)*(y/h)*(1-y/h)
 factor=1-192*h/(math.pi**5*w)*sum(math.tanh(n*math.pi*w/(2*h))/n**5 for n in range(1,160,2))
 total=0
 for n in range(1,160,2):
  ratio=(math.exp(-n*math.pi*z/h)+math.exp(-n*math.pi*(w-z)/h))/(1+math.exp(-n*math.pi*w/h))
  total+=(1-ratio)*math.sin(n*math.pi*y/h)/n**3
 return 48*q/(math.pi**3*w*h*factor)*total

def audit_case(folder):
 folder=Path(folder);card,c=accepted_case(folder)
 fields=field_checks(folder/"native_fields.txt",c["dimension"],c["operation"]["feed_EG_mol_m3"],c)
 if c.get("flow_model")!="fully-developed-poiseuille":raise InputError("analytic review only applies to fully developed family")
 errors=[]
 for a in field_rows(folder/"native_fields.txt",c["dimension"]):
  d=c["dimension"];exact=analytic_velocity(c,a[1],a[2] if d==3 else 0);velocity=a[d+4]
  errors.append(abs(velocity-exact))
 g=c["geometry"];q=c["operation"]["flow_m3_s"];mu=c["properties"]["viscosity_Pa_s"]
 corr=1 if c["dimension"]==2 else 1-192*g["gap_m"]/(math.pi**5*g["width_m"])*sum(math.tanh(n*math.pi*g["width_m"]/(2*g["gap_m"]))/n**5 for n in range(1,160,2))
 dp=12*mu*g["length_m"]*q/(g["width_m"]*g["gap_m"]**3*corr)
 velocity_error=max(errors)/(q/(g["gap_m"]*g["width_m"]));pressure_error=abs(card["metrics"]["pressure_drop_Pa"]-dp)/dp
 passed=velocity_error<=.01 and pressure_error<=.005
 return {"status":"independent-analytic-check-passed" if passed else "independent-analytic-check-failed","case":str(folder.resolve()),"fields":fields,"velocity_max_abs_error_m_s":max(errors),"velocity_max_error_relative_to_mean":velocity_error,"pressure_drop_analytic_Pa":dp,"pressure_drop_relative_error":pressure_error,"declared_numerical_limits":{"velocity_error_relative_to_mean":.01,"pressure_drop_relative_error":.005},"analytic_basis":"2D Poiseuille parabola or 3D rectangular-duct Fourier series, odd modes below 160; independent of native scalar PDE discretization","engineering_acceptance":"not-assessed"}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("command",choices=("analytic","mesh"));p.add_argument("--case");p.add_argument("--coarse");p.add_argument("--fine");p.add_argument("--out",required=True);p.add_argument("--project-root",default=None);a=p.parse_args()
 try:
  if a.command=="analytic":
   if not a.case:raise InputError("--case required")
   r=audit_case(a.case)
  else:
   if not a.coarse or not a.fine:raise InputError("--coarse and --fine required")
   r=compare(a.coarse,a.fine)
  write(a.out,r,a.project_root);print(json.dumps(r,ensure_ascii=False))
  if r["status"] not in ("mesh-comparison-within-declared-limits","independent-analytic-check-passed"):raise SystemExit(2)
 except (InputError,KeyError,ValueError) as e:print(str(e),file=sys.stderr);raise SystemExit(2)
if __name__=="__main__":main()
