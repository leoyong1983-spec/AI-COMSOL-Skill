"""Independent analytic velocity checks and controlled mesh comparison for spatial cases."""
from pathlib import Path
import argparse,json,math,sys
from reaction_core import InputError,owned_path
from spatial_reactor import field_checks,SPECIES
from comsol_native import decode

def read(path):return json.loads(Path(path).read_text(encoding="utf-8-sig"))
def write(path,obj,project_root=None):
 path=owned_path(project_root or Path.cwd(),path,"07_handoff")
 if path.exists():raise InputError("review output exists")
 path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")

def physics_key(c):
 return {k:c[k] for k in ("project_id","dimension","flow_model","geometry","operation","properties","kinetics","faradaic_efficiency")}
def compare(coarse,fine,limit=.02):
 a,b=read(Path(coarse)/"result_card.json"),read(Path(fine)/"result_card.json")
 ca,cb=read(Path(coarse)/"config_used.json"),read(Path(fine)/"config_used.json")
 if physics_key(ca)!=physics_key(cb):raise InputError("mesh comparison requires identical physical inputs and dimension")
 if ca["mesh"]==cb["mesh"]:raise InputError("mesh comparison needs different meshes")
 names=("I_anode_A","EG_surface_min_mol_m3","j_anode_max_A_m2","GLY_generation_mol_s","pressure_drop_Pa")
 differences={k:abs(a["metrics"][k]-b["metrics"][k])/max(abs(b["metrics"][k]),1e-15) for k in names}
 return {"status":"mesh-comparison-within-declared-limits" if max(differences.values())<=limit else "further-refinement-required","relative_limit":limit,"differences":differences,"coarse_mesh":ca["mesh"],"fine_mesh":cb["mesh"],"coarse":str(Path(coarse).resolve()),"fine":str(Path(fine).resolve()),"scope":"numerical convergence only; not experimental validation"}
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
 folder=Path(folder);c=read(folder/"config_used.json");card=read(folder/"result_card.json")
 fields=field_checks(folder/"native_fields.txt",c["dimension"],c["operation"]["feed_EG_mol_m3"])
 if c.get("flow_model")!="fully-developed-poiseuille":raise InputError("analytic review only applies to fully developed family")
 errors=[]
 for row in decode(folder/"native_fields.txt").splitlines():
  if not row.strip() or row.lstrip().startswith("%"):continue
  try:a=[float(v) for v in row.split()]
  except ValueError:continue
  d=c["dimension"];exact=analytic_velocity(c,a[1],a[2] if d==3 else 0);velocity=a[d+4]
  errors.append(abs(velocity-exact))
 g=c["geometry"];q=c["operation"]["flow_m3_s"];mu=c["properties"]["viscosity_Pa_s"]
 corr=1 if c["dimension"]==2 else 1-192*g["gap_m"]/(math.pi**5*g["width_m"])*sum(math.tanh(n*math.pi*g["width_m"]/(2*g["gap_m"]))/n**5 for n in range(1,160,2))
 dp=12*mu*g["length_m"]*q/(g["width_m"]*g["gap_m"]**3*corr)
 return {"status":"independent-analytic-check-computed","case":str(folder.resolve()),"fields":fields,"velocity_max_abs_error_m_s":max(errors),"velocity_max_error_relative_to_mean":max(errors)/(q/(g["gap_m"]*g["width_m"])),"pressure_drop_analytic_Pa":dp,"pressure_drop_relative_error":abs(card["metrics"]["pressure_drop_Pa"]-dp)/dp,"analytic_basis":"2D Poiseuille parabola or3D rectangular-duct Fourier series,160terms limit; independent of native scalar PDE discretization","engineering_acceptance":"not-assessed"}

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
 except (InputError,KeyError,ValueError) as e:print(str(e),file=sys.stderr);raise SystemExit(2)
if __name__=="__main__":main()
