"""Run the installed DWSIM neutral reference-property engine with confined artifacts."""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path
from reaction_core import owned_path, InputError
from workspace import dwsim_install

def process(raw_file,output):
    raw=json.loads(Path(raw_file).read_text(encoding="utf-8-sig"))
    if raw.get("status")!="calculated-with-property-level-status":raise InputError("DWSIM run failed: "+str(raw.get("error")))
    records=[]
    fields={"density":("density_kg_m3","kg/m3"),"viscosity":("viscosity_Pa_s","Pa*s"),"thermalConductivity":("thermal_conductivity_W_m_K","W/(m*K)"),"heatCapacityCp":("heat_capacity_J_kg_K","J/(kg*K)")}
    for source in raw["rows"]:
        row={k:source[k] for k in ("temperature_K","pressure_Pa_abs","EG_mole_fraction")}
        status={}
        for key,(name,unit) in fields.items():
            a=source.get(key)
            value=a[0] if isinstance(a,list) and len(a)==1 else None
            if not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
                row[name]=None;status[key]="rejected-zero-negative-or-missing-native-value"
            else:row[name]=value;status[key]="native-calculated-positive-units-CAPE-OPEN-SI"
        row["property_status"]=status;records.append(row)
    report={"schema_version":1,"engine":"DWSIM CalculatorInterface native .NET","version":raw.get("dwsim_thermodynamics_version"),"source_sha256":hashlib.sha256(Path(raw_file).read_bytes()).hexdigest(),"property_package":"Raoult's Law","composition_basis":"mole fraction EG in neutral water/EG","pressure_basis":"absolute Pa","validity":"reference properties only; no electrolyte salts or activities; not KOH electrolyte validation","rows":records,"excluded":"native heatCapacityCp returned zero in tested interface; kept null, not used in reactor heat balance","engineering_acceptance":"not-assessed"}
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    (output/"dwsim-property-reference.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    keys=["temperature_K","pressure_Pa_abs","EG_mole_fraction"]+[name for name,unit in fields.values()]
    with (output/"dwsim-property-reference.csv").open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore");writer.writeheader();writer.writerows(records)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--out",required=True);p.add_argument("--dwsim-root",default=None);p.add_argument("--project-root",default=None)
    a=p.parse_args();root=Path(a.project_root or Path.cwd()).resolve();out=owned_path(root,a.out,"07_handoff")
    if out.exists():raise InputError("choose a new output directory")
    install=dwsim_install(a.dwsim_root);work=owned_path(root,root/"AGENTCODING/04_work/dwsim_properties"/out.name);qa=owned_path(root,root/"AGENTCODING/05_qa/dwsim_properties"/out.name)
    if work.exists() or qa.exists():raise InputError("run ID already exists")
    work.mkdir(parents=True);qa.mkdir(parents=True);temp=work/"temp";temp.mkdir()
    env=os.environ.copy()
    for key in ("TEMP","TMP","TMPDIR"):env[key]=str(temp)
    csc=Path(r"C:\Windows\Microsoft.NET\Framework64\v4.0.30319\csc.exe");exe=work/"DwsimProperties.exe"
    args=[str(csc),"/nologo","/target:exe","/platform:x64","/r:Microsoft.CSharp.dll","/r:System.Web.Extensions.dll","/out:"+str(exe),str(Path(__file__).with_name("DwsimProperties.cs"))]
    with (qa/"compile.log").open("wb") as log:
        subprocess.run(args,env=env,cwd=work,stdout=log,stderr=subprocess.STDOUT,timeout=90,check=True)
    with (qa/"run.log").open("wb") as log:
        subprocess.run([str(exe),str(install),str(qa)],env=env,cwd=work,stdout=log,stderr=subprocess.STDOUT,timeout=240,check=True)
    report=process(qa/"dwsim-properties-raw.json",out)
    print(json.dumps({"status":"property-level-checked","rows":len(report["rows"]),"version":report["version"],"out":str(out),"Cp":"rejected-native-zero"},ensure_ascii=False))

if __name__=="__main__":main()
