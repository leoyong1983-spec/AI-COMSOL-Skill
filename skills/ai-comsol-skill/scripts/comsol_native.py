"""Generate a native COMSOL Global ODE/DAE model and verify saved-model readback."""
from __future__ import annotations
import csv
import json
import math
import os
import subprocess
from pathlib import Path
from reaction_core import Reactor, InputError, F, R, SIDES, pools, pool_for, owned_path

def jstr(s):return json.dumps(str(s),ensure_ascii=True)
def jarray(values):return "new String[]{"+",".join(jstr(s) for s in values)+"}"
def sumexpr(terms):return "("+("+".join(terms) if terms else "0")+")"

def generate(reactor:Reactor,output:Path):
    c=reactor.c
    params={"Cref":"1[mol/m^3]","Nref":"1[mol]","tref":"1[s]","Iref":"1[A]","Vref":"1[V]","Icmd":f'{c["current_A"]}[A]',"Tiso":f'{c["temperature_K"]}[K]',"Fnum":f'{F}[C/mol]',"Rnum":f'{R}[J/(mol*K)]'}
    if reactor.mem_names:params["Amem"]=f'{c["membrane"]["area_m2"]}[m^2]'
    var={}
    eq=[]
    names=[reactor.symbol[k] for k in reactor.keys]
    initial=[str(x) for x in reactor.initial]
    labels=["time_s"]+["|".join(k) for k in reactor.keys]
    expressions=["t"]+names
    initial_eta=reactor.current_densities(reactor.initial)[2]
    for side in SIDES:
        tag="c" if side=="cathode" else "a"
        params["A"+tag]=f'{c["electrode_area_m2"][side]}[m^2]'
    for side in pools(c):
        tag="s" if side=="cell" else "c" if side=="cathode" else "a"
        comp=c["compartments"][side]
        params["V"+tag]=f'{comp["volume_m3"]}[m^3]'
        params["Q"+tag]=f'{comp["flow_m3_s"]}[m^3/s]'
        for species in reactor.liquid:params["feed_"+tag+"_"+species]=f'{comp["feed_mol_m3"][species]}[mol/m^3]'
    for r in c["reactions"]:
        rid=r["id"];tag="c" if r["side"]=="cathode" else "a"
        if c["kinetic_mode"]=="tafel":
            k=r["kinetics"]
            params["jref_"+rid]=f'{k["j_ref_A_m2"]}[A/m^2]'
            params["csref_"+rid]=f'{k["substrate_ref_mol_m3"]}[mol/m^3]'
            params["alpha_"+rid]=str(k["alpha_eff"])
            params["offset_"+rid]=f'{k["eta_offset_V"]}[V]'
            sub=r.get("substrate")
            factor=f'({reactor.symbol[("c",pool_for(c,r["side"]),sub)]}*Cref/csref_{rid})^{k["order"]}' if sub else "1"
            var["j_"+rid]=f'jref_{rid}*{factor}*exp(alpha_{rid}*Fnum*(eta_{tag}*Vref-offset_{rid})/(Rnum*Tiso))'
        else:
            params["fe_"+rid]=str(r["faradaic_efficiency"])
            var["j_"+rid]=f'Icmd*fe_{rid}/A{tag}'
        var["rate_"+rid]=f'j_{rid}*A{tag}/({r["electron_number"]}*Fnum)'
        labels.extend(["j_A_m2|"+rid,"rate_mol_s|"+rid])
        expressions.extend(["j_"+rid,"rate_"+rid])
    if reactor.mem_names:
        carrier=c["membrane"]["carrier_species"]
        var["transfer_"+carrier]="Icmd/Fnum"
    for species,permeability in c["membrane"].get("neutral_permeability_m_s",{}).items():
        params["perm_"+species]=f'{permeability}[m/s]'
        var["transfer_"+species]=f'perm_{species}*Amem*Cref*({reactor.symbol[("c","cathode",species)]}-{reactor.symbol[("c","anode",species)]})'
    for key in reactor.keys:
        symbol=reactor.symbol[key]
        if key[0]=="c":
            _,side,species=key;tag="s" if side=="cell" else "c" if side=="cathode" else "a"
            flux=[f'Q{tag}*(feed_{tag}_{species}-{symbol}*Cref)']
            flux.extend(f'({r["stoichiometry"].get(species,0)})*rate_{r["id"]}' for r in c["reactions"] if pool_for(c,r["side"])==side and r["stoichiometry"].get(species,0))
            if species in reactor.mem_names:flux.append(("-" if side=="cathode" else "")+"transfer_"+species)
            rhs=sumexpr(flux)+f'/(V{tag}*Cref)'
        elif key[0]=="out":
            _,side,species=key;tag="s" if side=="cell" else "c" if side=="cathode" else "a"
            rhs=f'Q{tag}*Cref*{reactor.symbol[("c",side,species)]}/Nref'
        elif key[0]=="gas":
            _,side,species=key
            rhs=sumexpr([f'({r["stoichiometry"][species]})*rate_{r["id"]}' for r in c["reactions"] if pool_for(c,r["side"])==side and r["stoichiometry"].get(species,0)])+"/Nref"
        elif key[0]=="extent":rhs="rate_"+key[1]+"/Nref"
        else:rhs="transfer_"+key[1]+"/Nref"
        eq.append(f'tref*({symbol}t-({rhs}))')
    if c["kinetic_mode"]=="tafel":
        for side in SIDES:
            tag="c" if side=="cathode" else "a"
            names.append("eta_"+tag);initial.append(str(initial_eta[side]))
            eq.append(f'({sumexpr(["j_"+r["id"] for r in c["reactions"] if r["side"]==side])}*A{tag}-Icmd)/Iref')
            labels.append("eta_V|"+side);expressions.append("eta_"+tag)
    tlist=" ".join(format(t,".15g") for t in reactor.times())
    lines=['import com.comsol.model.*;','import com.comsol.model.util.*;','public class PairedReactor {',' public static Model run() {','  Model m=ModelUtil.create("PairedReaction");',f'  m.label({jstr(c["case_id"]+" - "+c["data_status"])});',f'  m.modelPath({jstr(output)});','  m.component().create("comp1",true);','  m.component("comp1").physics().create("ge","GlobalEquations");']
    for key,value in params.items():lines.append(f'  m.param().set({jstr(key)},{jstr(value)});')
    lines.append('  m.component("comp1").variable().create("v1");')
    for key,value in var.items():lines.append(f'  m.component("comp1").variable("v1").set({jstr(key)},{jstr(value)});')
    for prop,val in [('name',names),('equation',eq),('initialValueU',initial)]:
        lines.append(f'  m.component("comp1").physics("ge").feature("ge1").set({jstr(prop)},{jarray(val)});')
    lines+=['  m.study().create("std1");','  m.study("std1").create("time","Transient");',f'  m.study("std1").feature("time").set("tlist",{jstr(tlist)});',f'  m.study("std1").feature("time").set("rtol",{jstr(c["solver"]["relative_tolerance"])});','  m.study("std1").createAutoSequences("all");','  m.sol("sol1").feature("t1").set("consistent","off");','  m.sol("sol1").runAll();','  export(m,"ev","tbl",'+jstr(output/"native_raw.csv")+');','  try {m.save('+jstr(output/"reactor_solved.mph")+');}catch(Exception e){throw new RuntimeException(e);}','  ModelUtil.remove("PairedReaction");','  Model reopened=ModelUtil.load("Reopened",'+jstr(output/"reactor_solved.mph")+');','  export(reopened,"ev2","tbl2",'+jstr(output/"native_reopened_raw.csv")+');','  System.out.println("PAIRED_REACTOR_SOLVED_SAVED_REOPENED");','  return reopened;',' }',' static void export(Model m,String ev,String tbl,String path) {','  m.result().numerical().create(ev,"EvalGlobal");','  m.result().numerical(ev).set("expr",'+jarray(expressions)+');','  m.result().numerical(ev).set("descr",'+jarray(labels)+');','  m.result().table().create(tbl,"Table");','  m.result().numerical(ev).set("table",tbl);','  m.result().numerical(ev).setResult();','  try {m.result().table(tbl).save(path);}catch(Exception e){throw new RuntimeException(e);}',' }',' public static void main(String[] args){run();}','}']
    at='  m.sol("sol1").feature("t1").set("consistent","off");'
    controls='\n'.join([at,'  m.sol("sol1").feature("t1").set("rtol",'+jstr(c["solver"]["relative_tolerance"])+');','  m.sol("sol1").feature("t1").set("atolglobalvaluemethod","manual");','  m.sol("sol1").feature("t1").set("atolglobalmethod","unscaled");','  m.sol("sol1").feature("t1").set("atolglobal","1e-10");'])
    java="\n".join(lines).replace(at,controls).replace("public static Model run() {","public static Model run() throws Exception {").replace("public static void main(String[] args){run();}","public static void main(String[] args) throws Exception {run();}")
    return java,labels

def decode(path):
    raw=Path(path).read_bytes()
    for encoding in ("utf-8-sig","utf-16","gb18030"):
        try:return raw.decode(encoding)
        except UnicodeError:pass
    return raw.decode("utf-8",errors="replace")

def parse_table(path,labels):
    result=[]
    for row in csv.reader(decode(path).splitlines()):
        if not row or row[0].lstrip().startswith("%"):continue
        try:values=[float(v) for v in row]
        except ValueError:continue
        if len(values)==len(labels)+1:values=values[1:] # COMSOL adds its own time column.
        if len(values)!=len(labels):raise InputError("unexpected native table width")
        if not all(math.isfinite(v) for v in values):raise InputError("nonfinite native result")
        result.append(values)
    if not result:raise InputError("native table contains no numerical rows")
    return result

def run_native(reactor,output,install,timeout=300):
    root=Path(reactor.c["project_root"]).resolve();output=owned_path(root,output,"07_handoff")
    if output.exists():raise InputError("native output directory already exists; choose a new run directory")
    output.mkdir(parents=True)
    (output/"config_used.json").write_text(json.dumps(reactor.c,ensure_ascii=False,indent=2),encoding="utf-8")
    run_id=output.name
    work=owned_path(root,root/"AGENTCODING/04_work/comsol_runtime"/run_id)
    qa=owned_path(root,root/"AGENTCODING/05_qa/comsol_runtime"/run_id)
    if work.exists() or qa.exists():raise InputError("runtime run_id already exists")
    work.mkdir(parents=True);qa.mkdir(parents=True)
    temp=work/"temp";temp.mkdir()
    env=os.environ.copy()
    for key in ("TEMP","TMP","TMPDIR"):env[key]=str(temp)
    install=Path(install)
    javac=install/"java/win64/jre/bin/javac.exe";batch=install/"bin/win64/comsolbatch.exe"
    if not javac.is_file() or not batch.is_file():raise InputError("COMSOL native executables not found")
    java,labels=generate(reactor,output)
    source=work/"PairedReactor.java";source.write_text(java,encoding="utf-8")
    (qa/"column-map.json").write_text(json.dumps(labels,ensure_ascii=False,indent=2),encoding="utf-8")
    compile_args=[str(javac),"-J-Duser.language=en","-J-Dfile.encoding=UTF-8","-proc:none","-encoding","UTF-8","-cp",str(install/"plugins/*"),"-d",str(work),str(source)]
    with (qa/"compile.log").open("wb") as log:
        r=subprocess.run(compile_args,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=90)
    if r.returncode:raise InputError("Java compilation failed; see "+str(qa/"compile.log"))
    args=[str(batch),"-np","2","-prefsdir",str(work/"prefs"),"-tmpdir",str(temp),"-configuration",str(work/"configuration"),"-data",str(work/"workspace"),"-recoverydir",str(work/"recovery"),"-inputfile",str(work/"PairedReactor.class"),"-nosave","-batchlog",str(qa/"batch.log")]
    (qa/"command.json").write_text(json.dumps(args,ensure_ascii=False,indent=2),encoding="utf-8")
    with (qa/"launcher.log").open("wb") as log:
        r=subprocess.run(args,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    logs="\n".join(decode(path) for path in (qa/"launcher.log",qa/"batch.log") if path.exists())
    (qa/"combined-log-utf8.txt").write_text(logs,encoding="utf-8")
    needed=[output/"reactor_solved.mph",output/"native_raw.csv",output/"native_reopened_raw.csv"]
    if r.returncode or any(not p.is_file() or p.stat().st_size==0 for p in needed):
        raise InputError("COMSOL did not produce all solved/readback artifacts; exit status alone is insufficient: "+str(qa))
    native=parse_table(needed[1],labels);reopened=parse_table(needed[2],labels)
    if len(native)!=len(reopened):raise InputError("save/reopen changed row count")
    reopen_error=max(abs(a-b)/max(1,abs(a)) for ra,rb in zip(native,reopened) for a,b in zip(ra,rb))
    if reopen_error>1e-8:raise InputError("save/reopen numerical mismatch")
    return native,labels,{"comsol_exit_code":r.returncode,"save_reopen_max_scaled_difference":reopen_error,"evidence_directory":str(qa),"native_model":str(needed[0]),"generated_java":str(source),"install":str(install)}
