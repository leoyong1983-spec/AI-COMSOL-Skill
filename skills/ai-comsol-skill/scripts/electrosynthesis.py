"""Project-scoped COMSOL reaction CLI. Use py -3.12 -B."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import sys
from workspace import comsol_install
from datetime import datetime, timezone
from pathlib import Path
from reaction_core import Reactor, InputError, SIDES, pools, owned_path
from comsol_native import run_native

def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")

def load(path, project_root=None):
    from workspace import bind_workspace
    c=json.loads(Path(path).read_text(encoding="utf-8-sig"))
    return Reactor(bind_workspace(c, project_root))

def csv_write(path,header,rows):
    with Path(path).open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.writer(f);writer.writerow(header);writer.writerows(rows)

def publish(reactor,rows,output,native_evidence=None,agreement=None):
    c=reactor.c
    checks=[reactor.checks(t,y) for t,y in rows]
    envelope={
      "species_balance_max_abs_mol":max(x["species_balance_max_abs_mol"] for x in checks),
      "element_balance_max_abs_mol_atoms":max(abs(v) for x in checks for v in x["element_balance_mol_atoms"].values()),
      "bulk_charge_max_abs_mol_charge_m3":max(abs(v) for x in checks for v in x["bulk_charge_mol_charge_m3"].values()),
      "electrode_current_max_abs_residual_A":max(abs(v) for x in checks for v in x["electrode_current_residual_A"].values()),
      "cumulative_electrode_charge_max_relative_residual":max(v for x in checks for v in x["cumulative_electrode_charge_relative_residual"].values()),
      "minimum_liquid_concentration_mol_m3":min(x["minimum_liquid_concentration_mol_m3"] for x in checks)}
    tolerances={"species_balance_max_abs_mol":2e-6,"element_balance_max_abs_mol_atoms":1e-5,"bulk_charge_max_abs_mol_charge_m3":1e-3,"electrode_current_max_abs_residual_A":1e-8,"cumulative_electrode_charge_max_relative_residual":1e-7}
    if any(envelope[k]>limit for k,limit in tolerances.items()) or envelope["minimum_liquid_concentration_mol_m3"] < -1e-7:
        write_json(output/"failed-balance-checks.json",{"checks":envelope,"limits":tolerances})
        raise InputError("mass/charge/current/nonnegativity check failed")
    csv_write(output/"state_history.csv",["time_s"]+["|".join(k) for k in reactor.keys],[[t]+list(y) for t,y in rows])
    csv_write(output/"reaction_history.csv",["time_s"]+[r["id"]+"|rate_mol_s" for r in c["reactions"]]+[r["id"]+"|FE" for r in c["reactions"]],[[t]+[x["reaction_rate_mol_s"][r["id"]] for r in c["reactions"]]+[x["branch_faradaic_efficiency"][r["id"]] for r in c["reactions"]] for (t,y),x in zip(rows,checks)])
    last_t,last=rows[-1]
    streams=[]
    for pool in pools(c):
        comp=c["compartments"][pool]
        streams.append({"stream_id":pool+"_liquid","state":"inventory-and-net-outlet","temperature_K":c["temperature_K"],"volume_m3":comp["volume_m3"],"net_outlet_m3_s":comp["flow_m3_s"],"true_species_mol_m3":{n:reactor.value(last,("c",pool,n)) for n in reactor.liquid},"cumulative_outlet_mol":{n:reactor.value(last,("out",pool,n)) for n in reactor.liquid},"cumulative_gas_mol":{n:reactor.value(last,("gas",pool,n)) for n in reactor.gas}})
    summary={"schema_version":1,"case_id":c["case_id"],"run_id":output.name,"project_id":c["project_id"],"created_utc":datetime.now(timezone.utc).isoformat(),"data_status":c["data_status"],"purpose":c["purpose"],"solver":"COMSOL native Global ODE/DAE" if native_evidence else "independent Python RK4 reference","model_scope":"0D isothermal, fully mixed, forward Tafel or prescribed FE; no resolved spatial fields","project_experimental_calibration":"not-performed-no-project-experimental-data","engineering_acceptance":"not-assessed","elapsed_model_time_s":last_t,"checks":envelope,"check_limits_for_software_validation":tolerances,"last_point":checks[-1],"streams":streams,"native_evidence":native_evidence,"independent_reference_agreement":agreement,"provenance":c.get("provenance"),"assumptions":c["assumptions"],"not_available":["local current/flow/temperature fields","validated lifetime/deactivation","plant net energy/TEA/LCA","product purity after separation"]}
    write_json(output/"config_used.json",c)
    config_hash=hashlib.sha256((output/"config_used.json").read_bytes()).hexdigest()
    summary["config_sha256"]=config_hash
    summary["performance"]=reactor.performance(last_t,last)
    if native_evidence:summary["native_mph_sha256"]=hashlib.sha256(Path(native_evidence["native_model"]).read_bytes()).hexdigest()
    write_json(output/"result_card.json",summary)
    md=["# 成对电合成反应计算结果","",f'- 工况：{c["case_id"]}。',f'- 参数性质：{c["data_status"]}；{c["purpose"]}。',f'- 执行：{summary["solver"]}；模型终点 {last_t:g} s。',"- 范围：0D等温完全混合；未校准项目实验数据，未作工程验收。","", "## 软件验证", "",f'- 元素闭合最大绝对差：{envelope["element_balance_max_abs_mol_atoms"]:.6g} mol原子。',f'- 物种库存/累计流出闭合差：{envelope["species_balance_max_abs_mol"]:.6g} mol。',f'- 体相电中性最大残差：{envelope["bulk_charge_max_abs_mol_charge_m3"]:.6g} mol电荷/m³。',f'- 两电极支路电流闭合最大差：{envelope["electrode_current_max_abs_residual_A"]:.6g} A。']
    if agreement:md += [f'- COMSOL与独立积分最大缩放差：{agreement["state_max_scaled_difference"]:.6g}。',f'- 保存重开最大缩放差：{native_evidence["save_reopen_max_scaled_difference"]:.6g}。']
    perf=summary["performance"]
    md += ["","## 累计反应量与转化","",f'- 乙二醇净反应消耗：{perf["eg_net_consumed_mol"]:.8g} mol。',f'- 乙醇酸根净生成：{perf["net_reaction_generation_mol"].get("GLY",0):.8g} mol。']
    if perf["eg_cumulative_conversion_initial_plus_feed_basis"] is not None:
        md.append(f'- 初始库存加累计进料基准的乙二醇累计转化率：{100*perf["eg_cumulative_conversion_initial_plus_feed_basis"]:.6g}%。该口径不是连续稳态单程转化率。')
    if perf["carbon_selectivity_net_generation_per_eg_consumed"] is not None:
        md.append('- 净生成量基准碳选择性：'+ '；'.join(f'{name} {100*value:.6g}%' for name,value in perf["carbon_selectivity_net_generation_per_eg_consumed"].items())+'。')
    md += ["","## 终点液体组成","", "| 槽室 | 物种 | 浓度 mol/m³ |","| --- | --- | --- |"]
    for stream in streams:
        for name,value in stream["true_species_mol_m3"].items():
            if value>1e-9 and name not in ("H2O","K","OH"):md.append(f'| {stream["stream_id"]} | {reactor.species[name]["label"]} | {value:.8g} |')
    md += ["","## 使用边界",""]+["- "+s for s in c["assumptions"]]+["","电功仅在已提供同一参考下的电极动力学/电势与欧姆损失时计算。反应热、分离产品纯度和50 kW/千吨级性能不能由本次0D结果自动得出。"]
    (output/"结果说明.md").write_text("\n".join(md),encoding="utf-8")
    return summary

def run(config,out,install,native,project_root=None):
    model=load(config,project_root);out=owned_path(model.c["project_root"],out,"07_handoff")
    reference=model.integrate()
    if native:
        table,labels,evidence=run_native(model,out,install)
        if len(table)!=len(reference):raise InputError("native/reference time-grid mismatch")
        rows=[];error=0.0;jerror=0.0;state_ratio=0.0;current_ratio=0.0
        for data,(time,yref) in zip(table,reference):
            if abs(data[0]-time)>1e-6:raise InputError("native time grid differs")
            y=data[1:1+len(model.keys)];rows.append((data[0],y))
            error=max(error,max(abs(a-b)/max(1,abs(b)) for a,b in zip(y,yref)))
            state_ratio=max(state_ratio,max(abs(a-b)/(1e-9+5e-5*max(abs(a),abs(b))) for a,b in zip(y,yref)))
            actual_j=model.current_densities(y)[0]
            for r in model.c["reactions"]:
                col=labels.index("j_A_m2|"+r["id"])
                jerror=max(jerror,abs(data[col]-actual_j[r["id"]])/max(1,abs(actual_j[r["id"]])))
                current_ratio=max(current_ratio,abs(data[col]-actual_j[r["id"]])/(1e-8+5e-5*max(abs(data[col]),abs(actual_j[r["id"]]))))
        agreement={"state_max_scaled_difference":error,"native_current_max_scaled_difference":jerror,"relative_tolerance":5e-5,"state_absolute_tolerance_SI":1e-9,"current_density_absolute_tolerance_A_m2":1e-8,"state_max_tolerance_ratio":state_ratio,"current_max_tolerance_ratio":current_ratio}
        if state_ratio>1 or current_ratio>1:raise InputError("native solver disagrees with independent current roots/RK4")
        return publish(model,rows,out,evidence,agreement)
    if out.exists():raise InputError("output directory already exists")
    out.mkdir(parents=True)
    write_json(out/"config_used.json",model.c)
    return publish(model,reference,out)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=["validate","reference","native"])
    parser.add_argument("--config",required=True)
    parser.add_argument("--out")
    parser.add_argument("--project-root", default=None, help="Selected workspace; defaults to current directory")
    parser.add_argument("--comsol-root",default=None)
    args=parser.parse_args()
    try:
        if args.command=="validate":
            m=load(args.config,args.project_root)
            print(json.dumps({"status":"input-structure-and-stoichiometry-passed","states":len(m.keys),"data_status":m.c["data_status"],"experimental_calibration":False},ensure_ascii=False));return 0
        if not args.out:raise InputError("--out is required")
        report=run(args.config,args.out,comsol_install(args.comsol_root) if args.command=="native" else None,args.command=="native",args.project_root)
        print(json.dumps({"status":"software-checks-passed","case_id":report["case_id"],"output":str(Path(args.out).resolve()),"checks":report["checks"],"agreement":report["independent_reference_agreement"],"engineering_acceptance":"not-assessed"},ensure_ascii=False))
        return 0
    except (InputError,ValueError,KeyError) as ex:
        print(json.dumps({"status":"rejected","reason":str(ex)},ensure_ascii=False),file=sys.stderr);return 2

if __name__=="__main__":raise SystemExit(main())
