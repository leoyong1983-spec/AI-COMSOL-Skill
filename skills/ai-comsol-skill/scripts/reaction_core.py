"""Balanced shared-electrolyte or two-compartment reactor; SI inputs, stdlib only."""
from __future__ import annotations
import copy
import math
import re
from pathlib import Path

F = 96485.33212
R = 8.314462618
SIDES = ("cathode", "anode")

def pools(c):
    return ("cell",) if c.get("cell_topology")=="shared-electrolyte" else SIDES

def pool_for(c,side):
    return "cell" if c.get("cell_topology")=="shared-electrolyte" else side

class InputError(ValueError):
    pass

def finite_tree(value, path="input"):
    """Reject non-standard JSON numbers before comparisons can hide NaN."""
    if isinstance(value,float) and not math.isfinite(value):
        raise InputError(path+": nonfinite number")
    if isinstance(value,dict):
        for k,v in value.items():finite_tree(v,path+"."+str(k))
    elif isinstance(value,list):
        for i,v in enumerate(value):finite_tree(v,path+"["+str(i)+"]")

def number(value, path, errors, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        errors.append(path + ": requires a finite SI number")
        return
    if positive and value <= 0:
        errors.append(path + ": must be positive")
    if nonnegative and value < 0:
        errors.append(path + ": must be nonnegative")

def validate(c):
    if not isinstance(c,dict):raise InputError("configuration must be an object")
    finite_tree(c)
    for key in ("compartments","electrode_area_m2","membrane","time","solver"):
        if not isinstance(c.get(key),dict):raise InputError(key+": requires an object")
    for key in ("species","reactions"):
        if not isinstance(c.get(key),list) or any(not isinstance(v,dict) for v in c[key]):
            raise InputError(key+": requires a list of objects")
    if set(c["compartments"])!=set(pools(c)):
        raise InputError("compartments must exactly match the selected topology; unused inventories are not allowed")
    if any(not isinstance(p,dict) for p in c["compartments"].values()):raise InputError("compartments require objects")
    errors=[]
    if c.get("cell_topology") not in (None,"two-compartment","shared-electrolyte"):
        errors.append("unsupported cell_topology")
    if c.get("schema_version") != 1: errors.append("schema_version must be 1")
    if c.get("data_status") not in ("demonstration","literature-scenario","project-research"):
        errors.append("data_status must identify demonstration, literature-scenario or project-research")
    if c.get("data_status") != "demonstration":
        if not c.get("parameter_source"): errors.append("parameter_source is missing")
        for key in ("reaction_pair","reactor_geometry","electrolyte_speciation","kinetics_and_side_reactions","membrane","operating_conditions","electrode_potential_reference"):
            if not c.get("provenance",{}).get(key): errors.append("provenance."+key+" is missing")
    for key in ("current_A","temperature_K"):
        number(c.get(key),key,errors,positive=True)
    if c.get("kinetic_mode")=="tafel":
        for key in ("anode_reference_potential_V","cathode_reference_potential_V"):
            number(c.get(key),key,errors)
        number(c.get("ohmic_resistance_ohm"),"ohmic_resistance_ohm",errors,nonnegative=True)
    if c.get("kinetic_mode") not in ("tafel","prescribed-fe"): errors.append("unsupported kinetic_mode")
    species=c.get("species",[])
    if not species: errors.append("species is empty")
    sm={}
    for s in species:
        name=s.get("id","")
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*",name): errors.append("invalid species ID: "+name)
        if name in sm: errors.append("duplicate species ID: "+name)
        sm[name]=s
        if s.get("phase") not in ("liquid","gas","electronic"): errors.append("unsupported phase: "+name)
        atoms=s.get("atoms")
        if not isinstance(atoms,dict):errors.append("atoms missing: "+name)
        elif any(not isinstance(el,str) or not re.fullmatch(r"[A-Z][a-z]?",el) or isinstance(v,bool) or not isinstance(v,(int,float)) or v<=0 or v!=int(v) for el,v in atoms.items()):
            errors.append("invalid element/atom count: "+name)
        elif name!="electron" and not atoms:errors.append("non-electron species needs atoms: "+name)
        charge=s.get("charge")
        if isinstance(charge,bool) or not isinstance(charge,(int,float)) or charge!=int(charge):errors.append("integer species charge required: "+name)
        if s.get("phase")=="electronic" and name!="electron":errors.append("only the electron electronic species is supported")
        if s.get("phase")=="gas" and charge!=0:errors.append("charged gas is outside this model")
    if errors:raise InputError("\n".join(errors))
    if sm.get("electron",{}).get("charge") != -1: errors.append("electron species must have charge -1")
    if sm.get("electron",{}).get("phase")!="electronic" or sm.get("electron",{}).get("atoms")!={}:errors.append("electron must be electronic with no atoms")
    for side in SIDES:
        number(c.get("electrode_area_m2",{}).get(side),"electrode_area_m2."+side,errors,positive=True)
    for side in pools(c):
        comp=c.get("compartments",{}).get(side,{})
        number(comp.get("volume_m3"),side+".volume_m3",errors,positive=True)
        number(comp.get("flow_m3_s"),side+".flow_m3_s",errors,nonnegative=True)
        for field in ("feed_mol_m3","initial_mol_m3"):
            vals=comp.get(field)
            if not isinstance(vals,dict):
                errors.append(side+"."+field+" is missing"); continue
            for name,s in sm.items():
                if s.get("phase")=="liquid": number(vals.get(name),side+"."+field+"."+name,errors,nonnegative=True)
            if any(k not in sm or sm[k].get("phase")!="liquid" for k in vals): errors.append(side+"."+field+" has unknown/nonliquid species")
            if all(isinstance(vals.get(k),(int,float)) for k,s in sm.items() if s.get("phase")=="liquid"):
                charge=sum(vals[k]*s.get("charge",0) for k,s in sm.items() if s.get("phase")=="liquid")
                if abs(charge)>1e-7: errors.append(side+"."+field+" is not electroneutral; residual mol-charge/m3="+str(charge))
    reactions=c.get("reactions",[])
    if not reactions: errors.append("reactions is empty")
    seen=set()
    for r in reactions:
        rid=r.get("id","")
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*",rid) or rid in seen: errors.append("invalid/duplicate reaction ID: "+rid)
        seen.add(rid)
        side=r.get("side")
        if side not in SIDES: errors.append(rid+": invalid side")
        n=r.get("electron_number")
        number(n,rid+".electron_number",errors,positive=True)
        stoich=r.get("stoichiometry",{})
        if not isinstance(stoich,dict):raise InputError(rid+": stoichiometry must be an object")
        balance={}
        charge=0
        for name,coef in stoich.items():
            if name not in sm: errors.append(rid+": unknown species "+name); continue
            if isinstance(coef,bool) or not isinstance(coef,(int,float)): errors.append(rid+": invalid stoichiometry"); continue
            if sm[name].get("phase")=="gas" and coef<0:errors.append(rid+": gas consumption requires a gas-feed/inventory model")
            for el,num in sm[name].get("atoms",{}).items(): balance[el]=balance.get(el,0)+coef*num
            charge+=coef*sm[name].get("charge",0)
        if any(abs(v)>1e-10 for v in balance.values()): errors.append(rid+": element imbalance "+str(balance))
        if abs(charge)>1e-10: errors.append(rid+": charge imbalance "+str(charge))
        expected=n if side=="anode" else -n if isinstance(n,(int,float)) else None
        if stoich.get("electron") != expected: errors.append(rid+": electron count/sign disagrees with side")
        if c.get("kinetic_mode")=="tafel":
            k=r.get("kinetics",{})
            for key in ("j_ref_A_m2","alpha_eff","substrate_ref_mol_m3"):
                number(k.get(key),rid+"."+key,errors,positive=True)
            number(k.get("eta_offset_V"),rid+".eta_offset_V",errors)
            number(k.get("order"),rid+".order",errors,nonnegative=True)
            if not k.get("source"): errors.append(rid+": kinetics source is missing")
            sub=r.get("substrate")
            if sub and (sub not in sm or sm[sub].get("phase")!="liquid"): errors.append(rid+": invalid substrate")
        else:
            number(r.get("faradaic_efficiency"),rid+".faradaic_efficiency",errors,nonnegative=True)
    for side in SIDES:
        rs=[r for r in reactions if r.get("side")==side]
        if not rs: errors.append(side+": no reaction branches")
        if c.get("kinetic_mode")=="prescribed-fe" and all(isinstance(r.get("faradaic_efficiency"),(int,float)) for r in rs):
            if abs(sum(r["faradaic_efficiency"] for r in rs)-1)>1e-9:
                errors.append(side+": declared branch FE must sum to 1; unresolved current cannot be silently allocated")
    mem=c.get("membrane",{})
    carrier=mem.get("carrier_species")
    if c.get("cell_topology")=="shared-electrolyte" and mem.get("neutral_permeability_m_s"):
        errors.append("shared-electrolyte has no membrane permeability")
    if c.get("cell_topology")!="shared-electrolyte":
        if mem.get("carrier_direction")!="cathode-to-anode": errors.append("only explicit cathode-to-anode carrier supported in this model")
        if mem.get("carrier_mol_per_mol_e")!=1 or sm.get(carrier,{}).get("charge")!=-1 or sm.get(carrier,{}).get("phase")!="liquid":
            errors.append("current membrane closure requires a monovalent anion carrying 1 mol per mol electron")
        number(mem.get("area_m2"),"membrane.area_m2",errors,positive=True)
    for name,permeability in mem.get("neutral_permeability_m_s",{}).items():
        if sm.get(name,{}).get("phase")!="liquid" or sm.get(name,{}).get("charge")!=0: errors.append("neutral permeability requires neutral liquid species: "+name)
        number(permeability,"membrane.permeability."+name,errors,nonnegative=True)
    t=c.get("time",{})
    for key in ("stop_s","output_step_s","reference_max_step_s"): number(t.get(key),"time."+key,errors,positive=True)
    number(c.get("solver",{}).get("relative_tolerance"),"solver.relative_tolerance",errors,positive=True)
    if errors: raise InputError("\n".join(errors))
    # This is a bulk-electroneutral model, not a charge-accumulating double layer.
    # Each half-reaction must generate/remove exactly the charge its electrons carry.
    return c

class Reactor:
    def __init__(self,config):
        self.c=validate(copy.deepcopy(config))
        self.species={s["id"]:s for s in self.c["species"]}
        self.liquid=[n for n,s in self.species.items() if s["phase"]=="liquid"]
        self.gas=[n for n,s in self.species.items() if s["phase"]=="gas"]
        self.keys=[]
        self.symbol={}
        self.initial=[]
        for side in pools(self.c):
            for name in self.liquid:
                self.add(("c",side,name),self.c["compartments"][side]["initial_mol_m3"][name])
                self.add(("out",side,name),0)
            for name in self.gas:self.add(("gas",side,name),0)
        for r in self.c["reactions"]:self.add(("extent",r["id"]),0)
        self.mem_names=[] if self.c.get("cell_topology")=="shared-electrolyte" else [self.c["membrane"]["carrier_species"]]+list(self.c["membrane"].get("neutral_permeability_m_s",{}))
        for name in self.mem_names:self.add(("mem",name),0)
        self.index={k:i for i,k in enumerate(self.keys)}

    def add(self,key,value):
        self.keys.append(key); self.symbol[key]="x"+str(len(self.keys)-1);self.initial.append(float(value))

    def value(self,y,key):return y[self.index[key]]

    def current_densities(self,y):
        if len(y)!=len(self.keys) or not all(math.isfinite(v) for v in y):raise InputError("nonfinite or wrong-size reactor state")
        c=self.c
        currents={}
        eta={}
        for side in SIDES:
            rs=[r for r in c["reactions"] if r["side"]==side]
            if c["kinetic_mode"]=="prescribed-fe":
                eta[side]=None
                for r in rs:currents[r["id"]]=c["current_A"]*r["faradaic_efficiency"]/c["electrode_area_m2"][side]
                continue
            def log_branch(r,e):
                k=r["kinetics"]
                sub=r.get("substrate")
                conc=self.value(y,("c",pool_for(c,side),sub)) if sub else k["substrate_ref_mol_m3"]
                if conc<0:raise InputError("negative substrate concentration; model is outside its physical domain")
                logfactor=0.0
                if sub and k["order"]:
                    if conc==0:return -math.inf
                    logfactor=k["order"]*(math.log(conc)-math.log(k["substrate_ref_mol_m3"]))
                exponent=k["alpha_eff"]*F*(e-k["eta_offset_V"])/(R*c["temperature_K"])
                value=math.log(k["j_ref_A_m2"])+logfactor+exponent
                if not math.isfinite(value):raise InputError("Tafel exponent outside finite numerical domain")
                return value
            target=c["current_A"]/c["electrode_area_m2"][side]
            if not math.isfinite(target):raise InputError("nonfinite current density")
            logtarget=math.log(target)
            def logtotal(e):
                values=[log_branch(r,e) for r in rs];largest=max(values)
                return largest if largest==-math.inf else largest+math.log(sum(math.exp(v-largest) for v in values))
            low,high=-2.0,2.0
            if logtotal(low)>logtarget or logtotal(high)<logtarget:
                raise InputError("Tafel root is not bracketed in [-2,2] V; check inputs")
            for _ in range(50):
                mid=(low+high)/2
                if logtotal(mid)>logtarget:high=mid
                else:low=mid
            eta[side]=(low+high)/2
            for r in rs:currents[r["id"]]=math.exp(log_branch(r,eta[side]))
        rates={r["id"]:currents[r["id"]]*c["electrode_area_m2"][r["side"]]/(r["electron_number"]*F) for r in c["reactions"]}
        transfers={} if c.get("cell_topology")=="shared-electrolyte" else {c["membrane"]["carrier_species"]:c["current_A"]/F}
        for name,permeability in c["membrane"].get("neutral_permeability_m_s",{}).items():
            transfers[name]=permeability*c["membrane"]["area_m2"]*(self.value(y,("c","cathode",name))-self.value(y,("c","anode",name)))
        return currents,rates,eta,transfers

    def rhs(self,y):
        c=self.c; _,rates,_,transfers=self.current_densities(y); dy=[0.0]*len(y)
        for side in pools(c):
            comp=c["compartments"][side]; q=comp["flow_m3_s"];v=comp["volume_m3"]
            for name in self.liquid:
                concentration=self.value(y,("c",side,name))
                flux=q*(comp["feed_mol_m3"][name]-concentration)
                flux+=sum(r["stoichiometry"].get(name,0)*rates[r["id"]] for r in c["reactions"] if pool_for(c,r["side"])==side)
                flux+=(-1 if side=="cathode" else 1)*transfers.get(name,0)
                dy[self.index[("c",side,name)]]=flux/v
                dy[self.index[("out",side,name)]]=q*concentration
            for name in self.gas:dy[self.index[("gas",side,name)]]=sum(r["stoichiometry"].get(name,0)*rates[r["id"]] for r in c["reactions"] if pool_for(c,r["side"])==side)
        for r in c["reactions"]:dy[self.index[("extent",r["id"])]]=rates[r["id"]]
        for name,flux in transfers.items():dy[self.index[("mem",name)]]=flux
        return dy

    def times(self):
        stop=self.c["time"]["stop_s"];step=self.c["time"]["output_step_s"]
        result=[step*i for i in range(int(stop/step)+1)]
        if result[-1]<stop:result.append(stop)
        return result

    def integrate(self):
        y=list(self.initial);time=0;rows=[(0,list(y))]
        for target in self.times()[1:]:
            while time<target:
                h=min(self.c["time"]["reference_max_step_s"],target-time)
                k1=self.rhs(y)
                k2=self.rhs([v+h*d/2 for v,d in zip(y,k1)])
                k3=self.rhs([v+h*d/2 for v,d in zip(y,k2)])
                k4=self.rhs([v+h*d for v,d in zip(y,k3)])
                y=[v+h*(a+2*b+2*d+e)/6 for v,a,b,d,e in zip(y,k1,k2,k3,k4)]
                if time+h==time:raise InputError("reference time step cannot advance floating-point time")
                time+=h
                if not all(math.isfinite(v) for v in y):raise InputError("nonfinite reference state")
                if any(self.value(y,("c",s,n)) < -1e-8 for s in pools(self.c) for n in self.liquid):raise InputError("negative concentration: reduce current or supply valid kinetics/feed; result rejected")
            rows.append((target,list(y)))
        return rows

    def checks(self,time,y):
        currents,rates,eta,transfers=self.current_densities(y);c=self.c
        species_residual=0.0
        element_residual={}
        for side in pools(c):
            comp=c["compartments"][side];q=comp["flow_m3_s"];v=comp["volume_m3"]
            for name in self.liquid:
                lhs=v*self.value(y,("c",side,name))+self.value(y,("out",side,name))
                rhs=v*comp["initial_mol_m3"][name]+q*comp["feed_mol_m3"][name]*time
                rhs+=sum(r["stoichiometry"].get(name,0)*self.value(y,("extent",r["id"])) for r in c["reactions"] if pool_for(c,r["side"])==side)
                if name in self.mem_names:rhs+=(-1 if side=="cathode" else 1)*self.value(y,("mem",name))
                species_residual=max(species_residual,abs(lhs-rhs))
            for name,s in self.species.items():
                if s["phase"]=="electronic":continue
                amount=self.value(y,("gas",side,name)) if s["phase"]=="gas" else v*(self.value(y,("c",side,name))-comp["initial_mol_m3"][name])+self.value(y,("out",side,name))-q*comp["feed_mol_m3"][name]*time
                for el,n in s["atoms"].items():element_residual[el]=element_residual.get(el,0)+amount*n
        charge={side:sum(self.value(y,("c",side,n))*self.species[n]["charge"] for n in self.liquid) for side in pools(c)}
        current_error={side:sum(currents[r["id"]]*c["electrode_area_m2"][side] for r in c["reactions"] if r["side"]==side)-c["current_A"] for side in SIDES}
        fes={r["id"]:currents[r["id"]]*c["electrode_area_m2"][r["side"]]/c["current_A"] for r in c["reactions"]}
        vcell=None if c["kinetic_mode"]!="tafel" else c["anode_reference_potential_V"]-c["cathode_reference_potential_V"]+eta["anode"]+eta["cathode"]+c["current_A"]*c["ohmic_resistance_ohm"]
        charge_time=c["current_A"]*time
        cumulative_error={side:abs(sum(r["electron_number"]*F*self.value(y,("extent",r["id"])) for r in c["reactions"] if r["side"]==side)-charge_time)/max(abs(charge_time),1e-12) for side in SIDES}
        return {"species_balance_max_abs_mol":species_residual,"element_balance_mol_atoms":element_residual,"bulk_charge_mol_charge_m3":charge,"electrode_current_residual_A":current_error,"cumulative_electrode_charge_relative_residual":cumulative_error,"branch_faradaic_efficiency":fes,"reaction_rate_mol_s":rates,"overpotential_V":eta,"cell_voltage_V":vcell,"dc_power_W":None if vcell is None else vcell*c["current_A"],"heat_duty_W":None,"heat_status":"not-computed-no-common-enthalpy-basis","minimum_liquid_concentration_mol_m3":min(self.value(y,("c",s,n)) for s in pools(c) for n in self.liquid)}

    def performance(self,time,y):
        c=self.c
        net={name:sum(r["stoichiometry"].get(name,0)*self.value(y,("extent",r["id"])) for r in c["reactions"]) for name in self.species if name!="electron"}
        active=[c["compartments"][name] for name in pools(c)]
        supplied=sum(p["volume_m3"]*p["initial_mol_m3"].get("EG",0)+p["flow_m3_s"]*p["feed_mol_m3"].get("EG",0)*time for p in active)
        consumed=-net.get("EG",0)
        fe={r["id"]:r["electron_number"]*F*self.value(y,("extent",r["id"]))/ (c["current_A"]*time) if time>0 else None for r in c["reactions"]}
        # Net product generation includes subsequent consumption; no double counting of intermediates.
        other_carbon_feed=any(s["atoms"].get("C",0)>0 and name!="EG" and any(p[field].get(name,0)>0 for p in c["compartments"].values() for field in ("initial_mol_m3","feed_mol_m3")) for name,s in self.species.items())
        carbon_selectivity={name:net.get(name,0)*self.species[name]["atoms"].get("C",0)/(2*consumed) for name in ("GLY","FMA","OX") if name in self.species} if consumed>0 and not other_carbon_feed else None
        return {"net_reaction_generation_mol":net,"eg_net_consumed_mol":consumed,"eg_initial_plus_cumulative_feed_mol":supplied,"eg_cumulative_conversion_initial_plus_feed_basis":consumed/supplied if supplied>0 else None,"carbon_selectivity_net_generation_per_eg_consumed":carbon_selectivity,"carbon_selectivity_status":"EG-only carbon feed, net product generation basis" if carbon_selectivity is not None else "not-applicable-zero-conversion-or-other-carbon-feed","cumulative_branch_faradaic_efficiency":fe,"warning":"scenario calculation; no experimental calibration or separated-product purity claim"}

def owned_path(project_root,path,subfolder=None):
    root=Path(project_root).resolve();base=(root/"AGENTCODING").resolve()
    if not base.is_relative_to(root):raise InputError("AGENTCODING must remain within the selected project")
    p=Path(path);p=(p if p.is_absolute() else root/p).resolve()
    if not p.is_relative_to(base) or p==base:raise InputError("output must remain under the selected project's AGENTCODING")
    if subfolder and not p.is_relative_to((base/subfolder).resolve()):raise InputError("this artifact belongs in AGENTCODING/"+subfolder)
    return p
