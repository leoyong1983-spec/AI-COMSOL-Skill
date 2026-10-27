"""Physical applicability and conservation checks for the declared spatial family."""
import math
from reaction_core import InputError, F, finite_tree

SPECIES=("EG","GLY","FMA","OX")
BALANCE_LIMIT=.02


def applicability(c):
    g,o,p,k=c["geometry"],c["operation"],c["properties"],c["kinetics"]
    mean=o["flow_m3_s"]/(g["gap_m"]*g["width_m"])
    dh=2*g["gap_m"]*g["width_m"]/(g["gap_m"]+g["width_m"])
    reynolds=p["density_kg_m3"]*mean*dh/p["viscosity_Pa_s"]
    entrance=.05*reynolds*dh
    if reynolds>=2000:raise InputError("laminar applicability screen failed: Reynolds number >= 2000")
    if c["flow_model"]=="fully-developed-poiseuille" and g["upstream_straight_length_m"]<entrance:
        raise InputError("assumed upstream length is shorter than the 0.05 Re Dh entrance-length estimate")
    if abs(o["temperature_K"]-k["reference_temperature_K"])>1e-6:
        raise InputError("temperature differs from kinetic reference; supply matched kinetic/transport data, not an unmodeled temperature extrapolation")
    return {"mean_velocity_m_s":mean,"hydraulic_diameter_m":dh,"Reynolds_number":reynolds,
            "laminar_screen_limit":2000,"entrance_length_estimate_m":entrance,
            "assumed_upstream_length_m":g["upstream_straight_length_m"],
            "kinetic_reference_temperature_K":k["reference_temperature_K"],
            "anode_EG_reference_mol_m3":k["anode_EG_reference_mol_m3"],
            "scope":"engineering applicability screen; entrance correlation is approximate, not geometric/experimental validation",
            "flow_mode_validation":"declared-fully-developed-family" if c["flow_model"]=="fully-developed-poiseuille" else "experimental-Navier-Stokes-unvalidated"}


def validate_spatial(c):
    finite_tree(c)
    if c.get("schema_version")!=1:raise InputError("schema_version must be 1")
    if c.get("dimension") not in (2,3):raise InputError("dimension must be 2 or 3")
    if c.get("flow_model") not in ("fully-developed-poiseuille","navier-stokes"):raise InputError("explicit supported flow_model required")
    required={
        "geometry":("length_m","gap_m","width_m","upstream_straight_length_m"),
        "operation":("temperature_K","flow_m3_s","feed_EG_mol_m3","support_electrolyte_mol_m3","anode_metal_potential_V_RHE","cathode_metal_potential_V_RHE"),
        "properties":("density_kg_m3","viscosity_Pa_s","conductivity_S_m")+tuple(n+"_diffusivity_m2_s" for n in SPECIES),
        "kinetics":("anode_anchor_current_A_m2","anode_anchor_potential_V_RHE","anode_tafel_V_dec","anode_EG_reference_mol_m3","EG_reaction_order","cathode_anchor_current_A_m2","cathode_anchor_overpotential_V","cathode_tafel_V_dec","reference_temperature_K"),
        "mesh":("axial_elements","gap_elements","span_elements","gap_ratio"),
        "solver":("relative_tolerance",)}
    signed={"anode_metal_potential_V_RHE","cathode_metal_potential_V_RHE","anode_anchor_potential_V_RHE"}
    for group,keys in required.items():
        values=c.get(group)
        if not isinstance(values,dict):raise InputError("missing "+group)
        for key in keys:
            value=values.get(key)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise InputError("nonfinite or missing numeric input "+group+"."+key)
            if key=="EG_reaction_order":
                if value<0:raise InputError("negative reaction order is outside this kinetic family")
            elif key not in signed and value<=0:raise InputError("nonpositive "+key)
    if c["geometry"]["gap_m"]>=min(c["geometry"]["length_m"],c["geometry"]["width_m"]):raise InputError("unsupported non-channel geometry")
    if any(c["mesh"][key]!=int(c["mesh"][key]) or c["mesh"][key]<2 for key in ("axial_elements","gap_elements","span_elements")):raise InputError("mesh counts must be integers >= 2")
    if c["solver"]["relative_tolerance"]>=1:raise InputError("relative_tolerance must be below 1")
    fe=c.get("faradaic_efficiency")
    if not isinstance(fe,dict) or set(fe)!=set(SPECIES[1:]):raise InputError("all carbon branch FE must be explicit")
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in fe.values()) or abs(sum(fe.values())-1)>1e-10:raise InputError("branch FE must be numeric, nonnegative and sum to 1")
    if not isinstance(c.get("assumptions"),list) or not c["assumptions"] or not isinstance(c.get("provenance"),dict) or not c["provenance"]:raise InputError("assumption and provenance records required")
    if not c.get("data_status") or not c.get("case_id"):raise InputError("data_status and case_id required")
    applicability(c)
    return c


def anode_current_density(c, concentration, electrolyte_potential):
    """Independent pointwise oracle; reference concentration never follows inlet."""
    k=c["kinetics"];o=c["operation"]
    return k["anode_anchor_current_A_m2"]*(concentration/k["anode_EG_reference_mol_m3"])**k["EG_reaction_order"]*10**((o["anode_metal_potential_V_RHE"]-electrolyte_potential-k["anode_anchor_potential_V_RHE"])/k["anode_tafel_V_dec"])


def conservation_checks(c,v):
    finite_tree(v,"native metrics")
    required=("Qin_m3_s","Qout_m3_s","I_anode_A","I_cathode_A","EG_surface_min_mol_m3","pressure_drop_Pa","EG_consumption_mol_s")+tuple(n+"_generation_mol_s" for n in SPECIES[1:])+tuple(n+"_net_boundary_out_mol_s" for n in SPECIES)
    if any(isinstance(v.get(k),bool) or not isinstance(v.get(k),(int,float)) for k in required):raise InputError("missing or nonnumeric native metric")
    ia=v["I_anode_A"];q=c["operation"]["flow_m3_s"]
    if min(ia,v["Qin_m3_s"],v["Qout_m3_s"],v["EG_surface_min_mol_m3"],v["pressure_drop_Pa"],v["EG_consumption_mol_s"])<=0:raise InputError("nonphysical current/flow/surface concentration/pressure/consumption")
    fe=c["faradaic_efficiency"]
    expected={"GLY":ia*fe["GLY"]/(4*F),"FMA":2*ia*fe["FMA"]/(6*F),"OX":ia*fe["OX"]/(8*F)}
    expected["EG"]=-(expected["GLY"]+expected["FMA"]/2+expected["OX"])
    checks={"flow_relative_residual":abs(v["Qin_m3_s"]-v["Qout_m3_s"])/q,
            "flow_command_relative_residual":max(abs(v["Qin_m3_s"]-q),abs(v["Qout_m3_s"]-q))/q,
            "current_relative_residual":abs(ia+v["I_cathode_A"])/ia}
    for n in SPECIES:
        generated=-v["EG_consumption_mol_s"] if n=="EG" else v[n+"_generation_mol_s"]
        checks[n+"_faradaic_relative_residual"]=abs(generated-expected[n])/max(abs(expected[n]),1e-15)
        checks[n+"_mass_relative_residual"]=abs(v[n+"_net_boundary_out_mol_s"]-generated)/max(abs(generated),1e-15)
        if abs(generated-expected[n])>1e-8*abs(expected[n])+1e-14:raise InputError("Faraday/electron stoichiometry check failed: "+n)
        if abs(v[n+"_net_boundary_out_mol_s"]-generated)>BALANCE_LIMIT*abs(generated)+1e-14:raise InputError("species flux conservation failed: "+n)
    checks["carbon_relative_residual"]=abs(2*v["EG_net_boundary_out_mol_s"]+2*v["GLY_net_boundary_out_mol_s"]+v["FMA_net_boundary_out_mol_s"]+2*v["OX_net_boundary_out_mol_s"])/(2*v["EG_consumption_mol_s"])
    if any(checks[k]>1e-3 for k in ("flow_relative_residual","flow_command_relative_residual","current_relative_residual")):raise InputError("flow/current/command conservation limits exceeded")
    if checks["carbon_relative_residual"]>BALANCE_LIMIT:raise InputError("carbon conservation limit exceeded")
    return checks
