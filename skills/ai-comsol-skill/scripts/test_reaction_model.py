"""Behavioral and numerical invariants; no licensed solver needed for this suite."""
import copy
import json
import math
import unittest
from pathlib import Path
from reaction_core import Reactor, InputError, F, validate, owned_path

ASSETS=Path(__file__).resolve().parents[1]/"assets"
def config(name="demo_eg_her.json"):
    return json.loads((ASSETS/name).read_text(encoding="utf-8"))

class ReactionModelTests(unittest.TestCase):
    def test_all_active_half_reactions_conserve_elements_and_charge(self):
        validate(config()); validate(config("literature_eg_equivalent.json"))

    def test_wrong_electron_count_is_rejected(self):
        c=config();c["reactions"][1]["electron_number"]=2
        with self.assertRaisesRegex(InputError,"electron count"):validate(c)

    def test_unknown_topology_and_membrane_in_shared_liquid_are_rejected(self):
        c=config();c["cell_topology"]="unknown"
        with self.assertRaisesRegex(InputError,"cell_topology"):validate(c)
        c=config("literature_eg_equivalent.json");c["membrane"]["neutral_permeability_m_s"]={"EG":1e-7}
        with self.assertRaisesRegex(InputError,"no membrane permeability"):validate(c)

    def test_cumulative_conversion_and_net_carbon_selectivity(self):
        c=config("literature_eg_equivalent.json");m=Reactor(c);t,y=m.integrate()[-1];p=m.performance(t,y)
        self.assertAlmostEqual(p["eg_net_consumed_mol"],.475e-3,places=12)
        self.assertAlmostEqual(p["eg_cumulative_conversion_initial_plus_feed_basis"],.0095,places=12)
        self.assertAlmostEqual(p["carbon_selectivity_net_generation_per_eg_consumed"]["GLY"],1,places=12)
        self.assertAlmostEqual(p["cumulative_branch_faradaic_efficiency"]["a_eg"],.92,places=12)

    def test_wrong_carbon_stoichiometry_is_rejected(self):
        c=config();c["reactions"][2]["stoichiometry"]["FMA"]=1
        with self.assertRaisesRegex(InputError,"element imbalance"):validate(c)

    def test_fe_is_checked_per_electrode_not_paired_total(self):
        c=config("literature_eg_equivalent.json")
        self.assertAlmostEqual(sum(r["faradaic_efficiency"] for r in c["reactions"]),2)
        validate(c)
        c["reactions"][1]["faradaic_efficiency"]-=.01
        with self.assertRaisesRegex(InputError,"unresolved current"):validate(c)

    def test_shared_batch_has_no_artificial_fresh_feed_and_matches_faradaic_amount(self):
        c=config("literature_eg_equivalent.json");m=Reactor(c);t,y=m.integrate()[-1]
        expected=c["current_A"]*.92*t/(4*F)
        self.assertAlmostEqual(m.value(y,("c","cell","GLY"))*c["compartments"]["cell"]["volume_m3"],expected,places=12)
        self.assertAlmostEqual(expected,.475e-3,places=12)
        self.assertEqual(m.value(y,("out","cell","EG")),0)
        self.assertLess(max(abs(x) for x in m.checks(t,y)["element_balance_mol_atoms"].values()),1e-10)

    def test_asymmetric_2e_4e_reaction_products_share_current(self):
        c=config("literature_eg_equivalent.json")
        c["reactions"]=[r for r in c["reactions"] if r["id"]!="a_oer"]
        for r in c["reactions"]:r["faradaic_efficiency"]=1
        m=Reactor(c);_,rates,_,_=m.current_densities(m.initial)
        self.assertAlmostEqual(rates["c_her"]/rates["a_eg"],2,places=12)

    def test_cstr_analytic_solution_for_fixed_fe(self):
        c=config("demo_kinetic.json");c["kinetic_mode"]="prescribed-fe"
        c["membrane"]["neutral_permeability_m_s"]={}
        c["time"].update(stop_s=60,output_step_s=60,reference_max_step_s=.2)
        m=Reactor(c);t,y=m.integrate()[-1]
        p=c["compartments"]["cathode"];rate=c["current_A"]*.9/(2*F)
        expected=p["feed_mol_m3"]["FF"]-rate/p["flow_m3_s"]*(1-math.exp(-p["flow_m3_s"]*t/p["volume_m3"]))
        self.assertAlmostEqual(m.value(y,("c","cathode","FF")),expected,places=8)

    def test_concentration_and_current_closure_in_kinetic_model(self):
        c=config();c["time"].update(stop_s=20,output_step_s=10,reference_max_step_s=.5)
        m=Reactor(c);t,y=m.integrate()[-1];checks=m.checks(t,y)
        self.assertLess(max(abs(v) for v in checks["electrode_current_residual_A"].values()),1e-8)
        self.assertLess(checks["species_balance_max_abs_mol"],1e-10)
        self.assertLess(max(abs(v) for v in checks["bulk_charge_mol_charge_m3"].values()),1e-8)

    def test_unbalanced_feed_is_not_silently_fixed(self):
        c=config();c["compartments"]["anode"]["feed_mol_m3"]["K"]=0
        with self.assertRaisesRegex(InputError,"not electroneutral"):validate(c)

    def test_unfilled_project_input_is_blocked(self):
        with self.assertRaises(InputError):validate(config("project_input_template.json"))

    def test_no_external_artifact_path(self):
        c=config();root=Path.cwd().resolve()
        with self.assertRaisesRegex(InputError,"AGENTCODING"):owned_path(root,root/"outside-candidate.json")
        with self.assertRaises(InputError):owned_path(root,root/"AGENTCODING"/".."/"escape.json")

    def test_substrate_depletion_is_rejected_instead_of_clipped(self):
        c=config("literature_eg_equivalent.json");c["current_A"]=500;c["time"].update(stop_s=600,output_step_s=60,reference_max_step_s=1)
        with self.assertRaisesRegex(InputError,"negative concentration"):Reactor(c).integrate()

if __name__=="__main__":unittest.main(verbosity=2)
