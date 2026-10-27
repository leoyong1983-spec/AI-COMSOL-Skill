"""Regression tests for reproduced scientific false positives and parameter binding."""
import copy
import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path
from test_reaction_model import config as reaction_config
from test_spatial_inputs import config as spatial_config
from reaction_core import InputError, Reactor, F, validate
from spatial_contract import validate_spatial, anode_current_density, conservation_checks
from spatial_review import compare_data, accepted_case, main as review_main
from spatial_reactor import field_checks, generate
from numerical_checks import verify_readback


def balanced_metrics(c):
    i=.073
    fe=c['faradaic_efficiency'];q=c['operation']['flow_m3_s']
    rates={'GLY':i*fe['GLY']/(4*F),'FMA':2*i*fe['FMA']/(6*F),'OX':i*fe['OX']/(8*F)}
    v={'I_anode_A':i,'I_cathode_A':-i,'Qin_m3_s':q,'Qout_m3_s':q,'EG_surface_min_mol_m3':880,'pressure_drop_Pa':11.2,'j_anode_max_A_m2':740}
    v['EG_consumption_mol_s']=rates['GLY']+rates['FMA']/2+rates['OX']
    v['EG_net_boundary_out_mol_s']=-v['EG_consumption_mol_s']
    for n,r in rates.items():v[n+'_generation_mol_s']=r;v[n+'_net_boundary_out_mol_s']=r
    return v


class ModelAuditTests(unittest.TestCase):
    def test_nonfinite_stoichiometry_and_atoms_rejected(self):
        for key in ('stoichiometry','atoms'):
            c=reaction_config()
            if key=='stoichiometry':c['reactions'][0][key]['electron']=float('nan')
            else:c['species'][0][key]={'H':float('nan')}
            with self.assertRaisesRegex(InputError,'nonfinite'):validate(c)

    def test_unused_inventory_cannot_dilute_reported_conversion(self):
        c=reaction_config('literature_eg_equivalent.json')
        c['compartments']['unused']=copy.deepcopy(c['compartments']['cell'])
        with self.assertRaisesRegex(InputError,'unused inventories'):validate(c)

    def test_charge_time_integral_matches_branch_electron_extents(self):
        c=reaction_config('literature_eg_equivalent.json');m=Reactor(c);t,y=m.integrate()[-1]
        self.assertLess(max(m.checks(t,y)['cumulative_electrode_charge_relative_residual'].values()),1e-10)
        y[m.index[('extent','a_eg')]]*=.5
        self.assertGreater(m.checks(t,y)['cumulative_electrode_charge_relative_residual']['anode'],.4)

    def test_high_tafel_slope_does_not_overflow_root_bracket(self):
        c=reaction_config('demo_kinetic.json')
        for r in c['reactions']:r['kinetics']['alpha_eff']=10
        m=Reactor(c)
        self.assertLess(max(abs(x) for x in m.checks(0,m.initial)['electrode_current_residual_A'].values()),1e-8)

    def test_missing_trace_product_is_not_a_successful_readback(self):
        with self.assertRaisesRegex(InputError,'save/reopen mismatch'):verify_readback([[4.75e-9]],[[0]])

    def test_readback_allows_roundoff_and_rejects_nonfinite(self):
        self.assertLess(verify_readback([[0,1000,4.75e-9]],[[1e-15,1000+1e-7,4.75e-9+1e-16]])['save_reopen_max_tolerance_ratio'],1)
        for after in ([[float('nan')]], [[float('inf')]], []):
            with self.assertRaises(InputError):verify_readback([[1]],after)

    def test_feed_change_preserves_literature_kinetic_reference(self):
        c=spatial_config();half=copy.deepcopy(c);half['operation']['feed_EG_mol_m3']=500
        self.assertAlmostEqual(anode_current_density(c,800,.01),anode_current_density(half,800,.01))
        self.assertAlmostEqual(anode_current_density(c,500,0)/anode_current_density(c,1000,0),.5**.3)
        java,_=generate(half,Path('case-output'))
        self.assertIn('cEG/cEGref',java)
        self.assertNotIn('cEG/cin',java)

    def test_export_catches_old_feed_normalized_law(self):
        c=spatial_config();c['operation']['feed_EG_mol_m3']=500
        # The old law would give 1000 A/m2 at 500 mol/m3 and zero electrolyte potential.
        with patch('spatial_reactor.decode',return_value='0 0 500 0 0 0 0 0 0 1000'):
            with self.assertRaisesRegex(InputError,'fixed-reference'):field_checks('unused',2,500,c)

    def test_temperature_and_flow_applicability_are_enforced(self):
        c=spatial_config();c['operation']['temperature_K']=333.15
        with self.assertRaisesRegex(InputError,'temperature'):validate_spatial(c)
        c=spatial_config();c['operation']['flow_m3_s']=1e-4
        with self.assertRaisesRegex(InputError,'Reynolds'):validate_spatial(c)
        c=spatial_config();c['geometry']['upstream_straight_length_m']=.001
        with self.assertRaisesRegex(InputError,'upstream'):validate_spatial(c)

    def test_declared_zero_order_and_zero_product_branch_are_valid(self):
        c=spatial_config();c['kinetics']['EG_reaction_order']=0;c['faradaic_efficiency']={'GLY':1,'FMA':0,'OX':0}
        validate_spatial(c);conservation_checks(c,balanced_metrics(c))

    def test_equal_in_out_wrong_flow_is_rejected(self):
        c=spatial_config();v=balanced_metrics(c);v['Qin_m3_s']*=2;v['Qout_m3_s']*=2
        with self.assertRaisesRegex(InputError,'command'):conservation_checks(c,v)

    def test_formate_two_carbon_split_is_checked_independently(self):
        c=spatial_config();v=balanced_metrics(c);v['FMA_generation_mol_s']/=2;v['FMA_net_boundary_out_mol_s']/=2
        with self.assertRaisesRegex(InputError,'Faraday'):conservation_checks(c,v)

    def test_mesh_rejects_failed_nonfinite_and_inactive_refinement(self):
        c=spatial_config();fine=copy.deepcopy(c);fine['mesh']['axial_elements']*=2
        card={'status':'native-spatial-solve-and-readback-passed','metrics':balanced_metrics(c)}
        bad=copy.deepcopy(card);bad['status']='rejected'
        with self.assertRaisesRegex(InputError,'rejected'):compare_data(bad,card,c,fine)
        bad=copy.deepcopy(card);bad['metrics']['I_anode_A']=float('nan')
        with self.assertRaisesRegex(InputError,'nonfinite'):compare_data(bad,card,c,fine)
        unused=copy.deepcopy(c);unused['mesh']['span_elements']*=2
        with self.assertRaisesRegex(InputError,'active dimensions'):compare_data(card,card,c,unused)
        with self.assertRaisesRegex(InputError,'active dimensions'):compare_data(card,card,fine,c)

    def test_mesh_fails_when_key_metric_changes_above_limit(self):
        c=spatial_config();fine=copy.deepcopy(c);fine['mesh']['gap_elements']*=2
        card={'status':'native-spatial-solve-and-readback-passed','metrics':balanced_metrics(c)}
        bad=copy.deepcopy(card);bad['metrics']['I_anode_A']*=1.1
        self.assertEqual(compare_data(card,bad,c,fine)['status'],'further-refinement-required')

    def test_case_review_rejects_modified_config_artifact(self):
        c=spatial_config();card={'status':'native-spatial-solve-and-readback-passed','interface_version':'spatial-2','config_sha256':'wrong'}
        with patch('spatial_review.read',side_effect=[card,c]),patch('spatial_review.Path.is_file',return_value=True),patch('spatial_review.Path.read_bytes',return_value=b'edited'):
            with self.assertRaisesRegex(InputError,'hash mismatch'):accepted_case('unused')

    def test_review_cli_returns_failure_and_retains_rejected_review(self):
        args=SimpleNamespace(command='mesh',coarse='a',fine='b',case=None,out='review.json',project_root=None)
        result={'status':'further-refinement-required'}
        with patch('spatial_review.argparse.ArgumentParser.parse_args',return_value=args),patch('spatial_review.compare',return_value=result),patch('spatial_review.write') as writer,patch('builtins.print'):
            with self.assertRaises(SystemExit) as ex:review_main()
            self.assertEqual(ex.exception.code,2)
            writer.assert_called_once_with('review.json',result,None)


if __name__=='__main__':unittest.main(verbosity=2)
