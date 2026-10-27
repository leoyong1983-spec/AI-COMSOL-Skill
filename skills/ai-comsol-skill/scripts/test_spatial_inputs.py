"""Behavioral guards for scientific inputs and false-positive spatial acceptance."""
import unittest,copy,json
from pathlib import Path
from unittest.mock import patch
from spatial_reactor import load_config,field_checks
from spatial_review import compare
from reaction_core import InputError
CFG=Path(__file__).resolve().parents[1]/"assets/spatial_pdnimoo4.json"
def config():return json.loads(CFG.read_text(encoding="utf-8"))
def validate(c):
 with patch("spatial_reactor.Path.read_text",return_value=json.dumps(c)):return load_config(CFG)
class SpatialGuards(unittest.TestCase):
 def test_nonfinite_transport_is_rejected(self):
  c=config();c["properties"]["EG_diffusivity_m2_s"]=float("nan")
  with self.assertRaisesRegex(InputError,"nonfinite"):validate(c)
 def test_unassigned_branch_current_is_rejected(self):
  c=config();c["faradaic_efficiency"]["OX"]=0
  with self.assertRaisesRegex(InputError,"sum"):validate(c)
 def test_unknown_flow_physics_is_rejected(self):
  c=config();c["flow_model"]="unspecified"
  with self.assertRaisesRegex(InputError,"flow_model"):validate(c)
 def test_wrong_project_cannot_be_silently_used(self):
  c=config();c["project_root"]=str(Path.cwd().parent/"different-workspace")
  with self.assertRaisesRegex(InputError,"ownership"):validate(c)
 def test_overshooting_substrate_field_cannot_be_accepted(self):
  with patch("spatial_reactor.decode",return_value="0 0 1200 0 0 0 0 0 0 0"):
   with self.assertRaisesRegex(InputError,"maximum"):field_checks("unused",2,1000)
 def test_negative_product_is_rejected_without_clipping(self):
  with patch("spatial_reactor.decode",return_value="0 0 900 -1 0 0 0 0 0 0"):
   with self.assertRaisesRegex(InputError,"negative"):field_checks("unused",2,1000)
 def test_mesh_comparison_rejects_changed_operating_condition(self):
  a=config();b=copy.deepcopy(a);b["operation"]["flow_m3_s"]*=2
  card={"metrics":{}}
  with patch("spatial_review.read",side_effect=[card,card,a,b]):
   with self.assertRaisesRegex(InputError,"identical physical"):compare("coarse","fine")
if __name__=="__main__":unittest.main(verbosity=2)
