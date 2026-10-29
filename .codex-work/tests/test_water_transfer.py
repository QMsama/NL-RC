import sys,unittest
from pathlib import Path
from datetime import date
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from water_balance import rm_transfer_thresholds
class Transfer(unittest.TestCase):
 def test_high_water_and_last_day(self):
  original=(30.,50.,100.,False)
  self.assertEqual(rm_transfer_thresholds(date(2024,8,1),original,21),(80.,90.,100.,False))
  self.assertEqual(rm_transfer_thresholds(date(2024,8,10),original,21),(80.,90.,100.,False))
  self.assertEqual(rm_transfer_thresholds(date(2024,8,11),original,21),original)
  self.assertEqual(rm_transfer_thresholds(date(2024,8,21),original,21),(50.,50.,100.,False))
  self.assertEqual(rm_transfer_thresholds(date(2024,8,22),original,21),(0.,50.,100.,False))
  self.assertEqual(rm_transfer_thresholds(date(2024,9,1),(0.,0.,40.,False),21),(0.,0.,40.,False))
 def test_controls_and_drying(self):
  for day in [21,24]:
   self.assertEqual(rm_transfer_thresholds(date(2024,7,3),(0.,0.,0.,True),day),(0.,0.,0.,True))
   self.assertEqual(rm_transfer_thresholds(date(2024,8,day+1),(30.,50.,100.,False),day)[0],0.)
  with self.assertRaises(ValueError):rm_transfer_thresholds(date(2024,8,1),(30.,50.,100.,False),25)
if __name__=='__main__':unittest.main()
