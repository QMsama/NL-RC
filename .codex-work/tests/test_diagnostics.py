import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import regional_model as m
from water_balance import regional_ditch_step
class Diagnostics(unittest.TestCase):
 def test_expected_requires_both_treatments(self):
  actual={'RM':{'mean_kg_ha':2.,'total_t':3.},'RCC':{'mean_kg_ha':4.,'total_t':5.}}
  self.assertIsNone(m.check_expected(actual,None)[1])
  for expected in [{},{'RM':{},'RCC':{}},{'RM':{'mean_kg_ha':2.}},{'RM':{'mean_kg_ha':2.},'RCC':{'total_t':5.}}]:
   with self.subTest(expected=expected),self.assertRaises(ValueError):m.check_expected(actual,expected)
  self.assertTrue(m.check_expected(actual,actual)[1])
 def test_float32_ditch_diagnostic(self):
  v=np.array([1000.,6000.]);rain=np.array([56.63,96.],dtype=np.float32);loss=np.array([3.4,4.8],dtype=np.float32);dr=np.array([10.,15.]);mask=np.array([True,True]);length=m.DITCH_L;area=m.A_TOTAL
  water,outlet=regional_ditch_step(v,rain,loss,dr,mask,length,area,m.V_MAX)
  residual=m.ditch_closure_residual(v,rain,loss,dr,water,outlet)
  self.assertLess(np.max(np.abs(residual)),1e-8)
if __name__=='__main__':unittest.main()
