import sys,unittest,copy
from pathlib import Path
from datetime import date
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import water_quality as wq

def config():
 return {t:[dict(date='2024-06-01',response='2024-06-02',stage='basal',frac=.58),dict(date='2024-06-18',response='2024-06-19',stage='tillering',frac=.18),dict(date='2024-07-06',response='2024-07-07',stage='tillering',frac=.24)] for t in ['RM','RCC']}
class Events(unittest.TestCase):
 def load(self,c):return wq.validate_regional_events(c,date(2024,6,1),date(2024,9,23))
 def test_allocation_conservation(self):
  ev=self.load(config())
  for t in ev:
   np.testing.assert_allclose(sum(np.array([0.,100.,221.])*e['frac'] for e in ev[t]),[0,100,221],rtol=0,atol=1e-12)
 def test_reject_invalid(self):
  for mode in ['sum','negative','nan','missing','duplicate','reverse','outside','early_response','late_response','stage']:
   c=config()
   if mode=='sum':c['RM'][0]['frac']=.6
   if mode=='negative':c['RM'][0]['frac']=-.58
   if mode=='nan':c['RM'][0]['frac']=float('nan')
   if mode=='missing':del c['RCC']
   if mode=='duplicate':c['RM'][1]['date']='2024-06-01'
   if mode=='reverse':c['RM'].reverse()
   if mode=='outside':c['RM'][0]['date']='2024-05-31'
   if mode=='early_response':c['RM'][0]['response']='2024-05-31'
   if mode=='late_response':c['RM'][-1]['response']='2024-09-24'
   if mode=='stage':c['RM'][0]['stage']='unknown'
   with self.subTest(mode=mode), self.assertRaises(ValueError):self.load(c)
 def test_response_boundary_and_event_switch(self):
  ev=self.load(config());peaks={(t,j):np.array([10.+j]) for t in ev for j in range(3)}
  before=wq.regional_concentration('RM',date(2024,6,1),1,ev,peaks)
  first=wq.regional_concentration('RM',date(2024,6,2),1,ev,peaks)
  next_event=wq.regional_concentration('RM',date(2024,6,19),1,ev,peaks)
  self.assertAlmostEqual(before[0],wq.BACKGROUND)
  self.assertAlmostEqual(first[0],wq.BACKGROUND+10)
  self.assertAlmostEqual(next_event[0],wq.BACKGROUND+11)
if __name__=='__main__':unittest.main()
