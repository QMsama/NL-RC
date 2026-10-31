import importlib.util
from pathlib import Path
import unittest
import numpy as np
from scipy.optimize import brentq

spec=importlib.util.spec_from_file_location('source_curve',Path(__file__).resolve().parents[1]/'exports/identify_source_areas.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class CurveTests(unittest.TestCase):
    def test_known_concave_curve_and_slope_roots(self):
        x=np.linspace(0,1,2001)
        y=(1-np.exp(-3*x))/(1-np.exp(-3))
        fit,info=module.fit_curve(x,y,24)
        self.assertLess(info['max_abs_error_percentage_points'],0.01)
        for slope in [1,2]:
            actual=brentq(lambda v:fit.derivative()(v)-slope,0,1)
            expected=np.log(3/(slope*(1-np.exp(-3))))/3
            self.assertAlmostEqual(actual,expected,delta=0.002)

if __name__=='__main__':unittest.main()
