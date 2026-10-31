from pathlib import Path
import json, os
os.environ['QT_QPA_PLATFORM']='offscreen'
from qgis.core import QgsApplication,QgsRasterLayer
from qgis.analysis import QgsRasterCalculator,QgsRasterCalculatorEntry
from osgeo import gdal
import numpy as np
app=QgsApplication([],False);app.initQgis()
root=Path.cwd();base=root/'.codex-work/nl_rc_2024/upload_original/.codex-work/versions/v0.4.0'
fit=json.loads((base/'thresholds_and_fit.json').read_text(encoding='utf-8'))
results={}
for t in ['RM','RCC']:
    inp=root/f'模型/2026.9.3/面上推广/2024/point21_三期累计TN流失/TN_{t}_00_season_0601_0923_kg_ha.tif'
    lyr=QgsRasterLayer(str(inp),t);assert lyr.isValid()
    e=QgsRasterCalculatorEntry();e.ref=t+'@1';e.raster=lyr;e.bandNumber=1
    low=fit[t]['secondary_threshold_kg_ha'];high=fit[t]['critical_threshold_kg_ha']
    expression=f'1 + ("{t}@1" >= {low!r}) + ("{t}@1" >= {high!r})'
    out=root/f'.codex-work/nl_rc_2024/qgis_source_{t}.tif'
    calc=QgsRasterCalculator(expression,str(out),'GTiff',lyr.extent(),lyr.width(),lyr.height(),[e]);result=calc.processCalculation();assert int(result)==0
    ds=gdal.Open(str(out));arr=ds.ReadAsArray();expected=gdal.Open(str(root/f'模型/2026.9.3/面上推广/2024/point21_关键源区/TN_source_class_{t}_2024_point21.tif')).ReadAsArray()
    valid=expected>0
    assert np.array_equal(arr[valid],expected[valid])
    assert np.all(arr[~valid]==ds.GetRasterBand(1).GetNoDataValue())
    results[t]={'expression':expression,'valid_pixel_class_identical':True,'nodata_preserved':True}
(base/'qgis_verification.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(results,ensure_ascii=False))
