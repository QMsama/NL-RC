"""Export point21 TN loads by management period and verify the written rasters."""
from pathlib import Path
from datetime import date, timedelta
import argparse
import csv
import hashlib
import json
import sys
import shutil

sys.dont_write_bytecode = True
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
import regional_model as model
import numpy as np
from osgeo import gdal

PERIODS = [
    ('01_basal', '基肥后管理期', '2024-06-01', '2024-06-17'),
    ('02_tillering', '分蘖肥后管理期', '2024-06-18', '2024-07-05'),
    ('03_topdress_to_maturity', '追肥后至成熟期', '2024-07-06', '2024-09-23'),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--reference', required=True)
    args = ap.parse_args()
    root, output, reference = map(Path, (args.root, args.output, args.reference))
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Choose a new output directory to preserve existing results.')
    output.mkdir(parents=True, exist_ok=True)
    version = REPO / '.codex-work/versions/v0.3.1'
    version.mkdir(parents=True, exist_ok=True)
    periods = [(code, label, date.fromisoformat(start), date.fromisoformat(end)) for code,label,start,end in PERIODS]
    seen = []
    sums = {}
    def collect(day, rm, rcc):
        matches = [code for code, _, start, end in periods if start <= day <= end]
        assert len(matches) == 1, (day, matches)
        seen.append(day)
        code = matches[0]
        if code not in sums:
            sums[code] = {'RM': np.zeros_like(rm), 'RCC': np.zeros_like(rcc)}
        sums[code]['RM'] += rm
        sums[code]['RCC'] += rcc
    run = model.run('water_2024_point21', root, daily_load_observer=collect)
    assert seen == [date(2024,6,1)+timedelta(days=i) for i in range(115)]
    template = gdal.Open(str(root/'模型/2026.9.3/模型背景数据/rice_crayfish/Rice_2021.tif'))
    mask = template.ReadAsArray()
    indices = np.flatnonzero((mask == 1) | (mask == 3))
    gt, proj = template.GetGeoTransform(), template.GetProjection()
    area = abs(gt[1]*gt[5]-gt[2]*gt[4])/10000
    assert abs(area-model.PIXEL_HA) < 1e-8
    stats, checks, files = [], {}, []
    ledger = list(csv.DictReader((run/'逐日区域账本.csv').open(encoding='utf-8-sig')))
    for treatment, klass in [('RM',1),('RCC',3)]:
        sel = mask.ravel()[indices] == klass
        valid = mask == klass
        full_file = next(run.glob(f'TN_Load_{treatment}_*.tif'))
        baseline_file = next(reference.glob(f'TN_Load_{treatment}_*.tif'))
        actual = gdal.Open(str(full_file)).ReadAsArray()
        baseline = gdal.Open(str(baseline_file)).ReadAsArray()
        assert np.array_equal(actual,baseline)
        total = sum(sums[code][treatment] for code,_,_,_ in periods)
        assert np.allclose(total[sel],actual[valid],rtol=2e-7,atol=1e-7)
        checks[treatment] = {'season_reference_array_identical': True, 'pixels':int(valid.sum())}
        season_mass = float(total[sel].sum()*area)
        period_arrays = []
        for code,label,start,end in periods + [('00_season','完整稻季',date(2024,6,1),date(2024,9,23))]:
            values = total if code == '00_season' else sums[code][treatment]
            vals = values[sel]
            mass = float(vals.sum()*area)
            ledger_mass = sum(float(row[treatment.lower()+'_load_kg']) for row in ledger if start.isoformat() <= row['date'] <= end.isoformat())
            assert abs(mass-ledger_mass) < 1e-5
            stats.append(dict(treatment=treatment,period=label,start=str(start),end=str(end),days=(end-start).days+1,pixels=int(sel.sum()),mean_kg_ha=float(vals.mean()),total_t=mass/1000,season_percent=100*mass/season_mass,ledger_difference_kg=mass-ledger_mass))
            for unit, factor in [('kg_ha',1.0),('kg_pixel',area)]:
                dest=output/f'TN_{treatment}_{code}_{start:%m%d}_{end:%m%d}_{unit}.tif'
                arr=np.full(mask.shape,-9999,dtype=np.float32)
                arr.ravel()[indices[sel]]=vals*factor
                ds=gdal.GetDriverByName('GTiff').Create(str(dest),template.RasterXSize,template.RasterYSize,1,gdal.GDT_Float32,options=['COMPRESS=LZW'])
                ds.SetGeoTransform(gt);ds.SetProjection(proj)
                ds.SetMetadata({'scenario':'water_2024_point21','period':label,'start':str(start),'end':str(end),'aggregation':'sum of daily TN outlet loads','unit':unit,'last_irrigation_date':'2024-08-21','irrigation_cutoff':'2024-08-22'})
                band=ds.GetRasterBand(1);band.SetNoDataValue(-9999);band.SetUnitType('kg/ha' if unit=='kg_ha' else 'kg/pixel');band.WriteArray(arr);band.ComputeStatistics(False)
                band=None;ds=None
                read=gdal.Open(str(dest));written=read.ReadAsArray()
                assert read.GetGeoTransform()==gt and read.GetProjection()==proj
                assert np.array_equal(written,arr)
                assert np.array_equal(written!=-9999,valid)
                assert np.isfinite(written[valid]).all() and (written[valid]>=0).all()
                files.append({'file':dest.name,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'bytes':dest.stat().st_size})
                if unit=='kg_ha' and code!='00_season':period_arrays.append(written[valid].astype(np.float64))
                if unit=='kg_ha' and code=='00_season': season_written=written[valid]
        delta=np.sum(period_arrays,axis=0)-season_written
        assert np.allclose(np.sum(period_arrays,axis=0),season_written,rtol=3e-7,atol=2e-6)
        checks[treatment]['written_stage_sum_max_error_kg_ha']=float(np.max(np.abs(delta)))
    for treatment in ['RM', 'RCC']:
        drainage = next(run.glob(f'TN_Drain_{treatment}_*.tif'))
        ref_drainage = next(reference.glob(f'TN_Drain_{treatment}_*.tif'))
        assert np.array_equal(gdal.Open(str(drainage)).ReadAsArray(),gdal.Open(str(ref_drainage)).ReadAsArray())
        target = output / drainage.name
        shutil.copy2(drainage, target)
        assert target.read_bytes() == drainage.read_bytes()
        checks[treatment]['season_drainage_reference_array_identical'] = True
        files.append({'file':target.name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'bytes':target.stat().st_size})
    with (version/'stage_statistics.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=stats[0].keys());writer.writeheader();writer.writerows(stats)
    verification={'passed':True,'scenario':'water_2024_point21','dates':len(seen),'periods':PERIODS,'checks':checks,'run_id':run.name,'reference_run_id':reference.name,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'files':files,'output':str(output),'float32_sum_tolerance':'rtol 3e-7, atol 2e-6 kg/ha'}
    (version/'verification.json').write_text(json.dumps(verification,ensure_ascii=False,indent=2),encoding='utf-8')
    for f in ['summary.json','逐日区域账本.csv']:
        shutil.copy2(run/f,version/f)
    print(json.dumps({'passed':True,'output':str(output),'run':str(run),'checks':checks},ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
