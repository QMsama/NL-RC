"""Fit constrained load-area splines and export spatial TN source classes."""
from pathlib import Path
import argparse, csv, json, hashlib, sys
sys.dont_write_bytecode=True
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize, LinearConstraint, brentq
from osgeo import gdal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def write_csv(path, rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def fit_curve(x,y,count):
    inner=np.unique(np.r_[np.linspace(0,1,count+2)[1:-1],0.001,0.005,0.01,0.02,0.98,0.995])
    knots=np.r_[np.zeros(4),inner,np.ones(4)]
    n=len(knots)-4
    B=BSpline.design_matrix(x,knots,3).toarray()
    # Derivative control points are nonnegative and nonincreasing.
    D=np.zeros((n-1,n))
    for i in range(n-1):D[i,i:i+2]=[-3,3]/(knots[i+4]-knots[i+1])
    C=np.vstack([D,D[:-1]-D[1:]])
    init=np.array([knots[i+1:i+4].mean() for i in range(n)])
    bounds=[(0.,1.)]*n;bounds[0]=(0.,0.);bounds[-1]=(1.,1.)
    result=minimize(lambda c:np.mean((B@c-y)**2),init,jac=lambda c:2*B.T@(B@c-y)/len(y),bounds=bounds,constraints=[LinearConstraint(C,0,np.inf)],method='SLSQP',options={'maxiter':2000,'ftol':1e-13})
    if not result.success:raise RuntimeError(result.message)
    spline=BSpline(knots,result.x,3)
    grid=np.linspace(0,1,10001);der=spline.derivative()(grid)
    assert der.min()>-1e-7 and np.diff(der).max()<1e-7
    assert abs(spline(0))<1e-10 and abs(spline(1)-1)<1e-10
    err=spline(x)-y
    diagnostics={'curve_R2':float(1-np.sum(err**2)/np.sum((y-y.mean())**2)),'curve_RMSE_percentage_points':float(np.sqrt(np.mean(err**2))*100),'max_abs_error_percentage_points':float(np.max(np.abs(err))*100),'knots':knots.tolist(),'coefficients':result.x.tolist(),'degree':3,'optimizer':result.message}
    return spline,diagnostics

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--output',required=True);ap.add_argument('--records',required=True);args=ap.parse_args()
    source,out,records=map(Path,(args.input,args.output,args.records));out.mkdir(parents=True,exist_ok=True);records.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,'font.size':10})
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    mapfig,mapaxes=plt.subplots(1,2,figsize=(10,6),layout='constrained')
    summary={};stats=[];sens=[];audit={'inputs':[],'checks':{},'method':'constrained cubic B-spline least squares; slope roots 1 and 2; within-treatment normalization'}
    for col,treatment in enumerate(['RM','RCC']):
        path=source/f'TN_{treatment}_00_season_0601_0923_kg_ha.tif'
        masspath=source/f'TN_{treatment}_00_season_0601_0923_kg_pixel.tif'
        ds=gdal.Open(str(path));arr=ds.ReadAsArray();nd=ds.GetRasterBand(1).GetNoDataValue();valid=arr!=nd
        vals=arr[valid].astype(float);assert np.isfinite(vals).all() and (vals>=0).all()
        massds=gdal.Open(str(masspath));mass=massds.ReadAsArray()[valid].astype(float)
        gt=ds.GetGeoTransform();area=abs(gt[1]*gt[5]-gt[2]*gt[4])/10000
        assert ds.GetProjection()==massds.GetProjection() and gt==massds.GetGeoTransform()
        assert np.allclose(vals*area,mass,rtol=2e-7,atol=1e-6)
        # Authoritative load intensity determines ranking; mass derived consistently from area.
        order=np.argsort(-vals,kind='stable');sortedvals=vals[order];n=len(vals)
        xp=np.arange(n+1)/n;yp=np.r_[0,np.cumsum(sortedvals)/sortedvals.sum()]
        x=np.unique(np.r_[np.linspace(0,1,2001),0.001,0.005,0.01]);y=np.interp(x,xp,yp)
        fits={}
        for count in [16,24,32]:
            sp,diag=fit_curve(x,y,count);roots={};cuts={}
            for level in [1,2]:
                derivative=sp.derivative()
                if derivative(0)<=level:root=0.
                elif derivative(1)>=level:root=1.
                else:root=brentq(lambda z:float(derivative(z)-level),0,1,xtol=1e-13)
                h=int(np.clip(round(root*n),0,n))
                cut=float(sortedvals[h-1]) if h else float(np.nextafter(np.float32(sortedvals[0]),np.float32(np.inf)))
                roots[level]=root;cuts[level]=cut
            assert cuts[2]>=cuts[1]
            classes=np.ones(n,dtype=np.uint8);classes[vals>=cuts[1]]=2;classes[vals>=cuts[2]]=3
            fits[count]=(sp,diag,roots,cuts,classes)
            sel=classes==3
            sens.append({'treatment':treatment,'interior_knots':count,'critical_threshold_kg_ha':cuts[2],'secondary_threshold_kg_ha':cuts[1],'critical_area_pct':float(sel.mean()*100),'critical_load_pct':float(vals[sel].sum()/vals.sum()*100),'curve_R2':diag['curve_R2']})
        sp,diag,roots,cuts,classes=fits[24]
        full=np.zeros(arr.shape,dtype=np.uint8);full[valid]=classes
        dest=out/f'TN_source_class_{treatment}_2024_point21.tif'
        target=gdal.GetDriverByName('GTiff').Create(str(dest),ds.RasterXSize,ds.RasterYSize,1,gdal.GDT_Byte,['COMPRESS=LZW'])
        target.SetGeoTransform(gt);target.SetProjection(ds.GetProjection());target.SetMetadata({'class_1':'other valid crop area','class_2':'sub-critical source candidate','class_3':'critical source candidate','normalization':'within treatment','scenario':'water_2024_point21','threshold_1':repr(cuts[1]),'threshold_2':repr(cuts[2])})
        band=target.GetRasterBand(1);band.SetNoDataValue(0);band.WriteArray(full)
        colors=gdal.ColorTable()
        for k,c in [(0,(255,255,255,0)),(1,(223,234,240,255)),(2,(252,189,90,255)),(3,(190,46,45,255))]:colors.SetColorEntry(k,c)
        band.SetRasterColorTable(colors);band=None;target=None
        rd=gdal.Open(str(dest));assert np.array_equal(rd.ReadAsArray(),full);assert rd.GetGeoTransform()==gt and rd.GetProjection()==ds.GetProjection()
        for k,label in [(1,'其他有效种植区域'),(2,'次关键源区'),(3,'关键源区候选')]:
            selected=classes==k;stats.append({'treatment':treatment,'class':k,'label':label,'pixels':int(selected.sum()),'area_ha':float(selected.sum()*area),'area_pct':float(selected.mean()*100),'load_t':float(vals[selected].sum()*area/1000),'load_pct':float(vals[selected].sum()/vals.sum()*100)})
        boundaries=[]
        for level in [1,2]:
            selected=vals>=cuts[level];actual=vals[selected].sum()/vals.sum();pred=float(sp(roots[level]))
            boundaries.append({'slope':level,'fitted_area_pct':roots[level]*100,'fitted_load_pct':pred*100,'actual_area_pct':selected.mean()*100,'actual_load_pct':actual*100,'load_difference_percentage_points':(actual-pred)*100,'relative_load_difference_pct':(actual-pred)/pred*100})
        summary[treatment]={'pixels':n,'mean_kg_ha':float(vals.mean()),'critical_threshold_kg_ha':cuts[2],'secondary_threshold_kg_ha':cuts[1],'fit':diag,'boundaries':boundaries,'sensitivity_agreement_pct':{str(k):float((fits[k][4]==classes).mean()*100) for k in [16,32]}}
        write_csv(records/f'curve_{treatment}.csv',[{'area_pct':float(a*100),'actual_cumulative_load_pct':float(b*100),'fitted_cumulative_load_pct':float(sp(a)*100),'slope':float(sp.derivative()(a))} for a,b in zip(x,y)])
        ax=axes[0,col];ax.plot(x*100,y*100,label='经验累计曲线',color='#333333',lw=2);ax.plot(x*100,sp(x)*100,label='约束样条拟合',color='#227aab',ls='--');ax.plot([0,100],[0,100],color='#bbbbbb',lw=.7);ax.set(xlim=(0,100),ylim=(0,100),xlabel='累计面积（%）',ylabel='累计TN负荷（%）',title=treatment+' 全季负荷—面积曲线');ax.legend(loc='lower right')
        ax=axes[1,col];ax.plot(x*100,sp.derivative()(x),color='#227aab');ax.set(xlim=(0,100),xlabel='累计面积（%）',ylabel='拟合曲线斜率',title='拟合斜率及空间分级界限')
        for level in [1,2]:
            ax.axhline(level,color='#777777',ls=':',lw=.8);ax.axvline(roots[level]*100,color='#be2e2d',ls='--',lw=.8);ax.text(roots[level]*100,level,f' {roots[level]*100:.2f}%',fontsize=9)
        from matplotlib.colors import ListedColormap,BoundaryNorm
        cmap=ListedColormap(['#dfeaf0','#fcbd5a','#be2e2d']);norm=BoundaryNorm([.5,1.5,2.5,3.5],3)
        mapaxes[col].imshow(np.ma.masked_equal(full,0),cmap=cmap,norm=norm,interpolation='nearest');mapaxes[col].set_title(treatment+' 全季TN关键源区候选');mapaxes[col].axis('off')
        for p in [path,masspath]:audit['inputs'].append({'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
        audit['checks'][treatment]={'valid_pixels_preserved':True,'grid_crs_readback_pass':True,'nonnegative_loads':True,'class_partition_pass':int(sum((classes==k).sum() for k in [1,2,3]))==n,'class_mass_closure_kg':float(sum(vals[classes==k].sum() for k in [1,2,3])*area-vals.sum()*area),'output_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}
    from matplotlib.patches import Patch
    mapfig.legend(handles=[Patch(color=c,label=l) for c,l in [('#be2e2d','关键源区候选'),('#fcbd5a','次关键源区'),('#dfeaf0','其他有效种植区域')]],loc='lower center',ncol=3)
    fig.savefig(records/'负荷面积曲线与斜率.png',dpi=220);mapfig.savefig(records/'关键源区分布预览.png',dpi=220)
    write_csv(records/'class_statistics.csv',stats);write_csv(records/'sensitivity.csv',sens)
    (records/'thresholds_and_fit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    audit['passed']=all(v['class_partition_pass'] and abs(v['class_mass_closure_kg'])<1e-5 for v in audit['checks'].values());assert audit['passed']
    (records/'verification.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:{a:b for a,b in v.items() if a!='fit'} for k,v in summary.items()},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
