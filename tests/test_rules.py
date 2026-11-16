import re, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'..'); sys.path.insert(0,'.'); import nova_hunt as nh
nh.CFG.update(dmax_mpc=5, det_frac_min=0.25)
gal=pd.DataFrame(dict(name=['GalA','GalB'],ra=[10.,60.],dec=[-30.,-40.],D=[3.,4.],R1=[1500.,900.],R2=[1500.,900.],PA=[0.,0.],logM=[10.5,10.]))
gal.to_csv('g.csv',index=False)
mu={'GalA':5*np.log10(3e6)-5,'GalB':5*np.log10(4e6)-5}
rows=[]; did=[0]; img_of={}
CASES={}
def src(name,expect,gal_,dx,dy,epochs,sepg=9.,gg=19.,sepr=9.,refm=23.):
    """epochs: list of (mjd, filt, mag, rb, tel, image_id)"""
    g=gal.set_index('name').loc[gal_]
    ra=g.ra+dx/3600/np.cos(np.radians(g.dec)); dec=g.dec+dy/3600
    for (t,f,m,rb,tel,iid) in epochs:
        did[0]+=1
        rows.append(dict(gid=None,gname=gal_,DETECTION_ID=did[0],IMAGE_ID=iid,mjd=t,filt=f,field_id=1,telescop=tel,ra=ra,dec=dec,fnu=1.,
            mag=m,magerr=.05,snr=10,rb=rb,sep_near_gaia=sepg,mag_g_near_gaia=gg,sep_near_ref=sepr,mag_auto_near_ref=refm,sep_bright_gaia=99.,mag_g_bright_gaia=20.))
    CASES[name]=(expect,ra,dec)
iid=[1000]
def ep(t,f,m,rb=.9,tel='BG3'):
    iid[0]+=1; return (t,f,m,rb,tel,iid[0])
A='GalA'; pk=lambda M: M+mu[A]
# 1 classic nova, 4 nights fading
src('nova_multi',True,A,300,0,[ep(60000.8+d,'q',pk(-8)+0.15*d) for d in [0,2,5,9]])
src('too_bright_SN',False,A,-300,0,[ep(60000.8+d,'q',pk(-13)+0.05*d) for d in [0,3,6,10]])
src('on_nucleus',False,A,0.5,0.5,[ep(60000.8+d,'q',pk(-8)+0.15*d) for d in [0,2,5,9]])
src('single_det',False,A,0,300,[ep(60000.8,'q',pk(-8))])
src('low_mean_rb',False,A,0,-300,[ep(60000.8+d,'q',pk(-8)+0.15*d,rb=.3) for d in [0,2,5,9]])
src('single_4det_2tel',True,A,400,400,[ep(60010.80,'q',pk(-8),tel='BG2'),ep(60010.81,'q',pk(-8),tel='BG3'),ep(60010.85,'q',pk(-8),tel='BG2'),ep(60010.86,'q',pk(-8),tel='BG3')])
src('single_4det_1tel',False,A,-400,400,[ep(60011.80+0.003*k,'q',pk(-8)) for k in range(4)])
src('pair_1p25h',True,A,400,-400,[ep(60012.80,'q',pk(-7.8),rb=.87),ep(60012.852,'q',pk(-7.8),rb=.9)])
src('pair_5min',False,A,-400,-400,[ep(60013.80,'q',pk(-7.8)),ep(60013.8035,'q',pk(-7.8))])
src('pair_lowrb',False,A,600,0,[ep(60014.80,'q',pk(-7.8),rb=.6),ep(60014.852,'q',pk(-7.8),rb=.9)])
src('long_300d',False,A,-600,0,[ep(60000.8+d,'q',pk(-8)+0.002*d) for d in [0,100,200,300]])
src('no_fade_60d',False,A,0,600,[ep(60000.8+d,'q',pk(-8)) for d in [0,20,40,60]])
src('gaia_coincident',False,A,0,-600,[ep(60000.8+d,'q',pk(-8)+0.15*d) for d in [0,2,5,9]],sepg=0.3,gg=21.5)
src('gaia_faint_0p9',True,A,700,100,[ep(60000.8+d,'q',pk(-8)+0.15*d) for d in [0,2,5,9]],sepg=0.9,gg=pk(-8)+2.3)
src('gaia_bright_0p9',False,A,-700,100,[ep(60000.8+d,'q',pk(-8)+0.15*d) for d in [0,2,5,9]],sepg=0.9,gg=pk(-8)-0.5)
src('ref_counterpart',False,A,100,700,[ep(60000.8+d,'q',pk(-8)+0.15*d) for d in [0,2,5,9]],sepr=0.5,refm=pk(-8)+0.5)
for k,(dx,dy) in enumerate([(800,800),(-800,800),(800,-800)]):   # bad night: 3 pairs, same host+night
    src(f'badnight_{k}',False,A,dx,dy,[ep(60020.80,'q',pk(-7.8)),ep(60020.85,'q',pk(-7.8))])
src('stray_lowrb_peak',True,A,-100,-700,[ep(60000.8,'q',pk(-12.5),rb=.15)]+[ep(60001.8+d,'q',pk(-8)+0.15*d,rb=.95) for d in [0,2,5,9]])
src('flicker_1of8',False,'GalB',300,0,[ep(60030.8,'q',18.5),ep(60033.8,'q',18.6)])
det=pd.DataFrame(rows)
gl=nh.load_galaxies('g.csv'); det['gid']=det.gname.map(dict(zip(gl.name,gl.gid))); det=det.drop(columns='gname')
# images: every detection image + extra non-detection exposures for flicker_1of8 (7 extra each night)
def fake(sql,project,dry_run=False):
    if 'T_LMAG' not in sql: return det
    out=[]
    for s,ra,dec in re.findall(r'STRUCT\((\d+) AS src, ([-\d.]+) AS ra, ([-\d.]+) AS dec\)',sql):
        ra,dec=float(ra),float(dec)
        near=det[(abs(det.ra-ra)<1e-3)&(abs(det.dec-dec)<1e-3)]
        for _,x in near.iterrows():
            out.append(dict(src=int(s),image_id=x.IMAGE_ID,mjd=x.mjd,filt=x.filt,t_lmag=21.5,tqc='green',ra_c=ra+.1,dec_c=dec))
        if abs(ra-(60+300/3600/np.cos(np.radians(-40))))<1e-3:
            for n in [60030.8,60033.8]:
                for k in range(1,8): out.append(dict(src=int(s),image_id=900000+int(n)*10+k,mjd=n+0.01*k,filt='q',t_lmag=21.5,tqc='green',ra_c=ra+.1,dec_c=dec))
        out.append(dict(src=int(s),image_id=990000+int(s),mjd=59990.8,filt='q',t_lmag=22.5,tqc='green',ra_c=ra+.1,dec_c=dec))  # deep pre-image
    return pd.DataFrame(out)
nh.run_query=fake
F,d=nh.run('g.csv',tns=None)
got=set()
for _,r in F.iterrows():
    for n,(e,ra,dec) in CASES.items():
        if abs(r.ra-ra)<1e-4 and abs(r.dec-dec)<1e-4: got.add(n)
bad=0
for n,(e,_,_) in CASES.items():
    ok=(n in got)==e; bad+=not ok
    print(f"{'PASS' if ok else 'FAIL'}  {n:<20s} expected {'kept' if e else 'rejected'}, got {'kept' if n in got else 'rejected'}")
r=F[[abs(a-CASES['stray_lowrb_peak'][1])<1e-4 for a in F.ra]]
print('stray_lowrb_peak M_peak =',round(r.M_peak.iloc[0],2),'(should be -8, not -12.5)')
print('\nALL RULE TESTS PASSED' if bad==0 else f'\n{bad} TEST(S) FAILED')
