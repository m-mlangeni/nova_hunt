import sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'..'); sys.path.insert(0,'.'); import nova_hunt as nh
nh.CFG.update(dmax_mpc=5, det_frac_min=0.25)
gal=pd.DataFrame(dict(ra=[10.,60.,120.],dec=[-30.,-40.,-50.],D=[3.,4.,2.],R1=[900.]*3,R2=[900.]*3,PA=[0.]*3,logM=[10.]*3))
gal.to_csv('g_blank.csv',index=False)   # NO name column, like REGALADE v2 thin
gl=nh.load_galaxies('g_blank.csv'); print('names:',gl.name.tolist())
rows=[]; did=[0]; iid=[5000]
def pair(gidx,dx,night):
    g=gl.iloc[gidx]; ra=g.ra+dx/3600/np.cos(np.radians(g.dec))
    for k,t in enumerate([night+0.80,night+0.85]):
        did[0]+=1; iid[0]+=1
        rows.append(dict(gid=g.gid,DETECTION_ID=did[0],IMAGE_ID=iid[0],mjd=t,filt='q',field_id=1,telescop='BG3',ra=ra,dec=g.dec,fnu=1.,
          mag=g.dm-7.8,magerr=.05,snr=10,rb=.9,sep_near_gaia=9.,mag_g_near_gaia=19.,sep_near_ref=9.,mag_auto_near_ref=23.,sep_bright_gaia=99.,mag_g_bright_gaia=20.))
# same night 60100: 2 pairs in each of 3 galaxies (6 total, but only 2 per galaxy) -> all kept
for gi in range(3):
    for dx in (200,-200): pair(gi,dx,60100)
# night 60200: 3 pairs in galaxy 0 -> a real bad-night group -> all removed
for dx in (300,-300,500): pair(0,dx,60200)
det=pd.DataFrame(rows)
def fake(sql,project,dry_run=False):
    if 'T_LMAG' not in sql: return det
    import re
    out=[]
    for s,ra,dec in re.findall(r'STRUCT\((\d+) AS src, ([-\d.]+) AS ra, ([-\d.]+) AS dec\)',sql):
        for _,x in det[(abs(det.ra-float(ra))<1e-3)&(abs(det.dec-float(dec))<1e-3)].iterrows():
            out.append(dict(src=int(s),image_id=x.IMAGE_ID,mjd=x.mjd,filt='q',t_lmag=21.5,tqc='green',ra_c=float(ra)+.1,dec_c=float(dec)))
    return pd.DataFrame(out)
nh.run_query=fake
F,_=nh.run('g_blank.csv',tns=None)
n60100=(np.floor(F.first_mjd-0.5)==60099).sum()+(np.floor(F.first_mjd-0.5)==60100).sum()
kept_100=((F.first_mjd>60100)&(F.first_mjd<60101)).sum(); kept_200=((F.first_mjd>60200)&(F.first_mjd<60201)).sum()
print(f"{'PASS' if kept_100==6 else 'FAIL'}: 6 single-night sources spread over 3 unnamed galaxies on one night -> kept {kept_100}")
print(f"{'PASS' if kept_200==0 else 'FAIL'}: 3 single-night sources in ONE galaxy on one night (bad epoch) -> kept {kept_200}")
print('host column filled:', F.host.notna().all(), F.host.iloc[0])
