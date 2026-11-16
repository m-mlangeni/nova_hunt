"""Host naming: names file, NED lookup (faked), saved lookups, and --no-name-lookup."""
import sys, re, os, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'..'); sys.path.insert(0,'.'); import nova_hunt as nh
from astropy.table import Table
import astroquery.ipac.ned as nedmod
nh.CFG.update(dmax_mpc=5, det_frac_min=0.25)
pd.DataFrame(dict(ra=[10.,60.,120.],dec=[-30.,-40.,-50.],D=[3.,4.,2.],R1=[900.]*3,R2=[900.]*3,PA=[0.]*3,logM=[10.]*3)).to_csv('g_names.csv',index=False)
gl=nh.load_galaxies('g_names.csv')
rows=[]; did=[0]; iid=[5000]
for gi in range(3):                      # one clean multi-night nova-like source per galaxy
    g=gl.iloc[gi]; ra=g.ra+200/3600/np.cos(np.radians(g.dec))
    for d in [0,2,5,9]:
        did[0]+=1; iid[0]+=1
        rows.append(dict(gid=g.gid,DETECTION_ID=did[0],IMAGE_ID=iid[0],mjd=60100.8+d,filt='q',field_id=1,telescop='BG3',ra=ra,dec=g.dec,fnu=1.,
          mag=g.dm-8+0.15*d,magerr=.05,snr=10,rb=.9,sep_near_gaia=9.,mag_g_near_gaia=19.,sep_near_ref=9.,mag_auto_near_ref=23.,sep_bright_gaia=99.,mag_g_bright_gaia=20.))
det=pd.DataFrame(rows)
def fake(sql,project,dry_run=False):
    if 'T_LMAG' not in sql: return det
    out=[]
    for s,ra,dec in re.findall(r'STRUCT\((\d+) AS src, ([-\d.]+) AS ra, ([-\d.]+) AS dec\)',sql):
        for _,x in det[(abs(det.ra-float(ra))<1e-3)&(abs(det.dec-float(dec))<1e-3)].iterrows():
            out.append(dict(src=int(s),image_id=x.IMAGE_ID,mjd=x.mjd,filt='q',t_lmag=21.5,tqc='green',ra_c=float(ra)+.1,dec_c=float(dec)))
    return pd.DataFrame(out)
nh.run_query=fake
calls={'n':0}
class FakeNed:
    @staticmethod
    def query_region(c, radius):
        calls['n']+=1
        nm={60.0:'NGC 1234',120.0:'ESO 999-G099'}[round(c.ra.deg,1)]
        return Table({'Object Name':[nm,'SDSS star'],'Type':['G','*'],'Separation':[0.1,0.05]})
nedmod.Ned=FakeNed
def seed(): pd.DataFrame([dict(ra=10.,dec=-30.,name='PGC2789',common='',source='seed')]).to_csv('galaxy_names.csv',index=False)
seed()
F,_=nh.run('g_names.csv',tns=None)
want={'PGC2789 = NGC 253','J060.0000-40.0000 = NGC 1234','J120.0000-50.0000 = ESO 999-G099'}
print('PASS' if set(F.host)==want else 'FAIL','names from file + lookup:',sorted(F.host),'| lookups:',calls['n'])
n0=calls['n']; nh.run('g_names.csv',tns=None)
print('PASS' if calls['n']==n0 else 'FAIL','second run reuses galaxy_names.csv (no new lookups)')
nh.CFG['name_lookup']=False; seed()
F,_=nh.run('g_names.csv',tns=None)
print('PASS' if calls['n']==n0 and 'PGC2789 = NGC 253' in set(F.host) else 'FAIL','--no-name-lookup: names file still used, no lookups')

# --- retries, group-vs-galaxy preference, failures never saved, old "not found" rows retried
import time; time.sleep = lambda s: None
nh.CFG['name_lookup']=True
state={'fail':2}
class FlakyNed:
    @staticmethod
    def query_region(c, radius):
        calls['n']+=1
        if state['fail']>0:
            state['fail']-=1; raise ConnectionError('rate limited')
        if round(c.ra.deg,1)==60.0:   # group listed closer than the galaxy itself
            return Table({'Object Name':['NGC 1234 GROUP','NGC 1234'],'Type':['GGroup','G'],'Separation':[0.01,0.2]})
        return Table({'Object Name':['ESO 999-G099'],'Type':['G'],'Separation':[0.1]})
nedmod.Ned=FlakyNed
import astroquery.simbad as sm
class DeadSimbad:
    @staticmethod
    def query_region(c, radius): raise RuntimeError('CDS down')
sm.Simbad=DeadSimbad
# names file as left by the old version: two failed lookups saved as "not found"
pd.DataFrame([dict(ra=10.,dec=-30.,name='PGC2789',common='',source='seed'),
              dict(ra=60.,dec=-40.,name='J060.0000-40.0000',common='',source='not found'),
              dict(ra=120.,dec=-50.,name='J120.0000-50.0000',common='',source='not found')]).to_csv('galaxy_names.csv',index=False)
F,_=nh.run('g_names.csv',tns=None)
print('PASS' if set(F.host)==want else 'FAIL','old "not found" retried; NED fails twice then answers; galaxy preferred over group:',sorted(F.host))
g=pd.read_csv('galaxy_names.csv'); print('PASS' if not (g.source=='not found').any() and len(g)==3 else 'FAIL','names file cleaned, no duplicates:',len(g),'rows')
# total outage: nothing saved, retried next time
class DeadNed:
    @staticmethod
    def query_region(c, radius): raise ConnectionError('down')
nedmod.Ned=DeadNed
pd.DataFrame([dict(ra=10.,dec=-30.,name='PGC2789',common='',source='seed')]).to_csv('galaxy_names.csv',index=False)
F,_=nh.run('g_names.csv',tns=None)
g=pd.read_csv('galaxy_names.csv')
print('PASS' if len(g)==1 else 'FAIL','outage: no failed lookups saved (rows:',len(g),'); hosts fall back to coordinate names:',sorted(F.host))
