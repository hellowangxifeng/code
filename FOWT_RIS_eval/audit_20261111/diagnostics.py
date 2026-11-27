"""Independent read-only diagnostics; writes evidence only into this audit folder."""
from pathlib import Path
import sys, json, ast, hashlib
import numpy as np
from scipy import signal
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
import fx_optics as fo
import replay as rp
import control_study as cs
import fowt_data as fd

evidence = {}
d = fo.load_design()
ref = np.asarray(d['referenceXY'], float)
power = np.asarray(d['unitDirectPowerW'], float).copy()
power /= power.sum()
p0 = np.asarray(d['pRis'], float)
qin0 = fo.q_in(d, np.eye(3), p0, ref)
states = rp.design_setup(d, power, ref)
evidence['design'] = {k: np.asarray(d[k]).tolist() for k in ['pTx','pRis','RRis','rxPlaneZ','nWater','nGlass','nAir','B','maxCmdDeg','angleStepDeg','currentCenter','improvedCenter','improvedU','improvedQ']}
evidence['design']['qstar_recomputed'] = states['collimated']['qs'].tolist()
evidence['design']['qin_weighted'] = (power @ qin0).tolist()
evidence['raw_path'] = {'configured': fd.RAW_DIR, 'exists': Path(fd.RAW_DIR).exists(), 'actual_exists': (ROOT/'浮式平台六自由度实验数据'/'FOWT Platform Testing').exists(), 'catalog_count':len(fd.catalog())}
evidence['aperture'] = {}
for tag in ['current','improved']:
    x,y,I,c = [np.asarray(d[tag+k],float) for k in ['X','Y','I','Center']]
    X,Y=np.meshgrid(x,y)
    area=(x[1]-x[0])*(y[1]-y[0])
    ap=rp.aperture_map(d,tag,c,0.05)
    reported=rp.aperture_at(ap,0,0)
    direct=I[(X-c[0])**2+(Y-c[1])**2<=.05**2].sum()*area
    recentered=float(__import__('scipy').ndimage.map_coordinates(ap['corr'], [[(len(y)-1)//2],[(len(x)-1)//2]],order=1)[0])
    evidence['aperture'][tag]={'code_mw':reported*1000,'direct_mw':direct*1000,'corr_at_array_center_mw':recentered*1000,'map_midpoint':[(x[0]+x[-1])/2,(y[0]+y[-1])/2],'aperture_center':c.tolist(),'dxdy':[x[1]-x[0],y[1]-y[0]]}
# Off-center impulse gives an unambiguous coordinate counterexample.
x=np.arange(11,dtype=float); I=np.zeros((11,11)); I[3,8]=1
toy={'toyX':x,'toyY':x,'toyI':I}
ap=rp.aperture_map(toy,'toy',np.array([8.,3.]),.1)
evidence['aperture']['off_center_impulse']={'expected':1.,'code':rp.aperture_at(ap,0,0)}

Jold,_=cs.sensitivity(d,power,ref,states['point'])
_,qo,_=rp.point_targets(d,np.eye(3),p0,ref,qin0,np.r_[states['point']['receiver'],d['rxPlaneZ']])
bfrozen=qo-qin0
def ctr(rpy,dp):
    R=fo.euler_matrix(rpy); p=p0+dp
    qi=fo.q_in(d,R,p,ref)
    land,_=fo.unit_landing(d,R,p,ref,qi+bfrozen)
    return power @ land[:,:2]
J=np.zeros((2,6)); c0=ctr(np.zeros(3),np.zeros(3))
for j in range(6):
    eps=1e-4 if j<3 else 1e-5
    r=np.zeros(3); dp=np.zeros(3)
    (r if j<3 else dp)[j%3]=np.rad2deg(eps) if j<3 else eps
    J[:,j]=(ctr(r,dp)-ctr(-r,-dp))/(2*eps)
evidence['sensitivity']={'code':Jold.tolist(),'frozen_command_with_recomputed_incidence':J.tolist(),'norm_relative_error':float(np.linalg.norm(Jold-J)/np.linalg.norm(J))}
evidence['static_collimated_grid']={'max_gradient_vs_snapped_error':float(np.max(abs(states['collimated']['b']-d['nAir']*np.sin(np.deg2rad(fo.snap_commands(d,qin0,np.tile(states['collimated']['qs'],(len(ref),1)))[0])))))}
evidence['point_nominal']={}
for key in ['current_point','point']:
    receiver=np.r_[states[key]['receiver'],d['rxPlaneZ']]
    cmd,qo,valid=rp.point_targets(d,np.eye(3),p0,ref,qin0,receiver)
    land,dw=fo.unit_landing(d,np.eye(3),p0,ref,qo)
    delta=power @ land[:,:2]-receiver[:2]
    evidence['point_nominal'][key]={'centroid_error_mm':(delta*1000).tolist(),'radial_centroid_error_mm':float(np.linalg.norm(delta)*1000),'invalid_power':float(power[~valid].sum()),'rms_reference_rays_mm':float(np.sqrt(power @ np.sum((land[:,:2]-receiver[:2])**2,axis=1))*1000)}
evidence['units']={}
# mrad stored as degrees * 1000 instead of radians * 1000.
evidence['units']['stored_angle_over_true_mrad']=180/np.pi
# Causal predictor index test: current implementation reads pose[-1] on first quadratic prediction.
pose=np.arange(20,dtype=float)[:,None]; pred1=cs.predict(pose,2,2)
pose2=pose.copy(); pose2[-1]+=1000; pred2=cs.predict(pose2,2,2)
evidence['predictor_future_leak']={'first_affected_index':int(np.flatnonzero(pred1[:,0]!=pred2[:,0])[0]),'original':float(pred1[3,0]),'changed_future_last_sample':float(pred2[3,0])}
cmd,qo,valid=fo.snap_commands(d,np.array([[1.,.5]]),np.array([[.99,0.]]))
legal_other=np.array([[1.,.5]])+d['nAir']*np.sin(np.deg2rad([[-d['maxCmdDeg'],-d['maxCmdDeg']]]))
evidence['local_repair_counterexample']={'returned_valid':bool(valid[0]),'returned_q':qo.tolist(),'legal_alternative_q':legal_other.tolist(),'alternative_norm':float(np.linalg.norm(legal_other))}
qs,cov=fo.power_weighted_target(d,np.array([[1.26,.42]]),np.array([1.]))
evidence['target_disk_counterexample']={'returned_target':qs.tolist(),'norm':float(np.linalg.norm(qs)),'claimed_coverage':cov}
ctrlrows=list(__import__('csv').DictReader(open(ROOT/'results/control_study.csv',encoding='utf-8')))
rec=rp.load_motion()[8]; tt,poses=cs.pose_series(rec)
ctrl=[]
for delay in [2,5,10]:
    est=cs.predict(poses,delay,1); err=poses-est
    ctrl.append({'delay_ms':delay*10,'code_J_rms_mm':float(np.sqrt(np.mean(np.sum((err@Jold.T)**2,axis=1)))*1000),'recomputed_incidence_J_rms_mm':float(np.sqrt(np.mean(np.sum((err@J.T)**2,axis=1)))*1000)})
evidence['latency_linearised_comparison']=ctrl
evidence['python_compile']={}
inventory=[]
for f in sorted(ROOT.rglob('*')):
    if f.suffix not in ['.m','.py'] or 'audit_20261111' in str(f): continue
    b=f.read_bytes(); src=b.decode('utf-8-sig'); rel=str(f.relative_to(ROOT))
    inventory.append({'file':rel,'lines':len(src.splitlines()),'sha256':hashlib.sha256(b).hexdigest()})
    if f.suffix=='.py':
        try: ast.parse(src); evidence['python_compile'][rel]='PASS'
        except SyntaxError as ex: evidence['python_compile'][rel]=str(ex)
evidence['inventory']=inventory
evidence['counts']={'m':sum(x['file'].endswith('.m') for x in inventory),'py':sum(x['file'].endswith('.py') for x in inventory),'total_lines':sum(x['lines'] for x in inventory)}
# Inspect saved atlas without executing its authoring script.
atlas=sio.loadmat(ROOT/'figure_atlas_matlab/data/atlas_data.mat',squeeze_me=True)
evidence['atlas']={'focus_gain':float(atlas['focus_gain']),'focus_P5_mw':float(atlas['focus_P5']),'wind_U':atlas['wind_U'].tolist(),'wind_tilt':atlas['wind_tilt'].tolist(),'ap_focus_over_linear':(atlas['ap_P_focus']/atlas['ap_P_linear']).tolist(),'pcondition_factor':.7,'pcondition_dbm':float(30+10*np.log10(.7))}
out=Path(__file__).with_name('python_evidence.json')
out.write_text(json.dumps(evidence,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in evidence.items() if k not in ['inventory','python_compile','atlas']},ensure_ascii=False,indent=2))
print('COUNTS',evidence['counts'],'EVIDENCE',out)
