#!/usr/bin/env python3
"""Noise-resolved covariance and measurement-Fisher calculation.

This is a classical linearized covariance calculation around a drive-locked
orbit. It is not QFI and cannot support a quantum-advantage claim.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from reconstruction_core import ModelParameters, find_periodic_orbit, integrate_one_period, jacobian_phase, monodromy

OMEGA=np.array([[0,1,0,0,0,0],[-1,0,0,0,0,0],[0,0,0,1,0,0],[0,0,-1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,-1,0]],float)

def periodic_covariance(x0,p,n_th1=0.1,n_th2=0.1,detector_variance=0.01):
    orbit=integrate_one_period(x0,p,dense=True)
    if not orbit.success or orbit.sol is None: raise RuntimeError('orbit interpolation failed')
    n=6
    q=np.array([p.kappa/2,p.kappa/2,p.gamma1*(2*n_th1+1)/2,p.gamma1*(2*n_th1+1)/2,p.gamma2*(2*n_th2+1)/2,p.gamma2*(2*n_th2+1)/2],float)
    Q=np.diag(q)
    def fun(t,v):
        x=orbit.sol(t); A=jacobian_phase(x,p.drive_frequency*t,p); V=v.reshape(n,n,order='F'); return (A@V+V@A.T+Q).ravel(order='F')
    zero=solve_ivp(fun,(0,p.period),np.zeros(n*n),method='DOP853',rtol=1e-9,atol=1e-11,max_step=p.period/200)
    if not zero.success: raise RuntimeError(zero.message)
    qT=zero.y[:,-1]
    M=np.asarray(monodromy(x0,p)['monodromy'],float)
    K=np.eye(n*n)-np.kron(M,M)
    v0=np.linalg.solve(K,qT); V0=v0.reshape(n,n,order='F'); V0=0.5*(V0+V0.T)
    full=solve_ivp(fun,(0,p.period),V0.ravel(order='F'),method='DOP853',rtol=1e-9,atol=1e-11,dense_output=True,max_step=p.period/200)
    if not full.success or full.sol is None: raise RuntimeError(full.message)
    times=np.linspace(0,p.period,1001); Vs=np.array([0.5*(full.sol(t).reshape(n,n,order='F')+full.sol(t).reshape(n,n,order='F').T) for t in times])
    vmin=float(np.min([np.min(np.linalg.eigvalsh(V+0.5j*OMEGA).real) for V in Vs]))
    measured=np.asarray(Vs[:,0,0],float)+detector_variance; mean_variance=float(np.mean(measured)); mean_signal=float(np.mean(orbit.sol(times)[0]))
    return {'V0':V0.tolist(),'mean_signal':mean_signal,'mean_variance':mean_variance,'min_quantum_physicality_eigenvalue':vmin,'detector_variance':detector_variance,'thermal_occupations':[n_th1,n_th2]}

def one_theta(theta,eps,base,**noise):
    def obs(th):
        p=ModelParameters(theta=float(th)); orbit=find_periodic_orbit(np.zeros(6),p,max_iter=600)
        if orbit.get('status')!='PASS': raise RuntimeError(orbit.get('reason','periodic orbit failed'))
        return periodic_covariance(np.asarray(orbit['x0']),p,**noise), orbit
    plus,op=obs(theta+eps); minus,om=obs(theta-eps); center,oc=obs(theta)
    dm=(plus['mean_signal']-minus['mean_signal'])/(2*eps); dv=(plus['mean_variance']-minus['mean_variance'])/(2*eps); var=center['mean_variance']; fc=float(dm*dm/var+0.5*dv*dv/(var*var))
    return {'theta':float(theta),'status':'PASS' if min(center['min_quantum_physicality_eigenvalue'],plus['min_quantum_physicality_eigenvalue'],minus['min_quantum_physicality_eigenvalue'])>=-1e-8 and np.isfinite(fc) else 'FAIL','center':center,'finite_difference_eps':eps,'d_mean_d_theta':float(dm),'d_variance_d_theta':float(dv),'classical_fisher_information':fc,'interpretation':'Measurement Fisher information only; not QFI and not a gain over a reference.'}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--thetas',nargs='+',type=float,default=[0.0,np.pi/2,np.pi]); ap.add_argument('--eps',type=float,default=1e-4); ap.add_argument('--n-th1',type=float,default=0.1); ap.add_argument('--n-th2',type=float,default=0.1); ap.add_argument('--detector-variance',type=float,default=0.01); a=ap.parse_args()
    records=[]
    for th in a.thetas:
        try: records.append(one_theta(th,a.eps,None,n_th1=a.n_th1,n_th2=a.n_th2,detector_variance=a.detector_variance))
        except Exception as exc: records.append({'theta':th,'status':'FAIL','reason':f'{type(exc).__name__}: {exc}'})
    out={'status':'PASS' if records and all(r['status']=='PASS' for r in records) else 'FAIL','kind':'classical_measurement_fisher','recorded_at_utc':datetime.now(timezone.utc).isoformat(),'python':sys.version,'platform':platform.platform(),'measurement':{'observable':'cavity amplitude quadrature X_a = Re(alpha)','record':'y(t)=X_a(t)+nu_det(t)','detector_noise_variance':a.detector_variance,'thermal_occupations':[a.n_th1,a.n_th2]},'records':records,'interpretation':'This output is an operational classical Fisher calculation around stable periodic orbits. It is not QFI and contains no matched-baseline gain claim.'}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8'); print(json.dumps({'status':out['status'],'records':len(records),'output':str(a.output)},indent=2)); return 0 if out['status']=='PASS' else 1
if __name__=='__main__': raise SystemExit(main())
