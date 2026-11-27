"""Independent deterministic audit using a conditional-Gaussian identity.

Unlike bsf_ballistic.m, the source Gaussian product is integrated analytically
before cubature. This checks finite waist, incidence Jacobian and cell averaging.
"""
from pathlib import Path
import json
import h5py
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite_e import hermegauss

ROOT=Path(__file__).resolve().parent


def read_scalar(group,name):
    return float(np.asarray(group[name]).squeeze())


def calculate(group):
    config=group["config"]
    geometry=group["geometry"]
    power=read_scalar(config["source"],"power")
    waist=read_scalar(config["source"],"waist")
    angle=read_scalar(config["source"],"divergence")
    extinction=read_scalar(config["water"],"absorption")+read_scalar(config["water"],"scattering")
    RT=np.asarray(geometry["RTx"]).T
    RR=np.asarray(geometry["RRis"]).T
    pT=np.asarray(geometry["pTx"]).ravel()
    pR=np.asarray(geometry["pRis"]).ravel()
    normal=np.asarray(geometry["normal"]).ravel()
    normalT=RT.T@normal
    x=np.asarray(group["x"]).ravel()
    y=np.asarray(group["y"]).ravel()
    X,Y=np.meshgrid(x,y)
    dx,dy=x[1]-x[0],y[1]-y[0]
    nodes,weights=leggauss(6)
    gaussian_nodes,gaussian_weights=hermegauss(5)
    gaussian_weights=gaussian_weights/np.sqrt(2*np.pi)
    result=np.zeros_like(X)
    for i in range(6):
        for j in range(6):
            q=np.stack([X+dx/2*nodes[i],Y+dy/2*nodes[j],np.zeros_like(X)],axis=-1)@RR.T+pR
            r=(q-pT)@RT
            z=r[...,2]
            variance=(waist**2+angle**2*z*z)/4
            gaussian_density=np.exp(-np.sum(r[...,:2]**2,axis=-1)/(2*variance))/(2*np.pi*variance)
            conditional_mean=(waist**2/4/variance)[...,None]*r[...,:2]
            conditional_sd=np.sqrt(waist**2*angle**2*z*z/(16*variance))
            expectation=np.zeros_like(X)
            for k in range(5):
                for l in range(5):
                    u=conditional_mean[...,0]+conditional_sd*gaussian_nodes[k]
                    v=conditional_mean[...,1]+conditional_sd*gaussian_nodes[l]
                    rr=r-np.stack([u,v,np.zeros_like(X)],axis=-1)
                    ray_length=np.linalg.norm(rr,axis=-1)
                    h=np.maximum(rr@normalT,0)/z
                    expectation+=gaussian_weights[k]*gaussian_weights[l]*h*np.exp(-extinction*ray_length)
            result+=weights[i]*weights[j]/4*power*gaussian_density*expectation
    return result


def main():
    records=[]
    for path in sorted((ROOT/"results").glob("panel_r*_c*.mat")):
        with h5py.File(path,"r") as data:
            group=data["result"]
            if "transportMetrics" not in group:
                continue
            actual=np.asarray(group["Ib"])
            independent=calculate(group)
            records.append(dict(file=path.name,relativeL2=float(np.linalg.norm(actual-independent)/np.linalg.norm(independent)),
                maxAbsoluteOverPeak=float(np.max(np.abs(actual-independent))/independent.max())))
    audit=dict(method="conditional Gaussian marginalization + 5x5 normal cubature + 6x6 area quadrature",
        panels=records,maximumRelativeL2=max(r["relativeL2"] for r in records),
        pass_=all(r["relativeL2"]<1e-10 for r in records))
    audit["passed"]=audit.pop("pass_")
    (ROOT/"results"/"ballistic_independent_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print(json.dumps(audit,indent=2))


if __name__=="__main__":
    main()
