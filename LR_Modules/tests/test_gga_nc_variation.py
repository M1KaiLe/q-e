#!/usr/bin/env python3
"""Compare actual dgradcorr with finite variations of a discrete GGA energy.

Only the FFT interface and a quadratic, spin-exchange invariant GGA are
substituted. The noncollinear projection and response use production Fortran.
The independent reference differentiates the nonlinear SCF potential on the
complete commensurate supercell, rather than constructing an analytic kernel.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def derivative(value, cell_count=1):
    wave=np.fft.fftfreq(len(value), d=cell_count/len(value))
    return np.fft.ifft(1j*wave[:,None]*np.fft.fft(value,axis=0),axis=0).real

def potential(rho, a=.7, c=.23, cell_count=1):
    """Gradient of E=integral[a/2 (|grad n+|²+|grad n-|²)+c grad n+·grad n-]."""
    mag=np.linalg.norm(rho[:,1:],axis=1)
    spins=np.stack(((rho[:,0]+mag)/2,(rho[:,0]-mag)/2),axis=1)
    lap=derivative(derivative(spins,cell_count),cell_count)
    vsp=-a*lap-c*lap[:,::-1]
    return np.column_stack((vsp.sum(axis=1)/2,
        (vsp[:,0]-vsp[:,1])[:,None]*rho[:,1:]/(2*mag[:,None]))),vsp

def serial(value):
    value=np.asarray(value).ravel(order='F')
    if np.iscomplexobj(value):
        return ' '.join(f'({x.real:.17e},{x.imag:.17e})' for x in value)
    return ' '.join(f'{x:.17e}' for x in value)

def sample(n, textured=True):
    x=2*np.pi*np.arange(n)/n
    rho=np.column_stack((2+.1*np.cos(x),.29*np.sin(x),
                         .21*np.cos(2*x),.63+.13*np.sin(3*x)))
    if not textured:
        rho[:,1:3]=0
    rng=np.random.default_rng(197)
    change=np.zeros((n,4),complex)
    for mode in range(-2,3):
        change+=np.exp(1j*mode*x)[:,None]*(rng.normal(size=4)+1j*rng.normal(size=4))[None,:]/7
    return rho,change

def finite_variation(rho,change,q,h):
    if q==0:
        result=np.zeros_like(change)
        for direction,factor in ((change.real,1),(change.imag,1j)):
            result+=factor*(potential(rho+h*direction)[0]-potential(rho-h*direction)[0])/(2*h)
        return result
    copies=4
    assert q==.25
    n=len(rho)
    x=2*np.pi*np.arange(n*copies)/n
    base=np.tile(rho,(copies,1))
    direction=(np.tile(change,(copies,1))*np.exp(1j*q*x)[:,None]).real
    diff=(potential(base+h*direction,cell_count=copies)[0]
          -potential(base-h*direction,cell_count=copies)[0])/(2*h)
    # Recover the complex +q amplitude of the real cosine perturbation.
    return 2*np.mean((diff*np.exp(-1j*q*x)[:,None]).reshape(copies,n,4),axis=0)

def build(source,directory):
    exe=directory/'probe.x'
    command=[shutil.which('gfortran'),'-O0','-g','-fcheck=all','-ffree-line-length-none',
             str(HERE/'gga_nc_variation_stubs.f90'),str(source),'-o',str(exe)]
    subprocess.run(command,cwd=directory,check=True,capture_output=True,text=True)
    return exe

def run(exe,rho,change,q,background_gradient=False):
    n=len(rho)
    vsp=potential(rho)[1]
    gmag=np.zeros((3,n,4))
    if background_gradient:
        gmag[0]=derivative(rho)
    lines=[f'{n} {q} 0 0 .7 .23',serial(rho),serial(change),serial(gmag),
           serial(vsp[:,0]-vsp[:,1]),serial(np.ones(n))]
    run=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
    numbers=np.fromstring(run.stdout,sep=' ').reshape(-1,2)
    return (numbers[:,0]+1j*numbers[:,1]).reshape((n,4),order='F')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'LR_Modules/dgradcorr.f90')
    parser.add_argument('--report',type=Path)
    parser.add_argument('--require-pass',action='store_true')
    args=parser.parse_args()
    records=[]
    with tempfile.TemporaryDirectory(prefix='qe-gga-nc-variation-') as tmp:
        exe=build(args.source,Path(tmp))
        for n,textured,q in ((25,False,0),(65,True,0),(65,True,.25),(17,True,0)):
            rho,change=sample(n,textured)
            fd1=finite_variation(rho,change,q,2e-5)
            fd2=finite_variation(rho,change,q,1e-5)
            fd=(4*fd2-fd1)/3
            for background_gradient in (False,True):
                response=run(exe,rho,change,q,background_gradient)
                relative=float(np.linalg.norm(response-fd)/np.linalg.norm(fd))
                records.append(dict(n=n,textured=textured,q=q,
                    background_gradient_supplied=background_gradient,
                    finite_difference_halving_change=float(np.linalg.norm(fd2-fd1)/np.linalg.norm(fd)),
                    relative_error=relative,max_abs_error=float(np.max(np.abs(response-fd))),
                    pass_=bool(relative<1e-7)))
    report=dict(source=str(args.source.resolve()),
                source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                functional='a/2*(|grad n+|²+|grad n-|²)+c*grad n+·grad n-',a=.7,c=.23,
                reference='Nonlinear discrete SCF potential central differences h and h/2; Richardson array extrapolation',
                records=records,all_pass=all(item['pass_'] for item in records))
    if args.report:
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if args.require_pass and not report['all_pass']:
        raise SystemExit(1)

if __name__=='__main__':
    main()
