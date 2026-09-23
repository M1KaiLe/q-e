#!/usr/bin/env python3
"""Rebuild static dynamical matrices from COMPLETE saved PH perturbations.

Requires numpy and the corrected q2qstar.x. Reads source files only. A new output
folder is required. Existing exported matrices supply geometry/q labels only;
all force-constant entries come from PARTIAL_DYN and DISPLACEMENT_PATTERN.
"""
from pathlib import Path
import argparse, copy, hashlib, json, os, shutil, subprocess
import xml.etree.ElementTree as ET
import numpy as np

def complex_array(text):
 v=np.fromstring(text,sep=' ').reshape(-1,2);return v[:,0]+1j*v[:,1]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--dyn-dir',type=Path,required=True);p.add_argument('--phsave',type=Path,required=True)
 p.add_argument('--prefix',required=True);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--q2qstar',type=Path,default=Path(__file__).resolve().parents[1]/'q2qstar.x')
 p.add_argument('--unitary-magnetic',action='store_true')
 a=p.parse_args();a.output=a.output.resolve();a.q2qstar=a.q2qstar.resolve()
 a.output.mkdir(parents=True,exist_ok=False);seeds=a.output/'seeds';seeds.mkdir()
 dyn0=a.dyn_dir/(a.prefix+'.dyn0');shutil.copy2(dyn0,a.output/dyn0.name)
 lines=dyn0.read_text().splitlines();nq=int(lines[1].split()[0]);manifest={'q2qstar_sha256':digest(a.q2qstar),'unitary_magnetic':a.unitary_magnetic,'sources':{},'outputs':{}}
 env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',I_MPI_FABRICS='shm')
 for iq in range(1,nq+1):
  path=a.dyn_dir/f'{a.prefix}.dyn{iq}.xml';tree=ET.parse(path);root=tree.getroot()
  nat=int(root.findtext('GEOMETRY_INFO/NUMBER_OF_ATOMS'));nm=3*nat
  pat=a.phsave/f'patterns.{iq}.xml';pr=ET.parse(pat).getroot()
  U=np.array([complex_array(n.text) for n in pr.iter('DISPLACEMENT_PATTERN')]).T
  if U.shape!=(nm,nm):raise ValueError(f'incomplete pattern basis: {pat}')
  if np.max(abs(U.conj().T@U-np.eye(nm)))>1e-8:raise ValueError(f'non-orthonormal patterns: {pat}')
  nirr=int(pr.findtext('.//NUMBER_IRR_REP'));dyn=np.zeros((nm,nm),complex)
  sources=[path,pat]
  for irr in range(nirr+1):
   f=a.phsave/f'dynmat.{iq}.{irr}.xml';r=ET.parse(f).getroot()
   if r.findtext('.//DONE_IRR','').strip().lower()!='true':raise ValueError(f'incomplete perturbation: {f}')
   dyn+=complex_array(r.findtext('.//PARTIAL_DYN')).reshape(nm,nm,order='F');sources.append(f)
  dyn=U@dyn@U.conj().T
  if not np.all(np.isfinite(dyn)):raise ValueError(f'non-finite response at iq={iq}')
  # q2qstar performs QE's normal little-group symmetrization and writes fresh
  # frequencies/eigenvectors. Never reuse the contaminated exported matrix.
  block=root.find('DYNAMICAL_MAT_.1')
  for i in range(nat):
   for j in range(nat):
    z=dyn[3*i:3*i+3,3*j:3*j+3].ravel(order='F')
    block.find(f'PHI.{i+1}.{j+1}').text='\n'+''.join(f'{v.real:24.16e} {v.imag:24.16e}\n' for v in z)
  root.find('GEOMETRY_INFO/NUMBER_OF_Q').text='1'
  for node in list(root):
   if (node.tag.startswith('DYNAMICAL_MAT_') and node is not block) or node.tag=='FREQUENCIES_THZ_CMM1':root.remove(node)
  seed=seeds/path.name;ET.indent(tree);tree.write(seed,encoding='utf-8',xml_declaration=True)
  output=a.output/f'{a.prefix}.dyn{iq}'
  cmd=[str(a.q2qstar),str(seed),str(output)]
  if a.unitary_magnetic:cmd.append('--unitary-magnetic')
  with (a.output/f'q2qstar-{iq}.log').open('w') as log:
   subprocess.run(cmd,cwd=a.output,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
  for f in sources:manifest['sources'][str(f.resolve())]=digest(f)
  target=Path(str(output)+'.xml');manifest['outputs'][target.name]=digest(target)
  print('Regenerated',target.name,flush=True)
 (a.output/'regeneration.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
