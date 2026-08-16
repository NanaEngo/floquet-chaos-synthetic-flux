#!/usr/bin/env python3
"""Generate the validated pilot covariance/Fisher diagnostic figure."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--results',type=Path,required=True); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args()
    d=json.loads(a.results.read_text())
    if d.get('status')!='PASS' or not d.get('records'): raise SystemExit('ERROR: covariance/Fisher result must be PASS')
    if any(r.get('status')!='PASS' for r in d['records']): raise SystemExit('ERROR: all plotted records must be PASS')
    theta=np.array([r['theta'] for r in d['records']]); fisher=np.array([r['classical_fisher_information'] for r in d['records']]); phys=np.array([r['center']['min_quantum_physicality_eigenvalue'] for r in d['records']])
    a.output_dir.mkdir(parents=True,exist_ok=True); fig,ax=plt.subplots(1,2,figsize=(7,3),constrained_layout=True)
    ax[0].plot(theta,fisher,'o-',color='#1f77b4'); ax[0].set_yscale('symlog',linthresh=1e-25); ax[0].set_xlabel(r'$\theta$ (rad)'); ax[0].set_ylabel(r'$F_C$'); ax[0].set_title('Measured-record Fisher information')
    ax[1].plot(theta,phys,'o-',color='#2ca02c'); ax[1].set_xlabel(r'$\theta$ (rad)'); ax[1].set_ylabel('Minimum physicality eigenvalue'); ax[1].set_title('Covariance physicality check')
    stem=a.output_dir/'covariance_fisher_pilot'; fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(stem.with_suffix('.png'),dpi=300,bbox_inches='tight'); plt.close(fig)
    manifest={'status':'PASS','source':str(a.results.resolve()),'outputs':[str(stem.with_suffix('.pdf')),str(stem.with_suffix('.png'))],'no_mock_data':True,'interpretation':'Pilot classical Fisher/physicality diagnostic; not QFI and not a gain comparison.'}
    (a.output_dir/'covariance_fisher_pilot.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__': main()
