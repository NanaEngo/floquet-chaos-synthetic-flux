#!/usr/bin/env python3
"""Finalize a reconstruction run without changing any scientific values."""
from __future__ import annotations
import argparse, hashlib, json, platform, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def git_revision(path: Path):
    try: return subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip()
    except (OSError, subprocess.CalledProcessError): return None

def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument('--results',type=Path,required=True); p.add_argument('--workspace',type=Path,required=True); a=p.parse_args()
    results=a.results.resolve(); workspace=a.workspace.resolve()
    model=json.loads((results/'model_gate.json').read_text())
    floquet=json.loads((results/'floquet_lyapunov.json').read_text())
    env=json.loads((results/'environment_manifest.json').read_text())
    status='PASS' if model.get('status')=='PASS' and floquet.get('status')=='PASS' and env.get('status')=='PASS' else 'FAIL'
    params=floquet.get('parameters', model.get('parameters', {}))
    (results/'parameters.json').write_text(json.dumps(params,indent=2,sort_keys=True)+'\n')
    extended = {}
    if (results/'flux_grid.json').is_file():
        extended['flux_grid'] = json.loads((results/'flux_grid.json').read_text()).get('status')
    if (results/'covariance_fisher.json').is_file():
        extended['covariance_fisher'] = json.loads((results/'covariance_fisher.json').read_text()).get('status')
    extended_status = 'PASS' if extended and all(value == 'PASS' for value in extended.values()) else ('PARTIAL' if extended else 'NOT_RUN')
    status_record={'status':status,'model_gate':model.get('status'),'floquet_lyapunov_gate':floquet.get('status'),'environment_gate':env.get('status'),'extended_status':extended_status,'extended_gates':extended,'reason':floquet.get('reason')}
    (results/'status.json').write_text(json.dumps(status_record,indent=2,sort_keys=True)+'\n')
    input_results = ['model_gate.json', 'floquet_lyapunov.json']
    for candidate in ('flux_grid.json', 'covariance_fisher.json'):
        if (results / candidate).is_file():
            input_results.append(candidate)
    manifest={
      'status':status,
      'created_at_utc':datetime.now(timezone.utc).isoformat(),
      'python':sys.version,
      'platform':platform.platform(),
      'workspace':str(workspace),
      'workspace_revision':git_revision(workspace),
      'fl_qom_root':env.get('fl_qom_root'),
      'fl_qom_revision':env.get('fl_qom_revision'),
      'environment_manifest':'environment_manifest.json',
      'source_hashes':'source_hashes.json',
      'input_results':input_results,
      'parameters':'parameters.json',
      'status_file':'status.json',
      'no_mock_data':True,
      'extended_status':extended_status,
      'extended_gates':extended,
    }
    (results/'run_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    checksum_paths=sorted(p for p in results.iterdir() if p.is_file() and p.name!='checksums.sha256')
    lines=[f'{sha256(p)}  {p.name}' for p in checksum_paths]
    (results/'checksums.sha256').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'status':status,'files_hashed':len(lines),'results':str(results)},indent=2))
    return 0 if status=='PASS' else 1
if __name__=='__main__': raise SystemExit(main())
