from pathlib import Path
import hashlib,json,subprocess
import yaml
root=Path(__file__).resolve().parent
records=json.loads((root/'provenance.json').read_text())
for rec in records:
 p=root/rec['destination']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==rec['prepared_sha256'], str(p)
cases=json.loads((root/'cases.json').read_text())
assert len(cases)==24
for case in cases:
 d=root/case['name']; cfg=yaml.safe_load((d/'config.yaml').read_text())
 assert cfg['mesh']['domain']==[-15.,15.,-5.,5.]
 subprocess.run(['bash','-n',str(d/'run.sh')],check=True)
 if case['kind']!='spatial':
  assert cfg['source']['derivative']=='adolc'
  assert cfg['sampling']['seed']==20260914
  assert cfg['solver']['diagonal_shift']==0
  assert cfg['solver']['max_iter']==5000
  assert cfg['spectral']['requested_eigenpairs']==8
  assert cfg['krylov']['relative_tolerance']==0.01
  assert cfg['krylov']['max_iterations']==50
  expected=2000 if case['kind'].startswith('official:') else (100 if case['name'].startswith('branch_') else 200)
  assert cfg['sampling']['sample_count']==expected
  assert cfg['sampling']['rho_strata_edges']==[.70,.90,.98,1.02,1.10,1.30]
 assert 'current_path(ARTICLE_CASE_DIR)' in (d/'main.cpp').read_text()
subprocess.run(['bash','-n',str(root/'run_all.sh')],check=True)
print(f'Static checks passed: {len(records)} copied sources and {len(cases)} cases. No numerical program executed.')
