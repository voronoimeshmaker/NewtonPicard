from pathlib import Path
import csv,json,hashlib,math,collections
root=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader(p.open()))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks=[]
paths={}
for k in ('D0','D1','D2','D3','D4'):
    d=next(root.glob(k+'_*/'));o=d/(d/'Result/LATEST.txt').read_text().strip();paths[k]=o;m=json.loads((o/'manifest.txt').read_text())
    assert m['exit_code']==0 and not m['protected_files_changed'];assert m['executable_sha256']==sha(root/'build'/('diag_'+k));checks.append(k+': successful exit, executable hash matches, protected snapshot unchanged')
before=json.loads((paths['D0']/'protected_before.json').read_text());changed=[f for f,h in before.items() if not Path(f).exists() or sha(Path(f))!=h];assert not changed,changed
checks.append('All protected files from initial D0 snapshot remain byte-identical after D1-D4')
g=read(paths['D0']/'gate.csv');assert len(g)==60 and all(r['pass']=='1' for r in g)
a=read(paths['D1']/'d1_cases.csv');assert len(a)==3281 and len({(r['mesh_family'],r['case_id']) for r in a})==3281
assert all(int(r['nconv'])>=8 for r in a)
pairs=read(paths['D1']/'d1_eigenpairs.csv');counts=collections.Counter((r['mesh_family'],r['case_id']) for r in pairs)
assert all(counts[r['mesh_family'],r['case_id']]==int(r['nconv']) for r in a)
checks.append('D1: 3281 unique cases; all nconv>=8; every converged pair exported')
sel=json.loads((paths['D2']/'selection.json').read_text())['cases'];assert sum('reference_recovery' in r['selection_reason'] for r in sel)==46
assert {int(r['case_id']) for r in sel if 'random_anomaly' in r['selection_reason']}=={590,1828,1627,983,820,259,1016}
f=read(paths['D2']/'d2_floor.csv');fd=read(paths['D2']/'d2_directional.csv');sp=read(paths['D2']/'d2_spectrum.csv');assert len(f)==2*len(sel)==328;assert len(fd)==len(f)*8;assert len(sp)==sum(int(r['n_max_ties']) for r in f)
checks.append('D2: all mandatory cases, both states, eight directional steps/state and all tied-maximiser spectra exported')
meshes=read(paths['D3']/'d3_mesh_symmetry.csv');assert len(meshes)==12
expected={'L1':6247,'main':25039,'L3':100765,'L4':403259};assert all(int(r['cells'])==expected[r['level']] for r in meshes if r['mesh_family']=='hexagonal')
checks.append('D3: all twelve meshes and four expected hexagonal cell counts')
a=read(paths['D4']/'d4_cases.csv');top=read(paths['D4']/'d4_topcells.csv');assert len(a)==82 and len(top)==82*11
assert all(0<float(r['IPR'])<=1 and 1<=int(r['n90'])<=int(r['cells']) for r in a)
checks.append('D4: 82 states, ten dominant cells plus maximal row cell per state, IPR and participation counts valid')
for k,o in paths.items():
    for file in o.glob('*.csv'):
        for row in read(file):
            assert None not in row and all(v is not None for v in row.values()),str(file)
checks.append('All result CSV rows match their column headers')
(root/'DELIVERY_VALIDATION.json').write_text(json.dumps({'passed':True,'checks':checks},indent=2));print('\n'.join(checks))
