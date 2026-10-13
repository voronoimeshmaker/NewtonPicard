from pathlib import Path
import hashlib,json,re,subprocess,datetime

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def capture(root):
    files=set()
    linked=subprocess.check_output(['ldd',str(root/'build/diag_D0')],text=True)
    files.update(Path(x).resolve() for x in re.findall(r'(/[^\s()]+)',linked) if Path(x).is_file())
    for name in ('PETSc','SLEPc'):
        flags=subprocess.check_output(['pkg-config','--libs',name],text=True)
        dirs=re.findall(r'-L(\S+)',flags)
        for d in dirs:
            for pattern in ('libpetsc.so*','libslepc.so*'):
                files.update(x.resolve() for x in Path(d).glob(pattern) if x.is_file())
    for d in ('/usr/local/lib','/usr/lib/x86_64-linux-gnu'):
        for pat in ('libadolc.so*','libVoronoiMeshMaker*','libvmm*','libVoroGrad*'):
            files.update(x.resolve() for x in Path(d).glob(pat) if x.is_file())
    trees={}
    for d in ('/usr/local/include/VoronoiMeshMaker','/usr/local/include/VoroGradToolkit','/usr/include/adolc','/usr/include/CGAL'):
        base=Path(d)
        if base.exists():
            entries={str(x.relative_to(base)):digest(x) for x in sorted(base.rglob('*')) if x.is_file()}
            trees[d]={'files':len(entries),'tree_sha256':hashlib.sha256(json.dumps(entries,sort_keys=True).encode()).hexdigest()}
    code={str(x.relative_to(root)):digest(x) for x in root.rglob('*') if x.is_file() and not any(part in ('Result','build','__pycache__') for part in x.relative_to(root).parts)}
    return {'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'note':'Post-execution verification; original execution provenance is in manifest.txt. Header tree digest hashes the sorted JSON mapping of relative path to SHA256.','library_sha256':{str(x):digest(x) for x in sorted(files)},'header_trees':trees,'diagnostic_sources_sha256':code}
if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];data=capture(root)
    for folder in root.glob('D*_*'):
        latest=folder/'Result/LATEST.txt'
        if latest.exists():
            out=folder/latest.read_text().strip();(out/'provenance_verified.json').write_text(json.dumps(data,indent=2))
            mf=out/'manifest.txt';m=json.loads(mf.read_text());m['post_execution_provenance_file']='provenance_verified.json';m['post_execution_provenance_sha256']=digest(out/'provenance_verified.json');mf.write_text(json.dumps(m,indent=2))
