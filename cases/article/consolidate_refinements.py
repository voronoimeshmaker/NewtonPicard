from pathlib import Path
import csv,json,datetime,time,os,sys
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'consolidated_refinements'
OUT.mkdir(exist_ok=True)
def consolidate(name):
    selected={}
    for run in sorted((ROOT/name/'Result').glob('run_*')):
        if not (run/'summary.csv').exists() or not (run/'cases.csv').exists():continue
        with (run/'cases.csv').open() as f: rows=list(csv.DictReader(f))
        for mesh in ['cartesian','hexagonal','random']:
            group=[q for q in rows if q.get('method')=='WLSQR' and q.get('mesh_family')==mesh]
            if len(group)==100 and {q.get('case_id') for q in group}=={str(i) for i in range(1,101)}:
                selected[mesh]=(run,group)
    if len(selected)!=3:raise RuntimeError(name+': incomplete families '+str(list(selected)))
    rows=[];sources={};summary=[]
    for mesh,(run,group) in selected.items():
        rows.extend(group);sources[mesh]=str(run)
        summary.append(dict(mesh_family=mesh,cases=len(group),converged=sum(q['converged']=='1' for q in group),converged_to_reference=sum(q['converged_to_reference']=='1' for q in group),nonconverged=sum(q['converged']=='0' for q in group)))
    assert len({(q['method'],q['mesh_family'],q['case_id']) for q in rows})==300
    dest=OUT/name;dest.mkdir(exist_ok=True)
    for filename,data in [('cases.csv',rows),('counts.csv',summary)]:
        tmp=dest/(filename+'.tmp')
        with tmp.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
        os.replace(tmp,dest/filename)
    (dest/'sources.json').write_text(json.dumps(dict(sources=sources,note='Random mesh resumed with local exact-clipping fallback; other families preserved. Nonconvergence and case exceptions remain recorded scientific outcomes.'),indent=2)+'\n')
    print(name,'COMPLETE 300 unique cases',flush=True)
    return summary
if __name__=='__main__':
    try:
        consolidate('branch_L3_picard')
        if len(sys.argv)>1:
            pid=int(sys.argv[1])
            while Path('/proc') .joinpath(str(pid)).exists():time.sleep(10)
        consolidate('branch_L4_picard')
        (OUT/'STATUS.txt').write_text('COMPLETE: L3 and L4, 300 unique records each. Original results preserved. UTC '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'\n')
    except Exception as e:
        (OUT/'STATUS.txt').write_text('INCOMPLETE: '+str(e)+'\n');raise
