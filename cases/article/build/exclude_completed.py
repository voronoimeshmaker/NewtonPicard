from pathlib import Path
import csv,json,os,subprocess
r=Path('/home/jflavio/Programas/NewtonPicard/cases/article')
completed=['01_fadeev_spatial','02_newton','03_picard','04_initial_rho']
for name in completed:
 run=sorted((r/name/'Result').glob('run_*'))[-1]
 if name=='01_fadeev_spatial':
  assert 'run_state=completed' in (run/'manifest.txt').read_text()
 else:
  rows=list(csv.DictReader((run/'cases.csv').open())); assert len(rows)==12000
  assert len({(x['method'],x['mesh_family'],x['case_id']) for x in rows})==12000
  assert (run/'summary.csv').is_file()
  assert 'Cases recorded: 12000' in (run/'console.log').read_text()
p=r/'run_all.sh'; text=p.read_text()
anchor='if [[ "$mode" == list ]]; then'
block='''# Campanhas integralmente concluidas e conferidas em 2026-10-10.
# Excluidas inclusive de --case e dos wrappers individuais para evitar repeticao.
pending=()
for c in "${cases[@]}"; do
 case "$c" in
 01_fadeev_spatial|02_newton|03_picard|04_initial_rho)
  if [[ -n "$selected" ]]; then
   echo "Caso $c ja concluido integralmente; nova execucao bloqueada." >&2
   exit 2
  fi
  ;;
 *) pending+=("$c");;
 esac
done
cases=("${pending[@]}")
'''
assert '# Campanhas integralmente concluidas' not in text
text=text.replace(anchor,block+anchor,1)
tmp=p.with_suffix('.sh.tmp');tmp.write_text(text);tmp.chmod(p.stat().st_mode);subprocess.run(['bash','-n',str(tmp)],check=True);os.replace(tmp,p)
print('Excluded:',', '.join(completed))
subprocess.run(['bash',str(p),'--list'],check=True)
q=subprocess.run(['bash',str(p),'--case','02_newton'],capture_output=True,text=True);assert q.returncode==2;print(q.stderr.strip())