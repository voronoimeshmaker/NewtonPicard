#!/usr/bin/env bash
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
mode=run
selected=
while (($#)); do
 case "$1" in
 --build-only) mode=build; shift;;
 --list) mode=list; shift;;
 --case) selected="${2:?Informe o caso}"; shift 2;;
 *) echo "Uso: bash run_all.sh [--list|--build-only] [--case nome]" >&2; exit 2;;
 esac
done
mapfile -t cases < <(python3 -c 'import json,sys; print("\n".join(c["name"] for c in json.load(open(sys.argv[1]))))' "$root/cases.json")
if [[ -n "$selected" ]]; then
 found=0
 for c in "${cases[@]}"; do [[ "$c" == "$selected" ]] && found=1; done
 [[ "$found" == 1 ]] || { echo "Caso desconhecido: $selected" >&2; exit 2; }
 cases=("$selected")
fi
# Campanhas integralmente concluidas e conferidas em 2026-10-10.
# Excluidas inclusive de --case e dos wrappers individuais para evitar repeticao.
pending=()
for c in "${cases[@]}"; do
 case "$c" in
 01_fadeev_spatial|02_newton|03_picard|04_initial_rho|05_updated_rho|amplitude_1e-4_newton|amplitude_1e-4_picard|amplitude_1e-4_initial_rho|amplitude_1e-4_updated_rho|amplitude_1e-2_newton|amplitude_1e-2_picard|amplitude_1e-2_initial_rho|amplitude_1e-2_updated_rho|amplitude_1e-1_newton|amplitude_1e-1_picard|amplitude_1e-1_initial_rho|amplitude_1e-1_updated_rho|amplitude_1e0_newton|amplitude_1e0_picard|amplitude_1e0_initial_rho|amplitude_1e0_updated_rho|branch_L1_picard)
  if [[ -n "$selected" ]]; then
   echo "Caso $c ja concluido integralmente; nova execucao bloqueada." >&2
   exit 2
  fi
  ;;
 *) pending+=("$c");;
 esac
done
cases=("${pending[@]}")
mapfile -t cases < <(python3 "$root/pending_cases.py" "${cases[@]}")
if ((${#cases[@]} == 0)); then
 echo "Nenhuma campanha pendente." >&2
 exit 0
fi
if [[ "$mode" == list ]]; then printf '%s\n' "${cases[@]}"; exit 0; fi
mkdir -p "$root/build" "$root/build/tmp"
export TMPDIR="$root/build/tmp"
exec 9>"$root/build/runner.lock"
flock -n 9 || { echo "Outro launcher está ativo neste pacote" >&2; exit 2; }
python3 "$root/audit_static.py" > "$root/build/static_check.log" 2>&1 || { cat "$root/build/static_check.log"; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
# Use the same pkg-config environment for PETSc and SLEPc; do not mix installations.
cmake -S "$root" -B "$root/build" -DCMAKE_BUILD_TYPE=Release > "$root/build/configure.log" 2>&1 || { cat "$root/build/configure.log"; exit 1; }
for c in "${cases[@]}"; do
 mkdir -p "$root/$c/logs" "$root/$c/tmp"
 cmake --build "$root/build" --target "article_$c" --parallel "${BUILD_JOBS:-2}" > "$root/$c/logs/build.log" 2>&1 || { cat "$root/$c/logs/build.log"; exit 1; }
done
[[ "$mode" == build ]] && exit 0
# Serial campaigns: no numerical program is invoked by --build-only or --list.
command -v taskset >/dev/null || { echo "taskset não encontrado" >&2; exit 1; }
cpu="${ARTICLE_CPU:-$(python3 -c 'import os; print(min(os.sched_getaffinity(0)))')}"
failed=0
for c in "${cases[@]}"; do
 echo "Executando $c"
 stamp="$(date -u +%Y%m%dT%H%M%S)_$$"
 if (cd "$root/$c"; export TMPDIR="$root/$c/tmp"; taskset --cpu-list "$cpu" "$root/build/article_$c") 2>&1 | tee "$root/$c/logs/run_${stamp}.log"; then
  echo "Concluido: $c"
 else
  echo "FALHA: $c; consulte o log. Continuando as campanhas pendentes independentes." >&2
  failed=1
 fi
done

exit "$failed"
