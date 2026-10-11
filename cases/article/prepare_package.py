from pathlib import Path
import shutil, hashlib, json, re
src=Path('/home/jflavio/Programas/GSMaker/Old/GSMakerLib')
dst=Path('/home/jflavio/Programas/NewtonPicard/cases/article')
dst.mkdir(parents=True,exist_ok=True)
manifest=[]
def copy(a,b,transform=None):
 data=a.read_bytes(); b.parent.mkdir(parents=True,exist_ok=True)
 b.write_bytes(data if transform is None else transform(data.decode()).encode())
 manifest.append(dict(source=str(a),destination=str(b.relative_to(dst)),source_sha256=hashlib.sha256(data).hexdigest(),copied_sha256=hashlib.sha256(b.read_bytes()).hexdigest()))
for tree in ['include','src']:
 for a in (src/'GSMakerLib'/tree).rglob('*'):
  if a.is_file() and a.suffix in ['.hpp','.h','.tpp','.cpp','.cc']:
   copy(a,dst/'support/library'/tree/a.relative_to(src/'GSMakerLib'/tree))
for name in ['ConsoleTee.hpp','FadeevNonlinearStudyCommon.hpp','SparseJacobian.hpp','MatrixFreePicardOperator.hpp']:
 copy(src/'studies/Common'/name,dst/'support/common'/name)
# Relocate includes only; keep the two historical numerical drivers distinct.
for orig,name in [('NewtonOnly10000','official'),('CampanhasComplementares/Common','complementary')]:
 a=src/'papers/PaperIII'/orig/'OfficialFourStrategies10000Common.hpp'
 copy(a,dst/'support'/f'{name}.hpp',lambda s:re.sub(r'"(?:\.\./)+studies/Common/', '"common/',s))
copy(src/'papers/PaperIII/FadeevConvergence/pap_p3_fadeev_convergence.cpp',dst/'support/spatial.cpp',lambda s:s.replace('../../../studies/Common/','common/').replace('int main(', 'int spatial_main('))
strategies={'newton':'NewtonOnly','picard':'PicardOnly','initial_rho':'MixedRhoInitial','updated_rho':'MixedRhoEveryIteration'}
base=(src/'papers/PaperIII/NewtonOnly10000/config.yaml').read_text()
cases=[]
def case(name,kind,cfg,description):
 d=dst/name;d.mkdir(exist_ok=True)
 (d/'config.yaml').write_text(cfg)
 if kind=='spatial': call='spatial_main()';inc='spatial.cpp'
 else:
  header,method=kind.split(':');inc=header+'.hpp'
  call=f'p3official::run_program(p3official::Strategy::{strategies[method]}, "{name}")'
 (d/'main.cpp').write_text(f'#include "../support/{inc}"\nint main() {{\n  std::filesystem::current_path(ARTICLE_CASE_DIR);\n  return {call};\n}}\n')
 (d/'README.md').write_text(f'# {name}\n\n{description}\n\nConfiguração: `config.yaml`. Execute `bash run.sh` nesta pasta ou de qualquer diretório.\nSaídas numéricas: `Result/run_*/`, logs em `logs/`, temporários em `tmp/`, todos nesta pasta.\nNenhuma campanha foi executada durante a preparação. Consulte `../AUDIT.md`.\n')
 (d/'run.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\ncase_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"\nexec bash "$case_dir/../run_all.sh" --case "$(basename "$case_dir")"\n')
 cases.append(dict(name=name,kind=kind,description=description))
case('01_fadeev_spatial','spatial',(src/'papers/PaperIII/FadeevConvergence/config.yaml').read_text(),'Verificação espacial: fonte analítica prescrita, malhas 48x16, 96x32 e 192x64; três famílias. Não é um teste não linear.')
for i,m in enumerate(strategies,2):
 case(f'{i:02d}_{m}','official:'+m,base,f'Campanha principal {m}: 2000 amostras, três famílias, WLSQR e MGG. As tabelas principais do artigo selecionam WLSQR.')
for amp in ['1e-4','1e-2','1e-1','1e0']:
 cfg=(src/f'papers/PaperIII/CampanhasComplementares/A-PalpiteInicial/amp_{amp}/config.yaml').read_text()
 for m in strategies: case(f'amplitude_{amp}_{m}','complementary:'+m,cfg,f'Variação de amplitude {amp}: 200 amostras, WLSQR, três famílias, estratégia {m}.')
for level in ['L1','L3','L4']:
 cfg=(src/f'papers/PaperIII/CampanhasComplementares/B-EscalaDeMalha/{level}/config.yaml').read_text()
 # Article branch table concerns Picard; retain its three reported resolutions.
 case(f'branch_{level}_picard','complementary:picard',cfg,f'Refinamento de ramos, nível histórico {level}, 100 amostras, Picard/WLSQR. Correspondência: L1=128x43; L3=512x171; L4=1024x341.')
(dst/'cases.json').write_text(json.dumps(cases,indent=2))
(dst/'provenance.json').write_text(json.dumps(manifest,indent=2))
(dst/'CMakeLists.txt').write_text('''cmake_minimum_required(VERSION 3.20)
project(ArticleCases LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(POLICY CMP0167)
 cmake_policy(SET CMP0167 NEW)
endif()
find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_check_modules(SLEPC REQUIRED IMPORTED_TARGET SLEPc)
find_package(yaml-cpp REQUIRED)
find_package(OpenMP REQUIRED)
find_package(CGAL REQUIRED)
find_package(TBB REQUIRED)
find_package(VoronoiMeshMaker REQUIRED)
if(NOT TARGET VoronoiMeshMaker::VoronoiMeshMaker)
 find_library(LEGACY_VMM NAMES VoronoiMeshMaker REQUIRED)
 add_library(VoronoiMeshMaker::VoronoiMeshMaker INTERFACE IMPORTED)
 set_target_properties(VoronoiMeshMaker::VoronoiMeshMaker PROPERTIES INTERFACE_LINK_LIBRARIES "${LEGACY_VMM};VoronoiMeshMaker::vmm;TBB::tbb")
endif()
find_package(VoroGradToolkit REQUIRED)
find_library(ADOLC_LIB NAMES adolc REQUIRED)
find_path(ADOLC_INCLUDE NAMES adolc/adolc.h REQUIRED)
file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS support/library/src/*.cpp)
add_library(article_historical STATIC ${LIB_SOURCES})
target_include_directories(article_historical PUBLIC support/library/include ${ADOLC_INCLUDE})
target_compile_definitions(article_historical PUBLIC GSMAKER_USE_SLEPC=1 GSMAKER_USE_ADOLC=1 GSMAKER_USE_VOROGRADTOOLKIT=1 GSMAKER_USE_VORONOI_GRADIENT_KIT=1 GSMAKER_USE_CGAL=1)
target_link_libraries(article_historical PUBLIC PkgConfig::PETSC PkgConfig::SLEPC MPI::MPI_CXX yaml-cpp::yaml-cpp ${ADOLC_LIB} VoroGradToolkit::VoroGradToolkit VoronoiMeshMaker::VoronoiMeshMaker CGAL::CGAL TBB::tbb OpenMP::OpenMP_CXX)
file(GLOB CASE_MAINS CONFIGURE_DEPENDS */main.cpp)
foreach(main IN LISTS CASE_MAINS)
 get_filename_component(dir "${main}" DIRECTORY)
 get_filename_component(name "${dir}" NAME)
 add_executable(article_${name} "${main}")
 target_compile_definitions(article_${name} PRIVATE ARTICLE_CASE_DIR="${dir}")
 target_link_libraries(article_${name} PRIVATE article_historical)
endforeach()
''')
(dst/'run_all.sh').write_text('''#!/usr/bin/env bash
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
mapfile -t cases < <(python3 -c 'import json,sys; print("\\n".join(c["name"] for c in json.load(open(sys.argv[1]))))' "$root/cases.json")
if [[ -n "$selected" ]]; then
 found=0
 for c in "${cases[@]}"; do [[ "$c" == "$selected" ]] && found=1; done
 [[ "$found" == 1 ]] || { echo "Caso desconhecido: $selected" >&2; exit 2; }
 cases=("$selected")
fi
if [[ "$mode" == list ]]; then printf '%s\\n' "${cases[@]}"; exit 0; fi
mkdir -p "$root/build" "$root/build/tmp"
export TMPDIR="$root/build/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
# Use the same pkg-config environment for PETSc and SLEPc; do not mix installations.
cmake -S "$root" -B "$root/build" -DCMAKE_BUILD_TYPE=Release > "$root/build/configure.log" 2>&1 || { cat "$root/build/configure.log"; exit 1; }
for c in "${cases[@]}"; do
 mkdir -p "$root/$c/logs" "$root/$c/tmp"
 cmake --build "$root/build" --target "article_$c" --parallel "${BUILD_JOBS:-2}" > "$root/$c/logs/build.log" 2>&1 || { cat "$root/$c/logs/build.log"; exit 1; }
done
[[ "$mode" == build ]] && exit 0
# Serial campaigns: no numerical program is invoked by --build-only or --list.
for c in "${cases[@]}"; do
 echo "Executando $c"
 stamp="$(date -u +%Y%m%dT%H%M%S)_$$"
 (cd "$root/$c"; export TMPDIR="$root/$c/tmp"; "$root/build/article_$c") 2>&1 | tee "$root/$c/logs/run_${stamp}.log"
done
''')
print('Prepared',len(cases),'cases;',len(manifest),'source files')