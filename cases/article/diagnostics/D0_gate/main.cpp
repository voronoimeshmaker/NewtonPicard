#include "../../support/official.hpp"
using namespace p3official;
struct Case {std::string mesh;std::size_t id;Real alpha,theta,c,rho;};
int main(int argc,char**argv){
 try {
  if(argc!=4)throw std::runtime_error("usage: diag_D0 config.yaml selected.tsv output.csv");
  gss::SlepcContext slepc;
  auto settings=load_settings(argv[1]);auto config=gss::YAMLConfig::load(argv[1]);auto study=ps::load_study_settings(argv[1]);
  study.random_mesh_seed=settings.mesh_seed;study.requested_eigenpairs=settings.requested_eigenpairs;study.eps_tolerance=settings.eps_tolerance;study.eps_max_iterations=settings.eps_max_iterations;
  std::ifstream in(argv[2]);std::vector<Case> cases;Case q;while(in>>q.mesh>>q.id>>q.alpha>>q.theta>>q.c>>q.rho)cases.push_back(q);
  if(cases.size()!=60)throw std::runtime_error("D0 requires 60 selected cases");
  std::ofstream out(argv[3]);out<<std::setprecision(17)<<"mesh_family,case_id,c,rho_csv,rho_new,relative_error,nconv,eps_iterations,eigenpair_residual,seconds,pass\n";
  const vgt::gradients::polynomial::WLSQR2D method{vgt::gradients::WLSOptions2D{vgt::gradients::WLSWeighting2D::InverseDistanceSquared}};
  bool passed=true;
  for(auto mesh:settings.mesh_families){
   auto mc=ps::make_mesh_config(config.mesh,mesh,settings.mesh_seed,settings.nominal_nx,settings.nominal_ny);
   gss::Fadeev cf{settings.calibration_alpha,settings.calibration_theta_degrees,Real{1}};
   auto ws=ps::make_workspace(mc,cf,study.jacobian_stencil_radius);
   auto bs=evaluate_base_source(ws,ws.exact,cf,method,config.solver.source_regularisation_delta,study);
   auto comp=make_compensator(ws,ws.exact,bs,Real{1});
   ModifiedSource src{ws,cf,method,config.solver.source_regularisation_delta,study,Real{1},comp};
   SourceJacobianAssembler jac{ws,method,config.solver.source_regularisation_delta,study,settings.source_derivative};
   auto vals=jac.values(ws.exact,src);auto b=sparse::to_csr(ws.jacobian_pattern,vals);SpectralEvaluator eig{ws.poisson.matrix,b,settings};
   for(auto row:cases)if(row.mesh==mesh){
    auto start=std::chrono::steady_clock::now();auto local=study;local.alpha=row.alpha;local.theta_degrees=row.theta;
    gss::Fadeev f{row.alpha,row.theta,Real{1}};refresh_case_data(ws,f);
    auto base=evaluate_base_source(ws,ws.exact,f,method,config.solver.source_regularisation_delta,local);
    auto comp2=make_compensator(ws,ws.exact,base,row.c);ModifiedSource source{ws,f,method,config.solver.source_regularisation_delta,local,row.c,comp2};
    auto result=eig.evaluate(ws,ws.exact,source,jac);auto err=std::abs(result.rho-row.rho)/row.rho;bool ok=std::isfinite(err)&&err<=1e-10;passed=passed&&ok;
    out<<mesh<<','<<row.id<<','<<row.c<<','<<row.rho<<','<<result.rho<<','<<err<<','<<result.converged_eigenpairs<<','<<result.eps_iterations<<','<<result.max_relative_error<<','<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<','<<ok<<'\n';out.flush();
    std::cout<<mesh<<" case "<<row.id<<" relative_error="<<std::setprecision(17)<<err<<" pass="<<ok<<std::endl;
   }
  }
  return passed?0:2;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
}
