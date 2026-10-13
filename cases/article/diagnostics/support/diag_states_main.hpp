#include "diag_state.hpp"
using namespace p3official;
struct DiagCase{std::string mesh;size_t id;Real alpha,theta,c,rho;};
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("usage: diagnostic config selected.tsv output_directory");gss::SlepcContext slepc;auto settings=load_settings(argv[1]);auto config=gss::YAMLConfig::load(argv[1]);auto study=ps::load_study_settings(argv[1]);study.random_mesh_seed=settings.mesh_seed;double delta=config.solver.source_regularisation_delta;
 std::ifstream input(argv[2]);std::vector<DiagCase> cases;DiagCase q;while(input>>q.mesh>>q.id>>q.alpha>>q.theta>>q.c>>q.rho)cases.push_back(q);
 auto stream=[&](std::string name,std::string header){std::ofstream f(std::string(argv[3])+"/"+name);f<<std::setprecision(17)<<header<<'\n';return f;};
#if DIAG_KIND==2
 auto floorout=stream("d2_floor.csv","mesh_family,case_id,state,M2,tau,r_min,min_cell,min_x,min_y,eta_max,n_floor,n_max_ties");
 auto tiesout=stream("d2_maximisers.csv","mesh_family,case_id,state,max_cell,x,y,coordinate_within_first_ten");
 auto wout=stream("d2_w.csv","mesh_family,case_id,state,max_cell,column,w");
 auto fdout=stream("d2_directional.csv","mesh_family,case_id,state,direction,epsilon_relative,epsilon,err_fr,err_tot,n_floor,wv_forward,wv_backward");
 auto spout=stream("d2_spectrum.csv","mesh_family,case_id,state,max_cell,rho_reference,rho_fr,rho_tot,drho,drho_rel,first_order_dmu,nconv_fr,nconv_tot,residual_tot");
#else
 auto caseout=stream("d4_cases.csv","mesh_family,case_id,state,rho,mu_re,mu_im,IPR,n90,cells,nconv,iterations");
 auto topout=stream("d4_topcells.csv","mesh_family,case_id,state,kind,rank,cell,x,y,abs_v,r_P,area,aspect_ratio,stationary_distance,stationary_k,row_norm");
#endif
 const vgt::gradients::polynomial::WLSQR2D method{vgt::gradients::WLSOptions2D{vgt::gradients::WLSWeighting2D::InverseDistanceSquared}};
 for(auto mesh:settings.mesh_families){auto mc=ps::make_mesh_config(config.mesh,mesh,settings.mesh_seed,settings.nominal_nx,settings.nominal_ny);gss::Fadeev cf{settings.calibration_alpha,settings.calibration_theta_degrees,Real{1}};auto ws=ps::make_workspace(mc,cf,study.jacobian_stencil_radius);auto bs=evaluate_base_source(ws,ws.exact,cf,method,delta,study);auto comp=make_compensator(ws,ws.exact,bs,Real{1});ModifiedSource src{ws,cf,method,delta,study,Real{1},comp};DiagJacobian jac{ws,method,delta,study,settings.source_derivative};auto vals=jac.values(ws.exact,src);auto b=sparse::to_csr(ws.jacobian_pattern,vals);spectral::MatrixFreePicardOperator op{ws.poisson.matrix,b};diag::Shell shell{op};auto&top=ws.mesh.topology();
#if DIAG_KIND==4
 auto aspects=diag::aspect(ws);
#endif
 for(auto row:cases)if(row.mesh==mesh){auto local=study;local.alpha=row.alpha;local.theta_degrees=row.theta;gss::Fadeev f{row.alpha,row.theta,Real{1}};refresh_case_data(ws,f);auto base=evaluate_base_source(ws,ws.exact,f,method,delta,local);auto comp2=make_compensator(ws,ws.exact,base,row.c);ModifiedSource source{ws,f,method,delta,local,row.c,comp2};
 for(int st=0;st<2;st++){std::string state_name=st?"Psi0":"PsiF";auto state=st?perturbed_state(ws,ws.exact,1e-4):ws.exact;auto values=jac.values(state,source);op.update_b_values(values);shell.correction=false;auto spec=diag::spectrum(shell.mat);auto fl=diag::floor(jac,state,delta);auto prefix=mesh+","+std::to_string(row.id)+","+state_name+",";
#if DIAG_KIND==2
 floorout<<prefix<<fl.M2<<','<<fl.tau<<','<<fl.rmin<<','<<fl.min<<','<<top.coord_x[fl.min]<<','<<top.coord_y[fl.min]<<','<<fl.eta<<','<<fl.count<<','<<fl.ties.size()<<'\n';floorout.flush();
 auto z=diag::zcolumn(jac,state,row.c*f.nonlinearity_scale());std::vector<DenseVector> lines;
 for(size_t t=0;t<fl.ties.size();t++){size_t cell=fl.ties[t];lines.push_back(diag::wline(jac,fl,cell,delta));tiesout<<prefix<<cell<<','<<(t<10?diag::number(top.coord_x[cell]):"NA")<<','<<(t<10?diag::number(top.coord_y[cell]):"NA")<<','<<(t<10)<<'\n';for(size_t i=0;i<state.size();i++)if(lines.back()[i]!=0)wout<<prefix<<cell<<','<<i<<','<<lines.back()[i]<<'\n';
 shell.rankone(z,lines.back());auto total=diag::spectrum(shell.mat,false);double correction=NAN;if(!spec.u.empty()){auto az=diag::copyvec(shell.az);correction=diag::dot(spec.u,az)*diag::dot(lines.back(),spec.v)/diag::dot(spec.u,spec.v);}double drho=total.rho()-spec.rho();spout<<prefix<<cell<<','<<row.rho<<','<<spec.rho()<<','<<total.rho()<<','<<drho<<','<<drho/spec.rho()<<','<<diag::number(correction)<<','<<spec.nconv<<','<<total.nconv<<','<<total.pairs[0].error<<'\n';spout.flush();}
 DenseVector random(state.size());std::mt19937_64 rng(20261013+row.id);std::uniform_real_distribution<double>dist(-1,1);for(auto&x:random)x=dist(rng);double sn=0;for(auto x:state)sn=std::max(sn,std::abs(x));auto csr=sparse::to_csr(ws.jacobian_pattern,values);
 for(int di=0;di<2;di++){auto v=di?random:spec.v;diag::normalise(v);auto fr=csr.multiply(v);double wp=-INFINITY,wm=INFINITY;for(auto&w:lines){double vdot=diag::dot(w,v);wp=std::max(wp,vdot);wm=std::min(wm,vdot);}auto tot=fr;for(size_t i=0;i<tot.size();i++)tot[i]+=z[i]*(wp+wm)/2;
 for(double er:{1e-4,1e-5,1e-6,1e-7}){double eps=er*sn;auto plus=state,minus=state;for(size_t i=0;i<state.size();i++){plus[i]+=eps*v[i];minus[i]-=eps*v[i];}auto bp=source(plus),bm=source(minus);double nf=0,ef=0,et=0;for(size_t i=0;i<state.size();i++){double fd=(bp[i]-bm[i])/(2*eps);nf+=fd*fd;ef+=(fd-fr[i])*(fd-fr[i]);et+=(fd-tot[i])*(fd-tot[i]);}fdout<<prefix<<(di?"random":"dominant_real_part")<<','<<er<<','<<eps<<','<<diag::number(sqrt(ef/nf))<<','<<diag::number(sqrt(et/nf))<<','<<fl.count<<','<<wp<<','<<wm<<'\n';}fdout.flush();}
#else
 DenseVector mass(state.size());double sum=0,sum2=0;for(size_t i=0;i<mass.size();i++){mass[i]=spec.v[i]*spec.v[i]+spec.vi[i]*spec.vi[i];sum+=mass[i];sum2+=mass[i]*mass[i];}std::vector<size_t>order(mass.size());std::iota(order.begin(),order.end(),0);std::sort(order.begin(),order.end(),[&](auto a,auto b){return mass[a]>mass[b];});double accum=0;size_t n90=0;while(n90<order.size()&&accum<.9*sum)accum+=mass[order[n90++]];caseout<<prefix<<spec.rho()<<','<<spec.pairs[0].mu.real()<<','<<spec.pairs[0].mu.imag()<<','<<sum2/(sum*sum)<<','<<n90<<','<<mass.size()<<','<<spec.nconv<<','<<spec.iterations<<'\n';caseout.flush();auto rn=diag::rownorm(ws.jacobian_pattern,values);
 auto writecell=[&](size_t i,std::string kind,size_t rank){auto [distance,k]=diag::stationary(top.coord_x[i],top.coord_y[i],row.theta);topout<<prefix<<kind<<','<<rank<<','<<i<<','<<top.coord_x[i]<<','<<top.coord_y[i]<<','<<sqrt(mass[i])<<','<<sqrt(fl.g[i]/fl.M2)<<','<<ws.volumes[i]<<','<<diag::number(aspects[i])<<','<<distance<<','<<k<<','<<rn[i]<<'\n';};for(size_t k=0;k<10;k++)writecell(order[k],"top_abs_v",k+1);writecell(std::max_element(rn.begin(),rn.end())-rn.begin(),"max_row_norm",1);topout.flush();
#endif
 }
 std::cout<<mesh<<" case "<<row.id<<" complete"<<std::endl;
 }
 }
 return 0;}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
