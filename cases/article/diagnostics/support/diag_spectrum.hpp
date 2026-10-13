#pragma once
namespace diag {
using namespace p3official;
inline void ck(PetscErrorCode e){spectral::detail::check(e,"diagnostic PETSc/SLEPc call");}
struct Shell {
 spectral::MatrixFreePicardOperator& base; spectral::detail::ShellContext* original{}; Mat mat{}; Vec az{},w{}; bool correction=false;
 static PetscErrorCode mult(Mat m,Vec x,Vec y){Shell*s;PetscFunctionBeginUser;PetscCall(MatShellGetContext(m,&s));PetscCall(MatMult(s->base.mat(),x,y));if(s->correction){PetscScalar dot;PetscCall(VecDot(s->w,x,&dot));PetscCall(VecAXPY(y,dot,s->az));}PetscFunctionReturn(PETSC_SUCCESS);}
 static PetscErrorCode trans(Mat m,Vec x,Vec y){Shell*s;PetscFunctionBeginUser;PetscCall(MatShellGetContext(m,&s));PetscCall(KSPSolveTranspose(s->original->a_solver,x,s->original->work));PetscCall(MatMultTranspose(s->original->b_matrix,s->original->work,y));if(s->correction){PetscScalar dot;PetscCall(VecDot(s->az,x,&dot));PetscCall(VecAXPY(y,dot,s->w));}PetscFunctionReturn(PETSC_SUCCESS);}
 Shell(spectral::MatrixFreePicardOperator& b):base(b){ck(MatShellGetContext(base.mat(),&original));auto n=static_cast<PetscInt>(base.cells());ck(MatCreateShell(PETSC_COMM_SELF,n,n,n,n,this,&mat));ck(MatShellSetOperation(mat,MATOP_MULT,(void(*)(void))mult));ck(MatShellSetOperation(mat,MATOP_MULT_TRANSPOSE,(void(*)(void))trans));ck(MatCreateVecs(mat,&az,nullptr));ck(VecDuplicate(az,&w));}
 void rankone(const DenseVector& z,const DenseVector& line){PetscScalar* a;ck(VecGetArray(w,&a));std::copy(z.begin(),z.end(),a);ck(VecRestoreArray(w,&a));ck(KSPSolve(original->a_solver,w,az));ck(VecGetArray(w,&a));std::copy(line.begin(),line.end(),a);ck(VecRestoreArray(w,&a));correction=true;}
 ~Shell(){MatDestroy(&mat);VecDestroy(&az);VecDestroy(&w);}
};
struct Pair{std::complex<double> mu;double error;int index;};
struct Spectrum{std::vector<Pair> pairs;DenseVector v,vi,u;int nconv{},iterations{},reason{};double kappa=std::numeric_limits<double>::quiet_NaN();double rho()const{return std::abs(pairs.at(0).mu);}};
inline DenseVector copyvec(Vec v){const PetscScalar*a;ck(VecGetArrayRead(v,&a));PetscInt n;ck(VecGetSize(v,&n));DenseVector r(a,a+n);ck(VecRestoreArrayRead(v,&a));return r;}
inline Spectrum spectrum(Mat mat,bool two=true){
 EPS e{};ck(EPSCreate(PETSC_COMM_SELF,&e));ck(EPSSetOperators(e,mat,nullptr));ck(EPSSetProblemType(e,EPS_NHEP));ck(EPSSetType(e,EPSKRYLOVSCHUR));ck(EPSSetWhichEigenpairs(e,EPS_LARGEST_MAGNITUDE));ck(EPSSetDimensions(e,8,PETSC_DEFAULT,PETSC_DEFAULT));ck(EPSSetTolerances(e,1e-11,3000));ck(EPSSetTwoSided(e,two?PETSC_TRUE:PETSC_FALSE));ck(EPSSolve(e));Spectrum s;ck(EPSGetConverged(e,&s.nconv));ck(EPSGetIterationNumber(e,&s.iterations));EPSConvergedReason r;ck(EPSGetConvergedReason(e,&r));s.reason=r;
 for(int i=0;i<s.nconv;++i){PetscScalar re,im;PetscReal err;ck(EPSGetEigenvalue(e,i,&re,&im));ck(EPSComputeError(e,i,EPS_ERROR_RELATIVE,&err));s.pairs.push_back({{double(re),double(im)},double(err),i});}
 std::sort(s.pairs.begin(),s.pairs.end(),[](auto&a,auto&b){return abs(a.mu)>abs(b.mu);});if(s.pairs.empty()){EPSDestroy(&e);throw std::runtime_error("No converged eigenpairs");}
 Vec v,vi,u;ck(MatCreateVecs(mat,&v,nullptr));ck(VecDuplicate(v,&vi));ck(VecDuplicate(v,&u));ck(EPSGetEigenvector(e,s.pairs[0].index,v,vi));s.v=copyvec(v);s.vi=copyvec(vi);
 if(two&&std::abs(s.pairs[0].mu.imag())/s.rho()<=1e-10){ck(EPSGetLeftEigenvector(e,s.pairs[0].index,u,nullptr));s.u=copyvec(u);PetscReal vn,un;PetscScalar uv;ck(VecNorm(v,NORM_2,&vn));ck(VecNorm(u,NORM_2,&un));ck(VecDot(u,v,&uv));s.kappa=vn*un/std::abs(uv);}
 VecDestroy(&v);VecDestroy(&vi);VecDestroy(&u);EPSDestroy(&e);return s;
}
inline double dot(const DenseVector&a,const DenseVector&b){double r=0;for(size_t i=0;i<a.size();++i)r+=a[i]*b[i];return r;}
inline double norm(const DenseVector&a){return sqrt(dot(a,a));}
inline void normalise(DenseVector&a){double n=norm(a);if(n==0)throw std::runtime_error("Zero direction");for(auto&x:a)x/=n;}
inline std::string number(double x){if(!std::isfinite(x))return "NA";std::ostringstream o;o<<std::setprecision(17)<<x;return o.str();}
}
