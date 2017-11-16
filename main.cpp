#include <igl/viewer/Viewer.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include<Eigen/SparseQR>
#include <../../spectra-0.5.0/include/MatOp/SparseGenMatProd.h>
#include <../../spectra-0.5.0/include/MatOp/SparseCholesky.h>
#include <../../spectra-0.5.0/include/GenEigsSolver.h>
#include <../../spectra-0.5.0/include/SymGEigsSolver.h>

using namespace std;
using namespace Eigen;


struct Spring
{
  public:
  int first;
  int second;
  double k;
  double L;
  Spring(int p1_index, int p2_index, double length, double stiffness)
  {
    first = p1_index;
    second = p2_index;
    L = length;
    k = stiffness;
  }
};

class Simulation
{
public:
  //Time integration variables
  VectorXd f;
  VectorXd x;
  VectorXd v;
  SparseMatrix<double> gradf;
  SparseMatrix<double> Identity;

  //System Setup Variables
  SparseMatrix<double> InvM; //Inverse Mass Matrix
  SparseMatrix<double> M; //Mass Matrix
  MatrixXd V; //Eigenvectors
  VectorXd D; //Eigenvalues
  vector<Spring> springs;
  MatrixXd Verts;

  //Physics Constants
  double g = 9.8; //gravity (m/ s^2)
  double h = 0.01; //timestep (s)
  double t = 0; //total time


  Simulation();
  void update();
  void drawlinearmodes();

  void getEnergy(double& return_E);
  void getForceVector(VectorXd& return_f, VectorXd& u);
  void getStiffnessMatrix(SparseMatrix<double>& return_K, VectorXd& u);

  //Helper functions
  void VToX(MatrixXd& V, VectorXd& x);
  void XToV(MatrixXd& V, VectorXd& x);
  void generalizedEigenSolve( int reducedSize);
};

void Simulation::generalizedEigenSolve(int reducedSize)
{
  //DOCUMENTATION HERE
  //https://spectralib.org/doc/classspectra_1_1symgeigssolver_3_01scalar_00_01selectionrule_00_01optype_00_01boptype_00_01geigs__cholesky_01_4
  //and here
  //https://github.com/dilevin/GAUSS/blob/a7e924639966728d16702ae19f4149c7baeaa600/src/Core/include/UtilitiesEigen.h#L99

  this->VToX(this->Verts, this->x);
  this->getStiffnessMatrix(this->gradf, this->x);

  Spectra::SparseGenMatProd<double> Aop(this->gradf);
  Spectra::SparseCholesky<double> Bop(this->M);

  Spectra::SymGEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY>
  geigs(&Aop, &Bop, reducedSize, M.rows());

  geigs.init();
  int nconv = geigs.compute();

  if(geigs.info() == Spectra::SUCCESSFUL)
  {
      this->D = geigs.eigenvalues();
      this->V = geigs.eigenvectors();
  }
  else
  {
      cout<<"EIG SOLVE FAILED: "<<endl<<geigs.info()<<endl;
      exit(0);
  }

  cout<<this->D<<endl;
  // cout<<this->V<<endl;

}

void Simulation::getEnergy(double& return_E)
{
  // cout<<"GET ENERGY"<<endl;
  return_E = 0;

  //Potential Energy
  double potE = 0.0;
  for(int i =0; i<this->springs.size(); ++i)
  {
    //E = 0.5*k*(||p1 - p2|| - L)^2
    potE += 0.5 * this->springs[i].k * pow(((this->Verts.row(this->springs[i].first) - this->Verts.row(this->springs[i].second)).norm() - this->springs[i].L), 2);
  }
  cout<<"Potential Energy"<<endl;
  cout<<potE<<endl;

  //Kinetic Energy
  // 1/2 M v^2
  double kineticE = 0.5*this->v.transpose()*this->InvM*this->v;
  cout<<"Kinetic Energy"<<endl;
  cout<<kineticE<<endl;

  return_E += kineticE + potE;

  cout<<return_E<<endl;

  //Gravity E
  //Ignore that
}

void Simulation::getForceVector(VectorXd& return_f, VectorXd& u)
{
  // cout<<"GET FORCE VECTOR"<<endl;
  return_f.setZero();

  //Spring force
  for(int i=0; i<this->springs.size(); ++i)
  {
    Vector3d p1 = u.segment<3>(3*this->springs[i].first);
    Vector3d p2 = u.segment<3>(3*this->springs[i].second);
    double dist = (p1 - p2).norm();
    //local forces p1----p2
    Vector3d force_l = this->springs[i].k*(dist - this->springs[i].L)*(p1 - p2)/dist;

    return_f.segment<3>(3*this->springs[i].first) += -1*force_l;
    return_f.segment<3>(3*this->springs[i].second) += force_l;
  }

  //Gravity
  for(int i=0; i<this->v.size()/3; ++i)
  {
    double invm = this->InvM.coeff(3*i, 3*i);
    if(invm != 0 )
    {
      return_f(3*i + 1) += -1*this->g/invm;
    }
    else
    {
      return_f(3*i + 1) = 0;
    }
  }

}

void Simulation::getStiffnessMatrix(SparseMatrix<double>& return_K, VectorXd& u)
{
  return_K.setZero();
  // cout<<"GET STIFFNESS"<<endl;
  Matrix3d I;
  I<<1, 0, 0,
    0, 1, 0,
    0, 0, 1;

  Matrix3d gradF_l;
  for(int i=0; i<this->springs.size(); ++i)
  {
    Vector3d p1 = u.segment<3>(3*this->springs[i].first);
    Vector3d p2 = u.segment<3>(3*this->springs[i].second);
    double dist = (p1 - p2).norm();
    gradF_l = this->springs[i].k *(I - this->springs[i].L*((I/dist)-((p1-p2)*(p1-p2).transpose() /pow(dist, 3))));

    for(int j =0; j<3; ++j)
    {
      for(int k =0; k<3; ++k)
      {
        return_K.coeffRef(3*this->springs[i].first + j, 3*this->springs[i].first + k) -= gradF_l.coeff(j, k);
        return_K.coeffRef(3*this->springs[i].second +j, 3*this->springs[i].second +k) -= gradF_l.coeff(j, k);
        return_K.coeffRef(3*this->springs[i].first + j, 3*this->springs[i].second +k) += gradF_l.coeff(j, k);
        return_K.coeffRef(3*this->springs[i].second +j, 3*this->springs[i].first + k) += gradF_l.coeff(j, k);
      }
    }
  }
  //grad of gravity is 0

}

void Simulation::VToX(MatrixXd& V, VectorXd& x)
{
  for(int i =0; i<V.rows(); ++i)
  {
    x.segment<3>(3*i) = V.row(i);
  }
}

void Simulation::XToV(MatrixXd& V, VectorXd& x)
{
  for(int i =0; i<V.rows(); ++i)
  {
    V.row(i) = x.segment<3>(3*i);
  }
}

void Simulation::update()
{
  // cout<<"UPDATE"<<endl;
  double energy;


  // //Implicit Euler
  int NEWTONMAXITERS = 100; double error_bound = 1e-5;
  VectorXd x_k = this->x;
  for(int n =0; n< NEWTONMAXITERS; ++n)
  {
    this->getForceVector(this->f, x_k);
    this->getStiffnessMatrix(this->gradf, x_k);
    VectorXd g = x_k - (this->x + this->h*(this->v + this->h*this->InvM*this->f));
    SparseMatrix<double> grad_g = Identity - this->h*this->h*this->InvM*this->gradf;

    // //Solver dx = -1*inv(grad_g)*g
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		sqr.compute(grad_g);
		VectorXd dx = -1*sqr.solve(g);

    x_k += dx;

    if(x_k != x_k)
    {
      cout<<"NANs!!"<<endl;
      cout<<n<<endl;
      break;
    }
    if(g.squaredNorm()< 1e-7)
    {
      cout<<"error sq norm"<<endl;
      cout<<g.squaredNorm()<<endl;
      break;
    }
    if(n>90)
    {
      cout<<"NM not converging."<<endl;
      break;
    }
  }
  this->v = (x_k - this->x)/this->h;
  this->x = x_k;
  this->XToV(this->Verts, this->x);
  this->t +=1;


  // //Verlet
  // //x_t+1 = x_t + h*v_t
  // //v_t+1 = v_t + h*InvM*f(x_t+1)
  // this->x = this->x + this->h*this->v;
  // this->getEnergy(energy);
  // this->getForceVector(this->f, this->x);
  // this->v = this->v + this->h*this->InvM*this->f;
  // this->XToV(this->Verts, this->x);

  return;
}

Simulation::Simulation()
{
  // SET VERTEXES AND SPRINGS
  //-------------------------
  this->Verts = (MatrixXd(8, 3)<<
                    0.0,0.0,0.0, //0 fix
                    0.0,0.0,1.0, //1 fix
                    0.0,1.0,0.0, //2 fix
                    0.0,1.0,1.0, //3 fix
                    1.0,0.0,0.0, //4
                    1.0,0.0,1.0,
                    1.0,1.0,0.0,
                    1.0,1.0,1.0
                  ).finished();

  // this->Verts = (MatrixXd(2, 3)<<
  //                   0.0,0.0,0.0, //0
  //                   0.0,0.0,1.0 //1
  //                 ).finished();

  double k_s = 100;

  this->springs.push_back(Spring(0, 1, (this->Verts.row(0)- this->Verts.row(1)).norm(), k_s));
  this->springs.push_back(Spring(0, 2, (this->Verts.row(0)- this->Verts.row(2)).norm(), k_s));
  this->springs.push_back(Spring(0, 4, (this->Verts.row(0)- this->Verts.row(4)).norm(), k_s));
  this->springs.push_back(Spring(0, 7, (this->Verts.row(0)- this->Verts.row(7)).norm(), k_s));

  this->springs.push_back(Spring(7, 3, (this->Verts.row(7)- this->Verts.row(3)).norm(), k_s));
  this->springs.push_back(Spring(7, 5, (this->Verts.row(7)- this->Verts.row(5)).norm(), k_s));
  this->springs.push_back(Spring(7, 6, (this->Verts.row(7)- this->Verts.row(6)).norm(), k_s));

  this->springs.push_back(Spring(5, 2, (this->Verts.row(5)- this->Verts.row(2)).norm(), k_s));

  this->springs.push_back(Spring(6, 2, (this->Verts.row(6)- this->Verts.row(2)).norm(), k_s));
  this->springs.push_back(Spring(4, 5, (this->Verts.row(4)- this->Verts.row(5)).norm(), k_s));
  this->springs.push_back(Spring(4, 6, (this->Verts.row(4)- this->Verts.row(6)).norm(), k_s));
  this->springs.push_back(Spring(4, 3, (this->Verts.row(4)- this->Verts.row(3)).norm(), k_s));

  this->springs.push_back(Spring(3, 2, (this->Verts.row(3)- this->Verts.row(2)).norm(), k_s));
  this->springs.push_back(Spring(3, 1, (this->Verts.row(3)- this->Verts.row(1)).norm(), k_s));

  this->springs.push_back(Spring(1, 5, (this->Verts.row(1)- this->Verts.row(5)).norm(), k_s));
  this->springs.push_back(Spring(1, 6, (this->Verts.row(1)- this->Verts.row(6)).norm(), k_s));
  //-------------------------

  int numOfVerts = this->Verts.rows();

  //SET MASSES
  //-------------------------
  VectorXd massVector;
  massVector.resize(3*numOfVerts);
  massVector.setZero();
  for(int i=0; i<numOfVerts; ++i){
    double vertex_mass = 1.0;
    massVector(3*i) = vertex_mass;
    massVector(3*i+1) = vertex_mass;
    massVector(3*i+2) = vertex_mass;
  }
  // -Inv mass matrix
  this->InvM.resize(3*numOfVerts, 3*numOfVerts);
  this->InvM.setZero();
  this->M.resize(3*numOfVerts, 3*numOfVerts);
  this->M.setZero();
  for(int i =0; i<3*numOfVerts; ++i){
    this->InvM.coeffRef(i, i) = 1/massVector(i);
    this->M.coeffRef(i, i) = massVector(i);
  }

  //FIX VERTICES
  //TODO Make this a function
  //fix Vertex 0, 1, 2, 3
  InvM.coeffRef(0,0) = 0;
  InvM.coeffRef(1,1) = 0;
  InvM.coeffRef(2,2) = 0;
  InvM.coeffRef(3,3) = 0;
  InvM.coeffRef(4,4) = 0;
  InvM.coeffRef(5,5) = 0;
  InvM.coeffRef(6,6) = 0;
  InvM.coeffRef(7,7) = 0;
  InvM.coeffRef(8,8) = 0;
  InvM.coeffRef(9,9) = 0;
  InvM.coeffRef(10,10) = 0;
  InvM.coeffRef(11,11) = 0;
  //-------------------------

  //Initialize Variables
  //-------------------------
  this->x.resize(3*numOfVerts);
  this->x.setZero();
  this->VToX(this->Verts, this->x);

  this->v.resize(3*numOfVerts);
  this->v.setZero();

  this->f.resize(3*numOfVerts);
  this->f.setZero();

  this->gradf.resize(3*numOfVerts, 3*numOfVerts);
  this->gradf.setZero();

  this->Identity.resize(3*numOfVerts, 3*numOfVerts);
  this->Identity.setIdentity();
  //-------------------------

  //Do Linear Modal Analysis
  // K*x = lambda*M*x
  //-------------------------
  int getThisManyBasisVectors = 20;
  this->generalizedEigenSolve(getThisManyBasisVectors);
  // cout<<"Check that this is Identity"<<endl;
  // cout<<this->V.transpose()*this->M*this->V<<endl;
  //-------------------------

}

void Simulation::drawlinearmodes()
{
  int drawBasisVector = 9+4;
  VectorXd u = this->x + this->V.col(drawBasisVector)*sin(sqrt(fabs(this->D(drawBasisVector)))*this->t*this->h);
  this->XToV(this->Verts, u);
  // cout<<Sim.Verts<<endl;
}

//INITIALIZE SIMULATION
Simulation Sim = Simulation();

bool drawSim(igl::viewer::Viewer& viewer)
{
  viewer.data.clear();

  Sim.update();
  viewer.data.add_points(Sim.Verts, RowVector3d(1,0,0));
  for (int i=0; i<Sim.springs.size(); ++i)
  {
    viewer.data.add_edges(Sim.Verts.row(Sim.springs[i].first), Sim.Verts.row(Sim.springs[i].second), RowVector3d(0, 1, 0));
  }
}

bool drawBasis(igl::viewer::Viewer& viewer)
{
  viewer.data.clear();

  Sim.drawlinearmodes();

  Sim.t+=1;
  // cout<<Sim.t<<endl;

  viewer.data.add_points(Sim.Verts, RowVector3d(1,0,0));
  for (int i=0; i<Sim.springs.size(); ++i)
  {
    viewer.data.add_edges(Sim.Verts.row(Sim.springs[i].first), Sim.Verts.row(Sim.springs[i].second), RowVector3d(0, 1, 0));
  }

}

int main(int argc, char *argv[])
{
  int linearModalBasis = 1;

  igl::viewer::Viewer viewer;
  viewer.core.is_animating = true;
  if(linearModalBasis == 0)
  {
    // Simulate the mesh under gravity
    viewer.callback_pre_draw = &drawSim;

  }
  else{
    viewer.callback_pre_draw = &drawBasis;
  }
  //

  viewer.launch();
}
