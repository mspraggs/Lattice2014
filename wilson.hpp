#ifndef WILSON_HPP
#define WILSON_HPP

#include <Eigen/Dense>

#include <complex>

#include <omp.h>

#include <lattice.hpp>
#include <utils.hpp>
#include <linear_operators/linear_operator.hpp>
#include <linear_operators/hopping_term.hpp>

using namespace Eigen;
using namespace std;

class Wilson : public LinearOperator
{

public:
  Wilson(const double mass,
    
	 const vector<complex<double> >& boundaryConditions,
	 const Lattice* lattice);
  ~Wilson();

  VectorXcd apply(const VectorXcd& psi);
  VectorXcd applyHermitian(const VectorXcd& psi);
  VectorXcd makeHermitian(const VectorXcd& psi);

  VectorXcd applyEvenEvenInv(const VectorXcd& psi);
  VectorXcd applyOddOdd(const VectorXcd& psi);
  VectorXcd applyEvenOdd(const VectorXcd& psi);
  VectorXcd applyOddEven(const VectorXcd& psi);

private:
  // Member variables
  
  HoppingTerm* hoppingMatrix_;
  
  double mass_;
    
  const Lattice* lattice_;
  vector<vector<complex<double> > > boundaryConditions_;
};

#endif