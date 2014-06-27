#include <linear_operators/wilson.hpp>

Wilson::Wilson(
  const double mass, 
  const vector<complex<double> >& boundaryConditions,
  const Lattice* lattice) : LinearOperator::LinearOperator()
{
  // Class constructor - we set the fermion mass, create a pointer to the 
  // lattice and compute the frequently used spin structures used within the
  // Dirac operator.
  this->operatorSize_ 
    = 12 * int(pow(lattice->spatialExtent, 3)) * lattice->temporalExtent;
  this->lattice_ = lattice;

  
  this->mass_ = mass;
  this->hoppingMatrix_ = new HoppingTerm(boundaryConditions, lattice, 1);

  // These should be generated depending on whether there's a hopping matrix
  // available
  this->evenIndices_ = this->hoppingMatrix_->getEvenIndices();
  this->oddIndices_ = this->hoppingMatrix_->getOddIndices();
}



Wilson::~Wilson()
{
  // Need to determine which member_functions
  
  delete this->hoppingMatrix_;
  
}



VectorXcd Wilson::apply(const VectorXcd& psi)
{
  VectorXcd eta;
  
  eta = ((4 + this->mass_) * psi);
  eta -= (0.5 * this->hoppingMatrix_->apply(psi));
  return eta;
}



VectorXcd Wilson::applyHermitian(const VectorXcd& psi)
{
  
  return psi;
}



VectorXcd Wilson::makeHermitian(const VectorXcd& psi)
{
  
  return psi;
}



VectorXcd Wilson::applyEvenEvenInv(const VectorXcd& psi)
{
  // Invert the even diagonal piece
  
  return (psi / ((1 + (3 / this->lattice_->chi())) + this->mass_));
}



VectorXcd Wilson::applyOddOdd(const VectorXcd& psi)
{
  // Invert the even diagonal piece
  
  return psi;
}



VectorXcd Wilson::applyEvenOdd(const VectorXcd& psi)
{
  
  return psi;
}



VectorXcd Wilson::applyOddEven(const VectorXcd& psi)
{
  
  VectorXcd eta = 0.5 * psi;
  return eta;
}