import numpy as np

import pyQCD

@types(mass_="double", hoppingMatrix_='HoppingTerm*')
class Wilson(object):

    @types(mass="const double")
    def __init__(self, mass):

        self.mass_ = mass
        self.hoppingMatrix_ = HoppingTerm(1)

    @types(psi='VectorXcd', eta='VectorXcd')
    def apply(self, psi):

        eta = (4 + self.mass_) * psi
        eta -= 0.5 * self.hoppingMatrix_.apply(psi)
        return eta

    def make_hermitian(self, psi):
        return psi

    def apply_hermitian(self, psi):
        return psi

    def apply_even_even_inv(self, psi):
        return psi / (1 + 3 / self.lattice_.chi() + self.mass_)

    def apply_odd_odd(self, psi):
        return psi

    def apply_even_odd(self, psi):
        return psi

    def apply_odd_even(self, psi):
        RAW_CPP("VectorXcd eta = 0.5 * psi;")
        return eta
