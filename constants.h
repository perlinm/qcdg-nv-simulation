#pragma once

#include <math.h>
#include <complex.h>

const std::complex<double> j(0.0,1.0);
const double pi = M_PI;

// physical constants in SI
const double mu0_SI = 4*pi*1e-7; // magnetic constant (tesla meters / amp)
const double c_SI = 299792458; // speed of light (meters / second)
const double hbar_SI = 6.582119514e-16; // reduced Planck constant (eV / second)
const double qe_SI = 1.60217646e-19; // charge of electron (coulombs)
const double ge_SI = -1.760859708e11; // gyromagnetic ratio of NV electron (Hz/tesla)
const double gC13_SI = 67.28284e6; // gyromagnetic ratio of C-13 (Hz/tesla)
const double gN15_SI = -27.116e6; // gyromagnetic ratio of N-15 (Hz/tesla)

// physical constants in natural units: e_0 = mu_0 = c = hbar = 1; basic units are s and Hz
const double alpha = 1/137.035999074; // fine structure constant
const double qe = sqrt(4*pi*alpha); // unit electric charge
const double me = 510998.928/hbar_SI; // mass of electron (Hz)
const double ge = -qe/(2*me) * 2.0023193043617; // gyromagnetic ratio of electron (s)
const double gC13 = gC13_SI * ge/ge_SI; // gyromagnetic ratio of C-13 (s)
const double gN15 = gN15_SI * ge/ge_SI; // gyromagnetic ratio of N-15 (s)

// SI values in natural units
const double kHz = 1000; // one kilohertz (s)
const double meter = 1 / c_SI; // one meter (s)
const double nm = 1e-9*meter; // one nanometer (s)
const double coulomb = qe / qe_SI; // one coulomb
const double tesla = ge_SI / ge; // one tesla (Hz^2)
const double gauss = 1e-4*tesla; // one gauss (Hz^2)
const double volt = 1/(qe*hbar_SI); // one volt (Hz)

// NV system constants
const double NV_ZFS = 2*pi*2.87e9; // NV center zero field splitting energy (Hz)
const double a0 = 0.35668 * nm; // diamond lattice parameter at 300 K

