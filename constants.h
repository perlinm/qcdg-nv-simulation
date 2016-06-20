#pragma once

#include <math.h>
#include <complex.h>

const std::complex<double> j(0.0,1.0);
const double pi = M_PI;

// physical constants in SI units
const double mu0_SI = 4*pi*1e-7; // magnetic constant (tesla meters / amp)
const double c_SI = 299792458; // speed of light (meters / second)
const double hbar_SI = 6.582119514e-16; // reduced Planck constant (eV / second)
const double q_e_SI = 1.60217646e-19; // charge of electron (coulombs)
const double g_e_SI = -1.760859708e11; // gyromagnetic ratio of NV electron (Hz/tesla)
const double g_C13_SI = 67.28284e6; // gyromagnetic ratio of C13 (Hz/tesla)
const double g_N15_SI = -27.116e6; // gyromagnetic ratio of N15 (Hz/tesla)

// physical constants in natural units: e_0 = mu_0 = c = hbar = 1
const double alpha = 1/137.035999074; // fine structure constant
const double q_e = sqrt(4*pi*alpha); // unit electric charge
const double m_e = 510998.928/hbar_SI; // mass of electron (Hz)
const double g_e = -q_e/(2*m_e) * 2.0023193043617; // gyromagnetic ratio of electron (s)
const double g_C13 = g_C13_SI * g_e/g_e_SI; // gyromagnetic ratio of C-13 (s)
const double g_N15 = g_N15_SI * g_e/g_e_SI; // gyromagnetic ratio of N-15 (s)

// SI values in natural units
const double kHz = 1000; // one kilohertz (Hz)
const double meter = 1 / c_SI; // one meter (s)
const double nm = 1e-9*meter; // one nanometer (s)
const double coulomb = q_e / q_e_SI; // one coulomb
const double tesla = g_e_SI / g_e; // one tesla (Hz^2)
const double gauss = 1e-4*tesla; // one gauss (Hz^2)
const double volt = 1/(q_e*hbar_SI); // one volt (Hz)

// NV system constants
const double NV_ZFS = 2*pi*2.87e9; // NV center zero field splitting energy (Hz)
const double a0 = 0.35668 * nm; // diamond lattice parameter at 300 K

