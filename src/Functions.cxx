#include <math.h>
#include <iostream>

#include "ScatteringRadiation/Functions.h"

MoliereDCS::MoliereDCS(double foilThickness, double energy, double mass) {
    m_foilThickness = foilThickness;
    m_energy = energy;
    // Convert momentum and mass to energy:
    m_momentum = sqrt( pow(energy,2) - pow(mass,2)); // all in MeV

}

// Validated
double MoliereDCS::solve_B(double b) {

  // Use Newton Raphson method.
  // Want 3 decimal place accuracy
  double threshold = 0.001;
  double x0 = b;
  double xnew = x0; double diff = 10;
  while (abs(diff) > threshold) {
      x0 = xnew;
      double fun_x = x0 - log(x0) - b;  // f(x0)
      double fun_xx = 1 - 1/x0;  // f'(x0)
      xnew = x0 - (fun_x/fun_xx);
      diff = xnew - x0;
  }
  return xnew;

}


// Using 2.53, 2.59, 2.6 from JMFV
double MoliereDCS::operator() (double &theta) {

    // Define individual input quantities
    double density_Ta = 16.69; // units: g/cm3
    double Aweight_Ta = 180.94788; // units: grams per mole
    double NA = 6.02214076e23; // Avogadro's number, units: atoms per mole
    double cm3_to_fm3 = 1e39;
    double micron_to_fm = 1e9;
    double s = m_foilThickness/cos(theta); // units: microns
    double N = (NA * density_Ta)/(cm3_to_fm3 * Aweight_Ta); // units: 1/fm^3
    int Z = 73; // atomic number of tantalum
    int Zprime = -1; // electrons have charge -e
    double esquared = 1.4399764; // units MeV * fm
    double beta = m_momentum/m_energy; // Standard relativistic beta, pc/E
    double gamma = 0.5772; // Euler's constant
    std::cout << "beta: " << beta << std::endl;
    std::cout << "N: " << N << std::endl;
    
    // Mid level calculations
    double p_beta_c = m_momentum*m_momentum/m_energy; // Units are MeV. Validated.
    double AM = pow(0.510,2)/(4*pow(137,2)*pow(m_momentum,2)) 
                * pow((0.885*pow(Z,-1/3)),-2)
                * (1.13 + 3.76 * pow(Z/(beta*137),2)); // units: none
    std::cout << "AM: " << AM << std::endl;

    // Actual input variables
    std::cout << "p_beta_c = " << p_beta_c << std::endl;
    double chiC_2 = s * micron_to_fm * N * 4 * M_PI * pow(Z*Zprime*esquared,2)/pow(p_beta_c,2); // units: fm * 1/fm^3 * (MeV*fm)^2/(MeV^2)
    double b = log(chiC_2/AM) + 1 - 2*gamma;
    double B = solve_B(b);
    double v = theta/sqrt(chiC_2*B);

    // Cross checks ...
    std::cout << "b is: " << b << std::endl;
    std::cout << "B is: " << B << std::endl;
    std::cout << "chiC^2 is: " << chiC_2 << std::endl;
    std::cout << "v is: " << v << std::endl;

    // Our three equations (approximate)
    double f0 = 2*exp(- v*v);
    double f1 = 2/pow(v,4) * pow(1 - pow())

}
