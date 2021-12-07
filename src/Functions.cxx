#include <math.h>
#include <iostream>
#include <TGraph.h>
#include <TSpline.h>

#include "ScatteringRadiation/Functions.h"

MoliereDCS::MoliereDCS(double foilThickness, double energy, double mass) {

    m_foilThickness = foilThickness;
    // Convert momentum and mass to energy:
    double momentum = sqrt( pow(energy,2) - pow(mass,2)); // all in MeV

    // Set of constants; members when need to be accessed later
    double density_Ta = 16.69; // units: g/cm3
    double Aweight_Ta = 180.94788; // units: grams per mole
    double NA = 6.02214076e23; // Avogadro's number, units: atoms per mole
    double cm3_to_fm3 = 1e39;
    double micron_to_fm = 1e9;
    int Z = 73; // atomic number of tantalum
    int Zprime = -1; // electrons have charge -e
    double esquared = 1.4399764; // units MeV * fm    
    double beta = momentum/energy; // Standard relativistic beta, pc/E
    m_gamma = 0.5772; // Euler's constant

    // Define individual needed quantities
    double N = (NA * density_Ta)/(cm3_to_fm3 * Aweight_Ta); // units: 1/fm^3
    m_p_beta_c = momentum*momentum/energy; // Units are MeV. Validated.

    // Mid level calculations
    m_AM = pow(0.510,2)/(4*pow(137,2)*pow(momentum,2)) 
                * pow((0.885*pow(Z,-1/3)),-2)
                * (1.13 + 3.76 * pow(Z/(beta*137),2)); // units: none
    m_chiC_2_const = micron_to_fm * N * 4 * M_PI * pow(Z*Zprime*esquared,2)/pow(m_p_beta_c,2);

    // These splines are used for low values of v.
    // They use data from Bethe (Phys Rev 89 (6) 1953)
    std::vector<double> xvals {0.0, 0.2, 0.4, 0.6, 0.8, 
                               1.0, 1.2, 1.4, 1.6, 1.8, 
                               2.0, 2.2, 2.4, 2.6, 2.8,
                               3.0, 3.2, 3.4, 3.6, 3.8,
                               4.0, 4.5, 5.0, 5.5, 6.0};
    std::vector<double> y_f0 {2.0, 1.9216, 1.7214, 1.4094, 1.0546,
                              0.7338, 0.4738, 0.2817, 0.1546, 0.0783,
                              0.0366, 0.01581, 0.00630, 0.00232, 0.00079,
                              0.000250, 7.3e-5, 1.9e-5, 4.7e-6, 1.1e-6,
                              2.3e-7, 3e-9, 2e-11, 2e-13, 5e-16};
    std::vector<double> y_f1 {0.8456, 0.7038, 0.3437, -0.0777, -0.3981,
                              -0.5285, -0.4770, -0.3183, -0.1396, -0.0006,
                              0.0782, 0.1054, 0.1008, 0.08262, 0.06247,
                              0.04550, 0.03288, 0.02402, 0.01791, 0.01366,
                              10.638e-3, 6.140e-3, 3.831e-3, 2.527e-3, 1.739e-3};
    std::vector<double> y_f2 {2.4929, 2.0694, 1.0488, -0.0044, -0.6068,
                              -0.6359, -0.3086, 0.0525, 0.2423, 0.2386,
                              0.1316, 0.0196, -0.0467, -0.0649, -0.0546,
                              -0.03568, -0.01923, -0.00847, -0.00264, 0.00005, 
                              1.0741e-3, 1.2294e-3, 0.8326e-3, 0.5368e-3, 0.3495e-3};
    TGraph f0_graph(xvals.size(),&xvals[0], &y_f0[0]);
    m_lookup_f0 = new TSpline3("f0_spline",&f0_graph);
    TGraph f1_graph(xvals.size(),&xvals[0], &y_f1[0]);
    m_lookup_f1 = new TSpline3("f1_spline",&f1_graph);
    TGraph f2_graph(xvals.size(),&xvals[0], &y_f2[0]);
    m_lookup_f2 = new TSpline3("f2_spline",&f2_graph);

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
double MoliereDCS::operator() (double *x, double *par) {

    // Parsing parameters:
    // don't care about p.
    double theta = x[0];

    // Path length 
    double s = m_foilThickness/cos(theta); // units: microns

    // Actual input variables
    double chiC_2 = s * m_chiC_2_const; // units: fm * 1/fm^3 * (MeV*fm)^2/(MeV^2)
    double b = log(chiC_2/m_AM) + 1 - 2*m_gamma;
    double B = solve_B(b);
    double v = theta/sqrt(chiC_2*B);

    // Cross checks ...
    /*std::cout << "s is: " << s << std::endl;
    std::cout << "chiC^2 is: " << chiC_2 << std::endl;
    std::cout << "b is: " << b << std::endl;
    std::cout << "B is: " << B << std::endl;
    std::cout << "v is: " << v << std::endl; */

    // Our three equations (approximate)
    double f0, f1, f2;
    if (v < 6) {
        f0 = m_lookup_f0->Eval(v);
        f1 = m_lookup_f1->Eval(v);
        f2 = m_lookup_f2->Eval(v);
    } else {
        f0 = 2*exp(- v*v);
        f1 = 2/pow(v,4) * pow((1 - 5./pow(v,2)),-4./5.);
        f2 = 16./pow(v,6) * (log(v) + m_gamma - 3./2.) * (1 + 9./pow(v,2)) - 38./pow(v,8);
    }
    /*std::cout << "f0: " << f0 << std::endl;
    std::cout << "f1: " << f1 << std::endl;
    std::cout << "f2: " << f2 << std::endl; */

    // Moliere function
    double FM = 1./(2.*M_PI) * pow(theta/sin(theta),0.5) * exp(chiC_2*B/16.)/(chiC_2*B) * (f0 + f1/B + f2/pow(B,2));
    //std::cout << "Theta = " << theta << " gives MDCS = " << FM << std::endl;

    return par[0]*FM;
}

double MoliereDCS::estimate_width() {

    double radiation_length_Ta = 0.4094*10000.; // Initial value in cm; converting to microns
    double x_over_x0 = m_foilThickness/radiation_length_Ta;

    double width = 13.6/m_p_beta_c * sqrt(x_over_x0);// * (1 + 0.038*log(x_over_x0));
    
    return width;
}
