#ifndef __FUNCTIONS__
#define __FUNCTIONS__

class MoliereDCS {

 public:
   MoliereDCS(double foilThickness, double energy, double mass);
   ~MoliereDCS() {};

   double operator() (double &theta);

 private :
   double m_foilThickness;
   double m_momentum;
   double m_energy;

   double solve_B(double b);

};

#endif