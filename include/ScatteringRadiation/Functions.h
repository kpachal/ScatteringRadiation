#ifndef __FUNCTIONS__
#define __FUNCTIONS__

class TSpline3; 

class MoliereDCS {

 public:
   MoliereDCS(double foilThickness, double energy, double mass);
   ~MoliereDCS() {};

   double operator() (double *x, double *p);

 private :

   double m_foilThickness;
   double m_gamma;
   double m_p_beta_c;
   double m_AM;
   double m_chiC_2_const;

   TSpline3 * m_lookup_f0;
   TSpline3 * m_lookup_f1;
   TSpline3 * m_lookup_f2;

   double solve_B(double b);

};

#endif