#include <iostream>
#include <string>
#include <TF1.h>
#include <TFile.h>
#include <TVector.h>

#include "ScatteringRadiation/Functions.h"

int main(int argc, char* argv[]) {

   std::string outputfilename = "theory_curves.root";

   // Electron beam of 31 MeV on 1 micron foil
   MoliereDCS func_MDCS_1(1.,31.,0.510);
   // Electron beam of 31 MeV on 5 micron foil
   MoliereDCS func_MDCS_5(5.,31.,0.510);
   // Electron beam of 31 MeV on 10 micron foil
   MoliereDCS func_MDCS_10(10.,31.,0.510);

   // Make TF1 of Moliere distribution
   // Only going up to 1 radian - seems weird at larger angles
   TF1 * MDCS_1 = new TF1("MDCS_1micron",func_MDCS_1,0.001,1.0,1); 
   TF1 * MDCS_5 = new TF1("MDCS_5micron",func_MDCS_5,0.001,1.0,1); 
   TF1 * MDCS_10 = new TF1("MDCS_10micron",func_MDCS_10,0.001,1.0,1); 
   MDCS_1->SetNpx(300);    MDCS_5->SetNpx(300);    MDCS_10->SetNpx(300);

   // Integrate while they are properly defined
   //double int_MDCS_1 = MDCS_1->Integral(0.0,1.0);
   /*Int_t np = 3000;
   double * x_1=new double[np];
   double * w_1=new double[np];
   MDCS_1->CalcGaussLegendreSamplingPoints(np,x_1,w_1,1e-15); 
   double testintegral_1 = MDCS_1->IntegralFast(np,x_1,w_1,0.0,1.0); */
   double int_MDCS_1 = MDCS_1->Integral(0.0,1.0);
   double int_MDCS_5 = MDCS_5->Integral(0.0,1.0);
   double int_MDCS_10 = MDCS_10->Integral(0.0,1.0);
   TVectorD integrals_MDCS(3);
   integrals_MDCS[0] = int_MDCS_1;
   integrals_MDCS[1] = int_MDCS_5;
   integrals_MDCS[2] = int_MDCS_10;
   std::cout << "Integrals are " << int_MDCS_1 << " " << int_MDCS_5 << " " << int_MDCS_10 << std::endl;
   //std::cout << "Alt integrals are " << testintegral_1 << std::endl;

   // Save it.
   TFile outputfile(outputfilename.c_str(),"RECREATE");
   outputfile.cd();
   // Functions
   MDCS_1->Write();
   MDCS_5->Write();
   MDCS_10->Write();
   // Integrals
   integrals_MDCS.Write("integrals_MDCS");
   outputfile.Close();

}

