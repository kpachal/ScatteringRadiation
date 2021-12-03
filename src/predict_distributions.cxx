#include <iostream>
#include <string>

#include "ScatteringRadiation/Functions.h"

int main(int argc, char* argv[]) {

   std::string outputfilename = "theory_curves.root";

   // Electron beam of 31 MeV on 10 micron foil
   MoliereDCS func_MDCS(10.,31.,0.510);

   // Test Moliere
   double theta = 0.1;
   func_MDCS(theta);

   // Make TF1 of Moliere distribution

}

