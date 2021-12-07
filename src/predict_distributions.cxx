#include <iostream>
#include <string>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TVector.h>

#include "ScatteringRadiation/Functions.h"

int main(int argc, char* argv[]) {

   std::string outputfilename = "theory_curves.root";

   std::vector<int> thicknesses = {1,5,10};
   std::vector<TF1*> functions;
   std::vector<double> integrals;
   std::vector<double> est_widths;
   for (auto thickness : thicknesses) {

    // Moliere multiple scattering cross sections
    // Electron beam of 31 MeV on foil of given thickness
    MoliereDCS func_MDCS(double(thickness),31.,0.510);

    // Make TF1 of Moliere distribution
    // Only going up to 1 radian - seems weird at larger angles
    TF1 * MDCS = new TF1(Form("MDCS_%imicron",thickness),func_MDCS,0.001,1.0,1); 
    // Display purposes
    MDCS->SetNpx(300);

    // Now to get good normalisation, need to retrieve what this should be compared to
    // and fit normalisation parameter.
    TString infilename = Form("results_%imicron_1e8events.root",thickness);
    std::cout << "Opening input file " << infilename << std::endl;
    TFile infile(infilename,"READ");
    TH1D * inhist = (TH1D*) infile.Get("angle_polar_eminus_norm");
    inhist->SetDirectory(0);

    // Set initial normalising parameter to highest bin value
    std::cout << "Setting initial normalisation to " << inhist->GetMaximum() << std::endl;
    MDCS->SetParameter(0,inhist->GetMaximum());

    // And fit it
    inhist->Fit(MDCS,"","",0.01,1.0);
    
    // Get parameter value
    std::cout << "Fitted parameter value is " << MDCS->GetParameter(0) << std::endl;

    // Integrate while it is properly defined
    double int_MDCS = MDCS->Integral(0.0,1.0);
    integrals.push_back(int_MDCS);

    // And save the function
    functions.push_back(MDCS);

    // Approximation of width: taken from PDG, based on Gaussian core
    double width = func_MDCS.estimate_width();
    std::cout << "Width estimated from TDR formula: " << width << " rads" << std::endl;

    // Width estimated as half-max point for Gaussian core.
    // Peak is the leftmost value.
    double testpoint = 0.001;
    double normpar = MDCS->GetParameter(0);
    double maxval = func_MDCS(&testpoint,&normpar);
    double width_obs = MDCS->GetX(maxval/2.0, 0, 0.1);
    std::cout << "Half-width estimated from TF1: " << width_obs << " rads" << std::endl;

   }

   // Save it.
   TFile outputfile(outputfilename.c_str(),"RECREATE");
   TVectorD integrals_MDCS(3);
   for (unsigned int index=0; index<integrals.size(); index++) integrals_MDCS[index] = integrals.at(index);
   outputfile.cd();
   // Functions
   for (auto thisfunc : functions) thisfunc->Write();
   // Integrals
   integrals_MDCS.Write("integrals_MDCS");
   outputfile.Close();

}

