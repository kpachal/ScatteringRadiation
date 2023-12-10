#include <iostream>
#include <string>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TVector.h>

#include "ScatteringRadiation/Functions.h"

double nNeutrons(double thickness, double energy) {

    double radiation_length_Ta = 0.4094*10000.; // Initial value in cm; converting to microns
    double T = thickness/radiation_length_Ta; // target thickness in radiation lengths
    std::cout << "\nT is " << T << std::endl;
    double E0 = energy; // Beam energy in MeV
    double Z = 23;

    double Y = 8.0e-4 * (1 + 0.12*Z - 0.001*pow(Z,2)) * pow(T,2)/E0 * (1 + 0.04/T);
    return Y;

}

int main(int argc, char* argv[]) {

   std::string outputfilename = "theory_curves.root";

   std::vector<int> thicknesses = {1};
   std::vector<double> energies = {5.,10.,15.,20.,25.,30.,31.,32.};
   std::vector<TF1*> functions;
   std::vector<double> integrals;
   std::vector<std::pair<double,double> > est_widths;
   for (auto thickness : thicknesses) {
        for (auto energy : energies) {
            for (bool isDegrees : {true, false}) {

                // Moliere multiple scattering cross sections
                // Electron beam of specified energy in MeV on foil of given thickness
                MoliereDCS func_MDCS(double(thickness),energy,0.510);
                func_MDCS.setUseDegrees(isDegrees);

                // Make TF1 of Moliere distribution
                // Only going up to 1 radian/90 degrees - seems weird at larger angles
                double limit_low = isDegrees ? 0.5 : 0.001;
                double limit_high = isDegrees ? 45 : 1.0;
                TString funcname = Form("MDCS_%imicron_%iMeV",thickness,energy);
                if (isDegrees) funcname += "_deg";
                TF1 * MDCS = new TF1(funcname,func_MDCS,limit_low,limit_high,1); 
                // Display purposes
                MDCS->SetNpx(300);

                // Now to get good normalisation, need to retrieve what this should be compared to
                // and fit normalisation parameter.
                TString infilename = Form("results_%imicron_1e7events_%iMeV.root",thickness,int(energy));
                std::cout << "Opening input file " << infilename << std::endl;
                TFile infile(infilename,"READ");
                TString inhistname = isDegrees ? "angle_polar_deg_eminus_norm" : "angle_polar_eminus_norm";
                TH1D * inhist = (TH1D*) infile.Get(inhistname);
                inhist->SetDirectory(0);

                // Set initial normalising parameter to highest bin value
                MDCS->SetParameter(0,inhist->GetMaximum());

                // And fit it
                inhist->Fit(MDCS,"","",limit_low,limit_high);
                
                // And save the function
                functions.push_back(MDCS);

                // Don't continue if degrees - want width estimates in rads
                if (isDegrees) continue;

                // Integrate while it is properly defined
                double int_MDCS = MDCS->Integral(0.0,limit_high);
                integrals.push_back(int_MDCS);

                // Approximation of width: taken from PDG, based on Gaussian core
                double width = func_MDCS.estimate_width();

                // Width estimated as half-max point for Gaussian core.
                // Peak is the leftmost value.
                double testpoint = 0.001;
                double normpar = MDCS->GetParameter(0);
                double maxval = func_MDCS(&testpoint,&normpar);
                double width_obs = MDCS->GetX(maxval/2.0, 0, 0.1);

                est_widths.push_back(std::make_pair(width,width_obs));

            }

            // Now neutron numbers.
            double nNeut = nNeutrons(thickness,energy);
            std::cout << "For foil thickness " << thickness << " microns, Y = " << nNeut << " n per MeV per e" << std::endl;
            std::cout << "For " << energy << " MeV beam energy and 5e10 electrons, that is " << nNeut * 5e10 * energy << " neutrons." << std::endl;

        }
   }

   // Format TVectors
   TVectorD integrals_MDCS(3);
   for (unsigned int index=0; index<integrals.size(); index++) integrals_MDCS[index] = integrals.at(index);
   TVectorD widths_est(3);
   TVectorD widths_fromTF1(3);
   for (unsigned int index=0; index<est_widths.size(); index++) {
       std::cout << "Width " << thicknesses.at(index) << ":" << std::endl;
       std::cout << "From TDR forumla: " << est_widths.at(index).first << std::endl;
       std::cout << "From TF1 itself: " << est_widths.at(index).second << std::endl;
       widths_est[index] = est_widths.at(index).first;
       widths_fromTF1[index] = est_widths.at(index).second;
   }

   // Save it.
   TFile outputfile(outputfilename.c_str(),"RECREATE");

   outputfile.cd();
   // Functions
   for (auto thisfunc : functions) thisfunc->Write();
   // Integrals
   integrals_MDCS.Write("integrals_MDCS");
   widths_est.Write("Widths_fromTDRFormula");
   widths_fromTF1.Write("Widths_fromTF1");
   outputfile.Close();

}

