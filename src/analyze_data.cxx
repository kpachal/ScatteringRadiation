#include <iostream>
#include <string>
#include <list>
#include <algorithm>
#include <math.h>

#include <TH1I.h>
#include <TROOT.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TLorentzVector.h>
#include <ROOT/RDataFrame.hxx>
#include <TGraph.h>
#include <TCanvas.h>

// List of known processes (100M events)
// 1091 2002 2003 2005 2010 2011 2012 2013 2014 4121
std::vector<int> process_types = {1091, 2002, 2003, 2005, 2010, 2011, 2012, 2013, 2014, 4121};

// List of output particles (100M events)
// -11 11 22 2112 1000731790 1000731800 1000731810
// Last 3 are tantalum isotopes
// https://twiki.cern.ch/twiki/pub/Geant4/ExtendingFnalDb/Nuclei.txt
std::map<std::string,int> particle_types = { {"eplus",-11} ,{"eminus", 11}, {"gamma", 22}, {"neutron", 2112} };

using namespace ROOT; // RDataFrame's namespace

std::vector<TLorentzVector> build_particles(int desiredPDGID, std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> mass, std::vector<int> pdgIDs) {
    std::vector<TLorentzVector> particles;
    for (int i=0; i<px.size(); i++) {
        if (pdgIDs.at(i) != desiredPDGID) continue;
        TLorentzVector p4;
        p4.SetXYZM(px.at(i),py.at(i),pz.at(i),mass.at(i));
        particles.push_back(p4);
    }
    return particles;
}

int main(int argc, char* argv[]) {

    std::string input_filename = "output_byevent.root";
    std::string output_filename = "results.root";
    std::string tree_name = "Events";
    int ip=1;
    while (ip<argc) {

      if (std::string(argv[ip]).substr(0,2)=="--") {

          // Input file
          if (std::string(argv[ip])=="--input") {
            if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
              input_filename = argv[ip+1];
              ip+=2;
            } else {std::cout<<"\nNo input file name inserted."<<std::endl; break;}
          }

          // Tree name
          else if (std::string(argv[ip])=="--tree") {
            if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
              tree_name = argv[ip+1];
              ip+=2;
            } else {std::cout<<"\nNo tree name inserted"<<std::endl; break;}
          }
      } else { //if command does not start with "--"
          std::cout << "\nCommand '"<<std::string(argv[ip])<<"' unknown"<<std::endl;
          break;
      }//end if "--"

    }//end while loop

    // Just stop if don't have needed info
    if (input_filename.empty()) {
      std::cout << "No file name to run on!" << std::endl;
      exit(1);
    } else if (tree_name.empty()) {
      std::cout << "You need to specify a tree to use." << std::endl;
      exit(1);
    }

    std::cout << "Running over file: " << input_filename << std::endl;
    std::vector<std::string> files_to_use = {input_filename};

    // Make the RDataFrame!
    TFile *fin = TFile::Open(input_filename.c_str(), "READ");
    TTree *tree = (TTree*)(fin -> Get(tree_name.c_str()));
    RDataFrame frame(*tree);

    // Add columns
    auto full_types = frame.Define("fullType",[](std::vector<int> processType, std::vector<int> processSubType)
        {   std::vector<int> fullTypes;
            for (int i=0; i<processType.size(); i++) fullTypes.push_back(processType.at(i)*1000+processSubType.at(i));
            return fullTypes;
        }, {"processType","processSubType"} );
    auto n_interesting = full_types.Define("nInteresting", [](std::vector<int> fullType)
        {   int nInteresting = 0;
            for (auto thistype : fullType) if (thistype !=1091) nInteresting++;
            return nInteresting;
        }, {"fullType"});
    auto interesting = n_interesting.Define("interestingType",[](std::vector<int> fullType, int nInteresting)
        {   std::vector<int> interestingTypes;
            for (auto thistype : fullType) if (thistype !=1091) interestingTypes.push_back(thistype);
            return interestingTypes;
        }, {"fullType","nInteresting"} );
    auto nHad = interesting.Define("nHadronic",[](std::vector<int> interestingType)
        {   int nHadronicProcesses = 0;
            for (auto thistype : interestingType) if (thistype > 3999 && thistype < 5000) nHadronicProcesses++;
            return nHadronicProcesses;
        }, {"interestingType"});
    auto with_dominant = nHad.Define("dominantType",[](std::vector<int> fullType)
        {
            // Here we dictate which process to count an event as being in the case that
            // there are a bunch of different processes that occur.

            // Define order of importance: 1 = most important
            int dominantType = 1091;
            for (auto type : fullType) {
                // Decreasing order of interestingness
                if (type > 4999) {
                    dominantType = type;
                    break;
                } else if (type > 3999) {
                    dominantType = type;
                    break;
                } else if (type == 2014) {
                    dominantType = type;
                    break;
                } else if (type == 2005) {
                    dominantType = type;
                    break;
                } else if (type == 2012) {
                    dominantType = type;
                    break;
                } else if (type == 2004) {
                    dominantType = type;
                    break;
                } else if (type == 2003) {
                    dominantType = type;
                    break;
                } else if (type == 2002) {
                    dominantType = type;
                    break;
                } else if (type == 2001) {
                    dominantType = type;
                    break;
                } else if (type == 2011) {
                    dominantType = type;
                    break;
                } else if (type == 2010) {
                    dominantType = type;
                    break;
                } else if (type == 2013) {
                    dominantType = type;
                    break;
                } else if (type != 1091) {
                    dominantType = type;
                    break;
                }
            }
            return dominantType;
        }
        , {"fullType"});

    // Make final state particles and sort by type.
    // Having 3-vectors will allow us to determine angles
    int PDGid_eplus = particle_types["eplus"];
    auto withEplus = with_dominant.Define("eplus",[PDGid_eplus](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_eplus, px, py, pz, mass, pdgIDs); },
    {"momX", "momY", "momZ", "mass", "pdgIDFinal"});
    int PDGid_eminus = particle_types["eminus"];
    auto withEminus = withEplus.Define("eminus",[PDGid_eminus](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_eminus, px, py, pz, mass, pdgIDs); },
    {"momX", "momY", "momZ", "mass", "pdgIDFinal"});
    int PDGid_gamma = particle_types["gamma"];
     auto withGamma = withEminus.Define("gamma",[PDGid_gamma](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_gamma, px, py, pz, mass, pdgIDs); },
    {"momX", "momY", "momZ", "mass", "pdgIDFinal"});
    int PDGid_neutron = particle_types["neutron"];
    auto final_frame = withGamma.Define("neutron",[PDGid_neutron](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_neutron, px, py, pz, mass, pdgIDs); },
    {"momX", "momY", "momZ", "mass", "pdgIDFinal"});

    // Assess full set of processes present
    auto setOfAllProcesses = final_frame.Reduce([](std::vector<int> processes_c1, std::vector<int> processes_c2) 
        {   std::unordered_set <int> uset;
            for (auto t1 : processes_c1) uset.insert(t1);
            for (auto t2 : processes_c2) uset.insert(t2);
            std::vector<int> cleanvec(uset.begin(), uset.end());
            std::sort(cleanvec.begin(), cleanvec.end());
            return cleanvec;
         }, "fullType");

    // Assess full set of output particles
    auto setOfAllParticles = final_frame.Reduce([](std::vector<int> particles_c1, std::vector<int> particles_c2) 
        {   std::unordered_set <int> uset;
            for (auto t1 : particles_c1) uset.insert(t1);
            for (auto t2 : particles_c2) uset.insert(t2);
            std::vector<int> cleanvec(uset.begin(), uset.end());
            std::sort(cleanvec.begin(), cleanvec.end());
            return cleanvec;
         }, "pdgIDFinal");

    // Save histograms to evaluate at the end
    std::vector<ROOT::RDF::RResultPtr<TH1D> > outputs;

    // Make general histograms
    auto hist_fullTypes = final_frame.Histo1D("fullType"); outputs.push_back(hist_fullTypes);
    auto hist_nInteresting = final_frame.Histo1D("nInteresting"); outputs.push_back(hist_nInteresting);
    auto hist_interestingTypes = final_frame.Histo1D("interestingType"); outputs.push_back(hist_interestingTypes);
    auto hist_dominantTypes = final_frame.Histo1D("dominantType"); outputs.push_back(hist_dominantTypes);

    // Make separate plots by final state particle
    for (const auto& [name, pdgID] : particle_types) {

        // Just get events with these particles, in case I want to count them at some point
        auto frame_theseparticles = final_frame.Filter([](std::vector<TLorentzVector> particles)
       {return particles.size()>0; },{name});

        auto frame_withQuantities = frame_theseparticles
        .Define(("energy_" + name).c_str(),
            [](std::vector<TLorentzVector> particles)
            {   std::vector<double> energy;
                for (auto particle : particles) energy.push_back(particle.E());
                return energy; },{name})
        .Define(("momentum_" + name).c_str(),
            [](std::vector<TLorentzVector> particles)
            {   std::vector<double> momentum;
                for (auto particle : particles) momentum.push_back(particle.P());
                return momentum;}, {name})
        .Define(("angle_" + name).c_str(),
            [](std::vector<TLorentzVector> particles)
            {   std::vector<double> angles;
                // Unit vector we want direction relative to
                TVector3 beamDir(0,0,1);
                for (auto particle : particles) {
                    auto particleDir = particle.Vect();
                    double angle = particleDir.Angle(beamDir);
                    angles.push_back(angle);
                }
                return angles;}, {name});

        // Now make histograms from each.
        auto particle_energy = frame_withQuantities.Histo1D({("energy_" + name).c_str(),("energy_" + name).c_str(),70,0,35},("energy_" + name).c_str());
        outputs.push_back(particle_energy);


    }

    // Make separate plots by process

    // Now: evaluate.
    std::cout << "Full set of processes:" <<  std::endl;
    for (auto const &i: setOfAllProcesses) {
        std::cout << i << " ";
    }
    std::cout << std::endl;

    std::cout << "Full set of final state particle pdgIDs:" <<  std::endl;
    for (auto const &i: setOfAllParticles) {
        std::cout << i << " ";
    }
    std::cout << std::endl;

    // Make output file.
    TFile * output_file = new TFile(output_filename.c_str(),"RECREATE");
    output_file->cd();
    for (auto hist : outputs) hist.GetValue().Write();
    output_file->Close();

    // Make tables of data
    //listRepresented
}