#include <iostream>
#include <string>
#include <list>
#include <algorithm>
#include <math.h>
#include <regex>

#include <TH1I.h>
#include <TROOT.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TLorentzVector.h>
#include <TVector2.h>
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

struct particle {
    TLorentzVector four_vector;
    TVector3 position_vector;
    int pdgID;
};

// Using 2.53, 2.59, 2.6 from JMFV
double MoliereDCS(double theta, double foilThickness, double momentum, double energy) {

    // Define individual input quantities
    double density_Ta = 16.69; // units: g/cm3
    double Aweight_Ta = 180.94788; // units: none
    double daltons_g = 1.66053906660e-24; // units: g, converted from 1e-27kg
    double NA = 6.02214076e23; // Avogadro's number, units: none
    double cm3_to_fm3 = 1e39;
    double micron_to_fm = 1e9;

    // approximate path length from foil thickness and angle we end up at
    double s = foilThickness/cos(theta); // units: microns
    double N = (NA * density_Ta)/(cm3_to_fm3 * Aweight_Ta * daltons_g); // units: 1/fm^3
    int Z = 73; // atomic number of tantalum
    int Zprime = -1; // electrons have charge -e
    double esquared = 1.4399764; // units MeV * fm

    // Mid level calculations
    double p_beta_c = momentum*momentum/energy; // beta = pc/E; p is in MeV/c and E is in MeV. Units are MeV

    // Actual input variables
    double chiC_2 = s * micron_to_fm * N * 4 * M_PI * pow(Z*Zprime*esquared,2)/pow(p_beta_c,2); // units: none
    //double b = ;
    //double B = ; // solve from b somehow?

}

std::vector<particle> build_particles(int desiredPDGID, std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> mass, std::vector<int> pdgIDs) {
    std::vector<particle> particles;
    for (int i=0; i<px.size(); i++) {
        // Skip if wrong pdgID
        if (pdgIDs.at(i) != desiredPDGID) continue;
        particle newparticle;
        TLorentzVector p4;
        p4.SetXYZM(px.at(i),py.at(i),pz.at(i),mass.at(i));
        // Skip if magnitude of momentum is identically zero - somehow these are sneaking in. 
        // TODO understand why they are here....
        if (p4.P() == 0) continue;
        newparticle.four_vector = p4;
        TVector3 pos;
        pos.SetXYZ(x.at(i),y.at(i),z.at(i));
        newparticle.position_vector = pos;
        newparticle.pdgID = pdgIDs.at(i);
        particles.push_back(newparticle);
    }
    return particles;
}

int main(int argc, char* argv[]) {

    std::string input_filename = "output.root";
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

    std::string output_filename = std::regex_replace(input_filename, std::regex("output"),"results");

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
            int dominantType = 0;
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
                } else if (type == 2003) {
                    dominantType = type;
                    break;
                } else if (type == 2002) {
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
    auto withEplus = with_dominant.Define("eplus",[PDGid_eplus](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_eplus, px, py, pz, x, y, z, mass, pdgIDs); },
    {"momX", "momY", "momZ", "posX", "posY", "posZ", "mass", "pdgIDFinal"});
    int PDGid_eminus = particle_types["eminus"];
    auto withEminus = withEplus.Define("eminus",[PDGid_eminus](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_eminus, px, py, pz, x, y, z, mass, pdgIDs); },
    {"momX", "momY", "momZ", "posX", "posY", "posZ", "mass", "pdgIDFinal"});
    int PDGid_gamma = particle_types["gamma"];
     auto withGamma = withEminus.Define("gamma",[PDGid_gamma](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_gamma, px, py, pz, x, y, z, mass, pdgIDs); },
    {"momX", "momY", "momZ", "posX", "posY", "posZ", "mass", "pdgIDFinal"});
    int PDGid_neutron = particle_types["neutron"];
    auto final_frame = withGamma.Define("neutron",[PDGid_neutron](std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> mass, std::vector<int> pdgIDs)
        {   return build_particles(PDGid_neutron, px, py, pz, x, y, z, mass, pdgIDs); },
    {"momX", "momY", "momZ", "posX", "posY", "posZ", "mass", "pdgIDFinal"});

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
    std::vector<ROOT::RDF::RResultPtr<TH2D> > outputs_2D;

    // Make general histograms
    auto hist_fullTypes = final_frame.Histo1D("fullType"); outputs.push_back(hist_fullTypes);
    auto hist_nInteresting = final_frame.Histo1D("nInteresting"); outputs.push_back(hist_nInteresting);
    auto hist_interestingTypes = final_frame.Histo1D("interestingType"); outputs.push_back(hist_interestingTypes);
    auto hist_dominantTypes = final_frame.Histo1D("dominantType"); outputs.push_back(hist_dominantTypes);

    // Make separate plots by final state particle
    for (const auto& [name, pdgID] : particle_types) {

        // Just get events with these particles, in case I want to count them at some point
        //auto frame_theseparticles = final_frame.Filter([](std::vector<particle> particles)
       //{return particles.size()>0; },{name});

        auto frame_withQuantities = final_frame
        .Define(("energy_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> energy;
                for (auto p : particles) energy.push_back(p.four_vector.E());
                return energy; },{name})
        .Define(("momentum_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> momentum;
                for (auto p : particles) momentum.push_back(p.four_vector.P());
                return momentum;}, {name})
        // Angle of momentum away from initial beam direction
        .Define(("angle_polar_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> angles;
                // Unit vector we want direction relative to
                TVector3 beamDir(0,0,1);
                for (auto p : particles) {
                    auto particleDir = p.four_vector.Vect();
                    double angle = particleDir.Angle(beamDir);
                    angles.push_back(angle);
                }
                return angles;}, {name})
        // Individual momentum components
        .Define(("momentum_x_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> xpos;
                for (auto p : particles) xpos.push_back(p.four_vector.Px());
                return xpos;}, {name})
        .Define(("momentum_y_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> ypos;
                for (auto p : particles) ypos.push_back(p.four_vector.Py());
                return ypos;}, {name})
        .Define(("momentum_z_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> zpos;
                for (auto p : particles) zpos.push_back(p.four_vector.Pz());
                return zpos;}, {name})
        // Individual position components
        .Define(("position_x_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> xpos;
                for (auto p : particles) xpos.push_back(p.position_vector.X());
                return xpos;}, {name})
        .Define(("position_y_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> ypos;
                for (auto p : particles) ypos.push_back(p.position_vector.Y());
                return ypos;}, {name})
        .Define(("position_z_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> zpos;
                for (auto p : particles) zpos.push_back(p.position_vector.Z());
                return zpos;}, {name})
        // Directional angles to check gaussianity: x and y separately
        .Define(("angle_x_"+name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> angles;
                // Unit vector we want direction relative to
                TVector2 beamCenter(0,1);
                for (auto p : particles) {
                    // Project location into only x and z dimensions
                    TVector2 particleLocation(p.position_vector.X(),p.position_vector.Z());                    
                    double angle = beamCenter.DeltaPhi(particleLocation);
                    angles.push_back(angle);
                }
                return angles;}, {name})
        .Define(("angle_y_"+name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> angles;
                // Unit vector we want direction relative to
                TVector2
                 beamCenter(0,1);
                for (auto p : particles) {
                    // Project location into only x and z dimensions
                    TVector2 particleLocation(p.position_vector.Y(),p.position_vector.Z());                    
                    double angle = beamCenter.DeltaPhi(particleLocation);
                    angles.push_back(angle);
                }
                return angles;}, {name})
        // Angle in degrees
        .Define(("angle_polar_deg_" + name).c_str(),
            [](std::vector<particle> particles)
            {   std::vector<double> angles;
                // Unit vector we want direction relative to
                TVector3 beamDir(0,0,1);
                for (auto p : particles) {
                    double angle = p.position_vector.Angle(beamDir);
                    // Convert angle to degrees
                    angles.push_back(angle*57.2958);
                }
                return angles;}, {name});

        // Now make histograms from each.
        auto particle_energy = frame_withQuantities.Histo1D({("energy_" + name).c_str(),("energy_" + name).c_str(),200,0,100},("energy_" + name).c_str());
        outputs.push_back(particle_energy);
        auto particle_momentum = frame_withQuantities.Histo1D({("momentum_" + name).c_str(),("momentum_" + name).c_str(),400,0,100},("momentum_" + name).c_str());
        outputs.push_back(particle_momentum); 
        auto particle_mom_x = frame_withQuantities.Histo1D({("momentum_x_" + name).c_str(),("momentum_x_" + name).c_str(),400,-50,50},("momentum_x_" + name).c_str());      
        outputs.push_back(particle_mom_x);
        auto particle_mom_y = frame_withQuantities.Histo1D({("momentum_y_" + name).c_str(),("momentum_y_" + name).c_str(),400,-50,50},("momentum_y_" + name).c_str());      
        outputs.push_back(particle_mom_y);
        auto particle_mom_z = frame_withQuantities.Histo1D({("momentum_z_" + name).c_str(),("momentum_z_" + name).c_str(),400,-50,50},("momentum_z_" + name).c_str());      
        outputs.push_back(particle_mom_z);
        auto momentum_xy = frame_withQuantities.Histo2D({("momentum_xy_" + name).c_str(),("momentum_xy_" + name).c_str(),1200,-50,50,1200,-50,50},("momentum_x_" + name).c_str(),("momentum_y_" + name).c_str());
        outputs_2D.push_back(momentum_xy);
        auto particle_angle = frame_withQuantities.Histo1D({("angle_polar_" + name).c_str(),("angle_polar_" + name).c_str(),314,0,3.14},("angle_polar_" + name).c_str());
        outputs.push_back(particle_angle);
        auto particle_pos_x = frame_withQuantities.Histo1D({("position_x_" + name).c_str(),("position_x_" + name).c_str(),402,-2010,2010},("position_x_" + name).c_str());      
        outputs.push_back(particle_pos_x);
        auto particle_pos_y = frame_withQuantities.Histo1D({("position_y_" + name).c_str(),("position_y_" + name).c_str(),402,-2010,2010},("position_y_" + name).c_str());      
        outputs.push_back(particle_pos_y);
        auto particle_pos_z = frame_withQuantities.Histo1D({("position_z_" + name).c_str(),("position_z_" + name).c_str(),402,-2010,2010},("position_z_" + name).c_str());      
        outputs.push_back(particle_pos_z);
        auto position_xy = frame_withQuantities.Histo2D({("position_xy_" + name).c_str(),("position_xy_" + name).c_str(),1202,-2010,2010,1202,-2010,2010},("position_x_" + name).c_str(),("position_y_" + name).c_str());
        outputs_2D.push_back(position_xy);
        auto particle_angle_polar_deg = frame_withQuantities.Histo1D({("angle_polar_deg_" + name).c_str(),("angle_polar_deg_" + name).c_str(),314,0,3.14},("angle_polar_deg_" + name).c_str());
        outputs.push_back(particle_angle_polar_deg);
        auto particle_angle_x = frame_withQuantities.Histo1D({("angle_x_" + name).c_str(),("angle_x_" + name).c_str(),628,-3.14,3.14},("angle_x_" + name).c_str());
        outputs.push_back(particle_angle_x); 
        auto particle_angle_y = frame_withQuantities.Histo1D({("angle_y_" + name).c_str(),("angle_y_" + name).c_str(),628,-3.14,3.14},("angle_y_" + name).c_str());
        outputs.push_back(particle_angle_y);                 
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
    for (auto hist : outputs_2D) hist.GetValue().Write();
    output_file->Close();

    // Make tables of data
    // RMS of scattering
    //std::cout << "RMS of electron scattering is"
}