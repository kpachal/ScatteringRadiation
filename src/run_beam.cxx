// Goal: hit electron beam on a target and model all outgoing particles
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

// My includes
#include "ScatteringRadiation/Construction.h"
#include "ScatteringRadiation/Beam.h"
#include "ScatteringRadiation/Action.h"


// GEANT includes
#include <FTFP_BERT.hh>
#include <QBBC.hh>
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

std::map<std::string,int> particle_types = { {"eplus",-11} ,{"eminus", 11}, {"gamma", 22}, {"neutron", 2112} };

int main(int argc, char* argv[])
{

  // Parse run parameters.
  long nEvents = 10000;
  double target_thickness = 10; // thickness of foil target or diameter of wire in microns
  // Enum defined in Construction.h
  target_options target_type = foil; // use foil target as default
  std::string output_simple = "output";
  bool do_save_only = false;
  bool do_visuals = false;
  std::string save_only = "";
  int ip=1;
  int seed = 0;

  while (ip<argc) {

    if (std::string(argv[ip]).substr(0,2)=="--") {

        // Number of events
        if (std::string(argv[ip])=="--nEvents") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            nEvents = std::stol(argv[ip+1]);
            ip+=2;
          } else {std::cout<<"\nNo number of events specified."<<std::endl; return 1;}
        }

        // Foil thickness
        else if (std::string(argv[ip])=="--targetThickness") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            target_thickness = std::stod(argv[ip+1]);
            ip+=2;
          } else {std::cout<<"\nNo target thickness inserted."<<std::endl; return 1;}
        }

        // Diameter of a wire target
        else if (std::string(argv[ip])=="--wireTarget") {
            target_type = wire;
            ip+=1;
        }

        // Turn on visuals?
        // Only make this possible for small numbers of events.
        else if (std::string(argv[ip])=="--visuals") {
            if (nEvents <= 1000) do_visuals = true;
            else std::cout << "Don't request visuals with this many events!" << std::endl;
            ip+=1;
        }

        // Output file name
        else if (std::string(argv[ip])=="--output") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            output_simple = argv[ip+1];
            ip+=2;
          } else {std::cout<<"\nNo output name inserted."<<std::endl; return 1;}
        }

        // Keep just one particle type
        // Options: eplus, eminus, gamma, neutron
        else if (std::string(argv[ip])=="--saveOnly") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            do_save_only = true;
            save_only = argv[ip+1];
            ip+=2;
          } else {std::cout<<"\nNo particle type inserted."<<std::endl; return 1;}
        }

        // Random number seed for parallelisation
        else if (std::string(argv[ip])=="--seed") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            seed = std::stoi(argv[ip+1]);
            ip+=2;
          } else {std::cout<<"\nNo seed inserted."<<std::endl; return 1;}
        }

    } else { //if command does not start with "--"
        std::cout << "\nCommand '"<<std::string(argv[ip])<<"' unknown"<<std::endl;
        break;
    }//end if "--"

  }//end while loop

  // Format output name
  std::stringstream stream;
  stream << output_simple << "_" << std::setprecision(2) << target_thickness << "micron_1e" << log10(nEvents) << "events";
  //if (target_type != foil) stream << "_" << (target_type == wire ? "wire" : "gas");
  if (do_save_only) stream << "_" << save_only;
  if (seed) stream << "_seed" << seed;
  stream << ".root";
  std::string output_filename = stream.str();
  std::cout << "Creating output file " << output_filename << std::endl;

  // Set seed before initialising runmanager
  if (seed) {
    CLHEP::HepRandom::setTheSeed(seed); G4Random::setTheSeed(seed);
  }

	//Get instance of runmanager
	G4RunManager * runManager = new G4RunManager;

	// Detector construction
	runManager->SetUserInitialization(new Construction(target_thickness, target_type));

  // Physics list
  G4VModularPhysicsList* physicsList = new QBBC;
  runManager->SetUserInitialization(physicsList);

  // This includes setting up the beamline
  // Send false and a blank number if we keep everything (default)
  int PDG_only = 0;
  if (do_save_only) PDG_only = particle_types[save_only];
  runManager->SetUserInitialization(new Action(output_filename,do_save_only,PDG_only)); 

  // Initialize G4 kernel
	runManager->Initialize();

  std::cout << "Will save outputs to " << output_filename << std::endl;

  // Setting up the user interface to be nice.
  if (do_visuals) {
    G4UIExecutive * ui = new G4UIExecutive(argc, argv);
    G4VisManager * visManager = new G4VisExecutive();
    visManager->Initialize();
    G4UImanager * uiManager = G4UImanager::GetUIpointer();
    
    // Running commands to display what I'm doing
    uiManager->ApplyCommand("/run/initialize");
    uiManager->ApplyCommand("/vis/open OGL");
    //uiManager->ApplyCommand("/vis/open HepRepFile");
    uiManager->ApplyCommand("/vis/viewer/set/viewpointvector 1 1 1");
    uiManager->ApplyCommand("/vis/drawVolume");
    uiManager->ApplyCommand("/vis/viewer/set/autorefresh true");
    uiManager->ApplyCommand("/vis/scene/add/trajectories smooth");
    uiManager->ApplyCommand("/vis/scene/endOfEventAction accumulate");

    ui->SessionStart();
  }

  //Cause the run manager to generate events using the
	//primary generator action registered above.
  runManager->BeamOn(nEvents);


	//After the run is complete, free up the memory used by run 
	//manager and return 
	delete runManager;

  return 0;
}
