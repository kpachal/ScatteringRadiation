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
  double target_thickness = 10; // thickness of foil in microns
  std::string output_simple = "output";
  bool do_save_only = false;
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
          } else {std::cout<<"\nNo number of events specified."<<std::endl; break;}
        }

        // Foil thickness
        else if (std::string(argv[ip])=="--foilThickness") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            target_thickness = std::stod(argv[ip+1]);
            ip+=2;
          } else {std::cout<<"\nNo foil thickness inserted"<<std::endl; break;}
        }

        // Output file name
        else if (std::string(argv[ip])=="--output") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            output_simple = argv[ip+1];
            ip+=2;
          } else {std::cout<<"\nNo output name inserted"<<std::endl; break;}
        }

        // Keep just one particle type
        // Options: eplus, eminus, gamma, neutron
        else if (std::string(argv[ip])=="--saveOnly") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            do_save_only = true;
            save_only = argv[ip+1];
            ip+=2;
          } else {std::cout<<"\nNo particle type inserted"<<std::endl; break;}
        }

        // Random number seed for parallelisation
        if (std::string(argv[ip])=="--seed") {
          if (ip+1<argc && std::string(argv[ip+1]).substr(0,2)!="--") {
            seed = std::stoi(argv[ip+1]);
            ip+=2;
          } else {std::cout<<"\nNo seed inserted."<<std::endl; break;}
        }        

    } else { //if command does not start with "--"
        std::cout << "\nCommand '"<<std::string(argv[ip])<<"' unknown"<<std::endl;
        break;
    }//end if "--"

  }//end while loop

  // Format output name
  std::stringstream stream;
  stream << output_simple << "_" << std::setprecision(2) << target_thickness << "micron_1e" << log10(nEvents) << "events";
  if (do_save_only) stream << "_" << save_only;
  if (seed) stream << "_seed" << seed;
  stream << ".root";
  std::string output_filename = stream.str();

  // Set seed before initialising runmanager
  if (seed) {
    CLHEP::HepRandom::setTheSeed(seed); G4Random::setTheSeed(seed);
  }

	//Get instance of runmanager
	G4RunManager * runManager = new G4RunManager;

	// Detector construction
	runManager->SetUserInitialization(new Construction(target_thickness));  

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
  /*G4UIExecutive * ui = new G4UIExecutive(argc, argv);
  G4VisManager * visManager = new G4VisExecutive();
  visManager->Initialize();
  G4UImanager * uiManager = G4UImanager::GetUI--finter();
  
  // Running commands to display what I'm doing
  uiManager->ApplyCommand("/run/initialize");
  uiManager->ApplyCommand("/vis/open OGL");
  uiManager->ApplyCommand("/vis/viewer/set/viewpointvector 1 1 1");
  uiManager->ApplyCommand("/vis/drawVolume");
  uiManager->ApplyCommand("/vis/viewer/set/autorefresh true");
  uiManager->ApplyCommand("/vis/scene/add/trajectories smooth");
  uiManager->ApplyCommand("/vis/scene/endOfEventAction accumulate");

  ui->SessionStart(); */
  // End of user interface setup

  //Cause the run manager to generate a single event using the
	//primary generator action registered above.
  runManager->BeamOn(nEvents);

	//After the run is complete, free up the memory used by run 
	//manager and return 

	delete runManager;  

  return 0;
}
