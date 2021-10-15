// Goal: hit electron beam on a target and model all outgoing particles
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

// My includes
#include "ScatteringRadiation/Foil.h"
#include "ScatteringRadiation/Beam.h"

// GEANT includes
#include <FTFP_BERT.hh>
#include "G4RunManager.hh"


int main(int argc, char* argv[])
{
	//Get instance of runmanager
	G4RunManager * runManager = new G4RunManager;

	// Detector construction
	runManager->SetUserInitialization(new Foil());  

  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  // Primary generator action
  runManager->SetUserAction(new Beam());

  // Initialize G4 kernel
  
	runManager->Initialize();

  //Cause the run manager to generate a single event using the
	//primary generator action registered above.
  runManager->BeamOn(1);

	//After the run is complete, free up the memory used by run 
	//manager and return 

	delete runManager;  

  return 0;
}
