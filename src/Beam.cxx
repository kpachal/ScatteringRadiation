#include "ScatteringRadiation/Beam.h"
#include "G4ParticleGun.hh"
#include "G4Proton.hh"


void Beam::GeneratePrimaries(G4Event* anEvent){

	//Create particle gun
	G4ParticleGun* myGun = new G4ParticleGun();
	
	//Specify particle to be emitted
	myGun->SetParticleDefinition(G4Proton::ProtonDefinition());
	
	//Set particle  kinetic energy and direction of travel
	myGun->SetParticleEnergy(50.*CLHEP::keV); 
	myGun->SetParticleMomentumDirection(G4ThreeVector(1.0,0,0));
	
	//Generate one instance of specified particle
	myGun->GeneratePrimaryVertex(anEvent);
}
