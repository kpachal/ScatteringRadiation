#include "ScatteringRadiation/Beam.h"
//#include "G4Electron.hh"

Beam::Beam() {

	//Create particle gun
	// 1 primary vertex per event (is this where we add pileup?)
	//myGun = new G4ParticleGun(1);
	myGun = new G4GeneralParticleSource();
}

Beam::~Beam(){
	delete myGun;
}

void Beam::GeneratePrimaries(G4Event* anEvent){
	
/*	//Specify particle to be emitted
	myGun->SetParticleDefinition(G4Electron::ElectronDefinition());
	
	//Set particle  kinetic energy and direction of travel
	myGun->SetParticlePosition(G4ThreeVector(0,0,-2.0*CLHEP::m));
	myGun->SetParticleEnergy(31.*CLHEP::MeV); 
	myGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1.0));
*/
	//Generate one instance of specified particle
	myGun->GeneratePrimaryVertex(anEvent);
}
