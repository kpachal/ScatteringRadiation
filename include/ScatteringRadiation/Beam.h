#ifndef __BEAM__
#define __BEAM__
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class Beam : public G4VUserPrimaryGeneratorAction{
public:

	Beam();
	virtual ~Beam();

	void GeneratePrimaries(G4Event* anEvent);

private:
	G4ParticleGun* myGun;
};

#endif