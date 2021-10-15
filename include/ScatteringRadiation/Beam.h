#ifndef __BEAM__
#define __BEAM__
#include "G4VUserPrimaryGeneratorAction.hh"

class Beam : public G4VUserPrimaryGeneratorAction{
public:

	Beam() {};
	virtual ~Beam() {};

	void GeneratePrimaries(G4Event* anEvent);

};

#endif