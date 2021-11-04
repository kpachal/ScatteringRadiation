#ifndef __SURFACE__
#define __SURFACE__

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "g4root.hh"

class SensitiveSurface : public G4VSensitiveDetector
{
public:
    SensitiveSurface(G4String);
    ~SensitiveSurface();
    
private:
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
                        
};

#endif