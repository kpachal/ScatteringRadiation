#include "ScatteringRadiation/SensitiveSurface.h"

SensitiveSurface::SensitiveSurface(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveSurface::~SensitiveSurface() 
{}

G4bool SensitiveSurface::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{



}