#include "ScatteringRadiation/SensitiveSurface.h"

SensitiveSurface::SensitiveSurface(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveSurface::~SensitiveSurface() 
{}

G4bool SensitiveSurface::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{

    // Check event number
    G4Track *track = aStep->GetTrack();

    // These are the point where it enters the volume
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    // And here it is leaving
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    // Position and momentum as it enters
    G4ThreeVector posParticle = preStepPoint->GetPosition();
    G4ThreeVector momParticle = preStepPoint->GetMomentum();

    // Make sure it's not from an interaction inside of my
    // detector (even if it is just vacuum)
    if (!aStep->IsFirstStepInVolume()) return false;

    // This would keep only charged particles: not what we want
    /*if (aStep->GetTrack()->GetTrackID()!=1)
      return false;*/

    // Identity of particle
    auto pdgID = track->GetParticleDefinition()->GetPDGEncoding();
    auto particle_name = track->GetParticleDefinition()->GetParticleName();

    // Get what types of interactions the particle has undergone

    auto E=track->GetVertexKineticEnergy();
    auto m=track->GetParticleDefinition()->GetPDGMass();
    auto p=sqrt((E+m)*(E+m)-m*m);

    auto mom=track->GetVertexMomentumDirection();

    /*auto ip=atan2(mom.x(),mom.z())/CLHEP::degree;
    auto oop=atan2(mom.y(),mom.z())/CLHEP::degree;
    G4cout << "In-plane and out-of-plane: " << ip << " " << oop << G4endl;*/

    /*auto gpos=step->GetPreStepPoint()->GetPosition();
    auto pos=step->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(gpos); 
    int cpnr=step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(); */  

    // Could save outputs from here if desired

}