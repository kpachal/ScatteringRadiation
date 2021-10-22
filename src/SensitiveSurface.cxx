#include "ScatteringRadiation/SensitiveSurface.h"

SensitiveSurface::SensitiveSurface(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveSurface::~SensitiveSurface() 
{}

G4bool SensitiveSurface::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{

    G4Track *track = aStep->GetTrack();

    // These are the point where it enters the volume
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    // And here it is leaving
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    // Position and momentum as it enters
    G4ThreeVector posParticle = preStepPoint->GetPosition();
    G4ThreeVector momParticle = preStepPoint->GetMomentum();

    G4cout << "Position and momentum of particle: " << posParticle << " " << momParticle << G4endl;

    // Alternatively ....
    if (!aStep->IsFirstStepInVolume())
      return false;
    // Check these
    if (aStep->GetTrack()->GetTrackID()!=1)
      return false;

    auto E=track->GetVertexKineticEnergy();
    auto m=track->GetParticleDefinition()->GetPDGMass();
    auto p=sqrt((E+m)*(E+m)-m*m);

    auto mom=track->GetVertexMomentumDirection();

    auto ip=atan2(mom.x(),mom.z())/CLHEP::degree;
    auto oop=atan2(mom.y(),mom.z())/CLHEP::degree;

    /*auto gpos=step->GetPreStepPoint()->GetPosition();
    auto pos=step->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(gpos); 
    int cpnr=step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(); */    

}