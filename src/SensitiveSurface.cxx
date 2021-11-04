#include "ScatteringRadiation/SensitiveSurface.h"

SensitiveSurface::SensitiveSurface(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveSurface::~SensitiveSurface() 
{}

G4bool SensitiveSurface::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{

    // Check event number
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

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

    // Special print out if something that isn't an e comes out (photon, neutron ...)
    if (pdgID != 11) {
      G4cout << "Special!!! " << pdgID << G4endl;
      G4cout << "Particle definition & name are: " << pdgID << " " << particle_name << G4endl;
      G4cout << "Position and momentum of particle: " << posParticle << " " << momParticle << G4endl;
    }

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

    // Save what I want to my output root file
    G4AnalysisManager * manager = G4AnalysisManager::Instance();
    manager->FillNtupleIColumn(0, evt);
    manager->FillNtupleIColumn(1, pdgID);
    manager->FillNtupleDColumn(2, posParticle[0]);
    manager->FillNtupleDColumn(3, posParticle[1]);
    manager->FillNtupleDColumn(4, posParticle[2]);
    manager->FillNtupleDColumn(5, momParticle[0]);
    manager->FillNtupleDColumn(6, momParticle[1]);
    manager->FillNtupleDColumn(7, momParticle[2]);
    manager->AddNtupleRow(0);

}