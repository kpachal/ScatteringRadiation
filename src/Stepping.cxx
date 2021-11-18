#include "ScatteringRadiation/Stepping.h"
#include "G4String.hh"

// Allows us to access VProcess
#include "G4SteppingManager.hh"

//#include "G4RunManager.hh"

SteppingAction::SteppingAction(EventAction *eventAction) {
    m_EventAction = eventAction;
}

SteppingAction::~SteppingAction() {

}

void SteppingAction::UserSteppingAction(const G4Step* aStep) {

    // Want to get what process happened at this step and note it down.
    G4StepPoint * prePoint = aStep->GetPreStepPoint();
    G4StepPoint * postPoint = aStep->GetPostStepPoint();
   
    // Pre step process isn't defined when electron enters. So just look at post.
    auto postPointProcess = postPoint->GetProcessDefinedStep()->GetProcessName();
    auto postPointProcessType = postPoint->GetProcessDefinedStep()->GetProcessType();
    auto postPointProcessSubType = postPoint->GetProcessDefinedStep()->GetProcessSubType();

    //if (postPointProcessType != 1) G4cout << "Step process: " << postPointProcess << ", " << postPointProcessType << ", " << postPointProcessSubType <<  G4endl;

    // Save processes from this step into ntuple for this event.
    // Check event number
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); 
    int pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();   
    G4AnalysisManager * manager = G4AnalysisManager::Instance();
    manager->FillNtupleIColumn(1, 0, evt);
    manager->FillNtupleIColumn(1, 1, postPointProcessType);
    manager->FillNtupleIColumn(1, 2, postPointProcessSubType);
    manager->FillNtupleIColumn(1, 3, pdgID);
    manager->AddNtupleRow(1);

    // Fill process info into the event information
    m_EventAction->FillProcessInfo(postPointProcessType, postPointProcessSubType);

}