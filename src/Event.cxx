#include "ScatteringRadiation/Event.h"

EventAction::EventAction(RunAction* runAction) {
    m_RunAction = runAction;
    m_do_save_only = false;
    m_save_only = 0;
}

EventAction::~EventAction() {

}

void EventAction::BeginOfEventAction(const G4Event*) {
    // Reset vectors to be empty
    m_pdgIDFinal.clear();
    m_processType.clear();
    m_processSubType.clear();
    m_posX.clear();
    m_posY.clear();
    m_posZ.clear();
    m_momX.clear();
    m_momY.clear();
    m_momZ.clear();
    m_mass.clear();
}

void EventAction::FillProcessInfo(int processType, int processSubType) {
   m_processType.push_back(processType);
   m_processSubType.push_back(processSubType);
}

void EventAction::setDoOnly(bool doOnly, int PDGIDOnly) {
    m_do_save_only = doOnly;
    m_save_only = PDGIDOnly;
}


void EventAction::FillTrackInfo(int pdgID, 
                              double posX, double posY, double posZ, 
                              double momX, double momY, double momZ, double mass) 
{

   // If we requested to store only one type of particle,
   // return here if this is not it.
   if (m_do_save_only && pdgID != m_save_only) return;

   m_pdgIDFinal.push_back(pdgID);
   m_posX.push_back(posX);
   m_posY.push_back(posY);
   m_posZ.push_back(posZ);
   m_momX.push_back(momX);
   m_momY.push_back(momY);
   m_momZ.push_back(momZ);
   m_mass.push_back(mass);
   
}

void EventAction::EndOfEventAction(const G4Event*) {

    // Which event?
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); 

    // Count final state particles
    int nFinalState = m_pdgIDFinal.size();

    // If there aren't any final state particles 
    // because we are saving only something rare, don't write this event.
    if (nFinalState == 0) return;

    // Print to check set correctly
    m_RunAction->FillEventInfo(evt, nFinalState, &m_pdgIDFinal, 
                               &m_processType, &m_processSubType,
                               &m_posX, &m_posY, &m_posZ,
                               &m_momX, &m_momY, &m_momZ, &m_mass);

}