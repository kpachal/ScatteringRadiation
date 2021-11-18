#include "ScatteringRadiation/Event.h"

EventAction::EventAction(RunAction* runAction) {
    m_RunAction = runAction;
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
}

void EventAction::FillProcessInfo(int processType, int processSubType) {
   m_processType.push_back(processType);
   m_processSubType.push_back(processSubType);
}

void EventAction::FillTrackInfo(int pdgID, 
                              double posX, double posY, double posZ, 
                              double momX, double momY, double momZ) 
{
   m_pdgIDFinal.push_back(pdgID);
   m_posX.push_back(posX);
   m_posY.push_back(posY);
   m_posZ.push_back(posZ);
   m_momX.push_back(momX);
   m_momY.push_back(momY);
   m_momZ.push_back(momZ);
}

void EventAction::EndOfEventAction(const G4Event*) {

    // Which event?
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); 

    // Count final state particles
    int nFinalState = m_pdgIDFinal.size();

    // Print to check set correctly
/*    std::cout << evt << std::endl;
    std::cout << nFinalState << std::endl;    
    for (auto i: m_pdgIDFinal) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_processType) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_processSubType) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_posX) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_posY) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_posZ) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_momX) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_momY) std::cout << i << ' '; std::cout << std::endl;
    for (auto i: m_momZ) std::cout << i << ' '; std::cout << std::endl; */

    m_RunAction->FillEventInfo(evt, nFinalState, &m_pdgIDFinal, 
                               &m_processType, &m_processSubType,
                               &m_posX, &m_posY, &m_posZ,
                               &m_momX, &m_momY, &m_momZ);

}