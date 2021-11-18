#ifndef __EVENT__
#define __EVENT__

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "g4root.hh"
#include "G4RunManager.hh"

#include "ScatteringRadiation/Run.h"

class EventAction : public G4UserEventAction {

  public :
    EventAction(RunAction* runAction);
    ~EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void FillProcessInfo(int processType, int processSubType);
    void FillTrackInfo(int pdgID, double posX, double posY, double posZ, 
                              double momX, double momY, double momZ);
  private :

    RunAction * m_RunAction;

    // Quantities I want to track across the event:
    std::vector<int> m_pdgIDFinal;   ///< IDs of stable particles that leave event
    std::vector<int> m_processType;   ///< Type of each process that occurred
    std::vector<int> m_processSubType;   ///< Subtype of each process that occurred
    std::vector<double> m_posX;   ///< x position of each particle as it leaves sensitive volume (could switch to tracking action)
    std::vector<double> m_posY;   ///< y position of each particle as it leaves sensitive volume
    std::vector<double> m_posZ;   ///< z position of each particle as it leaves sensitive volume
    std::vector<double> m_momX;   ///< x momentum of each particle as it leaves sensitive volume (could switch to tracking action)
    std::vector<double> m_momY;   ///< y momentum of each particle as it leaves sensitive volume
    std::vector<double> m_momZ;   ///< z momentum of each particle as it leaves sensitive volume

};

#endif

