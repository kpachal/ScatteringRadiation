#ifndef __RUN__
#define __RUN__

#include "G4UserRunAction.hh"
#include "g4root.hh"
#include "TTree.h"

class RunAction : public G4UserRunAction {

  public :
    RunAction();
    ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    const void FillEventInfo(int &eventNumber, 
                       int &nFinalState,
                       std::vector<int> * pdgIDFinal,
                       std::vector<int> * processType,
                       std::vector<int> * processSubType,
                       std::vector<double> * posX,
                       std::vector<double> * posY,
                       std::vector<double> * posZ,
                       std::vector<double> * momX,
                       std::vector<double> * momY,
                       std::vector<double> * momZ
                       );


  private :

    // Tree and branches that I am managing myself
    TTree * m_customTree;
    int m_eventNumber;
    int m_nFinalState;
    std::vector<int> * m_pdgIDFinal;   ///< IDs of stable particles that leave event
    std::vector<int> * m_processType;   ///< Type of each process that occurred
    std::vector<int> * m_processSubType;   ///< Subtype of each process that occurred
    std::vector<double> * m_posX;   ///< x position of each particle as it leaves sensitive volume (could switch to tracking action)
    std::vector<double> * m_posY;   ///< y position of each particle as it leaves sensitive volume
    std::vector<double> * m_posZ;   ///< z position of each particle as it leaves sensitive volume
    std::vector<double> * m_momX;   ///< x momentum of each particle as it leaves sensitive volume (could switch to tracking action)
    std::vector<double> * m_momY;   ///< y momentum of each particle as it leaves sensitive volume
    std::vector<double> * m_momZ;   ///< z momentum of each particle as it leaves sensitive volume    

};

#endif

