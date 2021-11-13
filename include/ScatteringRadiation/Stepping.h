#ifndef __STEPPING__
#define __STEPPING__

#include "G4UserSteppingAction.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"

#include "ScatteringRadiation/Construction.h"
#include "ScatteringRadiation/Event.h"

class SteppingAction : public G4UserSteppingAction {

  public :
    SteppingAction();
    ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  private :

};

#endif

