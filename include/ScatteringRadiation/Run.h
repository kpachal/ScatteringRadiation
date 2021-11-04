#ifndef __RUN__
#define __RUN__

#include "G4UserRunAction.hh"
#include "g4root.hh"

class RunAction : public G4UserRunAction {

  public :
    RunAction();
    ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

  private :

};

#endif

