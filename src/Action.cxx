#include "ScatteringRadiation/Action.h"

void Action::Build() const {

Beam * beam = new Beam();
SetUserAction(beam);

RunAction * doOnRun = new RunAction();
SetUserAction(doOnRun);

EventAction * doOnEvent = new EventAction();
SetUserAction(doOnEvent);

SteppingAction * doOnStep = new SteppingAction();
SetUserAction(doOnStep);

}