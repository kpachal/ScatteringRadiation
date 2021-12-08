#include "ScatteringRadiation/Action.h"

void Action::Build() const {

    Beam * beam = new Beam();
    SetUserAction(beam);

    RunAction * doOnRun = new RunAction(m_outputFilename);
    SetUserAction(doOnRun);

    EventAction * doOnEvent = new EventAction(doOnRun);
    doOnEvent->setDoOnly(m_do_save_only, m_PDGID_only);
    SetUserAction(doOnEvent);

    TrackingAction * doOnTrack = new TrackingAction(doOnEvent);
    SetUserAction(doOnTrack);

    SteppingAction * doOnStep = new SteppingAction(doOnEvent);
    SetUserAction(doOnStep);

}