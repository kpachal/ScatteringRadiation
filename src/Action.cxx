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
    if (m_save_highangle_only) doOnTrack->setTrackHighAngleOnly(true, 20.0); // in degrees
    SetUserAction(doOnTrack);

    SteppingAction * doOnStep = new SteppingAction(doOnEvent);
    SetUserAction(doOnStep);

}