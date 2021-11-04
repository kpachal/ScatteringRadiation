#include "ScatteringRadiation/Event.h"

EventAction::EventAction() {

}

EventAction::~EventAction() {

}

void EventAction::BeginOfEventAction(const G4Event*) {
    // things to happen per event, if desired!

}

void EventAction::EndOfEventAction(const G4Event*) {

    // Write to ntuple
}