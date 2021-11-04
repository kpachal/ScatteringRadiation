#include "ScatteringRadiation/Run.h"

RunAction::RunAction() {

}

RunAction::~RunAction() {

}

void RunAction::BeginOfRunAction(const G4Run*) {
    G4AnalysisManager * manager = G4AnalysisManager::Instance();
    manager->OpenFile("output.root");

    manager->CreateNtuple("Particles","Particles");
    manager->CreateNtupleIColumn("eventNumber"); //0
    manager->CreateNtupleIColumn("pdgID"); //1
    manager->CreateNtupleDColumn("posX"); //2
    manager->CreateNtupleDColumn("posY"); //3
    manager->CreateNtupleDColumn("posZ"); //4
    manager->CreateNtupleDColumn("momX"); //5
    manager->CreateNtupleDColumn("momY"); //6
    manager->CreateNtupleDColumn("momZ"); //7
    manager->FinishNtuple(0);

}

void RunAction::EndOfRunAction(const G4Run*) {

    G4AnalysisManager * manager = G4AnalysisManager::Instance();
    manager->Write();
    manager->CloseFile("output.root");
}