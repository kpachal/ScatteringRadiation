#include "ScatteringRadiation/Run.h"
#include "TFile.h"

RunAction::RunAction() {

}

RunAction::~RunAction() {

}

void RunAction::BeginOfRunAction(const G4Run*) {
    
    // Examples of how to make vectors instead
    // https://github.com/jintonic/gears/blob/master/gears.cc

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

    manager->CreateNtuple("Processes","Processes");
    manager->CreateNtupleIColumn("eventNumber"); //0
    manager->CreateNtupleIColumn("processType"); //1
    manager->CreateNtupleIColumn("processSubType"); //2
    manager->CreateNtupleIColumn("outgoingParticlePDGID"); //3
    manager->FinishNtuple(1);

    std::vector<int> pdgIDFinal;   ///< IDs of stable particles that leave event
    std::vector<int> processType;   ///< Type of each process that occurred
    std::vector<int> processSubType;   ///< Subtype of each process that occurred
    std::vector<double> posX;   ///< x position of each particle as it leaves sensitive volume (could switch to tracking action)
    std::vector<double> posY;   ///< y position of each particle as it leaves sensitive volume
    std::vector<double> posZ;   ///< z position of each particle as it leaves sensitive volume
    std::vector<double> momX;   ///< x momentum of each particle as it leaves sensitive volume (could switch to tracking action)
    std::vector<double> momY;   ///< y momentum of each particle as it leaves sensitive volume
    std::vector<double> momZ;   ///< z momentum of each particle as it leaves sensitive volume
    manager->CreateNtuple("Events","Events");
    manager->CreateNtupleIColumn("eventNumber"); // 0  
    manager->CreateNtupleIColumn("nFinalState"); // 1
    manager->CreateNtupleIColumn("pdgIDFinal",pdgIDFinal); // 2
    manager->CreateNtupleIColumn("processType",processType); // 3
    manager->CreateNtupleIColumn("processSubType",processSubType); // 4
    manager->CreateNtupleDColumn("posX", posX); // 5
    manager->CreateNtupleDColumn("posY", posY); // 6
    manager->CreateNtupleDColumn("posZ", posZ); // 7
    manager->CreateNtupleDColumn("momX", momX); // 8
    manager->CreateNtupleDColumn("momY", momY); // 9
    manager->CreateNtupleDColumn("momZ", momZ); // 10  
    manager->FinishNtuple(2); 

    m_eventNumber = 0;
    m_nFinalState = 0;
    // initialise pointers to zero?
    m_customTree = new TTree("Events","Events");
    m_customTree->Branch("eventNumber",&m_eventNumber,"eventNumber/I");
    m_customTree->Branch("eventNumber",&m_nFinalState,"nFinalState/I");
    m_customTree->Branch("pdgIDFinal",&m_pdgIDFinal);
    m_customTree->Branch("processType",&m_processType);
    m_customTree->Branch("processSubType",&m_processSubType);    
    m_customTree->Branch("posX", &m_posX);    
    m_customTree->Branch("posY", &m_posY);    
    m_customTree->Branch("posZ", &m_posZ);    
    m_customTree->Branch("momX", &m_momX);    
    m_customTree->Branch("momY", &m_momY);
    m_customTree->Branch("momZ", &m_momZ);

}

const void RunAction::FillEventInfo(int &eventNumber, 
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
                       ) 
{
  m_eventNumber = eventNumber;
  m_nFinalState = nFinalState;
  m_pdgIDFinal = pdgIDFinal;
  m_processType = processType;
  m_processSubType = processSubType;
  m_posX = posX;
  m_posY = posY;
  m_posZ = posZ;
  m_momX = momX;
  m_momY = momY;
  m_momZ = momZ;
  m_customTree->Fill();
}

void RunAction::EndOfRunAction(const G4Run*) {

    G4AnalysisManager * manager = G4AnalysisManager::Instance();
    manager->Write();
    manager->CloseFile("output.root");

    TFile outfile2("output_byevent.root","RECREATE");
    outfile2.cd();
    m_customTree->Write();
    outfile2.Close();
}