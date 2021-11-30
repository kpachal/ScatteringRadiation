#include "ScatteringRadiation/Run.h"
#include "TFile.h"

RunAction::RunAction(std::string outputFilename) {
  m_outputFilename = outputFilename;
}

RunAction::~RunAction() {

}

void RunAction::BeginOfRunAction(const G4Run*) {

    m_eventNumber = 0;
    m_nFinalState = 0;
    // initialise pointers to zero?
    m_customTree = new TTree("Events","Events");
    m_customTree->Branch("eventNumber",&m_eventNumber,"eventNumber/I");
    m_customTree->Branch("nFinalState",&m_nFinalState,"nFinalState/I");
    m_customTree->Branch("pdgIDFinal",&m_pdgIDFinal);
    m_customTree->Branch("processType",&m_processType);
    m_customTree->Branch("processSubType",&m_processSubType);    
    m_customTree->Branch("posX", &m_posX);    
    m_customTree->Branch("posY", &m_posY);    
    m_customTree->Branch("posZ", &m_posZ);    
    m_customTree->Branch("momX", &m_momX);    
    m_customTree->Branch("momY", &m_momY);
    m_customTree->Branch("momZ", &m_momZ);
    m_customTree->Branch("mass", &m_mass);
    
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
                       std::vector<double> * momZ,
                       std::vector<double> * mass
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
  m_mass = mass;
  m_customTree->Fill();
}

void RunAction::EndOfRunAction(const G4Run*) {

    TFile outfile2(m_outputFilename.c_str(),"RECREATE");
    outfile2.cd();
    m_customTree->Write();
    outfile2.Close();
}