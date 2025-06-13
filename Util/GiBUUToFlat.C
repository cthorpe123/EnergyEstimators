#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

const double _EPSILON_ = 0.000001;
const double mmu = 0.10566; 
const double me =  0.00511; 

void GiBUUToFlat(){

  TFile* f_in = TFile::Open("EventOutput.Pert.00000001.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("RootTuple")) ;

  Double_t        weight;
  Double_t        lepIn_E;
  Double_t        lepIn_Px;
  Double_t        lepIn_Py;
  Double_t        lepIn_Pz;
  Double_t        lepOut_E;
  Double_t        lepOut_Px;
  Double_t        lepOut_Py;
  Double_t        lepOut_Pz;
  vector<int>     *barcode=0;
  vector<double>  *Px=0;
  vector<double>  *Py=0;
  vector<double>  *Pz=0;
  vector<double>  *E=0;

  t_in->SetBranchAddress("weight",&weight);
  t_in->SetBranchAddress("lepIn_E",&lepIn_E);
  t_in->SetBranchAddress("lepIn_Px",&lepIn_Px);
  t_in->SetBranchAddress("lepIn_Py",&lepIn_Py);
  t_in->SetBranchAddress("lepIn_Pz",&lepIn_Pz);
  t_in->SetBranchAddress("lepOut_E",&lepOut_E);
  t_in->SetBranchAddress("lepOut_Px",&lepOut_Px);
  t_in->SetBranchAddress("lepOut_Py",&lepOut_Py);
  t_in->SetBranchAddress("lepOut_Pz",&lepOut_Pz);
  t_in->SetBranchAddress("barcode",&barcode);
  t_in->SetBranchAddress("Px",&Px);
  t_in->SetBranchAddress("Py",&Py);
  t_in->SetBranchAddress("Pz",&Pz);
  t_in->SetBranchAddress("E",&E);

  TFile* f_out = new TFile("GiBUUEvents.root","RECREATE"); 
  TTree* t_out = new TTree("eventtree","eventtree"); 

  Double_t nu_e;
  Int_t ccnc;
  Int_t nu_pdg;  

  Int_t lepton_pdg;
  TLorentzVector lepton_p4;

  std::vector<int> pdg;
  std::vector<TLorentzVector> p4;

  t_out->Branch("nu_e",&nu_e);
  t_out->Branch("nu_pdg",&nu_pdg);
  t_out->Branch("ccnc",&ccnc);
  t_out->Branch("lepton_pdg",&lepton_pdg); 
  t_out->Branch("lepton_p4",&lepton_p4);
  t_out->Branch("pdg",&pdg);
  t_out->Branch("p4",&p4);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){
    //if(ievent > 10000) break;
    if(ievent % 10000 == 0) std::cout << "Event " << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    p4.clear();
    pdg.clear();

    lepton_p4 = TLorentzVector(lepOut_Px,lepOut_Py,lepOut_Pz,lepOut_E);

    // Extremely hacky way of getting the initial nu mass/pdg/ccnc status
    lepton_pdg = 0;
    if(abs(lepton_p4.M() - mmu) < _EPSILON_) lepton_pdg = 13;  
    else if(abs(lepton_p4.M() - me) < _EPSILON_) lepton_pdg = 11;  
    if(lepton_pdg == 0) ccnc = 0;
    else ccnc = 1;
    if(lepton_pdg == 13 && ccnc == 1) nu_pdg = 14;
    else if(lepton_pdg == 11 && ccnc == 1) nu_pdg = 12;

    nu_e = lepIn_E;

    for(size_t i_p=0;i_p<barcode->size();i_p++){
      pdg.push_back(barcode->at(i_p));
      p4.push_back(TLorentzVector(Px->at(i_p),Py->at(i_p),Pz->at(i_p),E->at(i_p)));
    }   

    t_out->Fill(); 

  }

  t_out->Write();
  f_out->Close();

  f_in->Close();

}

