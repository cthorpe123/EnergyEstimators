// Need to setup a different version of root to read the NEUT ntuple. You also need two terminals
// One the first one, we need a different root version
// unsetup root 
// setup root v6_28_12 -q e20:p3915:prof 
// Then run NEUTToFlat to with this setup
// The on a second one, use the normal setup with root v6_12_06a and run the FixP4 function below 
// This second function is needed as root v6_28_12 doesn't seem to like streaming vector<TLorentzVector>

#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

void NEUTToFlat(){

  TFile* f_in = TFile::Open("NEUT.flat.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("FlatTree_VARS")) ;

  Float_t         Enu_true;
  Int_t           PDGnu;
  Int_t           nfsp;
  Float_t         px[10];   //[nfsp]
  Float_t         py[10];   //[nfsp]
  Float_t         pz[10];   //[nfsp]
  Float_t         E[10];   //[nfsp]
  Int_t           in_pdg[10];   //[nfsp]

  t_in->SetBranchAddress("Enu_true",&Enu_true);
  t_in->SetBranchAddress("PDGnu",&PDGnu);
  t_in->SetBranchAddress("nfsp",&nfsp);
  t_in->SetBranchAddress("px",px);
  t_in->SetBranchAddress("py",py);
  t_in->SetBranchAddress("pz",pz);
  t_in->SetBranchAddress("E",E);
  t_in->SetBranchAddress("pdg",in_pdg);

  Double_t nu_e;
  Int_t ccnc;
  Int_t nu_pdg;  

  Int_t lepton_pdg;
  TLorentzVector lepton_p4;

  std::vector<int> pdg;
  //std::vector<TLorentzVector> p4;
  std::vector<double> p4_x;
  std::vector<double> p4_y;
  std::vector<double> p4_z;
  std::vector<double> p4_E;

  TFile* f_out = new TFile("NEUTEvents_tmp.root","RECREATE"); 
  TTree* t_out = new TTree("eventtree_tmp","eventtree_tmp"); 

  t_out->Branch("nu_e",&nu_e);
  t_out->Branch("nu_pdg",&nu_pdg);
  t_out->Branch("ccnc",&ccnc);
  t_out->Branch("lepton_pdg",&lepton_pdg); 
  t_out->Branch("lepton_p4",&lepton_p4);
  t_out->Branch("pdg",&pdg);
  t_out->Branch("p4_x",&p4_x);
  t_out->Branch("p4_y",&p4_y);
  t_out->Branch("p4_z",&p4_z);
  t_out->Branch("p4_E",&p4_E);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 10000) break;
    if(ievent % 10000 == 0) std::cout << "Event " << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    p4_x.clear();
    p4_y.clear();
    p4_z.clear();
    p4_E.clear();
    pdg.clear();

    if(nfsp == 0) continue;

    nu_e = Enu_true;
    nu_pdg = PDGnu;  

    lepton_pdg = in_pdg[0]; 
    lepton_p4 = TLorentzVector(px[0],py[0],pz[0],E[0]);

    // Hacky way to check if event is cc or nc
    ccnc = nu_pdg != lepton_pdg; 

    for(int i=1;i<nfsp;i++){
      pdg.push_back(in_pdg[i]);
      //p4.push_back(TLorentzVector(px[i],py[i],pz[i],E[i]));
      p4_x.push_back(px[i]);
      p4_y.push_back(py[i]);
      p4_z.push_back(pz[i]);
      p4_E.push_back(E[i]);
    }

    t_out->Fill(); 

  }


  t_out->Write();
  f_out->Close();
  f_in->Close();

}


void FixP4(){

  TFile* f_in = TFile::Open("NEUTEvents_tmp.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("eventtree_tmp")) ;

  Double_t nu_e;
  Int_t ccnc;
  Int_t nu_pdg;  

  Int_t lepton_pdg;
  TLorentzVector* in_lepton_p4=0;

  std::vector<int>* in_pdg=0;
  //std::vector<TLorentzVector> p4;
  std::vector<double>* in_p4_x=0;
  std::vector<double>* in_p4_y=0;
  std::vector<double>* in_p4_z=0;
  std::vector<double>* in_p4_E=0;

  t_in->SetBranchAddress("nu_e",&nu_e);
  t_in->SetBranchAddress("nu_pdg",&nu_pdg);
  t_in->SetBranchAddress("ccnc",&ccnc);
  t_in->SetBranchAddress("lepton_pdg",&lepton_pdg); 
  t_in->SetBranchAddress("lepton_p4",&in_lepton_p4);
  t_in->SetBranchAddress("pdg",&in_pdg);
  t_in->SetBranchAddress("p4_x",&in_p4_x);
  t_in->SetBranchAddress("p4_y",&in_p4_y);
  t_in->SetBranchAddress("p4_z",&in_p4_z);
  t_in->SetBranchAddress("p4_E",&in_p4_E);

  TFile* f_out = new TFile("NEUTEvents.root","RECREATE"); 
  TTree* t_out = new TTree("eventtree","eventtree"); 

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

    if(ievent % 10000 == 0) std::cout << "Event " << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    lepton_p4 = *in_lepton_p4;
    pdg = *in_pdg;

    for(size_t i=0;i<in_p4_x->size();i++) p4.push_back(TLorentzVector(in_p4_x->at(i),in_p4_y->at(i),in_p4_z->at(i),in_p4_E->at(i)));

    t_out->Fill();

  }

  t_out->Write();
  f_out->Close();
  f_in->Close();

}




