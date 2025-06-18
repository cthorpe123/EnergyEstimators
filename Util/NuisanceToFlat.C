#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

const std::string generator = "GiBUU";

// GiBUU ntuple is made of many files - scaling branch needs to be corrected
// for this, use this variable to do so
const int num_files = 1;

void NuisanceToFlat(){

  TFile* f_in = TFile::Open(("../rootfiles/Nuisance/"+generator+".flat.root").c_str());
  TTree* t_in = static_cast<TTree*>(f_in->Get("FlatTree_VARS")) ;

  Char_t          cc;
  Int_t           PDGnu;
  Float_t         Enu_true;
  Int_t           nfsp;
  Float_t         px[100];   //[nfsp]
  Float_t         py[100];   //[nfsp]
  Float_t         pz[100];   //[nfsp]
  Float_t         E[100];   //[nfsp]
  Int_t           pdg_in[100];   //[nfsp]
  Float_t         Weight;
  Double_t        fScaleFactor;

  t_in->SetBranchAddress("cc", &cc);
  t_in->SetBranchAddress("PDGnu", &PDGnu);
  t_in->SetBranchAddress("Enu_true", &Enu_true);
  t_in->SetBranchAddress("nfsp", &nfsp);
  t_in->SetBranchAddress("px", px);
  t_in->SetBranchAddress("py", py);
  t_in->SetBranchAddress("pz", pz);
  t_in->SetBranchAddress("E", E);
  t_in->SetBranchAddress("pdg", pdg_in);
  t_in->SetBranchAddress("Weight", &Weight);
  t_in->SetBranchAddress("fScaleFactor", &fScaleFactor);

  TFile* f_out = new TFile((generator+"Events_tmp.root").c_str(),"RECREATE"); 
  TTree* t_out = new TTree("eventtree_tmp","eventtree_tmp"); 
  
  Double_t scale;
  Double_t weight;
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

  t_out->Branch("scale",&scale);
  t_out->Branch("weight",&weight);
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

    scale = fScaleFactor; 
    weight = Weight;
    nu_e = Enu_true;
    ccnc = cc;
    nu_pdg = PDGnu;

    //lepton_pdg = fspl;
    //lepton_p4 = TLorentzVector(pxl,pyl,pzl,El);

    scale /= num_files; 

    int nlep = 0;
    for(size_t i_p=0;i_p<nfsp;i_p++){
      if(abs(pdg_in[i_p]) == 13 || abs(pdg_in[i_p]) == 11){
        lepton_pdg = pdg_in[i_p];
        lepton_p4 = TLorentzVector(px[i_p],py[i_p],pz[i_p],E[i_p]);
        nlep++;
      } 
      pdg.push_back(pdg_in[i_p]);
      p4_x.push_back(px[i_p]);
      p4_y.push_back(py[i_p]);
      p4_z.push_back(pz[i_p]);
      p4_E.push_back(E[i_p]);
    }   

    if(nlep > 1){
      std::cout << "Event has more than one FS lepton " << nlep << std::endl;
      continue;
    }

    t_out->Fill(); 

  }

  t_out->Write();
  f_out->Close();

  f_in->Close();

}

void FixP4(){

  TFile* f_in = TFile::Open((generator+"Events_tmp.root").c_str());
  TTree* t_in = static_cast<TTree*>(f_in->Get("eventtree_tmp")) ;

  Double_t scale;
  Double_t weight;
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

  t_in->SetBranchAddress("scale",&scale);
  t_in->SetBranchAddress("weight",&weight);
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

  TFile* f_out = new TFile((generator+"Events.root").c_str(),"RECREATE"); 
  TTree* t_out = new TTree("eventtree","eventtree"); 

  TLorentzVector lepton_p4;
  std::vector<int> pdg;
  std::vector<TLorentzVector> p4;

  t_out->Branch("scale",&scale);
  t_out->Branch("weight",&weight);
  t_out->Branch("nu_e",&nu_e);
  t_out->Branch("nu_pdg",&nu_pdg);
  t_out->Branch("ccnc",&ccnc);
  t_out->Branch("lepton_pdg",&lepton_pdg); 
  t_out->Branch("lepton_p4",&lepton_p4);
  t_out->Branch("pdg",&pdg);
  t_out->Branch("p4",&p4);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //std::cout << ievent << std::endl;
    if(ievent % 10000 == 0) std::cout << "Event " << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    p4.clear();

    lepton_p4 = *in_lepton_p4;
    pdg = *in_pdg;

    for(size_t i=0;i<in_p4_x->size();i++) p4.push_back(TLorentzVector(in_p4_x->at(i),in_p4_y->at(i),in_p4_z->at(i),in_p4_E->at(i)));

    t_out->Fill();

  }

  t_out->Write();
  f_out->Close();
  f_in->Close();

}

