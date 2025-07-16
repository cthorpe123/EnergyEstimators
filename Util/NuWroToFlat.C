R__LOAD_LIBRARY(event1.so);

#include "event1.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

void NuWroToFlat(){

  TFile* f_in = TFile::Open("/gluster/data/dune/cthorpe/DIS/NuWro/NuWroCard_CC_Ar_numu.prep.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("treeout")) ;

  // Because we've loaded event1.so we can read the entire event object on a single branch
  event *e = new event();
  t_in->SetBranchAddress("e",&e);

  TFile* f_out = new TFile("NuWroEvents.root","RECREATE"); 
  TTree* t_out = new TTree("eventtree","eventtree"); 

  Double_t scale = 1;
  Double_t weight; 
  Double_t nu_e;
  Int_t ccnc;
  Int_t nu_pdg;  
 
  Int_t lepton_pdg;
  TLorentzVector lepton_p4;

  std::vector<int> pdg;
  std::vector<TLorentzVector> p4;

  std::vector<int> pdg_nofsi;
  std::vector<TLorentzVector> p4_nofsi;
 
  t_out->Branch("scale",&scale);
  t_out->Branch("weight",&weight);
  t_out->Branch("nu_e",&nu_e);
  t_out->Branch("nu_pdg",&nu_pdg);
  t_out->Branch("ccnc",&ccnc);
  t_out->Branch("lepton_pdg",&lepton_pdg); 
  t_out->Branch("lepton_p4",&lepton_p4);
  t_out->Branch("pdg_nofsi",&pdg_nofsi);
  t_out->Branch("p4_nofsi",&p4_nofsi);
  t_out->Branch("pdg",&pdg);
  t_out->Branch("p4",&p4);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){
    //if(ievent > 10000) break;
    if(ievent % 10000 == 0) std::cout << "Event " << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);
    
    p4.clear();
    pdg.clear();
    p4_nofsi.clear();
    pdg_nofsi.clear();

    weight = e->weight;
    nu_e = e->in.at(0).t/1e3;
    ccnc = e->flag.cc;
    nu_pdg = e->in.at(0).pdg;

    lepton_pdg = e->post.at(0).pdg;
    lepton_p4 = TLorentzVector(e->post.at(0).x/1e3,e->post.at(0).y/1e3,e->post.at(0).z/1e3,e->post.at(0).t/1e3);

    for(size_t i_p=1;i_p<e->out.size();i_p++){
      pdg_nofsi.push_back(e->out.at(i_p).pdg);
      p4_nofsi.push_back(TLorentzVector(e->out.at(i_p).x/1e3,e->out.at(i_p).y/1e3,e->out.at(i_p).z/1e3,e->out.at(i_p).t/1e3));
    }   

    for(size_t i_p=1;i_p<e->post.size();i_p++){
      pdg.push_back(e->post.at(i_p).pdg);
      p4.push_back(TLorentzVector(e->post.at(i_p).x/1e3,e->post.at(i_p).y/1e3,e->post.at(i_p).z/1e3,e->post.at(i_p).t/1e3));
    }   

    t_out->Fill(); 
     
  }

  t_out->Write();
  f_out->Close();

  f_in->Close();

}

