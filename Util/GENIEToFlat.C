#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

void GENIEToFlat(){

  TFile* f_in = TFile::Open("../GENIE/gntp.0.gst.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("gst")) ;

  Double_t        wght;
  Int_t           neu;
  Bool_t          cc;
  Double_t        Ev;
  Int_t           fspl;
  Double_t        El;
  Double_t        pxl;
  Double_t        pyl;
  Double_t        pzl;
  Int_t           nf;
  Int_t           pdgf[44];   //[nf]
  Double_t        Ef[44];   //[nf]
  Double_t        pxf[44];   //[nf]
  Double_t        pyf[44];   //[nf]
  Double_t        pzf[44];   //[nf]

  t_in->SetBranchAddress("wght",&wght);
  t_in->SetBranchAddress("neu",&neu);
  t_in->SetBranchAddress("cc",&cc);
  t_in->SetBranchAddress("Ev",&Ev);
  t_in->SetBranchAddress("fspl",&fspl);
  t_in->SetBranchAddress("El",&El);
  t_in->SetBranchAddress("pxl",&pxl);
  t_in->SetBranchAddress("pyl",&pyl);
  t_in->SetBranchAddress("pzl",&pzl);
  t_in->SetBranchAddress("nf",&nf);
  t_in->SetBranchAddress("pdgf",pdgf);
  t_in->SetBranchAddress("Ef",Ef);
  t_in->SetBranchAddress("pxf",pxf);
  t_in->SetBranchAddress("pyf",pyf);
  t_in->SetBranchAddress("pzf",pzf);

  TFile* f_out = new TFile("GENIEEvents.root","RECREATE"); 
  TTree* t_out = new TTree("eventtree","eventtree"); 
  
  Double_t weight;
  Double_t nu_e;
  Int_t ccnc;
  Int_t nu_pdg;  

  Int_t lepton_pdg;
  TLorentzVector lepton_p4;

  std::vector<int> pdg;
  std::vector<TLorentzVector> p4;

  t_out->Branch("weight",&weight);
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

    weight = wght;
    nu_e = Ev;
    ccnc = cc;
    nu_pdg = neu;

    lepton_pdg = fspl;
    lepton_p4 = TLorentzVector(pxl,pyl,pzl,El);

    for(size_t i_p=0;i_p<nf;i_p++){
      pdg.push_back(pdgf[i_p]);
      p4.push_back(TLorentzVector(pxf[i_p],pyf[i_p],pzf[i_p],Ef[i_p]));
    }   

    t_out->Fill(); 

  }

  t_out->Write();
  f_out->Close();

  f_in->Close();

}

