#include "Funcs/Funcs.h"
#include "TLorentzVector.h"

// calculate true vs reco numu spectra for each generator,
// normalised to 1 KT X 1e21 POT, use in other osc fits

// in the case of numu - spectra is with uniform 100% survival prob
// in the case of nue - spectra is with uniform 100% appearance prob

bool nue_mode = false;

void RatesNoOsc(){

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();

  // h_flux is the flux in nu/m2/GeV/POT, this is the total flux nu/cm2/10^21 POT
  double total_flux = h_flux->Integral("width")*1e21/1e4; 
  std::cout << "Total Flux: " << total_flux << " nu/cm^2/10^21 POT" << std::endl; 

  std::vector<std::string> InputFiles_v;

  if(!nue_mode) InputFiles_v = {"rootfiles/GENIEEvents.root","rootfiles/NuWroEvents.root","rootfiles/NEUTEvents.root","rootfiles/GiBUUEvents.root"};
  else InputFiles_v = {"rootfiles/GENIE_NueEvents.root","rootfiles/NuWro_NueEvents.root","rootfiles/NEUT_NueEvents.root","rootfiles/GiBUU_NueEvents.root"};

  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<double> xsec(InputFiles_v.size(),0.0);
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy;
  std::vector<std::vector<TH1D*>> h_RecoEnergy;

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    TFile* f = TFile::Open(InputFiles_v.at(i_f).c_str()) ;
    TTree* t = static_cast<TTree*>(f->Get("eventtree")) ;

    Double_t scale;
    Double_t weight;
    Double_t nu_e;
    Int_t ccnc;
    Int_t nu_pdg;  
    Int_t lepton_pdg;
    TLorentzVector* lepton_p4=0;
    std::vector<int>* pdg=0;
    std::vector<TLorentzVector>* p4=0;

    t->SetBranchAddress("scale",&scale);
    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);

    h_TrueEnergy_RecoEnergy.push_back(std::vector<TH2D*>());

    for(std::string estimator : estimators_str){
      h_TrueEnergy_RecoEnergy.back().push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str(),"Events/KT/GeV^{2}/10^{21} POT;True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);",200,0.1,8.0,200,0.1,8.0));
    }

    int target_nu_pdg = 14;
    if(nue_mode) target_nu_pdg = 12;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(abs(nu_pdg) != target_nu_pdg || ccnc != 1) continue;

      double W = CalcW(pdg,p4);
      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);

      if(generator != "GiBUU") weight = 1;
      weight *= scale*1e38*40;

      if(nprot < 1) continue;

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,i_e);
        h_TrueEnergy_RecoEnergy.back().at(i_e)->Fill(nu_e,nu_e_reco,weight);
      }

    }

  }

  TFile* f_out = nue_mode ? new TFile("rootfiles/NueRates.root","RECREATE") : new TFile("rootfiles/NuMuRates.root","RECREATE");

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      for(int i=1;i<h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetNbinsX()+1;i++)
        for(int j=1;j<h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetNbinsY()+1;j++)
          h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->SetBinContent(i,j,Rate(total_flux,h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetBinContent(i,j)));
      h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->SetTitle("Events in 1 KT X 10^{21} POT");
      h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->Write();
    }
  }
 
 f_out->Close();

}
