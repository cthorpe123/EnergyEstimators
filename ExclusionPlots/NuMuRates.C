#include "Funcs/Funcs.h"
#include "Funcs/OscFitter.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

OscModel osc_model;

// Detector exposure to normalise plots to in KT x 10^21 POT
double fid_mass = 1; // active mass in KT
double POT = 1; // POT in 10^21

void NuMuRates(){

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  TH1D* h_flux_osc = static_cast<TH1D*>(f_flux->Get("numu_fluxosc"));
  h_flux->SetDirectory(0);
  h_flux_osc->SetDirectory(0);
  f_flux->Close();

  // h_flux is the flux in nu/m2/GeV/POT, this is the total flux nu/cm2/10^21 POT
  double total_flux = h_flux->Integral("width")*1e21/1e4; 
  std::cout << "Total Flux: " << total_flux << " nu/cm^2/10^21 POT" << std::endl; 

  std::vector<std::string> InputFiles_v = {"GENIEEvents.root","NuWroEvents.root","NEUTEvents.root","GiBUUEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSqCV;
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_Plus;
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_Minus;

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    h_RecoEnergy_DeltaMSqCV.push_back(std::vector<TH1D*>());
    h_RecoEnergy_DeltaMSq_Plus.push_back(std::vector<TH1D*>());
    h_RecoEnergy_DeltaMSq_Minus.push_back(std::vector<TH1D*>());

    int nbins = 100;
    for(std::string estimator : estimators_str){
      h_RecoEnergy_DeltaMSqCV.back().push_back(new TH1D((generator+"_RecoEnergy_DeltaMSqCV_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/10^{21} POT",nbins,0.1,8.0));
      h_RecoEnergy_DeltaMSq_Plus.back().push_back(new TH1D((generator+"_RecoEnergy_DeltaMSq_Plus_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/10^{21} POT",nbins,0.1,8.0));
      h_RecoEnergy_DeltaMSq_Minus.back().push_back(new TH1D((generator+"_RecoEnergy_DeltaMSq_Minus_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/10^{21} POT",nbins,0.1,8.0));
    }

    TFile* f = TFile::Open(("/gluster/data/dune/cthorpe/DIS/"+InputFiles_v.at(i_f)).c_str());
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

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 10000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      osc_model.SetDeltaMSq23(c_osc::deltamsq23);
      double osc_weight = osc_model.NuMuSurvProb(nu_e); 
      osc_model.SetDeltaMSq23(c_osc::deltamsq23+3*c_osc::deltamsq23_sigma);
      double osc_weight_deltam2_plus = osc_model.NuMuSurvProb(nu_e); 
      osc_model.SetDeltaMSq23(c_osc::deltamsq23-3*c_osc::deltamsq23_sigma);
      double osc_weight_deltam2_minus = osc_model.NuMuSurvProb(nu_e); 

      if(generator != "GiBUU") weight = 1.0;
      weight *= scale*1e38*40;

      //std::cout << weight << std::endl;

      if(nu_pdg != 14 || ccnc != 1) continue;

      double W = CalcW(pdg,p4);
      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);
      std::vector<TVector3> neutron_mom = GetNeutronMom(pdg,p4);

      if(nprot < 1) continue;

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,neutron_mom,i_e);
        h_RecoEnergy_DeltaMSqCV.back().at(i_e)->Fill(nu_e_reco,weight*osc_weight);
        h_RecoEnergy_DeltaMSq_Plus.back().at(i_e)->Fill(nu_e_reco,weight*osc_weight_deltam2_plus);
        h_RecoEnergy_DeltaMSq_Minus.back().at(i_e)->Fill(nu_e_reco,weight*osc_weight_deltam2_minus);
      }
    }

  }

  TFile* f_out = new TFile("rootfiles/NuMuRatesDeltaMSq.root","RECREATE");

  // Reco energy plots scaled to 1 KT x 1e21 POT of exposure
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      for(int i=1;i<h_RecoEnergy_DeltaMSqCV.at(i_f).at(i_e)->GetNbinsX()+1;i++){
        h_RecoEnergy_DeltaMSqCV.at(i_f).at(i_e)->SetBinContent(i,Rate(total_flux,h_RecoEnergy_DeltaMSqCV.at(i_f).at(i_e)->GetBinContent(i)));
        h_RecoEnergy_DeltaMSq_Plus.at(i_f).at(i_e)->SetBinContent(i,Rate(total_flux,h_RecoEnergy_DeltaMSq_Plus.at(i_f).at(i_e)->GetBinContent(i)));
        h_RecoEnergy_DeltaMSq_Minus.at(i_f).at(i_e)->SetBinContent(i,Rate(total_flux,h_RecoEnergy_DeltaMSq_Minus.at(i_f).at(i_e)->GetBinContent(i)));
      }
      h_RecoEnergy_DeltaMSqCV.at(i_f).at(i_e)->Write();
      h_RecoEnergy_DeltaMSq_Plus.at(i_f).at(i_e)->Write();
      h_RecoEnergy_DeltaMSq_Minus.at(i_f).at(i_e)->Write();
    }
  }   

  f_out->Close();

}

