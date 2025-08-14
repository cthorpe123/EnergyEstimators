#include "Funcs/Funcs.h"
#include "Funcs/EnergyEstimatorFuncs.h"
#include "Funcs/Smearing.h"
#include "TLorentzVector.h"

bool nue_mode = false;

void CalculateResponseMatrices(){

  std::vector<std::string> InputFiles_v;
  if(!nue_mode) InputFiles_v = {"GENIEEventsFiltered.root","NuWroEventsFiltered.root","NEUTEventsFiltered.root","GiBUUEventsFiltered.root"};
  else InputFiles_v = {"GENIE_NueEventsFiltered.root","NuWro_NueEventsFiltered.root","NEUT_NueEventsFiltered.root","GiBUU_NueEventsFiltered.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy(kMAX);
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_Smeared(kMAX);

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);    

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
    int nprot;
    double W;
    std::vector<double>* est_nu_e=0;

    t->SetBranchAddress("scale",&scale);
    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);
    t->SetBranchAddress("W",&W);
    t->SetBranchAddress("nprot",&nprot);
    t->SetBranchAddress("est_nu_e",&est_nu_e);

    for(int i_e=0;i_e<kMAX;i_e++){
      std::string est = estimators_str.at(i_e);
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH2D((gen+"_TrueEnergy_RecoEnergy_"+est).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);",200,0.2,8.0,200,0.2,8.0));
      h_TrueEnergy_RecoEnergy_Smeared.at(i_e).push_back(new TH2D((gen+"_TrueEnergy_RecoEnergy_Smeared_"+est).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);",200,0.2,8.0,200,0.2,8.0));
    }

    int target_nu_pdg = 14;
    if(nue_mode) target_nu_pdg = 12;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 100000) break;
      if(ievent % 20000 == 0) std::cout << gen << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(gen != "GiBUU") weight = 1.0;
      weight *= scale*1e38*40;

      if(abs(nu_pdg) != target_nu_pdg || ccnc != 1) continue;
      if(nprot < 1) continue;

      std::vector<double> energies = GetEnergyEst(lepton_p4,pdg,p4);

      // Calculate predictions with kinematic smearing 
      smearing::smear_mom(*lepton_p4,13);
      for(int i=0;i<p4->size();i++) smearing::smear_mom(p4->at(i),pdg->at(i));
      std::vector<double> energies_smeared = GetEnergyEst(lepton_p4,pdg,p4);
      energies_smeared.at(kTotalEDep) = energies.at(kTotalEDep)*smearing::rng->Gaus(1.0,smearing::resolutions.at(0));

      for(int i_e=0;i_e<kMAX;i_e++){
        h_TrueEnergy_RecoEnergy.at(i_e).back()->Fill(nu_e,energies.at(i_e),weight);
        h_TrueEnergy_RecoEnergy_Smeared.at(i_e).back()->Fill(nu_e,energies_smeared.at(i_e),weight);
      }

    }

  }

  gSystem->Exec("mkdir -p Plots/ResponsePlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);
  std::string name;

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      std::string est = estimators_str.at(i_e);
      std::string gen = Generators_v.at(i_f);
      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f);
      TH2D* h_smeared = h_TrueEnergy_RecoEnergy_Smeared.at(i_e).at(i_f);

      Normalise(h);
      Normalise(h_smeared);

      h->Draw("colz");
      h->SetStats(0);
      h->SetMinimum(0.0);
      h->SetMaximum(0.8);
      if(!nue_mode) name = "Plots/ResponsePlots/NuMu_TrueEnergy_RecoEnergy_" + est + "_" + gen + ".png";
      else name =  "Plots/ResponsePlots/Nue_TrueEnergy_RecoEnergy_" + est + "_" + gen + ".png";

      h->SetContour(1000);
      c->Print(name.c_str()); 
      c->Clear();

      h_smeared->Draw("colz");
      h_smeared->SetStats(0);
      h->SetMinimum(0.0);
      h->SetMaximum(0.8);
      if(!nue_mode) name = "Plots/ResponsePlots/Smeared_NuMu_TrueEnergy_RecoEnergy_" + est + "_" + gen + ".png";
      else name =  "Plots/ResponsePlots/Smeared_Nue_TrueEnergy_RecoEnergy_" + est + "_" + gen + ".png";

      h_smeared->SetContour(1000);
      c->Print(name.c_str()); 
      c->Clear();

    }
  }

  TFile* f_out = nue_mode ? new TFile("rootfiles/NueResponseMatrices.root","RECREATE") : new TFile("rootfiles/NuMuResponseMatrices.root","RECREATE");
  for(size_t i_e=0;i_e<estimators_str.size();i_e++)
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++)
      h_TrueEnergy_RecoEnergy.at(i_e).at(i_f)->Write();

  f_out->Close(); 

}
