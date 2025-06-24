#include "Funcs/Funcs.h"
#include "TLorentzVector.h"

void CalculateRates(){

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();

  // h_flux is the flux in nu/m2/GeV/POT, this is the total flux nu/cm2/10^21 POT
  double total_flux = h_flux->Integral("width")*1e21/1e4; 
  std::cout << "Total Flux: " << total_flux << " nu/cm^2/10^21 POT" << std::endl; 

  std::vector<std::string> InputFiles_v = {"rootfiles/GENIEEvents.root","rootfiles/NuWroEvents.root","rootfiles/NEUTEvents.root","rootfiles/GiBUUEvents.root"};
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
    h_RecoEnergy.push_back(std::vector<TH1D*>());

    for(std::string estimator : estimators_str){
      h_TrueEnergy_RecoEnergy.back().push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str(),"Events/KT/GeV^{2}/10^{21} POT;True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);",100,0.1,8.0,100,0.1,8.0));
      h_RecoEnergy.back().push_back(new TH1D((generator+"_RecoEnergy_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/10^{21} POT",100,0.1,8.0));
    }

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      if(ievent > 50000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != 14 || ccnc != 1) continue;

      double W = CalcW(pdg,p4);
      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);

      if(generator != "GiBUU") weight = 1;
      weight *= scale*1e38*40;

      if(nprot < 1) continue;

      xsec.at(i_f) += weight;

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,i_e);
        h_TrueEnergy_RecoEnergy.back().at(i_e)->Fill(nu_e,nu_e_reco,weight);
        h_RecoEnergy.back().at(i_e)->Fill(nu_e_reco,weight);
      }



    }

    for(size_t i_e=0;i_e<h_TrueEnergy_RecoEnergy.back().size();i_e++){
      DivideByBinWidth(h_RecoEnergy.back().at(i_e)); 
      DivideByBinWidth2D(h_TrueEnergy_RecoEnergy.back().at(i_e)); 
    }

  }

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++) std::cout << Generators_v.at(i_f) << "  Cross Section = " << xsec.at(i_f) << " 1e-38 cm^2  Rate = " << Rate(total_flux,xsec.at(i_f)) << " events/KT/10^21 POT" << std::endl;


  gSystem->Exec("mkdir -p Plots/RatePlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      for(int i=1;i<h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetNbinsX()+1;i++)
        for(int j=1;j<h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetNbinsY()+1;j++)
          h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->SetBinContent(i,j,Rate(total_flux,h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetBinContent(i,j)));
     
      h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->Draw("colz");
      h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->SetStats(0);
      c->Print(("Plots/RatePlots/Rate_TrueEnergy_RecoEnergy_" + estimators_str.at(i_e) + "_" + Generators_v.at(i_f) + ".png").c_str()); 
      c->Clear();
    }
  }
  

}
