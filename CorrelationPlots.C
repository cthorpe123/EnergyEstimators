#include "Funcs/Funcs.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

// Detector exposure to normalise plots to in KT x 10^21 POT
double exposure = 10*1; 

// Strength of oscillation effect
double osc_strength = 1.0;

void CorrelationPlots(){

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

  // Ratio of oscillated and unoscillated fluxes
  h_flux_osc->Divide(h_flux);

  TH1D* h_uniform = static_cast<TH1D*>(h_flux_osc->Clone("h_uniform"));
  h_uniform->Divide(h_uniform);

  std::vector<std::string> InputFiles_v = {"rootfiles/GENIEEvents.root","rootfiles/NuWroEvents.root","rootfiles/NEUTEvents.root","rootfiles/GiBUUEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH1D*>> h_RecoEnergy(kMAX,std::vector<TH1D*>());
  std::vector<std::vector<TH1D*>> h_RecoEnergy_Osc(kMAX,std::vector<TH1D*>());

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    for(int i_e=0;i_e<estimators_str.size();i_e++){
      std::string estimator = estimators_str.at(i_e);
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);Events/KT/GeV^{2}/8^{21} POT",100,0.1,8.0,100,0.1,8.0));
      h_RecoEnergy.at(i_e).push_back(new TH1D((generator+"_RecoEnergy_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",100,0.1,8.0));
      h_RecoEnergy_Osc.at(i_e).push_back(new TH1D((generator+"_RecoEnergy_Osc_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",100,0.1,8.0));
    }

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

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      double osc_weight = 1.0 - osc_strength*(1.0 - h_flux_osc->GetBinContent(h_flux_osc->FindBin(nu_e)));

      if(generator != "GiBUU") weight = 1.0;
      weight *= scale*1e38*40;

      //std::cout << weight << std::endl;

      if(nu_pdg != 14 || ccnc != 1) continue;

      double W = CalcW(pdg,p4);
      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);

      if(nprot < 1) continue;

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,i_e);
        h_TrueEnergy_RecoEnergy.at(i_e).back()->Fill(nu_e,nu_e_reco,weight);
        h_RecoEnergy.at(i_e).back()->Fill(nu_e_reco,weight);
        h_RecoEnergy_Osc.at(i_e).back()->Fill(nu_e_reco,weight*osc_weight);
      }
    }

  }


  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      for(int i=1;i<h_RecoEnergy.at(i_e).at(i_f)->GetNbinsX()+1;i++){
        h_RecoEnergy.at(i_e).at(i_f)->SetBinContent(i,exposure*Rate(total_flux,h_RecoEnergy.at(i_e).at(i_f)->GetBinContent(i)));
        h_RecoEnergy.at(i_e).at(i_f)->SetBinError(i,sqrt(h_RecoEnergy.at(i_e).at(i_f)->GetBinContent(i)));
        h_RecoEnergy_Osc.at(i_e).at(i_f)->SetBinContent(i,exposure*Rate(total_flux,h_RecoEnergy_Osc.at(i_e).at(i_f)->GetBinContent(i)));
        h_RecoEnergy_Osc.at(i_e).at(i_f)->SetBinError(i,sqrt(h_RecoEnergy_Osc.at(i_e).at(i_f)->GetBinContent(i)));
      }
    }
  }   

  gSystem->Exec("mkdir -p Plots/CorrelationPlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);


  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    std::string estimator = estimators_str.at(i_e);

    TMatrixDSym cov = MakeCovariance(h_RecoEnergy.at(i_e));
    TMatrixDSym cov_osc = MakeCovariance(h_RecoEnergy_Osc.at(i_e));

    int nbins = h_RecoEnergy.at(i_e).at(0)->GetNbinsX();
    int low = h_RecoEnergy.at(i_e).at(0)->GetBinLowEdge(1);
    int high = h_RecoEnergy.at(i_e).at(0)->GetBinLowEdge(nbins+1);

    TH2D* h_corr = new TH2D(("h_corr_"+estimator).c_str(),"Correlation;Estimated Neutrino Energy (GeV);Estimated Neutrino Energy (GeV)",nbins,low,high,nbins,low,high);
    TH2D* h_corr_osc = new TH2D(("h_corr_osc_"+estimator).c_str(),"Correlation;Estimated Neutrino Energy (GeV);Estimated Neutrino Energy (GeV)",nbins,low,high,nbins,low,high);

    for(int i=1;i<nbins+1;i++){
      for(int j=1;j<nbins+1;j++){
        h_corr->SetBinContent(i,j,cov[i-1][j-1]/sqrt(cov[i-1][i-1])/sqrt(cov[j-1][j-1]));
        h_corr_osc->SetBinContent(i,j,cov[i-1][j-1]/sqrt(cov[i-1][i-1])/sqrt(cov[j-1][j-1]));
        //std::cout << cov[i-1][j-1] << std::endl;
      }
    }

    h_corr->Draw("colz");
    h_corr->SetStats(0);
    c->Print(("Plots/CorrelationPlots/Correlation_"+estimator+".png").c_str());
    c->Clear(); 

    h_corr_osc->Draw("colz");
    h_corr_osc->SetStats(0);
    c->Print(("Plots/CorrelationPlots/Correlation_Osc_"+estimator+".png").c_str());
    c->Clear(); 

  } 

}

