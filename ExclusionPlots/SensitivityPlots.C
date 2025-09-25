#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"
#include "../Funcs/PlotSetup.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

//const int markers[4] = {8,21,22,34};
const int markers[4] = {4,25,26,28};

// Detector exposure to normalise plots to in KT x 10^21 POT
double fid_mass = 10; // active mass in KT
double POT = 1; // POT in 10^21

double scale_numu = 1;
double scale_nue = 1;

void SensitivityPlots(){

  PlotSetup();

  bool make_smear_plots = true;

  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_CV(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_Plus(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_Minus(kMAX);
  std::vector<std::vector<TGraph*>> g_chi2_deltam2_plus(kMAX);
  std::vector<std::vector<TGraph*>> g_chi2_deltam2_minus(kMAX);

  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_CV(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Plus(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Minus(kMAX);
  std::vector<std::vector<TGraph*>> g_chi2_deltacp_plus(kMAX);
  std::vector<std::vector<TGraph*>> g_chi2_deltacp_minus(kMAX);

  TH1D* h_axes = new TH1D("h_axes",";;#sqrt{X^{2}} Ratio",kMAX,-0.5,kMAX-0.5);

  // Load the histograms with predictions
  TFile* f_in_numu = TFile::Open("NuMuRatesDeltaMSq.root");
  TFile* f_in_nue = TFile::Open("NueRatesDeltaCP.root");

  // Reco energy plots scaled to 1 KT x 1e21 POT of exposure
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    std::string estimator = estimators_str.at(i_e);
    for(std::string generator : Generators_v){
      std::string name = "RecoEnergy";
      if(make_smear_plots) name += "_Smeared"; 
      h_RecoEnergy_DeltaMSq_CV.at(i_e).push_back(static_cast<TH1D*>(f_in_numu->Get((generator+"_"+name+"_DeltaMSq_CV_"+estimator).c_str())));
      h_RecoEnergy_DeltaMSq_Plus.at(i_e).push_back(static_cast<TH1D*>(f_in_numu->Get((generator+"_"+name+"_DeltaMSq_Plus_"+estimator).c_str())));
      h_RecoEnergy_DeltaMSq_Minus.at(i_e).push_back(static_cast<TH1D*>(f_in_numu->Get((generator+"_"+name+"_DeltaMSq_Minus_"+estimator).c_str())));
      h_RecoEnergy_DeltaCP_CV.at(i_e).push_back(static_cast<TH1D*>(f_in_nue->Get((generator+"_"+name+"_DeltaCP_CV_"+estimator).c_str())));
      h_RecoEnergy_DeltaCP_Plus.at(i_e).push_back(static_cast<TH1D*>(f_in_nue->Get((generator+"_"+name+"_DeltaCP_Plus_"+estimator).c_str())));
      h_RecoEnergy_DeltaCP_Minus.at(i_e).push_back(static_cast<TH1D*>(f_in_nue->Get((generator+"_"+name+"_DeltaCP_Minus_"+estimator).c_str())));
    }
  }  

  // Calculate chi2 w.r.t no delta CP as a function of exposure
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    std::string est = estimators_str.at(i_e);

    for(size_t i_f=0;i_f<Generators_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      std::cout << est << "  " << gen << std::endl;

      TH1D* h_deltacp_cv = h_RecoEnergy_DeltaCP_CV.at(i_e).at(i_f); 
      TH1D* h_deltacp_plus = h_RecoEnergy_DeltaCP_Plus.at(i_e).at(i_f); 
      TH1D* h_deltacp_minus = h_RecoEnergy_DeltaCP_Minus.at(i_e).at(i_f); 

      TH1D* h_deltam2_cv = h_RecoEnergy_DeltaMSq_CV.at(i_e).at(i_f); 
      TH1D* h_deltam2_plus = h_RecoEnergy_DeltaMSq_Plus.at(i_e).at(i_f); 
      TH1D* h_deltam2_minus = h_RecoEnergy_DeltaMSq_Minus.at(i_e).at(i_f); 

      // scale plots according to some exposure
      h_deltacp_cv->Scale(scale_nue);
      h_deltacp_plus->Scale(scale_nue);
      h_deltacp_minus->Scale(scale_nue);

      h_deltam2_cv->Scale(scale_numu);
      h_deltam2_plus->Scale(scale_numu);
      h_deltam2_minus->Scale(scale_numu);

      for(int i=1;i<h_deltacp_cv->GetNbinsX()+1;i++){
        h_deltacp_cv->SetBinError(i,sqrt(h_deltacp_cv->GetBinContent(i)));
        h_deltacp_plus->SetBinError(i,sqrt(h_deltacp_plus->GetBinContent(i)));
        h_deltacp_minus->SetBinError(i,sqrt(h_deltacp_minus->GetBinContent(i)));
      }

      for(int i=1;i<h_deltam2_cv->GetNbinsX()+1;i++){
        h_deltam2_cv->SetBinError(i,sqrt(h_deltam2_cv->GetBinContent(i)));
        h_deltam2_plus->SetBinError(i,sqrt(h_deltam2_plus->GetBinContent(i)));
        h_deltam2_minus->SetBinError(i,sqrt(h_deltam2_minus->GetBinContent(i)));
      }

      double deltam2_chi2_plus=0.0,deltam2_chi2_minus=0.0; 
      for(int i=1;i<h_deltam2_cv->GetNbinsX()+1;i++){
        deltam2_chi2_plus += pow((h_deltam2_plus->GetBinContent(i) - h_deltam2_cv->GetBinContent(i))/h_deltam2_plus->GetBinError(i),2); 
        deltam2_chi2_minus += pow((h_deltam2_minus->GetBinContent(i) - h_deltam2_cv->GetBinContent(i))/h_deltam2_minus->GetBinError(i),2); 
      }

      double deltacp_chi2_plus=0.0,deltacp_chi2_minus=0.0; 
      for(int i=1;i<h_deltacp_cv->GetNbinsX()+1;i++){
        deltacp_chi2_plus += pow((h_deltacp_plus->GetBinContent(i) - h_deltacp_cv->GetBinContent(i))/h_deltacp_plus->GetBinError(i),2); 
        deltacp_chi2_minus += pow((h_deltacp_minus->GetBinContent(i) - h_deltacp_cv->GetBinContent(i))/h_deltacp_minus->GetBinError(i),2); 
      }

      std::cout << "DeltaCP Plus: " << sqrt(deltacp_chi2_plus) << " Minus: " << sqrt(deltacp_chi2_minus) << std::endl;
      std::cout << "DeltaM2 Plus: " << sqrt(deltam2_chi2_plus) << " Minus: " << sqrt(deltam2_chi2_minus) << std::endl;

      std::vector<Double_t> x = {static_cast<Double_t>(i_e)};
      std::vector<Double_t> y_deltam2_plus = {sqrt(deltam2_chi2_plus)};
      std::vector<Double_t> y_deltam2_minus = {sqrt(deltam2_chi2_minus)};
      std::vector<Double_t> y_deltacp_plus = {sqrt(deltacp_chi2_plus)};
      std::vector<Double_t> y_deltacp_minus = {sqrt(deltacp_chi2_minus)};

      g_chi2_deltam2_plus.at(i_e).push_back(new TGraph(x.size(),&(x[0]),&(y_deltam2_plus[0]))); 
      g_chi2_deltam2_minus.at(i_e).push_back(new TGraph(x.size(),&(x[0]),&(y_deltam2_minus[0]))); 
      g_chi2_deltacp_plus.at(i_e).push_back(new TGraph(x.size(),&(x[0]),&(y_deltacp_plus[0]))); 
      g_chi2_deltacp_minus.at(i_e).push_back(new TGraph(x.size(),&(x[0]),&(y_deltacp_minus[0]))); 

      g_chi2_deltam2_plus.at(i_e).back()->SetMarkerColor(colors.at(i_e));
      g_chi2_deltam2_minus.at(i_e).back()->SetMarkerColor(colors.at(i_e));
      g_chi2_deltacp_plus.at(i_e).back()->SetMarkerColor(colors.at(i_e));
      g_chi2_deltacp_minus.at(i_e).back()->SetMarkerColor(colors.at(i_e));

      g_chi2_deltam2_plus.at(i_e).back()->SetMarkerStyle(markers[i_f]);
      g_chi2_deltam2_minus.at(i_e).back()->SetMarkerStyle(markers[i_f]);
      g_chi2_deltacp_plus.at(i_e).back()->SetMarkerStyle(markers[i_f]);
      g_chi2_deltacp_minus.at(i_e).back()->SetMarkerStyle(markers[i_f]);

      double size = 2.0;
      g_chi2_deltam2_plus.at(i_e).back()->SetMarkerSize(size);
      g_chi2_deltam2_minus.at(i_e).back()->SetMarkerSize(size);
      g_chi2_deltacp_plus.at(i_e).back()->SetMarkerSize(size);
      g_chi2_deltacp_minus.at(i_e).back()->SetMarkerSize(size);

    }

  }

  gSystem->Exec("mkdir -p Plots/");

  h_axes->Draw("HIST");
  SetAxisFontsH(h_axes); 
  h_axes->SetMinimum(0.45);
  h_axes->SetMaximum(1.25);
  h_axes->SetStats(0);
  for(int i_e=0;i_e<kMAX;i_e++) h_axes->GetXaxis()->SetBinLabel(i_e+1,estimators_leg.at(i_e).c_str());
  h_axes->GetXaxis()->SetLabelOffset(0.01);
  h_axes->GetXaxis()->SetLabelSize(0.065);

  std::vector<TGraph*> legs;
  for(size_t i_f=0;i_f<Generators_v.size();i_f++){
    std::string gen = Generators_v.at(i_f);
    legs.push_back(new TGraph());
    legs.back()->SetMarkerStyle(markers[i_f]);
    l->AddEntry(legs.back(),gen.c_str(),"P");
  }

  double ref = g_chi2_deltam2_plus.at(kMuonKinWNP).at(0)->GetY()[0];
  for(size_t i_f=0;i_f<Generators_v.size();i_f++){
    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      g_chi2_deltam2_plus.at(i_e).at(i_f)->GetY()[0] /= ref;
      g_chi2_deltam2_plus.at(i_e).at(i_f)->Draw("P");       
    }  
  }
  if(!make_smear_plots) c->Print("Plots/DeltaMSqPlusSen.pdf");
  else c->Print("Plots/SmearedDeltaMSqPlusSen.pdf");
  p_plot->Clear();

  ref = g_chi2_deltam2_minus.at(kMuonKinWNP).at(0)->GetY()[0];
  h_axes->Draw("HIST");
  SetAxisFontsH(h_axes); 
  h_axes->GetXaxis()->SetLabelOffset(0.01);
  h_axes->GetXaxis()->SetLabelSize(0.065);
  for(size_t i_f=0;i_f<Generators_v.size();i_f++){
    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      g_chi2_deltam2_minus.at(i_e).at(i_f)->GetY()[0] /= ref;
      g_chi2_deltam2_minus.at(i_e).at(i_f)->Draw("P");       
    }  
  }
  if(!make_smear_plots) c->Print("Plots/DeltaMSqMinusSen.pdf");
  else c->Print("Plots/SmearedDeltaMSqMinusSen.pdf");
  p_plot->Clear();

  ref = g_chi2_deltacp_plus.at(kMuonKinWNP).at(0)->GetY()[0];
  h_axes->Draw("HIST");
  SetAxisFontsH(h_axes); 
  h_axes->GetXaxis()->SetLabelOffset(0.01);
  h_axes->GetXaxis()->SetLabelSize(0.065);
  for(size_t i_f=0;i_f<Generators_v.size();i_f++){
    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      g_chi2_deltacp_plus.at(i_e).at(i_f)->GetY()[0] /= ref;
      g_chi2_deltacp_plus.at(i_e).at(i_f)->Draw("P");       
    }  
  }
  if(!make_smear_plots) c->Print("Plots/DeltaCPPlusSen.pdf");
  else c->Print("Plots/SmearedDeltaCPPlusSen.pdf");
  p_plot->Clear();

  h_axes->Draw("HIST");
  SetAxisFontsH(h_axes); 
  h_axes->GetXaxis()->SetLabelOffset(0.01);
  h_axes->GetXaxis()->SetLabelSize(0.065);
  ref = g_chi2_deltacp_minus.at(kMuonKinWNP).at(0)->GetY()[0];
  for(size_t i_f=0;i_f<Generators_v.size();i_f++){
    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      g_chi2_deltacp_minus.at(i_e).at(i_f)->GetY()[0] /= ref;
      g_chi2_deltacp_minus.at(i_e).at(i_f)->Draw("P");       
    }  
  }
  if(!make_smear_plots) c->Print("Plots/DeltaCPMinusSen.pdf");
  else c->Print("Plots/SmearedDeltaCPMinusSen.pdf");
  p_plot->Clear();

}

