#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

// Detector exposure to normalise plots to in KT x 10^21 POT
//double fid_mass = 10; // active mass in KT
//double POT = 1; // POT in 10^21

void DeltaCPExclusionPlots(){

  std::vector<std::string> estimators_tmp = { "MuonKin" , "MuonKinWNP" , "PeLEELike0Pi"  , "TotalEDep" , "SFMethod" };
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_CV(estimators_tmp.size());
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Plus(estimators_tmp.size());
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Minus(estimators_tmp.size());

  // Load the histograms with predictions
  TFile* f_in = TFile::Open("NueRatesDeltaCP.root");

  // Reco energy plots scaled to 1 KT x 1e21 POT of exposure
  for(size_t i_e=0;i_e<estimators_tmp.size();i_e++){
    std::string estimator = estimators_tmp.at(i_e);
    for(std::string generator : Generators_v){
      h_RecoEnergy_DeltaCP_CV.at(i_e).push_back(static_cast<TH1D*>(f_in->Get((generator+"_RecoEnergy_DeltaCP_CV_"+estimator).c_str())));
      h_RecoEnergy_DeltaCP_Plus.at(i_e).push_back(static_cast<TH1D*>(f_in->Get((generator+"_RecoEnergy_DeltaCP_Plus_"+estimator).c_str())));
      h_RecoEnergy_DeltaCP_Minus.at(i_e).push_back(static_cast<TH1D*>(f_in->Get((generator+"_RecoEnergy_DeltaCP_Minus_"+estimator).c_str())));
    }
  }   

  gSystem->Exec("mkdir -p Plots/ExclusionPlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  std::vector<TH1D*> h_bands_plus;
  std::vector<TH1D*> h_lines_plus;
  std::vector<TH1D*> h_bands_minus;
  std::vector<TH1D*> h_lines_minus;

  // Calculate chi2 w.r.t no delta CP as a function of exposure
  for(size_t i_e=0;i_e<estimators_tmp.size();i_e++){

    std::vector<TGraph*> g_sqrtchi2_plus;
    std::vector<TGraph*> g_sqrtchi2_minus;

    TMultiGraph* mg_plus = new TMultiGraph();
    TMultiGraph* mg_minus = new TMultiGraph();

    for(size_t i_f=0;i_f<Generators_v.size();i_f++){

      // Get the histograms with exposure of 1 KT * 10^21 POT
      const TH1D* h_unscaled_deltam2cv = h_RecoEnergy_DeltaCP_CV.at(i_e).at(i_f); 
      const TH1D* h_unscaled_plusdeltam2 = h_RecoEnergy_DeltaCP_Plus.at(i_e).at(i_f); 
      const TH1D* h_unscaled_minusdeltam2 = h_RecoEnergy_DeltaCP_Minus.at(i_e).at(i_f); 

      std::vector<double> exposure_v;
      std::vector<double> sqrtchi2_plus_v,sqrtchi2_minus_v; 

      // Exposure in KT * 1e21 POT
      double exposure = 1.0; 
      while(exposure < 500){

        TH1D* h_scaled_deltam2cv = static_cast<TH1D*>(h_unscaled_deltam2cv->Clone("h_scaled_deltam2cv"));
        TH1D* h_scaled_plusdeltam2 = static_cast<TH1D*>(h_unscaled_plusdeltam2->Clone("h_scaled_plusdeltam2"));
        TH1D* h_scaled_minusdeltam2 = static_cast<TH1D*>(h_unscaled_minusdeltam2->Clone("h_scaled_minusdeltam2"));

        // Set bin content and error to the exposure
        for(int i=1;i<h_scaled_deltam2cv->GetNbinsX()+1;i++){
          h_scaled_deltam2cv->SetBinContent(i,h_scaled_deltam2cv->GetBinContent(i)*exposure);
          h_scaled_deltam2cv->SetBinError(i,sqrt(h_scaled_deltam2cv->GetBinContent(i)));
          h_scaled_plusdeltam2->SetBinContent(i,h_scaled_plusdeltam2->GetBinContent(i)*exposure);
          h_scaled_plusdeltam2->SetBinError(i,sqrt(h_scaled_plusdeltam2->GetBinContent(i)));
          h_scaled_minusdeltam2->SetBinContent(i,h_scaled_minusdeltam2->GetBinContent(i)*exposure);
          h_scaled_minusdeltam2->SetBinError(i,sqrt(h_scaled_minusdeltam2->GetBinContent(i)));
        }

        // Calculate the chi2 between the zero deltaCP prediction and the two extreme deltaCP pred
        double chi2_plus=0.0,chi2_minus=0.0; 
        for(int i=1;i<h_scaled_deltam2cv->GetNbinsX()+1;i++){
          chi2_plus += pow((h_scaled_plusdeltam2->GetBinContent(i) - h_scaled_deltam2cv->GetBinContent(i))/h_scaled_plusdeltam2->GetBinError(i),2); 
          chi2_minus += pow((h_scaled_minusdeltam2->GetBinContent(i) - h_scaled_deltam2cv->GetBinContent(i))/h_scaled_minusdeltam2->GetBinError(i),2); 
        }

        chi2_plus /= h_scaled_deltam2cv->GetNbinsX();
        chi2_minus /= h_scaled_deltam2cv->GetNbinsX();

        exposure_v.push_back(exposure);
        sqrtchi2_plus_v.push_back(sqrt(chi2_plus)); 
        sqrtchi2_minus_v.push_back(sqrt(chi2_minus)); 

        delete h_scaled_deltam2cv;
        delete h_scaled_plusdeltam2;
        delete h_scaled_minusdeltam2;
        exposure += 1;
        
      }

      // Make graphs showing the exclusion curves for each generator/estimator combination
      g_sqrtchi2_plus.push_back(new TGraph(exposure_v.size(),&(exposure_v[0]),&(sqrtchi2_plus_v[0])));
      g_sqrtchi2_minus.push_back(new TGraph(exposure_v.size(),&(exposure_v[0]),&(sqrtchi2_minus_v[0])));

      g_sqrtchi2_plus.back()->SetLineStyle(i_f+1);
      g_sqrtchi2_plus.back()->SetLineColor(colors.at(i_e));
      g_sqrtchi2_plus.back()->SetLineWidth(2);
      g_sqrtchi2_minus.back()->SetLineColor(colors.at(i_e));
      g_sqrtchi2_minus.back()->SetLineStyle(i_f+1);
      g_sqrtchi2_minus.back()->SetLineWidth(2);
      mg_plus->Add(g_sqrtchi2_plus.back());    
      mg_minus->Add(g_sqrtchi2_minus.back());    

      l->AddEntry(g_sqrtchi2_plus.back(),Generators_v.at(i_f).c_str(),"L");

    }

    mg_plus->Draw("AL");
    mg_plus->SetTitle(";Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
    l->Draw();
    c->Print(("Plots/DeltaCPPlus_" + estimators_tmp.at(i_e) + ".pdf").c_str());
    c->Clear();

    mg_minus->Draw("AL");
    mg_minus->SetTitle(";Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
    l->Draw();
    c->Print(("Plots/DeltaCPMinus_" + estimators_tmp.at(i_e) + ".pdf").c_str());
    c->Clear();

    l->Clear();

    // Make band histograms
    int points = g_sqrtchi2_plus.at(0)->GetN();
    double low = g_sqrtchi2_plus.at(0)->GetX()[0];
    double high = g_sqrtchi2_plus.at(0)->GetX()[points-1];
    double width = (high-low)/points;        
    h_bands_plus.push_back(new TH1D(("h_bands_plus_"+estimators_tmp.at(i_e)).c_str(),"",points,low-width/2,high+width/2));
    h_lines_plus.push_back(new TH1D(("h_lines_plus_"+estimators_tmp.at(i_e)).c_str(),"",points,low-width/2,high+width/2));
    h_bands_minus.push_back(new TH1D(("h_bands_minus_"+estimators_tmp.at(i_e)).c_str(),"",points,low-width/2,high+width/2));
    h_lines_minus.push_back(new TH1D(("h_lines_minus_"+estimators_tmp.at(i_e)).c_str(),"",points,low-width/2,high+width/2));

    // Calculate the difference between the two extreme generator exclusion curves
    // DeltaCP Plus
    for(int i=0;i<points;i++){
      double min=1e10,max=-1e10;
      for(size_t i_f=0;i_f<g_sqrtchi2_plus.size();i_f++){
        min = std::min(min,g_sqrtchi2_plus.at(i_f)->GetY()[i]);
        max = std::max(max,g_sqrtchi2_plus.at(i_f)->GetY()[i]);
      }
      double center = (max + min)/2; 
      double width = (max - min)/2;
      h_bands_plus.back()->SetBinContent(i,center);
      h_lines_plus.back()->SetBinContent(i,center);
    }

    // DeltaCP Minus
    for(int i=0;i<points;i++){
      double min=1e10,max=-1e10;
      for(size_t i_f=0;i_f<g_sqrtchi2_minus.size();i_f++){
        min = std::min(min,g_sqrtchi2_minus.at(i_f)->GetY()[i]);
        max = std::max(max,g_sqrtchi2_minus.at(i_f)->GetY()[i]);
      }
      double center = (max + min)/2; 
      double width = (max - min)/2;
      h_bands_minus.back()->SetBinContent(i,center);
      h_lines_minus.back()->SetBinContent(i,center);
    }

  }  // i_e

  THStack* hs_bands_plus = new THStack("hs_bands_plus","#delta_{CP}=0 Exclusion at #delta_{CP}=#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
  THStack* hs_lines_plus = new THStack("hs_lines_plus","#delta_{CP}=0 Exclusion at #delta_{CP}=#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
  THStack* hs_bands_minus = new THStack("hs_bands_minus","#delta_{CP}=0 Exclusion at #delta_{CP}=-#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
  THStack* hs_lines_minus = new THStack("hs_lines_minus","#delta_{CP}=0 Exclusion at #delta_{CP}=-#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");

  for(size_t i_e=0;i_e<h_bands_plus.size();i_e++){
    //h_bands_plus.at(i_e)->SetFillColor(i_e+2);
    //h_bands_plus.at(i_e)->SetFillStyle(3344 + 110*i_e);
    h_lines_plus.at(i_e)->SetLineColor(colors.at(i_e));
    h_lines_plus.at(i_e)->SetLineWidth(3);
    hs_bands_plus->Add(h_bands_plus.at(i_e));
    hs_lines_plus->Add(h_lines_plus.at(i_e));
    l->AddEntry(h_lines_plus.at(i_e),estimators_tmp.at(i_e).c_str(),"L");
  }

  hs_bands_plus->Draw("nostack e4");
  //hs_lines_plus->Draw("nostack HIST same");
  for(size_t i_e=0;i_e<h_lines_plus.size();i_e++) h_lines_plus.at(i_e)->Draw("L same");
  l->Draw();
  c->Print("Plots/PlusDeltaCPZero_Exclusion.pdf"); 
  c->Clear();
  l->Clear();

  for(size_t i_e=0;i_e<h_bands_minus.size();i_e++){
    h_lines_minus.at(i_e)->SetLineColor(colors.at(i_e));
    h_lines_minus.at(i_e)->SetLineWidth(3);
    hs_bands_minus->Add(h_bands_minus.at(i_e));
    hs_lines_minus->Add(h_lines_minus.at(i_e));
    l->AddEntry(h_lines_minus.at(i_e),estimators_tmp.at(i_e).c_str(),"L");
  }

  hs_bands_minus->Draw("nostack e4");
  for(size_t i_e=0;i_e<h_lines_minus.size();i_e++) h_lines_minus.at(i_e)->Draw("L same");
  l->Draw();
  c->Print("Plots/MinusDeltaCPZero_Exclusion.pdf"); 
  c->Clear();
  l->Clear();
 
}

