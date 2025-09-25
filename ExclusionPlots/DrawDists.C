#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"
#include "../Funcs/PlotSetup.h"

double scale = 200;

void DrawDists(){

  PlotSetup();
  bool make_smear_plots = false;

  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_CV(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_Plus(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaMSq_Minus(kMAX);

  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_CV(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Plus(kMAX);
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Minus(kMAX);

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

  for(int i_e=0;i_e<kMAX;i_e++){
    for(size_t i_f=0;i_f<Generators_v.size();i_f++){
      std::string gen = Generators_v.at(i_f);    
      std::string est = estimators_str.at(i_e);

      TH1D* h_deltacp_cv = h_RecoEnergy_DeltaCP_CV.at(i_e).at(i_f); 
      TH1D* h_deltacp_plus = h_RecoEnergy_DeltaCP_Plus.at(i_e).at(i_f); 
      TH1D* h_deltacp_minus = h_RecoEnergy_DeltaCP_Minus.at(i_e).at(i_f); 

      TH1D* h_deltam2_cv = h_RecoEnergy_DeltaMSq_CV.at(i_e).at(i_f); 
      TH1D* h_deltam2_plus = h_RecoEnergy_DeltaMSq_Plus.at(i_e).at(i_f); 
      TH1D* h_deltam2_minus = h_RecoEnergy_DeltaMSq_Minus.at(i_e).at(i_f); 

      h_deltacp_cv->Scale(scale);
      h_deltacp_plus->Scale(scale);
      h_deltacp_minus->Scale(scale);

      h_deltam2_cv->Scale(scale);
      h_deltam2_plus->Scale(scale);
      h_deltam2_minus->Scale(scale);

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

      THStack* hs_deltam2 = new THStack("hs_deltam2",";E_{est} (GeV);Events");

      h_deltam2_cv->SetLineColor(colors.at(i_e));
      h_deltam2_cv->SetMarkerColor(colors.at(i_e));
      h_deltam2_cv->SetLineWidth(2);
      h_deltam2_cv->SetMarkerStyle(20);
      hs_deltam2->Add(h_deltam2_cv,"HIST");
      l->AddEntry(h_deltam2_cv,"CV","L");     

      h_deltam2_plus->SetLineColor(colors.at(i_e));
      h_deltam2_plus->SetMarkerColor(colors.at(i_e));
      h_deltam2_plus->SetMarkerStyle(21);
      hs_deltam2->Add(h_deltam2_plus,"e1");
      l->AddEntry(h_deltam2_plus,"+3#sigma","P");     

      h_deltam2_minus->SetLineColor(colors.at(i_e));
      h_deltam2_minus->SetMarkerColor(colors.at(i_e));
      h_deltam2_minus->SetMarkerStyle(22);
      hs_deltam2->Add(h_deltam2_minus,"e1");
      l->AddEntry(h_deltam2_minus,"+3#sigma","P");     

      hs_deltam2->Draw("nostack");
      SetAxisFonts(hs_deltam2);
      if(!make_smear_plots) c->Print(("Plots/NuMu_Dist_"+gen+"_"+est+".pdf").c_str());
      else c->Print(("Plots/NuMu_Dist_Smeared_"+gen+"_"+est+".pdf").c_str());
      p_plot->Clear();
      l->Clear();

      delete hs_deltam2;

      THStack* hs_deltacp = new THStack("hs_deltacp",";E_{est} (GeV);Events");

      h_deltacp_cv->SetLineColor(colors.at(i_e));
      h_deltacp_cv->SetMarkerColor(colors.at(i_e));
      h_deltacp_cv->SetLineWidth(2);
      hs_deltacp->Add(h_deltacp_cv,"HIST");
      l->AddEntry(h_deltacp_cv,"#delta_{CP} = 0","L");     

      h_deltacp_plus->SetLineColor(colors.at(i_e));
      h_deltacp_plus->SetMarkerColor(colors.at(i_e));
      h_deltacp_plus->SetMarkerStyle(21);
      hs_deltacp->Add(h_deltacp_plus,"e1");
      l->AddEntry(h_deltacp_plus,"#delta_{CP} = #pi/2","P");     

      h_deltacp_minus->SetLineColor(colors.at(i_e));
      h_deltacp_minus->SetMarkerColor(colors.at(i_e));
      h_deltacp_minus->SetMarkerStyle(22);
      hs_deltacp->Add(h_deltacp_minus,"e1");
      l->AddEntry(h_deltacp_minus,"#delta_{CP} = -#pi/2","P");     

      hs_deltacp->Draw("nostack");
      SetAxisFonts(hs_deltacp);
      if(!make_smear_plots) c->Print(("Plots/Nue_Dist_"+gen+"_"+est+".pdf").c_str());
      else c->Print(("Plots/Nue_Dist_Smeared_"+gen+"_"+est+".pdf").c_str());
      p_plot->Clear();
      l->Clear();

      delete hs_deltacp;

    }
  } 

}
