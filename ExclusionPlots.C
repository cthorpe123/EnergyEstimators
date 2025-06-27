#include "Funcs/Funcs.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

// Detector exposure to normalise plots to in KT x 10^21 POT

double fid_mass = 10; // active mass in KT
double POT = 1; // POT in 10^21

void ExclusionPlots(){

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

  std::vector<std::string> InputFiles_v = {"rootfiles/GENIE_NueEvents.root","rootfiles/NuWro_NueEvents.root","rootfiles/NEUT_NueEvents.root","rootfiles/GiBUU_NueEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy;
  std::vector<std::vector<TH1D*>> h_RecoEnergy_NoDeltaCP;
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Plus;
  std::vector<std::vector<TH1D*>> h_RecoEnergy_DeltaCP_Minus;

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    h_TrueEnergy_RecoEnergy.push_back(std::vector<TH2D*>());
    h_RecoEnergy_NoDeltaCP.push_back(std::vector<TH1D*>());
    h_RecoEnergy_DeltaCP_Plus.push_back(std::vector<TH1D*>());
    h_RecoEnergy_DeltaCP_Minus.push_back(std::vector<TH1D*>());

    int nbins = 100;
    for(std::string estimator : estimators_str){
      h_TrueEnergy_RecoEnergy.back().push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);Events/KT/GeV^{2}/8^{21} POT",nbins,0.1,8.0,nbins,0.1,8.0));
      h_RecoEnergy_NoDeltaCP.back().push_back(new TH1D((generator+"_RecoEnergy_NoDeltaCP_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",nbins,0.1,8.0));
      h_RecoEnergy_DeltaCP_Plus.back().push_back(new TH1D((generator+"_RecoEnergy_DeltaCP_Plus_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",nbins,0.1,8.0));
      h_RecoEnergy_DeltaCP_Minus.back().push_back(new TH1D((generator+"_RecoEnergy_DeltaCP_Minus_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",nbins,0.1,8.0));
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

      //if(ievent > 10000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      double osc_weight = nue_app_prob(nu_e); 
      double osc_weight_deltaCP_plus = nue_app_prob(nu_e,3.142/2); 
      double osc_weight_deltaCP_minus = nue_app_prob(nu_e,-3.142/2); 

      if(generator != "GiBUU") weight = 1.0;
      weight *= scale*1e38*40;

      //std::cout << weight << std::endl;

      if(nu_pdg != 12 || ccnc != 1) continue;

      double W = CalcW(pdg,p4);
      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);

      if(nprot < 1) continue;

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,i_e);
        h_TrueEnergy_RecoEnergy.back().at(i_e)->Fill(nu_e,nu_e_reco,weight);
        h_RecoEnergy_NoDeltaCP.back().at(i_e)->Fill(nu_e_reco,weight*osc_weight);
        h_RecoEnergy_DeltaCP_Plus.back().at(i_e)->Fill(nu_e_reco,weight*osc_weight_deltaCP_plus);
        h_RecoEnergy_DeltaCP_Minus.back().at(i_e)->Fill(nu_e_reco,weight*osc_weight_deltaCP_minus);
      }

    }

  }

  // Reco energy plots scaled to 1 KT x 1e21 POT of exposure
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      for(int i=1;i<h_RecoEnergy_NoDeltaCP.at(i_f).at(i_e)->GetNbinsX()+1;i++){
        h_RecoEnergy_NoDeltaCP.at(i_f).at(i_e)->SetBinContent(i,Rate(total_flux,h_RecoEnergy_NoDeltaCP.at(i_f).at(i_e)->GetBinContent(i)));
        h_RecoEnergy_DeltaCP_Plus.at(i_f).at(i_e)->SetBinContent(i,Rate(total_flux,h_RecoEnergy_DeltaCP_Plus.at(i_f).at(i_e)->GetBinContent(i)));
        h_RecoEnergy_DeltaCP_Minus.at(i_f).at(i_e)->SetBinContent(i,Rate(total_flux,h_RecoEnergy_DeltaCP_Minus.at(i_f).at(i_e)->GetBinContent(i)));
      }
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
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    std::vector<TGraph*> g_sqrtchi2_plus;
    std::vector<TGraph*> g_sqrtchi2_minus;

    TMultiGraph* mg_plus = new TMultiGraph();
    TMultiGraph* mg_minus = new TMultiGraph();

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::cout << estimators_str.at(i_e) << "  " << Generators_v.at(i_f) << std::endl;

      const TH1D* h_unscaled_zerodeltacp = h_RecoEnergy_NoDeltaCP.at(i_f).at(i_e); 
      const TH1D* h_unscaled_plusdeltacp = h_RecoEnergy_DeltaCP_Plus.at(i_f).at(i_e); 
      const TH1D* h_unscaled_minusdeltacp = h_RecoEnergy_DeltaCP_Minus.at(i_f).at(i_e); 

      std::vector<double> exposure_v;
      std::vector<double> sqrtchi2_plus_v,sqrtchi2_minus_v; 

      // Exposure in KT * 1e21 POT
      double exposure = 1.0; 
      while(exposure < 2000){

        TH1D* h_scaled_zerodeltacp = static_cast<TH1D*>(h_unscaled_zerodeltacp->Clone("h_scaled_zerodeltacp"));
        TH1D* h_scaled_plusdeltacp = static_cast<TH1D*>(h_unscaled_plusdeltacp->Clone("h_scaled_plusdeltacp"));
        TH1D* h_scaled_minusdeltacp = static_cast<TH1D*>(h_unscaled_minusdeltacp->Clone("h_scaled_minusdeltacp"));

        for(int i=1;i<h_scaled_zerodeltacp->GetNbinsX()+1;i++){
           h_scaled_zerodeltacp->SetBinContent(i,h_scaled_zerodeltacp->GetBinContent(i)*exposure);
           h_scaled_zerodeltacp->SetBinError(i,sqrt(h_scaled_zerodeltacp->GetBinContent(i)));
           h_scaled_plusdeltacp->SetBinContent(i,h_scaled_plusdeltacp->GetBinContent(i)*exposure);
           h_scaled_plusdeltacp->SetBinError(i,sqrt(h_scaled_plusdeltacp->GetBinContent(i)));
           h_scaled_minusdeltacp->SetBinContent(i,h_scaled_minusdeltacp->GetBinContent(i)*exposure);
           h_scaled_minusdeltacp->SetBinError(i,sqrt(h_scaled_minusdeltacp->GetBinContent(i)));
        }

        double chi2_plus=0.0,chi2_minus=0.0; 
        for(int i=1;i<h_scaled_zerodeltacp->GetNbinsX()+1;i++){
             chi2_plus += pow((h_scaled_plusdeltacp->GetBinContent(i) - h_scaled_zerodeltacp->GetBinContent(i))/h_scaled_plusdeltacp->GetBinError(i),2); 
             chi2_minus += pow((h_scaled_minusdeltacp->GetBinContent(i) - h_scaled_zerodeltacp->GetBinContent(i))/h_scaled_minusdeltacp->GetBinError(i),2); 
        }

        chi2_plus /= h_scaled_zerodeltacp->GetNbinsX();
        chi2_minus /= h_scaled_zerodeltacp->GetNbinsX();

        //std::cout << "exposure = " << exposure << "KT x 1e21 POT  sqrt(chi2_plus) = " << sqrt(chi2_plus) << "  sqrt(chi2_minus) = " << sqrt(chi2_minus) << std::endl;

        exposure_v.push_back(exposure);
        sqrtchi2_plus_v.push_back(sqrt(chi2_plus)); 
        sqrtchi2_minus_v.push_back(sqrt(chi2_minus)); 

        delete h_scaled_zerodeltacp;
        delete h_scaled_plusdeltacp;
        delete h_scaled_minusdeltacp;
        exposure += 1;
      }

        g_sqrtchi2_plus.push_back(new TGraph(exposure_v.size(),&(exposure_v[0]),&(sqrtchi2_plus_v[0])));
        g_sqrtchi2_minus.push_back(new TGraph(exposure_v.size(),&(exposure_v[0]),&(sqrtchi2_minus_v[0])));
        
        g_sqrtchi2_plus.back()->SetLineColor(i_f+2);
        g_sqrtchi2_plus.back()->SetLineWidth(2);
        g_sqrtchi2_minus.back()->SetLineColor(i_f+2);
        g_sqrtchi2_minus.back()->SetLineWidth(2);
        mg_plus->Add(g_sqrtchi2_plus.back());    
        mg_minus->Add(g_sqrtchi2_minus.back());    
        
        l->AddEntry(g_sqrtchi2_plus.back(),Generators_v.at(i_f).c_str(),"L");

    }

    mg_plus->Draw("AL");
    mg_plus->SetTitle(";Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
    l->Draw();
    c->Print(("Plots/ExclusionPlots/DeltaCPPlus_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();

    mg_minus->Draw("AL");
    mg_minus->SetTitle(";Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
    l->Draw();
    c->Print(("Plots/ExclusionPlots/DeltaCPMinus_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();
  
    l->Clear();

    int points = g_sqrtchi2_plus.at(0)->GetN();
    double low = g_sqrtchi2_plus.at(0)->GetX()[0];
    double high = g_sqrtchi2_plus.at(0)->GetX()[points-1];
    double width = (high-low)/points;        

    h_bands_plus.push_back(new TH1D(("h_bands_plus_"+estimators_str.at(i_e)).c_str(),"",points,low-width/2,high+width/2));
    h_lines_plus.push_back(new TH1D(("h_lines_plus_"+estimators_str.at(i_e)).c_str(),"",points,low-width/2,high+width/2));
    h_bands_minus.push_back(new TH1D(("h_bands_minus_"+estimators_str.at(i_e)).c_str(),"",points,low-width/2,high+width/2));
    h_lines_minus.push_back(new TH1D(("h_lines_minus_"+estimators_str.at(i_e)).c_str(),"",points,low-width/2,high+width/2));
   
    for(int i=0;i<points;i++){
        double min=1e10,max=-1e10;
        for(size_t i_f=0;i_f<g_sqrtchi2_plus.size();i_f++){
         min = std::min(min,g_sqrtchi2_plus.at(i_f)->GetY()[i]);
         max = std::max(max,g_sqrtchi2_plus.at(i_f)->GetY()[i]);
        }
        double center = (max + min)/2; 
        double width = (max - min)/2;
        h_bands_plus.back()->SetBinContent(i,center);
        h_bands_plus.back()->SetBinError(i,width);
        h_lines_plus.back()->SetBinContent(i,center);
    }
   
    for(int i=0;i<points;i++){
        double min=1e10,max=-1e10;
        for(size_t i_f=0;i_f<g_sqrtchi2_minus.size();i_f++){
         min = std::min(min,g_sqrtchi2_minus.at(i_f)->GetY()[i]);
         max = std::max(max,g_sqrtchi2_minus.at(i_f)->GetY()[i]);
        }
        double center = (max + min)/2; 
        double width = (max - min)/2;
        h_bands_minus.back()->SetBinContent(i,center);
        h_bands_minus.back()->SetBinError(i,width);
        h_lines_minus.back()->SetBinContent(i,center);
    }

  } 
      
  THStack* hs_bands_plus = new THStack("hs_bands_plus","#delta_{CP}=0 Exclusion at #delta_{CP}=#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
  THStack* hs_lines_plus = new THStack("hs_lines_plus","#delta_{CP}=0 Exclusion at #delta_{CP}=#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
  THStack* hs_bands_minus = new THStack("hs_bands_minus","#delta_{CP}=0 Exclusion at #delta_{CP}=-#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");
  THStack* hs_lines_minus = new THStack("hs_lines_minus","#delta_{CP}=0 Exclusion at #delta_{CP}=-#pi/2;Exposure (KT x 10^{21} POT);sqrt(#chi^{2}/ndof)");

  for(size_t i_e=0;i_e<h_bands_plus.size();i_e++){
    h_bands_plus.at(i_e)->SetFillColor(i_e+2);
    h_bands_plus.at(i_e)->SetFillStyle(3144);
    h_lines_plus.at(i_e)->SetLineColor(i_e+2);
    h_lines_plus.at(i_e)->SetLineWidth(2);
    hs_bands_plus->Add(h_bands_plus.at(i_e));
    hs_lines_plus->Add(h_lines_plus.at(i_e));
    l->AddEntry(h_lines_plus.at(i_e),estimators_str.at(i_e).c_str(),"L");
  }

   hs_bands_plus->Draw("nostack e4");
   hs_lines_plus->Draw("nostack HIST same");
   l->Draw();
   c->Print("Plots/ExclusionPlots/PlusDeltaCPZero_Exclusion.png"); 
   c->Clear();
   l->Clear();

  for(size_t i_e=0;i_e<h_bands_minus.size();i_e++){
    h_bands_minus.at(i_e)->SetFillColor(i_e+2);
    h_bands_minus.at(i_e)->SetFillStyle(3144);
    h_lines_minus.at(i_e)->SetLineColor(i_e+2);
    h_lines_minus.at(i_e)->SetLineWidth(2);
    hs_bands_minus->Add(h_bands_minus.at(i_e));
    hs_lines_minus->Add(h_lines_minus.at(i_e));
    l->AddEntry(h_lines_minus.at(i_e),estimators_str.at(i_e).c_str(),"L");
  }

   hs_bands_minus->Draw("nostack e4");
   hs_lines_minus->Draw("nostack HIST same");
   l->Draw();
   c->Print("Plots/ExclusionPlots/MinusDeltaCPZero_Exclusion.png"); 
   c->Clear();
   l->Clear();



}

