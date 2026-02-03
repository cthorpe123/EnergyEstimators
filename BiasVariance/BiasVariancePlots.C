#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/Smearing.h"
#include "../Funcs/PlotSetup.h"
#include "TLorentzVector.h"

void BiasVariancePlots(){

  PlotSetup();

  bool draw_smeared = true;
  bool rebin = false;

  std::vector<std::string> InputFiles_v = {"GENIEEvents.root"/*,"NEUTEvents.root","GiBUUEvents.root","NuWroEvents.root"*/};
  std::vector<std::string> Generators_v = {"GENIE","NEUT","GiBUU","NuWro"};

  std::vector<double> true_binning_v;
  true_binning_v.push_back(0.5);
  for(int i=0;i<10;i++) true_binning_v.push_back(true_binning_v.back()+0.25);
  for(int i=0;i<4;i++) true_binning_v.push_back(true_binning_v.back()+0.5);
  true_binning_v.push_back(true_binning_v.back()+1.0);
  int true_nbins = true_binning_v.size()-1;
  double* true_binning_a = &true_binning_v[0];

  std::vector<double> reco_binning_v;
  reco_binning_v.push_back(0.0);
  reco_binning_v.push_back(0.1);
  for(int i=0;i<249;i++) reco_binning_v.push_back(reco_binning_v.back()+0.05);
  int reco_nbins = reco_binning_v.size()-1;
  double* reco_binning_a = &reco_binning_v[0];

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_Smeared(kMAX+1,std::vector<TH2D*>());

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    for(size_t i_e=0;i_e<kMAX;i_e++){
      std::string estimator = estimators_str.at(i_e);
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str(),"",true_nbins,true_binning_a,200,-1.0,10.0));
      h_TrueEnergy_RecoEnergy_Smeared.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_Smeared_"+estimator).c_str(),"",true_nbins,true_binning_a,200,-1.0,10.0));
    }
    h_TrueEnergy_RecoEnergy_Smeared.at(kMAX).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_Smeared2_TotalEDep").c_str(),"",true_nbins,true_binning_a,200,-1.0,10.0));


    TFile* f = TFile::Open(("/gluster/data/dune/cthorpe/DIS/"+InputFiles_v.at(i_f)).c_str());
    TTree* t = static_cast<TTree*>(f->Get("eventtree")) ;

    Double_t weight;
    Double_t nu_e;
    Int_t ccnc;
    Int_t nu_pdg;  
    Int_t lepton_pdg;
    TLorentzVector* lepton_p4=0;
    std::vector<int>* pdg=0;
    std::vector<TLorentzVector>* p4=0;

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 100000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;
      //if(generator == "GiBUU" && weight > 1) continue;

      if(nu_pdg != 14 || ccnc != 1) continue;

      if(nu_e > 8) continue;
      if(lepton_p4->Vect().Mag() < 0.1) continue; 

      int nprot = GetNProt(pdg,p4);
      double W = CalcW(pdg,p4);
      //if(nprot < 1) continue;

      std::vector<double> energies =  GetEnergyEst(lepton_p4,pdg,p4);

      // Calculate predictions with kinematic smearing 
      smearing::smear_mom(*lepton_p4,13);
      for(int i=0;i<p4->size();i++) smearing::smear_mom(p4->at(i),pdg->at(i));
      int nprot_smeared = GetNProt(pdg,p4);
      double W_smeared = CalcW(pdg,p4);
      std::vector<double> energies_smeared = GetEnergyEst(lepton_p4,pdg,p4);
      //energies_smeared.at(kTotalEDep) = energies.at(kTotalEDep)*smearing::rng->Gaus(1.0,smearing::resolutions.at(0));

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = energies.at(i_e);
        double nu_e_reco_smeared = energies_smeared.at(i_e);
        h_TrueEnergy_RecoEnergy.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight);  
        h_TrueEnergy_RecoEnergy_Smeared.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_smeared,weight);  
      }

      // Flat 20% smearing on energy
      h_TrueEnergy_RecoEnergy_Smeared.at(kMAX).at(i_f)->Fill(nu_e,energies.at(kTotalEDep)*smearing::rng->Gaus(1.0,smearing::resolutions.at(0)),weight);  

    }

  }// i_f


  gSystem->Exec("mkdir -p Plots/");

  // First calculate the 1D bias/variance plots
  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);
    THStack* hs_Bias = new THStack(("hs_Bias"+gen).c_str(),";True Neutrino Energy (GeV);Frac. Bias");
    THStack* hs_Variance = new THStack(("hs_Variance"+gen).c_str(),";True Neutrino Energy (GeV);Frac. Variance");
    std::vector<TH1D*> h_Bias;
    std::vector<TH1D*> h_Variance;
    std::vector<TH1D*> h_Bias_Smeared;
    std::vector<TH1D*> h_Variance_Smeared;

    for(size_t i_e=0;i_e<kMAX;i_e++){

      std::string est = estimators_str.at(i_e);

      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f);
      h_Bias.push_back(new TH1D(("h_Bias_"+gen+"_"+est).c_str(),"",true_nbins,true_binning_a));
      h_Variance.push_back(new TH1D(("h_Variance_"+gen+"_"+est).c_str(),"",true_nbins,true_binning_a));
      GetBiasVariance(h,h_Bias.back(),h_Variance.back()); 

      TH2D* h_smeared = h_TrueEnergy_RecoEnergy_Smeared.at(i_e).at(i_f);
      h_Bias_Smeared.push_back(new TH1D(("h_Bias_Smeared_"+gen+"_"+est).c_str(),"",true_nbins,true_binning_a));
      h_Variance_Smeared.push_back(new TH1D(("h_Variance_Smeared_"+gen+"_"+est).c_str(),"",true_nbins,true_binning_a));
      GetBiasVariance(h_smeared,h_Bias_Smeared.back(),h_Variance_Smeared.back()); 

      h_Bias.back()->SetLineColor(colors.at(i_e));
      h_Bias.back()->SetLineWidth(2);
      hs_Bias->Add(h_Bias.back());

      if(draw_smeared){
        h_Bias_Smeared.back()->SetLineColor(colors.at(i_e));
        h_Bias_Smeared.back()->SetLineWidth(2);
        //h_Bias_Smeared.back()->SetLineStyle(2);
        //hs_Bias->Add(h_Bias_Smeared.back());
      }

      h_Variance.back()->SetLineColor(colors.at(i_e));
      h_Variance.back()->SetLineWidth(2);
      //hs_Variance->Add(h_Variance.back());

      if(draw_smeared){
        h_Variance_Smeared.back()->SetLineColor(colors.at(i_e));
        h_Variance_Smeared.back()->SetLineWidth(2);
        //h_Variance_Smeared.back()->SetLineStyle(2);
        hs_Variance->Add(h_Variance_Smeared.back());
      }

      l->AddEntry(h_Bias.back(),estimators_leg.at(i_e).c_str(),"L");

      //delete h;
      delete h_smeared;

    } 

    // Make the variance plot with the 20% smearing
    h_Bias_Smeared.push_back(new TH1D(("h_Bias_Smeared2_"+gen+"_TotalEDep").c_str(),"",true_nbins,true_binning_a));
    h_Variance_Smeared.push_back(new TH1D(("h_Variance_Smeared2_"+gen+"_TotalEDep").c_str(),"",true_nbins,true_binning_a));
    GetBiasVariance(h_TrueEnergy_RecoEnergy_Smeared.at(kMAX).at(i_f),h_Bias_Smeared.back(),h_Variance_Smeared.back()); 
    h_Variance_Smeared.back()->SetLineColor(colors.at(kTotalEDep));
    h_Variance_Smeared.back()->SetLineWidth(2);
    h_Variance_Smeared.back()->SetLineStyle(2);
    hs_Variance->Add(h_Variance_Smeared.back());

    p_plot->cd();
    hs_Bias->Draw("nostack HIST");
    SetAxisFonts(hs_Bias);
    hs_Bias->SetMinimum(-0.2);
    hs_Bias->SetMaximum(0.02);
    c->Print(("Plots/Bias_Energy_"+gen+".pdf").c_str());  
    p_plot->Clear();

    p_plot->cd();
    hs_Variance->Draw("nostack HIST");
    SetAxisFonts(hs_Variance);
    hs_Variance->SetMinimum(0.0);
    hs_Variance->SetMaximum(0.19);
    c->Print(("Plots/Variance_Energy_"+gen+".pdf").c_str());  
    p_plot->Clear();

    l->Clear();

    delete hs_Bias;
    delete hs_Variance;

    // Calculate the increase in variance between the not smeared and smeared calculations
    // First calculate the 1D bias/variance plots

    THStack* hs_Variance2 = new THStack(("hs_Variance2"+gen).c_str(),";True Neutrino Energy (GeV);Increase in Frac. Variance");

    for(size_t i_e=0;i_e<kMAX;i_e++){
      std::string est = estimators_str.at(i_e);
      h_Variance_Smeared.at(i_e)->Add(h_Variance.at(i_e),-1); 
      hs_Variance2->Add(h_Variance_Smeared.at(i_e));
      l->AddEntry(h_Variance_Smeared.at(i_e),estimators_leg.at(i_e).c_str(),"L");
    }

    h_Variance_Smeared.at(kMAX)->Add(h_Variance.at(kTotalEDep),-1); 
    hs_Variance2->Add(h_Variance_Smeared.at(kMAX));
    //l->AddEntry(h_Variance_Smeared.back(),estimators_leg.at(i_e).c_str(),"L");

    p_plot->cd();
    hs_Variance2->Draw("nostack HIST");
    SetAxisFonts(hs_Variance);
    c->Print(("Plots/Variance_Energy_Increase_"+gen+".pdf").c_str());    
    p_plot->Clear();
    l->Clear();

    delete hs_Variance2;

  }

  // Compare bias and variance between generators
  for(size_t i_e=0;i_e<kMAX;i_e++){

    std::string est = estimators_str.at(i_e);

    THStack* hs_Bias = new THStack(("hs_Bias_"+est).c_str(),";True Neutrino Energy (GeV);Frac. Bias");
    THStack* hs_Variance = new THStack(("hs_Variance_"+est).c_str(),";True Neutrino Energy (GeV);Frac. Variance");
    std::vector<TH1D*> h_bias(InputFiles_v.size());
    std::vector<TH1D*> h_variance(InputFiles_v.size());

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f);
      h_bias.at(i_f) = new TH1D(("h_Bias2_"+gen+"_"+est).c_str(),"",true_nbins,true_binning_a);
      h_variance.at(i_f) = new TH1D(("h_Variance2_"+gen+"_"+est).c_str(),"",true_nbins,true_binning_a);
      GetBiasVariance(h,h_bias.at(i_f),h_variance.at(i_f)); 

      h_bias.at(i_f)->SetLineColor(colors.at(i_e));
      h_bias.at(i_f)->SetLineStyle(i_f+1);
      h_bias.at(i_f)->SetLineWidth(2);
      hs_Bias->Add(h_bias.at(i_f));
      l->AddEntry(h_bias.at(i_f),gen.c_str(),"L");

      h_variance.at(i_f)->SetLineColor(colors.at(i_e));
      h_variance.at(i_f)->SetLineStyle(i_f+1);
      h_variance.at(i_f)->SetLineWidth(2);
      hs_Variance->Add(h_variance.at(i_f));

    } 

    p_plot->cd();
    hs_Bias->Draw("nostack HIST");
    SetAxisFonts(hs_Bias);
    hs_Bias->SetMaximum(0.02);
    hs_Bias->SetMinimum(-0.5);
    c->Print(("Plots/Bias_Bands_"+est+".pdf").c_str()); 
    p_plot->Clear();

    p_plot->cd();
    hs_Variance->Draw("nostack HIST");
    SetAxisFonts(hs_Variance);
    hs_Variance->SetMaximum(0.21);
    //hs_Variance->SetMinimum(0.0);
    c->Print(("Plots/Variance_Bands_"+est+".pdf").c_str()); 
    p_plot->Clear();

    l->Clear();

  }

}
