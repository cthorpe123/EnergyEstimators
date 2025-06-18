#include "Funcs/Funcs.h"
#include "TLorentzVector.h"

double osc_strength = 0.0;

void EnergyEstimatorPlots(){

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("Flux/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEFD_fastmc.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  TH1D* h_flux_osc = static_cast<TH1D*>(f_flux->Get("numu_fluxosc"));
  h_flux->SetDirectory(0);
  h_flux_osc->SetDirectory(0);
  f_flux->Close();

  // Ratio of oscillated and unoscillated fluxes
  h_flux_osc->Divide(h_flux);

  std::vector<std::string> InputFiles_v = {"rootfiles/GENIEEvents.root","rootfiles/NuWroEvents.root","rootfiles/NEUTEvents.root","rootfiles/GiBUUEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy;
  std::vector<std::vector<TH1D*>> h_RecoEnergy;
  std::vector<std::vector<TH1D*>> h_RecoEnergy_Osc;

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    h_TrueEnergy_RecoEnergy.push_back(std::vector<TH2D*>());
    h_RecoEnergy.push_back(std::vector<TH1D*>());
    h_RecoEnergy_Osc.push_back(std::vector<TH1D*>());
    for(std::string estimator : estimators_str){
      h_TrueEnergy_RecoEnergy.back().push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV)",50,0.1,7.0,50,0.1,7.0));
      h_RecoEnergy.back().push_back(new TH1D((generator+"_RecoEnergy_"+estimator).c_str(),";Estimated Neutrino Energy (GeV)",50,0.1,7.0));
      h_RecoEnergy_Osc.back().push_back(new TH1D((generator+"_RecoEnergy_Osc_"+estimator).c_str(),";Estimated Neutrino Energy (GeV)",50,0.1,7.0));
    }

    TFile* f = TFile::Open(InputFiles_v.at(i_f).c_str()) ;
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

      if(ievent > 10000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      double osc_weight = h_flux_osc->GetBinContent(h_flux_osc->FindBin(nu_e));

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != 14 || ccnc != 1) continue;

      double W = CalcW(pdg,p4);
      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);

      if(nprot < 1) continue;

      double Erec_MuonKin = T2KEnergy(lepton_p4);
      double Erec_MuonKinW = T2KEnergyW(lepton_p4,W);
      double Erec_MuonKinWNP = ubooneEnergy(lepton_p4,W,nprot);
      double Erec_PeLEELike = peleeEnergy(lepton_p4,proton_mom);
      double Erec_PeLEELike0Pi = peleeEnergy(lepton_p4,proton_mom);
      double Erec_TotalEDep = totaledepEnergy(lepton_p4,proton_mom,pion_mom,pizero_mom);

      h_TrueEnergy_RecoEnergy.back().at(kMuonKin)->Fill(nu_e,Erec_MuonKin,weight);
      h_TrueEnergy_RecoEnergy.back().at(kMuonKinW)->Fill(nu_e,Erec_MuonKinW,weight);
      h_TrueEnergy_RecoEnergy.back().at(kMuonKinWNP)->Fill(nu_e,Erec_MuonKinWNP,weight);
      h_TrueEnergy_RecoEnergy.back().at(kPeLEELike)->Fill(nu_e,Erec_PeLEELike,weight);
      if(!pion_mom.size() && !pizero_mom.size()) h_TrueEnergy_RecoEnergy.back().at(kPeLEELike0Pi)->Fill(nu_e,Erec_PeLEELike,weight);
      h_TrueEnergy_RecoEnergy.back().at(kTotalEDep)->Fill(nu_e,Erec_TotalEDep,weight);

      h_RecoEnergy.back().at(kMuonKin)->Fill(Erec_MuonKin,weight);
      h_RecoEnergy.back().at(kMuonKinW)->Fill(Erec_MuonKinW,weight);
      h_RecoEnergy.back().at(kMuonKinWNP)->Fill(Erec_MuonKinWNP,weight);
      h_RecoEnergy.back().at(kPeLEELike)->Fill(Erec_PeLEELike,weight);
      if(!pion_mom.size() && !pizero_mom.size()) h_RecoEnergy.back().at(kPeLEELike0Pi)->Fill(Erec_PeLEELike0Pi,weight);
      h_RecoEnergy.back().at(kTotalEDep)->Fill(Erec_TotalEDep,weight);

      h_RecoEnergy_Osc.back().at(kMuonKin)->Fill(Erec_MuonKin,weight*osc_weight);
      h_RecoEnergy_Osc.back().at(kMuonKinW)->Fill(Erec_MuonKinW,weight*osc_weight);
      h_RecoEnergy_Osc.back().at(kMuonKinWNP)->Fill(Erec_MuonKinWNP,weight*osc_weight);
      h_RecoEnergy_Osc.back().at(kPeLEELike)->Fill(Erec_PeLEELike,weight*osc_weight);
      if(!pion_mom.size() && !pizero_mom.size()) h_RecoEnergy_Osc.back().at(kPeLEELike0Pi)->Fill(Erec_PeLEELike0Pi,weight*osc_weight);
      h_RecoEnergy_Osc.back().at(kTotalEDep)->Fill(Erec_TotalEDep,weight*osc_weight);

    }

  }

  gSystem->Exec("mkdir -p Plots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      Normalise(h_TrueEnergy_RecoEnergy.at(i_f).at(i_e));
      h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->Draw("colz");
      h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->SetStats(0);
      c->Print(("Plots/TrueEnergy_RecoEnergy_" + estimators_str.at(i_e) + "_" + Generators_v.at(i_f) + ".png").c_str()); 
      c->Clear();
    }
  }

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    THStack* hs_Bias = new THStack("hs_Bias",";True Neutrino Energy (GeV);Frac. Bias");
    THStack* hs_Variance = new THStack("hs_Variance",";True Neutrino Energy (GeV);Frac. Variance");
    std::vector<TH1D*> h_Bias; 
    std::vector<TH1D*> h_Variance; 

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      h_Bias.push_back(nullptr);         
      h_Variance.push_back(nullptr);         
      GetBiasVariance(h_TrueEnergy_RecoEnergy.at(i_f).at(i_e),h_Bias.back(),h_Variance.back()); 
      l->AddEntry(h_Bias.back(),Generators_v.at(i_f).c_str(),"L");

      h_Bias.back()->SetLineWidth(2);
      h_Bias.back()->SetLineColor(i_f+1);
      hs_Bias->Add(h_Bias.back());  

      h_Variance.back()->SetLineWidth(2);
      h_Variance.back()->SetLineColor(i_f+1);
      hs_Variance->Add(h_Variance.back());  

    }

    hs_Bias->Draw("HIST nostack");
    hs_Bias->SetMaximum(0.8);
    hs_Bias->SetMinimum(-0.5);
    gPad->Modified();
    l->Draw();
    c->Print(("Plots/Bias_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();

    hs_Variance->Draw("HIST nostack");
    //hs_Variance->SetMaximum(0.25);
    //hs_Variance->SetMinimum(0.0);
    gPad->Modified();
    l->Draw();
    c->Print(("Plots/Variance_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();

    l->Clear();

    delete hs_Bias;
    delete hs_Variance;

  }

  // Oscillation plots

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    THStack* hs = new THStack("hs",";Estimated Neutrino Energy (GeV);Events");

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      double scale = 1.0/h_RecoEnergy.at(i_f).at(i_e)->Integral();
      h_RecoEnergy.at(i_f).at(i_e)->Scale(scale);
      h_RecoEnergy_Osc.at(i_f).at(i_e)->Scale(scale);

      h_RecoEnergy.at(i_f).at(i_e)->SetLineColor(i_f+1);
      h_RecoEnergy_Osc.at(i_f).at(i_e)->SetLineColor(i_f+1);
      h_RecoEnergy.at(i_f).at(i_e)->SetLineWidth(2);
      h_RecoEnergy_Osc.at(i_f).at(i_e)->SetLineWidth(2);
      h_RecoEnergy_Osc.at(i_f).at(i_e)->SetLineStyle(2);

      hs->Add(h_RecoEnergy.at(i_f).at(i_e));
      hs->Add(h_RecoEnergy_Osc.at(i_f).at(i_e));
      l->AddEntry(h_RecoEnergy.at(i_f).at(i_e),Generators_v.at(i_f).c_str(),"L");

    }

    hs->Draw("HIST nostack");
    gPad->Modified();
    l->Draw();
    c->Print(("Plots/Osc_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();
    l->Clear();

    delete hs;

  }


  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    THStack* hs_Ratio = new THStack("hs_Ratio",";Estimated Neutrino Energy (GeV);Osc/No Osc");

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      h_RecoEnergy_Osc.at(i_f).at(i_e)->Divide(h_RecoEnergy.at(i_f).at(i_e));

      h_RecoEnergy_Osc.at(i_f).at(i_e)->SetLineColor(i_f+1);
      h_RecoEnergy_Osc.at(i_f).at(i_e)->SetLineWidth(2);
      h_RecoEnergy_Osc.at(i_f).at(i_e)->SetLineStyle(1);

      hs_Ratio->Add(h_RecoEnergy_Osc.at(i_f).at(i_e));
      l->AddEntry(h_RecoEnergy_Osc.at(i_f).at(i_e),Generators_v.at(i_f).c_str(),"L");

    }

    hs_Ratio->Draw("HIST nostack");
    hs_Ratio->SetMinimum(0.0);
    hs_Ratio->SetMaximum(1.0);
    gPad->Modified();
    l->Draw();
    c->Print(("Plots/Osc_Ratio_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();
    l->Clear();

    delete hs_Ratio;

  }

  // Build covariance matrices 
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    std::vector<TH1D*> h_univ;
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++)
      h_univ.push_back(h_RecoEnergy.at(i_f).at(i_e));
        
    TH2D* h_Cov = MakeCovariance(h_univ);     
    h_Cov->SetName(("Covariance"+estimators_str.at(i_e)).c_str());

  } 

  

}
