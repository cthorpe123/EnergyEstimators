#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "TLorentzVector.h"

void FSIStudy(){

  std::vector<std::string> InputFiles_v = {"GENIEEvents2.root","NuWroEvents2.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro"};

  std::vector<double> binning_v;
  binning_v.push_back(0.35);
  for(int i=0;i<31;i++) binning_v.push_back(binning_v.back()+0.15);
  for(int i=0;i<10;i++) binning_v.push_back(binning_v.back()+0.3);
  int nbins = binning_v.size()-1;
  double* binning_a = &binning_v[0];

  std::vector<std::vector<TH1D*>> h_EnergyBias(kMAX,std::vector<TH1D*>());
  std::vector<std::vector<TH1D*>> h_EnergyBias_NoFSI(kMAX,std::vector<TH1D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_NoFSI(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_Dense(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_NoFSI_Dense(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH1D*>> h_Change(kMAX,std::vector<TH1D*>());
  std::vector<std::vector<TH1D*>> h_Change_Abs(kMAX,std::vector<TH1D*>());

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

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
    std::vector<int>* pdg_nofsi=0;
    std::vector<TLorentzVector>* p4_nofsi=0;

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);
    t->SetBranchAddress("pdg_nofsi",&pdg_nofsi);
    t->SetBranchAddress("p4_nofsi",&p4_nofsi);

    for(size_t i_e=0;i_e<kMAX;i_e++){
      std::string est = estimators_str.at(i_e);
      h_EnergyBias.at(i_e).push_back(new TH1D((generator+"_EnergyBias_"+est).c_str(),"",100,-1,1));
      h_EnergyBias_NoFSI.at(i_e).push_back(new TH1D((generator+"_EnergyBias_NoFSI_"+est).c_str(),"",100,-1,1));
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+est).c_str(),"",nbins,binning_a,200,0.0,10.0));
      h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_NoFSI_"+est).c_str(),"",nbins,binning_a,200,0.0,10.0));
      h_TrueEnergy_RecoEnergy_Dense.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_Dense_"+est).c_str(),"",200,0.2,8.0,200,0.2,8.0));
      h_TrueEnergy_RecoEnergy_NoFSI_Dense.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_NoFSI_Dense_"+est).c_str(),"",200,0.2,8.0,200,0.2,8.0));
      h_Change.at(i_e).push_back(new TH1D((generator+"_Change_"+est).c_str(),"",50,-0.3,0.3));
      h_Change_Abs.at(i_e).push_back(new TH1D((generator+"_Change_Abs_"+est).c_str(),"",50,0.0,0.3));
    }

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 100000) break;
      if(ievent % 50000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != 14 || ccnc != 1) continue;

      std::vector<double> energies = GetEnergyEst(lepton_p4,pdg,p4);
      std::vector<double> energies_nofsi = GetEnergyEst(lepton_p4,pdg_nofsi,p4_nofsi);

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = energies.at(i_e);
        double nu_e_reco_nofsi = energies_nofsi.at(i_e);
        if(nu_e_reco > 0){
          h_EnergyBias.at(i_e).at(i_f)->Fill((nu_e_reco-nu_e)/nu_e,weight);  
          h_TrueEnergy_RecoEnergy.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight); 
          h_TrueEnergy_RecoEnergy_Dense.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight); 
        }
        if(nu_e_reco_nofsi > 0){
          h_EnergyBias_NoFSI.at(i_e).at(i_f)->Fill((nu_e_reco_nofsi-nu_e)/nu_e,weight);  
          h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight); 
          h_TrueEnergy_RecoEnergy_NoFSI_Dense.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight); 
        }

        if(nu_e_reco > 0 && nu_e_reco_nofsi > 0){
          h_Change.at(i_e).at(i_f)->Fill((nu_e_reco - nu_e_reco_nofsi)/nu_e_reco_nofsi,weight);
          h_Change_Abs.at(i_e).at(i_f)->Fill(abs(nu_e_reco - nu_e_reco_nofsi)/nu_e_reco_nofsi,weight);
        }

      }

    }

  }// i_f


  gSystem->Exec("mkdir -p Plots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string est = estimators_str.at(i_e);
      std::string gen = Generators_v.at(i_f);

      THStack* hs = new THStack("hs",";(E_{est} - E_{true})/E_{true};");

      TH1D* h = h_EnergyBias.at(i_e).at(i_f);
      h->SetLineWidth(2);
      h->SetLineColor(1);
      l->AddEntry(h,"With FSI","L");
      hs->Add(h);

      TH1D* h_nofsi = h_EnergyBias_NoFSI.at(i_e).at(i_f);
      h_nofsi->SetLineWidth(2);
      h_nofsi->SetLineColor(2);
      l->AddEntry(h_nofsi,"No FSI","L");
      hs->Add(h_nofsi);

      hs->Draw("HIST nostack");
      l->Draw();
      c->Print(("Plots/BiasShape_"+ est + "_" + gen + ".pdf").c_str());
      c->Clear();
      l->Clear();

      delete hs;

    }
  }

  // calculate bias and variance afo neutrino energy with and without FSI

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);
    THStack* hs_Bias = new THStack(("hs_Bias"+gen).c_str(),";True Neutrino Energy (GeV);Frac. Bias");
    THStack* hs_Variance = new THStack(("hs_Variance"+gen).c_str(),";True Neutrino Energy (GeV);Frac. Variance");
    std::vector<TH1D*> h_Bias;
    std::vector<TH1D*> h_Variance;

    for(size_t i_e=0;i_e<kMAX;i_e++){

      std::string est = estimators_str.at(i_e);

      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f); 
      h_Bias.push_back(new TH1D(("h_Bias_"+gen+"_"+est).c_str(),"",nbins,binning_a));
      h_Variance.push_back(new TH1D(("h_Variance_"+gen+"_"+est).c_str(),"",nbins,binning_a));
      GetBiasVariance(h,h_Bias.back(),h_Variance.back()); 

      h_Bias.back()->SetLineWidth(2);
      h_Bias.back()->SetLineColor(colors.at(i_e));
      hs_Bias->Add(h_Bias.back());
      l->AddEntry(h_Bias.back(),est.c_str(),"L"); 

      h_Variance.back()->SetLineWidth(2);
      h_Variance.back()->SetLineColor(colors.at(i_e));
      hs_Variance->Add(h_Variance.back());

      TH2D* h_nofsi = h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f); 
      h_Bias.push_back(new TH1D(("h_Bias_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a));
      h_Variance.push_back(new TH1D(("h_Variance_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a));
      GetBiasVariance(h_nofsi,h_Bias.back(),h_Variance.back()); 

      h_Bias.back()->SetLineWidth(2);
      h_Bias.back()->SetLineColor(colors.at(i_e));
      h_Bias.back()->SetLineStyle(2);
      hs_Bias->Add(h_Bias.back());

      h_Variance.back()->SetLineWidth(2);
      h_Variance.back()->SetLineColor(colors.at(i_e));
      h_Variance.back()->SetLineStyle(2);
      hs_Variance->Add(h_Variance.back());

    }

    hs_Bias->Draw("HIST nostack");
    l->Draw();
    c->Print(("Plots/Bias_"+gen+".pdf").c_str());
    c->Clear();

    hs_Variance->Draw("HIST nostack");
    l->Draw();
    c->Print(("Plots/Variance_"+gen+".pdf").c_str());
    c->Clear();

    l->Clear();
    delete hs_Bias;
    delete hs_Variance;       

    h_Bias.clear();
    h_Variance.clear();

  }

  // Calculate the change in the estimated energy from adding FSI
  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);
    THStack* hs_Change = new THStack(("hs_Change"+gen).c_str(),";(E_{est}^{FSI} - E_{est}^{No FSI})/E_{est}^{No FSI};");

    for(size_t i_e=0;i_e<kMAX;i_e++){

      if(i_e == kMuonKin) continue;

      std::string est = estimators_str.at(i_e);
      TH1D* h = h_Change.at(i_e).at(i_f); 
      h->Scale(1.0/h->Integral());
      h->SetLineWidth(2);
      h->SetLineColor(colors.at(i_e));
      hs_Change->Add(h);
      l->AddEntry(h,est.c_str(),"L"); 
    }

    hs_Change->Draw("nostack HIST");
    l->Draw();
    c->Print(("Plots/Change_"+gen+".pdf").c_str()); 
    c->Clear();
    l->Clear();

  }


  // Calculate the change in the estimated energy from adding FSI
  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);
    THStack* hs_Change_Abs = new THStack(("hs_Change_Abs"+gen).c_str(),";|E_{est}^{FSI} - E_{est}^{No FSI}|/E_{est}^{No FSI};");

    for(size_t i_e=0;i_e<kMAX;i_e++){

      if(i_e == kMuonKin) continue;

      std::string est = estimators_str.at(i_e);
      TH1D* h = h_Change_Abs.at(i_e).at(i_f); 
      h->Scale(1.0/h->Integral());
      h->SetLineWidth(2);
      h->SetLineColor(colors.at(i_e));
      hs_Change_Abs->Add(h);
      l->AddEntry(h,est.c_str(),"L"); 
    }

    hs_Change_Abs->Draw("nostack HIST");
    l->Draw();
    c->Print(("Plots/Change_Abs_"+gen+".pdf").c_str()); 
    c->Clear();
    l->Clear();

  }


  std::vector<std::vector<TH1D*>> h_Ratio(kMAX);

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);

    for(size_t i_e=0;i_e<kMAX;i_e++){

      if(i_e != kMuonKinWNP && i_e != kTotalEDep) continue;

      std::string est = estimators_str.at(i_e);

      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f); 
      TH1D* h_bias = new TH1D(("h_Bias_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      TH1D* h_variance = new TH1D(("h_Variance_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetBiasVariance(h,h_bias,h_variance); 


      TH2D* h_nofsi = h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f); 
      TH1D* h_bias_nofsi = new TH1D(("h_Bias_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      TH1D* h_variance_nofsi = new TH1D(("h_Variance_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetBiasVariance(h_nofsi,h_bias_nofsi,h_variance_nofsi); 

      h_Ratio.at(i_e).push_back(static_cast<TH1D*>(h_bias->Clone(("h_Ratio_"+est).c_str())));
      h_Ratio.at(i_e).back()->Add(h_bias_nofsi,-1);

    }

  }

  THStack* hs = new THStack("hs",";True Neutrino Energy (GeV);Frac. Bias No FSI - Frac. Bias FSI");

  for(size_t i_e=0;i_e<kMAX;i_e++){

    if(i_e != kMuonKinWNP && i_e != kTotalEDep) continue;

    std::string est = estimators_str.at(i_e);

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      std::string gen = Generators_v.at(i_f);
      h_Ratio.at(i_e).at(i_f)->SetLineWidth(2);
      h_Ratio.at(i_e).at(i_f)->SetLineColor(colors.at(i_e));
      h_Ratio.at(i_e).at(i_f)->SetLineStyle(i_f+1);
      hs->Add(h_Ratio.at(i_e).at(i_f));
      //l->AddEntry(h_Ratio.at(i_e).at(i_f),Generators_v.at(i_f).c_str(),"L");
      l->AddEntry(h_Ratio.at(i_e).at(i_f),(est+" "+gen).c_str(),"L");
    }

  }

  hs->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/Ratio.pdf"); 
  c->Clear();
  l->Clear();


  // Compare the energy spectra before and after applying FSI for the two generators

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("../Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();
  h_flux->Scale(1.0/h_flux->Integral());


  for(size_t i_e=0;i_e<kMAX;i_e++){
    std::string est = estimators_str.at(i_e);
    if(i_e != kMuonKinWNP && i_e != kTotalEDep) continue;

    THStack* hs = new THStack("hs",";E_{est} (GeV);Pred");
    std::vector<TH1D*> h_v(InputFiles_v.size());
    std::vector<TH1D*> h_nofsi_v(InputFiles_v.size());

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      TH2D* h = static_cast<TH2D*>(h_TrueEnergy_RecoEnergy_Dense.at(i_e).at(i_f)->Clone(("res_"+gen+"_"+est).c_str())); 
      TH2D* h_nofsi = static_cast<TH2D*>(h_TrueEnergy_RecoEnergy_NoFSI_Dense.at(i_e).at(i_f)->Clone(("res_nofsi_"+gen+"_"+est).c_str())); 

      Normalise(h);
      Normalise(h_nofsi); 

      // Fold the 2D plot into a 1D plot
      h_v.at(i_f) = new TH1D(("h_reco_"+gen+"_"+est).c_str(),";E_{est} (GeV);Events",200,0.2,8.0); 
      for(int j=1;j<h_v.at(i_f)->GetNbinsX()+1;j++){
        double events = 0.0;
        for(int i=1;i<h->GetNbinsX()+1;i++){
          double flux = h_flux->GetBinContent(h_flux->FindBin(h->GetXaxis()->GetBinCenter(i)));
          double E = h->GetXaxis()->GetBinCenter(i);
          events += h->GetBinContent(i,j)*flux;
        }
        h_v.at(i_f)->SetBinContent(j,events);
      }

      h_nofsi_v.at(i_f) = new TH1D(("h_reco_nofsi_"+gen+"_"+est).c_str(),";E_{est} (GeV);Events",200,0.2,8.0); 
      for(int j=1;j<h_nofsi_v.at(i_f)->GetNbinsX()+1;j++){
        double events = 0.0;
        for(int i=1;i<h_nofsi->GetNbinsX()+1;i++){
          double flux = h_flux->GetBinContent(h_flux->FindBin(h_nofsi->GetXaxis()->GetBinCenter(i)));
          double E = h_nofsi->GetXaxis()->GetBinCenter(i);
          events += h_nofsi->GetBinContent(i,j)*flux;
        }
        h_nofsi_v.at(i_f)->SetBinContent(j,events);
      }

      h_v.at(i_f)->Rebin();
      h_nofsi_v.at(i_f)->Rebin();

      h_v.at(i_f)->SetLineWidth(2);
      h_v.at(i_f)->SetLineColor(i_f+1);

      h_nofsi_v.at(i_f)->SetLineWidth(2);
      h_nofsi_v.at(i_f)->SetLineColor(i_f+1);
      h_nofsi_v.at(i_f)->SetLineStyle(2);

      hs->Add(h_v.at(i_f));
      hs->Add(h_nofsi_v.at(i_f));

      l->AddEntry(h_v.at(i_f),(gen+"").c_str(),"L");
      l->AddEntry(h_nofsi_v.at(i_f),(gen+" No FSI").c_str(),"L");

    }


    hs->Draw("nostack HIST");
    l->Draw();

    c->Print(("Plots/Spectrum_"+est+".pdf").c_str()); 
    c->Clear();
    l->Clear();

    delete hs;


  }


}
