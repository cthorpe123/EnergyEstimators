#include "../Funcs/Funcs.h"
#include "TLorentzVector.h"

void BiasVariancePlots(){

  std::vector<std::string> InputFiles_v = {"GENIEEvents.root","NEUTEvents.root","GiBUUEvents.root","NuWroEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NEUT","GiBUU","NuWro"};

  std::vector<double> binning_v;
  binning_v.push_back(0.35);
  for(int i=0;i<31;i++) binning_v.push_back(binning_v.back()+0.15);
  for(int i=0;i<10;i++) binning_v.push_back(binning_v.back()+0.3);
  int nbins = binning_v.size()-1;
  double* binning_a = &binning_v[0];

  std::vector<double> w_binning_v;
  w_binning_v.push_back(0.8);
  for(int i=0;i<21;i++) w_binning_v.push_back(w_binning_v.back()+0.2); 
  int w_nbins = w_binning_v.size()-1;
  double* w_binning_a = &w_binning_v[0];

  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_W(kMAX,std::vector<TH3D*>());

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    for(size_t i_e=0;i_e<kMAX;i_e++) h_TrueEnergy_RecoEnergy_W.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_W_"+estimators_str.at(i_e)).c_str(),"",nbins,binning_a,nbins,binning_a,w_nbins,w_binning_a));

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

      //if(ievent > 5000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;
      //if(generator == "GiBUU" && weight > 1) continue;

      if(nu_pdg != 14 || ccnc != 1) continue;

      double W = CalcW(pdg,p4);
      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);
      std::vector<TVector3> neutron_mom = GetNeutronMom(pdg,p4);

      if(nprot < 1) continue;

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,neutron_mom,i_e);
        h_TrueEnergy_RecoEnergy_W.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,W,weight);  
      }

    }

  }// i_f


  gSystem->Exec("mkdir -p Plots/BiasVariancePlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  // First calculate the 1D bias/variance plots
  std::vector<THStack*> hs_Bias,hs_Variance; 
  std::vector<std::vector<TH1D*>> h_Bias(kMAX),h_Variance(kMAX);
  double min_bias = 1.0,max_bias = -1.0;
  double min_variance = 1.0,max_variance = -1.0;
  for(size_t i_e=0;i_e<kMAX;i_e++){

    hs_Bias.push_back(new THStack(("hs_Bias_"+estimators_str.at(i_e)).c_str(),";True Neutrino Energy (GeV);Frac. Bias"));
    hs_Variance.push_back(new THStack(("hs_Variance_"+estimators_str.at(i_e)).c_str(),";True Neutrino Energy (GeV);Frac. Variance"));

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      TH2D* h_proj = static_cast<TH2D*>(h_TrueEnergy_RecoEnergy_W.at(i_e).at(i_f)->Project3D("yx"));
      h_Bias.at(i_e).push_back(new TH1D(("h_Bias_"+Generators_v.at(i_f)+"_"+estimators_str.at(i_e)).c_str(),"",nbins,binning_a));
      h_Variance.at(i_e).push_back(new TH1D(("h_Variance_"+Generators_v.at(i_f)+"_"+estimators_str.at(i_e)).c_str(),"",nbins,binning_a));
      GetBiasVariance(h_proj,h_Bias.at(i_e).back(),h_Variance.at(i_e).back()); 

      h_Bias.at(i_e).back()->SetLineColor(i_f+2);
      h_Bias.at(i_e).back()->SetLineWidth(2);
      hs_Bias.back()->Add(h_Bias.at(i_e).back());
      min_bias = std::min(min_bias,h_Bias.at(i_e).back()->GetBinContent(h_Bias.at(i_e).back()->GetMinimumBin()));
      max_bias = std::max(max_bias,h_Bias.at(i_e).back()->GetBinContent(h_Bias.at(i_e).back()->GetMaximumBin()));

      h_Variance.at(i_e).back()->SetLineColor(i_f+2);
      h_Variance.at(i_e).back()->SetLineWidth(2);
      hs_Variance.back()->Add(h_Variance.at(i_e).back());
      min_variance = std::min(min_variance,h_Variance.at(i_e).back()->GetBinContent(h_Variance.at(i_e).back()->GetMinimumBin()));
      max_variance = std::max(max_variance,h_Variance.at(i_e).back()->GetBinContent(h_Variance.at(i_e).back()->GetMaximumBin()));

      if(i_e == 0) l->AddEntry(h_Bias.at(i_e).back(),Generators_v.at(i_f).c_str(),"L");

      delete h_proj;

    } 

  }

  for(size_t i_e=0;i_e<kMAX;i_e++){
    hs_Bias.at(i_e)->Draw("nostack HIST");
    hs_Bias.at(i_e)->SetMaximum(max_bias);
    hs_Bias.at(i_e)->SetMinimum(min_bias);
    l->Draw();
    c->Print(("Plots/BiasVariancePlots/Bias_Energy_"+estimators_str.at(i_e)+".png").c_str());  
    c->Clear();

    hs_Variance.at(i_e)->Draw("nostack HIST");
    hs_Variance.at(i_e)->SetMaximum(max_variance);
    hs_Variance.at(i_e)->SetMinimum(min_variance);
    l->Draw();
    c->Print(("Plots/BiasVariancePlots/Variance_Energy_"+estimators_str.at(i_e)+".png").c_str());  
    c->Clear();
  }

  // Calculate bias as a function of hadronic invariant mass as well  
  hs_Bias.clear();
  h_Bias.clear();
  h_Bias.resize(kMAX);
  min_bias = 1.0,max_bias = -1.0;
  for(size_t i_e=0;i_e<kMAX;i_e++){

    hs_Bias.push_back(new THStack(("hs_Bias_W_"+estimators_str.at(i_e)).c_str(),";Visible Hadronic Invariant Mass (GeV);Frac. Bias"));

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      //std::cout << Generators_v.at(i_f) << std::endl;
      const TH3D* h = h_TrueEnergy_RecoEnergy_W.at(i_e).at(i_f);
      h_Bias.at(i_e).push_back(new TH1D(("h_Bias_W_"+Generators_v.at(i_f)+"_"+estimators_str.at(i_e)).c_str(),"",w_nbins,w_binning_a));

      // Calculation of bias afo visible invariant mass
      for(int i_w=1;i_w<h->GetNbinsZ()+1;i_w++){
        double bias = 0.0;
        double events = 0.0;
        for(int i=1;i<h->GetNbinsX()+1;i++){
          for(int j=1;j<h->GetNbinsY()+1;j++){
            bias += (h->GetYaxis()->GetBinCenter(j) - h->GetXaxis()->GetBinCenter(i))*h->GetBinContent(i,j,i_w)/h->GetXaxis()->GetBinCenter(i);
            events += h->GetBinContent(i,j,i_w); 
          }
        }
        if(events > 0)
          h_Bias.at(i_e).back()->SetBinContent(i_w,bias/events);  
      } 

      h_Bias.at(i_e).back()->SetLineColor(i_f+2);
      h_Bias.at(i_e).back()->SetLineWidth(2);
      hs_Bias.at(i_e)->Add(h_Bias.at(i_e).back());

      min_bias = std::min(min_bias,h_Bias.at(i_e).back()->GetBinContent(h_Bias.at(i_e).back()->GetMinimumBin()));
      max_bias = std::max(max_bias,h_Bias.at(i_e).back()->GetBinContent(h_Bias.at(i_e).back()->GetMaximumBin()));

      //l->AddEntry(h_Bias.at(i_e).back(),Generators_v.at(i_f).c_str(),"L");

    }

  }

  for(size_t i_e=0;i_e<kMAX;i_e++){
    hs_Bias.at(i_e)->Draw("nostack HIST");
    hs_Bias.at(i_e)->SetMaximum(max_bias);
    hs_Bias.at(i_e)->SetMinimum(min_bias);
    l->Draw();
    c->Print(("Plots/BiasVariancePlots/Bias_W_"+estimators_str.at(i_e)+".png").c_str());  
    c->Clear();
  }


  // Calculate 2D bias plots for invariant mass and true neutrino energy 
  std::vector<std::vector<TH2D*>> h_Bias_2D(kMAX);
  min_bias = 1.0,max_bias = -1.0;
  for(size_t i_e=0;i_e<kMAX;i_e++){

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      const TH3D* h = h_TrueEnergy_RecoEnergy_W.at(i_e).at(i_f);
      h_Bias_2D.at(i_e).push_back(new TH2D(("h_Bias_2D_"+Generators_v.at(i_f)+"_"+estimators_str.at(i_e)).c_str(),";True Neutrino Energy (GeV);W_{vis} (GeV);Frac. Bias in Neutrino Energy",nbins,binning_a,w_nbins,w_binning_a));

      for(int i_w=1;i_w<h->GetNbinsZ()+1;i_w++){

        // Get a 2D slice of the 3D histogram at constant W
        TH2D* h_2d = static_cast<TH2D*>(h->Project3D("yx"));
        for(int i=1;i<h->GetNbinsX()+1;i++)
          for(int j=1;j<h->GetNbinsY()+1;j++)
            h_2d->SetBinContent(i,j,h->GetBinContent(i,j,i_w));

        // Then get the bias and variance histograms       
        TH1D* h_bias_tmp = new TH1D("h_bias_tmp","",nbins,binning_a);
        TH1D* h_variance_tmp = new TH1D("h_variance_tmp","",nbins,binning_a);
        GetBiasVariance(h_2d,h_bias_tmp,h_variance_tmp);

        // Set the content of the 2D bias/variance plots 
        for(int i=1;i<h->GetNbinsX()+1;i++)
          h_Bias_2D.at(i_e).back()->SetBinContent(i,i_w,h_bias_tmp->GetBinContent(i));

        delete h_bias_tmp;
        delete h_variance_tmp; 
        delete h_2d;

        min_bias = std::min(min_bias,h_Bias.at(i_e).back()->GetBinContent(h_Bias.at(i_e).back()->GetMinimumBin()));
        max_bias = std::max(max_bias,h_Bias.at(i_e).back()->GetBinContent(h_Bias.at(i_e).back()->GetMaximumBin()));

      }
    }
  }

  for(size_t i_e=0;i_e<kMAX;i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      h_Bias_2D.at(i_e).at(i_f)->Draw("colz");
      h_Bias_2D.at(i_e).at(i_f)->SetMaximum(max_bias);
      h_Bias_2D.at(i_e).at(i_f)->SetMinimum(min_bias);
      h_Bias_2D.at(i_e).at(i_f)->SetStats(0);  
      c->Print(("Plots/BiasVariancePlots/Bias_2D_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)+".png").c_str());  
      c->Clear();
    }
  }

}
