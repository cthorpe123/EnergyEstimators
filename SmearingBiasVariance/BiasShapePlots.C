#include "../Funcs/Funcs.h"
#include "TLorentzVector.h"

void BiasShapePlots(){

  int n_w_points = 50;

  std::vector<std::string> InputFiles_v = {"GENIEEvents.root","NEUTEvents.root","GiBUUEvents.root","NuWroEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NEUT","GiBUU","NuWro"};

  std::vector<std::vector<TH2D*>> h_EnergyBias_W(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_EnergyBias_W_HighBinning(kMAX,std::vector<TH2D*>());

  std::vector<std::vector<TH2D*>> h_EnergyBias_Energy(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_EnergyBias_Energy_HighBinning(kMAX,std::vector<TH2D*>());

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

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);

    for(size_t i_e=0;i_e<kMAX;i_e++) h_EnergyBias_W.at(i_e).push_back(new TH2D((generator+"_EnergyBias_W_"+estimators_str.at(i_e)).c_str(),"",100,-1,1,5,0.0,5.0));
    for(size_t i_e=0;i_e<kMAX;i_e++) h_EnergyBias_W_HighBinning.at(i_e).push_back(new TH2D((generator+"_EnergyBias_W__HighBinning"+estimators_str.at(i_e)).c_str(),"",10000,-1,1,n_w_points,0.93,5.0));

    for(size_t i_e=0;i_e<kMAX;i_e++) h_EnergyBias_Energy.at(i_e).push_back(new TH2D((generator+"_EnergyBias_Energy_"+estimators_str.at(i_e)).c_str(),"",100,-1,1,5,0.2,8.0));
    for(size_t i_e=0;i_e<kMAX;i_e++) h_EnergyBias_Energy_HighBinning.at(i_e).push_back(new TH2D((generator+"_EnergyBias_Energy__HighBinning"+estimators_str.at(i_e)).c_str(),"",10000,-1,1,n_w_points,0.2,8.0));

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 200000) break;
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
        if(nu_e_reco < 0) continue;
        h_EnergyBias_W.at(i_e).at(i_f)->Fill((nu_e_reco-nu_e)/nu_e,W,weight);  
        h_EnergyBias_W_HighBinning.at(i_e).at(i_f)->Fill((nu_e_reco-nu_e)/nu_e,W,weight);  

        h_EnergyBias_Energy.at(i_e).at(i_f)->Fill((nu_e_reco-nu_e)/nu_e,nu_e,weight);  
        h_EnergyBias_Energy_HighBinning.at(i_e).at(i_f)->Fill((nu_e_reco-nu_e)/nu_e,nu_e,weight);  
      }

    }

  }// i_f


  gSystem->Exec("mkdir -p Plots/BiasShapePlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      std::cout << estimators_str.at(i_e) << "  " << Generators_v.at(i_f) << std::endl;
      const TH2D* h = h_EnergyBias_W.at(i_e).at(i_f);
      THStack* hs = new THStack("hs",";(E_{e} - E_{t})/E_{t};");
      THStack* hs_cum = new THStack("hs",";(E_{e} - E_{t})/E_{t};Cumulative Dist.");
      TH1D* h_all = h->ProjectionX("h_bias_all");
      h_all->SetLineWidth(3);
      h_all->SetLineColor(1);
      h_all->Scale(1.0/h_all->Integral());
      TH1D* h_all_cum = static_cast<TH1D*>(h_all->GetCumulative("h_all_cum"));
      hs->Add(h_all);
      hs_cum->Add(h_all_cum);
      l->AddEntry(h_all,"All Events","L");
      std::vector<TH1D*> h_1D;
      std::vector<TH1D*> h_1D_cum;
      for(int i=1;i<h->GetNbinsY()+1;i++){
        double w_low = h->GetYaxis()->GetBinLowEdge(i);
        double w_high = h->GetYaxis()->GetBinLowEdge(i+1);

        h_1D.push_back(h->ProjectionX(("h_bias_"+std::to_string(i)).c_str(),i,i));
        h_1D.back()->SetLineColor(i+1);
        h_1D.back()->SetLineWidth(2);
        h_1D.back()->Scale(1.0/h_1D.back()->Integral());
        hs->Add(h_1D.back());  

        h_1D_cum.push_back(static_cast<TH1D*>(h_1D.back()->GetCumulative(("h_bias_cum_"+std::to_string(i)).c_str())));
        h_1D_cum.back()->SetLineColor(i+1);
        h_1D_cum.back()->SetLineWidth(2);
        hs_cum->Add(h_1D_cum.back());  

        std::string cap = to_string_with_precision(w_low) + " < W < " + to_string_with_precision(w_high) + " GeV";
        l->AddEntry(h_1D.back(),cap.c_str(),"L");
      }

      hs->Draw("HIST nostack");
      l->Draw();
      c->Print(("Plots/BiasShapePlots/Bias_W_"+estimators_str.at(i_e) + "_" + Generators_v.at(i_f) + ".png").c_str());
      c->Clear();

      hs_cum->Draw("HIST nostack");
      l->Draw();
      c->Print(("Plots/BiasShapePlots/Cum_Bias_W_"+estimators_str.at(i_e) + "_" + Generators_v.at(i_f) + ".png").c_str());
      c->Clear();

      l->Clear();

      delete hs;
      delete hs_cum;

    }
  }

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      std::cout << estimators_str.at(i_e) << "  " << Generators_v.at(i_f) << std::endl;
      const TH2D* h = h_EnergyBias_Energy.at(i_e).at(i_f);
      THStack* hs = new THStack("hs",";(E_{e} - E_{t})/E_{t};");
      TH1D* h_all = h->ProjectionX("h_bias_all");
      h_all->SetLineWidth(3);
      h_all->SetLineColor(1);
      h_all->Scale(1.0/h_all->Integral());
      hs->Add(h_all);
      l->AddEntry(h_all,"All Events","L");
      std::vector<TH1D*> h_1D;
      for(int i=1;i<h->GetNbinsY()+1;i++){
        double w_low = h->GetYaxis()->GetBinLowEdge(i);
        double w_high = h->GetYaxis()->GetBinLowEdge(i+1);
        h_1D.push_back(h->ProjectionX(("h_bias_"+std::to_string(i)).c_str(),i,i));
        h_1D.back()->SetLineColor(i+1);
        h_1D.back()->SetLineWidth(2);
        h_1D.back()->Scale(1.0/h_1D.back()->Integral());
        hs->Add(h_1D.back());  
        std::string cap = to_string_with_precision(w_low) + " < E_{t} < " + to_string_with_precision(w_high) + " GeV";
        l->AddEntry(h_1D.back(),cap.c_str(),"L");
      }
      hs->Draw("HIST nostack");
      l->Draw();
      c->Print(("Plots/BiasShapePlots/Bias_Energy_"+estimators_str.at(i_e) + "_" + Generators_v.at(i_f) + ".png").c_str());
      c->Clear();
      l->Clear();

      delete hs;

    }
  }

  // Pseudo KS score as a function of W
  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    THStack* hs = new THStack("hs","");
    std::vector<TH1D*> h_scores;
    for(size_t i_e=0;i_e<estimators_str.size();i_e++){

      if(i_e == kPeLEELike0Pi) continue;

      // easiest way to get the binning in W is to just project one of the 2d histograms
      TH1D* h = h_EnergyBias_W_HighBinning.at(0).at(0)->ProjectionY(estimators_str.at(i_e).c_str());
      h_scores.push_back(h);
      h_scores.back()->Reset();

      for(int i_w=1;i_w<h_scores.back()->GetNbinsX()+1;i_w++){
        TH1D* h_all = h_EnergyBias_W_HighBinning.at(i_e).at(i_f)->ProjectionX("h_all");
        TH1D* h_w = h_EnergyBias_W_HighBinning.at(i_e).at(i_f)->ProjectionX("h_w",i_w,i_w);
        if(h_w->Integral() == 0) continue;
        h_scores.back()->SetBinContent(i_w,PseudoKSTest(h_all,h_w));
        delete h_all;
        delete h_w;
      }
        hs->Add(h_scores.back());
        l->AddEntry(h_scores.back(),estimators_str.at(i_e).c_str(),"L");
        h_scores.back()->SetLineColor(i_e+1);
        h_scores.back()->SetLineWidth(2);
    }

    hs->Draw("nostack HIST");     
    l->Draw();
    c->Print(("Plots/BiasShapePlots/KS_" + Generators_v.at(i_f) + ".png").c_str());
    c->Clear();
    l->Clear();

  }



  // Pseudo KS score as a function of true neutrino energy
  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    THStack* hs = new THStack("hs","");
    std::vector<TH1D*> h_scores;
    for(size_t i_e=0;i_e<estimators_str.size();i_e++){

      //if(i_e == kPeLEELike0Pi) continue;

      // easiest way to get the binning in Energy is to just project one of the 2d histograms
      TH1D* h = h_EnergyBias_Energy_HighBinning.at(0).at(0)->ProjectionY(estimators_str.at(i_e).c_str());
      h_scores.push_back(h);
      h_scores.back()->Reset();

      for(int i_w=1;i_w<h_scores.back()->GetNbinsX()+1;i_w++){
        TH1D* h_all = h_EnergyBias_Energy_HighBinning.at(i_e).at(i_f)->ProjectionX("h_all");
        TH1D* h_w = h_EnergyBias_Energy_HighBinning.at(i_e).at(i_f)->ProjectionX("h_w",i_w,i_w);
        if(h_w->Integral() == 0) continue;
        h_scores.back()->SetBinContent(i_w,PseudoKSTest(h_all,h_w));
        delete h_all;
        delete h_w;
      }
        hs->Add(h_scores.back());
        l->AddEntry(h_scores.back(),estimators_str.at(i_e).c_str(),"L");
        h_scores.back()->SetLineColor(i_e+1);
        h_scores.back()->SetLineWidth(2);
    }

    hs->Draw("nostack HIST");     
    l->Draw();
    c->Print(("Plots/BiasShapePlots/KS_Energy_" + Generators_v.at(i_f) + ".png").c_str());
    c->Clear();
    l->Clear();

  }

}
