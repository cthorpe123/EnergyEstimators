#include "../../Funcs/Funcs.h"
#include "../../Funcs/EnergyEstimatorFuncs.h"
#include "../../Funcs/PlotSetup.h"
#include "TLorentzVector.h"

void GetWeightedBias(TH2D* h_all,TH2D* h_sub,TH1D* h_Bias){

  int nbins_x = h_sub->GetNbinsX();
  int nbins_y = h_sub->GetNbinsY();

  for(int i=1;i<nbins_x+1;i++){

    double mean = 0.0;
    double events = 0.0;
    double events_all = 0.0;

    for(int j=1;j<nbins_y+1;j++){
      mean += h_sub->GetYaxis()->GetBinCenter(j)*h_sub->GetBinContent(i,j);
      events += h_sub->GetBinContent(i,j);
      events_all += h_all->GetBinContent(i,j);
    }

    mean /= events;

    h_Bias->SetBinContent(i,(events/events_all)*(mean - h_Bias->GetBinCenter(i))/h_Bias->GetBinCenter(i));

    if(events == 0.0) h_Bias->SetBinContent(i,0);

  }

}

void InvestigatingNuWro(){

  PlotSetup();

  std::vector<std::string> InputFiles_v = {"NuWroEvents2.root"};
  std::vector<std::string> Generators_v = {"NuWro"};

  std::vector<double> true_binning_v;
  true_binning_v.push_back(0.5);
  for(int i=0;i<10;i++) true_binning_v.push_back(true_binning_v.back()+0.25);
  for(int i=0;i<4;i++) true_binning_v.push_back(true_binning_v.back()+0.5);
  true_binning_v.push_back(true_binning_v.back()+1.0);
  int nbins = true_binning_v.size()-1;
  double* binning_a = &true_binning_v[0];

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_NoFSI(kMAX,std::vector<TH2D*>());

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_1p0pi(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_1p0pi_NoFSI(kMAX,std::vector<TH2D*>());

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_Np0pi(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_Np0pi_NoFSI(kMAX,std::vector<TH2D*>());

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_1pNpi(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_1pNpi_NoFSI(kMAX,std::vector<TH2D*>());

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_NpNpi(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_NpNpi_NoFSI(kMAX,std::vector<TH2D*>());

  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_Not1p0pi(kMAX,std::vector<TH2D*>());
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_Not1p0pi_NoFSI(kMAX,std::vector<TH2D*>());




  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    for(size_t i_e=0;i_e<kMAX;i_e++){
      std::string est = estimators_str.at(i_e);
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));
      h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_NoFSI_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));

      h_TrueEnergy_RecoEnergy_1p0pi.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_1p0pi_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));
      h_TrueEnergy_RecoEnergy_1p0pi_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_1p0pi_NoFSI_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));

      h_TrueEnergy_RecoEnergy_Np0pi.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_Np0pi_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));
      h_TrueEnergy_RecoEnergy_Np0pi_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_Np0pi_NoFSI_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));

      h_TrueEnergy_RecoEnergy_1pNpi.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_1pNpi_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));
      h_TrueEnergy_RecoEnergy_1pNpi_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_1pNpi_NoFSI_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));

      h_TrueEnergy_RecoEnergy_NpNpi.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_NpNpi_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));
      h_TrueEnergy_RecoEnergy_NpNpi_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_NpNpi_NoFSI_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));

      h_TrueEnergy_RecoEnergy_Not1p0pi.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_Not1p0pi_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));
      h_TrueEnergy_RecoEnergy_Not1p0pi_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_Not1p0pi_NoFSI_"+est).c_str(),"",nbins,binning_a,200,-1,10.0));

    }

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

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 200000) break;
      if(ievent % 50000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != 14 || ccnc != 1) continue;

      //if(nu_e < min_e || nu_e > max_e) continue;

      std::vector<double> energies = GetEnergyEst(lepton_p4,pdg,p4);
      std::vector<double> energies_nofsi = GetEnergyEst(lepton_p4,pdg_nofsi,p4_nofsi);

      int np = 0;
      int npi = 0;
      for(size_t i_p=0;i_p<pdg->size();i_p++){
        if(pdg->at(i_p) == 2212 && p4->at(i_p).Vect().Mag() > 0.3) np++;
        if(abs(pdg->at(i_p)) == 211 && p4->at(i_p).Vect().Mag() > 0.1) npi++;
        if(pdg->at(i_p) == 111) npi++;
      }

      int np_nofsi = 0;
      int npi_nofsi = 0;
      for(size_t i_p=0;i_p<pdg_nofsi->size();i_p++){
        if(pdg_nofsi->at(i_p) == 2212 && p4_nofsi->at(i_p).Vect().Mag() > 0.3) np_nofsi++;
        if(abs(pdg_nofsi->at(i_p)) == 211 && p4_nofsi->at(i_p).Vect().Mag() > 0.1) npi_nofsi++;
        if(pdg_nofsi->at(i_p) == 111) npi_nofsi++;
      }

      for(int i_e=0;i_e<kMAX;i_e++){

        double nu_e_reco = energies.at(i_e);
        double nu_e_reco_nofsi = energies_nofsi.at(i_e);

        if(nu_e_reco > 0){
          h_TrueEnergy_RecoEnergy.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight); 
          if(np == 1 && npi == 0) h_TrueEnergy_RecoEnergy_1p0pi.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight); 
          else if(np > 1 && npi == 0) h_TrueEnergy_RecoEnergy_Np0pi.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight); 
          else if(np == 1 && npi > 0) h_TrueEnergy_RecoEnergy_1pNpi.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight); 
          else h_TrueEnergy_RecoEnergy_NpNpi.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight); 
          if(!(np == 1 && npi == 0)) h_TrueEnergy_RecoEnergy_Not1p0pi.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco,weight);
        }

        if(nu_e_reco_nofsi > 0){
          h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight); 
          if(np_nofsi == 1 && npi_nofsi == 0) h_TrueEnergy_RecoEnergy_1p0pi_NoFSI.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight); 
          else if(np_nofsi > 1 && npi_nofsi == 0) h_TrueEnergy_RecoEnergy_Np0pi_NoFSI.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight); 
          else if(np_nofsi == 1 && npi_nofsi > 0) h_TrueEnergy_RecoEnergy_1pNpi_NoFSI.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight); 
          else h_TrueEnergy_RecoEnergy_NpNpi_NoFSI.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight); 
          if(!(np_nofsi == 1 && npi_nofsi == 0)) h_TrueEnergy_RecoEnergy_Not1p0pi_NoFSI.at(i_e).at(i_f)->Fill(nu_e,nu_e_reco_nofsi,weight);
        }
      }

    }

  }// i_f


  //TCanvas* c = new TCanvas("c","c");

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);
    std::vector<TH1D*> h_Bias,h_Bias_NoFSI;
    THStack* hs_Bias_Change = new THStack(("hs_Bias_Change_"+gen).c_str(),";True Neutrino Energy (GeV);Change in Frac. Bias");

    for(size_t i_e=0;i_e<kMAX;i_e++){

      if(i_e != kMuonKinWNP && i_e != kSFMethod) continue;

      std::string est = estimators_str.at(i_e);

      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f); 
      h_Bias.push_back(new TH1D(("h_Bias_"+gen+"_"+est).c_str(),"",nbins,binning_a));
      GetWeightedBias(h,h,h_Bias.back()); 

      TH2D* h_nofsi = h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f); 
      h_Bias_NoFSI.push_back(new TH1D(("h_Bias_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a));
      GetWeightedBias(h_nofsi,h_nofsi,h_Bias_NoFSI.back()); 

      h_Bias.back()->SetLineWidth(2);
      h_Bias.back()->SetLineColor(colors.at(i_e));

      TH1D* h_bias = h_Bias.back();
      TH1D* h_bias_nofsi = h_Bias_NoFSI.back();
      h_bias->Add(h_bias_nofsi,-1);
      hs_Bias_Change->Add(h_bias);



    }
    hs_Bias_Change->Draw("HIST nostack");
    hs_Bias_Change->SetMinimum(-0.06);
    hs_Bias_Change->SetMaximum(0.022);
    c->Print(("Plots/Bias_Change_"+gen+".pdf").c_str());
    p_plot->Clear();
    //c->Clear();

    for(size_t i=0;i<h_Bias.size();i++){
      delete h_Bias.at(i);
      delete h_Bias_NoFSI.at(i);
    }

    delete hs_Bias_Change;
    l->Clear();

  }


  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);
    std::vector<TH1D*> h_Bias;
    THStack* hs_Bias_Change = new THStack(("hs_Bias_Change_"+gen).c_str(),";True Neutrino Energy (GeV);Change in Frac. Bias");

    for(size_t i_e=0;i_e<kMAX;i_e++){

      if(i_e != kMuonKinWNP && i_e != kSFMethod) continue;

      std::string est = estimators_str.at(i_e);
      std::string est_leg = estimators_leg.at(i_e);
      
      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f); 
      TH2D* h_1p0pi = h_TrueEnergy_RecoEnergy_1p0pi.at(i_e).at(i_f); 
      TH2D* h_Np0pi = h_TrueEnergy_RecoEnergy_Np0pi.at(i_e).at(i_f); 
      TH2D* h_1pNpi = h_TrueEnergy_RecoEnergy_1pNpi.at(i_e).at(i_f); 
      TH2D* h_NpNpi = h_TrueEnergy_RecoEnergy_NpNpi.at(i_e).at(i_f); 

      TH2D* h_NoFSI = h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f); 
      TH2D* h_1p0pi_NoFSI = h_TrueEnergy_RecoEnergy_1p0pi_NoFSI.at(i_e).at(i_f); 
      TH2D* h_Np0pi_NoFSI = h_TrueEnergy_RecoEnergy_Np0pi_NoFSI.at(i_e).at(i_f); 
      TH2D* h_1pNpi_NoFSI = h_TrueEnergy_RecoEnergy_1pNpi_NoFSI.at(i_e).at(i_f); 
      TH2D* h_NpNpi_NoFSI = h_TrueEnergy_RecoEnergy_NpNpi_NoFSI.at(i_e).at(i_f); 

      TH1D* h_Bias_1p0pi = new TH1D(("h_Bias_1p0pi_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h,h_1p0pi,h_Bias_1p0pi); 
      TH1D* h_Bias_1p0pi_NoFSI = new TH1D(("h_Bias_1p0pi_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h_NoFSI,h_1p0pi_NoFSI,h_Bias_1p0pi_NoFSI); 
      h_Bias_1p0pi->Add(h_Bias_1p0pi_NoFSI,-1);

      TH1D* h_Bias_Np0pi = new TH1D(("h_Bias_Np0pi_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h,h_Np0pi,h_Bias_Np0pi); 
      TH1D* h_Bias_Np0pi_NoFSI = new TH1D(("h_Bias_Np0pi_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h_NoFSI,h_Np0pi_NoFSI,h_Bias_Np0pi_NoFSI); 
      h_Bias_Np0pi->Add(h_Bias_Np0pi_NoFSI,-1);

      TH1D* h_Bias_1pNpi = new TH1D(("h_Bias_1pNpi_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h,h_1pNpi,h_Bias_1pNpi); 
      TH1D* h_Bias_1pNpi_NoFSI = new TH1D(("h_Bias_1pNpi_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h_NoFSI,h_1pNpi_NoFSI,h_Bias_1pNpi_NoFSI); 
      h_Bias_1pNpi->Add(h_Bias_1pNpi_NoFSI,-1);

      TH1D* h_Bias_NpNpi = new TH1D(("h_Bias_NpNpi_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h,h_NpNpi,h_Bias_NpNpi); 
      TH1D* h_Bias_NpNpi_NoFSI = new TH1D(("h_Bias_NpNpi_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h_NoFSI,h_NpNpi_NoFSI,h_Bias_NpNpi_NoFSI); 
      h_Bias_NpNpi->Add(h_Bias_NpNpi_NoFSI,-1);

      h_Bias_1p0pi->SetLineColor(colors.at(i_e)); 
      h_Bias_1p0pi->SetLineWidth(2);
      h_Bias_1p0pi->SetLineStyle(2);
      h_Bias.push_back(h_Bias_1p0pi);  

      h_Bias_Np0pi->SetLineColor(colors.at(i_e)); 
      h_Bias_Np0pi->SetLineWidth(2);
      h_Bias_Np0pi->SetLineStyle(3);
      h_Bias.push_back(h_Bias_Np0pi);  

      h_Bias_1pNpi->SetLineColor(colors.at(i_e)); 
      h_Bias_1pNpi->SetLineWidth(2);
      h_Bias_1pNpi->SetLineStyle(4);
      h_Bias.push_back(h_Bias_1pNpi);  

      h_Bias_NpNpi->SetLineColor(colors.at(i_e)); 
      h_Bias_NpNpi->SetLineWidth(2);
      h_Bias_NpNpi->SetLineStyle(5);
      h_Bias.push_back(h_Bias_NpNpi);  

      l->AddEntry(h_Bias_1p0pi,(est_leg+" 1p0pi").c_str(),"L");
      hs_Bias_Change->Add(h_Bias_1p0pi);      
      if(i_e != kSFMethod){
        l->AddEntry(h_Bias_Np0pi,(est_leg+" Np0pi").c_str(),"L");
        hs_Bias_Change->Add(h_Bias_Np0pi);      
        l->AddEntry(h_Bias_1pNpi,(est_leg+" 1pNpi").c_str(),"L");
        hs_Bias_Change->Add(h_Bias_1pNpi);      
        l->AddEntry(h_Bias_NpNpi,(est_leg+" NpNpi").c_str(),"L");
        hs_Bias_Change->Add(h_Bias_NpNpi);      
      }

    }

    hs_Bias_Change->Draw("HIST nostack");
    hs_Bias_Change->SetMinimum(-0.06);
    hs_Bias_Change->SetMaximum(0.022);
    c->Print(("Plots/Bias_Change_Comp_"+gen+".pdf").c_str());
    p_plot->Clear();

    for(TH1D* h : h_Bias) delete h;
    h_Bias.clear();

    l->Clear();

    delete hs_Bias_Change;

  }

  p_legend->cd();
  l->Clear();

  TLegend* l2 = new TLegend(0.1,0.0,0.9,0.95);

  l2->Draw();
  l2->SetNColumns(3);
  l2->SetBorderSize(0);   

  p_plot->cd();

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);
    std::vector<TH1D*> h_Bias;
    THStack* hs_Bias_Change = new THStack(("hs_Bias_Change_"+gen).c_str(),";True Neutrino Energy (GeV);Change in Frac. Bias");

    for(size_t i_e=0;i_e<kMAX;i_e++){

      if(i_e != kMuonKinWNP && i_e != kSFMethod) continue;

      std::string est = estimators_str.at(i_e);
      std::string est_leg = estimators_leg.at(i_e);

      TH2D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f); 
      TH2D* h_1p0pi = h_TrueEnergy_RecoEnergy_1p0pi.at(i_e).at(i_f); 
      TH2D* h_Not1p0pi = h_TrueEnergy_RecoEnergy_Not1p0pi.at(i_e).at(i_f); 

      TH2D* h_NoFSI = h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f); 
      TH2D* h_1p0pi_NoFSI = h_TrueEnergy_RecoEnergy_1p0pi_NoFSI.at(i_e).at(i_f); 
      TH2D* h_Not1p0pi_NoFSI = h_TrueEnergy_RecoEnergy_Not1p0pi_NoFSI.at(i_e).at(i_f); 

      TH1D* h_Bias_All = new TH1D(("h_Bias_All_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h,h,h_Bias_All); 

      TH1D* h_Bias_All_NoFSI = new TH1D(("h_Bias_All_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h_NoFSI,h_NoFSI,h_Bias_All_NoFSI); 

      TH1D* h_Bias_1p0pi = new TH1D(("h_Bias_1p0pi_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h,h_1p0pi,h_Bias_1p0pi); 

      TH1D* h_Bias_1p0pi_NoFSI = new TH1D(("h_Bias_1p0pi_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h_NoFSI,h_1p0pi_NoFSI,h_Bias_1p0pi_NoFSI); 

      TH1D* h_Bias_Not1p0pi = new TH1D(("h_Bias_Not1p0pi_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h,h_Not1p0pi,h_Bias_Not1p0pi); 

      TH1D* h_Bias_Not1p0pi_NoFSI = new TH1D(("h_Bias_Not1p0pi_NoFSI_"+gen+"_"+est).c_str(),"",nbins,binning_a);
      GetWeightedBias(h_NoFSI,h_Not1p0pi_NoFSI,h_Bias_Not1p0pi_NoFSI); 

      h_Bias_All->Add(h_Bias_All_NoFSI,-1);
      h_Bias_1p0pi->Add(h_Bias_1p0pi_NoFSI,-1);
      h_Bias_Not1p0pi->Add(h_Bias_Not1p0pi_NoFSI,-1);

      h_Bias_All->SetLineColor(colors.at(i_e)); 
      h_Bias_All->SetLineWidth(2);
      h_Bias.push_back(h_Bias_All);  

      h_Bias_1p0pi->SetLineColor(colors.at(i_e)); 
      h_Bias_1p0pi->SetLineWidth(2);
      h_Bias_1p0pi->SetLineStyle(2);
      h_Bias.push_back(h_Bias_1p0pi);  

      h_Bias_Not1p0pi->SetLineColor(colors.at(i_e)); 
      h_Bias_Not1p0pi->SetLineWidth(2);
      h_Bias_Not1p0pi->SetLineStyle(3);
      h_Bias.push_back(h_Bias_Not1p0pi);  

      l2->AddEntry(h_Bias_All,est_leg.c_str(),"L");
      hs_Bias_Change->Add(h_Bias_All);      

      if(i_e != kSFMethod){
        l2->AddEntry(h_Bias_1p0pi,(est_leg+" 1p0pi").c_str(),"L");
        hs_Bias_Change->Add(h_Bias_1p0pi);      
        l2->AddEntry(h_Bias_Not1p0pi,(est_leg+" NpXpi").c_str(),"L");
        hs_Bias_Change->Add(h_Bias_Not1p0pi);      
      }

    }

    hs_Bias_Change->Draw("HIST nostack");
    hs_Bias_Change->SetMinimum(-0.06);
    hs_Bias_Change->SetMaximum(0.01);
    c->Print(("Plots/Bias_Change_Comp2_"+gen+".pdf").c_str());
    p_plot->Clear();

    for(TH1D* h : h_Bias) delete h;
    h_Bias.clear();

  }

}
