#include "Funcs/Funcs.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;

// Detector exposure to normalise plots to in KT x 10^21 POT

double fid_mass = 10; // active mass in KT
double POT = 1; // POT in 10^21
double exposure = fid_mass*POT; 

void DeltaChi2Plots(){

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

    for(std::string estimator : estimators_str){
      h_TrueEnergy_RecoEnergy.back().push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);Events/KT/GeV^{2}/8^{21} POT",30,0.1,8.0,30,0.1,8.0));
      h_RecoEnergy_NoDeltaCP.back().push_back(new TH1D((generator+"_RecoEnergy_NoDeltaCP_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",30,0.1,8.0));
      h_RecoEnergy_DeltaCP_Plus.back().push_back(new TH1D((generator+"_RecoEnergy_DeltaCP_Plus_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",30,0.1,8.0));
      h_RecoEnergy_DeltaCP_Minus.back().push_back(new TH1D((generator+"_RecoEnergy_DeltaCP_Minus_"+estimator).c_str(),";Estimated Neutrino Energy (GeV);Events/KT/GeV/8^{21} POT",30,0.1,8.0));
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

      //if(ievent > 50000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;
      weight *= scale*1e38*40;

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
      }
    }

  }

  TFile* f_out = new TFile("rootfiles/DeltaChi2Plots.root","RECREATE");

  gSystem->Exec("mkdir -p Plots/NueAppPlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  // Oscillation plots
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      for(int i=1;i<h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetNbinsX()+1;i++)
        for(int j=1;j<h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetNbinsY()+1;j++)
          h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->SetBinContent(i,j,exposure*Rate(total_flux,h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->GetBinContent(i,j)));

      h_TrueEnergy_RecoEnergy.at(i_f).at(i_e)->Write();
      
    }
  }



  // Calculate the estimated spectrum for different deltaCP values
  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    // Use the first generator as the "model" we're going to compare the other generators to
    TH2D* h_true_reco_cv = h_TrueEnergy_RecoEnergy.at(0).at(i_e);
    int nbins_true = h_true_reco_cv->GetNbinsX();
    int nbins_reco = h_true_reco_cv->GetNbinsY();
    double min = h_true_reco_cv->GetYaxis()->GetBinLowEdge(1);
    double max = h_true_reco_cv->GetYaxis()->GetBinLowEdge(nbins_reco+1);

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      TH2D* h_delta_chi2_belt = new TH2D((Generators_v.at(i_f)+"_delta_chi2_belt_"+estimators_str.at(i_e)).c_str(),"#Delta#chi^{2} = #chi^{2}(0) - #chi^{2}(#delta_{CP}) ;CV #delta_{CP};Fake Data #delta_{CP};",100,-3.1415,3.1415,100,-3.1415,3.1415);

      std::cout << estimators_str.at(i_e) << "  " << Generators_v.at(i_f) << std::endl;

      TH2D* h_true_reco = h_TrueEnergy_RecoEnergy.at(i_f).at(i_e);

      TH1D* h_nodeltacp_cv = new TH1D("h_nodeltacp_cv","",nbins_reco,min,max);
      TH1D* h_nodeltacp = new TH1D("h_nodeltacp","",nbins_reco,min,max);

      // Calc the spectrum for zero deltaCP first
      for(int j=1;j<nbins_reco+1;j++){
        double events_cv = 0;
        double events = 0;
        for(int i=1;i<nbins_true+1;i++){
          events_cv += h_true_reco_cv->GetBinContent(i,j)*nue_app_prob(h_true_reco_cv->GetXaxis()->GetBinCenter(i));
          events += h_true_reco->GetBinContent(i,j)*nue_app_prob(h_true_reco->GetXaxis()->GetBinCenter(i));
        }

        h_nodeltacp_cv->SetBinContent(j,events_cv);
        h_nodeltacp_cv->SetBinError(j,sqrt(events_cv));

        h_nodeltacp->SetBinContent(j,events);
        h_nodeltacp->SetBinError(j,sqrt(events));

      } 

      for(int i_dcp_cv=1;i_dcp_cv<h_delta_chi2_belt->GetNbinsX()+1;i_dcp_cv++){
        for(int i_dcp=1;i_dcp<h_delta_chi2_belt->GetNbinsY()+1;i_dcp++){

          double deltacp_cv = h_delta_chi2_belt->GetXaxis()->GetBinCenter(i_dcp_cv);
          double deltacp = h_delta_chi2_belt->GetYaxis()->GetBinCenter(i_dcp);

          TH1D* h_deltacp_cv = new TH1D("h_deltacp_cv","",nbins_reco,min,max);
          TH1D* h_deltacp = new TH1D("h_deltacp","",nbins_reco,min,max);

          for(int j=1;j<nbins_reco+1;j++){
            double events_cv = 0;
            double events = 0;
            for(int i=1;i<nbins_true+1;i++){
              events += h_true_reco->GetBinContent(i,j)*nue_app_prob(h_true_reco->GetXaxis()->GetBinCenter(i),deltacp);
              events_cv += h_true_reco_cv->GetBinContent(i,j)*nue_app_prob(h_true_reco_cv->GetXaxis()->GetBinCenter(i),deltacp_cv);
            }
            h_deltacp->SetBinContent(j,events);
            h_deltacp->SetBinError(j,sqrt(events));

            h_deltacp_cv->SetBinContent(j,events_cv);
            h_deltacp_cv->SetBinError(j,sqrt(events_cv));
          } 

          // Calculate the chi2     
          double delta_chi2 = 0.0;
          for(int j=1;j<nbins_reco+1;j++){
            delta_chi2 +=  pow((h_deltacp->GetBinContent(j) - h_nodeltacp_cv->GetBinContent(j))/h_deltacp->GetBinError(j),2) 
              - pow((h_deltacp->GetBinContent(j) - h_deltacp_cv->GetBinContent(j))/h_deltacp->GetBinError(j),2);
          }

          h_delta_chi2_belt->SetBinContent(i_dcp_cv,i_dcp,delta_chi2/nbins_reco);

          delete h_deltacp;
          delete h_deltacp_cv;

        }
      }

      delete h_nodeltacp_cv;
      delete h_nodeltacp;


      h_delta_chi2_belt->Draw("colz");
      h_delta_chi2_belt->SetStats(0);
      c->Print(("Plots/NueAppPlots/DeltaChi2Belt_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)+".png").c_str());
      c->Clear();

      h_delta_chi2_belt->Write();

      //delete h_delta_chi2_belt;
    }

    delete h_true_reco_cv;

  }


  f_out->Close();

}

