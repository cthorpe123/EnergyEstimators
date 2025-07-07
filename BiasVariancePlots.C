#include "Funcs/Funcs.h"
#include "TLorentzVector.h"

void BiasVariancePlots(){

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("Flux/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEFD_fastmc.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  TH1D* h_flux_osc = static_cast<TH1D*>(f_flux->Get("numu_fluxosc"));
  h_flux->SetDirectory(0);
  h_flux_osc->SetDirectory(0);
  f_flux->Close();

  // Ratio of oscillated and unoscillated fluxes
  h_flux_osc->Divide(h_flux);

  std::vector<std::string> InputFiles_v = {"GENIEEvents.root","NEUTEvents.root","GiBUUEvents.root","NuWroEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NEUT","GiBUU","NuWro"};

  std::vector<double> binning_v;
  binning_v.push_back(0.2);
  for(int i=0;i<32;i++) binning_v.push_back(binning_v.back()+0.15);
  for(int i=0;i<10;i++) binning_v.push_back(binning_v.back()+0.3);
  int nbins = binning_v.size()-1;
  double* binning_a = &binning_v[0];

  std::vector<std::vector<TH1D*>> h_Bias_Energy(kMAX);
  std::vector<std::vector<TH1D*>> h_Variance_Energy(kMAX);
  std::vector<std::vector<TH1D*>> h_Bias_W(kMAX);
  std::vector<std::vector<TH1D*>> h_Variance_W(kMAX);

  std::vector<std::vector<TH2D*>> h_Bias_Energy_W(kMAX);
  std::vector<std::vector<TH2D*>> h_Variance_Energy_W(kMAX);

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    TFile* f = TFile::Open(("/exp/uboone/data/users/cthorpe/DIS/DUNE/rootfiles/"+InputFiles_v.at(i_f)).c_str());
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

    std::vector<double> true_e_v;
    std::vector<std::vector<double>> est_e_v(kMAX);
    std::vector<double> W_v;
    std::vector<double> weight_v;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 50000) break;
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

      true_e_v.push_back(nu_e);
      W_v.push_back(W);
      weight_v.push_back(weight);

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,i_e);
        est_e_v.at(i_e).push_back(nu_e_reco);    
      }

    }

    // Calcualte the bias and variance as a function of true neutrino energy 
    for(int i_e=0;i_e<kMAX;i_e++){
      std::cout << "Calculating bias and variance vs energy for estimator " << estimators_str.at(i_e) << std::endl;
      h_Bias_Energy.at(i_e).push_back(new TH1D(("h_Bias_Energy_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)).c_str(),"",nbins,binning_a));
      h_Variance_Energy.at(i_e).push_back(new TH1D(("h_Variance_Energy_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)).c_str(),"",nbins,binning_a));
      for(int i=1;i<h_Bias_Energy.at(i_e).back()->GetNbinsX()+1;i++){
        double min = h_Bias_Energy.at(i_e).back()->GetBinLowEdge(i); 
        double max = h_Bias_Energy.at(i_e).back()->GetBinLowEdge(i+1);
        double mean = 0.0;
        double bias = 0.0;
        double events = 0.0;      
        for(size_t i_ev=0;i_ev<true_e_v.size();i_ev++){
          if(true_e_v.at(i_ev) > min && true_e_v.at(i_ev) <= max && est_e_v.at(i_e).at(i_ev) > 0 && est_e_v.at(i_e).at(i_ev) < 8.0){
            mean += est_e_v.at(i_e).at(i_ev)*weight_v.at(i_ev);
            bias += (est_e_v.at(i_e).at(i_ev) - true_e_v.at(i_ev))*weight_v.at(i_ev);
            events += weight_v.at(i_ev);
        
          }
        }
        mean /= events;
        bias /= events;

        double var = 0.0;
        for(size_t i_ev=0;i_ev<true_e_v.size();i_ev++){
          if(true_e_v.at(i_ev) > min && true_e_v.at(i_ev) <= max && est_e_v.at(i_e).at(i_ev) > 0 && est_e_v.at(i_e).at(i_ev) < 8.0){
            var += pow(est_e_v.at(i_e).at(i_ev) - mean,2)*weight_v.at(i_ev);
          }
        }
        var /= events;       

        double center = (min+max)/2;
        if(events > 0){ 
          h_Bias_Energy.at(i_e).back()->SetBinContent(i,bias/center);
          h_Variance_Energy.at(i_e).back()->SetBinContent(i,var/center/center);
        }
      } 
    }


    // Calcualte the bias and variance as a function of visible invariant mass 
    for(int i_e=0;i_e<kMAX;i_e++){
      std::cout << "Calculating bias and variance vs W for estimator " << estimators_str.at(i_e) << std::endl;
      h_Bias_W.at(i_e).push_back(new TH1D(("h_Bias_W_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)).c_str(),"",50,0.8,5.0));
      h_Variance_W.at(i_e).push_back(new TH1D(("h_Variance_W_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)).c_str(),"",50,0.8,5.0));
      for(int i=1;i<h_Bias_W.at(i_e).back()->GetNbinsX()+1;i++){
        double min = h_Bias_W.at(i_e).back()->GetBinLowEdge(i); 
        double max = h_Bias_W.at(i_e).back()->GetBinLowEdge(i+1);
        double bias = 0.0;
        double events = 0.0;      
        for(size_t i_ev=0;i_ev<true_e_v.size();i_ev++){
          if(W_v.at(i_ev) > min && W_v.at(i_ev) <= max && est_e_v.at(i_e).at(i_ev) > 0 && est_e_v.at(i_e).at(i_ev) < 8.0){
            //std::cout << true_e_v.at(i_ev) << " " << est_e_v.at(i_e).at(i_ev) << std::endl;
            bias += ((est_e_v.at(i_e).at(i_ev) - true_e_v.at(i_ev))/true_e_v.at(i_ev))*weight_v.at(i_ev);
            events += weight_v.at(i_ev);
          }
        }
        bias /= events;

        if(events > 0){
          h_Bias_W.at(i_e).back()->SetBinContent(i,bias);
        }
      } 
    }

    // Calcualte the bias and variance as a function of true neutrino energy and W in 2D 
    for(int i_e=0;i_e<kMAX;i_e++){
      std::cout << "Calculating 2D bias and variance for estimator " << estimators_str.at(i_e) << std::endl;
      h_Bias_Energy_W.at(i_e).push_back(new TH2D(("h_Bias_Energy_W_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)).c_str(),"",nbins,binning_a,nbins,0.8,5.0));
      h_Variance_Energy_W.at(i_e).push_back(new TH2D(("h_Variance_Energy_W_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_f)).c_str(),"",nbins,binning_a,nbins,0.8,5.0));

      for(int i=1;i<nbins+1;i++){
        for(int j=1;j<nbins+1;j++){
          std::cout << "bin " << i << "/" << nbins+1 << " " << j  << "/" << nbins+1 << std::endl;

          double min_e = h_Bias_Energy_W.at(i_e).back()->GetXaxis()->GetBinLowEdge(i); 
          double max_e = h_Bias_Energy_W.at(i_e).back()->GetXaxis()->GetBinLowEdge(i+1);
          double min_w = h_Bias_Energy_W.at(i_e).back()->GetYaxis()->GetBinLowEdge(j); 
          double max_w = h_Bias_Energy_W.at(i_e).back()->GetYaxis()->GetBinLowEdge(j+1); 

          double mean = 0.0;
          double bias = 0.0;
          double events = 0.0;      

          for(size_t i_ev=0;i_ev<true_e_v.size();i_ev++){
            if(true_e_v.at(i_ev) > min_e && true_e_v.at(i_ev) <= max_e && W_v.at(i_ev) > min_w && W_v.at(i_ev) <= max_w && est_e_v.at(i_e).at(i_ev) > 0 && est_e_v.at(i_e).at(i_ev) < 8.0){
              mean += est_e_v.at(i_e).at(i_ev)*weight_v.at(i_ev);
              bias += (est_e_v.at(i_e).at(i_ev) - true_e_v.at(i_ev))*weight_v.at(i_ev);
              events += weight_v.at(i_ev);
            }
          }
          mean /= events;
          bias /= events;

          double var = 0.0;
          for(size_t i_ev=0;i_ev<true_e_v.size();i_ev++){
            if(true_e_v.at(i_ev) > min_e && true_e_v.at(i_ev) <= max_e && W_v.at(i_ev) > min_w && W_v.at(i_ev) <= max_w && est_e_v.at(i_e).at(i_ev) > 0 && est_e_v.at(i_e).at(i_ev) < 8.0){
              var += pow(est_e_v.at(i_e).at(i_ev) - mean,2)*weight_v.at(i_ev);
            }
          }
          var /= events;       


          double center_e = (min_e+max_e)/2;
          double center_w = (min_w+max_w)/2;
          if(events > 0){ 
            h_Bias_Energy_W.at(i_e).back()->SetBinContent(i,j,bias/center_e);
            h_Variance_Energy_W.at(i_e).back()->SetBinContent(i,j,var/center_e/center_e);
          }
        } 
      }

    }// i_e 

  }// i_f


  gSystem->Exec("mkdir -p Plots/BiasVariancePlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    THStack* hs_Bias = new THStack("hs_Bias",";True Neutrino Energy (GeV);Frac. Bias");
    THStack* hs_Variance = new THStack("hs_Variance",";True Neutrino Energy (GeV);Frac. Variance");
    THStack* hs_Bias_W = new THStack("hs_Bias_W",";W_{vis} (GeV);Frac. Bias");

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      l->AddEntry(h_Bias_Energy.at(i_e).at(i_f),Generators_v.at(i_f).c_str(),"L");

      h_Bias_Energy.at(i_e).at(i_f)->SetLineWidth(2);
      h_Bias_Energy.at(i_e).at(i_f)->SetLineColor(i_f+1);
      hs_Bias->Add(h_Bias_Energy.at(i_e).at(i_f));  

      h_Variance_Energy.at(i_e).at(i_f)->SetLineWidth(2);
      h_Variance_Energy.at(i_e).at(i_f)->SetLineColor(i_f+1);
      hs_Variance->Add(h_Variance_Energy.at(i_e).at(i_f));  

      h_Bias_W.at(i_e).at(i_f)->SetLineWidth(2);
      h_Bias_W.at(i_e).at(i_f)->SetLineColor(i_f+1);
      hs_Bias_W->Add(h_Bias_W.at(i_e).at(i_f));  

    }

    hs_Bias->Draw("HIST nostack");
    //hs_Bias->SetMaximum(0.2);
    //hs_Bias->SetMinimum(-0.5);
    gPad->Modified();
    l->Draw();
    c->Print(("Plots/BiasVariancePlots/Bias_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();

    hs_Bias_W->Draw("HIST nostack");
    //hs_Bias->SetMaximum(0.2);
    //hs_Bias->SetMinimum(-0.5);
    gPad->Modified();
    l->Draw();
    c->Print(("Plots/BiasVariancePlots/W_Bias_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();

    hs_Variance->Draw("HIST nostack");
    //hs_Variance->SetMaximum(0.25);
    //hs_Variance->SetMinimum(0.0);
    gPad->Modified();
    l->Draw();
    c->Print(("Plots/BiasVariancePlots/Variance_" + estimators_str.at(i_e) + ".png").c_str());
    c->Clear();

    l->Clear();

    delete hs_Bias;
    delete hs_Variance;
    delete hs_Bias_W;

  }


  for(size_t i_e=0;i_e<estimators_str.size();i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      h_Bias_Energy_W.at(i_e).at(i_f)->Draw("colz");
      h_Bias_Energy_W.at(i_e).at(i_f)->SetStats(0);
      h_Bias_Energy_W.at(i_e).at(i_f)->SetTitle("Fractional Bias;True Neutrino Energy (GeV);W_{vis} (GeV)");
      c->Print(("Plots/BiasVariancePlots/TwoD_Bias_" + estimators_str.at(i_e) + "_" + Generators_v.at(i_f) + ".png").c_str());      
      c->Clear();
      h_Variance_Energy_W.at(i_e).at(i_f)->Draw("colz");
      h_Variance_Energy_W.at(i_e).at(i_f)->SetStats(0);
      h_Bias_Energy_W.at(i_e).at(i_f)->SetTitle("Fractional Variance;True Neutrino Energy (GeV);W_{vis} (GeV)");
      c->Print(("Plots/BiasVariancePlots/TwoD_Variance_" + estimators_str.at(i_e) + "_" + Generators_v.at(i_f) + ".png").c_str());      
      c->Clear();
    }
  }

}
