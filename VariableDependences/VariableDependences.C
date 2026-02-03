#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/Funcs.h"
#include "../Funcs/PlotSetup.h"
#include "TLorentzVector.h"

void VariableDependences(){

  PlotSetup(); 

  bool rebin = false;

  std::vector<std::string> InputFiles_v = {"GENIEEventsFiltered.root","NEUTEventsFiltered.root","GiBUUEventsFiltered.root","NuWroEventsFiltered.root"};
  std::vector<std::string> Generators_v = {"GENIE","NEUT","GiBUU","NuWro"};
  std::vector<std::string> vars = {"W","Energy","Angle","MissingE","Neutrons","LeptonMom"};
  std::vector<std::string> axis_titles = {"W_{vis} (GeV)","E_{true} (GeV)","#theta_{lepton} (rad.)","Missing Hadronic Energy (GeV)","FS Neutrons","Muon Momentum (GeV)"};
  std::vector<std::string> dist_axis_titles = {"PDF (1/GeV)","PDF (1/GeV)","PDF","PDF (1/GeV)","PDF","PDF (1/GeV)"};

  std::vector<std::vector<std::vector<TH2D*>>> h_EnergyBias(kMAX,std::vector<std::vector<TH2D*>>(Generators_v.size(),std::vector<TH2D*>(vars.size())));
  std::vector<std::vector<std::vector<TH3D*>>> h_TrueEnergy_RecoEnergy(kMAX,std::vector<std::vector<TH3D*>>(Generators_v.size(),std::vector<TH3D*>(vars.size())));
  std::vector<std::vector<TH1D*>> h_var(Generators_v.size(),std::vector<TH1D*>(vars.size()));

  std::vector<int> points;
  std::vector<std::vector<double>> binning;
  std::vector<double*> binning_a;
  std::vector<std::pair<double,double>> bias_ranges;

  std::vector<double> bins;
  for(int i=0;i<201;i++) bins.push_back(0.2+i*(8.0-0.2)/200);
  int n_bins = bins.size()-1;
  double* bins_a = &bins[0];

  // W histogram setup
  std::vector<double> w_binning = {0.8,1.0,1.2,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0}; 
  int w_points = w_binning.size()-1;
  double* w_binning_a = &w_binning[0]; 
  points.push_back(w_points);
  binning.push_back(w_binning);
  binning_a.push_back(w_binning_a);
  bias_ranges.push_back(std::make_pair(-0.6,0.1));

  // Energy histogram setup
  std::vector<double> e_binning = {0.2,0.5,0.8,1.1,1.4,1.7,2.0,2.3,2.6,2.9,3.2,3.5,3.8,4.1,4.4,4.7,5.0,5.3,5.6,5.9,6.2,6.5,6.8,7.1,7.4,7.7,8.0}; 
  int e_points = e_binning.size()-1;
  double* e_binning_a = &e_binning[0]; 
  points.push_back(e_points);
  binning.push_back(e_binning);
  binning_a.push_back(e_binning_a);
  bias_ranges.push_back(std::make_pair(-0.45,0.1));

  // Angle histogram setup
  std::vector<double> ang_binning = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.142}; 
  int ang_points = ang_binning.size()-1;
  double* ang_binning_a = &ang_binning[0]; 
  points.push_back(ang_points);
  binning.push_back(ang_binning);
  binning_a.push_back(ang_binning_a);
  bias_ranges.push_back(std::make_pair(-0.65,0.1));

  // Missing energy histogram setup
  std::vector<double> miss_binning = {0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.2}; 
  int miss_points = miss_binning.size()-1;
  double* miss_binning_a = &miss_binning[0]; 
  points.push_back(miss_points);
  binning.push_back(miss_binning);
  binning_a.push_back(miss_binning_a);
  bias_ranges.push_back(std::make_pair(-0.6,0.0));

  // Neutrons histogram setup
  std::vector<double> n_binning = {-0.5,0.5,1.5,2.5,3.5,4.5,6.5}; 
  int n_points = n_binning.size()-1;
  double* n_binning_a = &n_binning[0]; 
  points.push_back(n_points);
  binning.push_back(n_binning);
  binning_a.push_back(n_binning_a);
  bias_ranges.push_back(std::make_pair(-0.65,0.0));

  // Muon Mom histogram setup
  std::vector<double> muon_p_binning = {0.3,0.5,0.8,1.1,1.4,1.7,2.0,2.3,2.6,2.9,3.2,3.5,3.8,4.1,4.4,4.7,5.0,5.3,5.6,5.9,6.2,6.5,6.8,7.1,7.4,7.7,8.0}; 
  int muon_p_points = muon_p_binning.size()-1;
  double* muon_p_binning_a = &muon_p_binning[0]; 
  points.push_back(muon_p_points);
  binning.push_back(muon_p_binning);
  binning_a.push_back(muon_p_binning_a);
  bias_ranges.push_back(std::make_pair(-1,1));

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);    

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
    int nprot;
    double W;
    std::vector<double>* est_nu_e=0;

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);
    t->SetBranchAddress("W",&W);
    t->SetBranchAddress("nprot",&nprot);
    t->SetBranchAddress("est_nu_e",&est_nu_e);


    for(size_t i_v=0;i_v<vars.size();i_v++){      
      std::string var = vars.at(i_v);
      for(size_t i_e=0;i_e<kMAX;i_e++){
        std::string estimator = estimators_str.at(i_e);
        h_EnergyBias.at(i_e).at(i_f).at(i_v) = new TH2D((gen+"_EnergyBias_"+var+"_"+estimator).c_str(),"",points.at(i_v),binning_a.at(i_v),100,0,2);
        h_TrueEnergy_RecoEnergy.at(i_e).at(i_f).at(i_v) = new TH3D((gen+"_TrueEnergy_RecoEnergy_"+var+"_"+estimator).c_str(),"",n_bins,bins_a,n_bins,bins_a,points.at(i_v),binning_a.at(i_v));
      }
      h_var.at(i_f).at(i_v) =  new TH1D((gen+"_"+var).c_str(),"",points.at(i_v),binning_a.at(i_v));
    }


    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 1000) break;
      if(ievent % 20000 == 0) std::cout << gen << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(gen != "GiBUU") weight = 1.0;
      //if(gen == "GiBUU" && weight > 1) continue;

      if(nu_pdg != 14 || ccnc != 1) continue;
      //if(nprot < 1) continue;

      double angle = acos(lepton_p4->CosTheta());
      double missing_e = GetMissingEnergy(pdg,p4); 
      int neutrons = 0;

      for(int i_p=0;i_p<pdg->size();i_p++)
        if(pdg->at(i_p) == 2112 && p4->at(i_p).Vect().Mag() > 0.3)
          neutrons++;

      std::vector<double> values = {W,nu_e,angle,missing_e,static_cast<double>(neutrons),lepton_p4->Vect().Mag()};

      for(size_t i_v=0;i_v<vars.size();i_v++){
        double val = values.at(i_v);
        h_var.at(i_f).at(i_v)->Fill(val,weight);
        for(int i_e=0;i_e<kMAX;i_e++){
          double nu_e_reco = est_nu_e->at(i_e);
           double bias = (nu_e_reco-nu_e)/nu_e;
          //double bias = nu_e_reco/nu_e;
          if(nu_e_reco < 0) continue;
          h_EnergyBias.at(i_e).at(i_f).at(i_v)->Fill(val,bias,weight);  
          h_TrueEnergy_RecoEnergy.at(i_e).at(i_f).at(i_v)->Fill(nu_e,nu_e_reco,val,weight);        
        } // i_e
      } // i_v

    } // ievent

  } // i_f

  gSystem->Exec("mkdir -p Plots/");

  // Bias AFO difference variables
  for(size_t i_v=0;i_v<vars.size();i_v++){

    std::string var = vars.at(i_v);
    std::string axis_title = axis_titles.at(i_v);

    gSystem->Exec(("mkdir -p Plots/" + var).c_str());

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);
      THStack* hs_Bias = new THStack(("hs_Bias"+gen+"_"+var).c_str(),(";"+axis_title+";Frac. Bias").c_str());
      THStack* hs_Variance = new THStack(("hs_Variance"+gen+"_"+var).c_str(),(";"+axis_title+";Frac. Variance").c_str());
      std::vector<TH1D*> h_Bias;
      std::vector<TH1D*> h_Variance;

      for(size_t i_e=0;i_e<kMAX;i_e++){

        if(var == "W" && i_e == kSFMethod) continue;
        std::string est = estimators_str.at(i_e);

        const TH3D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f).at(i_v);
        h_Bias.push_back(new TH1D(("h_Bias_W_"+gen+"_"+est+"_"+var).c_str(),"",points.at(i_v),binning_a.at(i_v)));
        h_Variance.push_back(new TH1D(("h_Variance_W_"+gen+"_"+est+"_"+var).c_str(),"",points.at(i_v),binning_a.at(i_v)));
        MakeBiasVarianceFrom3D(h,h_Bias.back(),h_Variance.back());
        h_Bias.back()->SetLineColor(colors.at(i_e));
        h_Bias.back()->SetLineWidth(2);
        hs_Bias->Add(h_Bias.back());

        h_Variance.back()->SetLineColor(colors.at(i_e));
        h_Variance.back()->SetLineWidth(2);
        hs_Variance->Add(h_Variance.back());

        l->AddEntry(h_Bias.back(),estimators_leg.at(i_e).c_str(),"L");

      }

      p_plot->cd();
      hs_Bias->Draw("nostack HIST");
      SetAxisFonts(hs_Bias);
      hs_Bias->SetMaximum(0.06);
      hs_Bias->SetMinimum(-0.5);
      c->Print(("Plots/"+var+"/Bias_"+gen+".pdf").c_str());  
      p_plot->Clear();

      hs_Variance->Draw("nostack HIST");
      SetAxisFonts(hs_Variance);
      hs_Variance->SetMaximum(0.27);
      hs_Variance->SetMinimum(0);
      c->Print(("Plots/"+var+"/Variance_"+gen+".pdf").c_str());  
      p_plot->Clear();
      l->Clear();

    }

  }

  // Distribution of each variable
  l->SetNColumns(4);
  for(size_t i_v=0;i_v<vars.size();i_v++){

    std::string axis_title = axis_titles.at(i_v);
    std::string y_axis_title = dist_axis_titles.at(i_v);
    std::string var = vars.at(i_v);

    THStack* hs = new THStack(("hs_"+var).c_str(),(";"+axis_title+";"+y_axis_title).c_str());

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      TH1D* h = h_var.at(i_f).at(i_v);
      h->Scale(1.0/h->Integral());
      DivideByBinWidth(h);
      h->SetLineColor(i_f+1);
      h->SetLineWidth(2);
      hs->Add(h);
      l->AddEntry(h,gen.c_str(),"L");

    }

    p_plot->cd();
    hs->Draw("HIST nostack");
    SetAxisFonts(hs);
    c->Print(("Plots/"+var+"/Dist.pdf").c_str());
    p_plot->Clear();

    l->Clear();

  }

  // Band and band width plots
  for(size_t i_v=0;i_v<vars.size();i_v++){

    std::string var = vars.at(i_v);
    std::string axis_title = axis_titles.at(i_v);

    // Calculate the difference in bias between the difference generators - calculate spread  
    std::vector<std::vector<TH1D*>> h_bias(kMAX,std::vector<TH1D*>(InputFiles_v.size()));

    for(size_t i_e=0;i_e<kMAX;i_e++){

      if(var == "W" && i_e == kSFMethod) continue;

      std::string est = estimators_str.at(i_e);

      THStack* hs_Bias = new THStack(("hs_Bias_"+var+"_"+est).c_str(),(";"+axis_title+";Frac. Bias").c_str());

      for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

        std::string gen = Generators_v.at(i_f);

        TH3D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f).at(i_v);
        h_bias.at(i_e).at(i_f) = new TH1D(("h_Bias2_"+gen+"_"+est+"_"+var).c_str(),"",points.at(i_v),binning_a.at(i_v));
        TH1D* h_variance = new TH1D(("h_Variance2_"+gen+"_"+est).c_str(),"",points.at(i_v),binning_a.at(i_v));
        MakeBiasVarianceFrom3D(h,h_bias.at(i_e).at(i_f),h_variance); 

        delete h;
        delete h_variance;

        h_bias.at(i_e).at(i_f)->SetLineColor(colors.at(i_e));
        h_bias.at(i_e).at(i_f)->SetLineStyle(i_f+1);
        h_bias.at(i_e).at(i_f)->SetLineWidth(2);
        hs_Bias->Add(h_bias.at(i_e).at(i_f));
        l->AddEntry(h_bias.at(i_e).at(i_f),gen.c_str(),"L");
               
      } 

      p_plot->cd();
      hs_Bias->Draw("nostack HIST");
      hs_Bias->SetMinimum(-0.62);
      hs_Bias->SetMaximum(0.12);
      SetAxisFonts(hs_Bias);
      c->Print(("Plots/"+var+"/Bias_Bands"+est+".pdf").c_str()); 
      p_plot->Clear();
      l->Clear();

    }

  }

}
