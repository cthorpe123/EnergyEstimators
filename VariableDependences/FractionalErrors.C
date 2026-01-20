#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/Funcs.h"
#include "../Funcs/PlotSetup.h"
#include "TLorentzVector.h"

void FractionalErrors(){

  PlotSetup(); 

  bool rebin = false;

  std::vector<std::string> InputFiles_v = {"GENIEEventsFiltered.root","NEUTEventsFiltered.root","GiBUUEventsFiltered.root","NuWroEventsFiltered.root"};
  std::vector<std::string> Generators_v = {"GENIE","NEUT","GiBUU","NuWro"};
  std::vector<std::string> vars = {"W","Angle","MissingE"};
  std::vector<std::string> axis_titles = {"W_{vis} (GeV)","#theta_{lepton} (rad.)","Missing Hadronic Energy (GeV)"};
  std::vector<std::string> legs = {"W_{vis}","#theta_{l}","E_{miss}"};
  std::vector<std::string> units = {"GeV","","GeV"};

  std::vector<std::vector<std::vector<TH2D*>>> h_err(Generators_v.size(),std::vector<std::vector<TH2D*>>(kMAX,std::vector<TH2D*>(vars.size())));

  std::vector<int> points;
  std::vector<std::vector<double>> binning;
  std::vector<double*> binning_a;
  std::vector<std::pair<double,double>> bias_ranges;

  // W histogram setup
  std::vector<double> w_binning = {0.8,1.0,1.5,2.0,2.5,3.0,4.0,5.0}; 
  int w_points = w_binning.size()-1;
  double* w_binning_a = &w_binning[0]; 
  points.push_back(w_points);
  binning.push_back(w_binning);
  binning_a.push_back(w_binning_a);
  bias_ranges.push_back(std::make_pair(-0.6,0.1));

  // Angle histogram setup
  std::vector<double> ang_binning = {0.0,0.2,0.4,0.6,0.9,1.2,1.8,3.142}; 
  int ang_points = ang_binning.size()-1;
  double* ang_binning_a = &ang_binning[0]; 
  points.push_back(ang_points);
  binning.push_back(ang_binning);
  binning_a.push_back(ang_binning_a);
  bias_ranges.push_back(std::make_pair(-0.65,0.1));

  // Missing energy histogram setup
  std::vector<double> miss_binning = {0.0,0.05,0.1,0.2,0.4,0.6,1.2}; 
  int miss_points = miss_binning.size()-1;
  double* miss_binning_a = &miss_binning[0]; 
  points.push_back(miss_points);
  binning.push_back(miss_binning);
  binning_a.push_back(miss_binning_a);
  bias_ranges.push_back(std::make_pair(-0.6,0.0));

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
        std::string est = estimators_str.at(i_e);
        h_err.at(i_f).at(i_e).at(i_v) = new TH2D((gen+"_frac_err_"+var+"_"+est).c_str(),"",100,-1,1,points.at(i_v),binning_a.at(i_v));
      }
    }


    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 20000 == 0) std::cout << gen << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(gen != "GiBUU") weight = 1.0;
      //if(gen == "GiBUU" && weight > 1) continue;

      if(nu_pdg != 14 || ccnc != 1) continue;
      if(nprot < 1) continue;

      double angle = acos(lepton_p4->CosTheta());
      double missing_e = GetMissingEnergy(pdg,p4); 
      int neutrons = 0;

      for(int i_p=0;i_p<pdg->size();i_p++)
        if(pdg->at(i_p) == 2112 && p4->at(i_p).Vect().Mag() > 0.3)
          neutrons++;

      std::vector<double> values = {W,angle,missing_e};

      for(size_t i_v=0;i_v<vars.size();i_v++){
        double val = values.at(i_v);
        for(int i_e=0;i_e<kMAX;i_e++){
          double nu_e_reco = est_nu_e->at(i_e);
          if(nu_e_reco < 0) continue;
          h_err.at(i_f).at(i_e).at(i_v)->Fill((nu_e_reco-nu_e)/nu_e,val,weight);  
        } // i_e
      } // i_v

    } // ievent

  } // i_f



  gSystem->Exec("mkdir -p Plots/");

  // Bias AFO difference variables

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
    std::string gen = Generators_v.at(i_f);
    for(size_t i_e=0;i_e<kMAX;i_e++){
      std::string est = estimators_str.at(i_e);
      for(size_t i_v=0;i_v<vars.size();i_v++){
        std::string var = vars.at(i_v);

        gSystem->Exec(("mkdir -p Plots/" + var).c_str());

        if(var == "W" && i_e == kSFMethod) continue;

        const TH2D* h = h_err.at(i_f).at(i_e).at(i_v);
        std::vector<TH1D*> h_slices(h->GetNbinsY());
        std::string name = h->GetName();

        l->SetNColumns(3);

        THStack* hs = new THStack("hs",";Fractional Bias;Events");

        // Cut the 2D hist into 1D slices
        for(int i=1;i<h->GetNbinsY()+1;i++){
          h_slices.at(i-1) = h->ProjectionX((name+"_"+std::to_string(i)).c_str(),i,i);
          h_slices.at(i-1)->SetLineWidth(2);    
          h_slices.at(i-1)->SetLineColor(i);
          h_slices.at(i-1)->Scale(1.0/h_slices.at(i-1)->Integral());
          std::string leg = to_string_with_precision(h->GetYaxis()->GetBinLowEdge(i),2) + " < " + legs.at(i_v) + " < " + to_string_with_precision(h->GetYaxis()->GetBinLowEdge(i),2) + " " + units.at(i_v);
          hs->Add(h_slices.at(i-1));
          l->AddEntry(h_slices.at(i-1),leg.c_str(),"L");
        }

        p_plot->cd();
        hs->Draw("nostack HIST");
        SetAxisFonts(hs);
        c->Print(("Plots/"+var+"/FracBias_"+gen+"_"+est+".pdf").c_str());
        p_plot->Clear();
        l->Clear();

        delete hs;

      }


    }

  }


}
