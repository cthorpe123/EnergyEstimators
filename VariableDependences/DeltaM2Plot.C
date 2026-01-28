#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"
#include "../Funcs/PlotSetup.h"

bool makeloadsofplots = true;

OscModel data_osc_model;
OscModel fit_osc_model;

std::pair<double,double> true_e_range = {0.8,5.0};

void DoFit(TH1D* h_flux,TH2D* h_data_true_reco,TH2D* h_model_true_reco,TH1D* h_data_reco,TH1D* h_model_reco,TH1D* h_model_reco_prefit,double true_deltamsq23,double& meas_deltamsq23){

  data_osc_model.SetDeltaMSq23(true_deltamsq23);

  Normalise(h_data_true_reco);
  Normalise(h_model_true_reco);

  int nbins = h_data_true_reco->GetNbinsY();
  double low = h_data_true_reco->GetYaxis()->GetBinLowEdge(1);
  double high = h_data_true_reco->GetYaxis()->GetBinLowEdge(nbins+1);

  // Fold the 2D data plot into a 1D plot
  for(int j=1;j<nbins+1;j++){
    double events = 0.0;
    for(int i=1;i<h_data_true_reco->GetNbinsX()+1;i++){
      if(h_data_true_reco->GetXaxis()->GetBinCenter(i) < true_e_range.first || h_data_true_reco->GetXaxis()->GetBinCenter(i) > true_e_range.second) continue;
      double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
      double E = h_data_true_reco->GetXaxis()->GetBinCenter(i);
      events += data_osc_model.NuMuSurvProb(E)*h_data_true_reco->GetBinContent(i,j)*flux;
    }
    h_data_reco->SetBinContent(j,events);
  }

  // Draw the model prediction with the input deltam2/amp 
  for(int j=1;j<nbins+1;j++){
    double events = 0.0;
    for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
      if(h_model_true_reco->GetXaxis()->GetBinCenter(i) < true_e_range.first || h_model_true_reco->GetXaxis()->GetBinCenter(i) > true_e_range.second) continue;
      double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
      double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
      events += data_osc_model.NuMuSurvProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
    }
    h_model_reco_prefit->SetBinContent(j,events);
  }

  // Set the function to be minimised in the fit
  ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff){

      fit_osc_model.SetDeltaMSq23(coeff[0]);
      fit_osc_model.SetNuMuDisAmp(coeff[1]);

      // Calculate the oscillated prediction
      for(int j=1;j<nbins+1;j++){
      double events = 0.0;
      for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
      if(h_model_true_reco->GetXaxis()->GetBinCenter(i) < true_e_range.first || h_model_true_reco->GetXaxis()->GetBinCenter(i) > true_e_range.second) continue;
      double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
      double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
      events += fit_osc_model.NuMuSurvProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
      }
      h_model_reco->SetBinContent(j,events);
      }

      h_model_reco->Scale(coeff[2]);

      double chi2 = 0.0;
      for(int j=1;j<nbins+1;j++){
        if(h_data_reco->GetBinContent(j) > 0)
          chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/sqrt(h_data_reco->GetBinContent(j)),2);         
      }

      //std::cout << coeff[0] << " " << coeff[1] << " chi2/nbins=" << chi2/nbins << std::endl;

      return chi2/nbins*1000000;

  }, 3);

  std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

  fMinimizer->SetVariable(0,"deltamsq23",true_deltamsq23,0.01*c_osc::deltamsq23);
  fMinimizer->SetVariableLimits(0,true_deltamsq23*0.8,true_deltamsq23*1.2);
  fMinimizer->SetVariable(1,"amp",data_osc_model.numu_dis_amp,0.01*data_osc_model.numu_dis_amp);
  fMinimizer->SetVariable(2,"scale",1.0,0.1);

  fMinimizer->SetFunction(min);
  fMinimizer->Minimize();

  meas_deltamsq23 = fMinimizer->X()[0];

}

void DeltaM2Plot(){

  PlotSetup(); 

  c->SetCanvasSize(800,420);

  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("../Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();
  h_flux->Scale(1.0/h_flux->Integral());

  std::string var = "W";

  TFile* f = TFile::Open("ResponseMatricesNuMu.root");

  double true_theta13 = c_osc::theta13;
  double true_theta23 = c_osc::theta23;

  double true_deltamsq23 = c_osc::deltamsq23;

  gSystem->Exec(("mkdir -p Plots/" + var).c_str());

  for(size_t i_f=0;i_f<Generators_v.size();i_f++){

    std::string gen = Generators_v.at(i_f);    

    double min_fit_ratio = 1;
    double max_fit_ratio = 1;

    //std::vector<TH1D*> h_fit_results;
    std::vector<TGraph*> g_fit_results(kMAX);

    std::string axis_title; 
    for(size_t i_e=0;i_e<estimators_str.size();i_e++){

      if(var == "W" && i_e == kSFMethod) continue;

      std::string est = estimators_str.at(i_e);
      std::cout << "Doing fitting with generator " << gen << " estimator " << est << std::endl;

      TH3D* h = static_cast<TH3D*>(f->Get((gen+"_TrueEnergy_RecoEnergy_"+var+"_"+est).c_str())); 

      int low_x_bin = h->FindBin(true_e_range.first);
      int high_x_bin = h->FindBin(true_e_range.second);
      h->GetXaxis()->SetRange(low_x_bin,high_x_bin);

      TH1D* h_fit_results = h->ProjectionZ();
      std::vector<Double_t> x_fit_results;
      std::vector<Double_t> y_fit_results;

      axis_title = h->GetZaxis()->GetTitle();

      for(int i_dm=1;i_dm<h_fit_results->GetNbinsX()+1;i_dm++){

        //if(i_e == kSFMethod /*|| i_e == kMuonKin*/) continue;

        // Pick a slice in the 3rd vairbale to fit
        h->GetZaxis()->SetRange(1,h->GetNbinsZ()); 
        TH2D* h_data_true_reco =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_data_true_reco")); 
        h->GetZaxis()->SetRange(1,i_dm); 
        TH2D* h_model_true_reco =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_model_true_reco")); 

        int nbins = h_data_true_reco->GetNbinsY();
        double low = h_data_true_reco->GetYaxis()->GetBinLowEdge(1);
        double high = h_data_true_reco->GetYaxis()->GetBinLowEdge(nbins+1);
        TH1D* h_data_reco = new TH1D("h_data_reco",";E_{est} (GeV);Events",nbins,low,high); 
        TH1D* h_model_reco = new TH1D("h_model_reco","E_{est} (GeV)",nbins,low,high);
        TH1D* h_model_reco_prefit = new TH1D("h_model_reco_prefit","E_{est} (GeV)",nbins,low,high);

        double meas_deltamsq23;
        DoFit(h_flux,h_data_true_reco,h_model_true_reco,h_data_reco,h_model_reco,h_model_reco_prefit,true_deltamsq23,meas_deltamsq23);

        if(makeloadsofplots){
          gSystem->Exec(("mkdir -p Plots/" + var + "/FitPlots/" + gen + "/" + est + "/").c_str());
          h_data_reco->Draw("HIST");
          h_data_reco->SetLineColor(1);
          h_data_reco->SetLineWidth(2);
          h_data_reco->SetStats(0);
          h_model_reco->Draw("HIST same");
          h_model_reco->SetLineColor(2);
          h_model_reco->SetLineWidth(2);
          h_model_reco_prefit->Draw("HIST same");
          h_model_reco_prefit->SetLineColor(3);
          h_model_reco_prefit->SetLineWidth(2);
          l->AddEntry(h_data_reco,(gen+" FD").c_str(),"L");
          l->AddEntry(h_model_reco_prefit,(gen+" Model, No Fit").c_str(),"L");
          l->AddEntry(h_model_reco,(gen+" Model, Fitted").c_str(),"L");
          l->AddEntry((TObject*)0,("Input #Delta m^{2}="+std::to_string(true_deltamsq23)+" Measured #Delta m^{2}="+std::to_string(meas_deltamsq23)).c_str(),"");
          std::string point;
          if(i_dm < 10) point = "00" + std::to_string(i_dm);
          else if(i_dm < 100) point = "0" + std::to_string(i_dm);
          else if(i_dm < 1000) point = std::to_string(i_dm);
          c->Print(("Plots/" + var + "/FitPlots/" + gen + "/" + est + "/" + "Point_" + point + "_" + gen + "_" + est +".png").c_str());
          p_plot->Clear();
          l->Clear();
        } 

        delete h_data_reco;
        delete h_model_reco;
        delete h_model_reco_prefit;

        x_fit_results.push_back(h_fit_results->GetBinCenter(i_dm));
        y_fit_results.push_back(meas_deltamsq23/true_deltamsq23);
        //h_fit_results.back()->SetBinContent(i_dm,meas_deltamsq23/true_deltamsq23);
        min_fit_ratio = std::min(min_fit_ratio,meas_deltamsq23/true_deltamsq23);
        max_fit_ratio = std::max(max_fit_ratio,meas_deltamsq23/true_deltamsq23);

      }

      g_fit_results.at(i_e) = new TGraph(x_fit_results.size(),&(x_fit_results[0]),&(y_fit_results[0]));
      delete h_fit_results;

    }

    std::string title = ";Cut on "+axis_title+"Measured #Delta m^{2}_{23}/Input #Delta m^{2}_{23}";
    //THStack* hs = new THStack("hs",title.c_str());     
    TMultiGraph* mg = new TMultiGraph("mg",title.c_str());

    TF1* f_line = new TF1("f_line","1",-1000,1000);
    f_line->SetLineColor(1);  
    f_line->SetLineWidth(2);
    f_line->SetLineStyle(9);

    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      if(var == "W" && i_e == kSFMethod) continue;
      g_fit_results.at(i_e)->SetLineColor(colors.at(i_e)); 
      g_fit_results.at(i_e)->SetLineWidth(2); 
      mg->Add(g_fit_results.at(i_e));
      l->AddEntry(g_fit_results.at(i_e),estimators_leg.at(i_e).c_str(),"L");
    }

    p_plot->cd();
    mg->Draw("AL"); 
    f_line->Draw("L same");
    SetAxisFontsMG(mg);
    mg->GetYaxis()->SetRangeUser(std::min(0.975,min_fit_ratio),std::max(1.025,max_fit_ratio));
    c->Print(("Plots/"+var+"/"+"DeltaM2FitResults_"+gen+".pdf").c_str());
    p_plot->Clear();
    l->Clear();

    delete mg;
    delete f_line;

  }

}

