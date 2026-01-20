#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"
#include "../Funcs/PlotSetup.h"

bool makeloadsofplots = true;

OscModel data_osc_model;
OscModel fit_osc_model;

std::pair<double,double> true_e_range = {0.8,5.0};

bool DoFit(TH1D* h_flux,TH2D* h_data_true_reco,TH2D* h_model_true_reco,TH1D* h_data_reco,TH1D* h_model_reco,TH1D* h_model_reco_prefit,double true_deltaCP,double& meas_deltaCP){

  data_osc_model.SetDeltaCP(true_deltaCP);

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
      events += data_osc_model.NueAppProb(E)*h_data_true_reco->GetBinContent(i,j)*flux;
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
      events += data_osc_model.NueAppProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
    }
    h_model_reco_prefit->SetBinContent(j,events);
  }

  // Set the function to be minimised in the fit
  ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff){

      fit_osc_model.SetDeltaCP(coeff[1]);
      //fit_osc_model.SetDeltaMSq13(coeff[2]);
      //fit_osc_model.SetTheta13(coeff[2]);

      // Calculate the oscillated prediction
      for(int j=1;j<nbins+1;j++){
      double events = 0.0;
      for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
      if(h_model_true_reco->GetXaxis()->GetBinCenter(i) < true_e_range.first || h_model_true_reco->GetXaxis()->GetBinCenter(i) > true_e_range.second) continue;
      double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
      double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
      events += fit_osc_model.NueAppProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
      }
      h_model_reco->SetBinContent(j,events);
      }

      h_model_reco->Scale(coeff[0]);

      double chi2 = 0.0;
      for(int j=1;j<nbins+1;j++){
        if(h_data_reco->GetBinContent(j) > 0)
          chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/sqrt(h_data_reco->GetBinContent(j)),2);         
      }

      //std::cout << coeff[0] << " " << coeff[1] << " chi2/nbins=" << chi2/nbins << std::endl;

      return chi2/nbins*1000000;

  }, 2);

  std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

  fMinimizer->SetMaxFunctionCalls(5000);

  fMinimizer->SetVariable(0,"scale",1.0,0.1);
  fMinimizer->SetVariable(1,"deltaCP",meas_deltaCP,0.01);
  //fMinimizer->SetVariable(2,"deltamsq13",c_osc::deltamsq13,0.01*c_osc::deltamsq13);
  //fMinimizer->SetVariable(2,"theta13",c_osc::theta13,0.01);

  double limit_low = true_deltaCP-3.14/2;
  double limit_up = true_deltaCP+3.14/2;

  fMinimizer->SetVariableLimits(1,limit_low,limit_up);

  fMinimizer->SetFunction(min);
  fMinimizer->Minimize();

  meas_deltaCP = fMinimizer->X()[1];

  if(abs(meas_deltaCP - limit_low) < 1e-5 || abs(meas_deltaCP - limit_up) < 1e-5) return false;
  else return true;

}

void DeltaCPPlot(){

  PlotSetup(); 

  c->SetCanvasSize(800,420);

  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("../Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();
  h_flux->Scale(1.0/h_flux->Integral());

  std::string var = "Angle";

  TFile* f = TFile::Open("ResponseMatricesNue.root");

  std::vector<double> true_deltaCP_v = {0,3.14/2,-3.14/2};
  std::vector<std::string> labels = {"Zero","Plus","Minus"};

  gSystem->Exec(("mkdir -p Plots/" + var).c_str());

  for(size_t i_dcp=0;i_dcp<true_deltaCP_v.size();i_dcp++){

    double true_deltaCP = true_deltaCP_v.at(i_dcp);

    for(size_t i_f=0;i_f<Generators_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);    

      double min_fit_ratio = 0;
      double max_fit_ratio = 0;

      //std::vector<TH1D*> h_fit_results;
      std::vector<TGraph*> g_fit_results(kMAX);

      std::string axis_title; 
      for(size_t i_e=0;i_e<estimators_str.size();i_e++){

        if(var == "W" && i_e == kSFMethod) continue;

        std::string est = estimators_str.at(i_e);
        std::cout << "Doing fitting with generator " << gen << " estimator " << est << std::endl;

        TH3D* h = static_cast<TH3D*>(f->Get((gen+"_TrueEnergy_RecoEnergy_"+var+"_"+est).c_str())); 
        //h->RebinX();
        //h->RebinY();
        //if(var != "Neutrons") h->RebinZ();

        TH1D* h_fit_results = h->ProjectionZ();
        std::vector<Double_t> x_fit_results;
        std::vector<Double_t> y_fit_results;

        //h_fit_results.push_back(h->ProjectionZ());
        //h_fit_results.back()->Reset();
        axis_title = h->GetZaxis()->GetTitle();

        //int low_x_bin = h->FindBin(true_e_range.first);
        //int high_x_bin = h->FindBin(true_e_range.second);
        //h->GetXaxis()->SetRange(low_x_bin,high_x_bin);

        double last_meas_deltaCP = true_deltaCP;
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

          double meas_deltaCP = last_meas_deltaCP;
          bool fit = DoFit(h_flux,h_data_true_reco,h_model_true_reco,h_data_reco,h_model_reco,h_model_reco_prefit,true_deltaCP,meas_deltaCP);

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
            //h_data_reco->SetTitle(("Input #delta_{CP}="+std::to_string(true_deltaCP)+" Measured #delta_{CP}="+std::to_string(meas_deltaCP)).c_str());
            l->AddEntry(h_data_reco,(gen+" FD").c_str(),"L");
            l->AddEntry(h_model_reco_prefit,(gen+" Model, No Fit").c_str(),"L");
            l->AddEntry(h_model_reco,(gen+" Model, Fitted").c_str(),"L");
            l->AddEntry((TObject*)0,("Input #delta_{CP}="+std::to_string(true_deltaCP)+" Measured #delta_{CP}="+std::to_string(meas_deltaCP)).c_str(),"");
            c->Print(("Plots/" + var + "/FitPlots/" + gen + "/" + est + "/" + "Point_" + std::to_string(i_dm) + "_" + gen + "_" + est +".pdf").c_str());
            p_plot->Clear();
            l->Clear();
          } 

          delete h_data_reco;
          delete h_model_reco;
          delete h_model_reco_prefit;

          if(fit){
            //h_fit_results.back()->SetBinContent(i_dm,meas_deltaCP - true_deltaCP);
            x_fit_results.push_back(h_fit_results->GetBinCenter(i_dm));
            y_fit_results.push_back(meas_deltaCP - true_deltaCP);
            min_fit_ratio = std::min(min_fit_ratio,meas_deltaCP - true_deltaCP);
            max_fit_ratio = std::max(max_fit_ratio,meas_deltaCP - true_deltaCP);
            last_meas_deltaCP = meas_deltaCP;
          }

        }

        g_fit_results.at(i_e) = new TGraph(x_fit_results.size(),&(x_fit_results[0]),&(y_fit_results[0]));

        delete h_fit_results;
        delete h;

      }

      std::string title = ";"+axis_title+"Measured #delta_{CP} - Input #delta_{CP} (rad)";
      //THStack* hs = new THStack("hs",title.c_str());     
      TMultiGraph* mg = new TMultiGraph("mg",title.c_str());

      TF1* f_line = new TF1("f_line","0",-1000,1000);
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

      mg->Draw("AL"); 
      mg->GetYaxis()->SetRangeUser(std::min(-0.05,min_fit_ratio),std::max(0.05,max_fit_ratio));
      f_line->Draw("L same");
      SetAxisFontsMG(mg);
      
      c->Print(("Plots/"+var+"/"+"DeltaCPFitResults_"+labels.at(i_dcp)+"_"+gen+".pdf").c_str());
      p_plot->Clear();
      l->Clear();

      delete mg;
      delete f_line;

    }

  }

}

