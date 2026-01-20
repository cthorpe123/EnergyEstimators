#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"
#include "../Funcs/PlotSetup.h"

OscModel data_osc_model;
OscModel fit_osc_model;

bool makeloadsofplots = true;

void DeltaCPPlot(){

  PlotSetup(); 

  // Load the histograms
  TFile* f = TFile::Open("ResponseMatricesNue.root");

  std::vector<std::string> Generators_v = {"GENIE","NuWro"};

  gSystem->Exec("mkdir -p Plots/");
  gSystem->Exec("mkdir -p Plots/FitPlots/");

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("../Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();
  h_flux->Scale(1.0/h_flux->Integral());

  for(size_t i_f=0;i_f<Generators_v.size();i_f++){
    std::string gen = Generators_v.at(i_f);    

    std::vector<TH1D*> h_fit_results(estimators_str.size());
    for(size_t i_e=0;i_e<estimators_str.size();i_e++)
      h_fit_results.at(i_e) = new TH1D(("h_fit_results_"+estimators_str.at(i_e)+"_"+gen).c_str(),";Input #delta_{CP} (rad);Measured #delta_{CP} - Input #delta_{CP} (rad)",20,-3.14/4,3.14/4);

    double min_fit_ratio = 0;
    double max_fit_ratio = 0;

    for(int i_dm=1;i_dm<h_fit_results.back()->GetNbinsX()+1;i_dm++){

      double true_deltaCP = h_fit_results.back()->GetBinCenter(i_dm);
      data_osc_model.SetDeltaCP(true_deltaCP);

      for(size_t i_e=0;i_e<estimators_str.size();i_e++){

        TH2D* h = static_cast<TH2D*>(f->Get((gen+"_TrueEnergy_RecoEnergy_"+estimators_str.at(i_e)).c_str())); 
        TH2D* h_nofsi = static_cast<TH2D*>(f->Get((gen+"_TrueEnergy_RecoEnergy_NoFSI_"+estimators_str.at(i_e)).c_str())); 

        if(i_dm == 1){
          Normalise(h);
          Normalise(h_nofsi);
        }

        if(i_e == kMuonKin) continue;

        const TH2D* h_data_true_reco = h;
        int nbins = h_data_true_reco->GetNbinsY();
        double low = h_data_true_reco->GetYaxis()->GetBinLowEdge(1);
        double high = h_data_true_reco->GetYaxis()->GetBinLowEdge(nbins+1);

        // Fold the 2D plot into a 1D plot
        TH1D* h_data_reco = new TH1D("h_data_reco",";E_{est} (GeV);Events",nbins,low,high); 
        for(int j=1;j<nbins+1;j++){
          double events = 0.0;
          for(int i=1;i<h_data_true_reco->GetNbinsX()+1;i++){
            double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
            double E = h_data_true_reco->GetXaxis()->GetBinCenter(i);
            if(E > 5.0 || E < 0.8) continue;
            events += data_osc_model.NueAppProb(E)*h_data_true_reco->GetBinContent(i,j)*flux;
          }
          h_data_reco->SetBinContent(j,events);
        }

        const TH2D* h_model_true_reco = h_nofsi;

        std::cout << "Fitting " << estimators_str.at(i_e) << " " << gen << std::endl;

        double meas_deltaCP = 0; 
        TH1D* h_model_reco = new TH1D("h_model_reco","",nbins,low,high);
        TH1D* h_model_reco_prefit = new TH1D("h_model_reco_prefit","",nbins,low,high);

        // Draw the model prediction with the input deltam2/amp 
        for(int j=1;j<nbins+1;j++){
          double events = 0.0;
          for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
            double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
            double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
            if(E > 5.0 || E < 0.8) continue;
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
            if(h_model_true_reco->GetXaxis()->GetBinCenter(i) < 0.8 || h_model_true_reco->GetXaxis()->GetBinCenter(i) > 5.0) continue;
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

        double limit_low = true_deltaCP-3.14/4;
        double limit_up = true_deltaCP+3.14/4;

        fMinimizer->SetVariableLimits(1,limit_low,limit_up);

        fMinimizer->SetFunction(min);
        fMinimizer->Minimize();

        meas_deltaCP = fMinimizer->X()[1];

        //if(abs(meas_deltaCP - limit_low) < 1e-5 || abs(meas_deltaCP - limit_up) < 1e-5) return false;
        //else return true;


/*

        // Set the function to be minimised in the fit
        ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff){

            fit_osc_model.SetDeltaMSq23(coeff[0]);
            fit_osc_model.SetNuMuDisAmp(coeff[1]);

            // Calculate the oscillated prediction
            for(int j=1;j<nbins+1;j++){
            double events = 0.0;
            for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
            double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
            double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
            if(E > 5.0 || E < 0.8) continue;
            events += fit_osc_model.NuMuSurvProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
            }
            h_model_reco->SetBinContent(j,events);
            }

            double chi2 = 0.0;
            for(int j=1;j<nbins+1;j++){
            if(h_data_reco->GetBinContent(j) > 0)
            chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/sqrt(h_data_reco->GetBinContent(j)),2);         
            }

            return chi2/nbins*100000;

        }, 2);

        std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
          ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

        fMinimizer->SetVariable(0,"deltamsq23",true_deltamsq23,0.01*c_osc::deltamsq23);
        fMinimizer->SetVariable(1,"amp",data_osc_model.numu_dis_amp,0.001*data_osc_model.numu_dis_amp);

        fMinimizer->SetFunction(min);
        fMinimizer->Minimize();

        meas_deltamsq23 = fMinimizer->X()[0];
*/







        if(makeloadsofplots){
          gSystem->Exec(("mkdir -p Plots/FitPlots/"+gen+"/"+estimators_str.at(i_e)).c_str());
          h_data_reco->Draw("HIST");
          h_data_reco->SetLineColor(1);
          h_data_reco->SetLineWidth(2);
          h_data_reco->SetStats(0);
          h_model_reco->Draw("HIST same");
          h_model_reco->SetLineColor(2);
          h_model_reco->SetLineWidth(2);
          h_model_reco_prefit->Draw("HIST same");
          h_model_reco_prefit->SetLineColor(2);
          h_model_reco_prefit->SetLineStyle(2);
          h_model_reco_prefit->SetLineWidth(2);
          l->AddEntry(h_data_reco,(gen+" FD").c_str(),"L");
          l->AddEntry(h_model_reco_prefit,(gen+" Model, No Fit").c_str(),"L");
          l->AddEntry(h_model_reco,(gen+" Model, Fitted").c_str(),"L");
          l->AddEntry((TObject*)0,("Input #delta_{CP}="+std::to_string(true_deltaCP)+" Measured #delta_{CP}="+std::to_string(meas_deltaCP)).c_str(),"");
          c->Print(("Plots/FitPlots/"+gen+"/"+estimators_str.at(i_e)+"/Point_" + std::to_string(i_dm) + "_" + gen + "_" + estimators_str.at(i_e) +".pdf").c_str());
          p_plot->Clear();
          l->Clear();
        } 

        h_fit_results.at(i_e)->SetBinContent(i_dm,meas_deltaCP-true_deltaCP);
        min_fit_ratio = std::min(min_fit_ratio,meas_deltaCP - true_deltaCP);
        max_fit_ratio = std::max(max_fit_ratio,meas_deltaCP - true_deltaCP);

        delete h_data_reco;
        delete h_model_reco;
        delete h_model_reco_prefit;

      }
    }

    TF1* f_line = new TF1("f_line","0",-1000,1000);
    f_line->SetLineColor(1);  
    f_line->SetLineWidth(2);
    f_line->SetLineStyle(9);

    THStack* hs = new THStack("hs",";Input #delta_{CP} (rad) ;Measured #delta_{CP} - Input #delta_{CP} (rad)");     

    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      if(i_e == kMuonKin) continue;
      h_fit_results.at(i_e)->SetLineColor(colors.at(i_e)); 
      h_fit_results.at(i_e)->SetLineWidth(2); 
      hs->Add(h_fit_results.at(i_e));
      l->AddEntry(h_fit_results.at(i_e),estimators_leg.at(i_e).c_str(),"L");
    }

    hs->Draw("HIST nostack"); 
    hs->SetMinimum(min_fit_ratio-0.1*(max_fit_ratio-min_fit_ratio));
    hs->SetMaximum(max_fit_ratio+0.1*(max_fit_ratio-min_fit_ratio));
    f_line->Draw("L same");
    hs->GetXaxis()->SetNdivisions(6);
    SetAxisFonts(hs);
    c->Print(("Plots/DeltaCPFitResults_"+gen+".pdf").c_str());
    p_plot->Clear();
    l->Clear();

    delete hs;
    delete f_line;

  }

}
