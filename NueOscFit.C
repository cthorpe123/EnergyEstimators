#include "Funcs/Funcs.h"

double app_prob(double E,double test_deltaCP,double test_theta12,double test_theta23,double test_theta13,double test_deltamsq13,double test_deltamsq12){

  double delta31 = 1.267*test_deltamsq13*L/E;
  double delta21 = 1.267*test_deltamsq12*L/E;

  double P = sin(test_theta23)*sin(test_theta23)*sin(2*test_theta13)*sin(2*test_theta13)*pow(sin(delta31 - aL),2)*delta31*delta31/pow(delta31 - aL,2)
    + sin(2*test_theta23)*sin(2*test_theta13)*sin(2*test_theta12)*sin(delta31 - aL)*delta31/(delta31 - aL)
    *(sin(aL)/aL)*delta21*cos(delta31 + test_deltaCP)
    + cos(test_theta23)*cos(test_theta23)*sin(2*test_theta12)*sin(2*test_theta12)*pow(sin(aL),2)/aL/aL*delta21*delta21; 

  return P;

}

// Use the first generator in the list as a model, and treat the
// other three as data

const double exposure = 5*40;

void NueOscFit(){

  std::vector<std::string> Generators_v = {"NuWro","NEUT","GiBUU","GENIE"};

  // Get all of the true vs reco spectra

  TFile* f_in = TFile::Open("rootfiles/NueRates.root");
  std::vector<std::vector<TH2D*>> h_energy_true_reco;
  for(std::string estimator : estimators_str){
    h_energy_true_reco.push_back(std::vector<TH2D*>());
    for(std::string generator : Generators_v){
      h_energy_true_reco.back().push_back(static_cast<TH2D*>(f_in->Get((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str())));
    }
  }

  gSystem->Exec("mkdir -p Plots/NueFitPlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  std::vector<std::vector<TH1D*>> h_fit_results(Generators_v.size(),std::vector<TH1D*>());
  std::vector<std::vector<TH1D*>> h_frac_error(Generators_v.size(),std::vector<TH1D*>());

  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    TH2D* h_model_true_reco = static_cast<TH2D*>(h_energy_true_reco.at(i_e).at(0)->Clone("h_model_true_reco"));

    for(size_t i_g=1;i_g<Generators_v.size();i_g++){

      std::cout << "Fitting " << estimators_str.at(i_e) << " " << Generators_v.at(i_g) << std::endl;

      TH2D* h_data_true_reco = static_cast<TH2D*>(h_energy_true_reco.at(i_e).at(i_g)->Clone("h_data_true_reco"));
      int nbins = h_data_true_reco->GetNbinsY();
      double low = h_data_true_reco->GetYaxis()->GetBinLowEdge(1);
      double high = h_data_true_reco->GetYaxis()->GetBinLowEdge(nbins+1);
      
      double true_deltaCP = 0.0;
      double true_deltamsq12 = deltamsq12;
      double true_deltamsq13 = deltamsq13;
      double true_theta12 = theta12;
      double true_theta13 = theta13;
      double true_theta23 = theta23;

      TH1D* h_data_reco = new TH1D("h_data_reco","",nbins,low,high);

      //std::cout << "Data Osc Parameters: true_deltamsq23=" << true_deltamsq23 << " true_amp=" << true_amp << std::endl;

      // Fold the 2D plot into a 1D plot
      for(int j=1;j<nbins+1;j++){
        double events = 0.0;
        for(int i=1;i<h_data_true_reco->GetNbinsX()+1;i++){
          double p = app_prob(h_data_true_reco->GetXaxis()->GetBinCenter(i),true_deltaCP,true_theta12,true_theta23,true_theta13,true_deltamsq13,true_deltamsq12);
          events += exposure*p*h_data_true_reco->GetBinContent(i,j);
          //std::cout << exposure << "  " << surv_prob(h_data_reco->GetXaxis()->GetBinCenter(i),true_deltamsq23,true_amp) << " " << h_data_true_reco->GetBinContent(i,j) << std::endl;
        }
        //std::cout << j << "  " << events << std::endl;
        h_data_reco->SetBinContent(j,events);
        h_data_reco->SetBinError(j,sqrt(events));
      }

      TH1D* h_model_reco = new TH1D("h_model_reco","",nbins,low,high);

      // Set the function to be minimised in the fit
      ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff ){

          // Calculate the oscillated prediction
          for(int j=1;j<nbins+1;j++){
            double events = 0.0;
            for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
            //double p = app_prob(h_model_true_reco->GetXaxis()->GetBinCenter(i),coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]);
            double p = app_prob(h_model_true_reco->GetXaxis()->GetBinCenter(i),coeff[0],true_theta12,coeff[1],true_theta13,true_deltamsq13,true_deltamsq12);
              events += exposure*p*h_model_true_reco->GetBinContent(i,j);
            }
            h_model_reco->SetBinContent(j,events);
          }

          double chi2 = 0.0;
          for(int j=1;j<nbins+1;j++){
            if(h_data_reco->GetBinContent(j) > 0)
              chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/h_data_reco->GetBinError(j),2);         
          }

          //std::cout << "chi2/nbins = " << chi2/nbins << std::endl; 
          return chi2/nbins;

          }, 2);

      std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
        ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

      fMinimizer->SetVariable(0,"deltaCP",0.0,0.01);
      fMinimizer->SetVariable(1,"theta23",0.0,0.01);
/*      fMinimizer->SetVariable(1,"theta12",theta12,0.01);
      fMinimizer->SetVariable(2,"theta23",theta23,0.01);
      fMinimizer->SetVariable(3,"theta13",theta13,0.01);
      fMinimizer->SetVariable(4,"deltamsq13",deltamsq13,0.01);
      fMinimizer->SetVariable(5,"deltamsq12",deltamsq12,0.01);*/

      fMinimizer->SetFunction(min);

      std::cout << "Fitting..." << std::endl;
      fMinimizer->Minimize();

      std::cout << "Best Fit" << std::endl;
      std::cout << "deltaCP = " << fMinimizer->X()[0] << " +/- " << sqrt(fMinimizer->CovMatrix(0,0)) << std::endl;
      std::cout << "theta23 = " << fMinimizer->X()[1] << " +/- " << sqrt(fMinimizer->CovMatrix(1,1)) << std::endl;

      h_data_reco->Draw("e1");
      h_data_reco->SetLineColor(1);
      h_data_reco->SetLineWidth(2);
      h_data_reco->SetStats(0);
      h_model_reco->Draw("HIST same");
      h_model_reco->SetLineColor(2);
      h_model_reco->SetLineWidth(2);
      c->Print(("Plots/NueFitPlots/"+Generators_v.at(i_g) + "_" + estimators_str.at(i_e) +".png").c_str());
      c->Clear();

      delete h_model_reco;
      delete h_data_reco;
/*
      h_fit_results.at(i_g).push_back(new TH1D((Generators_v.at(i_g) + "_" + estimators_str.at(i_e)).c_str(),"",2,0,2));
      h_fit_results.at(i_g).back()->SetBinContent(1,fMinimizer->X()[0]/true_deltamsq23);
      h_fit_results.at(i_g).back()->SetBinContent(2,fMinimizer->X()[1]/true_amp);
      h_fit_results.at(i_g).back()->SetBinError(1,sqrt(fMinimizer->CovMatrix(0,0))/true_deltamsq23);
      h_fit_results.at(i_g).back()->SetBinError(2,sqrt(fMinimizer->CovMatrix(1,1))/true_amp);
      h_fit_results.at(i_g).back()->SetLineColor(i_e+1); 
      h_fit_results.at(i_g).back()->SetLineWidth(2); 

      h_frac_error.at(i_g).push_back(new TH1D(("frac_error_"+Generators_v.at(i_g) + "_" + estimators_str.at(i_e)).c_str(),"",2,0,2));
      h_frac_error.at(i_g).back()->SetBinContent(1,sqrt(fMinimizer->CovMatrix(0,0))/fMinimizer->X()[0]);
      h_frac_error.at(i_g).back()->SetBinContent(2,sqrt(fMinimizer->CovMatrix(1,1))/fMinimizer->X()[1]);
      h_frac_error.at(i_g).back()->SetLineColor(i_e+1); 
      h_frac_error.at(i_g).back()->SetLineWidth(2); 
*/
    }

    delete h_model_true_reco;

  }


 /* 
  for(size_t i_g=1;i_g<Generators_v.size();i_g++){

    THStack* hs = new THStack("hs",";;Fit/Input");    

    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      hs->Add(h_fit_results.at(i_g).at(i_e));
      l->AddEntry(h_fit_results.at(i_g).at(i_e),estimators_str.at(i_e).c_str(),"L");
    }

    hs->Draw("e1 nostack");
    hs->GetXaxis()->SetBinLabel(1,"#Delta m^{2}_{23}");
    hs->GetXaxis()->SetBinLabel(2,"4cos^{2}(#theta_{13})sin^{2}(#theta_{23})(1-cos^{2}(#theta_{13})sin^{2}(#theta_{23}))");
    hs->SetMinimum(0.8);
    hs->SetMaximum(1.2);
    c->Update();

    l->Draw();
    c->Print(("Plots/NumuFitPlots/FitComp_"+Generators_v.at(i_g)+".png").c_str());
    c->Clear();
    l->Clear();

    delete hs;
  } 


  for(size_t i_g=1;i_g<Generators_v.size();i_g++){

    THStack* hs = new THStack("hs",";;Frac Error In Fit");    

    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      hs->Add(h_frac_error.at(i_g).at(i_e));
      l->AddEntry(h_frac_error.at(i_g).at(i_e),estimators_str.at(i_e).c_str(),"L");
    }

    hs->Draw("HIST nostack");
    hs->GetXaxis()->SetBinLabel(1,"#Delta m^{2}_{23}");
    hs->GetXaxis()->SetBinLabel(2,"4cos^{2}(#theta_{13})sin^{2}(#theta_{23})(1-cos^{2}(#theta_{13})sin^{2}(#theta_{23}))");
    l->Draw();

    c->Print(("Plots/NumuFitPlots/FracError_"+Generators_v.at(i_g)+".png").c_str());
    c->Clear();
    l->Clear();

    delete hs;

  } 
*/

}
