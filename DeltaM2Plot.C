#include "Funcs/Funcs.h"

// NuMu survival prob

////////////////////////////////////////////////////////////////////////////////
// numu disappearance probability
// equation 1 from https://www.nature.com/articles/s41586-020-2177-0 
// Patrick Dunne (Patrick DUNE?) tells me this is good approximation at DUNE's baseline/matter effect

// sintheta23 and costheta13 are degenerate in this situation, just fit a single "amplitude" parameter
double amp = 4*cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23)*(1-cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23));
double surv_prob(double E,double test_deltamsq23=deltamsq23,double test_amp=amp){
  double delta23 = 1.267*test_deltamsq23*L/E;
  return 1 - test_amp*sin(delta23)*sin(delta23);
}

// Use the first generator in the list as a model of the mapping
// from true neutrino energy into reco neutrino energy
// then use the responses from the other generators to make the reco
// distribution and fit deltam2 to each

void DeltaM2Plot(){

  gSystem->Exec("mkdir -p Plots/DeltaM2Plots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  std::vector<std::string> Generators_v = {"GiBUU","GENIE","NuWro","NEUT"};

  // Get all of the true vs reco spectra

  TFile* f_in = TFile::Open("rootfiles/NuMuResponseMatrices.root");
  std::vector<std::vector<TH2D*>> h_energy_true_reco;
  for(std::string estimator : estimators_str){
    h_energy_true_reco.push_back(std::vector<TH2D*>());
    for(std::string generator : Generators_v){
      h_energy_true_reco.back().push_back(static_cast<TH2D*>(f_in->Get((generator+"_TrueEnergy_RecoEnergy_"+estimator).c_str())));
    }
  }

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();
  h_flux->Scale(1.0/h_flux->Integral());

  double true_theta13 = theta13;
  double true_theta23 = theta23;
  double true_amp = 4*cos(true_theta13)*cos(true_theta13)*sin(true_theta23)*sin(true_theta23)*(1-cos(true_theta13)*cos(true_theta13)*sin(true_theta23)*sin(true_theta23));

  std::vector<std::vector<TH1D*>> h_fit_results(estimators_str.size(),std::vector<TH1D*>());
  for(size_t i_e=0;i_e<estimators_str.size();i_e++)
    for(size_t i_g=0;i_g<Generators_v.size();i_g++)
      h_fit_results.at(i_e).push_back(new TH1D(("h_fit_results_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_g)).c_str(),";Input #Delta m^{2}_{23};Measured #Delta m^{2}_{23}/Input #Delta m^{2}_{23}",50,0.9*deltamsq23,1.1*deltamsq23));

  for(int i_dm=1;i_dm<h_fit_results.back().back()->GetNbinsX()+1;i_dm++){

    double true_deltamsq23 = h_fit_results.back().back()->GetBinCenter(i_dm);

    for(size_t i_e=0;i_e<estimators_str.size();i_e++){

      TH2D* h_data_true_reco = static_cast<TH2D*>(h_energy_true_reco.at(i_e).at(0)->Clone("h_data_true_reco"));
      int nbins = h_data_true_reco->GetNbinsY();
      double low = h_data_true_reco->GetYaxis()->GetBinLowEdge(1);
      double high = h_data_true_reco->GetYaxis()->GetBinLowEdge(nbins+1);


      TH1D* h_data_reco = new TH1D("h_data_reco",";Estimated Neutrino Energy;Events",nbins,low,high); 
      // Fold the 2D plot into a 1D plot
      for(int j=1;j<nbins+1;j++){
        double events = 0.0;
        for(int i=1;i<h_data_true_reco->GetNbinsX()+1;i++){
          double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
          events += surv_prob(h_data_true_reco->GetXaxis()->GetBinCenter(i),true_deltamsq23,true_amp)*h_data_true_reco->GetBinContent(i,j)*flux;
        }
        //std::cout << j << "  " << events << std::endl;
        h_data_reco->SetBinContent(j,events);
        //h_data_reco->SetBinError(j,sqrt(events));
      }

      h_data_reco->Scale(1.0/h_data_reco->Integral());

      for(size_t i_g=0;i_g<Generators_v.size();i_g++){

        std::cout << "Fitting " << estimators_str.at(i_e) << " " << Generators_v.at(i_g) << std::endl;

        TH2D* h_model_true_reco = static_cast<TH2D*>(h_energy_true_reco.at(i_e).at(i_g)->Clone("h_model_true_reco"));
        TH1D* h_model_reco = new TH1D("h_model_reco","",nbins,low,high);

        // Set the function to be minimised in the fit
        ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff ){

            // Calculate the oscillated prediction
            for(int j=1;j<nbins+1;j++){
            double events = 0.0;
            for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
            double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
            events += surv_prob(h_model_true_reco->GetXaxis()->GetBinCenter(i),coeff[0],true_amp)*h_model_true_reco->GetBinContent(i,j)*flux;
            }
            h_model_reco->SetBinContent(j,events);
            }
            h_model_reco->Scale(1.0/h_model_reco->Integral());

            double chi2 = 0.0;
            for(int j=1;j<nbins+1;j++){
            if(h_data_reco->GetBinContent(j) > 0)
            chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/*/sqrt(h_data_reco->GetBinContent(j))*/,2);         
            }

            //if(i_g == 0) std::cout << "chi2/nbins = " << chi2/nbins << std::endl; 
            return chi2/nbins*10000;

        }, 1);

        std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
          ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

        fMinimizer->SetVariable(0,"deltamsq23",deltamsq23,0.001*deltamsq23);
        fMinimizer->SetVariableLimits(0,0.5*deltamsq23,1.5*deltamsq23); 
 
        //fMinimizer->SetVariable(1,"amp",true_amp,0.01*true_amp);
        //fMinimizer->SetVariable(2,"theta23",theta23,0.01*theta23);
        fMinimizer->SetFunction(min);
        fMinimizer->Minimize();

        //std::cout << "Data Osc Parameters: true_deltamsq23=" << true_deltamsq23 << " true_amp=" << true_amp << std::endl;
        //std::cout << "Fit  Osc Parameters:  fit_deltamsq23=" << fMinimizer->X()[0] << "  fit_amp=" << fMinimizer->X()[1] << std::endl;
        /*
           h_data_reco->Draw("HIST");
           h_data_reco->SetLineColor(1);
           h_data_reco->SetLineWidth(2);
           h_data_reco->SetStats(0);
           h_model_reco->Draw("HIST same");
           h_model_reco->SetLineColor(2);
           h_model_reco->SetLineWidth(2);
           c->Print(("Plots/DeltaM2Plots/"+Generators_v.at(i_g) + "_" + estimators_str.at(i_e) +".png").c_str());
           c->Clear();
           */
        /*

           h_frac_error.at(i_g).push_back(new TH1D(("frac_error_"+Generators_v.at(i_g) + "_" + estimators_str.at(i_e)).c_str(),"",2,0,2));
           h_frac_error.at(i_g).back()->SetBinContent(1,sqrt(fMinimizer->CovMatrix(0,0))/fMinimizer->X()[0]);
           h_frac_error.at(i_g).back()->SetBinContent(2,sqrt(fMinimizer->CovMatrix(1,1))/fMinimizer->X()[1]);
           h_frac_error.at(i_g).back()->SetLineColor(i_e+1); 
           h_frac_error.at(i_g).back()->SetLineWidth(2); 
           */

        if(i_g == 0) std::cout << fMinimizer->X()[0]/true_deltamsq23 << std::endl;

        h_fit_results.at(i_e).at(i_g)->SetBinContent(i_dm,fMinimizer->X()[0]/true_deltamsq23);

        delete h_model_reco;
        delete h_model_true_reco;
      }

      delete h_data_reco;
      delete h_data_true_reco;

    }

  }


  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    THStack* hs = new THStack("hs",";Input #Delta m^{2}_{23};Measured #Delta m^{2}_{23}/Input #Delta m^{2}_{23}");     

    for(size_t i_g=0;i_g<Generators_v.size();i_g++){
      h_fit_results.at(i_e).at(i_g)->SetLineColor(i_g+1); 
      h_fit_results.at(i_e).at(i_g)->SetLineWidth(2); 
      hs->Add(h_fit_results.at(i_e).at(i_g));
      l->AddEntry(h_fit_results.at(i_e).at(i_g),Generators_v.at(i_g).c_str(),"L");
    }

    hs->Draw("HIST nostack"); 
    hs->SetMinimum(0.9);
    hs->SetMaximum(1.1);
    l->Draw();
    c->Print(("Plots/DeltaM2Plots/FitResults_" + estimators_str.at(i_e) +".png").c_str());
    c->Clear();
    l->Clear();
}


}
