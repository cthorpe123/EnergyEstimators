#include "Funcs/Funcs.h"
#include "Funcs/OscFitter.h"

OscModel data_osc_model;
OscModel fit_osc_model;

// Use the first generator in the list as a model of the mapping
// from true neutrino energy into reco neutrino energy
// then use the responses from the other generators to make the reco
// distribution and fit deltam2 to each

bool makeloadsofplots = true;

void DeltaCPPlot(){

  gSystem->Exec("mkdir -p Plots/DeltaCPPlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  std::vector<std::string> Generators_v = {"GENIE","GiBUU","NEUT","NuWro"};

  // Get all of the true vs reco spectra

  TFile* f_in = TFile::Open("rootfiles/NueResponseMatrices.root");
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

  std::vector<std::vector<TH1D*>> h_fit_results(estimators_str.size(),std::vector<TH1D*>());
  for(size_t i_e=0;i_e<estimators_str.size();i_e++)
    for(size_t i_g=0;i_g<Generators_v.size();i_g++)
      h_fit_results.at(i_e).push_back(new TH1D(("h_fit_results_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_g)).c_str(),";Input #delta_{CP};Measured #delta_{CP}/Input #delta_{CP}",50,-3.142/2,3.142/2));

  for(int i_dm=1;i_dm<h_fit_results.back().back()->GetNbinsX()+1;i_dm++){

    double true_deltaCP = h_fit_results.back().back()->GetBinCenter(i_dm);
    data_osc_model.SetDeltaCP(true_deltaCP);

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
          double E = h_data_true_reco->GetXaxis()->GetBinCenter(i);
          events += data_osc_model.NueAppProb(E)*h_data_true_reco->GetBinContent(i,j)*flux;
        }
        h_data_reco->SetBinContent(j,events);
      }

      for(size_t i_g=0;i_g<Generators_v.size();i_g++){

        std::cout << "Fitting " << estimators_str.at(i_e) << " " << Generators_v.at(i_g) << std::endl;

        TH2D* h_model_true_reco = static_cast<TH2D*>(h_energy_true_reco.at(i_e).at(i_g)->Clone("h_model_true_reco"));
        TH1D* h_model_reco = new TH1D("h_model_reco","",nbins,low,high);

        double meas_deltaCP = 0; 

        // Set the function to be minimised in the fit
        ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff){

          fit_osc_model.SetDeltaCP(coeff[0]);
          fit_osc_model.SetDeltaMSq13(coeff[1]);
          fit_osc_model.SetTheta13(coeff[2]);

        // Calculate the oscillated prediction
        for(int j=1;j<nbins+1;j++){
        double events = 0.0;
        for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
        double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
        double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
        events += fit_osc_model.NueAppProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
        }
        h_model_reco->SetBinContent(j,events);
        }

        double chi2 = 0.0;
        for(int j=1;j<nbins+1;j++){
        if(h_data_reco->GetBinContent(j) > 0)
        chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/sqrt(h_data_reco->GetBinContent(j)),2);         
        }

        return chi2/nbins*100000;

        }, 3);

        std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
        ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

        fMinimizer->SetVariable(0,"deltaCP",true_deltaCP,0.01);
        fMinimizer->SetVariable(1,"deltamsq13",c_osc::deltamsq13,0.001*c_osc::deltamsq13);
        //fMinimizer->SetVariable(2,"deltamsq12",c_osc::deltamsq12,0.001*c_osc::deltamsq12);
        //fMinimizer->SetVariable(3,"theta23",c_osc::theta23,0.01*c_osc::theta23);
        //fMinimizer->SetVariable(4,"theta12",c_osc::theta12,0.01*c_osc::theta12);
        fMinimizer->SetVariable(2,"theta13",c_osc::theta13,0.01*c_osc::theta13);

        fMinimizer->SetVariableLimits(0,-3.1415/2,3.1415/2); 
        //fMinimizer->SetVariableLimits(1,0.5*deltamsq23,1.5*deltamsq23); 

        fMinimizer->SetFunction(min);
        fMinimizer->Minimize();

        //std::cout << "Data Osc Parameters: true_deltamsq23=" << true_deltamsq23 << " true_amp=" << true_amp << std::endl;
        //std::cout << "Fit  Osc Parameters:  fit_deltamsq23=" << fMinimizer->X()[0] << "  fit_amp=" << fMinimizer->X()[1] << std::endl;

        meas_deltaCP = fMinimizer->X()[0];

        if(makeloadsofplots){
          gSystem->Exec(("mkdir -p Plots/DeltaCPPlots/Point_"+std::to_string(i_dm)).c_str()); 
          h_data_reco->Draw("HIST");
          h_data_reco->SetLineColor(1);
          h_data_reco->SetLineWidth(2);
          h_data_reco->SetStats(0);
          h_model_reco->Draw("HIST same");
          h_model_reco->SetLineColor(2);
          h_model_reco->SetLineWidth(2);
          h_data_reco->SetTitle(("Input #delta_{CP}="+std::to_string(true_deltaCP)+" Measured #delta_{CP}="+std::to_string(meas_deltaCP)).c_str());
          l->AddEntry(h_data_reco,(Generators_v.at(0)+" FD").c_str(),"L");
          l->AddEntry(h_model_reco,(Generators_v.at(i_g)+" Model").c_str(),"L");
          l->Draw();
          c->Print(("Plots/DeltaCPPlots/Point_" + std::to_string(i_dm) +"/" + Generators_v.at(i_g) + "_" + estimators_str.at(i_e) +".png").c_str());
          c->Clear();
          l->Clear();
        } 

        h_fit_results.at(i_e).at(i_g)->SetBinContent(i_dm,meas_deltaCP-true_deltaCP);

        delete h_model_reco;
        delete h_model_true_reco;
      }

      delete h_data_reco;
      delete h_data_true_reco;

    }

  }


  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    THStack* hs = new THStack("hs",";Input #delta_{CP};Measured #delta_{CP} - Input #delta_{CP}");

    for(size_t i_g=0;i_g<Generators_v.size();i_g++){
      h_fit_results.at(i_e).at(i_g)->SetLineColor(i_g+1); 
      h_fit_results.at(i_e).at(i_g)->SetLineWidth(2); 
      hs->Add(h_fit_results.at(i_e).at(i_g));
      l->AddEntry(h_fit_results.at(i_e).at(i_g),Generators_v.at(i_g).c_str(),"L");
    }

    hs->Draw("HIST nostack"); 
    //hs->SetMinimum(-0.8);
    //hs->SetMaximum(0.8);
    l->Draw();
    c->Print(("Plots/DeltaCPPlots/FitResults_" + estimators_str.at(i_e) +".png").c_str());
    c->Clear();
    l->Clear();
  }


}
