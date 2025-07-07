#include "Funcs/Funcs.h"
#include "Funcs/OscFitter.h"

OscModel data_osc_model;
OscModel fit_osc_model;

// Fit with grid search
double gs_deltam2_low = 0.5*c_osc::deltamsq23;
double gs_deltam2_high = 1.5*c_osc::deltamsq23;
void FitGridSearch(const TH1D* h_data_reco,TH1D* h_flux,const TH2D* h_model_true_reco,TH1D* h_model_reco,double& best_deltam2){

  int nbins = h_model_true_reco->GetNbinsX();
  double deltam2 = gs_deltam2_low;
  double best_chi2 = 1e10;

  while(deltam2 <= gs_deltam2_high){ 

    fit_osc_model.SetDeltaMSq23(deltam2);

    // Calculate the oscillated prediction
    for(int j=1;j<nbins+1;j++){
      double events = 0.0;
      for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
        double flux = h_flux->GetBinContent(h_flux->FindBin(h_model_true_reco->GetXaxis()->GetBinCenter(i)));
        double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
        events += fit_osc_model.NuMuSurvProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
      }
      h_model_reco->SetBinContent(j,events);
    }

    double chi2 = 0.0;
    for(int j=1;j<nbins+1;j++){
      if(h_data_reco->GetBinContent(j) > 0)
        chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/sqrt(h_data_reco->GetBinContent(j)),2);         
    }

   if(chi2 < best_chi2){
       best_chi2 = chi2;
       best_deltam2 = deltam2;
   }

   deltam2 += (gs_deltam2_high - gs_deltam2_low)/500;

  }

  fit_osc_model.SetDeltaMSq23(best_deltam2);

  for(int j=1;j<nbins+1;j++){
    double events = 0.0;
    for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
      double flux = h_flux->GetBinContent(h_flux->FindBin(h_model_true_reco->GetXaxis()->GetBinCenter(i)));
      double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
      events += fit_osc_model.NuMuSurvProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
    }
    h_model_reco->SetBinContent(j,events);
  }
  
}

// Use the first generator in the list as a model of the mapping
// from true neutrino energy into reco neutrino energy
// then use the responses from the other generators to make the reco
// distribution and fit deltam2 to each

bool makeloadsofplots = true;

void DeltaM2Plot(){

  gSystem->Exec("mkdir -p Plots/DeltaM2Plots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT"};

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

  double true_theta13 = c_osc::theta13;
  double true_theta23 = c_osc::theta23;

  std::vector<std::vector<TH1D*>> h_fit_results(estimators_str.size(),std::vector<TH1D*>());
  for(size_t i_e=0;i_e<estimators_str.size();i_e++)
    for(size_t i_g=0;i_g<Generators_v.size();i_g++)
      h_fit_results.at(i_e).push_back(new TH1D(("h_fit_results_"+estimators_str.at(i_e)+"_"+Generators_v.at(i_g)).c_str(),";Input #Delta m^{2}_{23};Measured #Delta m^{2}_{23}/Input #Delta m^{2}_{23}",100,0.8*c_osc::deltamsq23,1.2*c_osc::deltamsq23));
 

  double min_fit_ratio = 1;
  double max_fit_ratio = 1;
 
  for(int i_dm=1;i_dm<h_fit_results.back().back()->GetNbinsX()+1;i_dm++){

    double true_deltamsq23 = h_fit_results.back().back()->GetBinCenter(i_dm);
    data_osc_model.SetDeltaMSq23(true_deltamsq23);

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
          events += data_osc_model.NuMuSurvProb(E)*h_data_true_reco->GetBinContent(i,j)*flux;
        }
        h_data_reco->SetBinContent(j,events);
      }

      for(size_t i_g=0;i_g<Generators_v.size();i_g++){

        std::cout << "Fitting " << estimators_str.at(i_e) << " " << Generators_v.at(i_g) << std::endl;

        TH2D* h_model_true_reco = static_cast<TH2D*>(h_energy_true_reco.at(i_e).at(i_g)->Clone("h_model_true_reco"));
        TH1D* h_model_reco = new TH1D("h_model_reco","",nbins,low,high);

        double meas_deltamsq23 = 0; 
        //FitGridSearch(h_data_reco,h_flux,h_model_true_reco,h_model_reco,meas_deltamsq23);

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

        if(makeloadsofplots){
          gSystem->Exec(("mkdir -p Plots/DeltaM2Plots/Point_"+std::to_string(i_dm)).c_str()); 
          h_data_reco->Draw("HIST");
          h_data_reco->SetLineColor(1);
          h_data_reco->SetLineWidth(2);
          h_data_reco->SetStats(0);
          h_model_reco->Draw("HIST same");
          h_model_reco->SetLineColor(2);
          h_model_reco->SetLineWidth(2);
          //h_data_reco->SetTitle(("Input #Delta m^{2}="+std::to_string(true_deltamsq23)+" Measured #Delta m^{2}="+std::to_string(fMinimizer->X()[0])).c_str());
          h_data_reco->SetTitle(("Input #Delta m^{2}="+std::to_string(true_deltamsq23)+" Measured #Delta m^{2}="+std::to_string(meas_deltamsq23)).c_str());
          l->AddEntry(h_data_reco,(Generators_v.at(0)+" FD").c_str(),"L");
          l->AddEntry(h_model_reco,(Generators_v.at(i_g)+" Model").c_str(),"L");
          l->Draw();
          c->Print(("Plots/DeltaM2Plots/Point_" + std::to_string(i_dm) +"/" + Generators_v.at(i_g) + "_" + estimators_str.at(i_e) +".png").c_str());
          c->Clear();
          l->Clear();
        } 

        h_fit_results.at(i_e).at(i_g)->SetBinContent(i_dm,meas_deltamsq23/true_deltamsq23);

        min_fit_ratio = std::min(min_fit_ratio,meas_deltamsq23/true_deltamsq23);
        max_fit_ratio = std::max(max_fit_ratio,meas_deltamsq23/true_deltamsq23);

        delete h_model_reco;
        delete h_model_true_reco;
      }

      delete h_data_reco;
      delete h_data_true_reco;

    }

  }


  for(size_t i_e=0;i_e<estimators_str.size();i_e++){

    THStack* hs = new THStack("hs",";Input #Delta m^{2}_{23} (eV^{2});Measured #Delta m^{2}_{23}/Input #Delta m^{2}_{23}");     

    for(size_t i_g=0;i_g<Generators_v.size();i_g++){
      h_fit_results.at(i_e).at(i_g)->SetLineColor(i_g+1); 
      h_fit_results.at(i_e).at(i_g)->SetLineWidth(2); 
      hs->Add(h_fit_results.at(i_e).at(i_g));
      l->AddEntry(h_fit_results.at(i_e).at(i_g),Generators_v.at(i_g).c_str(),"L");
    }

    hs->Draw("HIST nostack"); 
    hs->SetMinimum(min_fit_ratio-0.1*(max_fit_ratio-min_fit_ratio));
    hs->SetMaximum(max_fit_ratio+0.1*(max_fit_ratio-min_fit_ratio));
    l->Draw();
    c->Print(("Plots/DeltaM2Plots/FitResults_" + estimators_str.at(i_e) +".png").c_str());
    c->Clear();
    l->Clear();
}


}
