#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"

OscModel data_osc_model;
OscModel fit_osc_model;

bool makeloadsofplots = true;

void DeltaM2Plot(){

  // Load the histograms
  TFile* f = TFile::Open("ResponseMatricesNuMu.root");

  std::vector<std::string> Generators_v = {"GiBUU","NEUT","GENIE","NuWro"};

  gSystem->Exec("mkdir -p Plots/");
  gSystem->Exec("mkdir -p Plots/FitPlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("../Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  h_flux->SetDirectory(0);
  f_flux->Close();
  h_flux->Scale(1.0/h_flux->Integral());

  for(size_t i_f=1;i_f<Generators_v.size();i_f++){
    std::string gen = Generators_v.at(i_f);    

    std::vector<TH1D*> h_fit_results(estimators_str.size());
    for(size_t i_e=0;i_e<estimators_str.size();i_e++)
      h_fit_results.at(i_e) = new TH1D(("h_fit_results_"+estimators_str.at(i_e)+"_"+gen).c_str(),";Input #Delta m^{2}_{23};Measured #Delta m^{2}_{23}/Input #Delta m^{2}_{23}",50,0.8*c_osc::deltamsq23,1.2*c_osc::deltamsq23);

    double true_theta13 = c_osc::theta13;
    double true_theta23 = c_osc::theta23;

    double min_fit_ratio = 1;
    double max_fit_ratio = 1;

    for(int i_dm=1;i_dm<h_fit_results.back()->GetNbinsX()+1;i_dm++){

      double true_deltamsq23 = h_fit_results.back()->GetBinCenter(i_dm);
      data_osc_model.SetDeltaMSq23(true_deltamsq23);

      for(size_t i_e=0;i_e<estimators_str.size();i_e++){

        TH2D* h = static_cast<TH2D*>(f->Get((Generators_v.at(0)+"_TrueEnergy_RecoEnergy_"+estimators_str.at(i_e)).c_str())); 
        TH2D* h_nofsi = static_cast<TH2D*>(f->Get((gen+"_TrueEnergy_RecoEnergy_"+estimators_str.at(i_e)).c_str())); 

        if(i_dm == 1){
          Normalise(h);
          Normalise(h_nofsi);
        }

        //if(i_e == kMuonKin) continue;

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
            if(E < 0.8 || E > 5) continue;
            events += data_osc_model.NuMuSurvProb(E)*h_data_true_reco->GetBinContent(i,j)*flux;
          }
          h_data_reco->SetBinContent(j,events);
        }

        const TH2D* h_model_true_reco = h_nofsi;

        std::cout << "Fitting " << estimators_str.at(i_e) << " " << gen << std::endl;

        double meas_deltamsq23 = 0; 
        TH1D* h_model_reco = new TH1D("h_model_reco","",nbins,low,high);
        TH1D* h_model_reco_prefit = new TH1D("h_model_reco_prefit","",nbins,low,high);

        // Draw the model prediction with the input deltam2/amp 
        for(int j=1;j<nbins+1;j++){
          double events = 0.0;
          for(int i=1;i<h_model_true_reco->GetNbinsX()+1;i++){
            double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
            double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
            if(E < 0.8 || E > 5) continue;
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
            double flux = h_flux->GetBinContent(h_flux->FindBin(h_data_true_reco->GetXaxis()->GetBinCenter(i)));
            double E = h_model_true_reco->GetXaxis()->GetBinCenter(i);
            if(E < 0.8 || E > 5) continue;
            events += fit_osc_model.NuMuSurvProb(E)*h_model_true_reco->GetBinContent(i,j)*flux;
            }
            h_model_reco->SetBinContent(j,events);
            }

            h_model_reco->Scale(coeff[2]);

            double chi2 = 0.0;
            for(int j=1;j<nbins+1;j++){
            if(h_data_reco->GetBinContent(j) > 0)
            chi2 += pow((h_model_reco->GetBinContent(j) - h_data_reco->GetBinContent(j))/*/sqrt(h_data_reco->GetBinContent(j))*/,2);         
            }

            return chi2/nbins*100000;

        }, 2);

        std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
          ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

        fMinimizer->SetVariable(0,"deltamsq23",true_deltamsq23,0.01*c_osc::deltamsq23);
        fMinimizer->SetVariable(1,"amp",data_osc_model.numu_dis_amp,0.001*data_osc_model.numu_dis_amp);
        fMinimizer->SetVariable(2,"scale",1.0,0.01);

        fMinimizer->SetFunction(min);
        fMinimizer->Minimize();

        meas_deltamsq23 = fMinimizer->X()[0];

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
          h_model_reco_prefit->SetLineColor(3);
          h_model_reco_prefit->SetLineWidth(2);
          h_data_reco->SetTitle(("Input #Delta m^{2}="+std::to_string(true_deltamsq23)+" Measured #Delta m^{2}="+std::to_string(meas_deltamsq23)).c_str());
          l->AddEntry(h_data_reco,(gen+" FD").c_str(),"L");
          l->AddEntry(h_model_reco_prefit,(gen+" Model, No Fit").c_str(),"L");
          l->AddEntry(h_model_reco,(gen+" Model, Fitted").c_str(),"L");
          l->Draw();
          //c->Print(("Plots/DeltaM2Plots/" + gen + "_" + estimators_str.at(i_e) +".pdf").c_str());
          c->Print(("Plots/FitPlots/"+gen+"/"+estimators_str.at(i_e)+"/Point_" + std::to_string(i_dm) + "_" + gen + "_" + estimators_str.at(i_e) +".pdf").c_str());
          c->Clear();
          l->Clear();


        } 

        h_fit_results.at(i_e)->SetBinContent(i_dm,meas_deltamsq23/true_deltamsq23);
        min_fit_ratio = std::min(min_fit_ratio,meas_deltamsq23/true_deltamsq23);
        max_fit_ratio = std::max(max_fit_ratio,meas_deltamsq23/true_deltamsq23);

        delete h_data_reco;
        delete h_model_reco;
        delete h_model_reco_prefit;

      }
    }

    THStack* hs = new THStack("hs",";Input #Delta m^{2}_{23} (eV^{2});Measured #Delta m^{2}_{23}/Input #Delta m^{2}_{23}");     

    for(size_t i_e=0;i_e<estimators_str.size();i_e++){
      //if(i_e == kMuonKin) continue;
      h_fit_results.at(i_e)->SetLineColor(colors.at(i_e)); 
      h_fit_results.at(i_e)->SetLineWidth(2); 
      hs->Add(h_fit_results.at(i_e));
      l->AddEntry(h_fit_results.at(i_e),estimators_str.at(i_e).c_str(),"L");
    }

    TF1* func = new TF1("func","1");
    func->SetLineColor(1);
    func->SetLineWidth(2);
    func->SetLineStyle(9);

    hs->Draw("HIST nostack"); 
    hs->SetMinimum(min_fit_ratio-0.1*(max_fit_ratio-min_fit_ratio));
    hs->SetMaximum(max_fit_ratio+0.1*(max_fit_ratio-min_fit_ratio));
    func->Draw("L same");
    l->Draw();
    c->Print(("Plots/FitResults_"+gen+".pdf").c_str());
    c->Clear();
    l->Clear();

  }

}
