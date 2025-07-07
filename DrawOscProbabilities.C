#include "Funcs/Funcs.h"
#include "Funcs/OscFitter.h"

void DrawOscProbabilities(){

  OscModel osc_model;

  // Load the numu flux histogram
  TFile* f_flux = TFile::Open("Flux/DUNE_FD_Flux.root");
  TH1D* h_flux = static_cast<TH1D*>(f_flux->Get("numu_flux"));
  TH1D* h_flux_osc = static_cast<TH1D*>(f_flux->Get("numu_fluxosc"));
  h_flux->SetDirectory(0);
  h_flux_osc->SetDirectory(0);
  f_flux->Close();

  h_flux_osc->Divide(h_flux);

  std::vector<double> E,nue_app,numu_surv,nue_app_ppi2,nue_app_mpi2;

  double e = 0.2;
  for(int i=0;i<1000;i++){
    E.push_back(e);

    osc_model.SetDeltaCP(0);
    nue_app.push_back(osc_model.NueAppProb(e));

    osc_model.SetDeltaCP(3.141/2);
    nue_app_ppi2.push_back(osc_model.NueAppProb(e));

    osc_model.SetDeltaCP(-3.141/2);
    nue_app_mpi2.push_back(osc_model.NueAppProb(e));

    osc_model.SetDeltaCP(0);
    numu_surv.push_back(osc_model.NuMuSurvProb(e));

    e += 0.01;
  }

  TGraph* g_nue_app = new TGraph(E.size(),&(E[0]),&(nue_app[0]));
  TGraph* g_nue_app_ppi2 = new TGraph(E.size(),&(E[0]),&(nue_app_ppi2[0]));
  TGraph* g_nue_app_mpi2 = new TGraph(E.size(),&(E[0]),&(nue_app_mpi2[0]));
  TGraph* g_numu_surv = new TGraph(E.size(),&(E[0]),&(numu_surv[0]));

  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.65,0.65,0.98,0.98);
  c->SetLogx();

  g_nue_app->Draw("AC");
  g_nue_app->SetTitle(";E_{#nu} (GeV);Probability");
  l->AddEntry(g_nue_app,"#nu_{e} Appearance Prob","L");
  g_nue_app->SetLineWidth(2);
  g_nue_app->SetLineColor(2);

  g_nue_app_ppi2->Draw("C same");
  l->AddEntry(g_nue_app_ppi2,"#nu_{e} Appearance Prob, #delta_{CP}=#pi/2","L");
  g_nue_app_ppi2->SetLineWidth(2);
  g_nue_app_ppi2->SetLineColor(3);

  g_nue_app_mpi2->Draw("C same");
  l->AddEntry(g_nue_app_mpi2,"#nu_{e} Appearance Prob, #delta_{CP}=-#pi/2","L");
  g_nue_app_mpi2->SetLineWidth(2);
  g_nue_app_mpi2->SetLineColor(4);

  g_numu_surv->Draw("C same");  
  l->AddEntry(g_numu_surv,"#nu_{#mu} Survival Prob","L");
  g_numu_surv->SetLineWidth(2);
  g_numu_surv->SetLineColor(1);

  h_flux_osc->Draw("HIST same");
  h_flux_osc->SetLineColor(2);  
  l->AddEntry(h_flux_osc,"DR #nu_{#mu} Survival Prob","L");

  l->Draw();

  g_nue_app->GetHistogram()->SetMaximum(1.0);
  g_nue_app->GetHistogram()->SetMinimum(0.0);
  c->Update();
  c->Print("Plots/P.png");

}
