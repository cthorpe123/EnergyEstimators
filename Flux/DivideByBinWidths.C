void Divide(TH1D* h){
  for(int i=1;i<h->GetNbinsX()+1;i++)
    h->SetBinContent(i,h->GetBinContent(i)/h->GetBinWidth(i));
}


void DivideByBinWidths(){

// Change title to have nus/GeV
std::string y_axis_title = "#nu/m^{2}/POT/GeV";

TFile* f_in = TFile::Open("histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEFD_fastmc.root");
TH1D* numu_flux = static_cast<TH1D*>(f_in->Get("numu_flux"));
TH1D* numu_fluxosc = static_cast<TH1D*>(f_in->Get("numu_fluxosc"));
numu_flux->SetDirectory(0);
numu_fluxosc->SetDirectory(0);
f_in->Close();
 
Divide(numu_flux);
Divide(numu_fluxosc);
numu_flux->GetYaxis()->SetTitle(y_axis_title.c_str());
numu_fluxosc->GetYaxis()->SetTitle(y_axis_title.c_str());

TFile* f_out = new TFile("DUNE_FD_Flux.root","RECREATE");
numu_flux->SetDirectory(f_out);
numu_fluxosc->SetDirectory(f_out);
numu_flux->Write();
numu_fluxosc->Write();
f_out->Close();


}
