#ifndef _PlotSetup_h_
#define _PlotSetup_h_

TCanvas* c = new TCanvas("c","c",800,600);
TPad *p_plot = new TPad("p_plot","p_plot",0,0,1,0.85);
TPad *p_legend = new TPad("p_legend","p_legend",0,0.85,1,1);
TLegend* l = new TLegend(0.1,0.0,0.9,1.0);

void PlotSetup(){

  p_legend->SetBottomMargin(0);
  p_legend->SetTopMargin(0.1);
  p_plot->SetTopMargin(0.01);
  p_plot->SetBottomMargin(0.13);
  p_plot->SetLeftMargin(0.1);

  l->SetBorderSize(0);
  l->SetNColumns(3);

  p_legend->Draw();
  p_legend->cd();
  l->Draw();
  c->cd();
  p_plot->Draw();
  p_plot->cd();

}

void SetAxisFonts(THStack* hs){

  hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->GetXaxis()->SetTitleOffset(1.05);
  hs->GetYaxis()->SetTitleOffset(1.03);
  hs->GetXaxis()->SetLabelSize(0.045); 
  hs->GetYaxis()->SetLabelSize(0.045); 

}

void SetAxisFontsMG(TMultiGraph* mg){

  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleOffset(1.05);
  mg->GetYaxis()->SetTitleOffset(1.03);
  mg->GetXaxis()->SetLabelSize(0.045); 
  mg->GetYaxis()->SetLabelSize(0.045); 

}

void SetAxisFontsH(TH1D* h){

  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleOffset(1.03);
  h->GetXaxis()->SetLabelSize(0.045); 
  h->GetYaxis()->SetLabelSize(0.045); 

}

#endif
