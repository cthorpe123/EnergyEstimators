#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/PlotSetup.h"
#include "../Funcs/OscFitter.h"

void BiasWithCuts(){

  PlotSetup(); 

  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<std::string> vars = {"W","Angle","MissingE","Neutrons"};
  std::vector<std::string> x_axis_titles = {"W_{vis} (GeV)","#theta (rad)","Missing Hadronic Energy (GeV)","N"};
  std::vector<std::string> y_axis_titles = {"B' - B","B' - B","B' - B","B' - B",};

  TFile* f = TFile::Open("ResponseMatricesNuMu.root");

  for(int i_v=0;i_v<vars.size();i_v++){

    std::string var = vars.at(i_v);

    gSystem->Exec(("mkdir -p Plots/" + var).c_str());

    for(size_t i_f=0;i_f<Generators_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);    

      std::string axis_title; 
      std::vector<TH1D*> h_Bias(kMAX);
      for(size_t i_e=0;i_e<estimators_str.size();i_e++){

        if(var == "W" && i_e == kSFMethod) continue;

        std::string est = estimators_str.at(i_e);

        TH3D* h = static_cast<TH3D*>(f->Get((gen+"_TrueEnergy_RecoEnergy_"+var+"_"+est).c_str())); 
        h_Bias.at(i_e) = static_cast<TH1D*>(h->ProjectionZ()->Clone(("h_Bias_"+est+"_"+gen).c_str()));

        for(int i_z=1;i_z<h->GetNbinsZ()+1;i_z++){
          h->GetZaxis()->SetRange(1,i_z); 
          TH2D* h_true_reco =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_true_reco")); 

          // Calculate the bias in this slice as single number
          double bias = 0.;
          double events = 0.;
          for(int i_x=1;i_x<h_true_reco->GetNbinsX()+1;i_x++){
            double e_true = h_true_reco->GetXaxis()->GetBinCenter(i_x);
            for(int i_y=1;i_y<h_true_reco->GetNbinsY()+1;i_y++){
              bias += (h_true_reco->GetYaxis()->GetBinCenter(i_y) - e_true)/e_true*h_true_reco->GetBinContent(i_x,i_y);
              events += h_true_reco->GetBinContent(i_x,i_y);
            } 
          }

          h_Bias.at(i_e)->SetBinContent(i_z,bias/events);
          delete h_true_reco;
        }

        for(size_t i_z=0;i_z<h_Bias.at(i_e)->GetNbinsX()+1;i_z++) 
          h_Bias.at(i_e)->SetBinContent(i_z,h_Bias.at(i_e)->GetBinContent(i_z) - h_Bias.at(i_e)->GetBinContent(h_Bias.at(i_e)->GetNbinsX()));

        axis_title = h->GetZaxis()->GetTitle();

      }


      TF1* f_line = new TF1("f_line","0",-1000,1000);
      f_line->SetLineColor(1);  
      f_line->SetLineWidth(2);
      f_line->SetLineStyle(9);

      std::string title = ";"+x_axis_titles.at(i_v)+";"+y_axis_titles.at(i_v);
      THStack* hs = new THStack("hs",title.c_str());
      for(int i_e=0;i_e<h_Bias.size();i_e++){
        if(var == "W" && i_e == kSFMethod) continue;
        std::string est = estimators_str.at(i_e);
        h_Bias.at(i_e)->SetLineWidth(2);
        h_Bias.at(i_e)->SetLineColor(colors.at(i_e));
        l->AddEntry(h_Bias.at(i_e),estimators_leg.at(i_e).c_str(),"L");
        hs->Add(h_Bias.at(i_e));
      }

      hs->Draw("nostack HIST");
      f_line->Draw("L same");
      SetAxisFonts(hs);
      hs->GetYaxis()->SetTitleOffset(0.95);
      c->Print(("Plots/"+var+"/BiasWithCuts_"+gen+".pdf").c_str());
      l->Clear();
      p_plot->Clear();

      delete f_line;
    }

  }

}

