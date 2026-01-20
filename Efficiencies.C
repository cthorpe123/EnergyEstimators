#include "Funcs/EnergyEstimatorFuncs.h"
#include "Funcs/Funcs.h"
#include "TLorentzVector.h"
#include "Funcs/PlotSetup.h"

double osc_strength = 0.0;
std::vector<int> styles = {1,2,7,9};

void Efficiencies(){

  PlotSetup();

  std::vector<std::string> InputFiles_v = {"GENIEEvents.root","NuWroEvents.root","NEUTEvents.root","GiBUUEvents.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<double> true_binning_v;
  true_binning_v.push_back(0.5);
  for(int i=0;i<10;i++) true_binning_v.push_back(true_binning_v.back()+0.25);
  for(int i=0;i<4;i++) true_binning_v.push_back(true_binning_v.back()+0.5);
  true_binning_v.push_back(true_binning_v.back()+1.0);
  int true_nbins = true_binning_v.size()-1;
  double* true_binning_a = &true_binning_v[0];

  std::vector<TH1D*> h_TrueEnergy;
  std::vector<TH1D*> h_TrueEnergy_Np;
  std::vector<TH1D*> h_TrueEnergy_Np0pi;
  std::vector<TH1D*> h_TrueEnergy_1p0pi;

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    h_TrueEnergy.push_back(new TH1D(("h_TrueEnergy_"+generator).c_str(),"",true_nbins,true_binning_a));
    h_TrueEnergy_Np.push_back(new TH1D(("h_TrueEnergy_Np"+generator).c_str(),"",true_nbins,true_binning_a));
    h_TrueEnergy_Np0pi.push_back(new TH1D(("h_TrueEnergy_Np0pi"+generator).c_str(),"",true_nbins,true_binning_a));
    h_TrueEnergy_1p0pi.push_back(new TH1D(("h_TrueEnergy_1P0pi"+generator).c_str(),"",true_nbins,true_binning_a));

    TFile* f = TFile::Open(("/gluster/data/dune/cthorpe/DIS/"+InputFiles_v.at(i_f)).c_str());
    TTree* t = static_cast<TTree*>(f->Get("eventtree")) ;

    Double_t weight;
    Double_t nu_e;
    Int_t ccnc;
    Int_t nu_pdg;  
    Int_t lepton_pdg;
    TLorentzVector* lepton_p4=0;
    std::vector<int>* pdg=0;
    std::vector<TLorentzVector>* p4=0;

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);

    double events = 0.0;
    double events_Np = 0.0;
    double events_Np0pi = 0.0;
    double events_1p0pi = 0.0;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != 14 || ccnc != 1) continue;

      int nprot = GetNProt(pdg,p4);
      std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
      std::vector<TVector3> neutron_mom = GetNeutronMom(pdg,p4);
      std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
      std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);

      h_TrueEnergy.back()->Fill(nu_e,weight);

      if(nprot < 1) continue;
      h_TrueEnergy_Np.back()->Fill(nu_e,weight);

      if(pion_mom.size() || pizero_mom.size()) continue;
      h_TrueEnergy_Np0pi.back()->Fill(nu_e,weight);

      if(nprot != 1) continue;
      h_TrueEnergy_1p0pi.back()->Fill(nu_e,weight);

    }

  }

  gSystem->Exec("mkdir -p Plots/EfficiencyPlots/");

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
    THStack* hs = new THStack("hs",";True Neutrino Energy (GeV);Fraction");

    h_TrueEnergy_Np.at(i_f)->Divide(h_TrueEnergy.at(i_f));
    h_TrueEnergy_Np.at(i_f)->SetLineWidth(2);
    h_TrueEnergy_Np.at(i_f)->SetLineColor(2);
    hs->Add(h_TrueEnergy_Np.at(i_f));
    l->AddEntry(h_TrueEnergy_Np.at(i_f),"Np","L");

    h_TrueEnergy_Np0pi.at(i_f)->Divide(h_TrueEnergy.at(i_f));
    h_TrueEnergy_Np0pi.at(i_f)->SetLineWidth(2);
    h_TrueEnergy_Np0pi.at(i_f)->SetLineColor(3);
    hs->Add(h_TrueEnergy_Np0pi.at(i_f));
    l->AddEntry(h_TrueEnergy_Np0pi.at(i_f),"Np0#pi","L");

    h_TrueEnergy_1p0pi.at(i_f)->Divide(h_TrueEnergy.at(i_f));
    h_TrueEnergy_1p0pi.at(i_f)->SetLineWidth(2);
    h_TrueEnergy_1p0pi.at(i_f)->SetLineColor(4);
    hs->Add(h_TrueEnergy_1p0pi.at(i_f));
    l->AddEntry(h_TrueEnergy_1p0pi.at(i_f),"1p0#pi","L");

    p_plot->cd();
    hs->Draw("HIST nostack");
    SetAxisFonts(hs);

    c->Print(("Plots/EfficiencyPlots/"+Generators_v.at(i_f)+".pdf").c_str());
    l->Clear();
    p_plot->Clear();

    delete hs;
  } 

  l->SetNColumns(InputFiles_v.size());

  THStack* hs_all = new THStack("hs_all",";True Neutrino Energy (GeV);Fraction");
  
  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
    h_TrueEnergy_Np.at(i_f)->SetLineColor(i_f+1);
    h_TrueEnergy_Np0pi.at(i_f)->SetLineColor(i_f+1); 
    h_TrueEnergy_1p0pi.at(i_f)->SetLineColor(i_f+1); 
    l->AddEntry(h_TrueEnergy_Np.at(i_f),(Generators_v.at(i_f)+" Np").c_str(),"L"); 
  }
   
  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
    h_TrueEnergy_Np.at(i_f)->SetLineStyle(1); 
    h_TrueEnergy_Np0pi.at(i_f)->SetLineStyle(2); 
    h_TrueEnergy_1p0pi.at(i_f)->SetLineStyle(3); 
    l->AddEntry(h_TrueEnergy_Np0pi.at(i_f),"Np0#pi","L"); 
  }

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
    hs_all->Add(h_TrueEnergy_Np.at(i_f));
    hs_all->Add(h_TrueEnergy_Np0pi.at(i_f));
    hs_all->Add(h_TrueEnergy_1p0pi.at(i_f));
    l->AddEntry(h_TrueEnergy_1p0pi.at(i_f),"1p0#pi","L"); 
  }
 
  p_plot->cd(); 
  hs_all->Draw("nostack HIST");
  SetAxisFonts(hs_all);

  c->Print("Plots/EfficiencyPlots/Efficiencies.pdf");
  c->Close();

}
