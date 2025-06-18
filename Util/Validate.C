#include "../Funcs/Funcs.h"
#include "TLorentzVector.h"

void Validate(){

  std::vector<std::string> InputFiles_v = {"../rootfiles/NuWroEvents.root","../rootfiles/NEUTEvents.root","../rootfiles/GENIEEvents.root","../rootfiles/GiBUUEvents.root"};
  std::vector<std::string> Generators_v = {"NuWro","NEUT","GENIE","GiBUU"};

  std::vector<TH1D*> h_NeutrinoEnergy;
  std::vector<TH1D*> h_W;
  std::vector<TH1D*> h_LeptonMom;
  std::vector<TH1D*> h_LeptonAngle;

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    TFile* f = TFile::Open(InputFiles_v.at(i_f).c_str()) ;
    TTree* t = static_cast<TTree*>(f->Get("eventtree")) ;
    std::string generator = Generators_v.at(i_f);    

    h_NeutrinoEnergy.push_back(new TH1D((generator+"_NeutrinoEnergy").c_str(),"",100,0.0,10.0));
    h_W.push_back(new TH1D((generator+"_W").c_str(),"",100,0.0,5.0));
    h_LeptonMom.push_back(new TH1D((generator+"_LeptonMom").c_str(),"",100,0.0,8.0));
    h_LeptonAngle.push_back(new TH1D((generator+"_LeptonAngle").c_str(),"",100,-1,1.0));

    Double_t scale;
    Double_t weight;
    Double_t nu_e;
    Int_t ccnc;
    Int_t nu_pdg;  
    Int_t lepton_pdg;
    TLorentzVector* lepton_p4=0;
    std::vector<int>* pdg=0;
    std::vector<TLorentzVector>* p4=0;

    t->SetBranchAddress("scale",&scale);
    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){
      t->GetEntry(ievent);

      //if(ievent > 1000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t->GetEntries() << std::endl;

      if(ccnc != 1 || nu_pdg != 14) continue;

      TLorentzVector p4_tot(0,0,0,0);
      int prot = 0;
      for(size_t i_p=0;i_p<pdg->size();i_p++){
        if(pdg->at(i_p) == 2212 && p4->at(i_p).Vect().Mag() > 0.3) p4_tot += p4->at(i_p); 
        if(abs(pdg->at(i_p)) == 211 && p4->at(i_p).Vect().Mag() > 0.1) p4_tot += p4->at(i_p); 
        if(pdg->at(i_p) == 111) p4_tot += p4->at(i_p); 
        if(pdg->at(i_p) == 2212 && p4->at(i_p).Vect().Mag() > 0.3) prot++;
      } 

      if(generator != "GiBUU") weight = 1;

      weight *= scale*1e38*40;

      h_NeutrinoEnergy.back()->Fill(nu_e,weight);
      h_LeptonMom.back()->Fill(lepton_p4->Vect().Mag(),weight); 
      h_LeptonAngle.back()->Fill(lepton_p4->Vect().CosTheta(),weight); 
      h_W.back()->Fill(p4_tot.M(),weight);
    }

  }

  gSystem->Exec("mkdir -p ValidationPlots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  THStack* hs_NeutrinoEnergy = new THStack("hs",";Neutrino Energy (GeV);Shape");
  THStack* hs_W = new THStack("hs",";W (GeV);d#sigma/dW (10^{-38} cm^{2}/GeV/Nucleus)");
  THStack* hs_LeptonMom = new THStack("hs",";Lepton Momentum (GeV);d#sigma/dP (10^{-38} cm^{2}/GeV/Nucleus)");
  THStack* hs_LeptonAngle = new THStack("hs",";Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-38} cm^{2}/Nucleus)");

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    l->AddEntry(h_NeutrinoEnergy.at(i_f),Generators_v.at(i_f).c_str(),"L");

    h_NeutrinoEnergy.at(i_f)->SetLineWidth(2);
    h_NeutrinoEnergy.at(i_f)->SetLineColor(i_f+2);
    hs_NeutrinoEnergy->Add(h_NeutrinoEnergy.at(i_f));

    Reweight(h_W.at(i_f));
    h_W.at(i_f)->SetLineWidth(2);
    h_W.at(i_f)->SetLineColor(i_f+2);
    hs_W->Add(h_W.at(i_f));

    Reweight(h_LeptonMom.at(i_f));
    h_LeptonMom.at(i_f)->SetLineWidth(2);
    h_LeptonMom.at(i_f)->SetLineColor(i_f+2);
    hs_LeptonMom->Add(h_LeptonMom.at(i_f));

    Reweight(h_LeptonAngle.at(i_f));
    h_LeptonAngle.at(i_f)->SetLineWidth(2);
    h_LeptonAngle.at(i_f)->SetLineColor(i_f+2);
    hs_LeptonAngle->Add(h_LeptonAngle.at(i_f));

  }

  hs_NeutrinoEnergy->Draw("HIST nostack");
  l->Draw();
  c->Print("ValidationPlots/NeutrinoEnergy.png");  
  c->Clear();

  hs_W->Draw("HIST nostack");
  l->Draw();
  c->Print("ValidationPlots/W.png");  
  c->Clear();

  hs_LeptonMom->Draw("HIST nostack");
  l->Draw();
  c->Print("ValidationPlots/LeptonMom.png");  
  c->Clear();

  hs_LeptonAngle->Draw("HIST nostack");
  l->Draw();
  c->Print("ValidationPlots/LeptonAngle.png");  
  c->Clear();

}
