#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"

void MakeResponseMatrices(){

  bool nue_mode = false;

  std::vector<std::string> InputFiles_v;
  if(!nue_mode) InputFiles_v = {"GENIEEventsFiltered.root","NuWroEventsFiltered.root","NEUTEventsFiltered.root","GiBUUEventsFiltered.root"};
  else InputFiles_v = {"GENIE_NueEventsFiltered.root","NuWro_NueEventsFiltered.root","NEUT_NueEventsFiltered.root","GiBUU_NueEventsFiltered.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  // First make the response matrices
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy(kMAX);

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

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
    int nprot;
    double W;
    std::vector<double>* est_nu_e=0;

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);
    t->SetBranchAddress("W",&W);
    t->SetBranchAddress("nprot",&nprot);
    t->SetBranchAddress("est_nu_e",&est_nu_e);

    for(size_t i_e=0;i_e<kMAX;i_e++){
      std::string est = estimators_str.at(i_e);
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);",200,0.0,8.0,200,0.0,8.0));
      h_TrueEnergy_RecoEnergy.at(i_e).back()->SetDirectory(0);
    }

    int target_nu_pdg = 14;
    if(nue_mode) target_nu_pdg = 12;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != target_nu_pdg || ccnc != 1) continue;

      for(int i_e=0;i_e<kMAX;i_e++){
        h_TrueEnergy_RecoEnergy.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),weight); 
      }

    }

    f->Close();

  }

  std::string filename = nue_mode ? "ResponseMatricesNue.root" : "ResponseMatricesNuMu.root";
  TFile* f_out = new TFile(filename.c_str(),"RECREATE");

  for(size_t i_e=0;i_e<kMAX;i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      h_TrueEnergy_RecoEnergy.at(i_e).at(i_f)->Write();
    }
  }

  f_out->Close(); 

}

