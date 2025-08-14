#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"

void MakeResponseMatrices(){

  gSystem->Exec("mkdir -p Plots/DeltaM2Plots/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.95,0.95);

  std::vector<std::string> InputFiles_v = {"GENIEEvents2.root","NuWroEvents2.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro"};

  // First make the response matrices
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy(kMAX);
  std::vector<std::vector<TH2D*>> h_TrueEnergy_RecoEnergy_NoFSI(kMAX);

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    
    std::string inputfile = InputFiles_v.at(i_f);

    TFile* f = TFile::Open(("/gluster/data/dune/cthorpe/DIS/"+inputfile).c_str());
    TTree* t = static_cast<TTree*>(f->Get("eventtree")) ;

    Double_t weight;
    Double_t nu_e;
    Int_t ccnc;
    Int_t nu_pdg;  
    Int_t lepton_pdg;
    TLorentzVector* lepton_p4=0;
    std::vector<int>* pdg=0;
    std::vector<TLorentzVector>* p4=0;
    std::vector<int>* pdg_nofsi=0;
    std::vector<TLorentzVector>* p4_nofsi=0;

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nu_e",&nu_e);
    t->SetBranchAddress("nu_pdg",&nu_pdg);
    t->SetBranchAddress("ccnc",&ccnc);
    t->SetBranchAddress("lepton_pdg",&lepton_pdg); 
    t->SetBranchAddress("lepton_p4",&lepton_p4);
    t->SetBranchAddress("pdg",&pdg);
    t->SetBranchAddress("p4",&p4);
    t->SetBranchAddress("pdg_nofsi",&pdg_nofsi);
    t->SetBranchAddress("p4_nofsi",&p4_nofsi);

    for(size_t i_e=0;i_e<kMAX;i_e++){
      std::string est = estimators_str.at(i_e);
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_"+est).c_str(),"",200,0.2,8.0,200,0.2,8.0));
      h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).push_back(new TH2D((generator+"_TrueEnergy_RecoEnergy_NoFSI_"+est).c_str(),"",200,0.2,8.0,200,0.2,8.0));
      h_TrueEnergy_RecoEnergy.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).back()->SetDirectory(0);
    }

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 200000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != 14 || ccnc != 1) continue;

      std::vector<double> energies = GetEnergyEst(lepton_p4,pdg,p4);
      std::vector<double> energies_nofsi = GetEnergyEst(lepton_p4,pdg_nofsi,p4_nofsi);

      for(int i_e=0;i_e<kMAX;i_e++){
        double nu_e_reco = energies.at(i_e);
        double nu_e_reco_nofsi = energies_nofsi.at(i_e);
        if(nu_e_reco > 0){
          h_TrueEnergy_RecoEnergy.at(i_e).back()->Fill(nu_e,nu_e_reco,weight); 
        }
        if(nu_e_reco_nofsi > 0){
          h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).back()->Fill(nu_e,nu_e_reco_nofsi,weight); 
        }

      }

    }

   f->Close();

  }

  TFile* f_out = new TFile("ResponseMatricesNuMu.root","RECREATE");
  for(size_t i_e=0;i_e<kMAX;i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      h_TrueEnergy_RecoEnergy.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_NoFSI.at(i_e).at(i_f)->Write();
    }
  }
  f_out->Close(); 

}
