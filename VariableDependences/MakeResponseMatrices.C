#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"

void MakeResponseMatrices(){

  bool nue_mode = true;

  std::vector<std::string> InputFiles_v;
  if(!nue_mode) InputFiles_v = {"GENIEEventsFiltered.root","NuWroEventsFiltered.root","NEUTEventsFiltered.root","GiBUUEventsFiltered.root"};
  else InputFiles_v = {"GENIE_NueEventsFiltered.root","NuWro_NueEventsFiltered.root","NEUT_NueEventsFiltered.root","GiBUU_NueEventsFiltered.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<double> bins;
  for(int i=0;i<201;i++) bins.push_back(0.2+i*(8.0-0.2)/200);
  int n_bins = bins.size()-1;
  double* bins_a = &bins[0];

  std::vector<double> w_binning = {0.8,1.0,1.2,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0}; 
  int w_points = w_binning.size()-1;
  double* w_binning_a = &w_binning[0]; 

  std::vector<double> ang_binning = {0.0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.142}; 
  int ang_points = ang_binning.size()-1;
  double* ang_binning_a = &ang_binning[0]; 

  std::vector<double> miss_binning = {0.0,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5}; 
  int miss_points = miss_binning.size()-1;
  double* miss_binning_a = &miss_binning[0]; 

  std::vector<double> n_binning = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5}; 
  int n_points = n_binning.size()-1;
  double* n_binning_a = &n_binning[0]; 

  // Muon Mom histogram setup
  std::vector<double> lepton_p_binning = {0.3,0.5,0.8,1.1,1.4,1.7,2.0,2.3,2.6,2.9,3.2,3.5,3.8,4.1,4.4,4.7,5.0,5.3,5.6,5.9,6.2,6.5,6.8,7.1,7.4,7.7,8.0}; 
  int lepton_p_points = lepton_p_binning.size()-1;
  double* lepton_p_binning_a = &lepton_p_binning[0]; 

  // First make the response matrices
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_W(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_Angle(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_MissingE(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_Neutrons(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_LeptonMom(kMAX);

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
      h_TrueEnergy_RecoEnergy_W.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_W_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);W_{vis} (GeV);",n_bins,bins_a,n_bins,bins_a,w_points,w_binning_a));
      h_TrueEnergy_RecoEnergy_Angle.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_Angle_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);#theta (rad);",n_bins,bins_a,n_bins,bins_a,ang_points,ang_binning_a));
      h_TrueEnergy_RecoEnergy_MissingE.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_MissingE_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);Missing Hadronic Energy (GeV);",n_bins,bins_a,n_bins,bins_a,miss_points,miss_binning_a));
      h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_Neutrons_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);FS Neutrons;",n_bins,bins_a,n_bins,bins_a,n_points,n_binning_a));
      h_TrueEnergy_RecoEnergy_LeptonMom.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_LeptonMom_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);p_{l} (GeV);",n_bins,bins_a,n_bins,bins_a,lepton_p_points,lepton_p_binning_a));
      h_TrueEnergy_RecoEnergy_W.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_Angle.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_MissingE.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_LeptonMom.at(i_e).back()->SetDirectory(0);
    }

    int target_nu_pdg = 14;
    if(nue_mode) target_nu_pdg = 12;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != target_nu_pdg || ccnc != 1) continue;

      double missing_e = GetMissingEnergy(pdg,p4); 
      double angle = acos(lepton_p4->CosTheta());
      int neutrons = 0;
      for(int i_p=0;i_p<pdg->size();i_p++)
        if(pdg->at(i_p) == 2112 && p4->at(i_p).Vect().Mag() > 0.3)
          neutrons++;
      double lepton_mom = lepton_p4->Vect().Mag();

      for(int i_e=0;i_e<kMAX;i_e++){
        h_TrueEnergy_RecoEnergy_W.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),W,weight); 
        h_TrueEnergy_RecoEnergy_MissingE.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),missing_e,weight); 
        h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),neutrons,weight); 
        h_TrueEnergy_RecoEnergy_Angle.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),angle,weight); 
        h_TrueEnergy_RecoEnergy_LeptonMom.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),lepton_mom,weight); 
      }

    }

    f->Close();

  }

  std::string filename = nue_mode ? "ResponseMatricesNue.root" : "ResponseMatricesNuMu.root";
  TFile* f_out = new TFile(filename.c_str(),"RECREATE");

  for(size_t i_e=0;i_e<kMAX;i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      h_TrueEnergy_RecoEnergy_W.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_Angle.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_MissingE.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_LeptonMom.at(i_e).at(i_f)->Write();
    }
  }

  f_out->Close(); 

}

