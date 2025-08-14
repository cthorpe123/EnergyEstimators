#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/OscFitter.h"
#include "../Funcs/Smearing.h"

void MakeResponseMatrices(){

  bool nue_mode = false;

  std::vector<std::string> InputFiles_v;
  if(!nue_mode) InputFiles_v = {"GENIEEventsFiltered.root"/*,"NuWroEventsFiltered.root","NEUTEventsFiltered.root","GiBUUEventsFiltered.root"*/};
  else InputFiles_v = {"GENIE_NueEventsFiltered.root","NuWro_NueEventsFiltered.root","NEUT_NueEventsFiltered.root","GiBUU_NueEventsFiltered.root"};
  std::vector<std::string> Generators_v = {"GENIE","NuWro","NEUT","GiBUU"};

  std::vector<double> bins;
  for(int i=0;i<201;i++) bins.push_back(0.2+i*(8.0-0.2)/200);
  int n_bins = bins.size()-1;
  double* bins_a = &bins[0];

  std::vector<double> n_binning = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5}; 
  int n_points = n_binning.size()-1;
  double* n_binning_a = &n_binning[0]; 

  std::vector<double> miss_binning = {0.0,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5}; 
  int miss_points = miss_binning.size()-1;
  double* miss_binning_a = &miss_binning[0]; 

  // First make the response matrices
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_Neutrons(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_Neutrons_Perfect(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_Neutrons_Smeared(kMAX);

  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_MissingE(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_MissingE_Perfect(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_MissingE_Smeared(kMAX);

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
      h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_Neutrons_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);FS Neutrons;",n_bins,bins_a,n_bins,bins_a,n_points,n_binning_a));
      h_TrueEnergy_RecoEnergy_Neutrons_Perfect.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_Neutrons_Perfect_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);FS Neutrons;",n_bins,bins_a,n_bins,bins_a,n_points,n_binning_a));
      h_TrueEnergy_RecoEnergy_Neutrons_Smeared.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_Neutrons_Smeared_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);FS Neutrons;",n_bins,bins_a,n_bins,bins_a,n_points,n_binning_a));
      h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_Neutrons_Perfect.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_Neutrons_Smeared.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_MissingE.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_MissingE_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);Missing Hadronic Energy (GeV);",n_bins,bins_a,n_bins,bins_a,miss_points,miss_binning_a));
      h_TrueEnergy_RecoEnergy_MissingE_Perfect.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_MissingE_Perfect_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);Missing Hadronic Energy (GeV);",n_bins,bins_a,n_bins,bins_a,miss_points,miss_binning_a));
      h_TrueEnergy_RecoEnergy_MissingE_Smeared.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_MissingE_Smeared_"+est).c_str(),";E_{true} (GeV);E_{est} (GeV);Missing Hadronic Energy (GeV);",n_bins,bins_a,n_bins,bins_a,miss_points,miss_binning_a));
      h_TrueEnergy_RecoEnergy_MissingE.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_MissingE_Perfect.at(i_e).back()->SetDirectory(0);
      h_TrueEnergy_RecoEnergy_MissingE_Smeared.at(i_e).back()->SetDirectory(0);
    }

    int target_nu_pdg = 14;
    if(nue_mode) target_nu_pdg = 12;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;

      if(nu_pdg != target_nu_pdg || ccnc != 1) continue;

      int neutrons = 0;
      for(int i_p=0;i_p<pdg->size();i_p++)
        if(pdg->at(i_p) == 2112 && p4->at(i_p).Vect().Mag() > 0.3)
          neutrons++;
      double missing_e = GetMissingEnergy(pdg,p4); 

      // Cheat way of adding neutrons - pretend they're protons
      std::vector<int> pdg_w_n = *pdg;
      for(int i_p=0;i_p<pdg_w_n.size();i_p++){
        if(pdg_w_n.at(i_p) == 2112) pdg_w_n.at(i_p) = 2212;
      }

      std::vector<double> energies_neutrons =  GetEnergyEst(lepton_p4,&pdg_w_n,p4);

      std::vector<int> pdg_w_n_smear = *pdg;
      std::vector<TLorentzVector> p4_w_n_smear = *p4; 

      for(int i_p=0;i_p<pdg_w_n_smear.size();i_p++){
        if(pdg_w_n_smear.at(i_p) == 2112){
          TLorentzVector& mom = p4_w_n_smear.at(i_p);
          TVector3 mom3 = mom.Vect().Unit()*(mom.Vect().Mag() - 0.3);
          double mass = mom.M();
          double r = std::max(smearing::rng->Gaus(1.0,smearing::resolutions.at(abs(2112))),-1.0);
          mom3 *= r; 
          mom3 += mom.Vect(); 
          mom = TLorentzVector(mom3,sqrt(mass*mass+mom3.Mag()*mom3.Mag()));
          pdg_w_n_smear.at(i_p) = 2212;
        } 
      }

      std::vector<double> energies_neutrons_smear =  GetEnergyEst(lepton_p4,&pdg_w_n_smear,&p4_w_n_smear);

      for(int i_e=0;i_e<kMAX;i_e++){
        h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),neutrons,weight); 
        h_TrueEnergy_RecoEnergy_Neutrons_Perfect.at(i_e).back()->Fill(nu_e,energies_neutrons.at(i_e),neutrons,weight); 
        h_TrueEnergy_RecoEnergy_Neutrons_Smeared.at(i_e).back()->Fill(nu_e,energies_neutrons_smear.at(i_e),neutrons,weight); 
        h_TrueEnergy_RecoEnergy_MissingE.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),missing_e,weight); 
        h_TrueEnergy_RecoEnergy_MissingE_Perfect.at(i_e).back()->Fill(nu_e,energies_neutrons.at(i_e),missing_e,weight); 
        h_TrueEnergy_RecoEnergy_MissingE_Smeared.at(i_e).back()->Fill(nu_e,energies_neutrons_smear.at(i_e),missing_e,weight); 
      }
    }

    f->Close();

  }

  std::string filename = nue_mode ? "ResponseMatricesNue.root" : "ResponseMatricesNuMu.root";
  TFile* f_out = new TFile(filename.c_str(),"RECREATE");

  for(size_t i_e=0;i_e<kMAX;i_e++){
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){
      h_TrueEnergy_RecoEnergy_Neutrons.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_Neutrons_Perfect.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_Neutrons_Smeared.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_MissingE.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_MissingE_Perfect.at(i_e).at(i_f)->Write();
      h_TrueEnergy_RecoEnergy_MissingE_Smeared.at(i_e).at(i_f)->Write();
    }
  }

  f_out->Close(); 

}

