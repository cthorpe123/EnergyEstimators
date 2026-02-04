#include "../Funcs/EnergyEstimatorFuncs.h"
#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;
#include "../Funcs/Smearing.h"

const double _EPSILON_ = 0.01;

const std::string generator = "NEUT_Nue";
const int target_nu_pdg = 12;

void Filter(){

  TFile* f_in = TFile::Open(("/gluster/data/dune/cthorpe/DIS/"+generator+"Events.root").c_str());
  TTree* t_in = static_cast<TTree*>(f_in->Get("eventtree")) ;

  TFile* f_out = TFile::Open(("/gluster/data/dune/cthorpe/DIS/"+generator+"EventsFiltered.root").c_str(),"RECREATE");
  f_out->cd();

  Int_t ccnc;
  Int_t nu_pdg;  
  std::vector<int>* pdg=0;
  std::vector<TLorentzVector>* p4=0;
  Int_t lepton_pdg;
  TLorentzVector* lepton_p4=0;

  t_in->SetBranchAddress("nu_pdg",&nu_pdg);
  t_in->SetBranchAddress("ccnc",&ccnc);
  t_in->SetBranchAddress("pdg",&pdg);
  t_in->SetBranchAddress("p4",&p4);
  t_in->SetBranchAddress("lepton_pdg",&lepton_pdg); 
  t_in->SetBranchAddress("lepton_p4",&lepton_p4);

  TTree* t_out = t_in->CloneTree(0);

  int nprot;
  double W;
  std::vector<double> est_nu_e(kMAX);
  std::vector<double> smeared_est_nu_e(kMAX);
  t_out->Branch("W",&W);
  t_out->Branch("nprot",&nprot);
  t_out->Branch("est_nu_e",&est_nu_e);
  t_out->Branch("smeared_est_nu_e",&smeared_est_nu_e);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    t_in->GetEntry(ievent);
    //if(ievent > 100000) break;
    if(ievent % 50000 == 0) std::cout << "Event " << ievent << "/" << t_in->GetEntries() << std::endl;

    if(abs(nu_pdg) != target_nu_pdg || ccnc != 1) continue;
    //if(GetNProt(pdg,p4) < 1) continue;

    bool bad_event = false;
    for(size_t i_p=0;i_p<pdg->size();i_p++){
      if(masses.find(abs(pdg->at(i_p))) != masses.end()){
        if(abs(p4->at(i_p).M() - masses.at(abs(pdg->at(i_p)))) > _EPSILON_){
          //std::cout << pdg->at(i_p) << "  " << p4->at(i_p).M() << std::endl;
          bad_event = true;
        }
      }
    }
       
    if(bad_event) continue; 

    W = CalcW(pdg,p4);
    nprot = GetNProt(pdg,p4);
    est_nu_e = GetEnergyEst(lepton_p4,pdg,p4);

/*
    std::vector<TVector3> proton_mom = GetParticleMom(pdg,p4,2212);
    std::vector<TVector3> neutron_mom = GetNeutronMom(pdg,p4);
    std::vector<TVector3> pion_mom = GetParticleMom(pdg,p4,211);
    std::vector<TVector3> pizero_mom = GetParticleMom(pdg,p4,111);

    for(int i_e=0;i_e<kMAX;i_e++)
      est_nu_e.at(i_e) = GetEnergy(lepton_p4,W,nprot,proton_mom,pion_mom,pizero_mom,neutron_mom,i_e);
*/

    // Calculate predictions with kinematic smearing 
    TLorentzVector s_lepton_p4 = *lepton_p4;
    smearing::smear_mom(s_lepton_p4,13);
    std::vector<TLorentzVector> s_p4 = *p4;
    for(int i=0;i<p4->size();i++) smearing::smear_mom(s_p4.at(i),pdg->at(i));
    smeared_est_nu_e = GetEnergyEst(&s_lepton_p4,pdg,&s_p4);
 
    t_out->Fill();    

  } 

  std::cout << t_out->GetEntries() << std::endl;

  t_out->Write();
  f_out->Close();
  
  f_in->Close();

}
