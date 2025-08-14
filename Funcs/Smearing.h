#ifndef _Smearing_h_
#define _Smearing_h_

#include "EnergyEstimatorFuncs.h"

namespace smearing {

TRandom2* rng = new TRandom2();

// List of momentum smearing resolutions (fractions)
const std::map<int,double> resolutions = {
  {0,0.20},   // 20% for calorimetric energy  
  {13,0.10},   // 10% for muons, from https://arxiv.org/pdf/1703.06187
  {12,0.15},   // 15% for electrons
  {2212,0.15}, // 15% for protons 
  {211,0.25},  // 25% for charged pions 
  {111,0.25},   // 25% for pi0
  {2112,1.0}   // 100% for neutrons
};


////////////////////////////////////////////////////////////////////////////////
// Smear the magnitude of momentum
// Calculate the surplus momentum above the detection threshold and smear this

void smear_mom(TLorentzVector& mom,int pdg){
  
  if(resolutions.find(abs(pdg)) == resolutions.end()) return;
  TVector3 mom3 = mom.Vect();
  double mass = mom.M();
  double r = std::max(rng->Gaus(1.0,resolutions.at(abs(pdg))),-1.0);
  mom3 *= r; 
  mom = TLorentzVector(mom3,sqrt(mass*mass+mom3.Mag()*mom3.Mag()));

}

void smear_ang(TLorentzVector& mom,int pdg){

  double mass = mom.M();

  TVector3 dir = mom.Vect().Unit();
  double ang = rng->Uniform(-3.14,3.14);
  double mix = rng->Gaus(0.0,0.02);
  TVector3 orth = dir.Orthogonal(); 
  orth.Rotate(ang,dir); 
  TVector3 dir2 = dir + mix*orth; 
  dir2 = dir2.Unit();  

  TVector3 p = mom.Vect().Mag()*dir2; 
  mom = TLorentzVector(p,sqrt(mass*mass+p.Mag()*p.Mag()));

}

////////////////////////////////////////////////////////////////////////////////
// Very handwavy efficiency effect 

// List of momentum smearing resolutions (fractions)
const std::map<int,double> efficiencies = {
  {2212,0.5}, // 75% for protons
  {211,1.0},  // 75% for charged pions 
  {111,1.0}   // 75% for pi0
};

void apply_eff(std::vector<int>& pdg,std::vector<TLorentzVector>& mom){

  std::vector<int> pdg_tmp;
  std::vector<TLorentzVector> mom_tmp;

  for(int i_p=0;i_p<pdg.size();i_p++){

    if(efficiencies.find(abs(pdg.at(i_p))) == efficiencies.end()) continue;;

    double r = rng->Uniform(0.0,1.0); 
    if(r < efficiencies.at(abs(pdg.at(i_p)))){
      pdg_tmp.push_back(pdg.at(i_p));
      mom_tmp.push_back(mom.at(i_p));
    } 

  }

  pdg = pdg_tmp;
  mom = mom_tmp;

}


}

#endif
