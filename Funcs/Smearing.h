#ifndef _Smearing_h_
#define _Smearing_h_

namespace smearing {

TRandom2* rng = new TRandom2();

// List of momentum smearing resolutions (fractions)
const std::map<int,double> resolutions = {
  {13,0.1},   // 10% for muons
  {2212,0.20}, // 20% for protons
  {211,0.25},  // 25% for charged pions 
  {111,0.25}   // 25% for pi0
};


////////////////////////////////////////////////////////////////////////////////
// Smear the magnitude of momentum
// Calculate the surplus momentum above the detection threshold and smear this
void smear_mom(TLorentzVector& mom,int pdg){

  // If particle below threshold or not in map, ignore
  if(thresholds.find(abs(pdg)) == thresholds.end() || mom.Vect().Mag() < thresholds.at(abs(pdg)).first) return;
  if(resolutions.find(abs(pdg)) == resolutions.end()) return;

  double threshold = thresholds.at(abs(pdg)).first; 
  TVector3 mom3 = mom.Vect();
  double mass = mom.M();

  double p = mom3.Mag() - threshold;

  double r = std::max(rng->Gaus(0.0,resolutions.at(abs(pdg))),-1.0);
  p += p*r; 

  double scale = (p + threshold)/mom3.Mag();
  mom3 *= scale;

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

}

#endif
