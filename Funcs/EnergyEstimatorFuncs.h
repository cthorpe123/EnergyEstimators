#ifndef _EnergyEstimatorFuncs_h_
#define _EnergyEstimatorFuncs_h_

////////////////////////////////////////////////////////////////////////////////
// Assumed detection thresholds

const std::map<int,std::pair<double,double>> thresholds = {
  {13,{0.1,5.0}},
  {2212,{0.3,2.0}},
  {211,{0.1,5.0}},
  {111,{0.0,5.0}}
  //{321,{0.2,5.0}}
};

////////////////////////////////////////////////////////////////////////////////
// Table of particle masses 

const std::map<int,double> masses = {
  {13,0.106},
  {2112,0.938},
  {2212,0.935},
  {211,0.140},
  {111,0.135}
};

const double Mp = 0.935;
const double Mn = 0.938;
const double Eb = 0.03;
const double ml = 0.106;
const double mpi = 0.140;
const double mpi0 = 0.135;
const double MA = 22*Mn + 18*Mp - 0.34381;
const double MA1 = MA - Mn;

////////////////////////////////////////////////////////////////////////////////
// Visible hadronic invariant mass 

double CalcW(const std::vector<int>* pdg_v,const std::vector<TLorentzVector>* p4_v){

  TLorentzVector p4_tot(0,0,0,0); 
  for(int i_p=0;i_p<pdg_v->size();i_p++){
    const int& pdg = pdg_v->at(i_p);
    const double p = p4_v->at(i_p).Vect().Mag();   
    if(thresholds.find(abs(pdg)) != thresholds.end() && p > thresholds.at(abs(pdg)).first && p < thresholds.at(abs(pdg)).second)
      p4_tot += p4_v->at(i_p);
  }

  return p4_tot.M();

}

////////////////////////////////////////////////////////////////////////////////
// Get number of visible protons 

int GetNProt(const std::vector<int>* pdg_v,const std::vector<TLorentzVector>* p4_v){

  int nprot = 0;

  for(int i_p=0;i_p<pdg_v->size();i_p++){
    const int& pdg = pdg_v->at(i_p);
    const double p = p4_v->at(i_p).Vect().Mag();   
    if(pdg != 2212) continue;
    if(p > thresholds.at(2212).first && p < thresholds.at(2212).second)
      nprot++; 
  }

  return nprot;
}

////////////////////////////////////////////////////////////////////////////////
// Get the 3 momenta in descending order of magniture of a given particle type

std::vector<TVector3> GetParticleMom(const std::vector<int>* pdg_v,const std::vector<TLorentzVector>* p4_v,const int target_pdg){

  std::vector<TVector3> mom;
  for(int i_p=0;i_p<pdg_v->size();i_p++){
    if(abs(pdg_v->at(i_p)) != target_pdg) continue;
    const double p = p4_v->at(i_p).Vect().Mag();   
    if(p > thresholds.at(target_pdg).first && p < thresholds.at(target_pdg).second)
      mom.push_back(p4_v->at(i_p).Vect());
  } 

  std::sort(mom.begin(),mom.end(),[](const TVector3& a,const TVector3& b){ return a.Mag() > b.Mag(); });

  return mom;

}

////////////////////////////////////////////////////////////////////////////////
// Get the 3 momenta of neutrons with momentum above 0.4 GeV 
// Threshold is estimated (very badly) from Fig. 4b of https://arxiv.org/pdf/2406.10583
std::vector<TVector3> GetNeutronMom(const std::vector<int>* pdg_v,const std::vector<TLorentzVector>* p4_v,double neutron_thresh = 0.3){

  std::vector<TVector3> mom;
  for(int i_p=0;i_p<pdg_v->size();i_p++){
    if(abs(pdg_v->at(i_p)) != 2112 || p4_v->at(i_p).Vect().Mag() < neutron_thresh) continue;
    mom.push_back(p4_v->at(i_p).Vect());
  } 

  std::sort(mom.begin(),mom.end(),[](const TVector3& a,const TVector3& b){ return a.Mag() > b.Mag(); });

  return mom;

}

////////////////////////////////////////////////////////////////////////////////
// Different neutrino energy calculators 

enum estimators { kMuonKin , kMuonKinWNP , kPeLEELike0Pi , kTotalEDep , kSFMethod , kMAX };
const std::vector<std::string> estimators_str = { "MuonKin" , "MuonKinWNP" , "PeLEELike0Pi"  , "TotalEDep" , "SFMethod" };

double T2KEnergy(const TLorentzVector* plepton){
  return (Mp*Mp - (Mn - Eb)*(Mn - Eb) - plepton->M()*plepton->M() + 2*(Mn - Eb)*plepton->E())/(2*(Mn - Eb - plepton->E() + plepton->P()*plepton->Vect().CosTheta()));
}

double T2KEnergyW(const TLorentzVector* plepton,const double& W){
  return (W*W - (Mn - Eb)*(Mn - Eb) - plepton->M()*plepton->M() + 2*(Mn - Eb)*plepton->E())/(2*(Mn - Eb - plepton->E() + plepton->P()*plepton->Vect().CosTheta()));
}

double ubooneEnergy(const TLorentzVector* plepton,const double& W,const int& nprot){
  return (W*W - nprot*nprot*(Mn - Eb)*(Mn - Eb) - plepton->M()*plepton->M() + nprot*2*(Mn - Eb)*plepton->E())/(2*(nprot*Mn - nprot*Eb - plepton->E() + plepton->P()*plepton->Vect().CosTheta()));
}

double peleeEnergy(const TLorentzVector* plepton,const std::vector<TVector3>& pprotons){
  double e = plepton->E();
  for(TVector3 p : pprotons) e += sqrt(p.Mag()*p.Mag() + Mp*Mp) - Mp;
  return e; 
}

double totaledepEnergy(const TLorentzVector* plepton,const std::vector<TVector3>& pprotons,const std::vector<TVector3>& ppions,const std::vector<TVector3>& ppizeros){
  double e = plepton->E();
  for(TVector3 p : pprotons) e += sqrt(p.Mag()*p.Mag() + Mp*Mp) - Mp;
  for(TVector3 p : ppions) e += sqrt(p.Mag()*p.Mag() + mpi*mpi);
  for(TVector3 p : ppizeros) e += sqrt(p.Mag()*p.Mag() + mpi0*mpi0);
  return e;
}

double sfmethodEnergy(const TLorentzVector* plepton,const std::vector<TVector3>& pprotons){
  const TVector3& prot = pprotons.at(0);
  const double Ep = sqrt(prot.Mag()*prot.Mag() + Mp*Mp); 
  double pt2 = (prot + plepton->Vect()).Perp2();
  double pL = (pow((MA + plepton->Pz() + prot.Z() - plepton->E() - Ep),2) - pt2 - MA1*MA1)/2/(MA + plepton->Pz() + prot.Z() - plepton->E() - Ep); 

  return plepton->Pz() + prot.Z() - pL;

}

// W Method with detected neutrons
double MuonKinWNPNeutronEnergy(const TLorentzVector* plepton,const std::vector<TVector3>& pprotons,const std::vector<TVector3>& ppions,const std::vector<TVector3>& ppizeros,const std::vector<TVector3>& pneutrons){

  TLorentzVector p4_tot(0,0,0,0); 
  for(TVector3 proton : pprotons) p4_tot += TLorentzVector(proton,sqrt(Mp*Mp+proton.Mag()*proton.Mag()));
  for(TVector3 pion : ppions) p4_tot += TLorentzVector(pion,sqrt(mpi*mpi+pion.Mag()*pion.Mag()));
  for(TVector3 pizero : ppizeros) p4_tot += TLorentzVector(pizero,sqrt(mpi0*mpi0+pizero.Mag()*pizero.Mag()));
  for(TVector3 neutron : pneutrons) p4_tot += TLorentzVector(neutron,sqrt(Mn*Mn+neutron.Mag()*neutron.Mag()));
  double W = p4_tot.M();
  int n = pprotons.size() + pneutrons.size();
  return ubooneEnergy(plepton,W,n); 

}

double totaledepNeutronEnergy(const TLorentzVector* plepton,const std::vector<TVector3>& pprotons,const std::vector<TVector3>& ppions,const std::vector<TVector3>& ppizeros,const std::vector<TVector3>& pneutrons){
  double e = plepton->E();
  for(TVector3 p : pprotons) e += sqrt(p.Mag()*p.Mag() + Mp*Mp) - Mp;
  for(TVector3 p : ppions) e += sqrt(p.Mag()*p.Mag() + mpi*mpi);
  for(TVector3 p : ppizeros) e += sqrt(p.Mag()*p.Mag() + mpi0*mpi0);
  for(TVector3 p : pneutrons) e += sqrt(p.Mag()*p.Mag() + Mn*Mn) - Mn;

  return e;
}


double GetEnergy(const TLorentzVector* plepton,const double& W, const int& nprot,const std::vector<TVector3>& pprotons,const std::vector<TVector3>& ppions,const std::vector<TVector3>& ppizeros,const std::vector<TVector3>& pneutrons,int est){
  if(est < 0 || est >= kMAX) throw std::invalid_argument("GetEnergy: invalid estimator");

  switch(est){
    case kMuonKin: return T2KEnergy(plepton);
    case kMuonKinWNP: return ubooneEnergy(plepton,W,nprot);
    case kPeLEELike0Pi: if(!ppions.size() && !ppizeros.size()) return peleeEnergy(plepton,pprotons); else return -1;
    case kTotalEDep: return totaledepEnergy(plepton,pprotons,ppions,ppizeros);
    case kSFMethod: if(pprotons.size() == 1 && !ppions.size() && !ppizeros.size()) return sfmethodEnergy(plepton,pprotons); else return -1;
    default: return -1;
  }

  return -1;

}

#endif
