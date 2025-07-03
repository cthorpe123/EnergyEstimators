#ifndef _Funcs_h_
#define _Funcs_h_

////////////////////////////////////////////////////////////////////////////////
// Assumed detection thresholds

const std::map<int,std::pair<double,double>> thresholds = {
  {13,{0.1,5.0}},
  {2212,{0.3,1.0}},
  {211,{0.1,5.0}},
  {111,{0.0,5.0}},
  {321,{0.2,5.0}}
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
// Different neutrino energy calculators 

enum estimators { kMuonKin , kMuonKinW , kMuonKinWNP , kPeLEELike , kPeLEELike0Pi , kTotalEDep , kSFMethod , kMAX };
const std::vector<std::string> estimators_str = { "MuonKin" , "MuonKinW" , "MuonKinWNP" , "PeLEELike" , "PeLEELike0Pi"  , "TotalEDep" , "SFMethod" };

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

double GetEnergy(const TLorentzVector* plepton,const double& W, const int& nprot,const std::vector<TVector3>& pprotons,const std::vector<TVector3>& ppions,const std::vector<TVector3>& ppizeros,int est){
  if(est < 0 || est >= kMAX) throw std::invalid_argument("GetEnergy: invalid estimator");

  switch(est){
    case kMuonKin: return T2KEnergy(plepton);
    case kMuonKinW: return T2KEnergyW(plepton,W);
    case kMuonKinWNP: return ubooneEnergy(plepton,W,nprot);
    case kPeLEELike: return peleeEnergy(plepton,pprotons);
    case kPeLEELike0Pi: if(!ppions.size() && !ppizeros.size()) return peleeEnergy(plepton,pprotons); else return -1;
    case kTotalEDep: return totaledepEnergy(plepton,pprotons,ppions,ppizeros);
    case kSFMethod: if(pprotons.size() == 1 && !ppions.size() && !ppizeros.size()) return sfmethodEnergy(plepton,pprotons); else return -1;
    default: return -1;
  }

  return -1;

}

////////////////////////////////////////////////////////////////////////////////
// Take a 2D histogram and normalise each vertical column to 1 

void Normalise(TH2D* h){
  for(int i_x=1;i_x<h->GetNbinsX()+1;i_x++){
    double content = 0.0;
    for(int i_y=1;i_y<h->GetNbinsY()+1;i_y++) content += h->GetBinContent(i_x,i_y);
    if(content  == 0.0) continue;
    for(int i_y=1;i_y<h->GetNbinsY()+1;i_y++) h->SetBinContent(i_x,i_y,h->GetBinContent(i_x,i_y)/content);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Get the bias and variance of a 2D hist as a function of the X axis variable 

void GetBiasVariance(const TH2D* h_Data,TH1D*& h_Bias,TH1D*& h_Variance){

  double nbins = h_Data->GetNbinsX();
  double low = h_Data->GetXaxis()->GetBinLowEdge(1);
  double high = h_Data->GetXaxis()->GetBinLowEdge(nbins+1);

  std::string name = h_Data->GetName();
  h_Bias = new TH1D((name + "_Bias").c_str(),";;Frac Bias",nbins,low,high);
  h_Variance = new TH1D((name + "_Variance").c_str(),";;Frac Variance",nbins,low,high);

  for(int i=1;i<nbins+1;i++){
    double mean = 0.0;
    double events = 0.0;
    for(int j=1;j<nbins+1;j++){
      mean += h_Data->GetYaxis()->GetBinCenter(j)*h_Data->GetBinContent(i,j);
      events += h_Data->GetBinContent(i,j);
    }

    h_Bias->SetBinContent(i,(mean/events - h_Bias->GetBinCenter(i))/h_Bias->GetBinCenter(i));
    if(events == 0.0) h_Bias->SetBinContent(i,0);

    double var = 0.0; 
    for(int j=1;j<nbins+1;j++){
      var += (h_Data->GetYaxis()->GetBinCenter(j) - mean)*(h_Data->GetYaxis()->GetBinCenter(j) - mean)*h_Data->GetBinContent(i,j);
    }

    h_Variance->SetBinContent(i,var/mean/mean/events/events);
    if(events == 0.0) h_Variance->SetBinContent(i,0);

  }

}

////////////////////////////////////////////////////////////////////////////////
// Calculate the covariance matrix from systematic variation histograms 

TMatrixDSym MakeCovariance(const std::vector<TH1D*> h_universes){

  int nbins = h_universes.at(0)->GetNbinsX();

  TMatrixDSym cov_mat(nbins);

  for(int i_b1=1;i_b1<nbins+1;i_b1++){
    for(int i_b2=1;i_b2<nbins+1;i_b2++){

      double mean1 = 0.0;
      for(size_t i_u=0;i_u<h_universes.size();i_u++){
        mean1 += h_universes.at(i_u)->GetBinContent(i_b1)/h_universes.size();
      }

      double mean2 = 0.0;
      for(size_t i_u=0;i_u<h_universes.size();i_u++){
        mean2 += h_universes.at(i_u)->GetBinContent(i_b2)/h_universes.size();
      }

      double cov = 0.0;
      for(size_t i_u=0;i_u<h_universes.size();i_u++){
        cov += (h_universes.at(i_u)->GetBinContent(i_b1) - mean1)*(h_universes.at(i_u)->GetBinContent(i_b2) - mean2);
      }

      cov /= h_universes.size();
      cov_mat[i_b1-1][i_b2-1] = cov;

    }
  }

  return cov_mat;

}

////////////////////////////////////////////////////////////////////////////////
// Take a true vs reco plot and draw the reco distribution, reweighted using
// a distribution in the true variable 

TH1D* MakeRewHist(const TH2D* h_true_reco,TH1D* h_true_ratio,std::string name){

  double nbins = h_true_reco->GetNbinsY();
  double low = h_true_reco->GetYaxis()->GetBinLowEdge(1);
  double high = h_true_reco->GetYaxis()->GetBinLowEdge(nbins+1);

  TH1D* h_reco = new TH1D(("h_rew_reco_"+name).c_str(),"",nbins,low,high);

  for(int i_br=1;i_br<h_true_reco->GetNbinsY()+1;i_br++){
    double events = 0.0;
    for(int i_bt=1;i_bt<h_true_reco->GetNbinsX()+1;i_bt++){
      double weight = h_true_ratio->GetBinContent(h_true_ratio->FindBin(h_true_reco->GetXaxis()->GetBinCenter(i_bt)));
      events += h_true_reco->GetBinContent(i_bt,i_br)*weight;
    }
    h_reco->SetBinContent(i_br,events);
  }

  return h_reco;

}

////////////////////////////////////////////////////////////////////////////////
// Divide histogram by its bin widths - convert events into diff cross sections

void DivideByBinWidth(TH1D* h) {
  int NBins = h->GetXaxis()->GetNbins();
  for (int i=1;i<NBins+1;i++){
    double CurrentEntry = h->GetBinContent(i);
    double NewEntry = CurrentEntry / h->GetBinWidth(i);
    double CurrentError = h->GetBinError(i);
    double NewError = CurrentError / h->GetBinWidth(i);
    h->SetBinContent(i,NewEntry); 
    h->SetBinError(i,NewError); 
  }
}

void DivideByBinWidth2D(TH2D* h) {
  int NBinsX = h->GetNbinsX();
  int NBinsY = h->GetNbinsY();
  for (int i=1;i<NBinsX+1;i++){
    for (int j=1;j<NBinsY+1;j++){
      double CurrentEntry = h->GetBinContent(i,j);
      double NewEntry = CurrentEntry / h->GetXaxis()->GetBinWidth(i) / h->GetYaxis()->GetBinWidth(j);
      double CurrentError = h->GetBinError(i,j);
      double NewError = CurrentError / h->GetXaxis()->GetBinWidth(i) / h->GetYaxis()->GetBinWidth(j);
      h->SetBinContent(i,j,NewEntry); 
      h->SetBinError(i,j,NewError); 
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Calculate rates from cross sections 

// Avagadro's number
const double NA = 6.022e23;

// number of mol of Ar40 in 1 ton
const double molperton = 1000/0.004;

// Gives rate per 1000 tons of detector mass per 10^21 POT
// Inputs are nu/cm2/10^21 POT
// XSec in 10^-38 cm^2
double Rate(double flux,double xsec){
  return 1e3*flux*xsec*NA*molperton*1e-38; 
}

#endif
