#ifndef _Funcs_h_
#define _Funcs_h_

const std::map<int,std::pair<double,double>> thresholds = {
  {13,{0.1,5.0}},
  {2212,{0.3,1.0}},
  {211,{0.1,5.0}},
  {111,{0.0,5.0}},
  {321,{0.2,5.0}}
};

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

enum estimators { kMuonKin , kMuonKinW , kMuonKinWNP , kPeLEELike , kPeLEELike0Pi , kTotalEDep , kMAX };
const std::vector<std::string> estimators_str = { "MuonKin" , "MuonKinW" , "MuonKinWNP" , "PeLEELike" , "PeLEELike0Pi"  , "TotalEDep" };

double T2KEnergy(const TLorentzVector* plepton){
  return (Mp*Mp - (Mn - Eb)*(Mn - Eb) - ml*ml + 2*(Mn - Eb)*plepton->E())/(2*(Mn - Eb - plepton->E() + plepton->P()*plepton->Vect().CosTheta()));
}

double T2KEnergyW(const TLorentzVector* plepton,const double& W){
  return (W*W - (Mn - Eb)*(Mn - Eb) - ml*ml + 2*(Mn - Eb)*plepton->E())/(2*(Mn - Eb - plepton->E() + plepton->P()*plepton->Vect().CosTheta()));
}

double ubooneEnergy(const TLorentzVector* plepton,const double& W,const int& nprot){
  return (W*W - nprot*nprot*(Mn - Eb)*(Mn - Eb) - ml*ml + nprot*2*(Mn - Eb)*plepton->E())/(2*(nprot*Mn - nprot*Eb - plepton->E() + plepton->P()*plepton->Vect().CosTheta()));
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

void Normalise(TH2D* h){
  for(int i_x=1;i_x<h->GetNbinsX()+1;i_x++){
    double content = 0.0;
    for(int i_y=1;i_y<h->GetNbinsY()+1;i_y++) content += h->GetBinContent(i_x,i_y);
    if(content  == 0.0) continue;
    for(int i_y=1;i_y<h->GetNbinsY()+1;i_y++) h->SetBinContent(i_x,i_y,h->GetBinContent(i_x,i_y)/content);
  }
}

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

#endif
