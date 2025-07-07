#ifndef _OscFitter_h_
#define _OscFitter_h_

// table of osc model constants
namespace c_osc {

  const double theta23 = 0.699834;
  const double theta13 = 0.160875;
  const double theta12 = 0.594371;

  const double sin2theta23 = sin(theta23)*sin(theta23); // sin^2(theta23)
  const double sin2theta13 = sin(theta13)*sin(theta13); // sin^2(theta13)
  const double sin2theta12 = sin(theta12)*sin(theta12); // sin^2(theta23)

  const double sin22theta23 = sin(2*theta23)*sin(2*theta23); // sin^2(2theta23)
  const double sin22theta13 = sin(2*theta13)*sin(2*theta13); // sin^2(2theta13)
  const double sin22theta12 = sin(2*theta12)*sin(2*theta12); // sin^2(2theta12)
  const double cos2theta23 = cos(theta23)*cos(theta23); // cos^2(theta23)
  const double cos2theta13 = cos(theta13)*cos(theta13); // cos^2(theta13)

  const double deltamsq12 = 0.759e-4; // delta m2 12 in eV^2
  const double deltamsq13 = 23.2e-4; // delta m2 13 in eV^2
  const double deltamsq23 = deltamsq13; // delta m2 23 in eV^2 (assuming normal hierarchy)

  const double L = 1285; // baseline in km
  const double GF = 1.166e-5; // Fermi constant in GeV^-2
  const double Ne = 2.9805e23; // Mean electrons/cm3 in earth crust from  https://arxiv.org/pdf/1707.02322 
  const double rho = 2.7; // mean density of earths crust in g/cm3
  const double aL = (rho/3.0)*L/3500; // matter effect term

}

class OscModel {

  public: 

  OscModel();

  double theta23 = c_osc::theta23; 
  double theta13 = c_osc::theta13; 
  double theta12 = c_osc::theta12; 

  double sin2theta23 = sin(theta23)*sin(theta23); // sin^2(theta23)
  double sin2theta13 = sin(theta13)*sin(theta13); // sin^2(theta13)
  double sin2theta12 = sin(theta12)*sin(theta12); // sin^2(theta23)
  double sin22theta23 = sin(2*theta23)*sin(2*theta23); // sin^2(2theta23)
  double sin22theta13 = sin(2*theta13)*sin(2*theta13); // sin^2(2theta13)
  double sin22theta12 = sin(2*theta12)*sin(2*theta12); // sin^2(2theta12)
  double cos2theta23 = cos(theta23)*cos(theta23); // cos^2(theta23)
  double cos2theta13 = cos(theta13)*cos(theta13); // cos^2(theta13)
  double cos2theta12 = cos(theta12)*cos(theta12); // cos^2(theta12)

  double deltamsq12 = c_osc::deltamsq12; // delta m2 12 in eV^2
  double deltamsq13 = c_osc::deltamsq13; // delta m2 13 in eV^2
  double deltamsq23 = c_osc::deltamsq23; // delta m2 23 in eV^2 (assuming normal hierarchy)

  double deltaCP = 0.0;

  double L = c_osc::L; // baseline in km
  double GF = c_osc::GF; // Fermi constant in GeV^-2
  double Ne = c_osc::Ne; // Mean electrons/cm3 in earth crust from  https://arxiv.org/pdf/1707.02322 
  double rho = c_osc::rho; // mean density of earths crust in g/cm3
  double aL = c_osc::aL; // matter effect term

  double numu_dis_amp = 4*cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23)*(1-cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23));

  void SetTheta23(double in_theta23);
  void SetTheta13(double in_theta13);
  void SetTheta12(double in_theta12);
  void SetDeltaMSq12(double in_deltamsq12);
  void SetDeltaMSq13(double in_deltamsq13);
  void SetDeltaMSq23(double in_deltamsq23);
  void SetDeltaCP(double in_deltaCP);
  void SetNuMuDisAmp(double amp);

  double NueAppProb(double E) const;
  double NuMuSurvProb(double E) const;

};

OscModel::OscModel(){

}

void OscModel::SetTheta23(double in_theta23){
  theta23 = in_theta23;
  sin2theta23 = sin(theta23)*sin(theta23);
  sin22theta23 = sin(2*theta23)*sin(2*theta23);
  cos2theta23 = cos(theta23)*cos(theta23);
  numu_dis_amp = 4*cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23)*(1-cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23));
}

void OscModel::SetTheta13(double in_theta13){
  theta13 = in_theta13;
  sin2theta13 = sin(theta13)*sin(theta13);
  sin22theta13 = sin(2*theta13)*sin(2*theta13);
  cos2theta13 = cos(theta13)*cos(theta13);
  numu_dis_amp = 4*cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23)*(1-cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23));
}

void OscModel::SetTheta12(double in_theta12){
  theta12 = in_theta12;
  sin2theta12 = sin(theta12)*sin(theta12);
  sin22theta12 = sin(2*theta12)*sin(2*theta12);
  cos2theta12 = cos(theta12)*cos(theta12);
  numu_dis_amp = 4*cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23)*(1-cos(theta13)*cos(theta13)*sin(theta23)*sin(theta23));
}

void OscModel::SetDeltaMSq12(double in_deltamsq12){
  deltamsq12 = in_deltamsq12;
}

void OscModel::SetDeltaMSq13(double in_deltamsq13){
  deltamsq13 = in_deltamsq13;
}

void OscModel::SetDeltaMSq23(double in_deltamsq23){
  deltamsq23 = in_deltamsq23;
}

void OscModel::SetDeltaCP(double in_deltaCP){
  deltaCP = in_deltaCP;
}

void OscModel::SetNuMuDisAmp(double amp){
  numu_dis_amp = amp;
}

////////////////////////////////////////////////////////////////////////////////
// nue appearance probability
// equation 1 from https://arxiv.org/pdf/2006.16043

double OscModel::NueAppProb(double E) const {

  double delta31 = 1.267*deltamsq13*L/E;
  double delta21 = 1.267*deltamsq12*L/E;

  double P = sin2theta23*sin22theta13*pow(sin(delta31 - aL),2)*delta31*delta31/pow(delta31 - aL,2)
    + sqrt(sin22theta23)*sqrt(sin22theta13)*sqrt(sin22theta12)*sin(delta31 - aL)*delta31/(delta31 - aL)
    *(sin(aL)/aL)*delta21*cos(delta31 + deltaCP)
    + cos2theta23*sin22theta12*pow(sin(aL),2)/aL/aL*delta21*delta21; 

  return P;

}

////////////////////////////////////////////////////////////////////////////////
// numu disappearance probability
// equation 1 from https://www.nature.com/articles/s41586-020-2177-0 
// Patrick Dunne (Patrick DUNE?) tells me this is good approximation at DUNE's baseline/matter effect

double OscModel::NuMuSurvProb(double E) const {
  double delta23 = 1.267*deltamsq23*L/E;
  return 1 - numu_dis_amp*sin(delta23)*sin(delta23);
}

#endif
