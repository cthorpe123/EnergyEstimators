#ifndef _Funcs_h_
#define _Funcs_h_

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

  int nbins_x = h_Data->GetNbinsX();
  int nbins_y = h_Data->GetNbinsY();

  for(int i=1;i<nbins_x+1;i++){
    
    double mean = 0.0;
    double events = 0.0;
    for(int j=1;j<nbins_y+1;j++){
      mean += h_Data->GetYaxis()->GetBinCenter(j)*h_Data->GetBinContent(i,j);
      events += h_Data->GetBinContent(i,j);
    }

    mean /= events;

    h_Bias->SetBinContent(i,(mean - h_Bias->GetBinCenter(i))/h_Bias->GetBinCenter(i));

    if(events == 0.0) h_Bias->SetBinContent(i,0);

    double var = 0.0; 
    for(int j=1;j<nbins_y+1;j++){
      var += (h_Data->GetYaxis()->GetBinCenter(j) - mean)*(h_Data->GetYaxis()->GetBinCenter(j) - mean)*h_Data->GetBinContent(i,j);
    }

    var /= events;    
    
    h_Variance->SetBinContent(i,var/mean/mean);
    if(events == 0.0) h_Variance->SetBinContent(i,0);

  }

}

////////////////////////////////////////////////////////////////////////////////
// 3D hist should have x and y axes as true and reconstructed energy
// z axis is variable of interest
void MakeBiasVarianceFrom3D(const TH3D* h, TH1D* h_bias,TH1D* h_variance){

  // Calculate the mean estimated energy afo true neutrino energy and third variable
  TH2D* h_energy = static_cast<TH2D*>(h->Project3D("zx"));  
  h_energy->Reset(); 
  for(int i=1;i<h->GetNbinsX()+1;i++){
    for(int k=1;k<h->GetNbinsZ()+1;k++){
     double mean = 0.0;
     double events = 0.0; 

     for(int j=1;j<h->GetNbinsY()+1;j++){
       mean += h->GetYaxis()->GetBinCenter(j)*h->GetBinContent(i,j,k);
       events += h->GetBinContent(i,j,k); 
     }
       
     h_energy->SetBinContent(i,k,mean/events);

    }
  } 

  // Calculate the bias and variance as function of third variable 
  double total_bias = 0.0;
  double total_variance = 0.0;
  double total_events = 0.0;
  for(int k=1;k<h->GetNbinsZ()+1;k++){
    double bias = 0.0;
    double events = 0.0;
    double variance = 0.0;
    for(int i=1;i<h->GetNbinsX()+1;i++){
      for(int j=1;j<h->GetNbinsY()+1;j++){
        if(h->GetBinContent(i,j,k) > 0){
          total_bias += (h->GetYaxis()->GetBinCenter(j) - h->GetXaxis()->GetBinCenter(i))*h->GetBinContent(i,j,k)/h->GetXaxis()->GetBinCenter(i);
          total_variance += h->GetBinContent(i,j,k)*pow((h->GetYaxis()->GetBinCenter(j) - h_energy->GetBinContent(i,k))/h_energy->GetBinContent(i,k),2);
          total_events += h->GetBinContent(i,j,k); 
          bias += (h->GetYaxis()->GetBinCenter(j) - h->GetXaxis()->GetBinCenter(i))*h->GetBinContent(i,j,k)/h->GetXaxis()->GetBinCenter(i);
          variance += h->GetBinContent(i,j,k)*pow((h->GetYaxis()->GetBinCenter(j) - h_energy->GetBinContent(i,k))/h_energy->GetBinContent(i,k),2);
          events += h->GetBinContent(i,j,k); 
        }
      }
    }
    //std::cout << "k=" << k << " events=" << events << " bias/events=" << bias/events << " variance/events=" << variance/events << std::endl;
    if(events > 0){
      h_bias->SetBinContent(k,bias/events);  
      h_variance->SetBinContent(k,variance/events);  
    }
  } 

  //std::cout << "total_bias=" << total_bias/total_events << " total_variance=" <<  total_variance/total_events << " total_events=" << total_events << std::endl;
  h_bias->SetBinContent(0,total_bias/total_events);
  h_variance->SetBinContent(0,total_variance/total_events);

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
//  Pseudo KS test used to evaulate level of similatity in the shapes of two 
//  histograms in a way that doesn't depend on specific shapes etc. 

double PseudoKSTest(const TH1D* h_model,const TH1D* h_data){

  TH1D* h_model_norm = static_cast<TH1D*>(h_model->Clone("h_model_norm"));
  TH1D* h_data_norm = static_cast<TH1D*>(h_data->Clone("h_data_norm"));
  h_model_norm->Scale(1.0/h_model->Integral());
  h_data_norm->Scale(1.0/h_data->Integral());

  TH1D* h_model_cum = static_cast<TH1D*>(h_model_norm->GetCumulative());
  TH1D* h_data_cum = static_cast<TH1D*>(h_data_norm->GetCumulative());

  double max_diff = 0.0;
  for(int i=1;i<h_model_cum->GetNbinsX()+1;i++)
    max_diff = std::max(abs(h_model_cum->GetBinContent(i) - h_data_cum->GetBinContent(i)),max_diff);
  

  delete h_model_norm;
  delete h_data_norm;
  delete h_model_cum;
  delete h_data_cum;
  
  return max_diff;

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


////////////////////////////////////////////////////////////////////////////////
// Convert double to string with fixed precision

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 1)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
} 

#endif
