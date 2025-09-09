#include "../Funcs/Funcs.h"
#include "../Funcs/EnergyEstimatorFuncs.h"
#include "../Funcs/Smearing.h"

#include "TLorentzVector.h"

bool nue_mode = false;

void FunWithNeutrons(){

  std::vector<std::string> InputFiles_v = {"GENIEEventsFiltered.root","NEUTEventsFiltered.root","GiBUUEventsFiltered.root","NuWroEventsFiltered.root"};
  std::vector<std::string> Generators_v = {"GENIE","NEUT","GiBUU","NuWro"};

  std::vector<double> binning_v;
  binning_v.push_back(0.35);
  for(int i=0;i<31;i++) binning_v.push_back(binning_v.back()+0.15);
  for(int i=0;i<10;i++) binning_v.push_back(binning_v.back()+0.3);
  int nbins = binning_v.size()-1;
  double* binning_a = &binning_v[0];

  double n_binning_a[] = {-0.5,0.5,1.5,2.5};

  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy(kMAX);
  std::vector<std::vector<TH3D*>> h_TrueEnergy_RecoEnergy_Coarse(kMAX);
  std::vector<std::vector<TH1D*>> h_Bias_NoNeutrons(kMAX);
  std::vector<std::vector<TH1D*>> h_Bias_PerfectNeutrons(kMAX);
  std::vector<std::vector<TH1D*>> h_Bias_SmearedNeutrons(kMAX);
  std::vector<std::vector<TH1D*>> h_Change_PerfectNeutrons(kMAX);
  std::vector<std::vector<TH1D*>> h_Change_SmearedNeutrons(kMAX);

  for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

    std::string generator = Generators_v.at(i_f);    

    TFile* f = TFile::Open(("/gluster/data/dune/cthorpe/DIS/"+InputFiles_v.at(i_f)).c_str());
    TTree* t = static_cast<TTree*>(f->Get("eventtree")) ;

    Double_t scale;
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

    t->SetBranchAddress("scale",&scale);
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

    for(int i_e=0;i_e<estimators_str.size();i_e++){
      std::string est = estimators_str.at(i_e);
      h_TrueEnergy_RecoEnergy.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_"+est).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);",200,0.2,8.0,200,0.2,8.0,3,-0.5,2.5));
      h_TrueEnergy_RecoEnergy_Coarse.at(i_e).push_back(new TH3D((generator+"_TrueEnergy_RecoEnergy_Coarse_"+est).c_str(),";True Neutrino Energy (GeV);Estimated Neutrino Energy (GeV);",nbins,binning_a,nbins,binning_a,3,n_binning_a));
      h_Bias_NoNeutrons.at(i_e).push_back(new TH1D((generator+"_Bias_NoNeutrons_"+est).c_str(),"",200,-0.5,0.5));
      h_Bias_PerfectNeutrons.at(i_e).push_back(new TH1D((generator+"_Bias_PerfectNeutrons_"+est).c_str(),"",200,-0.5,0.5));
      h_Bias_SmearedNeutrons.at(i_e).push_back(new TH1D((generator+"_Bias_SmearedNeutrons_"+est).c_str(),"",200,-0.5,0.5));

      h_Change_PerfectNeutrons.at(i_e).push_back(new TH1D((generator+"_Change_PerfectNeutrons_"+est).c_str(),"",100,-0.5,0.8));
      h_Change_SmearedNeutrons.at(i_e).push_back(new TH1D((generator+"_Change_SmearedNeutrons_"+est).c_str(),"",100,-0.5,0.8));
    }

    int target_nu_pdg = 14;

    double no_neutrons = 0;
    double neutrons = 0;

    for(Long64_t ievent=0;ievent<t->GetEntries();ievent++){

      //if(ievent > 100000) break;
      if(ievent % 20000 == 0) std::cout << generator << " Event " << ievent << "/" << t->GetEntries() << std::endl;
      t->GetEntry(ievent);

      if(generator != "GiBUU") weight = 1.0;
      weight *= scale*1e38*40;

      if(abs(nu_pdg) != target_nu_pdg || ccnc != 1) continue;
      if(nprot < 1) continue;

      bool has_neutron = false;
      int n_neutron = 0;
      for(int i_p=0;i_p<pdg->size();i_p++)
        if(pdg->at(i_p) == 2112 && p4->at(i_p).Vect().Mag() > 0.3){
          n_neutron++;
        }
      has_neutron = n_neutron > 0;       

      if(has_neutron) neutrons += weight;
      else no_neutrons += weight;

      // Cheat way of adding neutrons - pretend they're protons
      std::vector<int> pdg_w_n = *pdg;
      for(int i_p=0;i_p<pdg_w_n.size();i_p++){
        if(pdg_w_n.at(i_p) == 2112) pdg_w_n.at(i_p) = 2212;
      }

      std::vector<double> energies_neutrons = GetEnergyEst(lepton_p4,&pdg_w_n,p4);

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

        h_TrueEnergy_RecoEnergy.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),0.0,weight);
        h_TrueEnergy_RecoEnergy.at(i_e).back()->Fill(nu_e,energies_neutrons.at(i_e),1.0,weight);
        h_TrueEnergy_RecoEnergy.at(i_e).back()->Fill(nu_e,energies_neutrons_smear.at(i_e),2.0,weight);

        h_TrueEnergy_RecoEnergy_Coarse.at(i_e).back()->Fill(nu_e,est_nu_e->at(i_e),0.0,weight);
        h_TrueEnergy_RecoEnergy_Coarse.at(i_e).back()->Fill(nu_e,energies_neutrons.at(i_e),1.0,weight);
        h_TrueEnergy_RecoEnergy_Coarse.at(i_e).back()->Fill(nu_e,energies_neutrons_smear.at(i_e),2.0,weight);

        h_Bias_NoNeutrons.at(i_e).back()->Fill((est_nu_e->at(i_e)-nu_e)/nu_e,weight);
        h_Bias_PerfectNeutrons.at(i_e).back()->Fill((energies_neutrons.at(i_e)-nu_e)/nu_e,weight);
        h_Bias_SmearedNeutrons.at(i_e).back()->Fill((energies_neutrons_smear.at(i_e)-nu_e)/nu_e,weight);

        //if(abs(energies_neutrons.at(i_e) - est_nu_e->at(i_e)) > 0){
        if(has_neutron && est_nu_e->at(i_e) > 0 && energies_neutrons.at(i_e) > 0)
          h_Change_PerfectNeutrons.at(i_e).back()->Fill((energies_neutrons.at(i_e) - est_nu_e->at(i_e))/est_nu_e->at(i_e),weight);


        if(has_neutron && est_nu_e->at(i_e) > 0 && energies_neutrons_smear.at(i_e) > 0)
          //if(abs(energies_neutrons_smear.at(i_e) - est_nu_e->at(i_e)) > 0)
          h_Change_SmearedNeutrons.at(i_e).back()->Fill((energies_neutrons_smear.at(i_e) - est_nu_e->at(i_e))/est_nu_e->at(i_e),weight);

      }

      }

      std::cout << "Fraction of events with neutrons = " << neutrons/(neutrons+no_neutrons) << std::endl;
      std::cout << "Fraction of events without neutrons = " << no_neutrons/(neutrons+no_neutrons) << std::endl;

    }

    gSystem->Exec("mkdir -p Plots/");
    TCanvas* c = new TCanvas("c","c");
    TLegend* l = new TLegend(0.73,0.73,0.98,0.98);

    // Overall bias and variance as a function of true neutrino energy
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);
      THStack* hs_Bias = new THStack(("hs_Bias"+gen).c_str(),";;Frac. Bias");
      THStack* hs_Variance = new THStack(("hs_Variance"+gen).c_str(),";;Frac. Variance");
      std::vector<TH1D*> h_Bias;
      std::vector<TH1D*> h_Variance;

      for(size_t i_e=0;i_e<kMAX;i_e++){

        if(i_e == kMuonKin || i_e == kSFMethod) continue;

        std::string est = estimators_str.at(i_e);

        const TH3D* h = h_TrueEnergy_RecoEnergy.at(i_e).at(i_f);
        h_Bias.push_back(new TH1D(("h_Bias_"+gen+"_"+est).c_str(),"",3,-0.5,2.5));
        h_Variance.push_back(new TH1D(("h_Variance_"+gen+"_"+est).c_str(),"",3,-0.5,2.5));

        MakeBiasVarianceFrom3D(h,h_Bias.back(),h_Variance.back());

        h_Bias.back()->SetLineColor(colors.at(i_e));
        h_Bias.back()->SetLineWidth(2);
        hs_Bias->Add(h_Bias.back());

        h_Variance.back()->SetLineColor(colors.at(i_e));
        h_Variance.back()->SetLineWidth(2);
        hs_Variance->Add(h_Variance.back());

        l->AddEntry(h_Bias.back(),est.c_str(),"L");

      }

      hs_Bias->Draw("nostack HIST");
      hs_Bias->GetXaxis()->SetBinLabel(1,"No Neutrons");
      hs_Bias->GetXaxis()->SetBinLabel(2,"Visible Neutrons");
      hs_Bias->GetXaxis()->SetBinLabel(3,"Smeared Neutrons");
      l->Draw();
      c->Print(("Plots/Bias_"+gen+".pdf").c_str());  
      c->Clear();

      hs_Variance->Draw("nostack HIST");
      hs_Variance->GetXaxis()->SetBinLabel(1,"No Neutrons");
      hs_Variance->GetXaxis()->SetBinLabel(2,"Visible Neutrons");
      hs_Variance->GetXaxis()->SetBinLabel(3,"Smeared Neutrons");
      l->Draw();
      c->Print(("Plots/Variance_"+gen+".pdf").c_str());  
      c->Clear();
      l->Clear();

    }

    // Bias shape
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      THStack* hs = new THStack("hs","");

      for(size_t i_e=0;i_e<kMAX;i_e++){

        if(i_e == kMuonKin || i_e == kSFMethod) continue;

        std::string est = estimators_str.at(i_e);

        TH1D* h_orig = h_Bias_NoNeutrons.at(i_e).at(i_f);
        TH1D* h_perfect = h_Bias_PerfectNeutrons.at(i_e).at(i_f);
        TH1D* h_smeared = h_Bias_SmearedNeutrons.at(i_e).at(i_f);
        h_orig->Scale(1.0/h_orig->Integral());
        h_perfect->Scale(1.0/h_perfect->Integral());
        h_smeared->Scale(1.0/h_smeared->Integral());

        l->AddEntry(h_orig,est.c_str(),"L");

        h_orig->SetLineColor(colors.at(i_e));
        h_orig->SetLineWidth(2);
        hs->Add(h_orig);

        h_perfect->SetLineColor(colors.at(i_e)+2);
        h_perfect->SetLineWidth(2);
        hs->Add(h_perfect);

        h_smeared->SetLineColor(colors.at(i_e)+2);
        h_smeared->SetLineWidth(2);
        hs->Add(h_smeared);

      }

      hs->Draw("nostack HIST");
      l->Draw();
      c->Print(("Plots/BiasShape_"+gen+".pdf").c_str());  
      c->Clear();
      l->Clear();

      delete hs;

    }

    // Bias shape as one plot per estimator/generator
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      for(size_t i_e=0;i_e<kMAX;i_e++){

        if(i_e == kMuonKin || i_e == kSFMethod) continue;

        THStack* hs = new THStack("hs",";(E_{est} - E_{true})/E_{true};Events");

        std::string est = estimators_str.at(i_e);

        TH1D* h_orig = h_Bias_NoNeutrons.at(i_e).at(i_f);
        TH1D* h_perfect = h_Bias_PerfectNeutrons.at(i_e).at(i_f);
        TH1D* h_smeared = h_Bias_SmearedNeutrons.at(i_e).at(i_f);
        h_orig->Scale(1.0/h_orig->Integral());
        h_perfect->Scale(1.0/h_perfect->Integral());
        h_smeared->Scale(1.0/h_smeared->Integral());

        h_orig->SetLineColor(colors.at(i_e));
        h_orig->SetLineWidth(2);
        hs->Add(h_orig);
        l->AddEntry(h_orig,"Original","L");

        h_perfect->SetLineColor(colors.at(i_e));
        h_perfect->SetLineWidth(2);
        h_perfect->SetLineStyle(2);
        hs->Add(h_perfect);
        l->AddEntry(h_perfect,"Perfect","L");

        h_smeared->SetLineColor(colors.at(i_e));
        h_smeared->SetLineStyle(3);
        h_smeared->SetLineWidth(2);
        hs->Add(h_smeared);
        l->AddEntry(h_smeared,"Smeared","L");

        hs->Draw("nostack HIST");
        l->Draw();
        c->Print(("Plots/BiasShape_"+est+"_"+gen+".pdf").c_str());  
        c->Clear();
        l->Clear();

        delete hs;

      }

    }

    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      for(size_t i_e=0;i_e<kMAX;i_e++){

        THStack* hs = new THStack("hs",";(E_{est}^{n} - E_{est})/E_{est};Events");

        if(i_e == kMuonKin || i_e == kSFMethod) continue;
        std::string est = estimators_str.at(i_e);

        TH1D* h_perfect = h_Change_PerfectNeutrons.at(i_e).at(i_f);
        TH1D* h_smeared = h_Change_SmearedNeutrons.at(i_e).at(i_f);
        h_perfect->Scale(1.0/h_perfect->Integral());
        h_smeared->Scale(1.0/h_smeared->Integral());

        h_perfect->SetLineColor(colors.at(i_e));
        h_perfect->SetLineStyle(2);
        h_perfect->SetLineWidth(2);
        hs->Add(h_perfect);
        l->AddEntry(h_perfect,(est+" Perfect").c_str(),"L");

        h_smeared->SetLineColor(colors.at(i_e));
        h_smeared->SetLineStyle(3);
        h_smeared->SetLineWidth(2);
        hs->Add(h_smeared);
        l->AddEntry(h_smeared,(est+" Smeared").c_str(),"L");

        hs->Draw("nostack HIST");
        l->Draw();
        c->Print(("Plots/Change_"+est+"_"+gen+".pdf").c_str());  
        c->Clear();
        l->Clear();

        delete hs;

      }

    }

    // Bias and variance as a function of neutrino energy for all three methods
    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      THStack* hs_Bias = new THStack("hs",";True Neutrino Energy (GeV);Frac. Bias");
      THStack* hs_Variance = new THStack("hs",";True Neutrino Energy (GeV);Frac. Variance");
      std::vector<TH1D*> h_Bias,h_Variance;

      for(size_t i_e=0;i_e<kMAX;i_e++){

        if(i_e == kMuonKin || i_e == kSFMethod) continue;
        std::string est = estimators_str.at(i_e);

        TH3D* h = h_TrueEnergy_RecoEnergy_Coarse.at(i_e).at(i_f);

        h->GetZaxis()->SetRange(1,1); 
        TH2D* h_orig =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_orig")); 
        h_Bias.push_back(new TH1D(("h_Bias_Orig_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        h_Variance.push_back(new TH1D(("h_Variance_Orig_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        GetBiasVariance(h_orig,h_Bias.back(),h_Variance.back()); 
        hs_Bias->Add(h_Bias.back());    
        hs_Variance->Add(h_Variance.back());    
        h_Bias.back()->SetLineColor(colors.at(i_e));
        h_Bias.back()->SetLineWidth(2);
        h_Variance.back()->SetLineColor(colors.at(i_e));
        h_Variance.back()->SetLineWidth(2);

        l->AddEntry(h_Bias.back(),est.c_str(),"L");

        h->GetZaxis()->SetRange(2,2); 
        TH2D* h_perfect =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_perfect")); 
        h_Bias.push_back(new TH1D(("h_Bias_Perfect_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        h_Variance.push_back(new TH1D(("h_Variance_Perfect_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        GetBiasVariance(h_perfect,h_Bias.back(),h_Variance.back()); 
        hs_Bias->Add(h_Bias.back());    
        hs_Variance->Add(h_Variance.back());    
        h_Bias.back()->SetLineColor(colors.at(i_e));
        h_Bias.back()->SetLineStyle(2);
        h_Bias.back()->SetLineWidth(2);
        h_Variance.back()->SetLineColor(colors.at(i_e));
        h_Variance.back()->SetLineStyle(2);
        h_Variance.back()->SetLineWidth(2);

        h->GetZaxis()->SetRange(3,3); 
        TH2D* h_smeared =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_smeared")); 
        h_Bias.push_back(new TH1D(("h_Bias_Smeared_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        h_Variance.push_back(new TH1D(("h_Variance_Smeared_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        GetBiasVariance(h_smeared,h_Bias.back(),h_Variance.back()); 
        hs_Bias->Add(h_Bias.back());    
        hs_Variance->Add(h_Variance.back());    
        h_Bias.back()->SetLineColor(colors.at(i_e));
        h_Bias.back()->SetLineStyle(3);
        h_Bias.back()->SetLineWidth(2);
        h_Variance.back()->SetLineColor(colors.at(i_e));
        h_Variance.back()->SetLineStyle(3);
        h_Variance.back()->SetLineWidth(2);

      }

      hs_Bias->Draw("nostack HIST"); 
      l->Draw();
      c->Print(("Plots/Bias_Energy_"+gen+".pdf").c_str());  
      c->Clear();

      hs_Variance->Draw("nostack HIST"); 
      l->Draw();
      c->Print(("Plots/Variance_Energy_"+gen+".pdf").c_str());  
      c->Clear();

      l->Clear();

    }


    for(size_t i_f=0;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      THStack* hs_Bias_Ratio  = new THStack("hs",";True Neutrino Energy (GeV);Frac. Bias With N - Frac. Bias Without N)");
      THStack* hs_Variance_Ratio  = new THStack("hs",";True Neutrino Energy (GeV);Frac. Variance With N - Frac. Variance Without N)");
      std::vector<TH1D*> h_Bias_Ratio;
      std::vector<TH1D*> h_Variance_Ratio;

      for(size_t i_e=0;i_e<kMAX;i_e++){

        if(i_e == kMuonKin || i_e == kSFMethod) continue;
        std::string est = estimators_str.at(i_e);

        TH3D* h = h_TrueEnergy_RecoEnergy_Coarse.at(i_e).at(i_f);

        h->GetZaxis()->SetRange(1,1); 
        TH2D* h_orig =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_orig2")); 
        TH1D* h_bias_orig = new TH1D(("h_Bias_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        TH1D* h_variance_orig = new TH1D(("h_Variance_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        GetBiasVariance(h_orig,h_bias_orig,h_variance_orig); 

        h->GetZaxis()->SetRange(2,2); 
        TH2D* h_perfect =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_perfect2")); 
        TH1D* h_bias_perfect = new TH1D(("h_Bias_Perfect2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        TH1D* h_variance_perfect = new TH1D(("h_Variance_Perfect2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        GetBiasVariance(h_perfect,h_bias_perfect,h_variance_perfect); 

        h->GetZaxis()->SetRange(3,3); 
        TH2D* h_smeared =  static_cast<TH2D*>(h->Project3D("yx")->Clone("h_smeared2")); 
        TH1D* h_bias_smeared = new TH1D(("h_Bias_Smeared2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        TH1D* h_variance_smeared = new TH1D(("h_Variance_Smeared2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        GetBiasVariance(h_smeared,h_bias_smeared,h_variance_smeared); 

        h_Bias_Ratio.push_back(static_cast<TH1D*>(h_bias_perfect->Clone("h_bias_ratio_perfect")));
        h_Bias_Ratio.back()->Add(h_bias_orig,-1);
        h_Bias_Ratio.back()->SetLineColor(colors.at(i_e));
        h_Bias_Ratio.back()->SetLineStyle(2);
        h_Bias_Ratio.back()->SetLineWidth(2);
        hs_Bias_Ratio->Add(h_Bias_Ratio.back()); 
        l->AddEntry(h_Bias_Ratio.back(),(est+" Perfect").c_str(),"L");

        h_Bias_Ratio.push_back(static_cast<TH1D*>(h_bias_smeared->Clone("h_bias_ratio_smeared")));
        h_Bias_Ratio.back()->Add(h_bias_orig,-1);
        h_Bias_Ratio.back()->SetLineColor(colors.at(i_e));
        h_Bias_Ratio.back()->SetLineStyle(3);
        h_Bias_Ratio.back()->SetLineWidth(2);
        hs_Bias_Ratio->Add(h_Bias_Ratio.back()); 
        l->AddEntry(h_Bias_Ratio.back(),(est+" Smeared").c_str(),"L");

        h_Variance_Ratio.push_back(static_cast<TH1D*>(h_variance_perfect->Clone("h_variance_ratio_perfect")));
        h_Variance_Ratio.back()->Add(h_variance_orig,-1);
        h_Variance_Ratio.back()->SetLineColor(colors.at(i_e));
        h_Variance_Ratio.back()->SetLineStyle(2);
        h_Variance_Ratio.back()->SetLineWidth(2);
        hs_Variance_Ratio->Add(h_Variance_Ratio.back()); 

        h_Variance_Ratio.push_back(static_cast<TH1D*>(h_variance_smeared->Clone("h_variance_ratio_smeared")));
        h_Variance_Ratio.back()->Add(h_variance_orig,-1);
        h_Variance_Ratio.back()->SetLineColor(colors.at(i_e));
        h_Variance_Ratio.back()->SetLineStyle(3);
        h_Variance_Ratio.back()->SetLineWidth(2);
        hs_Variance_Ratio->Add(h_Variance_Ratio.back()); 

        delete h_orig;
        delete h_bias_orig;
        delete h_variance_orig;
        delete h_perfect;
        delete h_bias_perfect;
        delete h_variance_perfect;
        delete h_smeared;
        delete h_bias_smeared;
        delete h_variance_smeared;

      }

      hs_Bias_Ratio->Draw("nostack HIST");
      l->Draw();
      c->Print(("Plots/Bias_Ratio_Energy_"+gen+".pdf").c_str());  
      c->Clear();

      hs_Variance_Ratio->Draw("nostack HIST");
      l->Draw();
      c->Print(("Plots/Variance_Ratio_Energy_"+gen+".pdf").c_str());  
      c->Clear();

      l->Clear();

    }

    // Ratio of bias of other generators to GENIE 
    for(size_t i_f=1;i_f<InputFiles_v.size();i_f++){

      std::string gen = Generators_v.at(i_f);

      THStack* hs_Bias_Ratio  = new THStack("hs",";True Neutrino Energy (GeV);Frac. Bias Generator - Frac. Bias GENIE");
      THStack* hs_Variance_Ratio  = new THStack("hs",";True Neutrino Energy (GeV);Frac. Bias Generator - Frac. Bias GENIE");
      std::vector<TH1D*> h_Bias_Ratio;
      std::vector<TH1D*> h_Variance_Ratio;

      for(size_t i_e=0;i_e<kMAX;i_e++){

        if(i_e == kMuonKin || i_e == kSFMethod) continue;
        std::string est = estimators_str.at(i_e);

        TH3D* h_genie = h_TrueEnergy_RecoEnergy_Coarse.at(i_e).at(0);

        h_genie->GetZaxis()->SetRange(1,1); 
        TH2D* h_genie_orig =  static_cast<TH2D*>(h_genie->Project3D("yx")->Clone("h_genie_orig2")); 
        TH1D* h_genie_bias_orig = new TH1D(("h_genie_Bias_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        TH1D* h_genie_variance_orig = new TH1D(("h_genie_Variance_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        GetBiasVariance(h_genie_orig,h_genie_bias_orig,h_genie_variance_orig); 

        h_genie->GetZaxis()->SetRange(2,2); 
        TH2D* h_genie_perfect =  static_cast<TH2D*>(h_genie->Project3D("yx")->Clone("h_genie_perfect2")); 
        TH1D* h_genie_bias_perfect = new TH1D(("h_genie_Bias_Perfect2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        TH1D* h_genie_variance_perfect = new TH1D(("h_genie_Variance_Perfect2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        GetBiasVariance(h_genie_perfect,h_genie_bias_perfect,h_genie_variance_perfect); 

        h_genie->GetZaxis()->SetRange(3,3); 
        TH2D* h_genie_smeared =  static_cast<TH2D*>(h_genie->Project3D("yx")->Clone("h_genie_smeared2")); 
        TH1D* h_genie_bias_smeared = new TH1D(("h_genie_Bias_Smeared2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        TH1D* h_genie_variance_smeared = new TH1D(("h_genie_Variance_Smeared2_"+gen+"_"+est).c_str(),"",nbins,binning_a);
        GetBiasVariance(h_genie_smeared,h_genie_bias_smeared,h_genie_variance_smeared); 

        TH3D* h_notgenie = h_TrueEnergy_RecoEnergy_Coarse.at(i_e).at(i_f);

        h_notgenie->GetZaxis()->SetRange(1,1); 
        TH2D* h_notgenie_orig =  static_cast<TH2D*>(h_notgenie->Project3D("yx")->Clone("h_notgenie_orig2")); 
        h_Bias_Ratio.push_back(new TH1D(("h_notgenie_Bias_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        h_Variance_Ratio.push_back(new TH1D(("h_notgenie_Variance_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        GetBiasVariance(h_notgenie_orig,h_Bias_Ratio.back(),h_Variance_Ratio.back()); 

        h_Bias_Ratio.back()->Add(h_genie_bias_orig,-1); 
        h_Variance_Ratio.back()->Add(h_genie_variance_orig,-1); 

        h_Bias_Ratio.back()->SetLineColor(colors.at(i_e)); 
        h_Bias_Ratio.back()->SetLineWidth(2); 
        h_Variance_Ratio.back()->SetLineColor(colors.at(i_e)); 
        h_Variance_Ratio.back()->SetLineWidth(2); 

        hs_Bias_Ratio->Add(h_Bias_Ratio.back());
        hs_Variance_Ratio->Add(h_Variance_Ratio.back());

        l->AddEntry(h_Bias_Ratio.back(),(est+" Original").c_str(),"L");


        h_notgenie->GetZaxis()->SetRange(2,2); 
        TH2D* h_notgenie_perfect =  static_cast<TH2D*>(h_notgenie->Project3D("yx")->Clone("h_notgenie_perfect2")); 
        h_Bias_Ratio.push_back(new TH1D(("h_notgenie_Bias_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        h_Variance_Ratio.push_back(new TH1D(("h_notgenie_Variance_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        GetBiasVariance(h_notgenie_perfect,h_Bias_Ratio.back(),h_Variance_Ratio.back()); 

        h_Bias_Ratio.back()->Add(h_genie_bias_perfect,-1); 
        h_Variance_Ratio.back()->Add(h_genie_variance_perfect,-1); 

        h_Bias_Ratio.back()->SetLineColor(colors.at(i_e)); 
        h_Bias_Ratio.back()->SetLineWidth(2); 
        h_Bias_Ratio.back()->SetLineStyle(2); 
        h_Variance_Ratio.back()->SetLineColor(colors.at(i_e)); 
        h_Variance_Ratio.back()->SetLineStyle(2); 
        h_Variance_Ratio.back()->SetLineWidth(2); 

        hs_Bias_Ratio->Add(h_Bias_Ratio.back());
        hs_Variance_Ratio->Add(h_Variance_Ratio.back());

        l->AddEntry(h_Bias_Ratio.back(),(est+" Perfect").c_str(),"L");

        h_notgenie->GetZaxis()->SetRange(3,3); 
        TH2D* h_notgenie_smeared =  static_cast<TH2D*>(h_notgenie->Project3D("yx")->Clone("h_notgenie_smeared2")); 
        h_Bias_Ratio.push_back(new TH1D(("h_notgenie_Bias_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        h_Variance_Ratio.push_back(new TH1D(("h_notgenie_Variance_Orig2_"+gen+"_"+est).c_str(),"",nbins,binning_a));
        GetBiasVariance(h_notgenie_smeared,h_Bias_Ratio.back(),h_Variance_Ratio.back()); 

        h_Bias_Ratio.back()->Add(h_genie_bias_smeared,-1); 
        h_Variance_Ratio.back()->Add(h_genie_variance_smeared,-1); 

        h_Bias_Ratio.back()->SetLineColor(colors.at(i_e)); 
        h_Bias_Ratio.back()->SetLineWidth(2); 
        h_Bias_Ratio.back()->SetLineStyle(3); 
        h_Variance_Ratio.back()->SetLineColor(colors.at(i_e)); 
        h_Variance_Ratio.back()->SetLineWidth(2); 
        h_Variance_Ratio.back()->SetLineStyle(3); 

        hs_Bias_Ratio->Add(h_Bias_Ratio.back());
        hs_Variance_Ratio->Add(h_Variance_Ratio.back());

        l->AddEntry(h_Bias_Ratio.back(),(est+" Smeared").c_str(),"L");

      }

      hs_Bias_Ratio->Draw("nostack HIST");
      l->Draw();
      c->Print(("Plots/Bias_Generator_Ratio_Energy_"+gen+".pdf").c_str());  
      c->Clear();

      hs_Variance_Ratio->Draw("nostack HIST");
      l->Draw();
      c->Print(("Plots/Variance_Generator_Ratio_Energy_"+gen+".pdf").c_str());  
      c->Clear();

      l->Clear();

    }


  }
