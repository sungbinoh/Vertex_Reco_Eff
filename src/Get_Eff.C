#include "canvas_margin.h"
#include "mylib.h"

void Get_Eff(int xxx=0){

  setTDRStyle();

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2(true);
  TH1::AddDirectory(kFALSE);

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString ENV_PLOT_PATH = getenv("PLOT_PATH");
  
  TString year = "2017";
  TString SampleName = "NvtxProfile_DYJets.root";
  
  TString base_filepath = WORKING_DIR+"/rootfiles/"+year+"/"+SampleName;
  TString base_plotpath = ENV_PLOT_PATH+"/"+year+"/";
  
  TFile *in = new TFile((base_filepath));
  
  vector<double> vector_vtx_eff;
  vector<double> vector_N_fake;
  vector_vtx_eff.clear();
  vector_N_fake.clear();

  for(int i = 0; i < 100; i++){
    cout << "Running for Ntrue vtx : " << i << endl;
    TString str_i = TString::Itoa(i, 10);
    map_Nvtx_hist["Nvtx_" + str_i] = (TH1D*)gDirectory -> Get("Nvtx_N_true_vtx_" + str_i);
    
    vector<double> eff_and_fake = Get_reco_eff_and_fake_prob(map_Nvtx_hist["Nvtx_" + str_i], i);
    cout << "eff_and_fake.at(0) : " << eff_and_fake.at(0) << ", eff_and_fake.at(1) : " << eff_and_fake.at(1) << endl;
    vector_vtx_eff.push_back(eff_and_fake.at(0));
    vector_N_fake.push_back(eff_and_fake.at(1));
    
    map_expected_hist["Nvtx_" + str_i + "_expected"] = MakeDistribution_with_eff(eff_and_fake.at(0), eff_and_fake.at(1), i);
    //map_expected_hist["Nvtx_" + str_i + "_expected"] = MakeDistribution_with_eff(0.7, 30, i);

  }
  
  //double integral_expected = expected_Nvtx -> Integral();
  //cout << "integral_expected : " << integral_expected << endl;
  
  if( !gSystem->mkdir(base_plotpath, kTRUE) ){
    cout
      << "###################################################" << endl
      << "Directoy " << base_plotpath << " is created" << endl
      << "###################################################" << endl
      << endl;
  }
  TFile *out = new TFile(base_plotpath+"/Output_plots.root","RECREATE");
  for(std::map< TString, TH1D* >::iterator mapit = map_Nvtx_hist.begin(); mapit!=map_Nvtx_hist.end(); mapit++){
    mapit->second->Scale(1. / mapit->second->Integral());
    mapit->second->Write();
  }
  for(std::map< TString, TH1D* >::iterator mapit = map_expected_hist.begin(); mapit!=map_expected_hist.end(); mapit++){
    mapit->second->Write();
  }
  for(std::map< TString, TH2D* >::iterator mapit = map_chi2.begin(); mapit!=map_chi2.end(); mapit++){
    mapit->second->Write();
  }
  /*
  for(std::map< TString, TH2D* >::iterator mapit = map_mean_diff.begin(); mapit!=map_mean_diff.end(); mapit++){
    mapit->second->Write();
  }
  */
  for(std::map< TString, TH2D* >::iterator mapit = map_Kolmogorov_Smirnov.begin(); mapit!=map_Kolmogorov_Smirnov.end(); mapit++){
    TCanvas *c = new TCanvas("","",800,800);
    gStyle -> SetOptStat(0);
    mapit->second->Draw("colz");
    TString current_name = mapit->second->GetName();
    c -> SaveAs(base_plotpath + "/" + current_name +".pdf");
    mapit->second->Write();
  }
  
  in->Close();
  out->Close();
    

}
