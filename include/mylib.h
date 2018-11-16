#ifndef mylib_h
#define mylib_h

map<TString, TH2D*> map_chi2; 
map<TString, TH2D*> map_mean_diff;
map<TString, TH2D*> map_Kolmogorov_Smirnov;
map<TString, TH1D*> map_Nvtx_hist;
map<TString, TH1D*> map_expected_hist;


vector<double> Get_reco_eff_and_fake_prob(TH1D* hist, int N_true_vertex){
  TString str_N_true_vertex = TString::Itoa(N_true_vertex, 10);
  

  vector<double> vertex_reco_eff_and_fake_prob;
  double vertex_reco_eff = 1.;
  double fake_prob = 0.;
  
  int N_reco_eff = 100;
  double reco_eff_interval = 1. / ( N_reco_eff + 0. );
  int N_fake_lambda = 100;
  double fake_prob_interval = 1.;
  
  TH1D *current_hist = (TH1D*)hist -> Clone("currrent_hist");
  current_hist -> Scale(1. / current_hist -> Integral());
  
  map_chi2["chi2_map" + str_N_true_vertex] = new TH2D("Nvtx_" + str_N_true_vertex + "_chi2map", "Nvtx_" + str_N_true_vertex + "_chi2map", N_reco_eff, 0., 1., N_fake_lambda, 0., (N_fake_lambda + 0.) * fake_prob_interval);
  map_mean_diff["mean_diff_map" + str_N_true_vertex] = new TH2D("Nvtx_" + str_N_true_vertex + "_meanDiffmap", "Nvtx_" + str_N_true_vertex + "_meanDiffmap", N_reco_eff, 0., 1., N_fake_lambda, 0., (N_fake_lambda + 0.) * fake_prob_interval); 
  map_Kolmogorov_Smirnov["Kolmogorov_Smirnov_map" + str_N_true_vertex] = new TH2D("Nvtx_" + str_N_true_vertex + "_Kolmogorov_Smirnovmap", "Nvtx_" + str_N_true_vertex + "_Kolmogorov_Smirnovmap", N_reco_eff, 0., 1., N_fake_lambda, 0., (N_fake_lambda + 0.) * fake_prob_interval);
  double min_chi_sqaure = 99999.;
  for(int i_prob = 0; i_prob < N_reco_eff; i_prob++){
    for(int i_fake = 0; i_fake < N_fake_lambda; i_fake++){
      double current_p = 0. + reco_eff_interval * (i_prob + 0.); // current vtx reco eff
      double current_lambda = fake_prob_interval + fake_prob_interval * (i_fake + 0.); // current lambda for Poisson distrubution for fake Nvtx
      
      //if(i_prob % 100 == 0 || i_fake % 100 == 0) cout << "i_prob : " << i_prob << ", i_fake : " << i_fake << endl;
      
      double current_chi2_sum = 0.;
      
      //////////////////////////////////////////////////////////////////////////////
      // Loop to get chi2 sum of current p (vtx reco eff) and q (fake prob)
      //////////////////////////////////////////////////////////////////////////////
      
      TH1D * h_for_chi2 = new TH1D("", "", 100, 0., 100.);
      
      for(int i_vtx = 0; i_vtx < 100; i_vtx++){
	
	double current_y = current_hist -> GetBinContent(i_vtx + 1);
	double expected_y = 0.;
	for(int i_k = 0; i_k < i_vtx; i_k++){ // i_k is number of reconstructed vertex among true vertices
	  if(i_k > N_true_vertex) break;
	  double expected_reco_nvtx = pow(current_p, i_k) * pow(1 - current_p, N_true_vertex - i_k) * TMath::Factorial(N_true_vertex) / (TMath::Factorial(N_true_vertex - i_k) * TMath::Factorial(i_k) );
	  double expected_fake_nvtx = TMath::Exp(-1. * current_lambda) * pow(current_lambda, i_vtx - i_k) / TMath::Factorial(i_vtx - i_k);
	  //if(i_vtx - i_k == 0) expected_fake_nvtx = 1. - current_q;
	  
	  expected_y = expected_y + expected_reco_nvtx * expected_fake_nvtx;
	}// end for i_k
	h_for_chi2 -> SetBinContent(i_vtx + 1, expected_y);
	
      }// end for i_vtx
      
      h_for_chi2 -> Scale(1. / h_for_chi2 -> Integral());
      for(int i_vtx = 0; i_vtx < 100; i_vtx++){
	double current_chi2 = 0.;
	if(current_hist -> GetBinContent(i_vtx + 1) + h_for_chi2 -> GetBinContent(i_vtx + 1) > 0)
	  current_chi2 = pow(current_hist -> GetBinContent(i_vtx + 1) - h_for_chi2 -> GetBinContent(i_vtx + 1), 2) / pow(current_hist -> GetBinContent(i_vtx + 1) + h_for_chi2 -> GetBinContent(i_vtx + 1), 2);
	current_chi2_sum = current_chi2_sum + current_chi2;
      }
      current_chi2_sum = current_chi2_sum / 100.;
      //cout << "current_chi2_sum : " << current_chi2_sum << endl;
      
      map_chi2["chi2_map" + str_N_true_vertex] -> Fill(current_p + 0.001, current_lambda, current_chi2_sum);
      map_mean_diff["mean_diff_map" + str_N_true_vertex] -> Fill(current_p + 0.001, current_lambda, fabs(current_hist -> GetMean() - h_for_chi2 -> GetMean()) );
      Double_t ks = current_hist->KolmogorovTest(h_for_chi2, "M");
      map_Kolmogorov_Smirnov["Kolmogorov_Smirnov_map" + str_N_true_vertex] -> Fill(current_p + 0.001, current_lambda, ks);
      //cout << "ks : " << ks << endl;
      if(current_chi2_sum < min_chi_sqaure){
	vertex_reco_eff = current_p;
	fake_prob = current_lambda;
	min_chi_sqaure = current_chi2_sum;
      }
      
    }// end for i_fake
  }// end for i_prob
  cout << "[Get_reco_eff_and_fake_prob] min_chi_sqaure : " << min_chi_sqaure << endl;
  
  vertex_reco_eff_and_fake_prob = {vertex_reco_eff, fake_prob,};
  
  return vertex_reco_eff_and_fake_prob;

}

TH1D *MakeDistribution_with_eff(double p, double lambda, int N){
  
  TH1D *out = new TH1D("expected_Nvtx", "expected_Nvtx", 100, 0., 100);
  
  for(int i_vtx = 0; i_vtx < 100; i_vtx++){

    double expected_y = 0.;
    for(int i_k = 0; i_k < i_vtx; i_k++){ // i_k is number of reconstructed vertex among true vertices                                                                                                                                                                      
      if(i_k > N) break;
      double expected_reco_nvtx = pow(p, i_k) * pow(1 - p, N - i_k) * TMath::Factorial(N) / (TMath::Factorial(N - i_k) * TMath::Factorial(i_k) );
      double expected_fake_nvtx = TMath::Exp(-1. * lambda) * pow(lambda, i_vtx - i_k) / TMath::Factorial(i_vtx - i_k);
      
      expected_y = expected_y + expected_reco_nvtx * expected_fake_nvtx;
    }// end for i_k                                                                                                                                                                                                                                                         
    out -> SetBinContent(i_vtx + 1, expected_y);
  }// end for i_vtx    
  
  out -> Scale(1. / out->Integral() );
  
  return out;

}


#endif
