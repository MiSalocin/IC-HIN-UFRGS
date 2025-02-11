void vn_extractor(){

  TFile * f = TFile::Open("/home/pedrolunardi/IC/collaboration/HydjetPlusPlus/AzimuthalCorrelation/corrfunc1D_cent_0_5.root");

  TH1D* h = (TH1D*) f->Get("signal_dphi_projection_0");

  // Cria a função de fit

  TF1 *fit = new TF1("fit", "[0]*(1 + 2*[1]*cos(x)+2*[2]*cos(2*x) + 2*[3]*cos(3*x))", -1.4, 4.6);

  fit->SetParameters(60, 0.1, 0.1, 0.1);

  fit->SetParLimits(1, 0.000001, 10);

  fit->SetParLimits(2, 0.000001, 10);

  fit->SetParLimits(3, 0.000001, 10);

  // Cria um canvas para desenhar o histograma e o fit

  TCanvas *c = new TCanvas("c", "c", 800, 600);

  c->cd();

  // Desenha o histograma no canvas

  h->Draw("p*Error");

  // Faz o fit no histograma do dphi_signal

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

  h->Fit(fit, "R Q");

  h->Fit(fit, "R Q");

  h->Fit(fit, "R M");

  // Desenha o fit no canvas

  fit->Draw("same");

  // Salva o canvas como um pdf

  c->Draw("SURF1");
  //c->SetTitle("30 - 60 %");
  //c->SetTitle("10 - 30 %");
  c->SetTitle("0 - 5 %");
  c->SaveAs("FourierFit_0_5_opendata.pdf");

}
