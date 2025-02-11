///*** To run it, please, do: <<root -l -b -q make_corrfunc.C>>


void make_corrfunc(){

int nbins=1;
float aptassmin[]={1.0};
float aptassmax[]={1.5};
float apttrigmin[]={3.0};
float apttrigmax[]={3.5};

Double_t c=0.;
Double_t c1,c2=1.;
Double_t c3=30508; //IMPORTANT (temporary): should be equal to the total number of events generated  

TFile* hfile = new TFile("corr_cent_0_5.root");

TFile *output = new TFile("corrfunc1D_cent_0_5.root","recreate");

int ptassmin, ptassmax, pttrigmin, pttrigmax;

TCanvas *canvas[nbins];
TString canvasNames[nbins];
canvasNames[0]=Form("signal_trig_%d_%d_ass_%d_%d",pttrigmin,pttrigmax,ptassmin,ptassmax);

TH2D *signal[nbins];
TH2D *background[nbins];
TH1D *dphi_signal[nbins];
TH1D *dphi_bg[nbins];

for (int i = 0; i < nbins; i++) {
   ptassmin = aptassmin[i]*10;
   ptassmax = aptassmax[i]*10;
   pttrigmin = apttrigmin[i]*10;
   pttrigmax = apttrigmax[i]*10;
   signal[i] = (TH2D*) hfile->Get(Form("signal_trig_%d_%d_ass_%d_%d",pttrigmin,pttrigmax,ptassmin,ptassmax));
   background[i] = (TH2D*) hfile->Get(Form("background_trig_%d_%d_ass_%d_%d",pttrigmin,pttrigmax,ptassmin,ptassmax));
}

for (int i =0; i < nbins; i++) {

   dphi_signal[i] = signal[i]->ProjectionY(Form("signal_dphi_projection_%d",i),4,10); //projection from -3.75 till -1.65 in deltaEta. Since was symmetrized, no need to get positive part
 
   dphi_bg[i] = background[i]->ProjectionY(Form("background_dphi_projection_%d",i),4,10); //projection from -3.75 till -1.65 in deltaEta. Since was symmetrized, no need to get positive part

     	
   //dphi_signal[i]->Sumw2();
   //dphi_bg[i]->Sumw2();
   dphi_signal[i]->Divide(dphi_bg[i]);
 
     	
   c=background[i]->GetBinContent(background[i]->FindBin(0,0));
   c1=background[i]->GetYaxis()->GetBinWidth(1);
   c2=background[i]->GetXaxis()->GetBinWidth(1);

   std::cout<<"c : " << c <<"; c1 : "<< c1 <<"; c2 : "<<c2<<std::endl;
   //std::cout<<"c : " << c<<std::endl;

   //dphi_signal[i]->Scale(c/(4.*c1*c2));
   dphi_signal[i]->Scale(c/(4.*c1*c2*c3));

   dphi_signal[i]->Write();		
   dphi_bg[i]->Write();			


   signal[i]->Divide(background[i]);
   //signal[i]->Scale(c/(4.*c1*c2));
   signal[i]->Scale(c/(4.*c1*c2*c3));
   signal[i]->Write();


   canvas[i] = new TCanvas("tc_"+canvasNames[i],"tc_"+canvasNames[i],500,500);
   canvas[i]->cd(1);
   canvas[i]->cd(1)->SetTopMargin(0.06);
   canvas[i]->cd(1)->SetLeftMargin(0.18);
   canvas[i]->cd(1)->SetRightMargin(0.043);
   canvas[i]->cd(1)->SetBottomMargin(0.145);
   signal[i]->Draw("SURF1");
   //signal[i]->SetTitle("30 - 60 %");
   //signal[i]->SetTitle("10 - 30 %");
   signal[i]->SetTitle("0 - 5 %");
   signal[i]->GetXaxis()->CenterTitle(true);
   signal[i]->GetYaxis()->CenterTitle(true);
   signal[i]->GetZaxis()->CenterTitle(true);
   signal[i]->GetXaxis()->SetTitleOffset(2.);
   signal[i]->GetYaxis()->SetTitleOffset(2.);
   signal[i]->GetZaxis()->SetTitleOffset(2.);
   signal[i]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}");
   gStyle->SetOptStat(0);
   canvas[i]->SaveAs(Form("fig_2Dcorrelation_signalOverBackg_cent_0_5_trig_%d_%d_ass_%d_%d.pdf",pttrigmin,pttrigmax,ptassmin,ptassmax));

}


}
