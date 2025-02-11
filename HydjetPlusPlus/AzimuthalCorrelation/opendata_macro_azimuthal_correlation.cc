///*** To run it, please, do
///*** 1) compile it by doing: << make opendata_macro_azimuthal_correlation >>
///*** 2) run the created executable: << ./opendata_macro_azimuthal_correlation  >>

#include<cstdio>
#include<cstdlib>
#include<iostream>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "Math/Vector4D.h"

///Define functions globally

void func_setAndSaveCanvas_(TCanvas *c1, TString figName);

void func_setAndDrawHistograms_(TH1D *h, TString titleX, TString titleY);

double deltaEta(const TVector3 a, const TVector3 b);

double deltaPhi(const TVector3 a, const TVector3 b);


///main function: access tree and loop over event and particles in the event
int main(){

std::cout<<"Starting..."<<std::endl;

const int particle_type = 0; // 0 for all (pions, kaons, protons); 1 for pions; 2 for kaons; 3 for protons

///access Hydjet++ generated events in ROOT  TTree format
std::unique_ptr<TFile> f( TFile::Open("/home/pedrolunardi/IC/Events/HiForestAOD_DATA_FlowAnalysis_630k.root") );
//std::unique_ptr<TFile> f( TFile::Open("/Users/cesarbernardes/Dropbox/Ubuntu_1204/AltasEnergias/ProfessorUFRGS/OrientacaoDeAlunos/Tools/CMSOpenData/cms_open_data_work_hin10/CMSSW_3_9_2_patch5/src/HiForest/HiForestProducer/HiForestAOD_DATA_FlowAnalysis_200k.root") );

const int NEventsArraySize = 30508; //IMPORTANT: total number of events generated

///some global variables -- after can think in add it in a ".h" file
bool doSymmetrisation_=true; //symmetrise the correlation function in each quadrant
bool doBackground_=true; //if "true" will do event mixing correlation for background - IMPORTANT: it will take a lot of time to finish

///File to store histograms - after we can use a variable to set centrality bins - for now I am just using one single bin
TFile *output = new TFile(Form("corr_cent_%d_%d.root",0,5), "recreate");

TVector3 particle;
TVector3 trigger;
TVector3 associated;

const int ptTbins_=1, ptNbins_=1; //number of pT ranges for Trigger & Associated particles 
std::vector<TVector3> allevents[NEventsArraySize][ptTbins_+ptNbins_];
//Max and Min boundaries for Trigger & Associated particles range
double ptTrigMax_[ptTbins_];
double ptTrigMin_[ptTbins_];
double ptAssMax_[ptNbins_];
double ptAssMin_[ptNbins_];
ptTrigMax_[0]=3.5; //in GeV/c
ptTrigMin_[0]=3.0;
ptAssMax_[0]=1.5;
ptAssMin_[0]=1.0;


//Get trees
auto input_tree = f->Get<TTree>("demo/Flow");

// Disable everything...
input_tree->SetBranchStatus("*", false);

// ...but the branch we need
input_tree->SetBranchStatus("Ntrk", true);
input_tree->SetBranchStatus("HFsumET", true);
input_tree->SetBranchStatus("trkPt", true);
input_tree->SetBranchStatus("trkEta", true);
input_tree->SetBranchStatus("trkPhi", true);

int Ntrk;
float HFsumET;

float trkPt[2000];
float trkPhi[2000];
float trkEta[2000];

input_tree->SetBranchAddress("Ntrk",&Ntrk);
input_tree->SetBranchAddress("HFsumET",&HFsumET);
input_tree->SetBranchAddress("trkPt",&trkPt);
input_tree->SetBranchAddress("trkEta",&trkEta);
input_tree->SetBranchAddress("trkPhi",&trkPhi);

auto hist_title_label = "";
if      (particle_type == 0){hist_title_label = "Dados CMS Open Data";}
else if (particle_type == 1){hist_title_label = "Pions";}
else if (particle_type == 2){hist_title_label = "Kaons";}
else if (particle_type == 3){hist_title_label = "Protons";}
else{}

///Histograms

// 1D - From particle four-vector
auto h_trkPt = new TH1D("pt",hist_title_label,200,0,20);
auto h_trkEta = new TH1D("eta",hist_title_label,200,-2.5,2.5);
auto h_trkPhi = new TH1D("phi",hist_title_label,100,-3.5,3.5);

// Histograms for correlation
auto h_trkPt_trigger_all = new TH1D("h_trkPt_trigger_all","All Trigger Particles",200,0,20);

TH2D* signal[ptTbins_][ptNbins_];
TH2D* background[ptTbins_][ptNbins_];
for (int j = 0; j < ptTbins_; j++) {
   for (int i = 0; i < ptNbins_; i++) {
    signal[j][i] = new TH2D(Form("signal_trig_%d_%d_ass_%d_%d",int(ptTrigMin_[j]*10), int(ptTrigMax_[j]*10), int(ptAssMin_[i]*10), int(ptAssMax_[i]*10) ) , ";#Delta#eta;#Delta#phi" , 33, -4.8 - 4.8/32.0, 4.8 + 4.8/32.0,  31, -0.5*TMath::Pi()+TMath::Pi()/32, 1.5*TMath::Pi()-TMath::Pi()/32);
    background[j][i] = new TH2D(Form("background_trig_%d_%d_ass_%d_%d",int(ptTrigMin_[j]*10), int(ptTrigMax_[j]*10), int(ptAssMin_[i]*10), int(ptAssMax_[i]*10) ) , ";#Delta#eta;#Delta#phi" , 33, -4.8 - 4.8/32.0, 4.8 + 4.8/32.0,  31, -0.5*TMath::Pi()+TMath::Pi()/32, 1.5*TMath::Pi()-TMath::Pi()/32);
   }
}

//std::cout << "HFsumET: " << HFsumET << std::endl;

for (int i = 0; input_tree->LoadTree(i) >= 0; ++i) {

   // Load the data for the given tree entry (i.e. a collision event)
   input_tree->GetEntry(i);

   //if (300 < HFsumET && HFsumET < 1210) continue; //30 - 60%

   //if (1210 < HFsumET && HFsumET < 2775) continue; // 10 - 30%
						   
   if (3390 < HFsumET && HFsumET < 4500) continue; // 0 - 5%

   // Loop over all particles in the event
   for (int j = 0; j < Ntrk; j++) {

      if (!(trkPt[j] > 0.3 && fabs(trkEta[j]) < 2.4 && fabs(trkPhi[j]) < 3.2)) continue;

      h_trkPt->Fill(trkPt[j]);
      h_trkEta->Fill(trkEta[j]);
      h_trkPhi->Fill(trkPhi[j]);

      particle.SetPtEtaPhi(trkPt[j],trkEta[j],trkPhi[j]);

      //store Associated particles information
      for (int ii = 0; ii < ptNbins_; ii++)
        if (trkPt[j] > ptAssMin_[ii] && trkPt[j] < ptAssMax_[ii] ) (allevents[i][ptTbins_ + ii]).push_back(particle);

      //store Trigger particles information  
      for (int jj = 0; jj < ptTbins_; jj++)
         if ( trkPt[j] > ptTrigMin_[jj] && trkPt[j] < ptTrigMax_[jj] ) (allevents[i][jj]).push_back(particle);

   }//end particle loop

   int Ntrig=0, Nass=0; //number of trigger & associated particles in each event
   double dEta=0.0, dPhi=0.0; //to build the correlation
   for (int i_trigbin = 0; i_trigbin < ptTbins_; i_trigbin++) {
      Ntrig = (allevents[i][i_trigbin]).size(); 
      for (int w = 0; w < Ntrig; w++) {
	 trigger = (allevents[i][i_trigbin])[w]; 
         h_trkPt_trigger_all->Fill(trigger.Pt());
	 for (int i_assbin = 0; i_assbin < ptNbins_; i_assbin++) {
	    Nass = (allevents[i][ptTbins_+i_assbin]).size();
	    for (int z = 0; z < Nass; z++) {
	       associated = (allevents[i][ptTbins_+i_assbin])[z];
	       dEta=deltaEta(trigger,associated); 
	       dPhi=deltaPhi(trigger,associated);
	       //std::cout<<"dEta : "<<dEta<<";  dPhi : "<<dPhi<<std::endl;
	       if (dEta != 0 && dPhi != 0) signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1./Ntrig); //See Eq. 2 of https://arxiv.org/pdf/1201.3158.pdf
	       if (dEta != 0  && dPhi != 0 && doSymmetrisation_){
                  dEta=deltaEta(trigger,associated); dPhi=deltaPhi(associated,trigger);
                  signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1./Ntrig); //See Eq. 2 of https://arxiv.org/pdf/1201.3158.pdf
                  dEta=deltaEta(associated,trigger); dPhi=deltaPhi(associated,trigger);
                  signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1./Ntrig); //See Eq. 2 of https://arxiv.org/pdf/1201.3158.pdf
                  dEta=deltaEta(associated,trigger); dPhi=deltaPhi(trigger,associated);
                  signal[i_trigbin][i_assbin]->Fill(dEta,dPhi,1./Ntrig); //See Eq. 2 of https://arxiv.org/pdf/1201.3158.pdf
	       }
            }		    
	 }	 
      }	      
   }
}

///Background part
if(doBackground_){

   int Ntrig_bg, Nass_bg;
   int NevtMix;
   double dEta_bg, dPhi_bg;

   NevtMix=1;
   if ( NEventsArraySize<11 ) NevtMix = NEventsArraySize-1;   

   for (int i_ev=0; i_ev<NEventsArraySize; i_ev++) {
      int i_mix =0;	   
      for (int j_ev=i_ev+1; j_ev<NEventsArraySize; j_ev++){
        if(i_mix>=NevtMix)break;

        for (int i_trigbin=0; i_trigbin<ptTbins_; i_trigbin++) {
           Ntrig_bg = (allevents[i_ev][i_trigbin]).size();
	   for (int i_assbin=0; i_assbin<ptNbins_; i_assbin++) {
              Nass_bg = (allevents[j_ev][ptTbins_+i_assbin]).size();
	      for (int j=0; j<Ntrig_bg; j++) {
		 trigger = (allevents[i_ev][i_trigbin])[j];
		 for (int i=0; i<Nass_bg; i++) {
		    associated = (allevents[j_ev][ptTbins_+i_assbin])[i];
		    dEta_bg=deltaEta(trigger,associated); dPhi_bg=deltaPhi(trigger,associated);
		    if (dEta_bg != 0  && dPhi_bg != 0) background[i_trigbin][i_assbin]->Fill(dEta_bg,dPhi_bg,1./Ntrig_bg); //See Eq. 3 of https://arxiv.org/pdf/1201.3158.pdf
		    if (dEta_bg != 0  && dPhi_bg != 0 && doSymmetrisation_){
                       dEta_bg=deltaEta(trigger,associated); dPhi_bg=deltaPhi(associated,trigger);
		       background[i_trigbin][i_assbin]->Fill(dEta_bg,dPhi_bg,1./Ntrig_bg); //See Eq. 3 of https://arxiv.org/pdf/1201.3158.pdf
                       dEta_bg=deltaEta(associated,trigger); dPhi_bg=deltaPhi(associated,trigger);
		       background[i_trigbin][i_assbin]->Fill(dEta_bg,dPhi_bg,1./Ntrig_bg); //See Eq. 3 of https://arxiv.org/pdf/1201.3158.pdf
                       dEta_bg=deltaEta(associated,trigger); dPhi_bg=deltaPhi(trigger,associated);
		       background[i_trigbin][i_assbin]->Fill(dEta_bg,dPhi_bg,1./Ntrig_bg); //See Eq. 3 of https://arxiv.org/pdf/1201.3158.pdf
                    }
		 }	 
	      }
           }		   
	}	
      	i_mix++;      
      }	      
   }	   
}

///Plot histograms

auto c_pt = new TCanvas("c_pt", "pT distribution", 500, 500);
c_pt->cd(0);
func_setAndDrawHistograms_(h_trkPt,"p_{T} [GeV/c]","Number of entries per Bin");
gPad->SetLogy();
func_setAndSaveCanvas_(c_pt,"h_pt");

auto c_eta = new TCanvas("c_eta", "eta distribution", 500, 500);
c_eta->cd(0);
func_setAndDrawHistograms_(h_trkEta,"#eta","Number of entries per Bin");
gPad->SetLogy(0);
func_setAndSaveCanvas_(c_eta,"h_eta");

auto c_phi = new TCanvas("c_phi", "phi distribution", 500, 500);
c_phi->cd(0);
func_setAndDrawHistograms_(h_trkPhi,"#phi","Number of entries per Bin");
gPad->SetLogy(0);
func_setAndSaveCanvas_(c_phi,"h_phi");

//For correlations
auto c_pt_trigger_all = new TCanvas("c_pt_trigger_all", "pT distribution", 500, 500);
c_pt_trigger_all->cd(0);
func_setAndDrawHistograms_(h_trkPt_trigger_all,"p_{T} [GeV/c]","Number of entries per Bin");
gPad->SetLogy();
func_setAndSaveCanvas_(c_pt_trigger_all,"h_pt_trigger_all");

TCanvas *canvas[ptTbins_][ptNbins_];
for (int j = 0; j < ptTbins_; j++) {
   for (int i = 0; i < ptNbins_; i++) {
    canvas[j][i] = new TCanvas(Form("c_signal_trig_%d_%d_ass_%d_%d",int(ptTrigMin_[j]*10), int(ptTrigMax_[j]*10), int(ptAssMin_[i]*10), int(ptAssMax_[i]*10) ), "eta vs phi correlation", 500, 500);
    canvas[j][i]->cd(0);
    signal[j][i]->Draw("COLZ");
    canvas[j][i]->SaveAs(Form("Figures/deta_dphi_signal_trig_%d_%d_ass_%d_%d.pdf",int(ptTrigMin_[j]*10), int(ptTrigMax_[j]*10), int(ptAssMin_[i]*10), int(ptAssMax_[i]*10)));
   }
}

TCanvas *canvas_bg[ptTbins_][ptNbins_];
for (int j = 0; j < ptTbins_; j++) {
   for (int i = 0; i < ptNbins_; i++) {
    canvas_bg[j][i] = new TCanvas(Form("c_background_trig_%d_%d_ass_%d_%d",int(ptTrigMin_[j]*10), int(ptTrigMax_[j]*10), int(ptAssMin_[i]*10), int(ptAssMax_[i]*10) ), "eta vs phi correlation", 500, 500);
    canvas_bg[j][i]->cd(0);
    background[j][i]->Draw("COLZ");
    canvas_bg[j][i]->SaveAs(Form("Figures/deta_dphi_background_trig_%d_%d_ass_%d_%d.pdf",int(ptTrigMin_[j]*10), int(ptTrigMax_[j]*10), int(ptAssMin_[i]*10), int(ptAssMax_[i]*10)));
   }
}

///Write the historgrams to the "output" file 
h_trkPt->Write();
h_trkEta->Write();
h_trkPhi->Write();
for (int j = 0; j < ptTbins_; j++) {
   for (int i = 0; i < ptNbins_; i++) {
      signal[j][i]->Write();
      background[j][i]->Write();
   }
}


std::cout<<"Done..."<<std::endl;

return 0;

}


///Functions: just for nomenclature always name the function using "func_name_"


//function to set canvas
void func_setAndSaveCanvas_(TCanvas *c1, TString figName){

   c1->Range(0.04000706,-0.4080692,0.7405774,1.786744);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);  
   c1->SetLeftMargin(0.1650854);
   c1->SetRightMargin(0.0341556);
   c1->SetTopMargin(0.08508404);
   c1->SetBottomMargin(0.1859244);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SaveAs("Figures/"+figName+".pdf");

}


//function to set histograms
void func_setAndDrawHistograms_(TH1D *h, TString titleX, TString titleY){

   h->Draw();	
   h->GetXaxis()->SetTitle(titleX);
   h->GetYaxis()->SetTitle(titleY);
   h->SetLineColor(2);
   h->SetLineWidth(2);

}

//function to calculate delta eta between two particles
double deltaEta(const TVector3 a, const TVector3 b){

return a.Eta() - b.Eta();

}

//function to calculate delta phi between two particles
double deltaPhi(const TVector3 a, const TVector3 b){

  double Dphi;
  Dphi = a.Phi() - b.Phi();
  while (Dphi < -0.5*TMath::Pi()) Dphi = Dphi + 2.0*TMath::Pi();
  while (Dphi > 1.5*TMath::Pi()) Dphi = Dphi - 2.0*TMath::Pi();
  return Dphi;

}



