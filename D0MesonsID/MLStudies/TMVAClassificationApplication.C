#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

// RUN: root -l -b -q TMVAClassificationApplication.C\(\"BDT\"\)

void TMVAClassificationApplication( TString myMethodList = "" )
{
        
    //---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
    // Cut optimisation
    Use["Cuts"]            = 1;
    Use["CutsD"]           = 1;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    //
    // 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 1;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 1;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 1;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 1; // k-nearest neighbour method
    //
    // Linear Discriminant Analysis
    Use["LD"]              = 1; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // Function Discriminant analysis
    Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    //
    // Support Vector Machine
    Use["SVM"]             = 1;
    //
    // Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    //
    // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 1;
    // ---------------------------------------------------------------
    Use["Plugin"]          = 0;
    Use["Category"]        = 0;
    Use["SVM_Gauss"]       = 0;
    Use["SVM_Poly"]        = 0;
    Use["SVM_Lin"]         = 0;
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;

    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
        std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
        for (UInt_t i=0; i<mlist.size(); i++) {
            std::string regMethod(mlist[i]);
            if (Use.find(regMethod) == Use.end()) {
                std::cout << "Method \"" << regMethod
                        << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
                std::cout << it->first << " ";
                }
                std::cout << std::endl;
                return;
            }
            Use[regMethod] = 1;
        }
    }

    // --------------------------------------------------------------------------------------------------
    // Create the Reader object
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    Float_t D3DDecayLength, D3DDecayLengthSignificance, DTrk1Chi2n, DTrk1PtErr, DTrk2Chi2n, DVtxProb, D3DPointingAngle, DDca, DTtrk1Pt, DTrk2Pt, DTrk2PtErr, DTrk1Eta, DTrk2Eta, DxyDCASignificanceDaugther1, DxyDCASignificanceDaugther2, DzDCASignificanceDaugther1, DzDCASignificanceDaugther2, DMass, DPhi, DPt, DRapidity, DGen;
    reader->AddVariable("D3DDecayLength", &D3DDecayLength);
    reader->AddVariable("D3DDecayLengthSignificance", &D3DDecayLengthSignificance);
    reader->AddVariable("DTrk1Chi2n", &DTrk1Chi2n);
    reader->AddVariable("DTrk1PtErr", &DTrk1PtErr);
    reader->AddVariable("DTrk2Chi2n", &DTrk2Chi2n);
    reader->AddVariable("DVtxProb", &DVtxProb);
    reader->AddVariable("D3DPointingAngle", &D3DPointingAngle);
    reader->AddVariable("DDca", &DDca);
    reader->AddVariable("DTtrk1Pt", &DTtrk1Pt);
    reader->AddVariable("DTrk2Pt", &DTrk2Pt);
    reader->AddVariable("DTrk2PtErr", &DTrk2PtErr);
    reader->AddVariable("DTrk1Eta", &DTrk1Eta);
    reader->AddVariable("DTrk2Eta", &DTrk2Eta);
    reader->AddVariable("DxyDCASignificanceDaugther1", &DxyDCASignificanceDaugther1);
    reader->AddVariable("DxyDCASignificanceDaugther2", &DxyDCASignificanceDaugther2);
    reader->AddVariable("DzDCASignificanceDaugther1", &DzDCASignificanceDaugther1);
    reader->AddVariable("DzDCASignificanceDaugther2", &DzDCASignificanceDaugther2);
    reader->AddVariable("DMass", &DMass);
    reader->AddVariable("DPhi", &DPhi);
    reader->AddVariable("DPt", &DPt);
    reader->AddVariable("DRapidity", &DRapidity);
    reader->AddVariable("DGen", &DGen);

    // Spectator variables declared in the training have to be added to the reader, too
    Float_t DGenpt, DGenphi, DGeny, DGenIndex;
    reader->AddSpectator("DGenpt", &DGenpt);
    reader->AddSpectator("DGenphi", &DGenphi);
    reader->AddSpectator("DGeny", &DGeny);
    reader->AddSpectator("DGenIndex", &DGenIndex);

    // Book the MVA methods
    TString dir    = "/home/allanfgodoi/Desktop/";
    //TString dir    = "/Users/cesarbernardes/Dropbox/Ubuntu_1204/AltasEnergias/ProfessorUFRGS/OrientacaoDeAlunos/IC_TCC/ReposGit/IC-HIN-UFRGS/D0MesonsID/MLStudies/";
    TString prefix = "TMVAClassification";

    // Book method(s)
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
        if (it->second) {
            TString methodName = TString(it->first) + TString(" method");
            TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            reader->BookMVA( methodName, weightfile );
        }
    }

    // Book output histograms
    UInt_t nbin = 100;
    TH1F *histBdt_all(0);
    TH1F *histBdt_signal(0);
    TH1F *histBdt_background(0);

    TH1F *histBdt_raw_DMass(0);
    TH1F *histBdt_raw_DPhi(0);
    TH1F *histBdt_raw_DPt(0);
    TH1F *histBdt_raw_DRapidity(0);

    TH1F *histBdt_DMass_cut(0);

    TH1D* h_bdtSignal(0);
    TH1D* h_bdtBackground(0);

    if (Use["BDT"]) {
        histBdt_all = new TH1F("MVA_BDT_all", "MVA_BDT_all", 100, -0.45, 0.25 );
        histBdt_signal = new TH1F("MVA_BDT_signal", "MVA_BDT_signal", 100, -0.45, 0.25 );
        histBdt_background = new TH1F("MVA_BDT_background", "MVA_BDT_background", 100, -0.45, 0.25 );
        
        histBdt_raw_DMass = new TH1F("raw_DMass", "raw_DMass", 100, 1.70, 2.05);
        histBdt_raw_DPhi = new TH1F("raw_DPhi", "raw_DPhi", 100, -3.5, 3.5);
        histBdt_raw_DPt = new TH1F("raw_DPt", "raw_DPt", 100, 3.4, 4.3);
        histBdt_raw_DRapidity = new TH1F("raw_DRapidity", "raw_DRapidity", 100, -0.9, 0.9);

        histBdt_DMass_cut = new TH1F("DMass_(>0.0)", "DMass_(>0.0)", 100, 1.70, 2.05);

        h_bdtSignal = new TH1D("BDT_Signal", "BDT_Signal", 100, -0.5, 0.2);
        h_bdtBackground = new TH1D("BDT_Background", "BDT_Background", 100, -0.5, 0.2);
    }

    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //
    TFile *input(0);
    TString fname = "/home/allanfgodoi/Desktop/tree_skim_MC_promptTest.root";
    //TString fname = "/Users/cesarbernardes/Dropbox/Ubuntu_1204/AltasEnergias/ProfessorUFRGS/OrientacaoDeAlunos/IC_TCC/TopicosDeEstudo/D0_Selection/D0_MC_SkimmedTrees/tree_skim_MC_prompt.root";
    if (!gSystem->AccessPathName( fname ))
        input = TFile::Open( fname ); // check if file in local directory exists
    else
        input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
    if (!input) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }
    std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

    // Event loop
    // Prepare the event tree
    // - Here the variable names have to corresponds to your tree
    // - You can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
    std::cout << "--- Select signal+background sample" << std::endl;
    TTree* theTree = (TTree*)input->Get("Dfinder/ntDkpi");
    std::vector<float> * vec_D3DDecayLength = 0; 
    std::vector<float> * vec_D3DDecayLengthSignificance = 0;
    std::vector<float> * vec_DTrk1Chi2n = 0;
    std::vector<float> * vec_DTrk1PtErr = 0;
    std::vector<float> * vec_DTrk2Chi2n = 0;
    std::vector<float> * vec_DVtxProb = 0;
    std::vector<float> * vec_D3DPointingAngle = 0;
    std::vector<float> * vec_DDca = 0;
    std::vector<float> * vec_DTtrk1Pt = 0;
    std::vector<float> * vec_DTrk2Pt = 0;
    std::vector<float> * vec_DTrk2PtErr = 0;
    std::vector<float> * vec_DTrk1Eta = 0;
    std::vector<float> * vec_DTrk2Eta = 0;
    std::vector<float> * vec_DxyDCASignificanceDaugther1 = 0;
    std::vector<float> * vec_DxyDCASignificanceDaugther2 = 0;
    std::vector<float> * vec_DzDCASignificanceDaugther1 = 0;
    std::vector<float> * vec_DzDCASignificanceDaugther2 = 0;
    std::vector<float> * vec_DMass = 0;
    std::vector<float> * vec_DPhi = 0;
    std::vector<float> * vec_DPt = 0;
    std::vector<float> * vec_DRapidity = 0;
    std::vector<float> * vec_DGen = 0;
    
    theTree->SetBranchAddress("D3DDecayLength", &vec_D3DDecayLength);
    theTree->SetBranchAddress("D3DDecayLengthSignificance", &vec_D3DDecayLengthSignificance);
    theTree->SetBranchAddress("DTrk1Chi2n", &vec_DTrk1Chi2n);
    theTree->SetBranchAddress("DTrk1PtErr", &vec_DTrk1PtErr);
    theTree->SetBranchAddress("DTrk2Chi2n", &vec_DTrk2Chi2n);
    theTree->SetBranchAddress("DVtxProb", &vec_DVtxProb);
    theTree->SetBranchAddress("D3DPointingAngle", &vec_D3DPointingAngle);
    theTree->SetBranchAddress("DDca", &vec_DDca);
    theTree->SetBranchAddress("DTtrk1Pt", &vec_DTtrk1Pt);
    theTree->SetBranchAddress("DTrk2Pt", &vec_DTrk2Pt);
    theTree->SetBranchAddress("DTrk2PtErr", &vec_DTrk2PtErr);
    theTree->SetBranchAddress("DTrk1Eta", &vec_DTrk1Eta);
    theTree->SetBranchAddress("DTrk2Eta", &vec_DTrk2Eta);
    theTree->SetBranchAddress("DxyDCASignificanceDaugther1", &vec_DxyDCASignificanceDaugther1);
    theTree->SetBranchAddress("DxyDCASignificanceDaugther2", &vec_DxyDCASignificanceDaugther2);
    theTree->SetBranchAddress("DzDCASignificanceDaugther1", &vec_DzDCASignificanceDaugther1);
    theTree->SetBranchAddress("DzDCASignificanceDaugther2", &vec_DzDCASignificanceDaugther2);
    theTree->SetBranchAddress("DMass", &vec_DMass);
    theTree->SetBranchAddress("DPhi", &vec_DPhi);
    theTree->SetBranchAddress("DPt", &vec_DPt);
    theTree->SetBranchAddress("DRapidity", &vec_DRapidity);
    theTree->SetBranchAddress("DGen", &vec_DGen);

    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    //this is a loop in the events (each collision)
    std::vector<float> aux_vec_all_trainingVariables;

    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

        if (ievt%10000 == 0){
            std::cout << "--- ... Processing event: " << ievt << std::endl;
        }

        theTree->GetEntry(ievt);

	//this is a loop in the D0 meson candidates in a given event
	//NB.: here it can include signal and background 
        for(int iD0=0; iD0<vec_D3DDecayLength->size(); iD0++){
            //use same sequence of variables as in the weight file

            aux_vec_all_trainingVariables.push_back((*vec_D3DDecayLength)[iD0]); 
	        if (!TMath::IsNaN((*vec_D3DDecayLengthSignificance)[iD0])) aux_vec_all_trainingVariables.push_back((*vec_D3DDecayLengthSignificance)[iD0]);
	        aux_vec_all_trainingVariables.push_back((*vec_DTrk1Chi2n)[iD0]);
	        aux_vec_all_trainingVariables.push_back((*vec_DTrk1PtErr)[iD0]);
	        aux_vec_all_trainingVariables.push_back((*vec_DTrk2Chi2n)[iD0]);
	        aux_vec_all_trainingVariables.push_back((*vec_DVtxProb)[iD0]);
	        aux_vec_all_trainingVariables.push_back((*vec_D3DPointingAngle)[iD0]);
	        aux_vec_all_trainingVariables.push_back((*vec_DDca)[iD0]);
            aux_vec_all_trainingVariables.push_back((*vec_DTtrk1Pt)[iD0]);
            aux_vec_all_trainingVariables.push_back((*vec_DTrk2Pt)[iD0]);
            aux_vec_all_trainingVariables.push_back((*vec_DTrk2PtErr)[iD0]);
            aux_vec_all_trainingVariables.push_back((*vec_DTrk1Eta)[iD0]);
            aux_vec_all_trainingVariables.push_back((*vec_DTrk2Eta)[iD0]);
	        if (!TMath::IsNaN((*vec_DxyDCASignificanceDaugther1)[iD0])) aux_vec_all_trainingVariables.push_back((*vec_DxyDCASignificanceDaugther1)[iD0]);
            if (!TMath::IsNaN((*vec_DxyDCASignificanceDaugther2)[iD0])) aux_vec_all_trainingVariables.push_back((*vec_DxyDCASignificanceDaugther2)[iD0]);
            if (!TMath::IsNaN((*vec_DzDCASignificanceDaugther1)[iD0])) aux_vec_all_trainingVariables.push_back((*vec_DzDCASignificanceDaugther1)[iD0]);
            if (!TMath::IsNaN((*vec_DzDCASignificanceDaugther2)[iD0])) aux_vec_all_trainingVariables.push_back((*vec_DzDCASignificanceDaugther2)[iD0]);

            // Return the MVA outputs and fill into histograms
	        // See method here: https://root.cern.ch/root/html608/Reader_8cxx_source.html#l00486

            if (Use["BDT"]) {
                float BDT_score = reader->EvaluateMVA(aux_vec_all_trainingVariables,"BDT method");

                histBdt_all->Fill(BDT_score);

                if ((*vec_DGen)[iD0]==23333 || (*vec_DGen)[iD0]==23344){              
                    histBdt_signal->Fill(BDT_score);
                }

                if ((*vec_DGen)[iD0]!=23333 && (*vec_DGen)[iD0]!=23344){
                    histBdt_background->Fill(BDT_score);
                }

                if (!TMath::IsNaN((*vec_DMass)[iD0])){
                    histBdt_raw_DMass->Fill((*vec_DMass)[iD0]);
                }

                if (!TMath::IsNaN((*vec_DPhi)[iD0])){
                    histBdt_raw_DPhi->Fill((*vec_DPhi)[iD0]);
                }

                if (!TMath::IsNaN((*vec_DPt)[iD0])){
                    histBdt_raw_DPt->Fill((*vec_DPt)[iD0]);
                }

                if (!TMath::IsNaN((*vec_DRapidity)[iD0])){
                    histBdt_raw_DRapidity->Fill((*vec_DRapidity)[iD0]);
                }

                if (BDT_score > 0){
                    // DMass cut
                    histBdt_DMass_cut->Fill((*vec_DMass)[iD0]);
                    // BDT>0.0 cut to calculate statistic significance
                    if ((*vec_DGen)[iD0]==23333 || (*vec_DGen)[iD0]==23344) h_bdtSignal->Fill(BDT_score);
                    else h_bdtBackground->Fill(BDT_score);
                }
            }
            //clean up the vectors to fill next D0 meson information
            aux_vec_all_trainingVariables.clear();
	} 
    }

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    // Create TCanvas and write histograms
    
    // All variables canvas
    TCanvas c_all_var("c_all_var","c_all_var", 1000, 1000);
    c_all_var.Divide(2, 2, 0.01, 0.01);
    c_all_var.cd(1);
        histBdt_all->Draw();
    c_all_var.cd(2);
        histBdt_signal->Draw();
    c_all_var.cd(3);
        histBdt_background->Draw();

    // Kinematic variables canvas
    TCanvas c_raw_kvar("c_raw_kvar","c_raw_kvar", 1000, 1000);
    c_raw_kvar.Divide(2, 2, 0.01, 0.01);
    c_raw_kvar.cd(1);
        histBdt_raw_DMass->Draw();
    c_raw_kvar.cd(2);
        histBdt_raw_DPhi->Draw();
    c_raw_kvar.cd(3);
        histBdt_raw_DPt->Draw();
    c_raw_kvar.cd(4);
        histBdt_raw_DRapidity->Draw();

    TCanvas c_DMass_cut("c_DMass_5","c_DMass_5", 1000, 1000);
        histBdt_DMass_cut->Draw();

    // Set and print a canvas containing the histograms of BDT>0.0 cut, for checking purposes
    TCanvas c_bdt_cut("BDT>0.0", "BDT>0.0", 2000, 1000);
    c_bdt_cut.Divide(2, 1, 0.01, 0.01);
    c_bdt_cut.cd(1); h_bdtSignal->Draw();
    c_bdt_cut.cd(2); h_bdtBackground->Draw();
    c_bdt_cut.Print("BDT>0.pdf");

    // Calculate and print the statistic significance
    float eff = h_bdtSignal->Integral()/(TMath::Sqrt(h_bdtSignal->Integral()+h_bdtBackground->Integral()));
    cout << "Statistic significance: " << eff << endl;

    // Write histograms in .root file
    TFile *target  = new TFile( "TMVApp.root","RECREATE" );
    if (Use["BDT"]) {
        histBdt_all->Write(); histBdt_signal->Write(); histBdt_background->Write();
        histBdt_raw_DMass->Write(); histBdt_raw_DPhi->Write(); histBdt_raw_DPt->Write(); histBdt_raw_DRapidity->Write();
        histBdt_DMass_cut->Write(); c_DMass_cut.Print("histBdt_DMass_cut.pdf");
        c_all_var.Write(); c_raw_kvar.Write();
        c_all_var.Print("c_all_var.pdf"); c_raw_kvar.Print("c_raw_kvar.pdf");
    }

    target->Close();
    std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}

int main( int argc, char** argv ) {
    TString methodList;
    for (int i=1; i<argc; i++) {
        TString regMethod(argv[i]);
        if(regMethod=="-b" || regMethod=="--batch") continue;
        if (!methodList.IsNull()) methodList += TString(",");
        methodList += regMethod;
    }
    TMVAClassificationApplication(methodList);
    return 0;
}
