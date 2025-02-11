#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include <random>
#include <TSystem.h>
#include <TTreeReaderArray.h>

using namespace std;

const string format = ".pdf";

struct Data{int nCol, nPart; double d;};

void drawTH2 (TH2 *th2, TCanvas *canvas, const string& xTitle, const string& yTitle, const string& saveAs){
    th2->GetXaxis()->SetTitle(xTitle.c_str());
    th2->GetYaxis()->SetTitle(yTitle.c_str());
    th2->SetContour(10000);
    if (gPad) gPad->SetRightMargin(0.12);
    if (gPad) gPad->SetLeftMargin(0.13);
    if (gPad) gPad->SetBottomMargin(0.12);
    th2->DrawClone("colz");
    canvas->SaveAs(saveAs.c_str());
    canvas->Clear();

}
void drawTH1 (TH1 *h, TCanvas *canvas, const string& xTitle, const string& yTitle, const string& saveAs){
    h->GetXaxis()->SetTitle(xTitle.c_str());
    h->GetYaxis()->SetTitle(yTitle.c_str());
    if (gPad) gPad->SetRightMargin(0.05);
    if (gPad) gPad->SetBottomMargin(0.12);
    h->DrawClone("");
    canvas->SaveAs(saveAs.c_str());
    canvas->Clear();
}

// Main code
void process (const char* dataFile = "./data.root",
              const string& location = "./graphs"){

    gSystem->mkdir(location.c_str());

    // Config the canvas

    auto *c = new TCanvas("canvas", "canvas", 700, 700);

    c->SetTicks();
    c->SetGrid();

    // ColorBlind Friendly palette
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Open the ROOT file containing generated data and get the “save” tree
    auto *file = new TFile(dataFile);
    auto *tree = (TTree*) file->Get("Data");

    // Create variables to hold the data
    Data data{};

    // Set up branches to read the data into the variables
    tree->SetBranchAddress("Collisions", &data);

    int nPartMax = 2500, nCollMax = 418, bMax = 18;

    auto partXEvent = new TH1I("ppe", "", 209, 0, nCollMax);
    auto colXEvent  = new TH1I("cpe", "", 209, 0, nPartMax);
    auto distXEvent = new TH1I("dpe", "", 209, 1e-4, bMax);
    auto distXPart  = new TH2I("dpp", "", 209,0,bMax,209,1e-4, nCollMax);
    auto distXCol   = new TH2I("dpc", "", 209,0,bMax,209,1e-4,nPartMax);
    auto colXPart   = new TH2I("cpp", "", 209,0,nCollMax,nPartMax,1e-4,nPartMax);

    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        partXEvent -> Fill(data.nPart);
        colXEvent  -> Fill(data.nCol);
        distXEvent -> Fill(data.d);
        distXPart  -> Fill(data.d, data.nPart);        distXCol   -> Fill(data.d, data.nCol);
        colXPart   -> Fill(data.nPart, data.nCol);
    }

    // Draw TH2 functions
    c->SetLogz();

    drawTH2(distXPart, c, "Impact parameter (fm)", "Number of participants", location + "/DistPart" + format);
    drawTH2(distXCol, c, "Impact parameter (fm)", "Number of binary collisions", location + "/DistCol" + format);
    drawTH2(colXPart, c, "Number of binary collisions", "Number of participants", location + "/ColPart" + format);


    // Draw TH1 Functions
    if (gPad) gPad->SetLeftMargin(0.14);
    drawTH1(distXEvent, c, "b (fm)", "Number of Events", location + "/DistEvent" + format);
    c->SetLogy();
    if (gPad) gPad->SetLeftMargin(0.13);
    drawTH1(partXEvent, c, "N_{part}", "Number of Events", location + "/PartEvent" + format);
    if (gPad) gPad->SetLeftMargin(0.13);
    drawTH1(colXEvent, c, "N_{col}", "Number of Events", location + "/ColEvent" + format);
}
