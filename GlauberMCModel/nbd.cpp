#include <TGraph2D.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TLine.h>
#include <TTree.h>
#include <TText.h>
#include <TH2F.h>
#include <TF1.h>
#include "TGaxis.h"
#include <TSystem.h>
#include "TLegend.h"
#include <TCanvas.h>
#include <TStyle.h>


using namespace std;

// Structures used to generate and save data.
struct Data {int nCol, nPart; double d;};

<<<<<<< HEAD
// Number of normalization values to test
constexpr int nDim = 30;

// Normalization range to test
constexpr double nMax = 1e-1;
constexpr double nMin = 5e-3;
constexpr double nDelta = (nMax - nMin) / (nDim - 1);
=======
void nbd(const char* dataFile = "./data.root",
         const char* compareFile = "./HFSumEt.root"){

    double hScale = 1./1000.;
    const int scanSize = 5;
    const int nDim = 10;
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437

// Number of each parameter value to test
int scanSize = 8;

// K and Mu intervals to test
double kMin = 0.6;
double kMax = 1.2;
double muMin = 1.4;
double muMax = 1.8;
const double kDelta = (kMax - kMin) / (scanSize - 1);
const double muDelta = (muMax - muMin) / (scanSize - 1);

// Bins to perform the chi2 test
constexpr int lowerBin = 10;
constexpr int upperBin = 100;

// Divisions to calculate the partitions
constexpr int minDiv = 5;
constexpr int nDiv = 20;

// Function to draw a line in a histogram
void drawLine (double x[], const int i, const double max, double energy) {
    const double yt = sqrt(log10(max));
    const double tx = 100-(static_cast<double>(minDiv)/nDiv
                              + static_cast<double>(i)/nDiv)*100;
    const double xt = (tx-100./nDiv>0.1) ? (x[i] + x[i+1])/2. : (x[i] + energy) / 2;
    auto *text = new TText(xt, yt,Form("%.0f-%.0f%%", tx, tx-100./nDiv));
    auto *line = new TLine(x[i],0,x[i],max);
    line->SetLineColor(kGreen + 2);
    line->SetLineStyle(9);
    line->Draw("SAME");
    text->SetTextSize(.03);
    text->SetTextAngle(90);
    text->SetTextAlign(22);
    if(tx-100./nDiv > 0.1 ? x[i+1] - x[i] > 0.12: energy -x[i]> 0.12)
        text->Draw("same");
}

void nbd(const char* dataFile = "./data.root",
         const char* compareFile = "./HFSumEt.root") {

    // Create the NBD function
    const auto NBD = new TF1("NBD", "TMath::Gamma(1000*x+[1])/(TMath::Gamma(1000*x+1)*TMath::Gamma([1])) *"
                             "(([0]/[1])**(1000*x))/(([0]/[1]+1)**(1000*x+[1]))", 0.0, 0.1);

    // Create the histogram to store the results
    TH2 *results = new TH2F("chiRes", "chiRes", scanSize, muMin - muDelta / 2,
                            muMax + muDelta / 2, scanSize, kMin - kDelta / 2, kMax + kDelta / 2);

    // Read experimental values
    cout << "\nReading experimental values: ";
    auto *expFile = new TFile(compareFile);
    auto *exp = (TH1F *) (expFile->Get("hfSum_new"));
    cout << "Done!\n";

    // Read simulated values
    cout << "Reading simulated values: ";
    auto *file = new TFile(dataFile);
    auto *tree = (TTree *) (file->Get("Data"));
    Data data{};
    tree->SetBranchAddress("Collisions", &data);
    vector<int> nCol;
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        nCol.emplace_back(data.nCol);
    }
    cout << "Done!\n";

<<<<<<< HEAD

    // Create variables to store the best values
    double bestK = 0, bestMu = 0, bestN = 0;
    double chi2best = INFINITY;

    // Run through the Mu values
    cout << "Iterating through values:\n";
    for (int iMu = 0; iMu < scanSize; iMu++) {
        const double mu = muMin + muDelta * iMu;
=======
    // Create graphing variables
    double bestK, bestMu;
    double chi2best = INFINITY;

    // Run through the Mu values
    for (int iMu = 0; iMu < scanSize; iMu++) {
        double mu = muMin + muDelta * iMu;
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437

        // Run through the K values
        for (int ik = 0; ik < scanSize; ik++) {
            double localChi2best = INFINITY;
<<<<<<< HEAD
            const double k = kMin + kDelta * ik;
            cout << "Mu value: " << mu << "   ";
            cout << "K value: " << k << "   ";
=======
            double k = kMin + kDelta * ik;
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437

            // Set the NBD parameters
            NBD->SetParameters(mu, k);
<<<<<<< HEAD
            const auto generated = new TH1F("Generated", "sumETf vs events", 100, 0, 5);
=======
            TH1F *generated = new TH1F("Generated", "sumETf vs events", 100, 0, 5);
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437

            // Generate the random data
            for (const int col: nCol) {
                double sumET = 0;
                for (int n = 0; n < col; n++) {
                    sumET += NBD->GetRandom();
                }
                generated->Fill(sumET);
            }

<<<<<<< HEAD
            // Calculate the chi2 of the given Mu and K
            // through a range of normalization values
            Double_t n = 0;
            for (int iNorm = 0; iNorm < nDim; iNorm++) {
                generated->Scale(nMin + nDelta * iNorm);
                Double_t chi2 = 0;
                for (int i = lowerBin; i < upperBin; i++) {
                    const Double_t Expected = exp->GetBinContent(i);
                    const Double_t Observed = generated->GetBinContent(i);
                    const Double_t diff = Observed - Expected;
                    if (diff == 0 || Expected == 0) continue;
                    chi2 += (diff * diff) / Expected;
                }
                if (chi2 == INFINITY) break;
                if (localChi2best > chi2) {
                    localChi2best = chi2;
=======
            // Calculate the Chi^2 ofr the given Mu and K
            // through a range of normalization values
            double chi2min;
            double n;
            for (int iNorm = 0; iNorm < nDim; iNorm++) {
                generated->Scale(nMin + nDelta * iNorm);
                chi2min = 0;
                double Observed, Expected, diff;
                for (int i = 20; i < 100; i++) {
                    Expected = (double)exp->GetBinContent(i);
                    Observed = (double)generated->GetBinContent(i);
                    diff = Observed - Expected;
                    chi2min += (diff*diff)/Expected;
                }
                if (localChi2best > chi2min) {
                    localChi2best = chi2min;
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437
                    n = nMin + nDelta * iNorm;
                }
                generated->Scale(1 / (nMin + nDelta * iNorm));
            }
<<<<<<< HEAD

            // Output the resulting values in the terminal
            results->SetBinContent(iMu + 1, ik + 1, localChi2best);
            cout
                    << "N = " << n << "   "
                    << "chi2 = " << localChi2best << "\n";

            // Saves the best K and Mu values
            if (localChi2best < chi2best) {
                bestK = k;
                bestMu = mu;
                bestN = n;
                chi2best = localChi2best;
            }
=======

            // Output the resulting values in the terminal
            results->SetBinContent(iMu + 1, ik + 1, localChi2best);
            cout
                    << "Chi2 = " << localChi2best
                    << " MU = " << mu
                    << " K = " << k
                    << " N = " << n << endl;

            // Saves the best K and Mu values
            if (localChi2best < chi2best) {bestK = k; bestMu = mu;}
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437
            localChi2best = INFINITY;
            delete generated;
        }
    }

    // Create an array to store the energy values
    const int len = static_cast<int>(nCol.size());
    float energy[len];

    // Generate a new histogram with the function configured with the best values
    cout << "Graphing results:\n";
    NBD->SetParameters(bestMu, bestK);
    const auto generated = new TH1F("Generated", "Generated energy;#Sigma E_{T} (TeV);Events per bin",
                                    100, 0, 5);


    // Generate and sort the energy values
    for (int i = 0; i <= len; i++) {
        double sumET = 0;
        for (int n = 0; n < nCol[i]; n++) {
            sumET += NBD->GetRandom();
        }
        energy[i] = static_cast<float>(sumET);
        generated->Fill(sumET);
    }
    sort(energy, energy + len);

    // Calculate the partitions
    double partitions[nDiv - minDiv];
    for (int i = minDiv; i < (nDiv); i++) {
        if ((len * i / nDiv) % 2 == 0)
            partitions[i - minDiv] = energy[len * i / nDiv];
        else
            partitions[i - minDiv] = (energy[len * i / nDiv] + energy[len * i / nDiv + 1]) / 2;
    }

    // Config the canvas
    auto *c = new TCanvas("canvas", "canvas", 1000, 743);
<<<<<<< HEAD
    TGaxis::SetMaxDigits(3);
=======
    c->SetLogy();
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437
    c->SetLogz();
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);

<<<<<<< HEAD
    // Personalize the resulting Chi2 histogram and draw it
    results->GetYaxis()->SetTitle("k (GeV)");
    results->GetXaxis()->SetTitle("#mu (GeV)");
    if (scanSize > 5){
        results->GetYaxis()->SetNdivisions(- scanSize - 600);
        results->GetXaxis()->SetNdivisions(- scanSize - 600);
    } else {
        results->GetYaxis()->SetNdivisions(-(2 * scanSize) - 600);
        results->GetXaxis()->SetNdivisions(-(2 * scanSize) - 600);
    }
=======
    // Personalize resulting Chi2h histogram and draw it
    results->GetYaxis()->SetTitle("k");
    results->GetXaxis()->SetTitle("#mu");
    results->GetYaxis()->SetNdivisions(2*scanSize+1,0,0);
    results->GetXaxis()->SetNdivisions(2*scanSize+1,0,0);
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437
    results->SetTitle("");
    results->DrawClone("colz");
    results->SetMarkerSize(1.5);
    results->Draw("same text");
<<<<<<< HEAD
    c->SaveAs("graph.pdf");

    // Personalize the comparative histogram with partitions and draw it
    c->SetLogy();
    generated->Scale(bestN);
    generated->SetLineColor(kMagenta);
    generated->SetLineWidth(2);
    generated->Draw("HIST");

    exp->SetLineColor(kBlack);
    exp->SetMarkerColor(kBlack);
    exp->SetMarkerStyle(20);
    exp->SetMarkerSize(.5);
    exp->Draw("same E0");

    auto *fitLeg = new TLegend(0.60, 0.80, 0.85, 0.88);
    fitLeg->SetTextSize(0.025);
    fitLeg->AddEntry(generated, "Simulation");
    fitLeg->AddEntry(exp, "CMS Open Data");
    fitLeg->Draw();


    c->SaveAs("resultswl.pdf");
    for (int i = 0; i < std::size(partitions);i++) {
        drawLine(partitions, i, generated->GetMaximum(), energy[len-1]);
    }
    c->SaveAs("results.pdf");
    cout << energy[len - 1] << ", " << energy[0];
    cout << "Done!\n";

    // Output the best values
    cout << "\nSmaller chi2: " << chi2best << ", Best mu: " << bestMu << ", Best K: " << bestK <<
         ", Best normalization: " << bestN << "\n\n";
}
=======
    c->SaveAs("ue.png");


    // Generate a new histogram with the best variable's values
    NBD->SetParameters(bestMu, bestK);
    TH1F *generated = new TH1F("Generated", "Energia gerada;#Sigma E_{T}", 100, 0, 5);
    for (int col : nCol) {
        double sumET = 0;
        for (int n = 0; n < col; n++) {
            sumET += NBD->GetRandom();
        }
        generated->Fill(sumET * hScale);
    }

    // Personalize the generated with the experimental histograms and draw them
    generated->Scale(0.021);
    generated->SetLineColor(kGreen);
    generated->Draw("HIST");
    exp->Draw("same E0");
    auto* fitLeg = new TLegend(0.60,0.80,0.85,0.88);
    fitLeg->SetTextSize(0.025);
    fitLeg->AddEntry(generated, "Simulacao");
    fitLeg->AddEntry(exp, "Dados experimentais");
    fitLeg->Draw();
    c->SaveAs("results.png");   

}
>>>>>>> 82e75058e8ca3415798da5809217a35bdbe01437
