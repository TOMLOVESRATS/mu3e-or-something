//Getting libraries that will be used 
#include "TCanvas.h"
//root canvas 
#include "TH1.h"
#include "TString.h"
//root 1d histogram
#include "TF1.h"
//root 1d function
#include "TF2.h"
//root 2d function
#include "TH2F.h"
//2D histogram
#include "TH3F.h"
//3D histogram.
#include "TRandom.h"
//ROOT’s random number generator
#include "TGraph2D.h"
// scatter plot
#include "TGraph2DErrors.h"
// with error bars
#include <TLegend.h>
//legend box
#include <TRandom3.h>
//root random generator again 
#include <TFile.h>
// used to save histogram
// These are normal c++ library stuff
#include <cmath>
//inculde math class
#include <random>
//Loads C++’s own random generator
#include <ctime>
//load time for some reason
#include <limits>
// constants for limit value 
#include "TSystem.h"
#include "TROOT.h"


// ----------------------------------------
// Global settings
// ----------------------------------------

double xMinPlot = 10.0 / 53.0;
double thetaMinPlot = 0;
double thetaMaxPlot = M_PI;

//root random generation 
TRandom3 randMomentum(0);
//empty stuff for storage
// Asymmetry histograms
TH1D* hAfb_vs_x      = nullptr;   // A_FB vs xp (reduced momentum/energy)
TH1D* hAfb_vs_costh  = nullptr;   // A_FB vs cos(theta_cut)
#include <vector>

// ----------------------------------------
// Michel spectrum function
// ----------------------------------------

double michel(double *x, double *par)
{
    double xp    = x[0];     // reduced energy
    double theta = x[1];     // emission angle

    // Michel parameters given in the 2d histogram part 
    double rho     = par[0];
    double eta     = par[1];
    double epsilon = par[2];
    double delta   = par[3];
    double P_u     = std::abs(par[4]);

    // PDG constant
    double x_0 = 9.67e-3;

    // Michel distribution
    double value =
        pow(xp, 2) *
        (3.0 * (1.0 - xp)
        + 2.0 * rho * (4.0 * xp / 3.0 - 1.0)
        + 3.0 * eta * x_0 * (1.0 - xp) / xp
        + std::abs(P_u) * epsilon * cos(theta)
          * (1.0 - xp + (2.0 / 3.0) * delta * (4.0 * xp - 3.0)));

    return value > 0 ? value : 0;
}

//fill A_FB vs xp and A_FB vs cos(theta_cut)

void fillAsymmetryHists(TH2D* fakeMCHistogram,
                        int   xBins,
                        double xMin,
                        double xMax)
{
    int thetaBins = fakeMCHistogram->GetNbinsY();

    // ---------- A_FB vs x_p ----------
    if (!hAfb_vs_x) {
        hAfb_vs_x = new TH1D("hAfb_vs_x",
                             "Forward-Backward Asymmetry vs x_{p};x_{p};A_{FB}",
                             xBins, xMin, xMax);
    }
    hAfb_vs_x->Reset();

    for (int ix = 1; ix <= xBins; ++ix) {
        double N_forward  = 0.0;
        double N_backward = 0.0;

        for (int it = 1; it <= thetaBins; ++it) {
            double thetaCenter = fakeMCHistogram->GetYaxis()->GetBinCenter(it);
            double costh = std::cos(thetaCenter);

            double content = fakeMCHistogram->GetBinContent(ix, it);
            if (content <= 0.0) continue;

            if (costh > 0.0)
                N_forward  += content;
            else if (costh < 0.0)
                N_backward += content;
        }

        double Afb = 0.0;
        if (N_forward + N_backward > 0.0)
            Afb = (N_forward - N_backward) / (N_forward + N_backward);

        hAfb_vs_x->SetBinContent(ix, Afb);
    }


// ---------- A_FB vs |cos(theta)|
if (!hAfb_vs_costh) {
    int nCosBins = 100; // |cosθ| from 0 to 1
    hAfb_vs_costh = new TH1D("hAfb_vs_costh",
                             "Forward-Backward Asymmetry vs |cos#theta|;|cos#theta|;A_{FB}",
                             nCosBins, 0.0, 1);
}
hAfb_vs_costh->Reset();

int nCosBins = hAfb_vs_costh->GetNbinsX();

// forward / backward counts per |cosθ| slice
std::vector<double> Nf_slice(nCosBins + 1, 0.0);
std::vector<double> Nb_slice(nCosBins + 1, 0.0);

// loop over theta bins
for (int it = 1; it <= thetaBins; ++it) {

    // use BIN CENTRE for theta
    double thetaCenter = fakeMCHistogram->GetYaxis()->GetBinCenter(it);
    double costh       = std::cos(thetaCenter);
    double absCos      = std::fabs(costh);

    int binAbs = hAfb_vs_costh->GetXaxis()->FindBin(absCos);

    // integrate over ALL xp for this θ-bin
    for (int ix = 1; ix <= xBins; ++ix) {

        double content = fakeMCHistogram->GetBinContent(ix, it);
        if (content <= 0.0) continue;

        if (costh > 0.0)
            Nf_slice[binAbs] += content;  // forward
        else if (costh < 0.0)
            Nb_slice[binAbs] += content;  // backward
        // costh == 0 contributes to neither
    }
}

// convert F/B counts into A_FB per |cosθ| bin
for (int b = 1; b <= nCosBins; ++b) {
    double F = Nf_slice[b];
    double B = Nb_slice[b];

    double Afb = 0.0;
    if (F + B > 0.0)
        Afb = (F - B) / (F + B);

    hAfb_vs_costh->SetBinContent(b, Afb);
}


}

// ----------------------------------------
// Create theoretical 2D histogram and fake Mc Simulation 
// ----------------------------------------
//compute for a given polarisaiton value
void createTheoreticalHistogram(long long numberPositrons,
                                double P_UChange,
                                int xBins, int thetaBins,
                                const char* theoryDir,
                                const char* fakeDir,
                                const char* symmetryDir)
{
    //setting x and theta for bining 
    double xMin = xMinPlot;
    double xMax = 1.0;
    double thetaMin = thetaMinPlot;
    double thetaMax = thetaMaxPlot;
// michelk parament packed into parSM 
    double rhoSM = 0.75, etaSM = 0, epsilonSM = 1, deltaSM = 0.75;
    double parSM[5] = {rhoSM, etaSM, epsilonSM, deltaSM, P_UChange};
//Computes the bin width in x and θ
    double binWidthX     = (xMax - xMin) / xBins;
    double binWidthTheta = (thetaMax - thetaMin) / thetaBins;
    TString hname_theory;
    hname_theory.Form("Theoretical_Histogram_SM_P_%+.1f", P_UChange);
//create the 2d histogram 
    TH2D *theoreticalSMHistogram =
        new TH2D(hname_theory,
                 "Theoretical Michel Decay; Reduced Energy; Theta; Events",
                 xBins, xMin, xMax,
                thetaBins, thetaMin, thetaMax);

//This will accumulate the total (integrated) value of the function across all bins.
    double totalEvents = 0;
//loop over bin
// i for x-bin j for theta
    for (int i = 1; i <= xBins; ++i) {
        for (int j = 1; j <= thetaBins; ++j) {
            // getting x0center and theta center 
            double xp    = theoreticalSMHistogram->GetXaxis()->GetBinCenter(i);
            double theta = theoreticalSMHistogram->GetYaxis()->GetBinCenter(j);
            //pack xp and theta to pass in to michel function 
            double input[2] = {xp, theta};
            double value    = michel(input, parSM);
            //value is from the michale spectrum
            //multipy by the bin area to get approximate number of events in that bin 
            double binContent = value * binWidthX * binWidthTheta;
            //update to total events
            totalEvents += binContent;
            //this histogram hold only the theoretical count 
            theoreticalSMHistogram->SetBinContent(i, j, binContent);
        }
    }

    // Normalise to requested total positrons
    if (totalEvents > 0)
        theoreticalSMHistogram->Scale(numberPositrons / totalEvents);

    theoreticalSMHistogram->SetStats(0);

    //fake monte carlo
    TString hname_fake;
    hname_fake.Form("FakeMC_Histogram_SM_P_%+.1f", P_UChange);
    TH2D *fakeMCHistogram =
        (TH2D*)theoreticalSMHistogram->Clone(hname_fake);
    fakeMCHistogram->SetTitle("Fake MC Michel Decay;Reduced Energy;Theta;Events"); 

    // Apply Gaussian fluctuations bin-by-bin
    for (int i = 1; i <= xBins; ++i) {
        for (int j = 1; j <= thetaBins; ++j) {

            double binContent = theoreticalSMHistogram->GetBinContent(i, j);
            double uncertainty = (binContent > 0.0) ? std::sqrt(binContent) : 0.0;
            double newBinContent = randMomentum.Gaus(binContent, uncertainty);
        // Force positivity (zero is OK)
        if (newBinContent < 0.0)
            newBinContent = 0.0;

            fakeMCHistogram->SetBinContent(i, j, newBinContent);
            fakeMCHistogram->SetBinError(i, j, uncertainty);
        }
    }
        // ---------------------------------------------------------
    // Create comparison histograms: difference and pull maps
    // ---------------------------------------------------------
    {
        // ----- Difference: FakeMC - Theory -----
        TString hname_diff;
        hname_diff.Form("Diff_FakeMinusTheory_SM_P_%+.1f", P_UChange);

        TH2D* hDiff = (TH2D*)fakeMCHistogram->Clone(hname_diff);
        hDiff->SetTitle("FakeMC - Theory; x_{p}; #theta; #Delta N");

        // Subtract theory contents bin-by-bin (errors stay from fakeMC)
        hDiff->Add(theoreticalSMHistogram, -1.0);

        // Draw and save
        TCanvas* cDiff = new TCanvas("cDiff", "FakeMC - Theory", 800, 600);
        cDiff->cd();
        hDiff->SetStats(0);
        hDiff->Draw("COLZ");

        std::string diffPdfName =
            std::string(symmetryDir) + "/Diff_FakeMinusTheory_SM_P_" + std::to_string(P_UChange) + ".pdf";
        cDiff->SaveAs(diffPdfName.c_str());

        std::string diffRootName =
            std::string(symmetryDir) + "/Diff_FakeMinusTheory_SM_P_" + std::to_string(P_UChange) + ".root";
        {
            TFile fDiff(diffRootName.c_str(), "RECREATE");
            hDiff->Write("diff");
        }
        delete cDiff;

        // ----- Pull map: (FakeMC - Theory) / sqrt(Theory) -----
        TString hname_pull;
        hname_pull.Form("Pull_SM_P_%+.1f", P_UChange);

        TH2D* hPull = (TH2D*)fakeMCHistogram->Clone(hname_pull);
        hPull->SetTitle("Pull: (FakeMC - Theory) / #sqrt{Theory}; x_{p}; #theta; pull");

        for (int i = 1; i <= xBins; ++i) {
            for (int j = 1; j <= thetaBins; ++j) {
                double nTh  = theoreticalSMHistogram->GetBinContent(i, j);
                double nMC  = fakeMCHistogram->GetBinContent(i, j);
                double diff = nMC - nTh;

                double sigma = (nTh > 0.0) ? std::sqrt(nTh) : 0.0;
                double pull  = 0.0;
                if (sigma > 0.0) pull = diff / sigma;

                hPull->SetBinContent(i, j, pull);
                hPull->SetBinError(i, j, 0.0); // often we don't use errors on pull bins
            }
        }

        TCanvas* cPull = new TCanvas("cPull", "Pull Map", 800, 600);
        cPull->cd();
        hPull->SetStats(0);
        hPull->Draw("COLZ");

        std::string pullPdfName =
            std::string(symmetryDir) + "/PullMap_SM_P_" + std::to_string(P_UChange) + ".pdf";
        cPull->SaveAs(pullPdfName.c_str());

        std::string pullRootName =
            std::string(symmetryDir) + "/PullMap_SM_P_" + std::to_string(P_UChange) + ".root";
        {
            TFile fPull(pullRootName.c_str(), "RECREATE");
            hPull->Write("pull");
        }
        delete cPull;
    }

    // Save plot
    //drawing the graph with c being the canvas name
    // --- Draw and save theory surface ---
    {
        TCanvas *c = new TCanvas("c_theory", "Theoretical Michel Decay SM", 800, 600);
        c->cd();

        TString title;
        title.Form("Theoretical Michel Decay (P = %.2f, N = 10^{15} e^{+});x_{p};#theta;Events",
                   P_UChange);
        theoreticalSMHistogram->SetTitle(title);
        theoreticalSMHistogram->Draw("SURF1");

        std::string pdfName =
            std::string(theoryDir) + "/Theoretical_Graph_SM_P_" + std::to_string(P_UChange) + ".pdf";
        c->SaveAs(pdfName.c_str());

        std::string rootName =
            std::string(theoryDir) + "/Theoretical_SM_Root_P_" + std::to_string(P_UChange) + ".root";
        TFile *f = new TFile(rootName.c_str(), "RECREATE");
        theoreticalSMHistogram->Write("theory");
        f->Close();
    }

    // --- Draw and save fake MC surface ---
    {
        TCanvas *c_fake = new TCanvas("c_fake", "Fake MC Michel Decay SM", 800, 600);
        c_fake->cd();

        TString title_fake;
        title_fake.Form("Fake MC Michel Decay (P = %.2f, N = 10^{15} e^{+});x_{p};#theta;Events",
                        P_UChange);
        fakeMCHistogram->SetTitle(title_fake);
        fakeMCHistogram->Draw("SURF1");

        std::string pdfNameFake =
            std::string(fakeDir) + "/FakeMC_Graph_SM_P_" + std::to_string(P_UChange) + ".pdf";
        c_fake->SaveAs(pdfNameFake.c_str());

        std::string rootNameFake =
            std::string(fakeDir) + "/FakeMC_SM_Root_P_" + std::to_string(P_UChange) + ".root";
        TFile *fFake = new TFile(rootNameFake.c_str(), "RECREATE");
        fakeMCHistogram->Write("fakeMC");
        fFake->Close();
    }

    // --- Fill asymmetry histograms from this fake MC ---
    fillAsymmetryHists(fakeMCHistogram, xBins, xMin, xMax);
// ------------------------------------------------------------
// Draw and save A_FB vs x  (momentum / xp)
// ------------------------------------------------------------
{
    TCanvas *c_afb_x = new TCanvas("c_afb_x", "A_FB vs x", 800, 600);
    c_afb_x->cd();

    TString titleAfbX;


    // Draw without error bars: use "P" (points only)
    hAfb_vs_x->Draw("P");

    // Save PDF in symmetry folder
    std::string pdfNameAfbX =
        std::string(symmetryDir) + "/Afb_vs_x_P_" + std::to_string(P_UChange) + ".pdf";
    c_afb_x->SaveAs(pdfNameAfbX.c_str());

    // Save ROOT in symmetry folder
    std::string rootNameAfbX =
        std::string(symmetryDir) + "/Afb_vs_x_P_" + std::to_string(P_UChange) + ".root";
    TFile *fAfbX = new TFile(rootNameAfbX.c_str(), "RECREATE");
    hAfb_vs_x->Write("Afb_vs_x");
    fAfbX->Close();
}
// ------------------------------------------------------------
// Draw and save A_FB vs cos(theta)
// ------------------------------------------------------------
{
    TCanvas *c_afb_cos = new TCanvas("c_afb_cos", "A_FB vs cos(theta)", 800, 600);
    c_afb_cos->cd();

    TString titleAfbCos;


    
    // Draw without error bars
    hAfb_vs_costh->Draw("P");
    hAfb_vs_costh->SetMarkerStyle(2);
    // Save PDF in symmetry folder
    std::string pdfNameAfbCos =
        std::string(symmetryDir) + "/Afb_vs_costh_P_" + std::to_string(P_UChange) + ".pdf";
    c_afb_cos->SaveAs(pdfNameAfbCos.c_str());

    // Save ROOT in symmetry folder
    std::string rootNameAfbCos =
        std::string(symmetryDir) + "/Afb_vs_costh_P_" + std::to_string(P_UChange) + ".root";
    TFile *fAfbCos = new TFile(rootNameAfbCos.c_str(), "RECREATE");
    hAfb_vs_costh->Write("Afb_vs_costh");
    fAfbCos->Close();
}


}
// ----------------------------------------
// Main
// ----------------------------------------
int main()
{
    const char* outDir = "/home/tom/Mu3e/Photos/Theory events photos";

    std::string theoryDir   = std::string(outDir) + "/Theory";
    std::string fakeDir     = std::string(outDir) + "/FakeMC";
    std::string symmetryDir = std::string(outDir) + "/Symmetry";

    ::gSystem->mkdir(theoryDir.c_str(),   true);
    ::gSystem->mkdir(fakeDir.c_str(),     true);
    ::gSystem->mkdir(symmetryDir.c_str(), true);

    long long numberOfPositrons = 1000000000000000LL; // 1e15
    int xBins = 100;
    int yBins = 100;

    // ======== LOOP OVER P VALUES ========
    for (double Ptest = -1.0; Ptest <= 1.0; Ptest += 0.1)
    {
        createTheoreticalHistogram(numberOfPositrons,
                                   Ptest,
                                   xBins, yBins,
                                   theoryDir.c_str(),
                                   fakeDir.c_str(),
                                   symmetryDir.c_str()); 
    }


    return 0;
}



