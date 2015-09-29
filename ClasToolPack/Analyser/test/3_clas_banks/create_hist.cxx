#if !defined(__CINT__)
#include "TIdentificator.h"
#include "TClasTool.h"
#include "TH2F.h"
#include "TProfile.h"
#endif


void create_hist(void)
{
    // Create 2D histograms and profiles of momentum vs betta for every kind
    // of particle in the given ROOT file.

#if defined(__CINT__)
    gROOT->Reset();
    gSystem->Load("libClasTool.so");
    gSystem->Load("libTIdentificator.so");
#endif

    TClasTool *input = new TClasTool();

    input->InitDSTReader("ROOTDSTR");
    input->AddFile("clas_42011_01_1.pass2.root");

    TIdentificator *t = new TIdentificator(input);

    TH2F *h1 = new TH2F("electron_2d",  "Betta vs Momentum", 1000, 0.0, 5.0, 50, 0.0, 1.2);
    TH2F *h2 = new TH2F("he_pion_2d",   "Betta vs Momentum", 1000, 0.0, 5.0, 50, 0.0, 1.2);
    TH2F *h3 = new TH2F("le_pion_2d",   "Betta vs Momentum", 1000, 0.0, 5.0, 50, 0.0, 1.2);
    TH2F *h4 = new TH2F("le_proton_2d", "Betta vs Momentum", 1000, 0.0, 5.0, 50, 0.0, 1.2);
    TH2F *h5 = new TH2F("positron_2d",  "Betta vs Momentum", 1000, 0.0, 5.0, 50, 0.0, 1.2);

    TProfile *p1 = new TProfile("electron_pro",  "Betta vs Momentum", 1000, 0.0, 5.0, 0.0, 1.2);
    TProfile *p2 = new TProfile("he_pion_pro",   "Betta vs Momentum", 1000, 0.0, 5.0, 0.0, 1.2);
    TProfile *p3 = new TProfile("le_pion_pro",   "Betta vs Momentum", 1000, 0.0, 5.0, 0.0, 1.2);
    TProfile *p4 = new TProfile("le_proton_pro", "Betta vs Momentum", 1000, 0.0, 5.0, 0.0, 1.2);
    TProfile *p5 = new TProfile("positron_pro",  "Betta vs Momentum", 1000, 0.0, 5.0, 0.0, 1.2);

    Int_t nEntries = input->GetEntries();

    for (Int_t k = 0; k < nEntries; k++) {
        input->Next();
        Int_t nRows = input->GetNRows("EVNT");
        for (Int_t i = 0; i < nRows; i++) {
            TString category = t->GetCategorization(i);
            if (category == "electron") {
                h1->Fill(t->Momentum(i), t->Betta(i));
                p1->Fill(t->Momentum(i), t->Betta(i));
            } else if (category == "high energy pion +") {
                h2->Fill(t->Momentum(i), t->Betta(i));
                p2->Fill(t->Momentum(i), t->Betta(i));
            } else if (category == "low energy pion +") {
                h3->Fill(t->Momentum(i), t->Betta(i));
                p3->Fill(t->Momentum(i), t->Betta(i));
            } else if (category == "low energy proton") {
                h4->Fill(t->Momentum(i), t->Betta(i));
                p4->Fill(t->Momentum(i), t->Betta(i));
            } else if (category == "positron") {
                h5->Fill(t->Momentum(i), t->Betta(i));
                p5->Fill(t->Momentum(i), t->Betta(i));
            }
        }
    };

    TFile *output_h = new TFile("particle_histograms.root","RECREATE","2D Histograms of particles");
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h5->Write();
    output_h->Close();

    TFile *output_p = new TFile("particle_profiles.root","RECREATE","Histogram Profiles of particles");
    p1->Write();
    p2->Write();
    p3->Write();
    p4->Write();
    p5->Write();
    output_p->Close();

    std::cout << "Done." << std::endl;
}



#if !defined(__CINT__)
int main(int argc, char **argv)
{
    create_hist();
    return 0;
}
#endif
