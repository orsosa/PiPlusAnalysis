#if !defined(__CINT__)
#include "TIdentificator.h"
#include "TClasTool.h"
#include "TH2F.h"
#include "TProfile.h"
#endif


void mass_distribution(void)
{
#if defined(__CINT__)
    gROOT->Reset();
    gSystem->Load("libClasTool.so");
    gSystem->Load("libTIdentificator.so");
#endif

    TClasTool *input = new TClasTool();

    input->InitDSTReader("ROOTDSTR");
    input->AddFile("clas_42011_01_1.pass2.root");

    TIdentificator *t = new TIdentificator(input);

    TH2F *h = new TH2F("mass_dist", "Mass Distribution", 1000, 0.0, 5.0, 50, 0.0, 1.2);
    TProfile *p = new TProfile("mass_dist_p", "Mass Distribution", 1000, 0.0, 5.0, 0.0, 1.2);

    Int_t nEntries = input->GetEntries();
    for (Int_t k = 0; k < nEntries; k++) {
        input->Next();
        Int_t nRows = input->GetNRows("EVNT");
        for (Int_t i = 0; i < nRows; i++)
            if (t->Charge(i) == 1 && t->StatSC(i) > 0) {
                h->Fill(t->Mass2(i), t->Momentum(i));
                p->Fill(t->Mass2(i), t->Momentum(i));
            }
    };

    TFile *output_m = new TFile("particle_mass_dist.root","RECREATE","Particle Mass Distribution");
    h->Write();
    p->Write();
    output_m->Close();

    cout << "Done." << endl;
}



#if !defined(__CINT__)
int main(int argc, char **argv)
{
    mass_distribution();
    return 0;
}
#endif
