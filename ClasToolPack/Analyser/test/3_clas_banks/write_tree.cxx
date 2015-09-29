#if !defined(__CINT__)
#include "TIdentificator.h"
#include "TClasTool.h"
#include "TH2F.h"
#endif


void write_tree()
{
    // Create a tree with a branch for betta, momentum and the kind of
    // particle (id), for every entry in the given ROOT file.

#if defined(__CINT__)
    gROOT->Reset();
    gSystem->Load("libClasTool.so");
    gSystem->Load("libTIdentificator.so");
#endif

    TClasTool *input = new TClasTool();

    input->InitDSTReader("ROOTDSTR");
    input->AddFile("clas_42011_01_1.pass2.root");

    TIdentificator *t = new TIdentificator(input);

    TFile *output = new TFile("particle_data.root", "RECREATE", "Data of particles");
    TTree *tree = new TTree("data", "Tree that holds the data");

    Double_t betta, moment;
    Int_t id;

    tree->Branch("betta", &betta, "betta/D");
    tree->Branch("moment", &moment, "moment/D");
    tree->Branch("particle", &id, "id/I");

    Int_t nEntries = input->GetEntries();

    for (Int_t k = 0; k < nEntries; k++) {
        input->Next();
        Int_t nRows = input->GetNRows("EVNT");
        for (Int_t i = 0; i < nRows; i++) {
            TString category = t->GetCategorization(i);
            if (category == "electron")
                id = 1;
            else if (category == "high energy pion +")
                id = 2;
            else if (category == "low energy pion +")
                id = 3;
            else if (category == "low energy proton")
                id = 4;
            else if (category == "positron")
                id = 5;
            else
                id = 0;
            moment = t->Momentum(i);
            betta = t->Betta(i);
            tree->Fill();
        }
    };

    output->Write();
    output->Close();
    cout << "Done." << endl;
}



#if !defined(__CINT__)
int main(int argc, char **argv)
{
    write_tree();
    return 0;
}
#endif
