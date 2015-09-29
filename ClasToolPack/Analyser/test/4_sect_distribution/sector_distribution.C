void sector_distribution()
{
    gSystem->Load("libTree.so");
    gSystem->Load("libClasTool.so");
    gSystem->Load("libTIdentificator.so");
    TClasTool *input = new TClasTool();
    input->InitDSTReader("ROOTDSTR");
    input->AddFile("clas_42011_01_1.pass2.root");
    TIdentificator *t = new TIdentificator(input);
    Int_t nEntries = input->GetEntries();

    TH1F *h[6];
    h[0]= new TH1F("Sector0","Distribucion Sector 0", 50,-40,0);
    h[1]= new TH1F("Sector1","Distribucion Sector 1", 50,-40,0);
    h[2]= new TH1F("Sector2","Distribucion Sector 2", 50,-40,0);
    h[3]= new TH1F("Sector3","Distribucion Sector 3", 50,-40,0);
    h[4]= new TH1F("Sector4","Distribucion Sector 4", 50,-40,0);
    h[5]= new TH1F("Sector5","Distribucion Sector 5", 50,-40,0);
                                    //32//20
    input->Next();
    for (int k = 0; k < nEntries; k++) {
        int nRows = input->GetNRows("EVNT");
        TString category = t->GetCategorization(0);//dejar numero de hit en 0 (primer hit)
        if (category == "electron") { //hacer con arreglo de TH2F
            int indice=t->Sector(0);
            h[indice]->Fill(t->Z(0));
        }
        input->Next();
    }
    TFile *output_m = new TFile("sect_dist.root","RECREATE","Sector Distribution");
    h[0]->Write();
    h[1]->Write();
    h[2]->Write();
    h[3]->Write();
    h[4]->Write();
    h[5]->Write();
    output_m->Close();
}
