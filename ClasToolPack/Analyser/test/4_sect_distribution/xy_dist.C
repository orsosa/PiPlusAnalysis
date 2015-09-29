void EquisvsY()
{
    gSystem->Load("libTree.so");
    gSystem->Load("libClasTool.so");
    gSystem->Load("libTIdentificator.so");
    TClasTool *input = new TClasTool();
    input->InitDSTReader("ROOTDSTR");
    input->AddFile("clas_42011_01_1.pass2.root");
    TIdentificator *t = new TIdentificator(input);
    Int_t nEntries = input->GetEntries();
    input->Next();
    TH2F *h = new TH2F("h","X vs Y",10,-20,20,1000,0.-40,40);
    TH2F *h1 = new TH2F("h1","X vs Y",10,-20,20,1000,0.-40,40);
    TH2F *h2 = new TH2F("h2","X vs Y",10,-20,20,1000,0.-40,40);
    TH2F *h3 = new TH2F("h3","X vs Y",10,-20,20,1000,0.-40,40);
    TH2F *h4 = new TH2F("h4","X vs Y",10,-20,20,1000,0.-40,40);
    TH2F *h5 = new TH2F("h5","X vs Y",10,-20,20,1000,0.-40,40);
    for (int k = 0; k < nEntries; k++) {
        //int nRows = input->GetNRows("EVNT");
        //for(int i = 0; i<nRows ;i++){
            TString category = t->GetCategorization(0);//dejar numero de hit en 0 (primer hit)
            if (category == "electron") { //hacer con arreglo de TH2F
                if(t->Sector(0)==0){
                    h->Fill(t->X(0),t->Y(0));
                }
                if(t->Sector(0)==1){
                    h1->Fill(t->X(0),t->Y(0));
                }
                if(t->Sector(0)==2){
                    h2->Fill(t->X(0),t->Y(0));
                }
                if(t->Sector(0)==3){
                    h3->Fill(t->X(0),t->Y(0));
                }
                if(t->Sector(0)==4){
                    h4->Fill(t->X(0),t->Y(0));
                }
                if(t->Sector(0)==5){
                    h5->Fill(t->X(0),t->Y(0));
                }
            }
        //}
        input->Next();
    }

    TFile *output_m = new TFile("XY.root","RECREATE","X vs Y");
    h->Write();
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h5->Write();
    output_m->Close();
}
