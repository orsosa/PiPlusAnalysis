{
  gROOT->Reset();
#include <fstream>
  Int_t bin;
  Double_t SigObs, SigBorn, Tail;

  TH1F *hist = new TH1F("hist","",12,-180.,180.);

  std::ifstream in;
  in.open("out.txt");
  while (in >> bin >> SigObs >> SigBorn >> Tail) {
    hist->SetBinContent(bin,SigObs);
  }

}
