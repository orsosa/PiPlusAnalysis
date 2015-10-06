#include <fstream>
#include <iostream>
#include <string>
#include "TString.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

const Int_t N_METAL = 4;

int main(int argc, char **argv)
{
	TFile *f;
	TFile *newf;
	TNtuple *ntuple;
	TNtuple *newntuple;
	TH1F *h1;
	TF1 *func;
	TString Metal;
	
	Double_t Q2, Xb, Zh, Pt, Phi, Val, Err;
	Double_t A, Ac, Acc;
	Double_t AErr, AcErr, AccErr;
	Double_t sumChi;

	Int_t nentries, empty;

	for(Int_t met = 0; met < N_METAL; met++){
		if(met == 0) Metal = "C";
		else if(met == 1) Metal = "Fe";
		else if(met == 2) Metal = "Pb";
		else if(met == 3) Metal = "D";
		f = new TFile("../Gen5DimData/" + Metal + "_5_dim_dist.root", "UPDATE");
		ntuple = (TNtuple*) f->Get("fit_data");

		nentries = ntuple->GetEntries();

		newntuple = new TNtuple("AAcAcc_data", "AAcAcc_data", "Q2:Xb:Zh:Pt:A:AErr:Ac:AcErr:Acc:AccErr:ChiSQ");
		func = new TF1("fit", "[0]+TMath::Cos(x*TMath::Pi()/180)*[1]+TMath::Cos(2*x*TMath::Pi()/180)*[2]");

		ntuple->SetBranchAddress("Xb", &Xb);
		ntuple->SetBranchAddress("Q2", &Q2);
		ntuple->SetBranchAddress("Xb", &Xb);
		ntuple->SetBranchAddress("Zh", &Zh);
		ntuple->SetBranchAddress("Pt", &Pt);
		ntuple->SetBranchAddress("Phi", &Phi);
		ntuple->SetBranchAddress("Val", &Val);
		ntuple->SetBranchAddress("Err", &Err);

		for(Int_t i = 0; i < nentries; i = i + 12){
			ntuple->GetEntry(i);
			h1 = new TH1F((const char*) Form("PhiDist Q2=%.3f Xb=%.3f Zh=%.3f Pt=%.3f", Q2, Xb, Zh, Pt), 
							(const char*) Form("PhiDist Q2=%.3f Xb=%.3f Zh=%.3f Pt=%.3f", Q2, Xb, Zh, Pt), 12, -180, 180);
			empty = 0;
			for(Int_t j = 1; j <= 12; j++){
				ntuple->GetEntry(i+j-1);
				h1->SetBinContent(j, Val);
				if(h1->GetBinContent(j) == 0)
					empty++;
				h1->SetBinError(j, Err*1.04);
			}
			if(empty != 12){
				h1->Fit(func);
				h1->Write();
				if(func->GetNDF() != 0){
					sumChi = func->GetChisquare();
					A = func->GetParameter(0);
					AErr = func->GetParError(0);
					Ac = func->GetParameter(1);
					AcErr = func->GetParError(1);
					Acc = func->GetParameter(2);
					AccErr = func->GetParError(2);
					newntuple->Fill(Q2, Xb, Zh, Pt, A, AErr, Ac, AcErr, Acc, AccErr, sumChi);
				}
			}
			delete h1;
		}
		f->Close();
		delete f;

		newf = new TFile(Metal + "newphihist.root", "RECREATE");
		newf->cd();
		newntuple->Write();
		newf->Close();
		delete newntuple;
		delete newf;		
	}
	return 0;	
}
