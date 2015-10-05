#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1.h"

TString dataLoc;
TString fDataExt;
TString fSimuExt;
TString elecExt;
TString pionExt;

Int_t nSimuFiles;

Double_t Q2_MIN;
Double_t Q2_MAX;
Double_t XB_MIN;
Double_t XB_MAX;
Double_t NU_MIN;
Double_t NU_MAX;
Double_t ZH_MIN;
Double_t ZH_MAX;
Double_t PT_MIN;
Double_t PT_MAX;
Double_t PHI_MIN;
Double_t PHI_MAX;

Int_t N_Q2;
Int_t N_XB;
Int_t N_NU;
Int_t N_ZH;
Int_t N_PT;
Int_t N_PHI;

Int_t XF_POS = 0;

Double_t delta_Q2;
Double_t delta_XB;
Double_t delta_NU;
Double_t delta_ZH;
Double_t delta_PT;
Double_t delta_PHI;

const Int_t N_METAL = 4;

TFile *plots;
TNtuple *fit_data;

void run_file(const Char_t Metal[], Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX, Double_t PT_MIN, Double_t PT_MAX) {
    TCut P_cut = Form("((T4>-0.55 && P<2.7) || P>2.7)");
    TCut Q2_cut = Form("Q2>%f && Q2<%f", Q2_MIN, Q2_MAX);
    TCut Xb_cut = Form("Xb>%f && Xb<%f", XB_MIN, XB_MAX);
    TCut Zh_cut = Form("Zh>%f && Zh<%f", ZH_MIN, ZH_MAX);
    TCut Pt_cut = Form("Pt>%f && Pt<%f", PT_MIN, PT_MAX);
	TCut Xf_cut;
	if(XF_POS == 0) Xf_cut = Form("Xf<0 || Xf>0");
	if(XF_POS == 1) Xf_cut = Form("Xf>0");
	if(XF_POS == -1) Xf_cut = Form("Xf<0");
	TCut Target_cut;
	if(Metal == "D")
    	Target_cut = Form("TargType==%d", 1);
	else
		Target_cut = Form("TargType==%d", 2);

    TCut cuts = P_cut && Q2_cut && Xb_cut && Zh_cut && Pt_cut && Xf_cut && Target_cut;
    TCut cuts_simul = P_cut && Q2_cut && Xb_cut && Zh_cut&& Pt_cut && Xf_cut;

    TChain *ntuple = new TChain("data_pion");
    ntuple->Add(dataLoc + Metal + fDataExt + pionExt);
    ntuple->Draw(">>list", cuts, "goff");
    ntuple->SetEventList((TEventList*) gDirectory->Get("list"));

    TChain *accept = new TChain("accept_pion");
	for(Int_t q = 0; q < nSimuFiles; q++)
		accept->Add(dataLoc + Metal + std::to_string(q+1) + fSimuExt);
    accept->Draw(">>list_acc",cuts_simul,"goff");
    accept->SetEventList((TEventList*) gDirectory->Get("list_acc"));

    TChain *thrown = new TChain("thrown_pion");
	for(Int_t q = 0; q < nSimuFiles; q++)
		thrown->Add(dataLoc + Metal + std::to_string(q+1) + fSimuExt);
    thrown->Draw(">>list_thr", cuts_simul, "goff");
    thrown->SetEventList((TEventList*)gDirectory->Get("list_thr"));

    ntuple->Draw((const char*)Form("PhiPQ>>htmp_data(%d,%f,%f)", N_PHI, PHI_MIN, PHI_MAX), "", "goff");
    TH1F *htmp_data = (TH1F*)gDirectory->GetList()->FindObject("htmp_data");
    htmp_data->Sumw2();
    accept->Draw((const char*)Form("PhiPQ>>htmp_acc(%d,%f,%f)", N_PHI, PHI_MIN, PHI_MAX), "", "goff");
    TH1F *htmp_acc = (TH1F*)gDirectory->GetList()->FindObject("htmp_acc");
    htmp_acc->Sumw2();
    thrown->Draw((const char*)Form("PhiPQ>>htmp_thr(%d,%f,%f)", N_PHI, PHI_MIN, PHI_MAX), "", "goff");
    TH1F *htmp_thr = (TH1F*)gDirectory->GetList()->FindObject("htmp_thr");
    htmp_thr->Sumw2();
    TH1F *htmp_acc_ratio = new TH1F("htmp_acc_ratio", "", N_PHI, PHI_MIN, PHI_MAX);
    TH1F *htmp_data_corr = new TH1F("htmp_data_corr", "", N_PHI, PHI_MIN, PHI_MAX);
    htmp_acc_ratio->Divide(htmp_acc,htmp_thr,1,1,"B");
    htmp_data_corr->Divide(htmp_data,htmp_acc_ratio,1,1);
    
    for(Int_t ii = 1; ii <= N_PHI; ii++){
        Double_t PhiVal = PHI_MIN + (ii-0.5)*delta_PHI;
        Double_t Val = htmp_data_corr->GetBinContent(ii);
        Double_t Err = htmp_data_corr->GetBinError(ii);
        fit_data->Fill((Q2_MIN+Q2_MAX)/2., (XB_MIN+XB_MAX)/2., (ZH_MIN+ZH_MAX)/2., (PT_MIN+PT_MAX)/2., PhiVal, Val, Err);
    }
    
    delete htmp_acc_ratio;
    delete htmp_data_corr;
    delete ntuple;
    delete accept;
    delete thrown;
}

int main(int argc, char **argv){
	TString Metal;
	TString rootFName;

	if(argc < 25){
		std::cout << "The number of arguments is incorrect" << std::endl;
	}

	Q2_MIN = (Double_t) std::stod(argv[1]);
	Q2_MAX = (Double_t) std::stod(argv[2]);
	N_Q2 = (Int_t) std::stoi(argv[3]);
	XB_MIN = (Double_t) std::stod(argv[4]);
	XB_MAX = (Double_t) std::stod(argv[5]);
	N_XB = (Int_t) std::stoi(argv[6]);
	NU_MIN = (Double_t) std::stod(argv[7]);
	NU_MAX = (Double_t) std::stod(argv[8]);
	N_NU = (Int_t) std::stoi(argv[9]);
	ZH_MIN = (Double_t) std::stod(argv[10]);
	ZH_MAX = (Double_t) std::stod(argv[11]);
	N_ZH = (Int_t) std::stoi(argv[12]); 
	PT_MIN = (Double_t) std::stod(argv[13]);
	PT_MAX = (Double_t) std::stod(argv[14]);
	N_PT = (Int_t) std::stoi(argv[15]);
	PHI_MIN = (Double_t) std::stod(argv[16]);
	PHI_MAX = (Double_t) std::stod(argv[17]);
	N_PHI = (Int_t) std::stoi(argv[18]);

	dataLoc = (TString) argv[19];
	fDataExt = (TString) argv[20];
	fSimuExt = (TString) argv[21];
	nSimuFiles = (Int_t) std::stoi(argv[22]);
	elecExt = (TString) argv[23];
	pionExt = (TString) argv[24];

	if(argc == 26) XF_POS = (Int_t) std::stoi(argv[25]);

	delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
	delta_XB = (XB_MAX-XB_MIN)/N_XB;
	delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
	delta_PT = (PT_MAX-PT_MIN)/N_PT;
	delta_PHI = (PHI_MAX-PHI_MIN)/N_PHI;
	delta_NU = (NU_MAX-NU_MIN)/N_NU;

	for(Int_t met = 0; met < N_METAL; met++){
		if(met == 0) Metal = "C";
		else if(met == 1) Metal = "Fe";
		else if(met == 2) Metal = "Pb";
		else if(met == 3) Metal = "D";
		rootFName = Metal + "_5_dim_dist.root";
		plots = new TFile(rootFName, "RECREATE");
		fit_data = new TNtuple("fit_data", "DATA FOR 5 DIM FIT", "Q2:Xb:Zh:Pt:Phi:Val:Err");
		for(Int_t i = 0; i < N_Q2; i++){
		    for(Int_t j = 0; j < N_XB; j++){
		        for(Int_t k = 0; k < N_ZH; k++){
		            for(Int_t l = 0; l < N_PT; l++){
		                run_file(Metal, Q2_MIN+i*delta_Q2, Q2_MIN+(i+1)*delta_Q2,
		                            XB_MIN+j*delta_XB, XB_MIN+(j+1)*delta_XB,
		                            ZH_MIN+k*delta_ZH, ZH_MIN+(k+1)*delta_ZH,
		                            PT_MIN+l*delta_PT, PT_MIN+(l+1)*delta_PT);
		                plots->cd();
		            }
		        }
		    }
		}
		plots->cd();
		fit_data->Write();
		plots->Close();

		delete plots;
		delete fit_data;
	}

    return 0;
}
