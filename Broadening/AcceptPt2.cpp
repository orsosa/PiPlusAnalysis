#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TAxis.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TCut.h"
#include "TEventList.h"

TString dataLoc;
TString fDataExt;
TString fSimuExt;
TString elecExt;
TString pionExt;

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

Double_t delta_Q2;
Double_t delta_XB;
Double_t delta_NU;
Double_t delta_ZH;
Double_t delta_PT;
Double_t delta_PHI;

Double_t *v_Q2;
Double_t *v_XB;
Double_t *v_NU;
Double_t *v_ZH;
Double_t *v_PT;
Double_t *v_PHI;

Int_t ii, jj, kk;

TCanvas *c1;
TH1F *h_data_corr;
TH1F *h_data_corr_RC;
TFile *plots = new TFile("pt2.root","RECREATE");

Float_t *rcFactorC;
Float_t *rcFactorFe;
Float_t *rcFactorPb;
Float_t *rcFactorDC;
Float_t *rcFactorDFe;
Float_t *rcFactorDPb;

const Double_t kMassProton = 0.938272;

void RCFactors();

void plots(Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX);

void run_file(TString Metal, Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX);

int main(int argc, char **argv){
	for(Int_t i = 0; i < argc; i++){
		std::cout << "Arg " << i << std::endl;
		std::cout << argv[i] << std::endl;
	}
	if(argc < 24){
		std::cout << "The number of arguments is incorrect" << std::endl;
		std::cout << argc << std::endl;
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
	
	delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
	delta_XB = (XB_MAX-XB_MIN)/N_XB;
	delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
	delta_PT = (PT_MAX-PT_MIN)/N_PT;
	delta_PHI = (PHI_MAX-PHI_MIN)/N_PHI;
	delta_NU = (NU_MAX-NU_MIN)/N_NU;
	v_Q2 = new Double_t[N_Q2+1];
	v_NU = new Double_t[N_NU+1];
	v_ZH = new Double_t[N_ZH+1];
	v_XB = new Double_t[N_XB+1];
	v_PT = new Double_t[N_PT+1];
	v_PHI = new Double_t[N_PHI+1];	

    for(Int_t i = 0; i < N_Q2+1; i++){
        if(i == 0) v_Q2[i] = Q2_MIN;
        else v_Q2[i] = v_Q2[i-1] + delta_Q2;
    }
    for(Int_t i = 0; i < N_XB+1; i++){    
        if(i == 0) v_XB[i] = XB_MIN;
        else v_XB[i] = v_XB[i-1] + delta_XB;    
    }
    for(Int_t i = 0; i < N_ZH+1; i++){
        if(i == 0) v_ZH[i] = ZH_MIN;
        else v_ZH[i] = v_ZH[i-1] + delta_ZH;
    }
    for(Int_t i = 0; i < N_PT+1; i++){
        if(i == 0) v_PT[i] = PT_MIN;
        else v_PT[i] = v_PT[i-1] + delta_PT;
    }    
    for(Int_t i = 0; i < N_PHI+1; i++){
        if(i == 0) v_PHI[i] = PHI_MIN;
        else v_PHI[i] = v_PHI[i-1] + delta_PHI;
    }
    
	Int_t NTot;

	if(RCOn){
		NTot = N_XB*N_Q2*N_ZH*N_PT*N_PHI;
		rcFactorC = new Float_t[NTot];
		rcFactorFe = new Float_t[NTot];
		rcFactorPb = new Float_t[NTot];
		rcFactorDC = new Float_t[NTot];
		rcFactorDFe = new Float_t[NTot];
		rcFactorDPb = new Float_t[NTot];
		RCFactors();
	}                 
    
	for(ii = 0; ii < N_Q2; ii++){
		for(jj = 0; jj < N_XB; jj++){
			for(kk = 0; kk < N_ZH; kk++){
				plots(v_Q2[ii], v_Q2[ii+1], v_XB[jj], v_XB[jj+1], v_ZH[kk], v_ZH[kk+1]);
			}
		}
	}
	plots->Close();
	return 0;
}

void plots(Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX){
	cout << "Running bin " << endl;
	cout << "Q2bin = " << ii << endl;
	cout << "Xbbin = " << jj << endl;
	cout << "Zhbin = " << kk << endl;
	h_data_corr = new TH1F("h_data_corr", "", N_PT, PT_MIN, PT_MAX);
	h_data_corr_RC = new TH1F("h_data_corr_RC", "", N_PT, PT_MIN, PT_MAX);
	run_file("C", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("Fe", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("Pb", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("D", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	plots->cd();
	h_data_corr->Write((const char*) Form("hist_%d%d%d", ii, jj, kk));
	h_data_corr_RC->Write((const char*) Form("hist_rc_%d%d%d", ii, jj, kk));
	delete h_data_corr;
	delete h_data_corr_RC;
}

void run_file(TString Metal, Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX){
	TCut Q2_cut = Form("Q2>%f && Q2<%f", Q2_MIN, Q2_MAX);
	TCut Xb_cut = Form("Xb>%f && Xb<%f", XB_MIN, XB_MAX);
	TCut Zh_cut = Form("Zh>%f && Zh<%f", ZH_MIN, ZH_MAX);
	TCut Target_cut;
	if(Metal == "D") Target_cut = Form("TargType==%d", 1);
	else Target_cut = Form("TargType==%d", 2);
	TCut cuts = Q2_cut && Xb_cut && Zh_cut && Target_cut;
	TCut cuts_simul = Q2_cut && Xb_cut && Zh_cut;
	TCut Phi_cut;

	TChain *ntuple = new TChain("data_pion");
	ntuple->Add(dataLoc + Metal + fDataExt + pionExt);
	
	ntuple->Draw("Pt*Pt>>htmp1(150, 0, 3.5)", cuts, "goff");
	TH1F *h1 = gDirectory->Get("htmp1");

	TF1 *func = new TF1("fit", "exp([0]+[1]*x)");
	Double_t upper_lim;
	Double_t mean = h1->GetMean();
	h1->Fit(func, "", "", mean, 3);
	Double_t B = func->GetParameter(0);
	Double_t K = func->GetParameter(1);

	if(fabs(K) > 0.001){
		upper_lim = -(B+1)/K;
	} 
	else{
		return 0;
	}

	TCut Pt2_cut = Form("Pt*Pt < %f", upper_lim);
	cuts  = cuts && Pt2_cut;
	cuts_simul = cuts && Pt2_cut;

	ntuple->Draw(">>list", cuts, "goff");
	ntuple->SetEventList((TEventList*)gDirectory->Get("list"));
	
	TChain *accept = new TChain("accept_pion");
	accept->Add(dataLoc + Metal + "1" + fSimuExt + pionExt);
	accept->Add(dataLoc + Metal + "2" + fSimuExt + pionExt);
	accept->Add(dataLoc + Metal + "3" + fSimuExt + pionExt);
	accept->Add(dataLoc + Metal + "4" + fSimuExt + pionExt);
	accept->Draw(">>list_acc", cuts_simul, "goff");
	accept->SetEventList((TEventList*) gDirectory->Get("list_acc"));

	TChain *thrown = new TChain("thrown_pion");
	thrown->Add(dataLoc + Metal + "1" + fSimuExt + pionExt);
	thrown->Add(dataLoc + Metal + "2" + fSimuExt + pionExt);
	thrown->Add(dataLoc + Metal + "3" + fSimuExt + pionExt);
	thrown->Add(dataLoc + Metal + "4" + fSimuExt + pionExt);
	thrown->Draw(">>list_thr",cuts_simul,"goff");
	thrown->SetEventList((TEventList*)gDirectory->Get("list_thr"));

	for(Int_t l = 0; l < N_PHI; l++){
		Phi_cut = Form("PhiPQ>%f && PhiPQ<%f", PHI_MIN+l*delta_PHI, PHI_MIN+(l+1)*delta_PHI);
		
		ntuple->Draw((const char*) Form("Pt*Pt>>htmp_data(%d,%f,%f)", N_PT, PT_MIN, PT_MAX), Phi_cut, "goff");
		TH1F *htmp_data = (TH1F*) gDirectory->GetList()->FindObject("htmp_data");
		htmp_data->Sumw2();
		
		accept->Draw((const char*)Form("Pt*Pt>>htmp_acc(%d,%f,%f)", N_PT, PT_MIN, PT_MAX), Phi_cut, "goff");
		TH1F *htmp_acc = (TH1F*) gDirectory->GetList()->FindObject("htmp_acc");
    	htmp_acc->Sumw2();

		thrown->Draw((const char*)Form("Pt*Pt>>htmp_thr(%d,%f,%f)", N_PT, PT_MIN, PT_MAX), Phi_cut, "goff");
    	TH1F *htmp_thr = (TH1F*) gDirectory->GetList()->FindObject("htmp_thr");
    	htmp_thr->Sumw2();
		TH1F *htmp_acc_ratio = new TH1F("htmp_acc_ratio", "", N_PT, PT_MIN, PT_MAX);
		TH1F *htmp_data_corr = new TH1F("htmp_data_corr", "", N_PT, PT_MIN, PT_MAX);
		htmp_acc_ratio->Divide(htmp_acc, htmp_thr, 1, 1, "B");
		htmp_data_corr->Divide(htmp_data, htmp_acc_ratio, 1, 1);
		h_data_corr->Add(htmp_data_corr,1);
		for(Int_t m=1; m<=htmp_data_corr->GetXaxis()->GetNbins(); m++) {
      		htmp_data_corr->SetBinContent(m, (htmp_data_corr->GetBinContent(m))*rcFactor[ii][jj][kk][m-1][l]);
      		htmp_data_corr->SetBinError(m, htmp_data_corr->GetBinError(m));
		}
		h_data_corr_RC->Add(htmp_data_corr,1);
		delete htmp_acc_ratio;
		delete htmp_data_corr;		
	}
	delete ntuple;
	delete accept;
	delete thrown;
}

void RCFactors(){
    TString Metal, Aux;
    ifstream inRC;
	Int_t index;
    Int_t Q2i, Xbi, Zhi, Pti, Phii;
    Float_t sigb, sigob, tail1, tail2, facno, fact; 
    for(Int_t met = 0; met < 6; met++){
        if(met == 0) Metal = "C";
        else if(met == 1) Metal = "Fe";
        else if(met == 2) Metal = "Pb";
        else if(met == 3) Metal = "D_C";
		else if(met == 4) Metal = "D_Fe";
		else if(met == 5) Metal = "D_Pb";
        inRC.open("RCFactor" + Metal + ".txt");
        if(inRC.is_open() == 0){
            std::cout << "File RCFactor" << Metal << ".txt not found " << std::endl;
            return;
        }        
        for(Int_t i = 0; i < 11; i++){
            inRC >> Aux;
        }
        for(Int_t i = 0; i < N_XB; i++){
            for(Int_t j = 0; j < N_Q2; j++){
                for(Int_t k = 0; k < N_ZH; k++){
                    for(Int_t l = 0; l < N_PT; l++){
                        for(Int_t m = 0; m < N_PHI; m++){
							index = i*N_Q2*N_ZH*N_PT*N_PHI+j*N_ZH*N_PT*N_PHI+k*N_PT*N_PHI+l*N_PHI+m;
                            inRC >> Q2i >> Xbi >> Zhi >> Pti >> Phii >> sigb >> sigob >> tail1 >> tail2 >> facno >> fact;
                            if(fact == 0){
                                if(met == 0) rcFactorFe[index] = 1;
                                else if(met == 1) rcFactorFe[index] = 1;
                                else if(met == 2) rcFactorPb[index] = 1;
                                else if(met == 3) rcFactorDC[index] = 1;
								else if(met == 4) rcFactorDFe[index] = 1;
								else if(met == 5) rcFactorDPb[index] = 1;
                            }
                            else{
                                if(met == 0) rcFactorFe[index] = fact;
                                else if(met == 1) rcFactorFe[index] = fact;
                                else if(met == 2) rcFactorPb[index] = fact;
                                else if(met == 3) rcFactorDC[index] = fact;
								else if(met == 4) rcFactorDFe[index] = fact;
								else if(met == 5) rcFactorDPb[index] = fact;
                            }
                        }
                    }
                }
            }
        }
        inRC.close();
    }
}
