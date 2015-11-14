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
#include "TF1.h"

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

//This should be changed later to pass as a input 
Int_t RCOn = 0;
Double_t PT2_MIN = 0;
Double_t PT2_MAX = 3;
Int_t N_PT2 = 30;
Double_t delta_PT2;
Double_t *v_PT2;

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

Int_t q2i, xbi, zhi;

TH1F *h_data_corr;
TH1F *h_data_corr_RC;
TFile *file;

Float_t *rcFactorC;
Float_t *rcFactorFe;
Float_t *rcFactorPb;
Float_t *rcFactorDC;
Float_t *rcFactorDFe;
Float_t *rcFactorDPb;

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
	elecExt = (TString) argv[22];
	pionExt = (TString) argv[23];
	
	RCOn = (Int_t) std::stoi(argv[24]);
	
	delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
	delta_XB = (XB_MAX-XB_MIN)/N_XB;
	delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
	delta_PT = (PT_MAX-PT_MIN)/N_PT;
	delta_PT2 = (PT2_MAX-PT2_MIN)/N_PT2;
	delta_PHI = (PHI_MAX-PHI_MIN)/N_PHI;
	delta_NU = (NU_MAX-NU_MIN)/N_NU;
	v_Q2 = new Double_t[N_Q2+1];
	v_NU = new Double_t[N_NU+1];
	v_ZH = new Double_t[N_ZH+1];
	v_XB = new Double_t[N_XB+1];
	v_PT = new Double_t[N_PT+1];
	v_PT2 = new Double_t[N_PT2+1];
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
		NTot = N_XB*N_Q2*N_ZH*N_PT2*N_PHI;
		rcFactorC = new Float_t[NTot];
		rcFactorFe = new Float_t[NTot];
		rcFactorPb = new Float_t[NTot];
		rcFactorDC = new Float_t[NTot];
		rcFactorDFe = new Float_t[NTot];
		rcFactorDPb = new Float_t[NTot];
		RCFactors();
	}
	
	file = new TFile("pt2.root","RECREATE");
    
	for(q2i = 0; q2i < N_Q2; q2i++){
		for(xbi = 0; xbi < N_XB; xbi++){
			for(zhi = 0; zhi < N_ZH; zhi++){
				plots(v_Q2[q2i], v_Q2[q2i+1], v_XB[xbi], v_XB[xbi+1], v_ZH[zhi], v_ZH[zhi+1]);
			}
		}
	}
	file->Close();
	
    delete v_Q2;
    delete v_XB;
    delete v_ZH;
    delete v_PT;
    delete v_PT2;
    delete v_PHI;
    
    delete rcFactorC;
    delete rcFactorFe;
    delete rcFactorPb;
    delete rcFactorDC;
    delete rcFactorDFe;
    delete rcFactorDPb;
	
	return 0;
}

void plots(Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX){
	std::cout << "Running bin " << std::endl;
	std::cout << "Q2bin = " << q2i << std::endl;
	std::cout << "Xbbin = " << xbi << std::endl;
	std::cout << "Zhbin = " << zhi << std::endl;
	run_file("C", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("Fe", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("Pb", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("DC", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("DFe", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
	run_file("DPb", Q2_MIN, Q2_MAX, XB_MIN, XB_MAX, ZH_MIN, ZH_MAX);
}

void run_file(TString Metal, Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX){
	h_data_corr = new TH1F("h_data_corr", "", N_PT2, PT2_MIN, PT2_MAX);
	h_data_corr_RC = new TH1F("h_data_corr_RC", "", N_PT2, PT2_MIN, PT2_MAX);

	Int_t index;
	TString solidMetal;
	TCut Target_cut;
	TCut Phi_cut;	
	TCut Q2_cut = Form("Q2>%f && Q2<%f", Q2_MIN, Q2_MAX);
	TCut Xb_cut = Form("Xb>%f && Xb<%f", XB_MIN, XB_MAX);
	TCut Zh_cut = Form("Zh>%f && Zh<%f", ZH_MIN, ZH_MAX);
	Float_t *rcFactor;
	
	if(Metal == "C") rcFactor = rcFactorC;
	else if(Metal == "Fe") rcFactor = rcFactorFe;
	else if(Metal == "Pb") rcFactor = rcFactorPb;
	else if(Metal == "DC") rcFactor = rcFactorDC;
	else if(Metal == "DFe") rcFactor = rcFactorDFe;
	else if(Metal == "DPb") rcFactor = rcFactorDPb;
	else{
		return;
	}
	
	if(Metal(0,1) == "D"){
		Target_cut = Form("TargType==%d", 1);
		solidMetal = Metal(2,2);
		Metal = "D";
	}
	else{
		Target_cut = Form("TargType==%d", 2);
		solidMetal = Metal;
	}
	
	TCut cuts = Q2_cut && Xb_cut && Zh_cut && Target_cut;
	TCut cuts_simul = Q2_cut && Xb_cut && Zh_cut;
	
	TChain *fntuple, *faccept, *fthrown;	

	fntuple = new TChain("data_pion");
	fntuple->Add(dataLoc + solidMetal + fDataExt + pionExt);
	fntuple->Draw("Pt*Pt>>htmp1(150, 0, 3.5)", cuts, "goff");
	
	TH1F *h1 = (TH1F*) gDirectory->Get("htmp1");

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
		return;
	}

	TCut Pt2_cut = Form("Pt*Pt < %f", upper_lim);
	cuts = cuts && Pt2_cut;
	cuts_simul = cuts && Pt2_cut;

	fntuple->Draw(">>list", cuts, "goff");
	fntuple->SetEventList((TEventList*)gDirectory->Get("list"));
	
	faccept = new TChain("accept_pion");
	faccept->Add(dataLoc + Metal + "1" + fSimuExt + pionExt);
	faccept->Add(dataLoc + Metal + "2" + fSimuExt + pionExt);
	faccept->Add(dataLoc + Metal + "3" + fSimuExt + pionExt);
	faccept->Add(dataLoc + Metal + "4" + fSimuExt + pionExt);
	faccept->Draw(">>list_acc", cuts_simul, "goff");
	faccept->SetEventList((TEventList*) gDirectory->Get("list_acc"));

	fthrown = new TChain("thrown_pion");
	fthrown->Add(dataLoc + Metal + "1" + fSimuExt + pionExt);
	fthrown->Add(dataLoc + Metal + "2" + fSimuExt + pionExt);
	fthrown->Add(dataLoc + Metal + "3" + fSimuExt + pionExt);
	fthrown->Add(dataLoc + Metal + "4" + fSimuExt + pionExt);
	fthrown->Draw(">>list_thr",cuts_simul,"goff");
	fthrown->SetEventList((TEventList*)gDirectory->Get("list_thr"));

	for(Int_t l = 0; l < N_PHI; l++){
		Phi_cut = Form("PhiPQ>%f && PhiPQ<%f", PHI_MIN+l*delta_PHI, PHI_MIN+(l+1)*delta_PHI);
		
		fntuple->Draw((const char*) Form("Pt*Pt>>htmp_data(%d,%f,%f)", N_PT2, PT2_MIN, PT2_MAX), Phi_cut, "goff");
		TH1F *htmp_data = (TH1F*) gDirectory->GetList()->FindObject("htmp_data");
		htmp_data->Sumw2();
		
		faccept->Draw((const char*)Form("Pt*Pt>>htmp_acc(%d,%f,%f)", N_PT2, PT2_MIN, PT2_MAX), Phi_cut, "goff");
		TH1F *htmp_acc = (TH1F*) gDirectory->GetList()->FindObject("htmp_acc");
    	htmp_acc->Sumw2();

		fthrown->Draw((const char*)Form("Pt*Pt>>htmp_thr(%d,%f,%f)", N_PT2, PT2_MIN, PT2_MAX), Phi_cut, "goff");
    	TH1F *htmp_thr = (TH1F*) gDirectory->GetList()->FindObject("htmp_thr");
    	htmp_thr->Sumw2();
		TH1F *htmp_acc_ratio = new TH1F("htmp_acc_ratio", "", N_PT2, PT2_MIN, PT2_MAX);
		TH1F *htmp_data_corr = new TH1F("htmp_data_corr", "", N_PT2, PT2_MIN, PT2_MAX);
		TH1F *htmp_data_corr_rc = new TH1F("htmp_data_corr_rc", "", N_PT2, PT2_MIN, PT2_MAX);
		htmp_acc_ratio->Divide(htmp_acc, htmp_thr, 1, 1, "B");
		htmp_data_corr->Divide(htmp_data, htmp_acc_ratio, 1, 1);
		h_data_corr->Add(htmp_data_corr,1);
		if(RCOn){
			for(Int_t m=1; m <= N_PT; m++){
				index = xbi*N_Q2*N_ZH*N_PT*N_PHI+q2i*N_ZH*N_PT*N_PHI+zhi*N_PT*N_PHI+(m-1)*N_PHI+l;
      			htmp_data_corr_rc->SetBinContent(m, (htmp_data_corr->GetBinContent(m))*rcFactor[index]);
      			htmp_data_corr_rc->SetBinError(m, htmp_data_corr->GetBinError(m));
			}
			h_data_corr_RC->Add(htmp_data_corr,1);
		}
		delete htmp_acc_ratio;
		delete htmp_data_corr;
		delete htmp_data_corr_rc;		
	}
	file->cd();
	h_data_corr->Write((const char*) Form("%s_hist_%d%d%d", (const char*)Metal, q2i, xbi, zhi));
	if(RCOn) h_data_corr_RC->Write((const char*) Form("%s_hist_rc_%d%d%d", (const char*)Metal, q2i, xbi, zhi));	
	delete h_data_corr;
	delete h_data_corr_RC;	
	
	delete fntuple;
	delete faccept;
	delete fthrown;
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
