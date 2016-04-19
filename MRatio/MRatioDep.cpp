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

TString depVar;

Int_t RCOn;
Int_t XF_POS = 0;

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

std::map<TString, Double_t> N_el;

Float_t *rcFactorC;
Float_t *rcFactorFe;
Float_t *rcFactorPb;
Float_t *rcFactorDC;
Float_t *rcFactorDFe;
Float_t *rcFactorDPb;

void RCFactors();

void runVar(TString depVar);

void runDataVar(TString depVar, TString Metal);

void runSimulVar(TString depVar, TString Metal);

void GetNel(TString Metal, Double_t Xb_min, Double_t Xb_max, Double_t Q2_min, Double_t Q2_max);

int main(int argc, char **argv){
	if(argc < 26){
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
	pionExt = (TString) argv[24];

	depVar = (TString) argv[25];

	RCOn = (Int_t) std::stoi(argv[26]);	

	if(argc == 28) XF_POS = (Int_t) std::stoi(argv[27]);
	
	std::cout << "The following settings are being used" << std::endl;
	std::cout << Q2_MIN << " < Q2 < " << Q2_MAX << ", N_Q2 = " << N_Q2 << std::endl;
	std::cout << XB_MIN << " < Xb < " << XB_MAX << ", N_XB = " << N_XB << std::endl;
	std::cout << ZH_MIN << " < Zh < " << ZH_MAX << ", N_ZH = " << N_ZH << std::endl;	
	std::cout << PT_MIN << " < Pt < " << PT_MAX << ", N_PT = " << N_PT << std::endl;	
	std::cout << PHI_MIN << " < Phi < " << PHI_MAX << ", N_PHI = " << N_PHI << std::endl;	
	std::cout << std::endl;
	std::cout << "Data Directory = " << dataLoc << std::endl;
	std::cout << "Calculating Multiplicity Ratio depending on " << depVar << std::endl;
	if(RCon)
		std::cout << "Applying Radiative Correction: Yes" << std::cout;
	else
		std::cout << "Applying Radiative Correction: No" << std::cout;
	if(XF_POS == 1)
		std::cout << "Applying Feynman Cut: Positive" << std::cout;
	else if(XF_POS == -1) 
		std::cout << "Applying Feynman Cut: Negative" << std::cout;	
		

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

	runVar(depVar);
	
	delete v_Q2;
	delete v_NU;
	delete v_ZH;
	delete v_XB;
	delete v_PT;
	delete v_PHI;
	
	delete rcFactorC;
	delete rcFactorFe;
	delete rcFactorPb;
	delete rcFactorDC;
	delete rcFactorDFe;
	delete rcFactorDPb;
}

void runVar(TString depVar){
	TString var1, var2, var3, var4;
	Int_t N_VAR;
	Int_t N_1, N_2, N_3, N_4;
	Int_t sys;
	Int_t index;
	Int_t q2i, xbi, zhi, pti, phii;

	if(depVar == "Q2"){
		N_VAR = N_Q2;
		N_1 = N_XB; N_2 = N_ZH;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Xb"; var2 = "Zh"; var3 = "Pt"; var4 = "PhiPQ";
	}
	else if(depVar == "Xb"){
		N_VAR = N_XB;
		N_1 = N_Q2; N_2 = N_ZH;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Zh"; var3 = "Pt"; var4 = "PhiPQ";
	}
	else if(depVar == "Zh"){
		N_VAR = N_ZH;
		N_1 = N_Q2; N_2 = N_XB;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Xb"; var3 = "Pt"; var4 = "PhiPQ";
	}
	else if(depVar == "Pt"){
		N_VAR = N_PT;
		N_1 = N_Q2; N_2 = N_XB;	N_3 = N_ZH;	N_4 = N_PHI;
		var1 = "Xb"; var2 = "Xb"; var3 = "Zh"; var4 = "PhiPQ";
	}
	else{
		std::cout << "Bad value of depVar" << std::endl;
		return;
	}

	TFile *fTemp;

    TH1F *hist_temp;
    TAxis *axis_temp;
    Int_t nbins;
    Double_t low;
    Double_t up;
	
    TH1F *hf_C_cor, *hf_C_cor_RC;
    TH1F *hf_D2C_cor, *hf_D2C_cor_RC; 
    TH1F *hf_Fe_cor, *hf_Fe_cor_RC;
    TH1F *hf_D2Fe_cor, *hf_D2Fe_cor_RC;
    TH1F *hf_Pb_cor, *hf_Pb_cor_RC;
	TH1F *hf_D2Pb_cor, *hf_D2Pb_cor_RC;

    TH1F *h_C_l, *h_Fe_l, *h_Pb_l;
    TH1F *h_C_s, *h_Fe_s, *h_Pb_s;
    TH1F *h_C_acc, *h_Fe_acc, *h_Pb_acc, *h_D_acc;
    TH1F *h_C_thr, *h_Fe_thr, *h_Pb_thr, *h_D_thr;

    TH1F *hsC, *hsFe, *hsPb, *hsD2;
    TH1F *hC_cor, *hFe_cor, *hPb_cor;
    TH1F *hC_cor_rc, *hFe_cor_rc, *hPb_cor_rc;
    TH1F *hD2C_cor, *hD2Fe_cor, *hD2Pb_cor;
    TH1F *hD2C_cor_rc, *hD2Fe_cor_rc, *hD2Pb_cor_rc;

    Double_t number_C[N_VAR];
    Double_t number_D2C[N_VAR];
    Double_t number_Fe[N_VAR];
    Double_t number_D2Fe[N_VAR];
    Double_t number_Pb[N_VAR];
    Double_t number_D2Pb[N_VAR];

    TH1F *hf_C, *hf_C_RC;
    TH1F *hf_Fe, *hf_Fe_RC;
    TH1F *hf_Pb, *hf_Pb_RC;

    sys = system("rm -f temp_" + depVar + ".root");
	if(sys){}

    runDataVar(depVar, "C");
    runSimulVar(depVar, "C");
    runDataVar(depVar, "Fe");
    runSimulVar(depVar, "Fe");
    runDataVar(depVar, "Pb");
    runSimulVar(depVar, "Pb");
    runSimulVar(depVar, "D");
	
	std::cout << "Finish Data and Simul Files" << std::endl;

	fTemp = new TFile("temp_" + depVar + ".root", "READ");
    hist_temp = (TH1F*) fTemp->Get("h_C_l0000"); 
    axis_temp = hist_temp->GetXaxis();       
    nbins = axis_temp->GetNbins();
    low = axis_temp->GetBinLowEdge(1);
    up = axis_temp->GetBinUpEdge(nbins); 

    hf_C_cor = new TH1F("hf_C_cor","",nbins,low,up);  
    hf_C_cor_RC = new TH1F("hf_C_cor_RC", "", nbins, low, up); 
    hf_D2C_cor = new TH1F("hf_D2C_cor","",nbins,low,up);
    hf_D2C_cor_RC = new TH1F("hf_D2C_cor_RC", "", nbins, low, up); 
    hf_Fe_cor = new TH1F("hf_Fe_cor","",nbins,low,up);
    hf_Fe_cor_RC = new TH1F("hf_Fe_cor_RC", "", nbins, low, up);    
    hf_D2Fe_cor = new TH1F("hf_D2Fe_cor","",nbins,low,up);
    hf_D2Fe_cor_RC = new TH1F("hf_D2Fe_cor_RC", "", nbins, low, up);        
    hf_Pb_cor = new TH1F("hf_Pb_cor","",nbins,low,up);
    hf_Pb_cor_RC = new TH1F("hf_Pb_cor_RC", "", nbins, low, up);       
    hf_D2Pb_cor = new TH1F("hf_D2Pb_cor","",nbins,low,up);
    hf_D2Pb_cor_RC = new TH1F("hf_D2Pb_cor_RC", "", nbins, low, up);  

    for(Int_t j = 0; j < N_1; j++){
        for(Int_t k = 0; k < N_2; k++){            
            for(Int_t l = 0; l < N_3; l++){
                for(Int_t m = 0; m < N_4; m++){ 
                    h_C_l = (TH1F*) fTemp->Get((const char*)Form("h_C_l%d%d%d%d", j, k, l, m));
                    h_C_l->Sumw2();
                    h_C_s = (TH1F*) fTemp->Get((const char*)Form("h_C_s%d%d%d%d", j, k, l, m));
                    h_C_s->Sumw2();
                    h_Fe_l = (TH1F*) fTemp->Get((const char*)Form("h_Fe_l%d%d%d%d", j, k, l, m));
                    h_Fe_l->Sumw2();
                    h_Fe_s = (TH1F*) fTemp->Get((const char*)Form("h_Fe_s%d%d%d%d", j, k, l, m));
                    h_Fe_s->Sumw2();
                    h_Pb_l = (TH1F*) fTemp->Get((const char*)Form("h_Pb_l%d%d%d%d", j, k, l, m));
                    h_Pb_l->Sumw2();
                    h_Pb_s = (TH1F*) fTemp->Get((const char*)Form("h_Pb_s%d%d%d%d", j, k, l, m));
                    h_Pb_s->Sumw2();
                    h_C_acc = (TH1F*) fTemp->Get((const char*)Form("h_C_acc%d%d%d%d", j, k, l, m));
                    h_C_acc->Sumw2();
                    h_C_thr = (TH1F*) fTemp->Get((const char*)Form("h_C_thr%d%d%d%d", j, k, l, m));
                    h_C_thr->Sumw2();
                    h_Fe_acc = (TH1F*) fTemp->Get((const char*)Form("h_Fe_acc%d%d%d%d", j, k, l, m));
                    h_Fe_acc->Sumw2();
                    h_Fe_thr = (TH1F*) fTemp->Get((const char*)Form("h_Fe_thr%d%d%d%d", j, k, l, m));
                    h_Fe_thr->Sumw2();
                    h_Pb_acc = (TH1F*) fTemp->Get((const char*)Form("h_Pb_acc%d%d%d%d", j, k, l, m));
                    h_Pb_acc->Sumw2();
                    h_Pb_thr = (TH1F*) fTemp->Get((const char*)Form("h_Pb_thr%d%d%d%d", j, k, l, m));
                    h_Pb_thr->Sumw2();
                    h_D_acc = (TH1F*) fTemp->Get((const char*)Form("h_D_acc%d%d%d%d", j, k, l, m));
                    h_D_acc->Sumw2();
                    h_D_thr = (TH1F*) fTemp->Get((const char*)Form("h_D_thr%d%d%d%d", j, k, l, m));
                    h_D_thr->Sumw2();
                    
                    //Applying acceptance correction to the temp root histogram file
                    hsC = new TH1F("hsC", "", nbins,low,up);
                    hC_cor = new TH1F("hC_cor","",nbins,low,up);
                    hC_cor_rc = new TH1F("hC_cor_rc", "", nbins, low, up);
                    hsC->Divide(h_C_acc,h_C_thr,1,1);
                    hC_cor->Divide(h_C_s,hsC,1,1);
                    hf_C_cor->Add(hC_cor,1);
                    
                    hsFe = new TH1F("hsFe","",nbins,low,up);
                    hFe_cor = new TH1F("hFe_cor","",nbins,low,up);
                    hFe_cor_rc = new TH1F("hFe_cor_rc", "", nbins, low, up);
                    hsFe->Divide(h_Fe_acc,h_Fe_thr,1,1);
                    hFe_cor->Divide(h_Fe_s,hsFe,1,1);
                    hf_Fe_cor->Add(hFe_cor,1);
                    
                    hsPb = new TH1F("hsPb","",nbins,low,up);
                    hPb_cor = new TH1F("hPb_cor","",nbins,low,up);
                    hPb_cor_rc = new TH1F("hPb_cor_rc", "", nbins, low, up);
                    hsPb->Divide(h_Pb_acc,h_Pb_thr,1,1);
                    hPb_cor->Divide(h_Pb_s,hsPb,1,1);
                    hf_Pb_cor->Add(hPb_cor,1);
                    
                    hsD2 = new TH1F("hsD2","",nbins,low,up);
                    hD2C_cor = new TH1F("hD2C_cor","",nbins,low,up);
                    hD2Fe_cor = new TH1F("hD2Fe_cor","",nbins,low,up);
                    hD2Pb_cor = new TH1F("hD2Pb_cor","",nbins,low,up);
					hD2C_cor_rc = new TH1F("hD2C_cor_rc", "", nbins, low, up);
                    hD2Fe_cor_rc = new TH1F("hD2Fe_cor_rc", "", nbins, low, up);
                    hD2Pb_cor_rc = new TH1F("hD2Pb_cor_rc", "", nbins, low, up);                    
                    hsD2->Divide(h_D_acc,h_D_thr,1,1);
                    hD2C_cor->Divide(h_C_l,hsD2,1,1);
                    hD2Fe_cor->Divide(h_Fe_l,hsD2,1,1);
                    hD2Pb_cor->Divide(h_Pb_l,hsD2,1,1);
                    hf_D2C_cor->Add(hD2C_cor,1);
                    hf_D2Fe_cor->Add(hD2Fe_cor,1);
                    hf_D2Pb_cor->Add(hD2Pb_cor,1);                    

					if(RCOn){
		                for(Int_t i = 0; i < N_VAR; i++){
							if(depVar == "Q2"){
								q2i = i; xbi = j; zhi = k; pti = l; phii = m;
							}
							else if(depVar == "Xb"){
								xbi = i; q2i = j; zhi = k; pti = l; phii = m;
							}
							else if(depVar == "Zh"){
								zhi = i; q2i = j; xbi = k; pti = l; phii = m;
							}
							else if(depVar == "Pt"){
								pti = i; q2i = j; xbi = k; zhi = l; phii = m;
							}
							else{
								std::cout << "Bad value of depVar" << std::endl;
								return;
							}
							index = xbi*N_Q2*N_ZH*N_PT*N_PHI+q2i*N_ZH*N_PT*N_PHI+zhi*N_PT*N_PHI+pti*N_PHI+phii;
		                    hC_cor_rc->SetBinContent(i+1, hC_cor->GetBinContent(i+1)*rcFactorC[index]);
		                    hC_cor_rc->SetBinError(i+1, hC_cor->GetBinError(i+1));
		                    hFe_cor_rc->SetBinContent(i+1, hFe_cor->GetBinContent(i+1)*rcFactorFe[index]);
	                        hFe_cor_rc->SetBinError(i+1, hFe_cor->GetBinError(i+1));
		                    hPb_cor_rc->SetBinContent(i+1, hPb_cor->GetBinContent(i+1)*rcFactorPb[index]);
		                    hPb_cor_rc->SetBinError(i+1, hPb_cor->GetBinError(i+1));
		                    hD2C_cor_rc->SetBinContent(i+1, hD2C_cor->GetBinContent(i+1)*rcFactorDC[index]);
		                    hD2Fe_cor_rc->SetBinContent(i+1, hD2Fe_cor->GetBinContent(i+1)*rcFactorDFe[index]);
		                    hD2Pb_cor_rc->SetBinContent(i+1, hD2Pb_cor->GetBinContent(i+1)*rcFactorDPb[index]);
		                    hD2C_cor_rc->SetBinError(i+1, hD2C_cor->GetBinError(i+1));
		                    hD2Fe_cor_rc->SetBinError(i+1, hD2Fe_cor->GetBinError(i+1));
		                    hD2Pb_cor_rc->SetBinError(i+1, hD2Pb_cor->GetBinError(i+1));
		                }
					}
					
                    hf_C_cor_RC->Add(hC_cor_rc,1);					
                    hf_Fe_cor_RC->Add(hFe_cor_rc,1);					
                    hf_Pb_cor_RC->Add(hPb_cor_rc,1);					
                    hf_D2C_cor_RC->Add(hD2C_cor_rc,1);
                    hf_D2Fe_cor_RC->Add(hD2Fe_cor_rc,1);
                    hf_D2Pb_cor_RC->Add(hD2Pb_cor_rc,1);					
					
					delete hsC;
					delete hC_cor;
					delete hC_cor_rc;
                    delete hsFe;
                    delete hFe_cor;
                    delete hFe_cor_rc;
                    delete hsPb;
                    delete hPb_cor;
                    delete hPb_cor_rc;
                    delete hsD2;
                    delete hD2C_cor;
                    delete hD2Fe_cor;
                    delete hD2Pb_cor;
                    delete hD2C_cor_rc;      
                    delete hD2Fe_cor_rc;        
                    delete hD2Pb_cor_rc;                      
                }
            }
        }
    }

	if(depVar == "Xb"){
		for(Int_t i = 0; i < N_XB; i++){
		    GetNel("C", v_XB[i], v_XB[i+1], Q2_MIN, Q2_MAX);
		    GetNel("Fe", v_XB[i], v_XB[i+1], Q2_MIN, Q2_MAX);
		    GetNel("Pb", v_XB[i], v_XB[i+1], Q2_MIN, Q2_MAX);
		    GetNel("D", v_XB[i], v_XB[i+1], Q2_MIN, Q2_MAX);    
		    if(N_el["C_acc"] == 0)
		        number_C[i] = 1;
		    else
		        number_C[i] = N_el["C_sld"]*N_el["C_thr"]/N_el["C_acc"];
		            
		    if(N_el["D_acc"] == 0){
		        number_D2C[i] = 1;
		        number_D2Fe[i] = 1;
		        number_D2Pb[i] = 1;
		    } else {
		        number_D2C[i] = N_el["C_liq"]*N_el["D_thr"]/N_el["D_acc"];
		        number_D2Fe[i] = N_el["Fe_liq"]*N_el["D_thr"]/N_el["D_acc"];
		        number_D2Pb[i] = N_el["Pb_liq"]*N_el["D_thr"]/N_el["D_acc"];
		    }
		    if(N_el["Fe_acc"] == 0)
		        number_Fe[i] = 1;
		    else
		        number_Fe[i] = N_el["Fe_sld"]*N_el["Fe_thr"]/N_el["Fe_acc"];
		    
		    if(N_el["Pb_acc"] == 0)
		        number_Pb[i] = 1;
		    else
		        number_Pb[i] = N_el["Pb_sld"]*N_el["Pb_thr"]/N_el["Pb_acc"]; 
		}
	}
	else if(depVar == "Q2"){
		for(Int_t j = 0; j < N_Q2; j++){
		    GetNel("C", XB_MIN, XB_MAX, v_Q2[j], v_Q2[j+1]);
		    GetNel("Fe", XB_MIN, XB_MAX, v_Q2[j], v_Q2[j+1]);
		    GetNel("Pb", XB_MIN, XB_MAX, v_Q2[j], v_Q2[j+1]);
		    GetNel("D", XB_MIN, XB_MAX, v_Q2[j], v_Q2[j+1]);
		    if(N_el["C_acc"] == 0)
		        number_C[j] = 1;
		    else
		        number_C[j] = N_el["C_sld"]*N_el["C_thr"]/N_el["C_acc"];
		            
		    if(N_el["D_acc"] == 0){
		        number_D2C[j] = 1;
		        number_D2Fe[j] = 1;
		        number_D2Pb[j] = 1;
		    } else {
		        number_D2C[j] = N_el["C_liq"]*N_el["D_thr"]/N_el["D_acc"];
		        number_D2Fe[j] = N_el["Fe_liq"]*N_el["D_thr"]/N_el["D_acc"];
		        number_D2Pb[j] = N_el["Pb_liq"]*N_el["D_thr"]/N_el["D_acc"];
		    }
		    if(N_el["Fe_acc"] == 0)
		        number_Fe[j] = 1;
		    else
		        number_Fe[j] = N_el["Fe_sld"]*N_el["Fe_thr"]/N_el["Fe_acc"];
		    
		    if(N_el["Pb_acc"] == 0)
		        number_Pb[j] = 1;
		    else
		        number_Pb[j] = N_el["Pb_sld"]*N_el["Pb_thr"]/N_el["Pb_acc"]; 
		}
	}
            
    hf_C = new TH1F("hf_C", "", nbins, low, up);
    hf_C->Divide(hf_C_cor, hf_D2C_cor, 1, 1);
    delete hf_C_cor;
    delete hf_D2C_cor;
    hf_C_RC = new TH1F("hf_C_RC", "", nbins, low, up);
    hf_C_RC->Divide(hf_C_cor_RC, hf_D2C_cor_RC, 1, 1); 
    delete hf_C_cor_RC;
    delete hf_D2C_cor_RC;
    
    hf_Fe = new TH1F("hf_Fe", "", nbins, low, up);
    hf_Fe->Divide(hf_Fe_cor, hf_D2Fe_cor, 1, 1);
    delete hf_Fe_cor;
    delete hf_D2Fe_cor;
    hf_Fe_RC = new TH1F("hf_Fe_RC", "", nbins, low, up);
    hf_Fe_RC->Divide(hf_Fe_cor_RC, hf_D2Fe_cor_RC, 1, 1);
    delete hf_Fe_cor_RC;
    delete hf_D2Fe_cor_RC;        
        
    hf_Pb = new TH1F("hf_Pb", "", nbins, low, up);
    hf_Pb->Divide(hf_Pb_cor, hf_D2Pb_cor, 1, 1);
    delete hf_Pb_cor;
    delete hf_D2Pb_cor;
    hf_Pb_RC = new TH1F("hf_Pb_RC", "", nbins, low, up);
    hf_Pb_RC->Divide(hf_Pb_cor_RC, hf_D2Pb_cor_RC, 1, 1);    
    delete hf_Pb_cor_RC;
    delete hf_D2Pb_cor_RC;    

	if(depVar == "Xb" || depVar == "Q2"){
		for(Int_t i = 0; i < N_VAR; i++){
			hf_C->SetBinContent(i+1, hf_C->GetBinContent(i+1)*number_D2C[i]/number_C[i]);
			hf_Fe->SetBinContent(i+1, hf_Fe->GetBinContent(i+1)*number_D2Fe[i]/number_Fe[i]);
			hf_Pb->SetBinContent(i+1, hf_Pb->GetBinContent(i+1)*number_D2Pb[i]/number_Pb[i]);		
		}
	}
        
    TString name_C = Form("C_hist");
    TString name_Fe = Form("Fe_hist");
    TString name_Pb = Form("Pb_hist");
    TString name_C_rc = Form("C_rc_hist");
    TString name_Fe_rc = Form("Fe_rc_hist");
    TString name_Pb_rc = Form("Pb_rc_hist");    
    
    TFile *fRatio = new TFile("ratio_" + depVar + ".root","UPDATE"); 
    fRatio->cd();
    hf_C->Write(name_C);
    hf_Fe->Write(name_Fe);
    hf_Pb->Write(name_Pb);
    hf_C_RC->Write(name_C_rc);
    hf_Fe_RC->Write(name_Fe_rc);
    hf_Pb_RC->Write(name_Pb_rc);       
	hf_C->Delete();
	hf_Fe->Delete();
	hf_Pb->Delete();
	hf_C_RC->Delete();
	hf_Fe_RC->Delete();
	hf_Pb_RC->Delete();
    fRatio->Close();  
    
    delete fRatio;
    fTemp->Close();
    delete fTemp;
    sys = system("rm -f temp_" + depVar +".root");
	if(sys){}
    return;                                 
}

void runDataVar(TString depVar, TString Metal){
    TCut liquid = "TargType==1";
    TCut solid = "TargType==2";    
    TCut cuts;
	TCut depVarcut;
	TCut var1cut, var2cut, var3cut, var4cut;
	TCut Xf_cut;
	if(XF_POS == 1) Xf_cut = Form("Xf>0");
	else if(XF_POS == -1) Xf_cut = Form("Xf<0");
	else Xf_cut = Form("Xf<0 || Xf>0");

	TString var1, var2, var3, var4;
	Int_t N_VAR;
	Int_t N_1, N_2, N_3, N_4;
	Double_t depV_MIN, depV_MAX;
	Double_t *lim1, *lim2, *lim3, *lim4;	

	if(depVar == "Q2"){
		N_VAR = N_Q2;
		N_1 = N_XB; N_2 = N_ZH;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Xb"; var2 = "Zh"; var3 = "Pt"; var4 = "PhiPQ";
		depV_MIN = Q2_MIN; depV_MAX = Q2_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_XB; lim2 = v_ZH, lim3 = v_PT; lim4 = v_PHI;
	}
	else if(depVar == "Xb"){
		N_VAR = N_XB;
		N_1 = N_Q2; N_2 = N_ZH;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Zh"; var3 = "Pt"; var4 = "PhiPQ";
		depV_MIN = XB_MIN; depV_MAX = XB_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_Q2; lim2 = v_ZH, lim3 = v_PT; lim4 = v_PHI;
	}
	else if(depVar == "Zh"){
		N_VAR = N_ZH;
		N_1 = N_Q2; N_2 = N_XB;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Xb"; var3 = "Pt"; var4 = "PhiPQ";
		depV_MIN = ZH_MIN; depV_MAX = ZH_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_Q2; lim2 = v_XB, lim3 = v_PT; lim4 = v_PHI;
	}
	else if(depVar == "Pt"){
		N_VAR = N_PT;
		N_1 = N_Q2; N_2 = N_XB;	N_3 = N_ZH;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Xb"; var3 = "Zh"; var4 = "PhiPQ";
		depV_MIN = PT_MIN; depV_MAX = PT_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_Q2; lim2 = v_XB, lim3 = v_ZH; lim4 = v_PHI;
	}
	else{
		std::cout << "Bad value of depVar" << std::endl;
		return;
	}
    
    TH1F *h_l[N_1][N_2][N_3][N_4];
    TH1F *h_s[N_1][N_2][N_3][N_4];
    
    TFile *fPion = new TFile(dataLoc + Metal + fDataExt + pionExt);
    TNtuple *ntuplePion = (TNtuple*) fPion->Get("data_pion");
    ntuplePion->Draw(">>list", depVarcut && Xf_cut, "goff");
    ntuplePion->SetEventList((TEventList*)gDirectory->Get("list"));

    for(Int_t j = 0; j < N_1; j++){
        var1cut = Form("%s>%f && %s<%f", (const char*)var1, lim1[j], (const char*)var1, lim1[j+1]);
        for(Int_t k = 0; k < N_2; k++){
            var2cut= Form("%s>%f && %s<%f", (const char*)var2, lim2[k], (const char*)var2, lim2[k+1]);
            for(Int_t l = 0; l < N_3; l++){
                var3cut = Form("%s>%f && %s<%f", (const char*)var3, lim3[l], (const char*)var3, lim3[l+1]);
                for(Int_t m = 0; m < N_4; m++){
                    std::cout << "Getting Data j = " << j << " k = " << k << " l = " << l << " m = " << m << std::endl;
                    var4cut = Form("%s>%f && %s<%f", (const char*)var4, lim4[m], (const char*)var4, lim4[m+1]);
                    cuts = var1cut && var2cut && var3cut && var4cut;
                    ntuplePion->Draw((const char*) Form("%s>>htmp_l(%d,%f,%f)", (const char*)depVar, N_VAR, depV_MIN, depV_MAX), liquid && cuts, "goff");
                    h_l[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_l");
                    h_l[j][k][l][m]->SetName((const char*) Form("h_%s_l%d%d%d%d", (const char*)Metal, j, k, l, m));            
                    ntuplePion->Draw((const char*) Form("%s>>htmp_s(%d,%f,%f)", (const char*)depVar, N_VAR, depV_MIN, depV_MAX), solid && cuts, "goff");
                    h_s[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_s");
                    h_s[j][k][l][m]->SetName((const char*) Form("h_%s_s%d%d%d%d", (const char*)Metal, j, k, l, m));
                }
            }                 
        }
    }

    TFile *fTemp = new TFile("temp_" + depVar +".root", "UPDATE");
    for(Int_t j = 0; j < N_1; j++){
        for(Int_t k = 0; k < N_2; k++){
            for(Int_t l = 0; l < N_3; l++){
                for(Int_t m = 0; m < N_4; m++){
                    h_l[j][k][l][m]->Write();
                    h_s[j][k][l][m]->Write();
                }
            }        
        }
    }
        
    fTemp->Close();
    delete fTemp;
    fPion->Close();
    delete fPion;              
}

void runSimulVar(TString depVar, TString Metal){
    TCut liquid = "TargType==1";
    TCut solid = "TargType==2";    
    TCut cuts;
	TCut depVarcut;
	TCut var1cut, var2cut, var3cut, var4cut;
	TCut Xf_cut;
	if(XF_POS == 1) Xf_cut = Form("Xf>0");
	else if(XF_POS == -1) Xf_cut = Form("Xf<0");
	else Xf_cut = Form("Xf<0 || Xf>0");	

	TString var1, var2, var3, var4;
	Int_t N_VAR;
	Int_t N_1, N_2, N_3, N_4;
    Int_t nentries;
	Double_t depV_MIN, depV_MAX;
	Double_t *lim1, *lim2, *lim3, *lim4;

	if(depVar == "Q2"){
		N_VAR = N_Q2;
		N_1 = N_XB; N_2 = N_ZH;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Xb"; var2 = "Zh"; var3 = "Pt"; var4 = "PhiPQ";
		depV_MIN = Q2_MIN; depV_MAX = Q2_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_XB; lim2 = v_ZH, lim3 = v_PT; lim4 = v_PHI;
	}
	else if(depVar == "Xb"){
		N_VAR = N_XB;
		N_1 = N_Q2; N_2 = N_ZH;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Zh"; var3 = "Pt"; var4 = "PhiPQ";
		depV_MIN = XB_MIN; depV_MAX = XB_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_Q2; lim2 = v_ZH, lim3 = v_PT; lim4 = v_PHI;
	}
	else if(depVar == "Zh"){
		N_VAR = N_ZH;
		N_1 = N_Q2; N_2 = N_XB;	N_3 = N_PT;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Xb"; var3 = "Pt"; var4 = "PhiPQ";
		depV_MIN = ZH_MIN; depV_MAX = ZH_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_Q2; lim2 = v_XB, lim3 = v_PT; lim4 = v_PHI;
	}
	else if(depVar == "Pt"){
		N_VAR = N_PT;
		N_1 = N_Q2; N_2 = N_XB;	N_3 = N_ZH;	N_4 = N_PHI;
		var1 = "Q2"; var2 = "Xb"; var3 = "Zh"; var4 = "PhiPQ";
		depV_MIN = PT_MIN; depV_MAX = PT_MAX;
		depVarcut = Form("%s>%f && %s<%f", (const char*)depVar, depV_MIN, (const char*)depVar, depV_MAX);
		lim1 = v_Q2; lim2 = v_XB, lim3 = v_ZH; lim4 = v_PHI;
	}
	else{
		std::cout << "Bad value of depVar" << std::endl;
		return;
	}
    
    TH1F *h_acc[N_1][N_2][N_3][N_4];
    TH1F *h_thr[N_1][N_2][N_3][N_4];
    
    TChain *accept = new TChain("accept_pion");
	for(Int_t q = 0; q < nSimuFiles; q++)
		accept->Add(dataLoc + Metal + std::to_string(q+1) + fSimuExt + pionExt);
    nentries = accept->GetEntries();
    accept->SetEstimate(nentries);
    accept->Draw(">>list_accept", depVarcut && Xf_cut, "goff");
    accept->SetEventList((TEventList*) gDirectory->Get("list_accept"));

    TChain *thrown = new TChain("thrown_pion");
	for(Int_t q = 0; q < nSimuFiles; q++)
		thrown->Add(dataLoc + Metal + std::to_string(q+1) + fSimuExt + pionExt);
    nentries = thrown->GetEntries();
    thrown->SetEstimate(nentries);
    thrown->Draw(">>list_thrown", depVarcut && Xf_cut, "goff");
    thrown->SetEventList((TEventList*) gDirectory->Get("list_thrown"));
    
    for(Int_t j = 0; j < N_1; j++){
        var1cut = Form("%s>%f && %s<%f", (const char*)var1, lim1[j], (const char*)var1, lim1[j+1]);
        for(Int_t k = 0; k < N_2; k++){
            var2cut = Form("%s>%f && %s<%f", (const char*)var2, lim2[k], (const char*)var2, lim2[k+1]);
            for(Int_t l = 0; l < N_3; l++){
                var3cut = Form("%s>%f && %s<%f", (const char*)var3, lim3[l], (const char*)var3, lim1[l+1]);
                for(Int_t m = 0; m < N_PHI; m++){
                    std::cout << "Getting Simul j = " << j << " k = " << k << " l = " << l << " m = " << m << std::endl;
                    var4cut = Form("%s>%f && %s<%f", (const char*)var4, lim4[m], (const char*)var4, lim4[m+1]);
                    cuts = var1cut && var2cut && var3cut && var4cut;
                    accept->Draw((const char*) Form("%s>>htmp_acc(%d,%f,%f)", (const char*)depVar, N_VAR, depV_MIN, depV_MAX), cuts,"goff");
                    h_acc[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_acc");
                    h_acc[j][k][l][m]->SetName((const char*) Form("h_%s_acc%d%d%d%d", (const char*)Metal, j, k, l, m));
                    thrown->Draw((const char*) Form("%s>>htmp_thr(%d,%f,%f)", (const char*)depVar, N_VAR, depV_MIN, depV_MAX), cuts,"goff");
                    h_thr[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_thr");
                    h_thr[j][k][l][m]->SetName((const char*) Form("h_%s_thr%d%d%d%d", (const char*)Metal, j, k, l, m));
                }
            }           
        }
    }
     

    TFile *fTemp = new TFile("temp_" + depVar + ".root", "UPDATE");
    for(Int_t j = 0; j < N_1; j++){
        for(Int_t k = 0; k < N_2; k++){
            for(Int_t l = 0; l < N_3; l++){
                for(Int_t m = 0; m < N_4; m++){
                    h_acc[j][k][l][m]->Write();
                    h_thr[j][k][l][m]->Write();
                }
            }          
        }
    }  
     
    fTemp->Close();
    delete fTemp;
    delete accept;
    delete thrown;  
}

void GetNel(TString Metal, Double_t Xb_min, Double_t Xb_max, Double_t Q2_min, Double_t Q2_max){
    TCut Xb_cut = Form("Xb>%f && Xb<%f", Xb_min, Xb_max);
    TCut Q2_cut = Form("Q2>%f && Q2<%f", Q2_min, Q2_max);
    TCut liquid = "TargType==1";
    TCut solid = "TargType==2";
    Int_t nentries;     

    if(Metal != "D"){
        TFile *fElec = new TFile(dataLoc + Metal + fDataExt + elecExt);
        TNtuple *ntupleElec = (TNtuple*) fElec->Get("data_elec");
        
        TEventList *eli;
        ntupleElec->Draw(">>list", Xb_cut && Q2_cut, "goff");
        ntupleElec->SetEventList((TEventList*)gDirectory->Get("list"));
        ntupleElec->Draw(">>eli", liquid, "goff");
        eli = (TEventList*) gDirectory->Get("eli");
        N_el[Metal + "_liq"] = (Double_t) eli->GetN();
        ntupleElec->Draw(">>eli", solid , "goff");
        eli = (TEventList*) gDirectory->Get("eli");
        N_el[Metal + "_sld"] = (Double_t) eli->GetN();  
        fElec->Close();    
        delete fElec;             
    }
    
    TChain *elAccept = new TChain("accept_elec");
    elAccept->Add(dataLoc + Metal + "1" + fSimuExt + elecExt);
    elAccept->Add(dataLoc + Metal + "2" + fSimuExt + elecExt);
    elAccept->Add(dataLoc + Metal + "3" + fSimuExt + elecExt);
    elAccept->Add(dataLoc + Metal + "4" + fSimuExt + elecExt);
    nentries = elAccept->GetEntries();        
    elAccept->SetEstimate(nentries);
    elAccept->Draw(">>list_accepted", Xb_cut && Q2_cut, "goff");
    elAccept->SetEventList((TEventList*) gDirectory->Get("list_accepted"));
    TEventList *elistAccepted = (TEventList*) gDirectory->Get("list_accepted");    
    N_el[Metal + "_acc"]= (Double_t) elistAccepted->GetN();
    
    TChain *elThrown = new TChain("thrown_elec");
    elThrown->Add(dataLoc + Metal + "1" + fSimuExt + elecExt);
    elThrown->Add(dataLoc + Metal + "2" + fSimuExt + elecExt);
    elThrown->Add(dataLoc + Metal + "3" + fSimuExt + elecExt);
    elThrown->Add(dataLoc + Metal + "4" + fSimuExt + elecExt);
    nentries = elThrown->GetEntries();
    elThrown->SetEstimate(nentries);
    elThrown->Draw(">>list_thrown", Xb_cut && Q2_cut, "goff");
    elThrown->SetEventList((TEventList*) gDirectory->Get("list_thrown"));
    TEventList *elistThrown = (TEventList*) gDirectory->Get("list_thrown");
    N_el[Metal + "_thr"] = (Double_t) elistThrown->GetN();        

    delete elAccept;
    delete elThrown;   
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
