#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <map>

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

std::map<std::string, Double_t> N_el;

Double_t *****rcFactorC;
Double_t *****rcFactorFe;
Double_t *****rcFactorPb;
Double_t *****rcFactorD;

const Double_t kMassProton = 0.938272;

int main(int argc, char **argv){
	if(argc != 25){
		std::cout << "The number of arguments is incorrect" << std::endl;
	}
	
	Q2_MIN = (Double_t) std::stod(argv[1]);
	Q2_MAX = (Double_t) std::stod(argv[2]);
	N_Q2 = (Int_t) std::stoi(argv[3]);
	XB_MIN = (Double_t) std::stod(argv[4]);
	XB_MAX = (Double_t) std::stod(argv[5]);
	N_XB = (Int_t) std::stoi(argv[6]);
	NU_MIN = (Double_t) std::std::stod(arg[7]);
	NU_MAX = (Double_t) std::std::stod(arg[8]);
	N_NU = (Int_t) std::stoi(argv[9]);
	ZH_MIN = (Double_t) std::std::stod(arg[10]);
	ZH_MAX = (Double_t) std::std::stod(arg[11]);
	N_ZH = (Int_t) std::stoi(argv[12]); 
	PT_MIN = (Double_t) std::std::stod(arg[13]);
	PT_MAX = (Double_t) std::std::stod(arg[14]);
	N_PT = (Int_t) std::stoi(argv[15]);
	PHI_MIN = (Double_t) std::std::stod(arg[16]);
	PHI_MAX = (Double_t) std::std::stod(arg[17]);
	N_PHI = (Int_t) std::stoi(argv[18]);
	
	dataLoc = (TString) argv[19];
	fDataExt = (TString) argv[20];
	fSimuExt = (TString) argv[21];
	nSimuFiles = (Int_t) std::stoi(argv[22]);
	elecExt = (TString) argv[23];
	pionExt = (TString) argv[24];

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
	rcFactorC = new Double_t[N_XB][N_Q2][N_ZH][N_PT][N_PHI];
	rcFactorFe = new Double_t[N_XB][N_Q2][N_ZH][N_PT][N_PHI];
	rcFactorPb = new Double_t[N_XB][N_Q2][N_ZH][N_PT][N_PHI];
	rcFactorD = new Double_t[N_XB][N_Q2][N_ZH][N_PT][N_PHI];

    RCFactors();
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
    runXb();
}

void RCFactors(){
    TString string;
    TString metal;
    ifstream inRC;
    Int_t Q2i, Xbi, Zhi, Pti, Phii;
    Float_t sigb, sigob, tail1, tail2, facno, fact; 
    for(Int_t r = 0; r < 4; r++){
        if(r == 0) metal = "C";
        else if(r == 1) metal = "Fe";
        else if(r == 2) metal = "Pb";
        else if(r == 3) metal = "D";
        inRC.open("RCFactor" + metal + ".txt");
        if(inRC.is_open() == 0){
            cout << "File RCFactor" << metal << ".txt not found " << endl;
            return;
        }        
        for(Int_t i = 0; i < 11; i++){
            inRC >> string;
        }
        for(Int_t i = 0; i < N_XB; i++){
            for(Int_t j = 0; j < N_Q2; j++){
                for(Int_t k = 0; k < N_ZH; k++){
                    for(Int_t l = 0; l < N_PT; l++){
                        for(Int_t m = 0; m < N_PHI; m++){
                            inRC >> Q2i >> Xbi >> Zhi >> Pti >> Phii >> sigb >> sigob >> tail1 >> tail2 >> facno >> fact;
                            if(fact == 0){
                                if(r == 0) rcFactorC[i][j][k][l][m] = 1;
                                else if(r == 1) rcFactorFe[i][j][k][l][m] = 1;
                                else if(r == 2) rcFactorPb[i][j][k][l][m] = 1;
                                else if(r == 3) rcFactorD[i][j][k][l][m] = 1;
                            }
                            else{
                                if(r == 0) rcFactorC[i][j][k][l][m] = fact;
                                else if(r == 1) rcFactorFe[i][j][k][l][m] = fact;
                                else if(r == 2) rcFactorPb[i][j][k][l][m] = fact;
                                else if(r == 3) rcFactorD[i][j][k][l][m] = fact;
                            }
                        }
                    }
                }
            }
        }
        inRC.close();
    }
}

void runXb(){    
    TFile *fTemp;
    
    TH1F *hist_temp;
    TAxis *axis_temp;
    Int_t nbins;
    Double_t low;
    Double_t up;
    
    TH1F *hf_C_cor;
    TH1F *hf_C_cor_RC;
    TH1F *hf_D2C_cor;
    TH1F *hf_D2C_cor_RC;  
    TH1F *hf_Fe_cor;
    TH1F *hf_Fe_cor_RC;
    TH1F *hf_D2Fe_cor;
    TH1F *hf_D2Fe_cor_RC;  
    TH1F *hf_Pb_cor;
    TH1F *hf_Pb_cor_RC;      
    TH1F *hf_D2Pb_cor;
    TH1F *hf_D2Pb_cor_RC;      
    
    TH1F *h_C_l; TH1F *h_Fe_l; TH1F *h_Pb_l;
    TH1F *h_C_s; TH1F *h_Fe_s; TH1F *h_Pb_s;
    TH1F *h_C_acc; TH1F *h_Fe_acc; TH1F *h_Pb_acc; TH1F *h_D_acc;
    TH1F *h_C_thr; TH1F *h_Fe_thr; TH1F *h_Pb_thr; TH1F *h_D_thr;
   
    TH1F *hsC; TH1F *hsFe; TH1F *hsPb; TH1F *hsD2;
    TH1F *hC_cor; TH1F *hFe_cor; TH1F *hPb_cor;
    TH1F *hC_cor_rc; TH1F *hFe_cor_rc; TH1F *hPb_cor_rc;
    TH1F *hD2C_cor; TH1F *hD2Fe_cor; TH1F *hD2Pb_cor;
    TH1F *hD2C_cor_rc; TH1F *hD2Fe_cor_rc; TH1F *hD2Pb_cor;
    
    Double_t number_C[N_XB];
    Double_t number_D2C[N_XB];
    Double_t number_Fe[N_XB];
    Double_t number_D2Fe[N_XB];
    Double_t number_Pb[N_XB];
    Double_t number_D2Pb[N_XB];
    
    TH1F *hf_C;
    TH1F *hf_C_RC;
    TH1F *hf_Fe;
    TH1F *hf_Fe_RC;
    TH1F *hf_Pb;
    TH1F *hf_Pb_RC; 
    
    TH1F *hff_C;
    TH1F *hff_C_RC;
    TH1F *hff_Fe;
    TH1F *hff_Fe_RC;
    TH1F *hff_Pb;
    TH1F *hff_Pb_RC;

    system("rm -f temp_Xb.root");        

    runDataXb("C");
    runSimulXb("C");
    runDataXb("Fe");
    runSimulXb("Fe");
    runDataXb("Pb");
    runSimulXb("Pb");
    runSimulXb("D");
            
    cout << "Finish Data and Simul Files" << endl;
            
    fTemp = new TFile("temp_Xb.root", "READ");    
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
    
    for(Int_t j = 0; j < N_Q2; j++){
        for(Int_t k = 0; k < N_ZH; k++){            
            for(Int_t l = 0; l < N_PT; l++){
                for(Int_t m = 0; m < N_PHI; m++){ 
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
                    for(Int_t i = 0; i < N_XB; i++){
                        hC_cor_rc->SetBinContent(i+1, hC_cor->GetBinContent(i+1)*rcFactorC[i][j][k][l][m]);
                        hC_cor_rc->SetBinError(i+1, hC_cor->GetBinError(i+1));
                    }
                    hf_C_cor->Add(hC_cor,1);
                    hf_C_cor_RC->Add(hC_cor_rc,1);              
                    delete hsC;
                    delete hC_cor;
                    delete hC_cor_rc;    
                    
                    hsFe = new TH1F("hsFe","",nbins,low,up);
                    hFe_cor = new TH1F("hFe_cor","",nbins,low,up);
                    hFe_cor_rc = new TH1F("hFe_cor_rc", "", nbins, low, up);
                    hsFe->Divide(h_Fe_acc,h_Fe_thr,1,1);
                    hFe_cor->Divide(h_Fe_s,hsFe,1,1);
                    for(Int_t i = 0; i < N_XB; i++){
                        hFe_cor_rc->SetBinContent(i+1, hFe_cor->GetBinContent(i+1)*rcFactorFe[i][j][k][l][m]);
                        hFe_cor_rc->SetBinError(i+1, hFe_cor->GetBinError(i+1));
                    }        
                    hf_Fe_cor->Add(hFe_cor,1);
                    hf_Fe_cor_RC->Add(hFe_cor_rc,1);
                    delete hsFe;
                    delete hFe_cor;
                    delete hFe_cor_rc;
                    
                    hsPb = new TH1F("hsPb","",nbins,low,up);
                    hPb_cor = new TH1F("hPb_cor","",nbins,low,up);
                    hPb_cor_rc = new TH1F("hPb_cor_rc", "", nbins, low, up);
                    hsPb->Divide(h_Pb_acc,h_Pb_thr,1,1);
                    hPb_cor->Divide(h_Pb_s,hsPb,1,1);
                    for(Int_t i = 0; i < N_XB; i++){
                        hPb_cor_rc->SetBinContent(i+1, hPb_cor->GetBinContent(i+1)*rcFactorPb[i][j][k][l][m]);
                        hPb_cor_rc->SetBinError(i+1, hPb_cor->GetBinError(i+1));
                    }          
                    hf_Pb_cor->Add(hPb_cor,1);
                    hf_Pb_cor_RC->Add(hPb_cor_rc,1);
                    delete hsPb;
                    delete hPb_cor;
                    delete hPb_cor_rc;  
                    
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
                    for(Int_t i = 0; i < N_XB; i++){
                        hD2C_cor_rc->SetBinContent(i+1, hD2C_cor->GetBinContent(i+1)*rcFactorD[i][j][k][l][m]);
                        hD2Fe_cor_rc->SetBinContent(i+1, hD2Fe_cor->GetBinContent(i+1)*rcFactorD[i][j][k][l][m]);
                        hD2Pb_cor_rc->SetBinContent(i+1, hD2Pb_cor->GetBinContent(i+1)*rcFactorD[i][j][k][l][m]);
                        hD2C_cor_rc->SetBinError(i+1, hD2C_cor->GetBinError(i+1));
                        hD2Fe_cor_rc->SetBinError(i+1, hD2Fe_cor->GetBinError(i+1));
                        hD2Pb_cor_rc->SetBinError(i+1, hD2Pb_cor->GetBinError(i+1));
                    }        
                    hf_D2C_cor->Add(hD2C_cor,1);
                    hf_D2Fe_cor->Add(hD2Fe_cor,1);
                    hf_D2Pb_cor->Add(hD2Pb_cor,1);
                    hf_D2C_cor_RC->Add(hD2C_cor_rc,1);
                    hf_D2Fe_cor_RC->Add(hD2Fe_cor_rc,1);
                    hf_D2Pb_cor_RC->Add(hD2Pb_cor_rc,1);
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
    
    for(Int_t i = 0; i < N_XB; i++){
        hf_C->SetBinContent(i+1, hf_C->GetBinContent(i+1)*number_D2C[i]/number_C[i]);
        hf_C_RC->SetBinContent(i+1, hf_C_RC->GetBinContent(i+1)*number_D2C[i]/number_C[i]);
        hf_Fe->SetBinContent(i+1, hf_Fe->GetBinContent(i+1)*number_D2Fe[i]/number_Fe[i]);
        hf_Fe_RC->SetBinContent(i+1, hf_Fe_RC->GetBinContent(i+1)*number_D2Fe[i]/number_Fe[i]);
        hf_Pb->SetBinContent(i+1, hf_Pb->GetBinContent(i+1)*number_D2Pb[i]/number_Pb[i]);
        hf_Pb_RC->SetBinContent(i+1, hf_Pb_RC->GetBinContent(i+1)*number_D2Pb[i]/number_Pb[i]);
    }   
    
    TString name_C = Form("C_hist");
    TString name_Fe = Form("Fe_hist");
    TString name_Pb = Form("Pb_hist");
    TString name_C_rc = Form("C_rc_hist");
    TString name_Fe_rc = Form("Fe_rc_hist");
    TString name_Pb_rc = Form("Pb_rc_hist");
    
    TFile *fRatio = new TFile("ratio_Xb.root","UPDATE"); 
    fRatio->cd();
    hf_C->Write(name_C);
    hf_Fe->Write(name_Fe);
    hf_Pb->Write(name_Pb);
    hf_C_RC->Write(name_C_rc);
    hf_Fe_RC->Write(name_Fe_rc);
    hf_Pb_RC->Write(name_Pb_rc);   
    fRatio->Close();  
    
    delete fRatio;
    delete hf_C;
    delete hf_Fe;
    delete hf_Pb;
    delete hf_C_RC;
    delete hf_Fe_RC;
    delete hf_Pb_RC;
    fTemp->Close();
    delete fTemp;
    system("rm -f temp_Xb.root");   
    return;                                 
}

void runDataXb(const Char_t Metal[]){
    Double_t Zh_min = 0.4;
    Double_t Zh_max = 0.8;
    TCut Xb_cut = Form("Xb>%f && Xb<%f", XB_MIN, XB_MAX);
    TCut Q2_cut;
    TCut Zh_cut = Form("Zh>%f && Zh<%f", Zh_min, Zh_max);
    TCut Pt_cut;
    TCut Phi_cut;
    TCut liquid = "TargType==1";
    TCut solid = "TargType==2";    
    TCut cuts;
    
    TH1F *h_l[N_Q2][N_ZH][N_PT][N_PHI];
    TH1F *h_s[N_Q2][N_ZH][N_PT][N_PHI];
    
    TFile *fPion = new TFile(dataLoc + Metal + fDataExt + pionExt);
    TNtuple *ntuplePion = (TNtuple*) fPion->Get("data_pion");
    ntuplePion->Draw(">>list", Zh_cut && Xb_cut, "goff");
    ntuplePion->SetEventList((TEventList*)gDirectory->Get("list"));

    for(Int_t j = 0; j < N_Q2; j++){
        Q2_cut = Form("Q2>%f && Q2<%f", v_Q2[j], v_Q2[j+1]);
        for(Int_t k = 0; k < N_ZH; k++){
            Zh_cut = Form("Zh>%f && Zh<%f", v_ZH[k], v_ZH[k+1]);
            for(Int_t l = 0; l < N_PT; l++){
                Pt_cut = Form("Pt>%f && Pt<%f", v_PT[l], v_PT[l+1]);
                for(Int_t m = 0; m < N_PHI; m++){
                    cout << "Getting Data j = " << j << " k = " << k << " l = " << l << " m = " << m << endl;
                    Phi_cut = Form("PhiPQ>%f && PhiPQ<%f", v_PHI[m], v_PHI[m+1]);
                    cuts = Q2_cut && Zh_cut && Pt_cut && Phi_cut;
                    ntuplePion->Draw((const char*) Form("Xb>>htmp_l(%d,%f,%f)", N_XB, XB_MIN, XB_MAX), liquid && cuts, "goff");
                    h_l[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_l");
                    h_l[j][k][l][m]->SetName((const char*) Form("h_%s_l%d%d%d%d", Metal, j, k, l, m));            
                    ntuplePion->Draw((const char*) Form("Xb>>htmp_s(%d,%f,%f)", N_XB, XB_MIN, XB_MAX), solid && cuts, "goff");
                    h_s[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_s");
                    h_s[j][k][l][m]->SetName((const char*) Form("h_%s_s%d%d%d%d", Metal, j, k, l, m));
                }
            }                 
        }
    }

    TFile *fTemp = new TFile("temp_Xb.root", "UPDATE");
    for(Int_t j = 0; j < N_Q2; j++){
        for(Int_t k = 0; k < N_ZH; k++){
            for(Int_t l = 0; l < N_PT; l++){
                for(Int_t m = 0; m < N_PHI; m++){
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

void runSimulXb(const Char_t Metal[]){
    Double_t Zh_min = 0.4;
    Double_t Zh_max = 0.8;
    TCut Xb_cut = Form("Xb>%f && Xb<%f", XB_MIN, XB_MAX);
    TCut Q2_cut;
    TCut Zh_cut = Form("Zh>%f && Zh<%f", Zh_min, Zh_max);
    TCut Pt_cut;
    TCut Phi_cut;   
    TCut liquid = "TargType==1";
    TCut solid = "TargType==2";    
    TCut cuts;
    
    Int_t nentries;   
    
    TH1F *h_acc[N_Q2][N_ZH][N_PT][N_PHI];
    TH1F *h_thr[N_Q2][N_ZH][N_PT][N_PHI];
    
    TChain *accept = new TChain("accept_pion");
	for(Int_t q = 0; q < nSimuFiles; q++){
		accept->Add(dataLoc + Metal + std::itoa(q+1) + fSimuExt + pionExt);
	}
    nentries = accept->GetEntries();
    accept->SetEstimate(nentries);
    accept->Draw(">>list_accept", Zh_cut && Xb_cut, "goff");
    accept->SetEventList((TEventList*) gDirectory->Get("list_accept"));

    TChain *thrown = new TChain("thrown_pion");
    thrown->Add(dataLoc + Metal + "1" + fSimuExt + pionExt);
    thrown->Add(dataLoc + Metal + "2" + fSimuExt + pionExt);
    thrown->Add(dataLoc + Metal + "3" + fSimuExt + pionExt);
    thrown->Add(dataLoc + Metal + "4" + fSimuExt + pionExt);         
    nentries = thrown->GetEntries();
    thrown->SetEstimate(nentries);
    thrown->Draw(">>list_thrown", Zh_cut && Xb_cut, "goff");
    thrown->SetEventList((TEventList*) gDirectory->Get("list_thrown"));
    
    for(Int_t j = 0; j < N_Q2; j++){
        Q2_cut = Form("Q2>%f && Q2<%f", v_Q2[j], v_Q2[j+1]);
        for(Int_t k = 0; k < N_ZH; k++){
            Zh_cut = Form("Zh>%f && Zh<%f", v_ZH[k], v_ZH[k+1]);
            for(Int_t l = 0; l < N_PT; l++){
                Pt_cut = Form("Pt>%f && Pt<%f", v_PT[l], v_PT[l+1]);
                for(Int_t m = 0; m < N_PHI; m++){
                    cout << "Getting Simul j = " << j << " k = " << k << " l = " << l << " m = " << m << endl;
                    Phi_cut = Form("PhiPQ>%f && PhiPQ<%f", v_PHI[m], v_PHI[m+1]);
                    cuts = Q2_cut && Zh_cut && Pt_cut && Phi_cut;                    
                    accept->Draw((const char*) Form("Xb>>htmp_acc(%d,%f,%f)", N_XB, XB_MIN, XB_MAX), cuts,"goff");
                    h_acc[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_acc");
                    h_acc[j][k][l][m]->SetName((const char*) Form("h_%s_acc%d%d%d%d", Metal, j, k, l, m));
                    thrown->Draw((const char*) Form("Xb>>htmp_thr(%d,%f,%f)", N_XB, XB_MIN, XB_MAX), cuts,"goff");
                    h_thr[j][k][l][m] = (TH1F*) gDirectory->GetList()->FindObject("htmp_thr");
                    h_thr[j][k][l][m]->SetName((const char*) Form("h_%s_thr%d%d%d%d", Metal, j, k, l, m));
                }
            }           
        }
    }
     

    TFile *fTemp = new TFile("temp_Xb.root", "UPDATE");
    for(Int_t j = 0; j < N_Q2; j++){
        for(Int_t k = 0; k < N_ZH; k++){
            for(Int_t l = 0; l < N_PT; l++){
                for(Int_t m = 0; m < N_PHI; m++){
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

void GetNel(const Char_t Metal[], Double_t Xb_min, Double_t Xb_max, Double_t Q2_min, Double_t Q2_max){
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
        N_el[string(Metal)+string("_liq")] = (Double_t) eli->GetN();
        ntupleElec->Draw(">>eli", solid , "goff");
        eli = (TEventList*) gDirectory->Get("eli");
        N_el[string(Metal)+string("_sld")] = (Double_t) eli->GetN();  
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
    N_el[string(Metal)+string("_acc")]= (Double_t) elistAccepted->GetN();
    
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
    N_el[string(Metal)+string("_thr")] = (Double_t) elistThrown->GetN();        

    delete elAccept;
    delete elThrown;   
}
