#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <map>

TString f_location = "/home/rodrigo/Data/PiData/";
TString fd_ext = "_data.root";
TString fs_ext = "_simul.root";

const Double_t Q2_MIN = 1.;
const Double_t Q2_MAX = 4.;
const Double_t XB_MIN = 0.1;
const Double_t XB_MAX = 0.55;
const Double_t ZH_MIN = 0.;
const Double_t ZH_MAX = 1.;
const Double_t PT_MIN = 0.;
const Double_t PT_MAX = 1;
const Double_t PHI_MIN = -180.;
const Double_t PHI_MAX = 180.;

const Double_t kMassProton = 0.938272;

const Double_t NU_MIN = Q2_MIN/(2*kMassProton*XB_MAX);
const Double_t NU_MAX = Q2_MAX/(2*kMassProton*XB_MIN);

const Int_t N_Q2 = 6;
const Int_t N_XB = 5;
const Int_t N_ZH = 10;
const Int_t N_PT = 30;    
const Int_t N_PHI = 12;
const Int_t N_NU = 3;

const Double_t delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
const Double_t delta_XB = (XB_MAX-XB_MIN)/N_XB;
const Double_t delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
const Double_t delta_PT = (PT_MAX-PT_MIN)/N_PT;
const Double_t delta_PHI = (PHI_MAX-PHI_MIN)/N_PHI;
const Double_t delta_NU = (NU_MAX-NU_MIN)/N_NU;

Double_t *v_Q2 = new Double_t[N_Q2+1];
Double_t *v_NU = new Double_t[N_NU+1];
Double_t *v_ZH = new Double_t[N_ZH+1];
Double_t *v_XB = new Double_t[N_XB+1];
Double_t *v_PT = new Double_t[N_PT+1];

Int_t ii, jj, kk;

TCanvas *c1;
TH1F *h_data_corr;
TH1F *h_data_corr_RC;
TFile *plots = new TFile("pt2.root","RECREATE");

Float_t rcFactor[6][5][10][30][12];

void AcceptPt2Script(){
    Float_t q2i, xbi, zhi, pti, phii, sigb, sigob, tail1, tail2, facno, fact;    
    TString string;
    ifstream in;
    in.open("RCFactorC.txt");
    for(Int_t i=0;i<11;i++){
        in >> string;
    }
    for(Int_t j=0;j<N_XB;j++){
        for(Int_t i=0;i<N_Q2;i++){
            for(Int_t k=0;k<N_ZH;k++){
                for(Int_t l=0;l<N_PT;l++){
                    for(Int_t m=0;m<N_PHI;m++){
                        in >> q2i >> xbi >> zhi >> pti >> phii >> sigb >> sigob >> tail1 >> tail2 >> facno >> fact;
                        if(fact == 0)
                        {
                            rcFactor[i][j][k][l][m] = 1;
                        }
                        else
                        {
                            rcFactor[i][j][k][l][m] = fact;
                        }
                    }
                }
            }        
        }    
    }
    in.close();

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
	plots->cd();
	h_data_corr->Write((const char*) Form("hist_%d%d%d", ii, jj, kk));
	h_data_corr_RC->Write((const char*) Form("hist_rc_%d%d%d", ii, jj, kk));
	delete h_data_corr;
	delete h_data_corr_RC;
}

void run_file(const Char_t Metal[], Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX){
	TCut Q2_cut = Form("Q2>%f && Q2<%f", Q2_MIN, Q2_MAX);
	TCut Xb_cut = Form("Xb>%f && Xb<%f", XB_MIN, XB_MAX);
	TCut Zh_cut = Form("Zh>%f && Zh<%f", ZH_MIN, ZH_MAX);
	TCut Target_cut = Form("TargType==%d", 2);
	TCut cuts = Q2_cut && Xb_cut && Zh_cut && Target_cut;
	TCut cuts_simul = Q2_cut && Xb_cut && Zh_cut;
	TCut Phi_cut;

	TChain *ntuple = new TChain("ntuple_data");
	ntuple->Add(f_location + Metal + fd_ext);
	
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
	
	TChain *accept = new TChain("ntuple_accept");
	accept->Add(f_location + Metal + "1" + fs_ext);
	accept->Add(f_location + Metal + "2" + fs_ext);
	accept->Add(f_location + Metal + "3" + fs_ext);
	accept->Add(f_location + Metal + "4" + fs_ext);	
	accept->Draw(">>list_acc", cuts_simul, "goff");
	accept->SetEventList((TEventList*) gDirectory->Get("list_acc"));

	TChain *thrown = new TChain("ntuple_thrown");
	thrown->Add(f_location + Metal + "1" + fs_ext);
	thrown->Add(f_location + Metal + "2" + fs_ext);
	thrown->Add(f_location + Metal + "3" + fs_ext);
	thrown->Add(f_location + Metal + "4" + fs_ext);
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
