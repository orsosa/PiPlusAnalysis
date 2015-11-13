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

	TFile *f = new TFile("pt2.root", "READ");

	Double_t mean_C, mean_D_C, mean_Fe, mean_D_Fe, mean_Pb, mean_D_Pb;
	Double_t mean_C_err, mean_D_C_err, mean_Fe_err, mean_D_Fe_err, mean_Pb_err, mean_D_Pb_err;
	Double_t pt2b_C, pt2b_Fe, pt2b_Pb;
	Double_t pt2b_C_err, pt2b_Fe_err, pt2b_Pb_err;

	Double_t mean_C_rc1, mean_D_C_rc1, mean_Fe_rc1, mean_D_Fe_rc1, mean_Pb_rc1, mean_D_Pb_rc1;
	Double_t mean_C_err_rc1, mean_D_C_err_rc1, mean_Fe_err_rc1, mean_D_Fe_err_rc1, mean_Pb_err_rc1, mean_D_Pb_err_rc1;
	Double_t pt2b_C_rc1, pt2b_Fe_rc1, pt2b_Pb_rc1;
	Double_t pt2b_C_err_rc1, pt2b_Fe_err_rc1, pt2b_Pb_err_rc1;

	Double_t mean_C_rc2, mean_D_C_rc2, mean_Fe_rc2, mean_D_Fe_rc2, mean_Pb_rc2, mean_D_Pb_rc2;
	Double_t mean_C_err_rc2, mean_D_C_err_rc2, mean_Fe_err_rc2, mean_D_Fe_err_rc2, mean_Pb_err_rc2, mean_D_Pb_err_rc2;
	Double_t pt2b_C_rc2, pt2b_Fe_rc2, pt2b_Pb_rc2;
	Double_t pt2b_C_err_rc2, pt2b_Fe_err_rc2, pt2b_Pb_err_rc2;

	ofstream out;
	out.open("Broadening.txt");
	out << "Transverse momentum broadening vs A1/3" << endl;
	out << "Zh\t\t\tdPt2 C\t\t\tErr C\t\t\tdPt2 Fe\t\t\tErr Fe\t\t\tdPt2 Pb\t\t\tErr Pb" << endl;

	for(q2i = 0; q2i < N_Q2; q2i++){
		for(xbi = 0; xbi < N_XB; xbi++){
			for(zhi = 0; zhi < N_ZH; zhi++){
				TH1F *h_C = (TH1F*) (f->Get((const char*)Form("hist_%d%d%d",i1,i2,i3)));
				mean_C = h_C->GetMean();
				mean_C_err = h_C->GetMeanError();
				TH1F *h_C_rc1 = dynamic_cast<TH1F*>(f_C->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_C_rc1 = h_C_rc1->GetMean();
				mean_C_err_rc1 = h_C_rc1->GetMeanError();
				f_C_rc->cd();
				TH1F *h_C_rc2 = dynamic_cast<TH1F*>(f_C_rc->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_C_rc2 = h_C_rc2->GetMean();
				mean_C_err_rc2 = h_C_rc2->GetMeanError();

				// Get D_C Statistic
				f_D_C->cd();
				TH1F *h_D_C = dynamic_cast<TH1F*>(f_D_C->Get((const char*)Form("hist_%d%d%d",i1,i2,i3)));
				mean_D_C = h_D_C->GetMean();
				mean_D_C_err = h_D_C->GetMeanError();
				TH1F *h_D_C_rc1 = dynamic_cast<TH1F*>(f_D_C->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_D_C_rc1 = h_D_C_rc1->GetMean();
				mean_D_C_err_rc1 = h_D_C_rc1->GetMeanError();
				f_D_C_rc->cd();
				TH1F *h_D_C_rc2 = dynamic_cast<TH1F*>(f_D_C_rc->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_D_C_rc2 = h_D_C_rc2->GetMean();
				mean_D_C_err_rc2 = h_D_C_rc2->GetMeanError();

				// Get Fe Statistic
				f_Fe->cd();
				TH1F *h_Fe = dynamic_cast<TH1F*>(f_Fe->Get((const char*)Form("hist_%d%d%d",i1,i2,i3)));
				mean_Fe = h_Fe->GetMean();
				mean_Fe_err = h_Fe->GetMeanError();
				TH1F *h_Fe_rc1 = dynamic_cast<TH1F*>(f_Fe->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_Fe_rc1 = h_Fe_rc1->GetMean();
				mean_Fe_err_rc1 = h_Fe_rc1->GetMeanError();
				f_Fe_rc->cd();
				TH1F *h_Fe_rc2 = dynamic_cast<TH1F*>(f_Fe_rc->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_Fe_rc2 = h_Fe_rc2->GetMean();
				mean_Fe_err_rc2 = h_Fe_rc2->GetMeanError();

				// Get D_Fe Statistic
				f_D_Fe->cd();
				TH1F *h_D_Fe = dynamic_cast<TH1F*>(f_D_Fe->Get((const char*)Form("hist_%d%d%d",i1,i2,i3)));
				mean_D_Fe = h_D_Fe->GetMean();
				mean_D_Fe_err = h_D_Fe->GetMeanError();
				TH1F *h_D_Fe_rc1 = dynamic_cast<TH1F*>(f_D_Fe->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_D_Fe_rc1 = h_D_Fe_rc1->GetMean();
				mean_D_Fe_err_rc1 = h_D_Fe_rc1->GetMeanError();
				f_D_Fe_rc->cd();
				TH1F *h_D_Fe_rc2 = dynamic_cast<TH1F*>(f_D_Fe_rc->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_D_Fe_rc2 = h_D_Fe_rc2->GetMean();
				mean_D_Fe_err_rc2 = h_D_Fe_rc2->GetMeanError();

				// Get Pb Statistic
				f_Pb->cd();
				TH1F *h_Pb = dynamic_cast<TH1F*>(f_Pb->Get((const char*)Form("hist_%d%d%d",i1,i2,i3)));
				mean_Pb = h_Pb->GetMean();
				mean_Pb_err = h_Pb->GetMeanError();
				TH1F *h_Pb_rc1 = dynamic_cast<TH1F*>(f_Pb->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_Pb_rc1 = h_Pb_rc1->GetMean();
				mean_Pb_err_rc1 = h_Pb_rc1->GetMeanError();
				f_Pb_rc->cd();
				TH1F *h_Pb_rc2 = dynamic_cast<TH1F*>(f_Pb_rc->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_Pb_rc2 = h_Pb_rc2->GetMean();
				mean_Pb_err_rc2 = h_Pb_rc2->GetMeanError();

				// Get D_Pb Statistic
				f_D_Pb->cd();
				TH1F *h_D_Pb = dynamic_cast<TH1F*>(f_D_Pb->Get((const char*)Form("hist_%d%d%d",i1,i2,i3)));
				mean_D_Pb = h_D_Pb->GetMean();
				mean_D_Pb_err = h_D_Pb->GetMeanError();
				TH1F *h_D_Pb_rc1 = dynamic_cast<TH1F*>(f_D_Pb->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_D_Pb_rc1 = h_D_Pb_rc1->GetMean();
				mean_D_Pb_err_rc1 = h_D_Pb_rc1->GetMeanError();
				f_D_Pb_rc->cd();
				TH1F *h_D_Pb_rc2 = dynamic_cast<TH1F*>(f_D_Pb_rc->Get((const char*)Form("hist_rc_%d%d%d",i1,i2,i3)));
				mean_D_Pb_rc2 = h_D_Pb_rc2->GetMean();
				mean_D_Pb_err_rc2 = h_D_Pb_rc2->GetMeanError();

				// Create canvas and pad for plots
				TCanvas *MyC = new TCanvas("MyC","transverse momuntum broadening of leading pions",6,21,1117,499);
				MyC->SetFillColor(42);
				MyC->SetGrid();
				MyC->GetFrame()->SetFillColor(21);
				MyC->GetFrame()->SetBorderSize(12);
				MyC->Update();
				TPad *S = new TPad("S", "transverse momuntum broadening of leading pions",0,0.160072,1,1);
				S->Draw();
				S->cd();
				S->SetFillColor(42);
				S->SetGrid();
				S->GetFrame()->SetFillColor(21);
				S->GetFrame()->SetBorderSize(12);
				S->Update();

				// Calculate Broadening for RC1
				pt2b_C_rc1 = mean_C_rc1 - mean_D_C_rc1;
				pt2b_Fe_rc1 = mean_Fe_rc1 - mean_D_Fe_rc1;
				pt2b_Pb_rc1 = mean_Pb_rc1 - mean_D_Pb_rc1;
				pt2b_C_err_rc1 = TMath::Sqrt(mean_C_err_rc1*mean_C_err_rc1 + mean_D_C_err_rc1*mean_D_C_err_rc1);
				pt2b_Fe_err_rc1 = TMath::Sqrt(mean_Fe_err_rc1*mean_Fe_err_rc1 + mean_D_Fe_err_rc1*mean_D_Fe_err_rc1);
				pt2b_Pb_err_rc1 = TMath::Sqrt(mean_Pb_err_rc1*mean_Pb_err_rc1 + mean_D_Pb_err_rc1*mean_D_Pb_err_rc1);
			  
				Double_t x25_rc1[3]  = {pow(12,0.333), pow(56,0.333), pow(208,0.333)}; 
				Double_t y25_rc1[3]  = {pt2b_C_rc1, pt2b_Fe_rc1, pt2b_Pb_rc1};
				Double_t ex25_rc1[3] = {0, 0, 0};
				Double_t ey25_rc1[3] = {pt2b_C_err_rc1, pt2b_Fe_err_rc1, pt2b_Pb_err_rc1}; 

				TGraphErrors *gr25_rc1 = new TGraphErrors(3,x25_rc1,y25_rc1,ex25_rc1,ey25_rc1); 

				// Calculate Broadening for RC2
				pt2b_C_rc2 = mean_C_rc2 - mean_D_C_rc2;
				pt2b_Fe_rc2 = mean_Fe_rc2 - mean_D_Fe_rc2;
				pt2b_Pb_rc2 = mean_Pb_rc2 - mean_D_Pb_rc2;
				pt2b_C_err_rc2 = TMath::Sqrt(mean_C_err_rc2*mean_C_err_rc2 + mean_D_C_err_rc2*mean_D_C_err_rc2);
				pt2b_Fe_err_rc2 = TMath::Sqrt(mean_Fe_err_rc2*mean_Fe_err_rc2 + mean_D_Fe_err_rc2*mean_D_Fe_err_rc2);
				pt2b_Pb_err_rc2 = TMath::Sqrt(mean_Pb_err_rc2*mean_Pb_err_rc2 + mean_D_Pb_err_rc2*mean_D_Pb_err_rc2);
			  
				Double_t x25_rc2[3]  = {pow(12,0.333), pow(56,0.333), pow(208,0.333)}; 
				Double_t y25_rc2[3]  = {pt2b_C_rc2, pt2b_Fe_rc2, pt2b_Pb_rc2};
				Double_t ex25_rc2[3] = {0, 0, 0};
				Double_t ey25_rc2[3] = {pt2b_C_err_rc2, pt2b_Fe_err_rc2, pt2b_Pb_err_rc2}; 

				TGraphErrors *gr25_rc2 = new TGraphErrors(3,x25_rc2,y25_rc2,ex25_rc2,ey25_rc2); 

				// Calculate Broadening for no RC
				pt2b_C = mean_C - mean_D_C;
				pt2b_Fe = mean_Fe - mean_D_Fe;
				pt2b_Pb = mean_Pb - mean_D_Pb;
				pt2b_C_err = TMath::Sqrt(mean_C_err*mean_C_err + mean_D_C_err*mean_D_C_err);
				pt2b_Fe_err = TMath::Sqrt(mean_Fe_err*mean_Fe_err + mean_D_Fe_err*mean_D_Fe_err);
				pt2b_Pb_err = TMath::Sqrt(mean_Pb_err*mean_Pb_err + mean_D_Pb_err*mean_D_Pb_err);
			  
				Double_t x25[3]  = {pow(12,0.333), pow(56,0.333), pow(208,0.333)}; 
				Double_t y25[3]  = {pt2b_C, pt2b_Fe, pt2b_Pb};
				Double_t ex25[3] = {0, 0, 0};
				Double_t ey25[3] = {pt2b_C_err, pt2b_Fe_err, pt2b_Pb_err}; 

				TGraphErrors *gr25 = new TGraphErrors(3,x25,y25,ex25,ey25); 
		
				gr25->SetMaximum(0.1);
				gr25->SetMinimum(0);
				TString title = Form("%.1f<Xb<%.1f %.1f<Q^{2}<%.1f %.1f<Z_{#pi^{+}}<%.1f   #void8   #pi^{+}", v_XB[i2], v_XB[i2+1], v_Q2[i1], v_Q2[i1+1], v_ZH[i3], v_ZH[i3+1]);
		
				TH1F *Graph4 = new TH1F("Graph4", title, 100, 1.95, 6.27715);
				Graph4->SetMinimum(0);
				Graph4->SetMaximum(0.1);
				Graph4->SetDirectory(0);
				Graph4->SetStats(0);
				Graph4->GetXaxis()->SetTitle("A^{1/3}");
				Graph4->GetYaxis()->SetTitle("#Delta<p_{t}^{2}> (Gev^{2}/c^{2})");
				gr25->SetHistogram(Graph4);
		
				gr25->Draw("AP"); 
				S->Update();
		
				S->GetFrame()->SetFillColor(21);
				S->GetFrame()->SetBorderSize(12);
				gr25->GetHistogram()->SetYTitle("#Delta<p_{t}^{2}> (Gev^{2}/c^{2})");
				gr25->GetHistogram()->SetXTitle("A^{1/3}");
				S->Modified();
				S->Update();
		
				TLatex *tex1 = new TLatex(2.75132,0.0455599,"CLAS PRELIMINARY");
				tex1->SetTextColor(22);
				tex1->SetTextSize(0.164659);
				tex1->SetLineWidth(2);
				tex1->Draw("same");
				S->Modified();
				S->Update();
		
				gr25->Draw("P:same");
				TMarker *m1 = new TMarker(x25[0],y25[0],21);
				m1->SetMarkerColor(2);
				m1->Draw("same");
				TMarker *m1_rc1 = new TMarker(x25_rc1[0], y25_rc1[0], 21);
				m1_rc1->SetMarkerColor(kBlue);
				m1_rc1->Draw("same");				
				TMarker *m1_rc2 = new TMarker(x25_rc2[0], y25_rc2[0], 21);
				m1_rc2->SetMarkerColor(kGreen);
				m1_rc2->Draw("same");	

				TMarker *m2 = new TMarker(x25[1],y25[1],23);
				m2->SetMarkerColor(2);
				m2->Draw("same");
				TMarker *m2_rc1 = new TMarker(x25_rc1[1], y25_rc1[1], 23);
				m2_rc1->SetMarkerColor(kBlue);
				m2_rc1->Draw("same");
				TMarker *m2_rc2 = new TMarker(x25_rc2[1], y25_rc2[1], 23);
				m2_rc2->SetMarkerColor(kGreen);
				m2_rc2->Draw("same");

				TMarker *m3 = new TMarker(x25[2],y25[2],8);
				m3->SetMarkerColor(2);
				m3->Draw("same");
				TMarker *m3_rc1 = new TMarker(x25_rc1[2], y25_rc1[2], 8);
				m3_rc1->SetMarkerColor(kBlue);
				m3_rc1->Draw("same");
				TMarker *m3_rc2 = new TMarker(x25_rc2[2], y25_rc2[2], 8);
				m3_rc2->SetMarkerColor(kGreen);
				m3_rc2->Draw("same");
			
				S->Update();
				TGraph *gr1 = new TGraph();
				TGraph *gr2 = new TGraph();
				TGraph *gr3 = new TGraph();
				gr1->SetMarkerStyle(21);
				gr2->SetMarkerStyle(23);
				gr3->SetMarkerStyle(8);
				gr1->SetMarkerColor(4);
				gr2->SetMarkerColor(4);
				gr3->SetMarkerColor(4);
				gr1->SetMarkerSize(2.5);
				gr2->SetMarkerSize(2.5);
				gr3->SetMarkerSize(2.5);
				TLegend *legend = new TLegend(0.902979,0.770192,0.993191,0.982151);
				legend->SetTextFont(72);
				legend->SetTextSize(0.0356984);
				legend->AddEntry(gr1," Carbon","p");
				legend->AddEntry(gr2," Iron","p");
				legend->AddEntry(gr3," Lead","p");
				legend->Draw();
				TString name=Form("source%d%d%d.cc",i1,i2,i3);
				TString name1=Form("%d%d%d.png",i1,i2,i3);
				TString name2=Form("%d%d%d.eps",i1,i2,i3);
				TString name3=Form("%d%d%d.ps",i1,i2,i3);
				MyC->cd();
				TString subtitle=Form("C_{m}=%.4f Fe_{m}=%.4f Pb_{m}=%.4f",pt2b_C,pt2b_Fe,pt2b_Pb);
				TLatex *tex = new TLatex(0.01,0.10241,subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(2);
				tex->Draw();
				subtitle=Form("C_{e}=%.4f Fe_{e}=%.4f Pb_{e}=%.4f",pt2b_C_err,pt2b_Fe_err,pt2b_Pb_err);
				tex = new TLatex(0.01,0.0281124,subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(2);
				tex->Draw();
				subtitle = Form("C_{m}=%.4f Fe_{m}=%.4f Pb_{m}=%.4f",pt2b_C_rc1,pt2b_Fe_rc1,pt2b_Pb_rc1);
				tex = new TLatex(0.34, 0.10241, subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(kBlue);
				tex->Draw();
				subtitle = Form("C_{e}=%.4f Fe_{e}=%.4f Pb_{e}=%.4f",pt2b_C_err_rc1,pt2b_Fe_err_rc1,pt2b_Pb_err_rc1);
				tex = new TLatex(0.34, 0.0281124, subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(kBlue);
				tex->Draw();
				subtitle = Form("C_{m}=%.4f Fe_{m}=%.4f Pb_{m}=%.4f",pt2b_C_rc2,0,pt2b_Pb_rc2);
				tex = new TLatex(0.67, 0.10241, subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(kGreen);
				tex->Draw();
				subtitle = Form("C_{e}=%.4f Fe_{e}=%.4f Pb_{e}=%.4f",pt2b_C_err_rc2,pt2b_Fe_err_rc2,pt2b_Pb_err_rc2);
				tex = new TLatex(0.67, 0.0281124, subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(kGreen);
				tex->Draw();
				//TLatex *tex1 = new TLatex(0.634894,0.10241,"Corrected data");
				//tex1->SetLineWidth(2);
				//tex1->Draw();
				MyC->Modified();
				MyC->Update();
				//MyC->SaveSource(name);
				//MyC->SetEditable(kFALSE);
				//MyC->SaveAs(name1);	

				out << pt2b_C_rc1 << "\t\t" << pt2b_C_err_rc1 << "\t\t" << pt2b_Fe_rc1 << "\t\t" << pt2b_Fe_err_rc1 << "\t\t" << pt2b_Pb_rc1 << "\t\t" << pt2b_Pb_err_rc1 << endl;
			}	
			out << endl;
		}
	}
	out.close();
}

