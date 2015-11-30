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
#include "TCanvas.h"
#include "TMarker.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"

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
Int_t N_METAL = 6;

TH1F *h_data_corr;
TH1F *h_data_corr_RC;
TFile *file;

int main(int argc, char **argv){
	Int_t met;
	TString Metal;
	
	for(Int_t i = 0; i < argc; i++){
		std::cout << "Arg " << i << std::endl;
		std::cout << argv[i] << std::endl;
	}
	if(argc < 24){
		std::cout << "The number of arguments is incorrect" << std::endl;
		std::cout << argc << std::endl;
		return 0;
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

	TFile *f = new TFile("/home/rodrigo/Data/GITData/PT2/Neg/pt2.root", "READ");
	
	Double_t mean[6];
	Double_t mean_err[6];
	Double_t pt2b[3];
	Double_t pt2b_err[3];
	
	Double_t mean_rc[6];
	Double_t mean_err_rc[6];
	Double_t pt2b_rc[3];
	Double_t pt2b_err_rc[3];

	ofstream out;
	out.open("Broadening.txt");
	out << "Transverse momentum broadening vs A1/3" << std::endl;
	out << "Zh\t\t\tdPt2 C\t\t\tErr C\t\t\tdPt2 Fe\t\t\tErr Fe\t\t\tdPt2 Pb\t\t\tErr Pb" << std::endl;
	
	TH1F *h;
	TH1F *h_rc;
	TCanvas *MyC;
	TPad *S;
	
	MyC = new TCanvas("MyC","transverse momuntum broadening of leading pions",6,21,1117,499);
	
	Int_t noHist = 0;

	for(q2i = 0; q2i < N_Q2; q2i++){
		for(xbi = 0; xbi < N_XB; xbi++){
			for(zhi = 0; zhi < N_ZH; zhi++){
				if(q2i == 3 && xbi == 2 && zhi == 3) continue;
				for(met = 0; met < N_METAL; met++){
					if(met == 0) Metal = "C"; 
					else if(met == 1) Metal = "Fe";
					else if(met == 2) Metal = "Pb";
					else if(met == 3) Metal = "DC";
					else if(met == 4) Metal = "DFe";
					else Metal = "DPb";
					
					std::cout << "Reading hist " << Metal << "_hist_" << q2i << xbi << zhi << std::endl;
					h = (TH1F*) f->Get((const char*)Form("%s_hist_%d%d%d", (const char*) Metal, q2i, xbi, zhi));
					if(h == 0){
						noHist = 1;
						break;
					}
					mean[met] = h->GetMean();
					mean_err[met] = h->GetMeanError();
					
					if(RCOn){
						h_rc = (TH1F*) (f->Get((const char*)Form("%s_hist_rc_%d%d%d", (const char*) Metal, q2i, xbi, zhi)));
						mean_rc[met] = h_rc->GetMean();
						mean_err_rc[met] = h_rc->GetMeanError();						
					}
				}
				
				if(noHist){
					noHist = 0;
					continue;
				}
				
				std::cout << "Finished getting the means" << std::endl;
				
				MyC->SetFillColor(42);
				MyC->SetGrid();
				MyC->GetFrame()->SetFillColor(21);
				MyC->GetFrame()->SetBorderSize(12);
				MyC->Update();
				
				std::cout << "Creating Pag" << std::endl;
				S = new TPad("S", "transverse momuntum broadening of leading pions",0,0.160072,1,1);
				S->Draw();
				S->cd();
				S->SetFillColor(42);
				S->SetGrid();
				S->GetFrame()->SetFillColor(21);
				S->GetFrame()->SetBorderSize(12);
				S->Update();

				for(met = 0; met < N_METAL/2; met++){
					pt2b[met] = mean[met] - mean[met+3];
					pt2b_err[met] = TMath::Sqrt(mean_err[met]*mean_err[met] + mean_err[met+3]*mean_err[met+3]);
					if(RCOn){
						pt2b_rc[met] = mean_rc[met] - mean_rc[met+3];
						pt2b_err_rc[met] = TMath::Sqrt(mean_err_rc[met]*mean_err_rc[met] + mean_err_rc[met+3]*mean_err_rc[met+3]);
					}
				}
							  
				Double_t x25[3]  = {pow(12,0.333), pow(56,0.333), pow(208,0.333)}; 
				Double_t y25[3]  = {pt2b[0], pt2b[1], pt2b[2]};
				Double_t ex25[3] = {0, 0, 0};
				Double_t ey25[3] = {pt2b_err[0], pt2b_err[1], pt2b_err[2]}; 

				TGraphErrors *gr25 = new TGraphErrors(3,x25,y25,ex25,ey25); 
			  
				Double_t x25_rc[3]  = {pow(12,0.333), pow(56,0.333), pow(208,0.333)}; 
				Double_t y25_rc[3]  = {pt2b_rc[0], pt2b_rc[1], pt2b_rc[2]};
				Double_t ex25_rc[3] = {0, 0, 0};
				Double_t ey25_rc[3] = {pt2b_err_rc[0], pt2b_err_rc[1], pt2b_err_rc[2]}; 

				TGraphErrors *gr25_rc = new TGraphErrors(3, x25_rc, y25_rc, ex25_rc, ey25_rc); 

				gr25->SetMaximum(0.1);
				gr25->SetMinimum(0);
				TString title = Form("%.1f<Xb<%.1f %.1f<Q^{2}<%.1f %.1f<Z_{#pi^{+}}<%.1f   #void8   #pi^{+}", v_XB[xbi], v_XB[xbi+1], v_Q2[q2i], v_Q2[q2i+1], v_ZH[zhi], v_ZH[zhi+1]);
		
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
				if(RCOn) gr25_rc->Draw("P:same");
								
				TMarker *m1 = new TMarker(x25[0],y25[0],21);
				m1->SetMarkerColor(2);
				m1->Draw("same");
				TMarker *m2 = new TMarker(x25[1],y25[1],23);
				m2->SetMarkerColor(2);
				m2->Draw("same");
				TMarker *m3 = new TMarker(x25[2],y25[2],8);
				m3->SetMarkerColor(2);
				m3->Draw("same");
				if(RCOn){
					TMarker *m1_rc = new TMarker(x25_rc[0], y25_rc[0], 21);
					m1_rc->SetMarkerColor(kBlue);
					m1_rc->Draw("same");								
					TMarker *m2_rc = new TMarker(x25_rc[1], y25_rc[1], 23);
					m2_rc->SetMarkerColor(kBlue);
					m2_rc->Draw("same");
					TMarker *m3_rc = new TMarker(x25_rc[2], y25_rc[2], 8);
					m3_rc->SetMarkerColor(kBlue);
					m3_rc->Draw("same");
				}				
			
				S->Update();
				TGraph *gr1 = new TGraph();
				TGraph *gr2 = new TGraph();
				TGraph *gr3 = new TGraph();
				gr1->SetMarkerStyle(21);
				gr2->SetMarkerStyle(23);
				gr3->SetMarkerStyle(8);
				gr1->SetMarkerColor(kRed);
				gr2->SetMarkerColor(kRed);
				gr3->SetMarkerColor(kRed);
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
				TString name1=Form("%d%d%d.png", q2i, xbi, zhi);
				MyC->cd();
				TString subtitle=Form("C_{m}=%.4f Fe_{m}=%.4f Pb_{m}=%.4f", pt2b[0], pt2b[1], pt2b[2]);
				TLatex *tex = new TLatex(0.01,0.10241,subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(2);
				tex->Draw();
				subtitle=Form("C_{e}=%.4f Fe_{e}=%.4f Pb_{e}=%.4f", pt2b_err[0], pt2b_err[1], pt2b_err[2]);
				tex = new TLatex(0.01,0.0281124,subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(2);
				tex->Draw();
				subtitle = Form("C_{m}=%.4f Fe_{m}=%.4f Pb_{m}=%.4f", pt2b_rc[0], pt2b_rc[1], pt2b_rc[2]);
				tex = new TLatex(0.34, 0.10241, subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(kBlue);
				tex->Draw();
				subtitle = Form("C_{e}=%.4f Fe_{e}=%.4f Pb_{e}=%.4f", pt2b_err_rc[0], pt2b_err_rc[1], pt2b_err_rc[2]);
				tex = new TLatex(0.34, 0.0281124, subtitle);
				tex->SetLineWidth(2);
				tex->SetTextColor(kBlue);
				tex->Draw();
				MyC->Modified();
				MyC->Update();
				MyC->SaveAs(name1);	
				
				MyC->Clear();
				
				std::cout << "Deleting Pointers" << std::endl;				
				std::cout << "After Deleting Pointers" << std::endl;		

				out << pt2b_rc[0] << "\t\t" << pt2b_err_rc[0] << "\t\t" << pt2b_rc[1] << "\t\t" << pt2b_err_rc[1] << "\t\t" << pt2b_rc[2] << "\t\t" << pt2b_err_rc[2] << std::endl;
			}	
			out << std::endl;
		}
	}
	out.close();
}

