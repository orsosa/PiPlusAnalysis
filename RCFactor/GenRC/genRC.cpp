#include <fstream>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TRadCor.h"
#include "TNtuple.h"
#include "TMath.h"
#include "haprad_constants.h"

TString dataLoc;
TString Metal;

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

const Int_t N_METAL = 4;

int main(int argc, char **argv)
{
	Double_t f1, f3;
	Double_t NAZ;
	Double_t xbmas, q2mas, zhmas, ptmas, phimas;
	Double_t m = TMath::Power((kMassNeutron + kMassPion), 2);
	Double_t ebeam = 5.015;
	Int_t sysReturn;

	TRadCor rc;

	TString Metal;

	std::ofstream out;

	Q2_MIN = (Double_t) std::stod(argv[1]);
	Q2_MAX = (Double_t) std::stod(argv[2]);
	N_Q2 = (Int_t) std::stoi(argv[3]);
	XB_MIN = (Double_t) std::stod(argv[4]);
	XB_MAX = (Double_t) std::stod(argv[5]);
	N_XB = (Int_t) std::stoi(argv[6]);
	ZH_MIN = (Double_t) std::stod(argv[7]);
	ZH_MAX = (Double_t) std::stod(argv[8]);
	N_ZH = (Int_t) std::stoi(argv[9]); 
	PT_MIN = (Double_t) std::stod(argv[10]);
	PT_MAX = (Double_t) std::stod(argv[11]);
	N_PT = (Int_t) std::stoi(argv[12]);
	PHI_MIN = (Double_t) std::stod(argv[13]);
	PHI_MAX = (Double_t) std::stod(argv[14]);
	N_PHI = (Int_t) std::stoi(argv[15]);
	
	dataLoc = argv[16];	
	Metal = argv[17];
	
	sysReturn = system("cp " + dataLoc + Metal + "newphihist.root .");
	if(sysReturn == 256){
		std::cout << "File " << dataLoc << Metal  << " not found" << std::endl;
		return 0;
	}
	sysReturn = system("mv " + Metal + "newphihist.root newphihist.root");
	
	out.open("RCFactor" + Metal + ".txt");
	out << "Q2\tXb\tZh\tPt\tPhi\tSigmaB\tSigmaOb\tTail1\tTaile2\tFact_noex\tFact_ex" << std::endl;
	for(Int_t xb_i = 1; xb_i <= N_XB; xb_i++){
		xbmas = XB_MIN + (xb_i - 0.5)*delta_XB;
		for(Int_t q2_i = 1; q2_i <= N_Q2; q2_i++) {
			q2mas = Q2_MIN + (q2_i - 0.5)*delta_Q2;
			for(Int_t zh_i = 1; zh_i <= N_ZH; zh_i++) {
				zhmas = ZH_MIN + (zh_i - 0.5)*delta_ZH;
				for(Int_t pt_i = 1; pt_i <= N_PT; pt_i++){
					ptmas = PT_MIN + (pt_i - 0.5)*delta_PT;
					for(Int_t phi_i = 1; phi_i <= N_PHI; phi_i++){
						phimas = PHI_MIN + (phi_i - 0.5)*delta_PHI;
						rc.CalculateRCFactor(ebeam, xbmas, q2mas, zhmas, ptmas, phimas, m, NAZ);
						f1 = rc.GetFactor1();
						f3 = rc.GetFactor3();
						
						out << q2_i << "\t" << xb_i << "\t" << zh_i << "\t" << pt_i << "\t" << phi_i;
						out << "\t0\t0\t0";
						out << "\t0\t" << f1 << "\t" << f3 << std::endl;
					}			
				}
			}
		}
	}
	out.close();
	
	sysReturn = system("rm newphihist.root");
	
	return 0;		
}
