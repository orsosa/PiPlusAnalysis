#include <fstream>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TRadCor.h"
#include "TNtuple.h"
#include "TMath.h"
#include "haprad_constants.h"
#include "BinMgr.hxx"

TString dataLoc;
TString Metal;

Int_t N_Q2, N_XB, N_ZH, N_PT, N_PHI;

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
	Double_t a1, a2, a3;
	Int_t sysReturn;
	Int_t dataSL;

	TRadCor rc;

	std::ofstream out;

	ifstream
	  
	dataLoc = argv[16];	
	Metal = argv[17];
	
	delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
	delta_XB = (XB_MAX-XB_MIN)/N_XB;
	delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
	delta_PT = (PT_MAX-PT_MIN)/N_PT;
	delta_PHI = (PHI_MAX-PHI_MIN)/N_PHI;
	
	dataSL = dataLoc.Length();
	if(dataLoc(dataSL-1, 1) != "/"){
		dataLoc = dataLoc + "/";
	}	
	
	a1 = 1;
	a2 = 0;
	a3 = a1/a2;
	
	std::cout << "The following settings are being used" << std::endl;
	std::cout << Q2_MIN << " < Q2 < " << Q2_MAX << ", N_Q2 = " << N_Q2 << std::endl;
	std::cout << XB_MIN << " < Xb < " << XB_MAX << ", N_XB = " << N_XB << std::endl;
	std::cout << ZH_MIN << " < Zh < " << ZH_MAX << ", N_ZH = " << N_ZH << std::endl;	
	std::cout << PT_MIN << " < Pt < " << PT_MAX << ", N_PT = " << N_PT << std::endl;	
	std::cout << PHI_MIN << " < Phi < " << PHI_MAX << ", N_PHI = " << N_PHI << std::endl;	
	std::cout << std::endl;
	std::cout << "Data Directory = " << dataLoc << std::endl;
	std::cout << "Running only Metal " << Metal << std::endl;
	
	sysReturn = system("cp " + dataLoc + Metal + "newphihist.root .");
	if(sysReturn == 256){
		std::cout << "File " << dataLoc << Metal  << " not found" << std::endl;
		return 0;
	}
	sysReturn = system("mv " + Metal + "newphihist.root newphihist.root");
	
	if(Metal == "C") NAZ = 0.5;
	else if(Metal == "Fe") NAZ = 0.5;
	else if(Metal == "Pb") NAZ = 82./208.;
	else if(Metal == "D_C") NAZ = 0.5;
	else if(Metal == "D_Fe") NAZ = 0.5;
	else if(Metal == "D_Pb") NAZ = 0.5;
	else NAZ = 0.5;
	
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
						rc.CalculateRCFactor(5.014, xbmas, q2mas, zhmas, ptmas, phimas, m, NAZ);
						f1 = rc.GetFactor1();
						f3 = rc.GetFactor3();
						if(TMath::IsNaN(f1) || f1 == a3) f1 = 0;
						if(TMath::IsNaN(f3) || f3 == a3) f3 = 0;
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
