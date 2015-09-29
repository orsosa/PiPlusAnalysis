#include <stdlib.h>
#include "haprad_constants.h"
#include "TRadCor.h"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
  Double_t NAZ=0.5; //NAZ=Z/A
  Double_t Eb=5.014;
  Double_t Q2=1.5;
  Double_t W=2.5;
  Double_t x=Q2/(W*W+Q2+kMassProton*kMassProton);
  Double_t z=0.5;
  Double_t pt=0.5;
  const Double_t phi_min = -180.;
  const Double_t phi_max = 180.;
  const Int_t phi_n = 12;
  
  std::ofstream out;
  system("rm -f out.txt");
  out.open("out.txt");	
  
  TRadCor rc;
  Double_t m = TMath::Power((kMassNeutron + kMassPion),2);

  for(Int_t phi_i=1; phi_i <= phi_n; phi_i++){
    Double_t phi = phi_min + (phi_i-0.5)*(phi_max-phi_min)/phi_n;
    //rc.CalculateRCFactor(Eb,x,Q2,z,pt,phi,m);
    //rc.CalculateRCFactor(Eb,0.2,1.5,0.5,0.5,phi,m);
    rc.CalculateRCFactor(Eb,0.145,1.25,0.45,0.7,phi,m,NAZ);
    Double_t factor = rc.GetFactor1();
    out << phi_i << "\t" << rc.GetSigBorn() << "\t" << rc.GetSigObs() << "\t" << factor << "\n";
    
    std::cout << std::endl << "Sigma Observed: " << rc.GetSigObs() << std::endl;
    std::cout << std::endl << "Sigma Born: " << rc.GetSigBorn() << std::endl;
    std::cout << std::endl << "SIDIS tail: " << rc.GetTail(0) << std::endl;
    std::cout << std::endl << "RC factor: " << factor << std::endl;
  }

  return 0;
}
