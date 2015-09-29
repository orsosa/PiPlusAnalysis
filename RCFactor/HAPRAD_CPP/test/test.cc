#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
  
Double_t Eb=5.014;
Double_t x=0.15; 
Double_t Q2=1.25;
Double_t z=0.05;
Double_t p_t=0.1;
const Double_t phi_min = -3.1415926;
const Double_t phi_max = 3.1415926;
const Int_t phi_n = 12;

int main(int argc, char *argv[]) {
  
  for(Int_t phi_i=1; phi_i <= phi_n; phi_i++) {
    Double_t phi = phi_min + (phi_i-0.5)*(phi_max-phi_min)/phi_n;
    
    std::cout << phi_i << "\t" << GetSigBorn(phi) << "\n";
    
  }
  
  return 0;
}

Int_t GetSigBorn(Double_t phi) {


Double_t N = kBarn*(kPi*kAlpha*kAlpha*Y*Sx*kMassProton)/(2*SqrtLq);


}
