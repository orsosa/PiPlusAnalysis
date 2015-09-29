/*
 * ===========================================================================
 *
 *       Filename:  Tuple_test.cxx
 *
 *    Description:  Create tuples from HapRad
 *
 *        Version:  1.0
 *        Created:  21/01/11 16:16:09
 *  Last Modified:  24/01/11 15:47:24
 *       Compiler:  gcc
 *
 *         Author:  Ricardo Oyarzun
 *          Email:  royarzun@alumnos.inf.utfsm.cl
 *
 * ===========================================================================
 */

#include <iostream>
#include <fstream>
#include <exception>
#include "TFile.h"
#include "TRadCor.h"
#include "TNtuple.h"
#include "haprad_constants.h"

void StudyBins()
{
    Double_t f1, f2, f3;
    TRadCor rc;
    Double_t m = TMath::Power((kMassNeutron + kMassPion), 2);

    const Double_t q2_min = 1.;
    const Double_t q2_max = 4.;
    const Double_t nu_min = 2.2;
    const Double_t nu_max = 4.0;
    const Double_t pt2_min = 0.0;
    const Double_t pt2_max = 1.5;
    const Double_t z_min = 0.4;
    const Double_t z_max = 0.7;
    const Double_t phi_min = -180.;
    const Double_t phi_max = 180.;
    const Int_t nu_n = 3;
    const Int_t z_n = 3;
    const Int_t pt2_n = 5;
    const Int_t q2_n = 3;
    const Int_t phi_n = 12;
    Double_t numas,xmas,q2mas,zmas,ptmas,phimas;

    std::ofstream out;
    system("rm -f out.txt");
    out.open("out.txt");	
    //out << "nu" << "\t" << "q2" << "\t" << "z" << "\t" << "pt" << "\t" << "phi" << "\t" << "f1" << "\t" << "f2" << "\t" << "f3" << "\t" << std::endl;	
    out << "nu" << "\t" << "q2" << "\t" << "z" << "\t" << "pt" << "\t" << "phi" << "\t" << "f" << "\t" << std::endl;
    for(Int_t nu_i = 1; nu_i <= nu_n; nu_i++) {
      numas = nu_min + (nu_i-0.5)*(nu_max-nu_min)/nu_n;
      for(Int_t q2_i = 1; q2_i <= q2_n; q2_i++) {
	q2mas = q2_min + (q2_i-0.5)*(q2_max-q2_min)/q2_n;
        xmas = q2mas/2./numas/kMassProton;
        for(Int_t z_i = 1; z_i <= z_n; z_i++) {
	  zmas = z_min+(z_i-0.5)*(z_max-z_min)/z_n;
	  for(Int_t pt2_i = 1; pt2_i <= pt2_n; pt2_i++) {
	    ptmas = sqrt(pt2_min+(pt2_i-0.5)*(pt2_max-pt2_min)/pt2_n);
	    std::cout << nu_i << "\t" << q2_i << "\t" << z_i << "\t" << pt2_i << "\t" << "$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
	    for(Int_t phi_i=1; phi_i <= phi_n; phi_i++){ 
	      phimas = phi_min + (phi_i-0.5)*(phi_max-phi_min)/phi_n;
	      rc.CalculateRCFactor(5.015, xmas, q2mas, zmas, ptmas, phimas, m);
	      f1 = rc.GetFactor1();
	      f2 = rc.GetFactor2();
	      f3 = rc.GetFactor3();
	      //out << nu_i << "\t" << q2_i << "\t" << z_i << "\t" << pt2_i << "\t" << phi_i << "\t" << f1 << "\t" << f2 << "\t" << f3 << "\n";
	      out << nu_i << "\t" << q2_i << "\t" << z_i << "\t" << pt2_i << "\t" << phi_i << "\t" << f3 << "\n";
	    }
	  }
	}
      }
    }
    
    return;

}

int main()
{

    StudyBins();
    return 0;

}

