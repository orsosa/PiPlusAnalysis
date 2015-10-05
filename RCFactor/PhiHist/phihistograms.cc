const Double_t Q2_MIN = 1.;
const Double_t Q2_MAX = 4.;
const Int_t N_Q2 = 6;
const Double_t ZH_MAX = 1.;
const Double_t ZH_MIN = 0.;
const Int_t N_ZH = 10;
const Double_t XB_MAX = 0.55;
const Double_t XB_MIN = 0.1;
const Int_t N_XB = 5;
const Double_t Phi_pq_MAX = 180.;
const Double_t Phi_pq_MIN = -180.;
const Int_t N_Phi_pq = 12;
const Double_t PT_MAX = 1;
const Double_t PT_MIN = 0.;
const Int_t N_PT = 5;
Int_t i,j,k,l,m;

phihistograms(){
    TString f_location = "/home/rodrigo/Dropbox/FisCCTVal2.0/AAcAcc/C/Data"; 
    TString f_name = "5_dim_dist.root";
    f_name = f_location + "/" + f_name;
	TFile *f = new TFile(f_name, "READ");
	TNtuple *ntuple = (TNtuple*) f->Get("fit_data");
	TFile *file = new TFile("phihist.root", "RECREATE");
	Float_t q2mas, xbmas, zhmas, ptmas, phimas, val, err;
	Float_t Amas, Acmas, Accmas;
	Float_t Aermas, Acermas, Accermas;
	Int_t nentries = ntuple->GetEntries();
	Int_t nhist = ntuple->GetEntries()/12;
	Float_t sumchi = 0;

	ntuple->SetBranchAddress("Q2", &q2mas);
	ntuple->SetBranchAddress("Xb", &xbmas);
	ntuple->SetBranchAddress("Zh", &zhmas);
	ntuple->SetBranchAddress("Pt", &ptmas);
	ntuple->SetBranchAddress("Phi", &phimas);
	ntuple->SetBranchAddress("Val", &val);
	ntuple->SetBranchAddress("Err", &err);

	TNtuple *newntuple = new TNtuple("AAcAcc_data", "AAcAcc_data", "Q2:Xb:Zh:Pt:A:Aer:Ac:Acer:Acc:Accer");
	TF1 *func = new TF1("fit", "[0]+TMath::Cos(x*TMath::Pi()/180)*[1]+TMath::Cos(2*x*TMath::Pi()/180)*[2]");

	for(Int_t i = 0; i < nentries; i = i + 12){
		ntuple->GetEntry(i);
		TH1F *h1 = new TH1F((const char*) Form("PhiDist Q2=%.3f Xb=%.3f Zh=%.3f Pt=%.3f", q2mas, xbmas, zhmas, ptmas), 
									(const char*) Form("PhiDist Q2=%.3f Xb=%.3f Zh=%.3f Pt=%.3f", q2mas, xbmas, zhmas, ptmas), 12, -180, 180);
		Int_t empty = 0;
		for(Int_t j = 1; j < 13; j++){
			ntuple->GetEntry(i+j-1);
			h1->SetBinContent(j, val);
			if(h1->GetBinContent(j) == 0)
				empty++;
			h1->SetBinError(j, err*1.04);
		}
		cout << "Empty: " << "\t" << empty << endl;
		
		if(empty != 12){
			h1->Fit(func);
			h1->Write();
			if(func->GetNDF() != 0){
                sumchi = sumchi + func->GetChisquare();
				Amas = func->GetParameter(0);
				Aermas = func->GetParError(0);
				Acmas = func->GetParameter(1);
				Acermas = func->GetParError(1);
				Accmas = func->GetParameter(2);
				Accermas = func->GetParError(2);
				newntuple->Fill(q2mas, xbmas, zhmas, ptmas, Amas, Aermas, Acmas, Acermas, Accmas, Accermas);
			}
		}
		else{
		    cout << "Empty histogram Q2=" << q2mas << " Xb=" << xbmas << " Zh=" << zhmas << " Pt=" << ptmas << endl;
		}
	}
	newntuple->Write();
	file->Close();
	cout << "ChiSquare = " << sumchi << endl;
}
