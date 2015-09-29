TString f_location = "/home/rodrigo/Data/PiData/"; // The location of data files with ntuples inside
TString fd_ext = "_data.root";
TString fs_ext = "_simul.root";
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
TH1F *htmp_data_corr;
TFile *plots = new TFile("5_dim_dist.root","RECREATE");
TNtuple *fit_data = new TNtuple("fit_data", "DATA FOR 5 DIM. FIT", "Q2:Xb:Zh:Pt:Phi:Val:Err");

pt_5dim_acc_corr(){
    const Double_t delta_Q2 = (Q2_MAX-Q2_MIN)/N_Q2;
    const Double_t delta_XB = (XB_MAX-XB_MIN)/N_XB;
    const Double_t delta_ZH = (ZH_MAX-ZH_MIN)/N_ZH;
    const Double_t delta_PT = (PT_MAX-PT_MIN)/N_PT;
    for(Int_t i=0;i<N_Q2;i++){
        for(Int_t j=0;j<N_XB;j++){
            for(Int_t k=0;k<N_ZH;k++){
                for(Int_t l=0;l<N_PT;l++){
                    cout << "Running nº: " << i << j << k << l << endl;
                    run_file("C",Q2_MIN+i*delta_Q2, Q2_MIN+(i+1)*delta_Q2,
                                XB_MIN+j*delta_XB, XB_MIN+(j+1)*delta_XB,
                                ZH_MIN+k*delta_ZH, ZH_MIN+(k+1)*delta_ZH,
                                PT_MIN+l*delta_PT, PT_MIN+(l+1)*delta_PT);
                    plots->cd();					
                }
            }
        }
    }
    plots->cd(); 
    fit_data->Write();
    plots->Close();
    return 0;
}

void run_file(const Char_t Metal[], Double_t Q2_MIN, Double_t Q2_MAX, Double_t XB_MIN, Double_t XB_MAX, Double_t ZH_MIN, Double_t ZH_MAX, Double_t PT_MIN, Double_t PT_MAX) {
    TCut p_cut = Form("((T4>-0.55 && P<2.7) || P>2.7)");
    TCut Q2_cut = Form("Q2>%f && Q2<%f", Q2_MIN, Q2_MAX);
    TCut Xb_cut = Form("Xb>%f && Xb<%f", XB_MIN, XB_MAX);
    TCut Zh_cut = Form("Zh>%f && Zh<%f", ZH_MIN, ZH_MAX);
    TCut Pt_cut = Form("Pt>%f && Pt<%f", PT_MIN, PT_MAX);
    TCut Xf_cut = Form("Xf<0");
    TCut Target_cut = Form("TargType==%d", 2);
    TCut cuts = p_cut && Q2_cut && Xb_cut && Zh_cut && Pt_cut && Xf_cut && Target_cut;
    TCut cuts_simul = p_cut && Q2_cut && Xb_cut && Zh_cut&& Pt_cut;
    TCut Phi_pq_cut;

    TChain *ntuple = new TChain("ntuple_data");
    ntuple->Add(f_location + Metal + fd_ext);
    ntuple->Draw(">>list", cuts, "goff");
    ntuple->SetEventList((TEventList*)gDirectory->Get("list"));

    TChain *accept = new TChain("ntuple_accept");
    accept->Add(f_location + Metal + "1" + fs_ext);
    accept->Add(f_location + Metal + "2" + fs_ext);
    accept->Add(f_location + Metal + "3" + fs_ext);
    accept->Add(f_location + Metal + "4" + fs_ext);
    accept->Draw(">>list_acc",cuts_simul,"goff");
    accept->SetEventList((TEventList*)gDirectory->Get("list_acc"));

    TChain *thrown = new TChain("ntuple_thrown");
    thrown->Add(f_location + Metal + "1" + fs_ext);
    thrown->Add(f_location + Metal + "2" + fs_ext);
    thrown->Add(f_location + Metal + "3" + fs_ext);
    thrown->Add(f_location + Metal + "4" + fs_ext);
    thrown->Draw(">>list_thr", cuts_simul, "goff");
    thrown->SetEventList((TEventList*)gDirectory->Get("list_thr"));

    ntuple->Draw((const char*)Form("PhiPQ>>htmp_data(%d,%f,%f)", N_Phi_pq, Phi_pq_MIN, Phi_pq_MAX), "", "goff");
    TH1F *htmp_data = (TH1F*)gDirectory->GetList()->FindObject("htmp_data");
    htmp_data->Sumw2();
    accept->Draw((const char*)Form("PhiPQ>>htmp_acc(%d,%f,%f)", N_Phi_pq, Phi_pq_MIN, Phi_pq_MAX), "", "goff");
    TH1F *htmp_acc = (TH1F*)gDirectory->GetList()->FindObject("htmp_acc");
    htmp_acc->Sumw2();
    thrown->Draw((const char*)Form("PhiPQ>>htmp_thr(%d,%f,%f)", N_Phi_pq, Phi_pq_MIN, Phi_pq_MAX), "", "goff");
    TH1F *htmp_thr = (TH1F*)gDirectory->GetList()->FindObject("htmp_thr");
    htmp_thr->Sumw2();
    TH1F *htmp_acc_ratio = new TH1F("htmp_acc_ratio", "", N_Phi_pq, Phi_pq_MIN, Phi_pq_MAX);
    htmp_data_corr = new TH1F("htmp_data_corr", "", N_Phi_pq, Phi_pq_MIN, Phi_pq_MAX);
    htmp_acc_ratio->Divide(htmp_acc,htmp_thr,1,1,"B");
    htmp_data_corr->Divide(htmp_data,htmp_acc_ratio,1,1);

    const Double_t delta_Phi = (Phi_pq_MAX-Phi_pq_MIN)/N_Phi_pq;
    
    for(Int_t ii=1; ii<=N_Phi_pq; ii++){
        Double_t PhiVal = Phi_pq_MIN+(ii-0.5)*delta_Phi;
        Double_t Val = htmp_data_corr->GetBinContent(ii);
        Double_t Err = htmp_data_corr->GetBinError(ii);
        fit_data->Fill((Q2_MIN+Q2_MAX)/2., (XB_MIN+XB_MAX)/2., (ZH_MIN+ZH_MAX)/2., (PT_MIN+PT_MAX)/2., PhiVal, Val, Err);
    }
    
    delete htmp_acc_ratio;
    delete htmp_data_corr;
    delete ntuple;
    delete accept;
    delete thrown;
}
