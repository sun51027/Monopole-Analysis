void FitSaturatedFraction_data(string year)
{
	TCanvas *c = new TCanvas("c","",800,600);
	TChain *tree = new TChain("monopoles");
        tree->Add("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/Data/Blind/*");
	
	int passHLT_Photon200;
	int nCandidates = 0;
	vector<float> *SubHits = 0;   
	vector<float> *SatSubHits = 0;

	tree->SetBranchAddress("passHLT_Photon200", &passHLT_Photon200);
	tree->SetBranchAddress("SatSubHits", &SatSubHits);
	tree->SetBranchAddress("SubHits", &SubHits);
	tree->SetBranchAddress("Point_N",&nCandidates);
	cout<<tree<<endl;
	cout<<"NEvents "<<tree->GetEntries()<<endl;	

	TH1D *h_fracSats = new TH1D("h_fracSats","",100,0,0.5);
	for(int ev=0; ev < tree->GetEntries();ev++){

		if(ev%100000==0)	cout<<ev<<"/"<<tree->GetEntries()<<endl;
		tree->GetEntry(ev);		
		for( int c = 0 ; c < nCandidates ; c++ ){

			double frac = (double)(*SatSubHits)[c]/(double)(*SubHits)[c];
			if (frac > 0 )	h_fracSats->Fill(frac);
		}
	}
	h_fracSats->Fit("gaus");
	
	h_fracSats->Draw();
	c->SaveAs("h_fracSats_gauss.pdf");
	
}
