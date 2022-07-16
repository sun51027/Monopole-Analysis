void FitSaturatedFraction(string year, string mass,bool matching_option)
{
	string matching;
	TChain *tree = new TChain("monopoles");
	tree->Add("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/Data/Blind/*");
        Bool_t passHLT_Photon200;
	Bool_t passHLT_Photon175;
	Bool_t passHLT_DoublePhoton70;
	Bool_t passHLT_DoublePhoton60;

        unsigned nCandidates;
	
        vector<double> *subHits=0;
        vector<double> *subSatHits=0;
        vector<double> *dEdXSig=0;
        tree->SetBranchAddress("passHLT_Photon200" , &passHLT_Photon200);
	tree->SetBranchAddress("passHLT_Photon175" , &passHLT_Photon175);
        tree->SetBranchAddress("passHLT_DoublePhoton70",&passHLT_DoublePhoton70);
        tree->SetBranchAddress("passHLT_DoublePhoton60",&passHLT_DoublePhoton60);
	tree->SetBranchAddress("cand_N",&nCandidates);
        tree->SetBranchAddress("cand_SubHits",&subHits);
        tree->SetBranchAddress("cand_SatSubHits",&subSatHits);
        tree->SetBranchAddress("cand_dEdXSig",&dEdXSig);
	const unsigned NEvents = tree->GetEntries();
	
	TH1D *h_fracSats = new TH1D("h_fracSats","",100,0,1);
        for(unsigned ev=0; ev<NEvents;ev++){

		if(ev%1000==0)	cout<<ev<<"/"<<NEvents<<endl;
		tree->GetEntry(ev);

		for( unsigned c = 0 ; c < nCandidates ; c++ ){
			h_fracSats->Fill((*subSatHits)[c]/(*subHits)[c]);
		}
	}
	h_fracSats->Draw();
}
