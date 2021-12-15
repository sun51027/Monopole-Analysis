////////////////////////
//
//	Build New tree without Signal region
//	Only background
//
//	 1 = Loose cut dEdx < 7 f51 < 0.6
//	 2 = tight cut dEdx < 9 f51 < 0.85
//	 3 = Ublind
//	
//	Built by Shih Lin 2021 10 27
//
/////////////////////////
void MakeNewTree_fromTrim(int option, string dataset)
{
	TFile *fin = new TFile(("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/Data/data_2018/Monopole2018Data_"+dataset+"_trimmed.root").c_str());
        TTree *tree = (TTree*)fin->Get("monopoles");

 	TFile *oFile = new TFile(("../Data/Blind/BlindedData_2018"+dataset+".root").c_str(),"recreate");
	TTree *NewTree = new TTree("monopoles","Blinded data");

        float passHLT_Photon200;
        float Run=0;
        float Event=0;
        float Dist=0;
        float HIso= 0;
        float XYPar0 = 0;
        float XYPar1 = 0;
        float XYPar2 = 0;
        float RZPar0 = 0;
        float RZPar1 = 0;
        float RZPar2 = 0;
        float Eta = 0;
        float SatSubHits = 0;
        float SubHits = 0;
        float seedFrac = 0;
        float e55 =  0;
	
	tree->SetBranchAddress("run", &Run);
	tree->SetBranchAddress("event", &Event);
	tree->SetBranchAddress("Dist", &Dist);
	tree->SetBranchAddress("HIso", &HIso);
	tree->SetBranchAddress("XYPar0", &XYPar0);
	tree->SetBranchAddress("XYPar1", &XYPar1);
	tree->SetBranchAddress("XYPar2", &XYPar2);
	tree->SetBranchAddress("RZPar0", &RZPar0);
	tree->SetBranchAddress("RZPar1", &RZPar1);
	tree->SetBranchAddress("RZPar2", &RZPar2);
	tree->SetBranchAddress("Eta", &Eta);
	tree->SetBranchAddress("SatSubHits", &SatSubHits);
	tree->SetBranchAddress("SubHits", &SubHits);
	tree->SetBranchAddress("seedFrac", &seedFrac);
	tree->SetBranchAddress("cand_e55", &e55);
	tree->SetBranchAddress("passHLT_Photon200", &passHLT_Photon200);

	long long LastEvent = -1;
	vector<float> m_SubHits;   
	vector<float> m_SatSubHits;
	vector<float> m_XYPar0;
	vector<float> m_XYPar1;
	vector<float> m_XYPar2;
	vector<float> m_RZPar0;
	vector<float> m_RZPar1;
	vector<float> m_RZPar2;
	vector<float> m_Eta;
	vector<float> m_seedFrac;
	vector<float> m_e55;
	vector<float> m_Dist;
	vector<float> m_HIso;
	vector<float> m_dEdxSig;   
	long long m_event;
	long m_run;
	int m_nPoint=0;
	float dEdxSig = 0;
	int passHLT_Photon200_;

	int nPoint = 0;


	NewTree->Branch("SubHits",&m_SubHits);
	NewTree->Branch("SatSubHits",&m_SatSubHits);
	NewTree->Branch("XYPar0",&m_XYPar0);
	NewTree->Branch("XYPar1",&m_XYPar1);
	NewTree->Branch("XYPar2",&m_XYPar2);
	NewTree->Branch("RZPar0",&m_RZPar0);
	NewTree->Branch("RZPar1",&m_RZPar1);
	NewTree->Branch("RZPar2",&m_RZPar2);
	NewTree->Branch("Eta",&m_Eta);
	NewTree->Branch("seedFrac",&m_seedFrac);
	NewTree->Branch("e55",&m_e55);
	NewTree->Branch("Dist",&m_Dist);
	NewTree->Branch("HIso",&m_HIso);
	NewTree->Branch("dEdxSig",&m_dEdxSig);
	NewTree->Branch("Point_N",&m_nPoint);
	NewTree->Branch("event",&m_event);
	NewTree->Branch("run",&m_run);
	NewTree->Branch("passHLT_Photon200", &passHLT_Photon200_);
        bool SignalFlag = 0;
	int totalp=0;
	for(unsigned ev=0; ev < tree->GetEntries();ev++){
		tree->GetEntry(ev);
		//Blind signal region
		dEdxSig = sqrt(-TMath::Log(TMath::BinomialI(0.07, SubHits, SatSubHits)));
		if(option == 1){
			if(dEdxSig > 7 && seedFrac > 0.6){
				SignalFlag = 1;
				continue;
			}
		}
		else if(option == 2){
			if(dEdxSig > 9 && seedFrac > 0.85){
				SignalFlag = 1;
				continue;
			}
		}
		if(Event!=LastEvent){//only diff events will enteri
		    if(LastEvent > -1 && SignalFlag != 1){	
			NewTree->Fill();
			totalp++;
		    }
			LastEvent = Event;
			
			m_SubHits.clear(); 
			m_SatSubHits.clear();
			m_XYPar0.clear();
			m_XYPar1.clear();
			m_XYPar2.clear();
			m_RZPar0.clear();
			m_RZPar1.clear();
			m_RZPar2.clear();
			m_Eta.clear();
			m_seedFrac.clear();
			m_e55.clear();
			m_Dist.clear();
			m_HIso.clear();
			m_dEdxSig.clear();   
			m_nPoint = 0;
			m_run = 0;
			m_event = 0;
			passHLT_Photon200_ = -1;
			SignalFlag = 0;
			nPoint = 0;
		}
			nPoint++;//number of candidate in an event
			m_SubHits.push_back(SubHits);
			m_SatSubHits.push_back(SatSubHits);
			m_XYPar0.push_back(XYPar0);
			m_XYPar1.push_back(XYPar1);
			m_XYPar2.push_back(XYPar2);
			m_RZPar0.push_back(RZPar0);
			m_RZPar1.push_back(RZPar1);
			m_RZPar2.push_back(RZPar2);
			m_Eta.push_back(Eta);
			m_seedFrac.push_back(seedFrac);
			m_e55.push_back(e55);
			m_Dist.push_back(Dist);
			m_HIso.push_back(HIso);
			m_dEdxSig.push_back(dEdxSig);   
			m_nPoint = nPoint;
			m_run = Run;
			m_event = Event;
			passHLT_Photon200_ = passHLT_Photon200;

		if(ev==tree->GetEntries()-1){
			NewTree->Fill();
		}
			
	   if(ev%1000000==0) cout << ev << " / " << tree->GetEntries() << endl;
	} // end of entry loop
	cout<<"finish option "<<option<<" dataset "<<dataset<<endl;
	oFile->cd();
	NewTree->Write();

	oFile->Close();

}
