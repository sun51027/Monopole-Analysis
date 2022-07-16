//
//	bkg_ACBD.cc
//	Created by  Shih Lin
//	
//	Just Draw some plots, not for do ABCD
//
#include <algorithm>
#include <string>
#include "../interface/Candidate_Data.h"
#include "../interface/PlotSet.h"
#include "../interface/MonoCuts_data.h"
using namespace std;

  void MonoCuts::doAnalysis_data(vector<MonoCandidate> &cand, unsigned nCandidates, bool TRG, unsigned ev)
  {
	Clear();

     for(unsigned c=0;c<nCandidates;c++){

	MonoCandidate &cands = cand[c];
	bool PreselectionCut = evalPreselection(cands);
	bool Region1Cut = evalRegion1(cands);
	bool Region2Cut = evalRegion2(cands);
	bool Region3Cut = evalRegion3(cands);
	bool Region4Cut = evalRegion4(cands);
	bool Region5Cut = evalRegion5(cands);
	bool Region6Cut = evalRegion6(cands);
	bool Region7Cut = evalRegion7(cands);
	bool Region8Cut = evalRegion8(cands);
	bool Region9Cut = evalRegion9(cands);

	if(TRG && PreselectionCut ) Preselection.push_back(cands); 

	if(TRG && PreselectionCut && Region1Cut) Region1.push_back(cands);
	if(TRG && PreselectionCut && Region2Cut) Region2.push_back(cands);
	if(TRG && PreselectionCut && Region3Cut) Region3.push_back(cands);
	if(TRG && PreselectionCut && Region4Cut) Region4.push_back(cands);
	if(TRG && PreselectionCut && Region5Cut) Region5.push_back(cands);
	if(TRG && PreselectionCut && Region6Cut) Region6.push_back(cands);
	if(TRG && PreselectionCut && Region7Cut) Region7.push_back(cands);
	if(TRG && PreselectionCut && Region8Cut) Region8.push_back(cands);
	if(TRG && PreselectionCut && Region9Cut) Region9.push_back(cands);

      }//for cand loop
    	
	//cut flow events calculating
        sort(Preselection.begin(),Preselection.begin()+Preselection.size());

        sort(Region1.begin(),Region1.begin()+Region1.size());
        sort(Region2.begin(),Region2.begin()+Region2.size());
        sort(Region3.begin(),Region3.begin()+Region3.size());
        sort(Region4.begin(),Region4.begin()+Region4.size());
        sort(Region5.begin(),Region5.begin()+Region5.size());
        sort(Region6.begin(),Region6.begin()+Region6.size());
        sort(Region7.begin(),Region7.begin()+Region7.size());
        sort(Region8.begin(),Region8.begin()+Region8.size());
        sort(Region9.begin(),Region9.begin()+Region9.size());

	if(Preselection.size()>0) 
	{
		FillFlowHistogram(0,Preselection);
		//Fill only highest dEdx in the plots
		//whatever this event have 1 or 2 or 3 candidate, we only plus 1
		count++;	
//		cout<<"Preselection.size() "<<Preselection.size()<<endl;
	}
  }
  void MonoCuts::FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand){
	PlotSet &z = CutFlow[n];
//	for(int i=0; i < CutFlowCand.size() ;i++){
	//Draw event
//	  z.GetPlot(FracSatVNstrips)->Fill(CutFlowCand[0].subHits_,CutFlowCand[0].subSatHits_/CutFlowCand[0].subHits_);
//	  z.GetPlot(DedXSig)->Fill(CutFlowCand[0].dEdXSig_);
//	  z.GetPlot(F51)->Fill(CutFlowCand[0].f51_);
          z.GetPlot(ABCD)->Fill(CutFlowCand[0].f51_,CutFlowCand[0].dEdXSig_);
//	}
  }
  void MonoCuts::Clear(){

	Preselection.clear();
	Region1.clear();
	Region2.clear();
	Region3.clear();
	Region4.clear();
	Region5.clear();
	Region6.clear();
	Region7.clear();
	Region8.clear();
	Region9.clear();
  }

  void MonoCuts::WritePlots(TFile *oFile){
	oFile->cd(trigName_.c_str());
	CutFlow[0].WritePlot();
//	Profile[0].WriteProfile();
	
  }
	void PRINTT(){
		cout<<"count "<<count<<endl;
	}
  //Cuts parameters
  const double MonoCuts::xyp0Cut_=0.6;
  const double MonoCuts::xyp2Cut_=1000;
  const double MonoCuts::rzp0Cut_=10;
  const double MonoCuts::rzp1Cut_=999;
  const double MonoCuts::rzp2Cut_=0.005;
  const double MonoCuts::distCut_ = 0.5;
  const double MonoCuts::hIsoCut_= 10;
  const double MonoCuts::e55Cut_ = 200;
  const double MonoCuts::dEdXSigCut_ = 9;//preselection cut
  const double MonoCuts::dEdXSig_looseCut_ = 7;//preselection cut
  const double MonoCuts::f51Cut_ = 0.85;
  const double MonoCuts::f51_looseCut_ = 0.75;

void bkg_ABCD()
{
	TFile *oFile = new TFile("output/BlindAnalysis_ABCD_MC.root","recreate");
	TChain *Tree = new TChain("monopoles");
	tree->Add(("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MCNtuple"+year+"/"+mass+"/*.root").c_str());
//	Tree->Add("../../Data/Blind/BlindedData_2018AB.root");
//	Tree->Add("../../Data/Blind/BlindedData_2018C.root");
//	Tree->Add("../../Data/Blind/BlindedData_2018D.root");
//	TFile *fin = new TFile("../../Data/Blind/BlindedData_loose_2018all.root");
//        TTree *Tree = (TTree*)fin->Get("monopoles");

	int LastEvent = -1;
	vector<float> *SubHits = 0;   
	vector<float> *SatSubHits = 0;
	vector<float> *XYPar0 = 0;
	vector<float> *XYPar1 = 0;
	vector<float> *XYPar2 = 0;
	vector<float> *RZPar0 = 0;
	vector<float> *RZPar1 = 0;
	vector<float> *RZPar2 = 0;
	vector<float> *Eta = 0;
	vector<float> *seedFrac = 0;
	vector<float> *e55 = 0;
	vector<float> *Dist = 0;
	vector<float> *HIso = 0;
	vector<float> *dEdxSig = 0;   
	long long event;
	int run;
	int nPoint=0;
	int passHLT_Photon200;


	Tree->SetBranchAddress("SubHits",&SubHits);
	Tree->SetBranchAddress("SatSubHits",&SatSubHits);
	Tree->SetBranchAddress("XYPar0",&XYPar0);
	Tree->SetBranchAddress("XYPar1",&XYPar1);
	Tree->SetBranchAddress("XYPar2",&XYPar2);
	Tree->SetBranchAddress("RZPar0",&RZPar0);
	Tree->SetBranchAddress("RZPar1",&RZPar1);
	Tree->SetBranchAddress("RZPar2",&RZPar2);
	Tree->SetBranchAddress("Eta",&Eta);
	Tree->SetBranchAddress("seedFrac",&seedFrac);
	Tree->SetBranchAddress("e55",&e55);
	Tree->SetBranchAddress("Dist",&Dist);
	Tree->SetBranchAddress("HIso",&HIso);
	Tree->SetBranchAddress("dEdxSig",&dEdxSig);
	Tree->SetBranchAddress("Point_N",&nPoint);
	Tree->SetBranchAddress("event",&event);
	Tree->SetBranchAddress("passHLT_Photon200", &passHLT_Photon200);
	
	vector<MonoCandidate> cand(10);
	MonoCuts TrgAnalysis("HLT_Photon",oFile);
	int nParticle = 0;
        for(int ev=0; ev<Tree->GetEntries();ev++){
		if(ev%1000000==0) cout<<ev<<"/"<<Tree->GetEntries()<<endl;
		Tree->GetEntry(ev);
		nParticle = nPoint;	
		if(nParticle > cand.size()) cand.resize(nParticle);
		for( int i=0;i < nParticle ;i++){
			if((*dEdxSig)[i]==-0) (*dEdxSig)[i]=0;
			cand[i] = MonoCandidate(
		        (*subHits)[i],
		        (*subSatHits)[i],
		        (*dist)[i],
		        (*hIso)[i],
		        (*xyp0)[i],
		        (*xyp1)[i],
		        (*xyp2)[i],
		        (*rzp0)[i],
		        (*rzp1)[i],
		        (*rzp2)[i],
		        (*e55)[i],
		        (*f51)[i],
		        (*dEdXSig)[i],
		        (*eta)[i],
			event
			);
			
		}
		TrgAnalysis.doAnalysis_data(cand,nParticle,passHLT_Photon200,ev);
	}
//	TrgAnalysis.SignalEff("HLT_Photon200",0);
	TrgAnalysis.PRINTT();
	TrgAnalysis.WritePlots(oFile);
	oFile->Close();

}
