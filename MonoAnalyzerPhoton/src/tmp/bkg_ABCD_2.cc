//
//	bkg_ACBD.cc
//	Created by  Shih Lin
//	Use updated data structure to calculate ABCD 
//	**Jack propose.... but it was wrong	
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
	int label = 0;
     for(unsigned c=0;c<nCandidates;c++){
	MonoCandidate &cands = cand[c];
	bool PreselectionCut = evalPreselection(cands);
	if(TRG && PreselectionCut ) Preselection.push_back(cands); 
      }
    	
        sort(Preselection.begin(),Preselection.begin()+Preselection.size());

	if(Preselection.size()>0) 
	{
		FillFlowHistogram(0,Preselection);
		for(int c= 0;c<Preselection.size();c++){
			if(Preselection[c].dEdXSig_ > 7 && Preselection[c].f51_ > 0.75){
				selectedCand.push_back(Preselection[c]);
			}
			else{
				otherCand.push_back(Preselection[c]);
			}
		}	
	}

	if(selectedCand.size()>0){
		sort(selectedCand.begin(),selectedCand.begin()+selectedCand.size());
		FillFlowHistogram(0,selectedCand);
		label++;
	}
	else if(otherCand.size()>0 && selectedCand.size() == 0){
		sort(otherCand.begin(),otherCand.begin()+otherCand.size());
		FillFlowHistogram(0,otherCand);
		label++;
	}
	if( label!=0)	cout<<"Every event only label 1: ev "<<ev<<" "<<label<<endl;
		
  }
  void MonoCuts::FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand){
	PlotSet &z = CutFlow[n];
        z.GetPlot(ABCD)->Fill(CutFlowCand[0].f51_,CutFlowCand[0].dEdXSig_);
  }
  void MonoCuts::Clear(){

	Preselection.clear();
	selectedCand.clear();
	otherCand.clear();
  }

  void MonoCuts::WritePlots(TFile *oFile){
	oFile->cd(trigName_.c_str());
	CutFlow[0].WritePlot();
//	Profile[0].WriteProfile();
	
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

void bkg_ABCD_2()
{
	TFile *oFile = new TFile("output/BlindAnalysis_ABCD_Jack.root","recreate");
	TChain *Tree = new TChain("monopoles");
	Tree->Add("../../Data/Blind/BlindedData_2018AB.root");
	Tree->Add("../../Data/Blind/BlindedData_2018C.root");
	Tree->Add("../../Data/Blind/BlindedData_2018D.root");
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
			(*SubHits)[i],
			(*SatSubHits)[i],
			(*Dist)[i],
			(*HIso)[i],
			(*XYPar0)[i],
			(*XYPar1)[i],
			(*XYPar2)[i],
			(*RZPar0)[i],
			(*RZPar1)[i],
			(*RZPar2)[i],
			(*e55)[i],
			(*seedFrac)[i],
			(*dEdxSig)[i], 
			(*Eta)[i],
			event
			);
			
		}
		TrgAnalysis.doAnalysis_data(cand,nParticle,passHLT_Photon200,ev);
	}
//	TrgAnalysis.SignalEff("HLT_Photon200",0);
	TrgAnalysis.WritePlots(oFile);
	oFile->Close();

}
