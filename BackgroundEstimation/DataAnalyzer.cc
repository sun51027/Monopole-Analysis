#include "iostream"
#include "TAttMarker.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "math.h"
#include "/wk_cms2/shihlin0314/CMSSW_8_0_29/src/BackgroundEstimation/Candidate_Data.h"
#include "/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MonoAnalyzerPhoton/PlotSet.h"

#include <algorithm>
#include <string>
using namespace std;
enum PlotName{
  FracSatVNstrips = 0, // fraction of saturated strips vs. number of strips
  DedXSig,             // dE/dX significance
  XYPar0,//d0
  XYPar1,//phi0
  XYPar2,//radius
  RZPar0,
  RZPar1,
  RZcurv,              // RZ curvature/uncertainty
  E55,
  F51,                 // frac 51
  HcalIso,             // Hcal Iso
  Dist,
  ABCD,
  Spike
};
static const unsigned nPlot = 20U;
class MonoCuts:public MonoCandidate
{
public:
  MonoCuts(){}
  MonoCuts(const string &trName,TFile *openFile):trigName_(trName){

        openFile->mkdir(trName.c_str());
        Plot.resize(1U);
        PlotSet &x = Plot[0];
        x.CreatPlot(FracSatVNstrips,new TH2D("FracSatVNstrips","",100,0,1000,100,0,1));
        x.CreatPlot(DedXSig,new TH1D("DedXSig","",100,0,30));
        x.CreatPlot(Dist,new TH1D("Dist","",100,-1,1));
        x.CreatPlot(HcalIso,new TH1D("HcalIso","",100,-1,10));
        x.CreatPlot(XYPar0,new TH1D("XYPar0","",100,-3,3));
        x.CreatPlot(XYPar1,new TH1D("XYPar1","",100,-10,10));
        x.CreatPlot(XYPar2,new TH1D("XYPar2","",100,-2000,2000));
        x.CreatPlot(RZPar0,new TH1D("RZPar0","",100,-20,20));
        x.CreatPlot(RZPar1,new TH1D("RZPar1","",100,-20,20));
        x.CreatPlot(RZcurv,new TH1D("RZcurv","",100,-0.01,0.01));
        x.CreatPlot(E55,new TH1D("E55","",100,-1,1200));
        x.CreatPlot(F51,new TH1D("F51","",100,0.2,1.1));
        x.CreatPlot(ABCD,new TH2D("ABCD","",100,0,1.1,100,0,30));

  }
  ~MonoCuts(){}
  void doAnalysis(vector<MonoCandidate> &cand, unsigned nCandidates,int TRG, unsigned ev)
  {
	CutFlowCand_TRG.clear();	
	CutFlowCand_Qual.clear();
	CutFlowCand_Energy.clear();
	CutFlowCand_F51.clear();
	CutFlowCand_Dedx.clear();
	N1CutCand_TRG.clear();	
//	N1CutCand_Qual.clear();
//	N1CutCand_Energy.clear();
//	N1CutCand_F51.clear();
//	N1CutCand_Dedx.clear();
	N1CutCand_HIso.clear();
	N1CutCand_Dist.clear();
	N1CutCand_XYPar0.clear();
	N1CutCand_XYPar1.clear();
	N1CutCand_XYPar2.clear();
	N1CutCand_RZPar0.clear();
	N1CutCand_RZPar1.clear();
	N1CutCand_RZPar2.clear();
	CutFlowCand_Dist.clear();
	CutFlowCand_HIso.clear();
	CutFlowCand_XYPar0.clear();
	CutFlowCand_XYPar2.clear();
	CutFlowCand_XYPar1.clear();
	CutFlowCand_RZPar0.clear();
	CutFlowCand_RZPar1.clear();
	CutFlowCand_RZPar2.clear();

	cout<<"ev "<<ev<<", cand number "<<nCandidates<<endl;
     for(unsigned c=0;c<nCandidates;c++){
	MonoCandidate &cands = cand[c];
	
	bool QualCut = evalQuality(cands);
	bool HIsoCut   = evalHIso(cands);
	bool DistCut   = evalDist(cands);
	bool XYPar0Cut = evalXYPar0(cands);
	bool XYPar1Cut = evalXYPar1(cands);
	bool XYPar2Cut = evalXYPar2(cands);
	bool RZPar0Cut = evalRZPar0(cands);
	bool RZPar1Cut = evalRZPar1(cands);
	bool RZPar2Cut = evalRZPar2(cands);
	bool ECut = evalE(cands);
	bool F51Cut = evalF51(cands);
	bool dEdXCut = evaldEdX(cands);

	//count for total events without TRG	


	if( DistCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar0Cut && RZPar1Cut && RZPar2Cut)N1CutCand_HIso.push_back(cands);  
	if( HIsoCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar0Cut && RZPar1Cut && RZPar2Cut)N1CutCand_Dist  .push_back(cands); 
	if( HIsoCut && DistCut && XYPar1Cut && XYPar2Cut && RZPar0Cut && RZPar1Cut && RZPar2Cut)  N1CutCand_XYPar0.push_back(cands);
	if( HIsoCut && DistCut && XYPar0Cut && XYPar2Cut && RZPar0Cut && RZPar1Cut && RZPar2Cut)  N1CutCand_XYPar1.push_back(cands);
	if( HIsoCut && DistCut && XYPar0Cut && XYPar1Cut && RZPar0Cut && RZPar1Cut && RZPar2Cut)  N1CutCand_XYPar2.push_back(cands);
	if( HIsoCut && DistCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar1Cut && RZPar2Cut)  N1CutCand_RZPar0.push_back(cands);
	if( HIsoCut && DistCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar0Cut && RZPar2Cut)  N1CutCand_RZPar1.push_back(cands);
	if( HIsoCut && DistCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar0Cut && RZPar1Cut)  N1CutCand_RZPar2.push_back(cands);
	//N-1 cut and relative efficiency
//        if( ECut && F51Cut && dEdXCut &&TRG) N1CutCand_Qual.push_back(cands); 
//        if( QualCut && F51Cut && dEdXCut &TRG) N1CutCand_Energy.push_back(cands);
//        if( QualCut && ECut && dEdXCut &&TRG ) N1CutCand_F51.push_back(cands);
//	if( QualCut && ECut && F51Cut &&TRG)  N1CutCand_Dedx.push_back(cands);
//	if( QualCut && ECut & F51Cut && dEdXCut) N1CutCand_TRG.push_back(cands);
	
	//-----------------------------------------------------------------	
	//---Cutflow histograms--------------------------------------------
	//-----------------------------------------------------------------

	//signal efficiency
	
	if(TRG) CutFlowCand_TRG.push_back(cands);

	if(TRG && DistCut ) 								CutFlowCand_Dist.push_back(cands); 
	if(TRG && DistCut && HIsoCut) 							CutFlowCand_HIso.push_back(cands); 
	if(TRG && DistCut && HIsoCut && XYPar0Cut) 					CutFlowCand_XYPar0.push_back(cands); 
	if(TRG && DistCut && HIsoCut && XYPar0Cut && XYPar1Cut) 			CutFlowCand_XYPar1.push_back(cands); 
	if(TRG && DistCut && HIsoCut && XYPar0Cut && XYPar1Cut && XYPar2Cut) 		CutFlowCand_XYPar2.push_back(cands); 
	if(TRG && DistCut && HIsoCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar0Cut) CutFlowCand_RZPar0.push_back(cands); 
	if(TRG && DistCut && HIsoCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar0Cut 
	       && RZPar1Cut) 								CutFlowCand_RZPar1.push_back(cands); 
	if(TRG && DistCut && HIsoCut && XYPar0Cut && XYPar1Cut && XYPar2Cut && RZPar0Cut 
	       && RZPar1Cut && RZPar2Cut) 						CutFlowCand_RZPar2.push_back(cands); 
	if(TRG && QualCut && ECut ) 							CutFlowCand_Energy.push_back(cands);
	if(TRG && QualCut && ECut && F51Cut) 						CutFlowCand_F51.push_back(cands);
	if(TRG && QualCut && ECut && F51Cut && dEdXCut) 				CutFlowCand_Dedx.push_back(cands);
	
      }//for cand loop

    	
        sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	
	if(CutFlowCand_TRG.size()>0) 
	{
		count++;
	}
        sort(CutFlowCand_Dist.begin(),CutFlowCand_Dist.begin()+CutFlowCand_Dist.size());
	
	if(CutFlowCand_Dist.size()>0) 
	{
		Dist_count++;
	}

        sort(CutFlowCand_HIso.begin(),CutFlowCand_HIso.begin()+CutFlowCand_HIso.size());
	if(CutFlowCand_HIso.size()>0) 
	{
		HIso_count++;	
	}
        sort(CutFlowCand_XYPar0.begin(),CutFlowCand_XYPar0.begin()+CutFlowCand_XYPar0.size());
	if(CutFlowCand_XYPar0.size()>0)
	{
		XYPar0_count++;	
	}
        sort(CutFlowCand_XYPar1.begin(),CutFlowCand_XYPar1.begin()+CutFlowCand_XYPar1.size());
	if(CutFlowCand_XYPar1.size()>0)
	{
		XYPar1_count++;	
	}
        sort(CutFlowCand_XYPar2.begin(),CutFlowCand_XYPar2.begin()+CutFlowCand_XYPar2.size());
	if(CutFlowCand_XYPar2.size()>0)
	{
		XYPar2_count++;	
	}
	
	sort(CutFlowCand_RZPar0.begin(),CutFlowCand_RZPar0.begin()+CutFlowCand_RZPar0.size());
	if(CutFlowCand_RZPar0.size()>0)
	{
		RZPar0_count++;	
	}
	sort(CutFlowCand_RZPar1.begin(),CutFlowCand_RZPar1.begin()+CutFlowCand_RZPar1.size());
	if(CutFlowCand_RZPar1.size()>0)
	{
		RZPar1_count++;	
	}
	sort(CutFlowCand_RZPar2.begin(),CutFlowCand_RZPar2.begin()+CutFlowCand_RZPar2.size());
	if(CutFlowCand_RZPar2.size()>0)
	{
		RZPar2_count++;	
	}
	sort(CutFlowCand_Energy.begin(),CutFlowCand_Energy.begin()+CutFlowCand_Energy.size());
	if(CutFlowCand_Energy.size()>0) E_count++;
	sort(CutFlowCand_F51.begin(),CutFlowCand_F51.begin()+CutFlowCand_F51.size());
	if(CutFlowCand_F51.size()>0){

		//PrintInfo();	
		 f51_count++;
	}
	sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());
	if(CutFlowCand_Dedx.size()>0)
	{
		//PrintInfo();	
		Dedx_count++;	
	}


  }
  bool operator<(const MonoCandidate &mc)const{
   if(dEdXSig_>mc.dEdXSig_) return true;
   else if(dEdXSig_==mc.dEdXSig_){
        if(f51_>mc.f51_) return true;
        else return false;
        }
    else return false;
  }
  void WritePlots(TFile *oFile){
        oFile->cd(trigName_.c_str());
	        Plot[0].WritePlot();
   }
  
  void PrintInfo(){	
		for(int i=0;i<CutFlowCand_F51.size();i++){
			cout<<i<<endl;
			cout<<"     eta "<<setprecision(5)<<CutFlowCand_F51[i].eta_<<endl;
			cout<<"      e55 "<<setprecision(5)<<CutFlowCand_F51[i].e55_<<endl;
			cout<<"     f51 "<<setprecision(5)<<CutFlowCand_F51[i].f51_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<CutFlowCand_F51[i].dEdXSig_<<endl;
			cout<<"----------------------------"<<endl;
		}
  }
  void SignalEff(const string trName)
  {
        //signal efficiency = no. of events after all selection cuts/all events
	cout<<trName<<" ================================="<<endl;
	cout<<"        TRG "<<count<<endl;
	cout<<"       Dist "<<Dist_count<<endl;
	cout<<"       HIso "<<HIso_count<<endl;
	cout<<"     XYPar0 "<<XYPar0_count<<endl;
	cout<<"     XYPar1 "<<XYPar1_count<<endl;
	cout<<"     XYPar2 "<<XYPar2_count<<endl;
	cout<<"     RZPar0 "<<RZPar0_count<<endl;
	cout<<"     RZPar1 "<<RZPar1_count<<endl;
	cout<<"     RZPar2 "<<RZPar2_count<<endl;
        cout<<"       ECut "<<E_count<<endl;
        cout<<"     F51Cut "<<f51_count<<endl;
        cout<<" dEdXSigCut "<<Dedx_count<<endl;
        cout<<"------------------------------------------"<<endl;
	cout<<"Dist "<<(double)Dedx_count/(double)Dist_count<<endl;
	cout<<"HIso "<<(double)Dedx_count/(double)(HIso_count)<<endl;
	cout<<"XYPar0 "<<(double)Dedx_count/(double)(XYPar0_count)<<endl;
	cout<<"XYPar1 "<<(double)Dedx_count/(double)(XYPar1_count)<<endl;
	cout<<"XYPar2 "<<(double)Dedx_count/(double)(XYPar2_count)<<endl;
	cout<<"RZPar0 "<<(double)Dedx_count/(double)(RZPar0_count)<<endl;
	cout<<"RZPar1 "<<(double)Dedx_count/(double)(RZPar1_count)<<endl;
	cout<<"RZPar2 "<<(double)Dedx_count/(double)(RZPar2_count)<<endl;
	cout<<"Energy "<<(double)Dedx_count/(double)E_count<<endl;
	cout<<"F51    "<<(double)Dedx_count/(double)f51_count<<endl;
        cout<<"Signal efficiency = "<<(double)Dedx_count/(double)count<<endl;
        cout<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"Relative efficiency"<<endl;
        cout<<"------------------------------------------"<<endl;
	cout<<"	   TRG "<<(double)RZPar2_count/(double)NoTRG<<endl;
	cout<<"	  Dist "<<(double)RZPar2_count/(double)NoDist<<endl;
	cout<<"	  HIso "<<(double)RZPar2_count/(double)NoHIso<<endl;
	cout<<"	XYPar0 "<<(double)RZPar2_count/(double)NoXYPar0<<endl;
	cout<<"	XYPar1 "<<(double)RZPar2_count/(double)NoXYPar1<<endl;
	cout<<"	XYPar2 "<<(double)RZPar2_count/(double)NoXYPar2<<endl;
	cout<<"	RZPar0 "<<(double)RZPar2_count/(double)NoRZPar0<<endl;
	cout<<"	RZPar1 "<<(double)RZPar2_count/(double)NoRZPar1<<endl;
	cout<<"	RZPar2 "<<(double)RZPar2_count/(double)NoRZPar2<<endl;
//        cout<<"        TRG "<<(double)dEdX_count/(double)NoTRG<<endl;
//        cout<<"QualityCuts "<<(double)dEdX_count/(double)NoQual<<endl;
//        cout<<"       ECut "<<(double)dEdX_count/(double)NoE<<endl;
//        cout<<"     F51Cut "<<(double)dEdX_count/(double)NoF51<<endl;
//        cout<<" dEdXSigCut "<<(double)dEdX_count/(double)NodEdXCut<<endl;
	  
      cout<<endl;
  }

  //Cuts parameters
  static const double xyp0Cut_ ;
  static const double xyp1Cut_ ;
  static const double xyp2Cut_ ;
  static const double rzp0Cut_ ;
  static const double rzp1Cut_ ;
  static const double rzp2Cut_ ;
  static const double distCut_ ;
  static const double hIsoCut_ ;
  static const double dEdXSigCut_ ;
  static const double e55Cut_ ;
  static const double f51Cut_ ;
  static const double photonCut_ ;
   // if you want to set parameter in the class, you should add constexpr
  //  ex. constexpr static const double x=1;

	
private:
  string trigName_;
  vector<PlotSet> Plot;
  //N-1 cut plot box
  vector<PlotSet> n_1Plot;
  //cutflow plot box
  vector<PlotSet> CutFlow;
  static const unsigned nCut = 9U;
  string cutName_[nCut];
  


  //cuts analysis
  bool evalQuality(MonoCandidate &mc) { return TMath::Abs(mc.xyp0_ )< xyp0Cut_ && TMath::Abs(mc.xyp1_ )< xyp1Cut_
		&& TMath::Abs(mc.xyp2_) > xyp2Cut_ 
		&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_  &&  
		TMath::Abs(mc.rzp2_) < rzp2Cut_ && TMath::Abs(mc.rzp1_)<rzp1Cut_ && TMath::Abs(mc.rzp0_) < rzp0Cut_;  }
  bool evalDist(MonoCandidate &mc) { return TMath::Abs(mc.dist_) < distCut_ ;}
  bool evalHIso(MonoCandidate &mc) { return mc.hIso_ < hIsoCut_ ;}
  bool evalXYPar0(MonoCandidate &mc ){ return TMath::Abs(mc.xyp0_) < xyp0Cut_ ;}
  bool evalXYPar1(MonoCandidate &mc ){ return TMath::Abs(mc.xyp1_) < xyp1Cut_ ;}

  bool evalXYPar2(MonoCandidate &mc ){ return TMath::Abs(mc.xyp2_) > xyp2Cut_ ;}

  bool evalRZPar0(MonoCandidate &mc ){ return TMath::Abs(mc.rzp0_) < rzp0Cut_ ;}
  bool evalRZPar1(MonoCandidate &mc ){ return TMath::Abs(mc.rzp1_) < rzp1Cut_ ;}
  bool evalRZPar2(MonoCandidate &mc) { return TMath::Abs(mc.rzp2_) < rzp2Cut_;  }
  bool evalE(MonoCandidate &mc) { return mc.e55_ > e55Cut_; }
  bool evalF51(MonoCandidate &mc) { return mc.f51_ > f51Cut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
 
  vector<MonoCandidate> CutFlowCand_TRG;
  vector<MonoCandidate> CutFlowCand_Qual;
  vector<MonoCandidate> CutFlowCand_Energy;
  vector<MonoCandidate> CutFlowCand_F51;
  vector<MonoCandidate> CutFlowCand_Dedx;

  vector<MonoCandidate> CutFlowCand_Dist;
  vector<MonoCandidate> CutFlowCand_HIso;
  vector<MonoCandidate> CutFlowCand_XYPar0;
  vector<MonoCandidate> CutFlowCand_XYPar2;
  vector<MonoCandidate> CutFlowCand_XYPar1;
  vector<MonoCandidate> CutFlowCand_RZPar0;
  vector<MonoCandidate> CutFlowCand_RZPar1;
  vector<MonoCandidate> CutFlowCand_RZPar2;

  vector<MonoCandidate> N1CutCand_Dist  ; 
  vector<MonoCandidate> N1CutCand_HIso  ; 
  vector<MonoCandidate> N1CutCand_XYPar0;
  vector<MonoCandidate> N1CutCand_XYPar1;
  vector<MonoCandidate> N1CutCand_XYPar2;
  vector<MonoCandidate> N1CutCand_RZPar0;
  vector<MonoCandidate> N1CutCand_RZPar1;
  vector<MonoCandidate> N1CutCand_RZPar2;


  vector<MonoCandidate> N1CutCand_TRG;
  vector<MonoCandidate> N1CutCand_Qual;
  vector<MonoCandidate> N1CutCand_Energy;
  vector<MonoCandidate> N1CutCand_F51;
  vector<MonoCandidate> N1CutCand_Dedx;
	
  //no. of every selections 
  int Qual=0;
  int E=0;
  int f51=0;
  int dEdX=0;
  //count for cutflow 
  int count=0;
  int Qual_count=0;
  int E_count=0;
  int f51_count=0;
  int Dedx_count=0; 
  int Dist_count =0 ;
  int HIso_count = 0;
  int XYPar0_count = 0;
  int XYPar1_count = 0;
  int XYPar2_count = 0;
  int RZPar0_count = 0;
  int RZPar1_count = 0;
  int RZPar2_count = 0;
 
  int Reco=0;
  int Gen=0;
  int photonLike=0;

  //relative eff without HLT
  int NoTRG=0;
  int NoQual=0;
  int NoE=0;
  int NoF51=0;
  int NodEdXCut=0;

  int NoDist  =0; 
  int NoHIso  =0;
  int NoXYPar0=0;
  int NoXYPar1=0;
  int NoXYPar2=0;
  int NoRZPar0=0;
  int NoRZPar1=0;
  int NoRZPar2=0;
  //Apply TRG 
  int count_TRG=0;
  int Qual_count_TRG=0;
  int E_count_TRG=0;
  int f51_count_TRG=0;
  int dEdX_count_TRG=0;
  //relative eff with HLT
  int No_TRGTRG=0;
  int No_TRGQual=0;
  int No_TRGE=0;
  int No_TRGF51=0;
  int No_TRGdEdXCut=0;

};

  //Cuts parameters
  const double MonoCuts::xyp0Cut_=0.6;//d0
  const double MonoCuts::xyp1Cut_=10;//phi0
  const double MonoCuts::xyp2Cut_=1000;//radius of the circle
  const double MonoCuts::rzp0Cut_=20;
  const double MonoCuts::rzp1Cut_=2.5;
  const double MonoCuts::rzp2Cut_=0.005;
  const double MonoCuts::distCut_ = 0.5;
  const double MonoCuts::hIsoCut_= 10;
  const double MonoCuts::dEdXSigCut_ = 0;
  const double MonoCuts::e55Cut_ = 200;
  const double MonoCuts::f51Cut_ = 0;

void DataAnalyzer()
{
	TFile *oFile = new TFile("BlindAnalysis.root","recreate");
	TFile *fin = new TFile("MonopoleData2018C_Blinded.root");
        TTree *Tree = (TTree*)fin->Get("monopoles");

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
	Tree->SetBranchAddress("run",&run);
	Tree->SetBranchAddress("passHLT_Photon200", &passHLT_Photon200);
	
	vector<MonoCandidate> cand(10);
	MonoCuts TrgAnalysis("HLT_Photon",oFile);
	int nParticle = 0;
        for(int ev=0; ev<Tree->GetEntries();ev++){
		Tree->GetEntry(ev);
		nParticle = dEdxSig->size();
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
		TrgAnalysis.doAnalysis(cand,nParticle,passHLT_Photon200,ev);
	}
	TrgAnalysis.SignalEff("HLT_Photon200");
	TrgAnalysis.WritePlots(oFile);
	oFile->Close();

}
