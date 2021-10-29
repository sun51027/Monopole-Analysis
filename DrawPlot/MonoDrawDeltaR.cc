//////////////////////////////////
//
//	File name: MonoDrawDeltaR.cc
//	author: Shih Lin
//	Purpose: To see the deltaR between candidate and gen
//		 decide the deltaR cut 
//		 (July 16 2021)
/////////////////////////////////////////////////////
#include "iostream"
#include "TAttMarker.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "math.h"
#include <algorithm>
#include <string>
using namespace std;
enum PlotName{
  FracSatVNstrips = 0, // fraction of saturated strips vs. number of strips
  DedXSig,             // dE/dX significance
  RZcurv,              // RZ curvature/uncertainty
//  Energy,
  E55,
  F51,                 // frac 51
  HcalIso,             // Hcal Iso
  ABCD,
  NoSelectionDeltaR,
  SelectedDeltaR
};
static const unsigned nPlot = 10U;
class PlotSet
{
public:
  PlotSet(){plots_.resize(nPlot);}
  ~PlotSet(){}
  void CreatPlot(const PlotName pn, TH1* h){ 
	plots_[pn] = h;
  }
  TH1 * GetPlot(const PlotName pn){ return plots_[pn]; }
  void Print(){ 
	cout<<plots_[0]<<endl;
	cout<<plots_[1]<<endl;}
  void WritePlot(){
	for(int pn=0;pn<nPlot;pn++){
	TH1 *h = plots_[pn];
	if(h){ 
	  h->Write();
	//  cout<<h<<"pass h->Write()"<<endl;
	  }
	}
  }
private:
  vector<TH1*> plots_;
};
class MonoCandidate
{
public:
  MonoCandidate(){}
  //This will be used in main function and to absord the data in root file
    MonoCandidate(double sh, double satsh, double dedxsig,double tiso, double xyp0, double xyp1, double xyp2,
    double rzp0, double rzp1, double rzp2,
    double dist, double f51, double f15,
    double e55, double hiso, double eta,
    double phi, double mono_eta, double mono_phi, double amon_eta, double amon_phi, 
    double event,double NPV):
  subHits_(sh),subSatHits_(satsh),dEdXSig_(dedxsig),tIso_(tiso),xyp0_(xyp0),
  xyp1_(xyp1),xyp2_(xyp2),rzp0_(rzp0),rzp1_(rzp1),rzp2_(rzp2),
  dist_(dist),f51_(f51),f15_(f15),e55_(e55),hIso_(hiso),
  eta_(eta),phi_(phi),mono_eta_(mono_eta), mono_phi_(mono_phi),
  amon_eta_(amon_eta), amon_phi_(amon_phi),event_(event),NPV_(NPV) { }
  //This will be used in comparing with cut
  MonoCandidate(const MonoCandidate &mc) : 
    subHits_(mc.subHits_),subSatHits_(mc.subSatHits_),dEdXSig_(mc.dEdXSig_),tIso_(mc.tIso_),
    xyp0_(mc.xyp0_),xyp1_(mc.xyp1_),xyp2_(mc.xyp2_),
    rzp0_(mc.rzp0_),rzp1_(mc.rzp1_),rzp2_(mc.rzp2_),
    dist_(mc.dist_),f51_(mc.f51_),f15_(mc.f15_),e55_(mc.e55_),
    hIso_(mc.hIso_),eta_(mc.eta_),phi_(mc.phi_),mono_eta_(mc.mono_eta_), mono_phi_(mc.mono_phi_),
  amon_eta_(mc.amon_eta_), amon_phi_(mc.amon_phi_),event_(mc.event_),NPV_(mc.NPV_) { } 
        
  ~MonoCandidate() {}
  void print(){
  } 
  bool operator<(const MonoCandidate &mc)const{
   if(dEdXSig_>mc.dEdXSig_) return true;
   else if(dEdXSig_==mc.dEdXSig_){
        if(f51_>mc.f51_) return true;
        else return false;
        }
    else return false;
  }
  //All candidates variable
  double subHits_;
  double subSatHits_;
  double dEdXSig_;
  double tIso_;
  double xyp0_;
  double xyp1_;
  double xyp2_;
  double rzp0_;
  double rzp1_;
  double rzp2_;

  double dist_;
  double f51_;
  double f15_;
  double e55_;
  double hIso_;
  double eta_;
  double phi_;
  double event_;
  double NPV_;
  double mono_eta_;
  double mono_phi_;
  double amon_eta_;
  double amon_phi_;

};	
class Photon
{
public:
  Photon(){}
  Photon(double pho_eta, double pho_phi, double pho_pt):pho_eta_(pho_eta),pho_phi_(pho_phi),pho_pt_(pho_pt){}
  Photon(const Photon &mc):pho_eta_(mc.pho_eta_),pho_phi_(mc.pho_phi_),pho_pt_(mc.pho_pt_){}
  double pho_eta_;
  double pho_phi_;
  double pho_pt_;

  ~Photon(){}

};
class MonoCuts:public MonoCandidate, public Photon
{
public:
  MonoCuts(){}

  MonoCuts(const string &trName,TFile *openFile):trigName_(trName){
	
	//create a new directory in the output file
	openFile->mkdir(trName.c_str());
	cout<<"create directory "<<trName.c_str()<<endl;	
	//create histrogram N-1 cutf for No trigger
	Plot.resize(1U);
	PlotSet &x = Plot[0];


        x.CreatPlot(NoSelectionDeltaR,new TH1D("NoSelectionDeltaR","",100,0,1));
        x.CreatPlot(SelectedDeltaR,new TH1D("SelectedDeltaR","",100,0,1));



//        x.CreatPlot(FracSatVNstrips,new TH2D("FracSatVNstrips","",100,0,1000,100,0,1));
//        x.CreatPlot(DedXSig,new TH1D("DedXSig","",100,0,30));
//        x.CreatPlot(RZcurv,new TH1D("RZcurv","",100,-0.01,0.01));
////	x.CreatPlot(Energy,new TH1D("Energy","",100,0,1400));
//	x.CreatPlot(E55,new TH1D("E55","",100,-1,1200));
//        x.CreatPlot(F51,new TH1D("F51","",100,0.2,1.1));
//        x.CreatPlot(HcalIso,new TH1D("HcalIso","",100,-1,10));
//        x.CreatPlot(ABCD,new TH2D("ABCD","",100,0,1.1,100,0,30));
	
	cout<<"the address of plots x "<<&x<<endl;

  }

  ~MonoCuts(){}
  void PrintFakeEvent(unsigned ev)
  {
	cout<<"Fake ev "<<ev<<endl;
  }
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev)
  {
//	cout<<"enter doAnalysis successfully"<<endl;
//	cout<<"nPhoton "<<nPhoton<<endl;
	CutFlowCand_TRG.clear();	
	CutFlowCand_Qual.clear();
	CutFlowCand_Energy.clear();
	CutFlowCand_F51.clear();
	CutFlowCand_Dedx.clear();
	N1CutCand_TRG.clear();	
	N1CutCand_Qual.clear();
	N1CutCand_Energy.clear();
	N1CutCand_F51.clear();
	N1CutCand_Dedx.clear();
	HighPtPhoton.clear();
	
//	cout<<"ev "<<ev<<", cand number "<<nCandidates<<endl;
     for(unsigned c=0;c<nCandidates;c++){
	MonoCandidate &cands = cand[c];
	bool QualCut = evalQuality(cands);
	bool ECut = evalE(cands);
	bool F51Cut = evalF51(cands);
	bool dEdXCut = evaldEdX(cands);

	//count for total events without TRG	

	if(TRG) CutFlowCand_TRG.push_back(cands);
	
	//N-1 cut and relative efficiency
        if( ECut && F51Cut && dEdXCut &&TRG) N1CutCand_Qual.push_back(cands); 
        if( QualCut && F51Cut && dEdXCut &TRG) N1CutCand_Energy.push_back(cands);
        if( QualCut && ECut && dEdXCut &&TRG ) N1CutCand_F51.push_back(cands);
	if( QualCut && ECut && F51Cut &&TRG)  N1CutCand_Dedx.push_back(cands);
	if( QualCut && ECut & F51Cut && dEdXCut) N1CutCand_TRG.push_back(cands);
	
	//-----------------------------------------------------------------	
	//---Cutflow histograms--------------------------------------------
	//-----------------------------------------------------------------

	//signal efficiency

	if(TRG && QualCut ) CutFlowCand_Qual.push_back(cands); 
	if(TRG && QualCut && ECut ) CutFlowCand_Energy.push_back(cands);
	if(TRG && QualCut && ECut && F51Cut) CutFlowCand_F51.push_back(cands);
	if(TRG && QualCut && ECut && F51Cut && dEdXCut) CutFlowCand_Dedx.push_back(cands);
	
      }//for cand loop

//	cout<<"entering HighPtPhoton loop"<<endl;	
    	
    if(nPhoton!=0){
    for(unsigned p=0;p<nPhoton;p++){
	Photon &photon = pho[p];
	bool PhotonPtCut = evalPhoton(photon);
	//HighPtPhoton.push_back(photon);
	if(PhotonPtCut) HighPtPhoton.push_back(photon);
    }
    }
	sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());

                PlotSet &x = Plot[0];
	   for(int i=0; i<CutFlowCand_TRG.size();i++){

		double m_deltaR=0;
		double am_deltaR=0;
		m_deltaR = sqrt(pow(CutFlowCand_TRG[i].eta_-CutFlowCand_TRG[0].mono_eta_,2)+
				pow(CutFlowCand_TRG[i].phi_-CutFlowCand_TRG[0].mono_phi_,2));
		am_deltaR= sqrt(pow(CutFlowCand_TRG[i].eta_-CutFlowCand_TRG[0].amon_eta_,2)+
                                pow(CutFlowCand_TRG[i].phi_-CutFlowCand_TRG[0].amon_phi_,2));


		x.GetPlot(NoSelectionDeltaR)->Fill(m_deltaR);
		x.GetPlot(NoSelectionDeltaR)->Fill(am_deltaR);
		
	}


	sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());
	if(CutFlowCand_Dedx.size()==1) 	    MonoNum1++;
	else if(CutFlowCand_Dedx.size()==2) MonoNum2++;
	else if(CutFlowCand_Dedx.size()==3) MonoNum3++;
	   for(int i=0; i<CutFlowCand_Dedx.size();i++){

		double m_deltaR=0;
		double am_deltaR=0;
		m_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].mono_eta_,2)+
				pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].mono_phi_,2));
		am_deltaR= sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].amon_eta_,2)+
                                pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].amon_phi_,2));


		x.GetPlot(SelectedDeltaR)->Fill(m_deltaR);
		x.GetPlot(SelectedDeltaR)->Fill(am_deltaR);
	}	
  }
	//Cut all
/*	if(CutFlowCand_Dedx.size()>0)
//	if(RealMono==true && RealAnti ==true)
	{
		        
		PlotSet &x = CutFlow[3];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_Dedx[0].subHits_,CutFlowCand_Dedx[0].subSatHits_/CutFlowCand_Dedx[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_Dedx[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_Dedx[0].rzp2_);
        	x.GetPlot(E55)->Fill(CutFlowCand_Dedx[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_Dedx[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_Dedx[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_Dedx[0].f51_,CutFlowCand_Dedx[0].dEdXSig_);
		dEdX_count++;	
	}

*/

  //}
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
/*	for(unsigned c=0; c<nCut; c++) n_1Plot[c].WritePlot();
        cout<<"n-1cut pass writeplots func"<<endl;

	for(unsigned c=0; c<nCut; c++) CutFlow[c].WritePlot();
	cout<<"cutflow pass writeplots func"<<endl;
*/	
  }
  
  void SignalEff(const string trName)
  {
        //signal efficiency = no. of events after all selection cuts/all events
	cout<<trName<<" ================================="<<endl;
//        cout<<"No TRG ==================================="<<endl;
        cout<<"        TRG "<<count<<endl;
        cout<<"QualityCuts "<<Qual_count<<endl;
        cout<<"       ECut "<<E_count<<endl;
        cout<<"     F51Cut "<<f51_count<<endl;
        cout<<" dEdXSigCut "<<dEdX_count<<endl;
	cout<<" Ecut lose% "<<(double)(Qual_count-E_count)/(double)Qual_count<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"Relative effciency count------------------"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"     No TRG "<<NoTRG<<endl;
        cout<<" No Quality "<<NoQual<<endl;
        cout<<"    No ECut "<< NoE <<endl;
        cout<<"  No F51Cut "<<NoF51<<endl;
        cout<<" No dEdXSig "<<NodEdXCut<<endl;
        cout<<"=========================================="<<endl;
        cout<<"Signal efficiency = "<<(double)dEdX_count/(double)TotalEvents<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"Relative efficiency-----------------------"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"        TRG "<<(double)dEdX_count/(double)NoTRG<<endl;
        cout<<"QualityCuts "<<(double)dEdX_count/(double)NoQual<<endl;
        cout<<"       ECut "<<(double)dEdX_count/(double)NoE<<endl;
        cout<<"     F51Cut "<<(double)dEdX_count/(double)NoF51<<endl;
        cout<<" dEdXSigCut "<<(double)dEdX_count/(double)NodEdXCut<<endl;
	cout<<"///////////////////////////////////////"<<endl;
//	cout<<count_test<<endl;
 	cout<<" 1 monopole event "<<MonoNum1<<endl;
 	cout<<" 2 monopole event "<<MonoNum2<<endl;
 	cout<<" 3 monopole event "<<MonoNum3<<endl;
	cout<<" ideal total mono "<<MonoNum1+2*MonoNum2<<endl;
        cout<<"Real monopole "<<RealNum<<endl;
        cout<<"Real antimono "<<RealAntiNum<<endl;
	cout<<"Fake monopole "<<FakeNum<<endl;
        cout<<"Total reconstructed event "<<Reco<<endl;
	cout<<"photonLike "<<photonLike<<endl;
	  
      cout<<endl;
  }

  //Cuts parameters
  static const double xyp0Cut_ ;
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
  static const unsigned nCut = 5U;
  string cutName_[nCut];
  


  //cuts analysis
  bool evalQuality(MonoCandidate &mc) { return mc.xyp0_ < xyp0Cut_&& mc.xyp2_ > xyp2Cut_ 
			&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_;  }
  bool evalE(MonoCandidate &mc) { return mc.e55_ > e55Cut_; }
  bool evalF51(MonoCandidate &mc) { return mc.f51_ > f51Cut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
  int TotalEvents = 15000;
 
  vector<MonoCandidate> CutFlowCand_TRG;
  vector<MonoCandidate> CutFlowCand_Qual;
  vector<MonoCandidate> CutFlowCand_Energy;
  vector<MonoCandidate> CutFlowCand_F51;
  vector<MonoCandidate> CutFlowCand_Dedx;

  vector<MonoCandidate> N1CutCand_TRG;
  vector<MonoCandidate> N1CutCand_Qual;
  vector<MonoCandidate> N1CutCand_Energy;
  vector<MonoCandidate> N1CutCand_F51;
  vector<MonoCandidate> N1CutCand_Dedx;
	
  vector<Photon> HighPtPhoton;
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
  int dEdX_count=0; 
  
  int MonoNum1=0;
  int MonoNum2=0;
  int MonoNum3=0;
  int FakeNum=0;
  int RealNum=0;
  int RealAntiNum=0;
  int RealEvNum=0;
  
  int Reco=0;
  int Gen=0;
  int photonLike=0;

  int count_test=0;
  int TRGminitest=0;
  int Qualminitest=0;
  //relative eff without HLT
  int NoTRG=0;
  int NoQual=0;
  int NoQualCand=0;
  int NoE=0;
  int NoECand=0;
  int NoF51=0;
  int NoF51Cand=0;
  int NodEdXCut=0;
  int NodEdXCutCand=0;
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
  const double MonoCuts::xyp0Cut_=0.6;
  const double MonoCuts::xyp2Cut_=1000;
  const double MonoCuts::rzp0Cut_=10;
  const double MonoCuts::rzp1Cut_=999;
  const double MonoCuts::rzp2Cut_=0.005;
  const double MonoCuts::distCut_ = 0.5;
  const double MonoCuts::hIsoCut_= 10;
//  const double MonoCuts::dEdXSigCut_ = 7;//loose
  const double MonoCuts::dEdXSigCut_ = 9;
//  const double MonoCuts::e55Cut_ = 0;
  const double MonoCuts::e55Cut_ = 200;
//  const double MonoCuts::f51Cut_ = 0.6;
  const double MonoCuts::f51Cut_ = 0.85;
  const double MonoCuts::photonCut_ = 200;

void MonoDrawDeltaR()
{
	cout<<"hello"<<endl;
	TFile *oFile = new TFile("Monopole_DeltaR.root","recreate");
	cout<<"new file line pass"<<endl;

	TFile *fin = new TFile("MonoNtuple2018_MC_1000.root");
	cout<<"open file success"<<endl;
        TTree *tree = (TTree*)fin->Get("monopoles");
	cout<<"open tree success"<<endl;
        Bool_t passHLT_Photon200;
	Bool_t passHLT_Photon175;
	Bool_t passHLT_DoublePhoton70;
        vector<bool> * trigResults = 0;
        vector<string> * trigNames = 0;

        unsigned nCandidates;
	

        vector<double> *subHits=0;
        vector<double> *subSatHits=0;
        vector<double> *dEdXSig=0;
        vector<double> * tIso = 0;
        vector<double> * xyp0 = 0;
        vector<double> * xyp1 = 0;
        vector<double> * xyp2 = 0;
        vector<double> * rzp0 = 0;
        vector<double> * rzp1 = 0;
        vector<double> * rzp2 = 0;

        vector<double> * dist = 0;
        vector<double> * f51 = 0;
        vector<double> * f15 = 0;
        vector<double> * e55 = 0;
        vector<double> * hIso = 0;
        vector<double> * eta = 0;
        vector<double> * phi = 0;
	vector<double> * pho_eta = 0;
	vector<double> * pho_phi = 0;
	vector<double> * pho_pt = 0;
        unsigned nPhoton;

	vector<double> * testE = 0;
	unsigned event;
	unsigned NPV;
	
        double mono_eta;
        double mono_phi;
        double amon_eta;
        double amon_phi;
//	double pho_eta;
//	double pho_phi;
//	double pho_pt;

	tree->SetBranchAddress("event",&event);
	tree->SetBranchAddress("trigResult",&trigResults);
        tree->SetBranchAddress("trigNames",&trigNames);
        tree->SetBranchAddress("passHLT_Photon200" , &passHLT_Photon200);
	tree->SetBranchAddress("passHLT_Photon175" , &passHLT_Photon175);
        tree->SetBranchAddress("passHLT_DoublePhoton70",&passHLT_DoublePhoton70);
	tree->SetBranchAddress("cand_N",&nCandidates);
        tree->SetBranchAddress("cand_SubHits",&subHits);
        tree->SetBranchAddress("cand_SatSubHits",&subSatHits);
        tree->SetBranchAddress("cand_dEdXSig",&dEdXSig);
        tree->SetBranchAddress("cand_TIso",&tIso);
        tree->SetBranchAddress("cand_f51",&f51);
        tree->SetBranchAddress("cand_f15",&f15);
        tree->SetBranchAddress("cand_e55",&e55);
        tree->SetBranchAddress("cand_HIso",&hIso);
        tree->SetBranchAddress("cand_XYPar0",&xyp0);
        tree->SetBranchAddress("cand_XYPar1",&xyp1);
        tree->SetBranchAddress("cand_XYPar2",&xyp2);
        tree->SetBranchAddress("cand_RZPar0",&rzp0);
        tree->SetBranchAddress("cand_RZPar1",&rzp1);
        tree->SetBranchAddress("cand_RZPar2",&rzp2);
        tree->SetBranchAddress("cand_eta",&eta);
        tree->SetBranchAddress("cand_phi",&phi);
        tree->SetBranchAddress("cand_dist",&dist);
        tree->SetBranchAddress("mono_eta",&mono_eta);
        tree->SetBranchAddress("mono_phi",&mono_phi);
        tree->SetBranchAddress("amon_eta",&amon_eta);
        tree->SetBranchAddress("amon_phi",&amon_phi);
        tree->SetBranchAddress("pho_N",&nPhoton);
	tree->SetBranchAddress("pho_eta",&pho_eta);
	tree->SetBranchAddress("pho_phi",&pho_phi);
	tree->SetBranchAddress("pho_pt",&pho_pt);


	double subHits_count=0;
	double subSatHits_count=0;
        

	const unsigned NEvents = tree->GetEntries();
	vector<double> Et;
	
	
	MonoCuts noTrgAnalysis("NO TRG",oFile);
	MonoCuts TrgAnalysis("HLT_Photon200",oFile);

	vector<MonoCandidate> cand(10);	
	vector<Photon> photon(0);
	
        for(unsigned ev=0; ev<NEvents;ev++){
             	cout<<"ev "<<ev<<endl;
		tree->GetEntry(ev);
		
		if(nCandidates>Et.size()) Et.resize(nCandidates);
		if(nCandidates>cand.size()) cand.resize(nCandidates);

		if(nPhoton>photon.size()) photon.resize(nPhoton);

		for(unsigned i=0;i<nCandidates;i++){
        	
			Et[i]= (*e55)[i]/(TMath::CosH(TMath::Abs((*eta)[i])));
/*			cout<<"Et pass"<<endl;
			cout<<"subHits "<<(*subHits)[i]<<endl;
			cout<<"subSatHits "<<(*subSatHits)[i]<<endl;
			cout<<"tIso "<<(*tIso)[i]<<endl;
			cout<<"xyp0 "<<(*xyp0)[i]<<endl;
			cout<<"xyp1 "<<(*xyp1)[i]<<endl;
			cout<<"xyp2 "<<(*xyp2)[i]<<endl;
			cout<<"rzp0 "<<(*rzp0)[i]<<endl;
			cout<<"rzp1 "<<(*rzp1)[i]<<endl;
			cout<<"rzp2 "<<(*rzp2)[i]<<endl;
			cout<<"dist "<<(*dist)[i]<<endl;
			cout<<"f51 " <<(*f51)[i]<<endl;
			cout<<"Et "  <<Et[i]<<endl;
			cout<<"hIso "<<(*hIso)[i]<<endl;
			cout<<"eta " <<(*eta)[i]<<endl;
			cout<<"phi " <<(*phi)[i]<<endl;
*/
		cand[i] = MonoCandidate(
	        (*subHits)[i],
	        (*subSatHits)[i],
	        (*dEdXSig)[i],
	        (*tIso)[i],
	        (*xyp0)[i],
	        (*xyp1)[i],
	        (*xyp2)[i],
	        (*rzp0)[i],
	        (*rzp1)[i],
	        (*rzp2)[i],
	        (*dist)[i],
	        (*f51)[i],
	        (*f15)[i],
		Et[i],
//	        (*e55)[i],
	        (*hIso)[i],
	        (*eta)[i],
	        (*phi)[i],
                mono_eta,
                mono_phi,
                amon_eta,
                amon_phi,
		event,
		NPV	
     		 );
		}
		if(nPhoton!=0){
		for(unsigned j=0;j<nPhoton;j++){
			(*pho_eta)[j],
			(*pho_phi)[j],
			(*pho_pt)[j]
			);
			}
		}
			noTrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,true,ev);			
//                        TrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon200,ev);
	}//for every event loop
//	TrgAnalysis.WritePlots(oFile);
//	TrgAnalysis.SignalEff("HLT_Photon200");

	noTrgAnalysis.WritePlots(oFile);
//	noTrgAnalysis.SignalEff("NO TRG");
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
