//////////////////////////////////
//
//	File name: MonoAnalyzerMETTrigger.cc
//	author: Shih Lin
//	Content: Matching
//		 Photon information(21/4/24 update) 
//		 MET / MET trigger path(21/5/11 update)
//		
//	purpose: Decide MET selection 
//		 Trigger efficiency 
//	Output: Monopole_2018_METTrigger_1000.root
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
  MET,
  E55,
  F51,                 // frac 51
  HcalIso,             // Hcal Iso
  ABCD,
};
static const unsigned nPlot = 40U;
class PlotSet
{
public:
  PlotSet(){plots_.resize(nPlot);}
  ~PlotSet(){}
  void CreatPlot(const PlotName pn, TH1* h){ 
	plots_[pn] = h;
  }
  TH1 * GetPlot(const PlotName pn){ return plots_[pn]; }
  void WritePlot(){
	for(int pn=0;pn<nPlot;pn++){
	TH1 *h = plots_[pn];
	if(h){ 
	  h->Write();
	  }
	}
  }
private:
  vector<TH1*> plots_;
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
class MonoCandidate
{
public:
  MonoCandidate(){}
  //This will be used in main function and to absord the data in root file
    MonoCandidate(double sh, double satsh, double dedxsig,double tiso, double xyp0, double xyp1, double xyp2,
    double rzp0, double rzp1, double rzp2,
    double dist, double f51, double f15,
    double e55, double hiso, double eta,
    double phi, double mono_eta, double mono_phi, double mono_px, double mono_py, double mono_p,
    double amon_eta, double amon_phi, double amon_px, double amon_py, double amon_p,
    double event,double NPV, double MET_pt, double MET_phi, double GenMET_px,double GenMET_py,double GenMET_pt,double GenMET_phi,
    double CaloMET_px,double CaloMET_py,double CaloMET_pt,double CaloMET_phi):
  subHits_(sh),subSatHits_(satsh),dEdXSig_(dedxsig),tIso_(tiso),xyp0_(xyp0),
  xyp1_(xyp1),xyp2_(xyp2),rzp0_(rzp0),rzp1_(rzp1),rzp2_(rzp2),
  dist_(dist),f51_(f51),f15_(f15),e55_(e55),hIso_(hiso),
  eta_(eta),phi_(phi),mono_eta_(mono_eta), mono_phi_(mono_phi),mono_px_(mono_px),mono_py_(mono_py),mono_p_(mono_p),
  amon_eta_(amon_eta), amon_phi_(amon_phi),amon_px_(amon_px),amon_py_(amon_py),amon_p_(amon_p),
  event_(event),NPV_(NPV),MET_pt_(MET_pt),MET_phi_(MET_phi), GenMET_px_(GenMET_px),GenMET_py_(GenMET_py),GenMET_pt_(GenMET_pt),GenMET_phi_(GenMET_phi),
  CaloMET_px_(CaloMET_px),CaloMET_py_(CaloMET_py),CaloMET_pt_(CaloMET_pt),CaloMET_phi_(CaloMET_phi) { }
  //This will be used in comparing with cut
  MonoCandidate(const MonoCandidate &mc) : 
    subHits_(mc.subHits_),subSatHits_(mc.subSatHits_),dEdXSig_(mc.dEdXSig_),tIso_(mc.tIso_),
    xyp0_(mc.xyp0_),xyp1_(mc.xyp1_),xyp2_(mc.xyp2_),
    rzp0_(mc.rzp0_),rzp1_(mc.rzp1_),rzp2_(mc.rzp2_),
    dist_(mc.dist_),f51_(mc.f51_),f15_(mc.f15_),e55_(mc.e55_),
    hIso_(mc.hIso_),eta_(mc.eta_),phi_(mc.phi_),
  mono_eta_(mc.mono_eta_), mono_phi_(mc.mono_phi_),mono_px_(mc.mono_px_),mono_py_(mc.mono_py_),mono_p_(mc.mono_p_),
  amon_eta_(mc.amon_eta_), amon_phi_(mc.amon_phi_),amon_px_(mc.amon_px_),amon_py_(mc.amon_py_),amon_p_(mc.amon_p_),
  event_(mc.event_),NPV_(mc.NPV_),MET_pt_(mc.MET_pt_),MET_phi_(mc.MET_phi_),
  GenMET_px_(mc.GenMET_px_),GenMET_py_(mc.GenMET_py_),GenMET_pt_(mc.GenMET_pt_),GenMET_phi_(mc.GenMET_phi_),
  CaloMET_px_(mc.CaloMET_px_),CaloMET_py_(mc.CaloMET_py_),CaloMET_pt_(mc.CaloMET_pt_),CaloMET_phi_(mc.CaloMET_phi_) { } 
        
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
  double mono_px_;
  double mono_py_;
  double mono_p_;
  double amon_eta_;
  double amon_phi_;
  double amon_px_;
  double amon_py_;
  double amon_p_;
  double MET_pt_;
  double MET_phi_;
  double GenMET_pt_;
  double GenMET_px_;
  double GenMET_py_;
  double GenMET_phi_;
  double CaloMET_pt_;
  double CaloMET_px_;
  double CaloMET_py_;
  double CaloMET_phi_;

};
class MonoCuts:public MonoCandidate, public Photon
{
public:
  MonoCuts(){}

  MonoCuts(const string &trName,TFile *openFile):trigName_(trName){
	
	//create a new directory in the output file
	openFile->mkdir(trName.c_str());
	//cout<<"create directory "<<trName.c_str()<<endl;	
	//create histrogram N-1 cutf for No trigger
	Plot.resize(1U);
	PlotSet &x = Plot[0];
        x.CreatPlot(FracSatVNstrips,new TH2D("FracSatVNstrips","",100,0,1000,100,0,1));
        x.CreatPlot(DedXSig,new TH1D("DedXSig","",30,0,30));
        x.CreatPlot(RZcurv,new TH1D("RZcurv","",30,-0.01,0.01));
	x.CreatPlot(MET,new TH1D("MET","",30,-1,1200));
	x.CreatPlot(E55,new TH1D("E55","",30,-1,1200));
        x.CreatPlot(F51,new TH1D("F51","",30,0.2,1.1));
        x.CreatPlot(HcalIso,new TH1D("HcalIso","",30,-1,10));
        x.CreatPlot(ABCD,new TH2D("ABCD","",100,0,1.1,100,0,30));
	//x.CreatPlot(Eff,new TProfile("Eff","",100,0,800,6,0,1.2));
	//Name the cut name
	 cutName_[0] = "Quality_";
//	 cutName_[1] = "Energy_";
	 cutName_[1] = "MET_";
	 cutName_[2] = "F51_";
	 cutName_[3] = "dEdxSig_";
	 cutName_[4] = "HLT_";
	
	//create plot with n-1 cut

        n_1Plot.resize(nCut);
        string trgn1name = "HLT_";
        string name = "N_1_"+trgn1name;         //  N_1_HLT_FracSatVNstrips

	for( int c = 0; c <nCut ;c++){

	string cutn1name = name + cutName_[c];	//N_1_HLT_Quality_FracSatVNstrips
	PlotSet &z = n_1Plot[c];//[1] = Quality,[2]energy, [3] f51, [4]dEdx
	//create the n-1 cut plots
	
	z.CreatPlot(FracSatVNstrips,new TH2D((cutn1name+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
	z.CreatPlot(DedXSig,new TH1D((cutn1name+"DedXSig").c_str(),"",30,0,30));
	z.CreatPlot(RZcurv,new TH1D((cutn1name+"RZcurv").c_str(),"",30,-0.01,0.01));
        z.CreatPlot(MET,new TH1D((cutn1name+"MET").c_str(),"",30,0,800));
        z.CreatPlot(E55,new TH1D((cutn1name+"E55").c_str(),"",30,-1,1200));
        z.CreatPlot(F51,new TH1D((cutn1name+"F51").c_str(),"",30,0.2,1.1));
        z.CreatPlot(HcalIso,new TH1D((cutn1name+"HcalIso").c_str(),"",30,-1,10));
        z.CreatPlot(ABCD,new TH2D((cutn1name+"ABCD").c_str(),"",100,0,1.1,100,0,30));
	}

	//create plot with cutflow
	CutFlow.resize(nCut); 
	for( int c = 0;c<nCut-1;c++){
	
	PlotSet &y = CutFlow[c];
	string cutflowName = "Flow_HLT_"+cutName_[c];

        y.CreatPlot(FracSatVNstrips,new TH2D((cutflowName+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
        y.CreatPlot(DedXSig,new TH1D((cutflowName+"DedXSig").c_str(),"",30,0,30));
        y.CreatPlot(RZcurv,new TH1D((cutflowName+"RZcurv").c_str(),"",30,-0.01,0.01));
        y.CreatPlot(MET,new TH1D((cutflowName+"MET").c_str(),"",30,0,800));
        y.CreatPlot(E55,new TH1D((cutflowName+"E55").c_str(),"",30,-1,1200));
        y.CreatPlot(F51,new TH1D((cutflowName+"F51").c_str(),"",30,0.2,1.1));
        y.CreatPlot(HcalIso,new TH1D((cutflowName+"HcalIso").c_str(),"",30,-1,10));
        y.CreatPlot(ABCD,new TH2D((cutflowName+"ABCD").c_str(),"",100,0,1.1,100,0,30));
	}
  	
	//create trigger efficiency
	
  }

  ~MonoCuts(){}
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev)
  {
	
//	cout<<"enter doAnalysis successfully"<<endl;	
	EcalShower.clear();
	EcalSpike.clear();
	Trk.clear();
	Quality.clear();	

        HighPtPhoton.clear();
	CutFlowCand_TRG.clear();	
	CutFlowCand_Qual.clear();
	CutFlowCand_MET.clear();
	CutFlowCand_F51.clear();
	CutFlowCand_Dedx.clear();
	N1CutCand_TRG.clear();	
	N1CutCand_Qual.clear();
	N1CutCand_MET.clear();
	N1CutCand_F51.clear();
	N1CutCand_Dedx.clear();
	WithoutMETCut.clear();
	N1CutWithMETCut.clear();
	
	EcalSpikeNum=0;
	EcalShowerNum=0;
	TrackerMono=0;

     for(unsigned c=0;c<nCandidates;c++){
	MonoCandidate &cands = cand[c];
	bool QualCut = evalQuality(cands);
	bool ECut = evalE(cands);
	bool METCut = evalMET(cands);
	bool F15Cut = evalF15(cands);
	bool F51Cut = evalF51(cands);
	bool dEdXCut = evaldEdX(cands);

//	cout<<"enter the cut loop"<<endl; 
	//no cut
	if(TRG) CutFlowCand_TRG.push_back(cands);
	//cut all 
//	if(TRG && QualCut ) CutFlowCand_Qual.push_back(cands); 
//	if(TRG && QualCut && ECut && METCut ) CutFlowCand_MET.push_back(cands);
//	if(TRG && QualCut && ECut && METCut && F51Cut) CutFlowCand_F51.push_back(cands);
//	if(TRG && QualCut && ECut && METCut && F51Cut && dEdXCut) CutFlowCand_Dedx.push_back(cands);

	if(TRG && QualCut ) CutFlowCand_Qual.push_back(cands); 
	if(TRG && QualCut && METCut ) CutFlowCand_MET.push_back(cands);
	if(TRG && QualCut && METCut && F51Cut) CutFlowCand_F51.push_back(cands);
	if(TRG && QualCut && METCut && F51Cut && dEdXCut) CutFlowCand_Dedx.push_back(cands);
        if(TRG && QualCut && F51Cut && dEdXCut) WithoutMETCut.push_back(cands);
	//N-1 cut 
//	if( ECut && METCut && F51Cut && dEdXCut &&TRG) N1CutCand_Qual.push_back(cands); 
//        if( QualCut && F51Cut && dEdXCut &TRG) N1CutCand_MET.push_back(cands);
//        if( QualCut && ECut && METCut && dEdXCut &&TRG ) N1CutCand_F51.push_back(cands);
//	if( QualCut && ECut && METCut && F51Cut &&TRG)  N1CutCand_Dedx.push_back(cands);
//	if( QualCut && ECut && METCut && F51Cut && dEdXCut) N1CutCand_TRG.push_back(cands);
		
	if(  METCut && F51Cut && dEdXCut && TRG) N1CutCand_Qual.push_back(cands); 
        if( QualCut && F51Cut && dEdXCut && TRG) N1CutCand_MET.push_back(cands);
        if( QualCut && METCut && dEdXCut && TRG) N1CutCand_F51.push_back(cands);
	if( QualCut && METCut && F51Cut  && TRG)  N1CutCand_Dedx.push_back(cands);
	if( QualCut && METCut && F51Cut  && dEdXCut) N1CutCand_TRG.push_back(cands);
	if( QualCut && F51Cut  && dEdXCut) N1CutWithMETCut.push_back(cands);
	
/*	//MET
	if(QualCut) 			 			Quality.push_back(cands); 
//	if(dEdXCut && ECut)   					EcalShower.push_back(cands);//need to add f51<0.85 in the following code)
//	if(dEdXCut && ECut && F51Cut)				EcalSpike.push_back(cands); 
//	if(dEdXCut)	         				Trk.push_back(cands);
	if(QualCut && dEdXCut)	         			Trk.push_back(cands);
	if(QualCut && dEdXCut && ECut && F51Cut)		EcalSpike.push_back(cands); 
	if(QualCut && dEdXCut && ECut ) 			EcalShower.push_back(cands); 
*/
     }//end of cand scan	

    if(nPhoton!=0){
    for(unsigned p=0;p<nPhoton;p++){
        Photon &photon = pho[p];
        bool PhotonPtCut = evalPhoton(photon);
        //HighPtPhoton.push_back(photon);
        if(PhotonPtCut) HighPtPhoton.push_back(photon);
    }
    }

	//cout<<"start to fill all plots"<<endl;
	sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	
	if(CutFlowCand_TRG.size()>0) 
	{
		//cout<<"start to Get CandFlowCand_TRG plot"<<endl;
	        PlotSet &x = Plot[0];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_TRG[0].subHits_,CutFlowCand_TRG[0].subSatHits_/CutFlowCand_TRG[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_TRG[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_TRG[0].rzp2_);
		x.GetPlot(MET)->Fill(CutFlowCand_TRG[0].CaloMET_pt_);
		x.GetPlot(E55)->Fill(CutFlowCand_TRG[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_TRG[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_TRG[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_TRG[0].f51_,CutFlowCand_TRG[0].dEdXSig_);
		count++;
	
	}
	sort(CutFlowCand_Qual.begin(),CutFlowCand_Qual.begin()+CutFlowCand_Qual.size());

	if(CutFlowCand_Qual.size()>0) 
	{
		//cout<<"start to Get CandFlowCand_Qual plot"<<endl;
		PlotSet &x = CutFlow[0];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_Qual[0].subHits_,CutFlowCand_Qual[0].subSatHits_/CutFlowCand_Qual[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_Qual[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_Qual[0].rzp2_);
		x.GetPlot(MET)->Fill(CutFlowCand_Qual[0].CaloMET_pt_);
        	x.GetPlot(E55)->Fill(CutFlowCand_Qual[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_Qual[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_Qual[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_Qual[0].f51_,CutFlowCand_Qual[0].dEdXSig_);
		Qual_count++;	

	}
        sort(CutFlowCand_MET.begin(),CutFlowCand_MET.begin()+CutFlowCand_MET.size());
	if(CutFlowCand_MET.size()>0)
	{
	        PlotSet &x = CutFlow[1];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_MET[0].subHits_,CutFlowCand_MET[0].subSatHits_/CutFlowCand_MET[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_MET[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_MET[0].rzp2_);
		x.GetPlot(MET)->Fill(CutFlowCand_MET[0].CaloMET_pt_);
        	x.GetPlot(E55)->Fill(CutFlowCand_MET[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_MET[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_MET[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_MET[0].f51_,CutFlowCand_MET[0].dEdXSig_);
		E_count++;	
	}
        sort(CutFlowCand_F51.begin(),CutFlowCand_F51.begin()+CutFlowCand_F51.size());
	if(CutFlowCand_F51.size()>0)
	{
	        PlotSet &x = CutFlow[2];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_F51[0].subHits_,CutFlowCand_F51[0].subSatHits_/CutFlowCand_F51[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_F51[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_F51[0].rzp2_);
		x.GetPlot(MET)->Fill(CutFlowCand_F51[0].CaloMET_pt_);
        	x.GetPlot(E55)->Fill(CutFlowCand_F51[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_F51[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_F51[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_F51[0].f51_,CutFlowCand_F51[0].dEdXSig_);
		f51_count++;	
	}
	
	sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());
	if(CutFlowCand_Dedx.size()>0)
	{
		        
		PlotSet &x = CutFlow[3];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_Dedx[0].subHits_,CutFlowCand_Dedx[0].subSatHits_/CutFlowCand_Dedx[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_Dedx[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_Dedx[0].rzp2_);
		x.GetPlot(MET)->Fill(CutFlowCand_Dedx[0].CaloMET_pt_);
        	x.GetPlot(E55)->Fill(CutFlowCand_Dedx[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_Dedx[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_Dedx[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_Dedx[0].f51_,CutFlowCand_Dedx[0].dEdXSig_);
		dEdX_count++;	
	}

	sort(WithoutMETCut.begin(),WithoutMETCut.begin()+WithoutMETCut.size());
	if(WithoutMETCut.size()>0)  CutAllexceptMET++;

	sort(N1CutWithMETCut.begin(),N1CutWithMETCut.begin()+N1CutWithMETCut.size());
	if(N1CutWithMETCut.size()>0)  CutAllexceptMETandTRG++;


	// N-1 Cut
        sort(N1CutCand_Qual.begin(),N1CutCand_Qual.begin()+N1CutCand_Qual.size());
	if(N1CutCand_Qual.size()>0) 
	{
	        PlotSet &z = n_1Plot[0];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_Qual[0].subHits_,N1CutCand_Qual[0].subSatHits_/N1CutCand_Qual[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_Qual[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_Qual[0].rzp2_);
       		z.GetPlot(MET)->Fill(N1CutCand_Qual[0].CaloMET_pt_);
       		z.GetPlot(E55)->Fill(N1CutCand_Qual[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_Qual[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_Qual[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_Qual[0].f51_,N1CutCand_Qual[0].dEdXSig_);
		NoQual++;

	}
        

	sort(N1CutCand_MET.begin(),N1CutCand_MET.begin()+N1CutCand_MET.size());
	if(N1CutCand_MET.size()>0)
	{
	        PlotSet &z = n_1Plot[1];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_MET[0].subHits_,N1CutCand_MET[0].subSatHits_/N1CutCand_MET[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_MET[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_MET[0].rzp2_);
       		z.GetPlot(MET)->Fill(N1CutCand_MET[0].CaloMET_pt_);
        	z.GetPlot(E55)->Fill(N1CutCand_MET[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_MET[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_MET[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_MET[0].f51_,N1CutCand_MET[0].dEdXSig_);
		NoE++;	
	}


        sort(N1CutCand_F51.begin(),N1CutCand_F51.begin()+N1CutCand_F51.size());
	if(N1CutCand_F51.size()>0)
	{
	        PlotSet &z = n_1Plot[2];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_F51[0].subHits_,N1CutCand_F51[0].subSatHits_/N1CutCand_F51[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_F51[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_F51[0].rzp2_);
       		z.GetPlot(MET)->Fill(N1CutCand_F51[0].CaloMET_pt_);
        	z.GetPlot(E55)->Fill(N1CutCand_F51[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_F51[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_F51[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_F51[0].f51_,N1CutCand_F51[0].dEdXSig_);
		NoF51++;	
	}
	

	sort(N1CutCand_Dedx.begin(),N1CutCand_Dedx.begin()+N1CutCand_Dedx.size());
	if(N1CutCand_Dedx.size()>0)
	{
	
	        PlotSet &z = n_1Plot[3];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_Dedx[0].subHits_,N1CutCand_Dedx[0].subSatHits_/N1CutCand_Dedx[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_Dedx[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_Dedx[0].rzp2_);
       		z.GetPlot(MET)->Fill(N1CutCand_Dedx[0].CaloMET_pt_);
        	z.GetPlot(E55)->Fill(N1CutCand_Dedx[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_Dedx[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_Dedx[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_Dedx[0].f51_,N1CutCand_Dedx[0].dEdXSig_);
		NodEdXCut++;	
	}

        sort(N1CutCand_TRG.begin(),N1CutCand_TRG.begin()+N1CutCand_TRG.size());
	if(N1CutCand_TRG.size()>0) 
	{
	        PlotSet &z = n_1Plot[4];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_TRG[0].subHits_,N1CutCand_TRG[0].subSatHits_/N1CutCand_TRG[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_TRG[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_TRG[0].rzp2_);
       		z.GetPlot(MET)->Fill(N1CutCand_TRG[0].CaloMET_pt_);
        	z.GetPlot(E55)->Fill(N1CutCand_TRG[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_TRG[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_TRG[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_TRG[0].f51_,N1CutCand_TRG[0].dEdXSig_);
		NoTRG++;
	}
		bool RealEvt=0;	
  
/*	cout<<"CutFlowCand_TRG size "<<CutFlowCand_TRG.size()<<endl;
	for(int i=0;i<CutFlowCand_TRG.size();i++){
//		cout<<"Trk size "<<Trk.size()<<endl;
                double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(CutFlowCand_TRG[i].eta_-CutFlowCand_TRG[0].mono_eta_,2)+
                                pow(CutFlowCand_TRG[i].phi_-CutFlowCand_TRG[0].mono_phi_,2));
                am_deltaR= sqrt(pow(CutFlowCand_TRG[i].eta_-CutFlowCand_TRG[0].amon_eta_,2)+
                                pow(CutFlowCand_TRG[i].phi_-CutFlowCand_TRG[0].amon_phi_,2));
//		cout<<"starting matching"<<endl;
                if(m_deltaR<0.15||am_deltaR<0.15)
                {
			cout<<i<<endl;
			cout<<"     eta "<<setprecision(5)<<CutFlowCand_TRG[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<CutFlowCand_TRG[i].phi_<<endl;
			cout<<"      Et "<<setprecision(5)<<CutFlowCand_TRG[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<CutFlowCand_TRG[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<CutFlowCand_TRG[i].f51_<<endl;
			cout<<"     f15 "<<setprecision(5)<<CutFlowCand_TRG[i].f15_<<endl;
			cout<<"     MET "<<setprecision(5)<<CutFlowCand_TRG[i].CaloMET_pt_<<endl;
			RealNum++;
			RealEvt=1;
		
			double mono_pt =0;
			mono_pt = sqrt(pow(CutFlowCand_TRG[0].mono_px_,2)+pow(CutFlowCand_TRG[0].mono_py_,2));
			double amon_pt =0;
			amon_pt = sqrt(pow(CutFlowCand_TRG[0].amon_px_,2)+pow(CutFlowCand_TRG[0].amon_py_,2));
              		cout<<"          candidate           monoGen         antimonoGen"<<endl;
                	cout<<"eta      "<<setprecision(5)<<CutFlowCand_TRG[i].eta_<<setw(20)<<CutFlowCand_TRG[i].mono_eta_
                        	         <<setw(20)<<CutFlowCand_TRG[i].amon_eta_<<endl;
                	cout<<"phi      "<<setprecision(5)<<CutFlowCand_TRG[i].phi_<<setw(20)<<CutFlowCand_TRG[i].mono_phi_
                       		         <<setw(20)<<CutFlowCand_TRG[i].amon_phi_<<endl;
			cout<<"pT	"<<setw(27)<<mono_pt<<setw(20)<<amon_pt<<endl;
                	cout<<"m deltaR "<<m_deltaR<<endl;
                	cout<<"a deltaR "<<am_deltaR<<endl;
			cout<<"--------------------------------------------"<<endl;
		}
	}*/

/*		cout<<"Print out Photon information----------"<<endl;        
                for(int j=0;j<HighPtPhoton.size();j++){
                        cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                        cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                        cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
                }
*/
                //to see monopole is photon-like or spike-like
                //define: photon-monopole angle<0.5 is photon-like
            /*    for(int j=0;j<HighPtPhoton.size();j++){
                        double photonMono_deltaR = 0;
                        photonMono_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-HighPtPhoton[j].pho_eta_,2)+
                                                 pow(CutFlowCand_Dedx[i].phi_-HighPtPhoton[j].pho_phi_,2));
                        if(photonMono_deltaR<0.15){
                                cout<<"deltaR "<<photonMono_deltaR<<endl;
                                cout<<"photon-like monopole"<<endl;
                                photonLike++;// otherwise not photon-like monopole
                        }
                }(*/


        	//cout<<"event tag"<<CutFlowCand_TRG[0].event_<<endl;
//		if(RealEvt==1) RealEvNum++;//If there is any monopole in event(where in tracker or ecal or number), the event is taken
//		RealEvt=0;
//			cout<<"=================================================="<<endl;		
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
	for(unsigned c=0; c<nCut; c++) n_1Plot[c].WritePlot();
  //       cout<<"n-1cut pass writeplots func"<<endl;

	for(unsigned c=0; c<nCut; c++) CutFlow[c].WritePlot();
//	cout<<"cutflow pass writeplots func"<<endl;
	
  }
  
  void SignalEff( const string trName,unsigned TotalEvents)
  {
        //signal efficiency = no. of events after all selection cuts/all events
	cout<<"3000GeV";	
	cout<<endl;
	cout<<trName<<" ================================="<<endl;
        cout<<"        TRG "<<count<<endl;
        cout<<"QualityCuts "<<Qual_count<<endl;
//        cout<<"       ECut "<<E_count<<endl;
        cout<<"     METCut "<<E_count<<endl;
        cout<<"     F51Cut "<<f51_count<<endl;
        cout<<" dEdXSigCut "<<dEdX_count<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"Relative effciency count------------------"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"     No TRG "<<NoTRG<<endl;//Quality+f51+dEdx
        cout<<" No Quality "<<NoQual<<endl;
        cout<<"  No METCut "<< NoE <<endl;
        //cout<<"    No ECut "<< NoE <<endl;
        cout<<"  No F51Cut "<<NoF51<<endl;
        cout<<" No dEdXSig "<<NodEdXCut<<endl;
        cout<<"=========================================="<<endl;
        cout<<"Signal efficiency = "<<(double)dEdX_count/(double)TotalEvents<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"Relative efficiency-----------------------"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"        TRG "<<(double)dEdX_count/(double)NoTRG<<endl;
        cout<<"QualityCuts "<<(double)dEdX_count/(double)NoQual<<endl;
        cout<<"     METCut "<<(double)dEdX_count/(double)NoE<<endl;
        cout<<"     F51Cut "<<(double)dEdX_count/(double)NoF51<<endl;
        cout<<" dEdXSigCut "<<(double)dEdX_count/(double)NodEdXCut<<endl;
	cout<<"=========================================="<<endl;

	double TrgEff1 =0 ,TrgEff2 = 0,TrgEff3 = 0;
	TrgEff1 = (double)count/(double)TotalEvents;
	TrgEff2 = (double)CutAllexceptMET/(double)CutAllexceptMETandTRG;
	TrgEff3 = (double)dEdX_count/(double)NoTRG;
	cout<<"TrgEff 1 "<<TrgEff1<<endl;	
	cout<<"TrgEff 2 "<<TrgEff2<<endl;	
	cout<<"TrgEff 3 "<<TrgEff3<<endl;	

	PlotSet &x = Plot[0];
	
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
  static const double METCut_ ;
  static const double f15Cut_ ;
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
  
  static const unsigned nCut = 6U;
  string cutName_[nCut];
  


  //cuts analysis
  bool evalQuality(MonoCandidate &mc) { return mc.xyp0_ < xyp0Cut_&& mc.xyp2_ > xyp2Cut_ 
			&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_;  }
  bool evalE(MonoCandidate &mc) { return mc.e55_ > e55Cut_; }
  bool evalMET(MonoCandidate &mc) { return mc.CaloMET_pt_ > METCut_; }
  bool evalF51(MonoCandidate &mc) { return mc.f51_ > f51Cut_ ; }
  bool evalF15(MonoCandidate &mc) { return mc.f15_ > f15Cut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
 
  vector<MonoCandidate> CutFlowCand_TRG;
  vector<MonoCandidate> CutFlowCand_Qual;
  vector<MonoCandidate> CutFlowCand_MET;
  vector<MonoCandidate> CutFlowCand_F51;
  vector<MonoCandidate> CutFlowCand_Dedx;
  vector<MonoCandidate> WithoutMETCut; 
  vector<MonoCandidate> N1CutWithMETCut;

  vector<MonoCandidate> N1CutCand_TRG;
  vector<MonoCandidate> N1CutCand_Qual;
  vector<MonoCandidate> N1CutCand_MET;
  vector<MonoCandidate> N1CutCand_F51;
  vector<MonoCandidate> N1CutCand_Dedx;
	
  vector<Photon> HighPtPhoton;

  vector<MonoCandidate> EcalShower;
  vector<MonoCandidate> EcalSpike;
  vector<MonoCandidate> Trk;
  vector<MonoCandidate> Quality;

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


  int TrackerMono=0;
  int EcalSpikeNum=0;
  int EcalShowerNum=0;
 
  int TwoSpikeLike=0;
  int TwoShowerLike=0;
  int OneShowerOneSpike=0;
  int OneShowerOneTrk=0;
  int OneSpikeOneTrk=0;
  int TwoTrk=0;
  int OneShowerLike=0;
  int OneSpikeLike=0;
  int OneTrk=0;
  int ZeroMono=0;		

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

  int CutAllexceptMET =0;
  int CutAllexceptMETandTRG=0;

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
//  const double MonoCuts::e55Cut_ = 20;
  const double MonoCuts::METCut_ = 300;
  const double MonoCuts::e55Cut_ = 0;
//  const double MonoCuts::f51Cut_ = 0.6;
  const double MonoCuts::f15Cut_ = 0.85;
  const double MonoCuts::f51Cut_ = 0.85;
  const double MonoCuts::photonCut_ = 0;
void MonoAnalyzerCaloMET_Trigger_3000()
{
	cout<<"hello"<<endl;
	TFile *oFile = new TFile("Monopole_CaloMETTrigger_2018_3000_.root","recreate");
	cout<<"new file line pass"<<endl;

	TFile *fin = new TFile("MonoNtuple2018_MC_3000_0722.root");
	cout<<"open file success"<<endl;
        TTree *tree = (TTree*)fin->Get("monopoles");
	cout<<"open tree success"<<endl;
        Bool_t passHLT_Photon200;
	Bool_t passHLT_Photon175;
	Bool_t passHLT_DoublePhoton70;
        Bool_t passHLT_PFMET300;
        Bool_t passHLT_PFMET170_HBHE_BeamHaloCleaned;
        Bool_t passHLT_PFMET140_PFMHT140_IDTight;
        Bool_t passHLT_PFMET250_HBHECleaned;
        Bool_t passHLT_PFMET300_HBHECleaned;
        Bool_t passHLT_PFMET200_HBHE_BeamHaloCleaned;
        Bool_t passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;
        Bool_t passHLT_CaloMET300_HBHECleaned;                

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
	double MET_phi;
	double MET_pt;
	double GenMET_pt;
	double GenMET_px;
	double GenMET_py;
	double GenMET_phi;
	double CaloMET_pt;
	double CaloMET_px;
	double CaloMET_py;
	double CaloMET_phi;



	double  mono_px;
        double  mono_py;
        double  mono_p;
        double  amon_px;
        double  amon_py;
        double  amon_p;

	tree->SetBranchAddress("event",&event);
	tree->SetBranchAddress("trigResult",&trigResults);
        tree->SetBranchAddress("trigNames",&trigNames);
        tree->SetBranchAddress("passHLT_Photon200" , &passHLT_Photon200);
	tree->SetBranchAddress("passHLT_Photon175" , &passHLT_Photon175);
        tree->SetBranchAddress("passHLT_DoublePhoton70",&passHLT_DoublePhoton70);
        tree->SetBranchAddress("passHLT_PFMET300",&passHLT_PFMET300);
        tree->SetBranchAddress("passHLT_PFMET170_HBHE_BeamHaloCleaned",&passHLT_PFMET170_HBHE_BeamHaloCleaned);
        tree->SetBranchAddress("passHLT_PFMET140_PFMHT140_IDTight",&passHLT_PFMET140_PFMHT140_IDTight);
        tree->SetBranchAddress("passHLT_PFMET250_HBHECleaned",&passHLT_PFMET250_HBHECleaned);
        tree->SetBranchAddress("passHLT_PFMET300_HBHECleaned",&passHLT_PFMET300_HBHECleaned);
        tree->SetBranchAddress("passHLT_PFMET200_HBHE_BeamHaloCleaned",&passHLT_PFMET200_HBHE_BeamHaloCleaned);
        tree->SetBranchAddress("passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",&passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned);
        tree->SetBranchAddress("passHLT_CaloMET300_HBHECleaned",&passHLT_CaloMET300_HBHECleaned);
	
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
        tree->SetBranchAddress("mono_px",&mono_px);
        tree->SetBranchAddress("mono_py",&mono_py);
        tree->SetBranchAddress("mono_p",&mono_p);
        tree->SetBranchAddress("amon_eta",&amon_eta);
        tree->SetBranchAddress("amon_phi",&amon_phi);
        tree->SetBranchAddress("amon_px",&amon_px);
        tree->SetBranchAddress("amon_py",&amon_py);
        tree->SetBranchAddress("amon_p",&amon_p);
        tree->SetBranchAddress("pho_N",&nPhoton);
	tree->SetBranchAddress("pho_eta",&pho_eta);
	tree->SetBranchAddress("pho_phi",&pho_phi);
	tree->SetBranchAddress("pho_pt",&pho_pt);
	tree->SetBranchAddress("mpt_phi",&MET_phi);
	tree->SetBranchAddress("mpt_pt",&MET_pt);
	tree->SetBranchAddress("GenMET_phi",&GenMET_phi);
	tree->SetBranchAddress("GenMET_px",&GenMET_px);
	tree->SetBranchAddress("GenMET_py",&GenMET_py);
	tree->SetBranchAddress("GenMET_pt",&GenMET_pt);
	tree->SetBranchAddress("CaloMET_px",&CaloMET_px);
	tree->SetBranchAddress("CaloMET_py",&CaloMET_py);
	tree->SetBranchAddress("CaloMET_pt",&CaloMET_pt);
	tree->SetBranchAddress("CaloMET_phi",&CaloMET_phi);
	

	double subHits_count=0;
	double subSatHits_count=0;
        

	const unsigned NEvents = tree->GetEntries();
	vector<double> Et;
	
	
	MonoCuts noTrgAnalysis("NOTRG",oFile);
//	MonoCuts TrgAnalysis_Photon                                 ("HLT_Photon200",oFile);
//	MonoCuts TrgAnalysis_PFMET300		                    ("HLT_PFMET300",oFile); //2016
//	MonoCuts TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned          ("HLT_PFMET170_HBHE_BeamHaloCleaned",oFile); //2016
//	MonoCuts TrgAnalysis_PFMET140_PFMHT140_IDTight              ("HLT_PFMET140_PFMHT140_IDTight",oFile);
//	MonoCuts TrgAnalysis_PFMET250_HBHECleaned                   ("HLT_PFMET250_HBHECleaned",oFile);
//	MonoCuts TrgAnalysis_PFMET300_HBHECleaned                   ("HLT_PFMET300_HBHECleaned",oFile);
//	MonoCuts TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned          ("HLT_PFMET200_HBHE_BeamHaloCleaned",oFile);
//	MonoCuts TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned   ("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",oFile);
	MonoCuts TrgAnalysis_CaloMET300_HBHECleaned                 ("HLT_CaloMET300_HBHECleaned",oFile);

	vector<MonoCandidate> cand(10);	
	vector<Photon> photon(0);
	
        for(unsigned ev=0; ev<NEvents;ev++){
             tree->GetEntry(ev);
		cout<<"ev "<<ev<<endl;	
		if(nCandidates>Et.size()) Et.resize(nCandidates);
		if(nCandidates>cand.size()) cand.resize(nCandidates);

		if(nPhoton>photon.size()) photon.resize(nPhoton);
	//	cout<<"Start to feed variables"<<endl;
		for(unsigned i=0;i<nCandidates;i++){
        	
			Et[i]= (*e55)[i]/(TMath::CosH(TMath::Abs((*eta)[i])));
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
//		Et[i],
	        (*e55)[i],
	        (*hIso)[i],
	        (*eta)[i],
	        (*phi)[i],
                mono_eta,
                mono_phi,
		mono_px,
                mono_py,
                mono_p,
                amon_eta,
                amon_phi,
                amon_px,
                amon_py,
		amon_p,
		event,
		NPV,
		//MET_px,
		//MET_py,
		MET_pt,
		MET_phi,
		GenMET_px,
		GenMET_py,
		GenMET_pt,
		GenMET_phi,
		CaloMET_px,
		CaloMET_py,
		CaloMET_pt,
		CaloMET_phi
			
     		 );
	//	cout<<"MET_pt "<<MET_pt<<endl;
		}
//		cout<<"enter photon loop"<<endl;
		if(nPhoton!=0){
		for(unsigned j=0;j<nPhoton;j++){
			photon[j] = Photon(
			(*pho_eta)[j],
			(*pho_phi)[j],
			(*pho_pt)[j]
			);
			}
		}
//		cout<<"start to do analysis"<<endl;
			noTrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,true,ev);
//		cout<<"do analysis :PFMET300"<<endl;	
  //       TrgAnalysis_PFMET300.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300,ev);
    //     TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET170_HBHE_BeamHaloCleaned,ev);
//         TrgAnalysis_PFMET140_PFMHT140_IDTight           .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET140_PFMHT140_IDTight,ev);
//         TrgAnalysis_PFMET250_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET250_HBHECleaned,ev);
        // TrgAnalysis_PFMET300_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300_HBHECleaned,ev);
  //       TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET200_HBHE_BeamHaloCleaned,ev);
    //     TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned,ev);
           TrgAnalysis_CaloMET300_HBHECleaned              .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_CaloMET300_HBHECleaned,ev);
//         TrgAnalysis_Photon		 		 .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon200,ev);
	}//for every event loop
	
// 			TrgAnalysis_Photon         			  .SignalEff("HLT_Photon200",NEvents); 
//		        TrgAnalysis_PFMET300		                  .SignalEff("HLT_PFMET300",NEvents);
//		        TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET170_HBHE_BeamHaloCleaned",NEvents);
//		        TrgAnalysis_PFMET140_PFMHT140_IDTight             .SignalEff("HLT_PFMET140_PFMHT140_IDTight",NEvents);
	//	        TrgAnalysis_PFMET250_HBHECleaned                  .SignalEff("HLT_PFMET250_HBHECleaned",NEvents);
      //		        TrgAnalysis_PFMET300_HBHECleaned                  .SignalEff("HLT_PFMET300_HBHECleaned",NEvents);
//		        TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET200_HBHE_BeamHaloCleaned",NEvents);
//		        TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned  .SignalEff("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",NEvents);
		        TrgAnalysis_CaloMET300_HBHECleaned                .SignalEff("HLT_CaloMET300_HBHECleaned",NEvents);                  


  //       TrgAnalysis_Photon                              .WritePlots(oFile);
//         TrgAnalysis_PFMET300		                 .WritePlots(oFile);       
//         TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .WritePlots(oFile);     
  //       TrgAnalysis_PFMET140_PFMHT140_IDTight           .WritePlots(oFile);
    //     TrgAnalysis_PFMET250_HBHECleaned                .WritePlots(oFile);   
         //TrgAnalysis_PFMET300_HBHECleaned                .WritePlots(oFile);
    //     TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .WritePlots(oFile);   
    //     TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.WritePlots(oFile);
         TrgAnalysis_CaloMET300_HBHECleaned              .WritePlots(oFile);   
	                                                     

	noTrgAnalysis.WritePlots(oFile);
	noTrgAnalysis.SignalEff("NOTRG",NEvents);
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
