//////////////////////////////////
//
//	File name: MonoAnalyzerPhoton.cc
//	author: Shih Lin
//	Content: **Mainly for analyzing HLT_Photon200 trigger**
//		   (the original analysis)
//		 the significance plots for every cuts
//		 count efficiency for every cuts
//		 Matching(21/2/25 update)
//		 Update file name (June 5 2021) 
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
  ABCD
};
static const unsigned nPlot = 8U;
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
        x.CreatPlot(FracSatVNstrips,new TH2D("FracSatVNstrips","",100,0,1000,100,0,1));
        x.CreatPlot(DedXSig,new TH1D("DedXSig","",100,0,30));
        x.CreatPlot(RZcurv,new TH1D("RZcurv","",100,-0.01,0.01));
//	x.CreatPlot(Energy,new TH1D("Energy","",100,0,1400));
	x.CreatPlot(E55,new TH1D("E55","",100,-1,1200));
        x.CreatPlot(F51,new TH1D("F51","",100,0.2,1.1));
        x.CreatPlot(HcalIso,new TH1D("HcalIso","",100,-1,10));
        x.CreatPlot(ABCD,new TH2D("ABCD","",100,0,1.1,100,0,30));
	
	cout<<"the address of plots x "<<&x<<endl;

	//Name the cut name
	 cutName_[0] = "Quality_";
	 cutName_[1] = "Energy_";
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
	z.CreatPlot(DedXSig,new TH1D((cutn1name+"DedXSig").c_str(),"",100,0,30));
	z.CreatPlot(RZcurv,new TH1D((cutn1name+"RZcurv").c_str(),"",100,-0.01,0.01));
  //      z.CreatPlot(Energy,new TH1D((cutn1name+"Energy").c_str(),"",100,0,1400));
        z.CreatPlot(E55,new TH1D((cutn1name+"E55").c_str(),"",100,-1,1200));
        z.CreatPlot(F51,new TH1D((cutn1name+"F51").c_str(),"",100,0.2,1.1));
        z.CreatPlot(HcalIso,new TH1D((cutn1name+"HcalIso").c_str(),"",100,-1,10));
        z.CreatPlot(ABCD,new TH2D((cutn1name+"ABCD").c_str(),"",100,0,1.1,100,0,30));
	cout<<"good "<<cutn1name<<endl;
	}

	//create plot with cutflow
	CutFlow.resize(nCut); 
	for( int c = 0;c<nCut-1;c++){
	
	PlotSet &y = CutFlow[c];
	string cutflowName = "Flow_HLT_"+cutName_[c];

        y.CreatPlot(FracSatVNstrips,new TH2D((cutflowName+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
        y.CreatPlot(DedXSig,new TH1D((cutflowName+"DedXSig").c_str(),"",100,0,30));
        y.CreatPlot(RZcurv,new TH1D((cutflowName+"RZcurv").c_str(),"",100,-0.01,0.01));
    //    y.CreatPlot(Energy,new TH1D((cutflowName+"Energy").c_str(),"",100,0,1400));
        y.CreatPlot(E55,new TH1D((cutflowName+"E55").c_str(),"",100,-1,1200));
        y.CreatPlot(F51,new TH1D((cutflowName+"F51").c_str(),"",100,0.2,1.1));
        y.CreatPlot(HcalIso,new TH1D((cutflowName+"HcalIso").c_str(),"",100,-1,10));
        y.CreatPlot(ABCD,new TH2D((cutflowName+"ABCD").c_str(),"",100,0,1.1,100,0,30));
        cout<<"good cutflow"<<cutflowName<<endl;
	}

  }

  ~MonoCuts(){}
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

//	cout<<"end of highphoton loop"<<endl;
	//cut flow events calculating
        sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	
	if(CutFlowCand_TRG.size()>0) 
	{
		if(CutFlowCand_TRG[0].dEdXSig_<999){	
	        PlotSet &x = Plot[0];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_TRG[0].subHits_,CutFlowCand_TRG[0].subSatHits_/CutFlowCand_TRG[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_TRG[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_TRG[0].rzp2_);
        	x.GetPlot(E55)->Fill(CutFlowCand_TRG[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_TRG[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_TRG[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_TRG[0].f51_,CutFlowCand_TRG[0].dEdXSig_);
		count++;
		}
		else{
		PlotSet &x = Plot[0];
                x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_TRG[1].subHits_,CutFlowCand_TRG[1].subSatHits_/CutFlowCand_TRG[1].subHits_);
                x.GetPlot(DedXSig)->Fill(CutFlowCand_TRG[1].dEdXSig_);
                x.GetPlot(RZcurv)->Fill(CutFlowCand_TRG[1].rzp2_);
                x.GetPlot(E55)->Fill(CutFlowCand_TRG[1].e55_);
                x.GetPlot(F51)->Fill(CutFlowCand_TRG[1].f51_);
                x.GetPlot(HcalIso)->Fill(CutFlowCand_TRG[1].hIso_);
                x.GetPlot(ABCD)->Fill(CutFlowCand_TRG[1].f51_,CutFlowCand_TRG[1].dEdXSig_);
                count++;
		}
		if(CutFlowCand_TRG[0].dEdXSig_>0&&CutFlowCand_TRG[0].dEdXSig_<5) TRGminitest++;
	
	}

	//count cutflow TRG+Qual and get plots
        sort(CutFlowCand_Qual.begin(),CutFlowCand_Qual.begin()+CutFlowCand_Qual.size());
	for(int i=0;i<CutFlowCand_Qual.size();i++){
		//cout<<"ev "<<ev<<", cutflow Quality order:"<<CutFlowCand_Qual[i].dEdXSig_<<endl;
	}
	if(CutFlowCand_Qual.size()>0) 
	{
	        if(CutFlowCand_Qual[0].dEdXSig_<999){
		PlotSet &x = CutFlow[0];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_Qual[0].subHits_,CutFlowCand_Qual[0].subSatHits_/CutFlowCand_Qual[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_Qual[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_Qual[0].rzp2_);
        	x.GetPlot(E55)->Fill(CutFlowCand_Qual[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_Qual[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_Qual[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_Qual[0].f51_,CutFlowCand_Qual[0].dEdXSig_);
		Qual_count++;	
		}
		else{
		
		PlotSet &x = CutFlow[0];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_Qual[1].subHits_,CutFlowCand_Qual[1].subSatHits_/CutFlowCand_Qual[1].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_Qual[1].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_Qual[1].rzp2_);
        	x.GetPlot(E55)->Fill(CutFlowCand_Qual[1].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_Qual[1].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_Qual[1].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_Qual[1].f51_,CutFlowCand_Qual[1].dEdXSig_);
		Qual_count++;
		}
                if(CutFlowCand_Qual[0].dEdXSig_>0&&CutFlowCand_Qual[0].dEdXSig_<5) Qualminitest++;

	}
	
	//count cutflow TRG+Qual+Energy and get plots
        sort(CutFlowCand_Energy.begin(),CutFlowCand_Energy.begin()+CutFlowCand_Energy.size());
	if(CutFlowCand_Energy.size()>0)
	{
	        PlotSet &x = CutFlow[1];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_Energy[0].subHits_,CutFlowCand_Energy[0].subSatHits_/CutFlowCand_Energy[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_Energy[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_Energy[0].rzp2_);
        	x.GetPlot(E55)->Fill(CutFlowCand_Energy[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_Energy[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_Energy[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_Energy[0].f51_,CutFlowCand_Energy[0].dEdXSig_);
		E_count++;	
	}

	//TRG+Qual+Energy+F51
        sort(CutFlowCand_F51.begin(),CutFlowCand_F51.begin()+CutFlowCand_F51.size());
	if(CutFlowCand_F51.size()>0)
	{
	        PlotSet &x = CutFlow[2];
	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_F51[0].subHits_,CutFlowCand_F51[0].subSatHits_/CutFlowCand_F51[0].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_F51[0].dEdXSig_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_F51[0].rzp2_);
        	x.GetPlot(E55)->Fill(CutFlowCand_F51[0].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_F51[0].f51_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_F51[0].hIso_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_F51[0].f51_,CutFlowCand_F51[0].dEdXSig_);
		f51_count++;	
	}
	
	sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());
	//Count no. of monopole signal
/*	if(CutFlowCand_Dedx.size()==1){ 	
	/*	cout<<"ev "<<ev<<" has monopole"<<endl;
		cout<<"event tag"<<CutFlowCand_Dedx[0].event_<<endl;
		cout<<"eta "<<CutFlowCand_Dedx[0].eta_<<endl;
		cout<<"phi "<<CutFlowCand_Dedx[0].phi_<<endl;
		
		cout<<"---------------------------------------"<<endl;
		MonoNum1++;
	}
	else if(CutFlowCand_Dedx.size()==2) MonoNum2++;
	else if(CutFlowCand_Dedx.size()==3) MonoNum3++;
	if(CutFlowCand_Dedx.size()>=1&&CutFlowCand_Dedx.size()<=3){
		cout<<"ev "<<ev<<"  ,no. of cand "<<CutFlowCand_Dedx.size()<<endl;
		bool RealMono=0;
		bool RealAnti=0;

	   //go through all candidate in this event
	   for(int i=0; i<CutFlowCand_Dedx.size();i++){

		double m_deltaR=0;
		double am_deltaR=0;
		m_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].mono_eta_,2)+
				pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].mono_phi_,2));
		am_deltaR= sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].amon_eta_,2)+
                                pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].amon_phi_,2));

		if(m_deltaR<0.15) 
		{
			cout<<i+1<<" monopole match "<<endl;
			RealNum++;
			Reco++;
			RealMono=true;
		}
		else if(am_deltaR<0.15) 
		{
			cout<<i+1<<" anti monopole match"<<endl;
			RealAntiNum++;
			Reco++;
			RealAnti=true;
		}
		else 
		{
			cout<<ev<<" ";
			cout<<i+1<<" FAKE"<<endl;
			FakeNum++;
			cout<<"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"<<endl;
		}

		cout<<"          candidate           monoGen         antimonoGen"<<endl;
		cout<<"eta      "<<setprecision(5)<<CutFlowCand_Dedx[i].eta_<<setw(20)<<CutFlowCand_Dedx[i].mono_eta_
				 <<setw(20)<<CutFlowCand_Dedx[i].amon_eta_<<endl;
		cout<<"phi      "<<setprecision(5)<<CutFlowCand_Dedx[i].phi_<<setw(20)<<CutFlowCand_Dedx[i].mono_phi_
				 <<setw(20)<<CutFlowCand_Dedx[i].amon_phi_<<endl;
		cout<<"m deltaR "<<m_deltaR<<endl;
		cout<<"a deltaR "<<am_deltaR<<endl;
                cout<<"ev "<<ev<<" has monopole"<<endl;
                cout<<"event tag"<<CutFlowCand_Dedx[0].event_<<endl;

		cout<<"----------------------------"<<endl;
	
		
		cout<<"Print out all candidate in ev "<<ev<<endl;
		for(int i=0;i<CutFlowCand_TRG.size();i++){
			cout<<i<<endl;
			cout<<"     eta "<<setprecision(5)<<CutFlowCand_TRG[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<CutFlowCand_TRG[i].phi_<<endl;
			cout<<"      Et "<<setprecision(5)<<CutFlowCand_TRG[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<CutFlowCand_TRG[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<CutFlowCand_TRG[i].f51_<<endl;

		}
		
		cout<<endl;	
		cout<<"Print out photon information"<<endl;
	//	cout<<"nPhoton "<<nPhoton<<endl;
	//	cout<<"Highphoton size "<<HighPtPhoton.size()<<endl;
		for(int j=0;j<HighPtPhoton.size();j++){
			cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                	cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                	cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
		}

		//to see monopole is photon-like or spike-like
		//define: photon-monopole angle<0.5 is photon-like
		for(int j=0;j<HighPtPhoton.size();j++){
			double photonMono_deltaR = 0;
			photonMono_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-HighPtPhoton[j].pho_eta_,2)+
                	               	         pow(CutFlowCand_Dedx[i].phi_-HighPtPhoton[j].pho_phi_,2));
			if(photonMono_deltaR<0.15){
				cout<<"deltaR "<<photonMono_deltaR<<endl;
				cout<<"photon-like monopole"<<endl;
				photonLike++;// otherwise not photon-like monopole
			}
		}
		cout<<endl;
		cout<<"=================================="<<endl;
		}
	}
*/	//Cut all
	if(CutFlowCand_Dedx.size()>0)
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




	///////////////////////////////////////////////////
	/////////  N1 Cut Plots and  Count   //////////////
	///////////////////////////////////////////////////


	//count n_1Plot TRG+Qual and get plots
        sort(N1CutCand_Qual.begin(),N1CutCand_Qual.begin()+N1CutCand_Qual.size());
	if(N1CutCand_Qual.size()>0) 
	{
	        PlotSet &z = n_1Plot[0];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_Qual[0].subHits_,N1CutCand_Qual[0].subSatHits_/N1CutCand_Qual[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_Qual[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_Qual[0].rzp2_);
        	z.GetPlot(E55)->Fill(N1CutCand_Qual[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_Qual[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_Qual[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_Qual[0].f51_,N1CutCand_Qual[0].dEdXSig_);
		NoQual++;

	}
        
	//except E
	sort(N1CutCand_Energy.begin(),N1CutCand_Energy.begin()+N1CutCand_Energy.size());
	if(N1CutCand_Energy.size()>0)
	{
	        PlotSet &z = n_1Plot[1];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_Energy[0].subHits_,N1CutCand_Energy[0].subSatHits_/N1CutCand_Energy[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_Energy[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_Energy[0].rzp2_);
        	z.GetPlot(E55)->Fill(N1CutCand_Energy[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_Energy[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_Energy[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_Energy[0].f51_,N1CutCand_Energy[0].dEdXSig_);
		NoE++;	
	}

	//except f51
        sort(N1CutCand_F51.begin(),N1CutCand_F51.begin()+N1CutCand_F51.size());
	if(N1CutCand_F51.size()>0)
	{
	        PlotSet &z = n_1Plot[2];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_F51[0].subHits_,N1CutCand_F51[0].subSatHits_/N1CutCand_F51[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_F51[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_F51[0].rzp2_);
        	z.GetPlot(E55)->Fill(N1CutCand_F51[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_F51[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_F51[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_F51[0].f51_,N1CutCand_F51[0].dEdXSig_);
		NoF51++;	
	}
	
	//except dedx
	sort(N1CutCand_Dedx.begin(),N1CutCand_Dedx.begin()+N1CutCand_Dedx.size());
	if(N1CutCand_Dedx.size()>0)
	{
	
	        PlotSet &z = n_1Plot[3];
	        z.GetPlot(FracSatVNstrips)->Fill(N1CutCand_Dedx[0].subHits_,N1CutCand_Dedx[0].subSatHits_/N1CutCand_Dedx[0].subHits_);
	        z.GetPlot(DedXSig)->Fill(N1CutCand_Dedx[0].dEdXSig_);
	        z.GetPlot(RZcurv)->Fill(N1CutCand_Dedx[0].rzp2_);
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
        	z.GetPlot(E55)->Fill(N1CutCand_TRG[0].e55_);
	        z.GetPlot(F51)->Fill(N1CutCand_TRG[0].f51_);
        	z.GetPlot(HcalIso)->Fill(N1CutCand_TRG[0].hIso_);
        	z.GetPlot(ABCD)->Fill(N1CutCand_TRG[0].f51_,N1CutCand_TRG[0].dEdXSig_);
		NoTRG++;
	}
	//cout<<"do analysis end"<<endl;	
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
        cout<<"n-1cut pass writeplots func"<<endl;

	for(unsigned c=0; c<nCut; c++) CutFlow[c].WritePlot();
	cout<<"cutflow pass writeplots func"<<endl;
	
  }
  
  void SignalEff(const string trName)
  {
        //signal efficiency = no. of events after all selection cuts/all events
	cout<<"2000GeV"<<endl;
	cout<<trName<<" ================================="<<endl;
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

void MonoAnalyzerPhoton_2000()
{
	cout<<"hello"<<endl;
	TFile *oFile = new TFile("Monopole_2018_Photon_2000.root","recreate");
	cout<<"new file line pass"<<endl;

	TFile *fin = new TFile("../MonoNtuple2018_MC_2000_0722.root");
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
//		cout<<"nPhoton "<<nPhoton<<endl;	
//		cout<<"photon size "<<photon.size()<<endl;

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
//		cout<<"enter photon loop"<<endl;
		if(nPhoton!=0){
		for(unsigned j=0;j<nPhoton;j++){
/*			cout<<"pho eta "<<(*pho_eta)[j]<<endl;
                        cout<<"pho phi"<<(*pho_phi)[j]<<endl;
                        cout<<"pho pt "<<(*pho_pt)[j]<<endl;
*/			photon[j] = Photon(
			(*pho_eta)[j],
			(*pho_phi)[j],
			(*pho_pt)[j]
			);
			}
		}
			noTrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,true,ev);			
                        TrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon200,ev);
	}//for every event loop
	TrgAnalysis.WritePlots(oFile);
	TrgAnalysis.SignalEff("HLT_Photon200");

	noTrgAnalysis.WritePlots(oFile);
	noTrgAnalysis.SignalEff("NO TRG");
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
