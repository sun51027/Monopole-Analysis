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
#include "Candidate.h"
#include "PlotSet.h"
#include <algorithm>
#include <string>
using namespace std;
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
        x.CreatPlot(Dist,new TH1D("Dist","",100,-1,1));
        x.CreatPlot(HcalIso,new TH1D("HcalIso","",100,-1,10));
        x.CreatPlot(XYPar0,new TH1D("XYPar0","",100,-3,3));
        x.CreatPlot(XYPar1,new TH1D("XYPar1","",100,-10,10));
        x.CreatPlot(XYPar2,new TH1D("XYPar2","",100,-20000,20000));
        x.CreatPlot(RZPar0,new TH1D("RZPar0","",100,-20,20));
        x.CreatPlot(RZPar1,new TH1D("RZPar1","",100,-20,20));
        x.CreatPlot(RZcurv,new TH1D("RZcurv","",100,-0.01,0.01));
	x.CreatPlot(E55,new TH1D("E55","",100,-1,1200));
        x.CreatPlot(F51,new TH1D("F51","",100,0.2,1.1));
        x.CreatPlot(ABCD,new TH2D("ABCD","",100,0,1.1,100,0,30));
	 
  	 cutName_[0] = "Dist_";
	 cutName_[1] = "HIso_";
	 cutName_[2] = "XYPar0_";
	 cutName_[3] = "XYPar1_";
	 cutName_[4] = "XYPar2_";
	 cutName_[5] = "RZPar0_";
	 cutName_[6] = "RZPar1_";
	 cutName_[7] = "RZcurv_";
	 cutName_[8] = "Energy_";
	 cutName_[9] = "e55_";
	 cutName_[10] = "Dedx_";
	
        n_1Plot.resize(nCut);
//        string trgn1name = "HLT_";
        string name = "N_1_";         //  N_1_HLT_FracSatVNstrips
	cout<<"nCut "<<nCut<<endl;
/*	for( int c = 0; c <nCut ;c++){

	string cutn1name = name + cutName_[c];	//N_1_HLT_Quality_FracSatVNstrips
	PlotSet &z = n_1Plot[c];//[1] = Quality,[2]energy, [3] f51, [4]dEdx
	
	z.CreatPlot(FracSatVNstrips,new TH2D((cutn1name+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
	z.CreatPlot(DedXSig,new TH1D((cutn1name+"DedXSig").c_str(),"",100,0,30));
        z.CreatPlot(Dist,new TH1D((cutn1name+"Dist").c_str(),"",100,-1,1));
        z.CreatPlot(HcalIso,new TH1D((cutn1name+"HcalIso").c_str(),"",100,-1,10));
        z.CreatPlot(XYPar0,new TH1D((cutn1name+"XYPar0").c_str(),"",100,-3,3));
        z.CreatPlot(XYPar1,new TH1D((cutn1name+"XYPar1").c_str(),"",100,-10,10));
        z.CreatPlot(XYPar2,new TH1D((cutn1name+"XYPar2").c_str(),"",100,-2000,2000));
        z.CreatPlot(RZPar0,new TH1D((cutn1name+"RZPar0").c_str(),"",100,-20,20));
        z.CreatPlot(RZPar1,new TH1D((cutn1name+"RZPar1").c_str(),"",100,-20,20));
	z.CreatPlot(RZcurv,new TH1D((cutn1name+"RZcurv").c_str(),"",100,-0.01,0.01));
        z.CreatPlot(E55,new TH1D((cutn1name+"E55").c_str(),"",100,-1,1200));
        z.CreatPlot(F51,new TH1D((cutn1name+"F51").c_str(),"",100,0.2,1.1));
        z.CreatPlot(ABCD,new TH2D((cutn1name+"ABCD").c_str(),"",100,0,1.1,100,0,30));
	}
	cout<<"create N_1 plots"<<endl;*/
	//create plot with cutflow
	CutFlow.resize(nCut); 
	cout<<"nCut "<<nCut<<endl;
	for( int c = 0;c<nCut;c++){
	
	PlotSet &y = CutFlow[c];
	string cutflowName = "Flow_"+cutName_[c];
	//[c=1] Flow_HLT_Energy_	
        y.CreatPlot(FracSatVNstrips,new TH2D((cutflowName+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
        y.CreatPlot(DedXSig,new TH1D((cutflowName+"DedXSig").c_str(),"",100,0,30));
        y.CreatPlot(Dist,new TH1D((cutflowName+"Dist").c_str(),"",100,-1,1));
        y.CreatPlot(HcalIso,new TH1D((cutflowName+"HcalIso").c_str(),"",100,-1,10));
        y.CreatPlot(XYPar0,new TH1D((cutflowName+"XYPar0").c_str(),"",100,-3,3));
        y.CreatPlot(XYPar1,new TH1D((cutflowName+"XYPar1").c_str(),"",100,-10,10));
        y.CreatPlot(XYPar2,new TH1D((cutflowName+"XYPar2").c_str(),"",100,-20000,20000));
        y.CreatPlot(RZPar0,new TH1D((cutflowName+"RZPar0").c_str(),"",100,-20,20));
        y.CreatPlot(RZPar1,new TH1D((cutflowName+"RZPar1").c_str(),"",100,-20,20));
        y.CreatPlot(RZcurv,new TH1D((cutflowName+"RZcurv").c_str(),"",100,-0.01,0.01));
        y.CreatPlot(E55,new TH1D((cutflowName+"E55").c_str(),"",100,-1,1200));
        y.CreatPlot(F51,new TH1D((cutflowName+"F51").c_str(),"",100,0.2,1.1));
        y.CreatPlot(Spike,new TH1D((cutflowName+"Spike").c_str(),"",100,-5,5));
        y.CreatPlot(ABCD,new TH2D((cutflowName+"ABCD").c_str(),"",100,0,1.1,100,0,30));
	}
	cout<<"create cutflow plots"<<endl;

  }

  ~MonoCuts(){}
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev)
  {
	Clear();
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
	        PlotSet &x = Plot[0];
		
           for(int i=0; i<CutFlowCand_TRG.size();i++){

                double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(CutFlowCand_TRG[i].eta_-CutFlowCand_TRG[0].mono_eta_,2)+
                                pow(CutFlowCand_TRG[i].phi_-CutFlowCand_TRG[0].mono_phi_,2));
                am_deltaR= sqrt(pow(CutFlowCand_TRG[i].eta_-CutFlowCand_TRG[0].amon_eta_,2)+
                                pow(CutFlowCand_TRG[i].phi_-CutFlowCand_TRG[0].amon_phi_,2));

                if(m_deltaR<0.15||am_deltaR<0.15){

	        x.GetPlot(FracSatVNstrips)->Fill(CutFlowCand_TRG[i].subHits_,CutFlowCand_TRG[i].subSatHits_/CutFlowCand_TRG[i].subHits_);
	        x.GetPlot(DedXSig)->Fill(CutFlowCand_TRG[i].dEdXSig_);
	        x.GetPlot(Dist)->Fill(CutFlowCand_TRG[i].dist_);
	        x.GetPlot(XYPar0)->Fill(CutFlowCand_TRG[i].xyp0_);
	        x.GetPlot(XYPar1)->Fill(CutFlowCand_TRG[i].xyp1_);
	        x.GetPlot(XYPar2)->Fill(CutFlowCand_TRG[i].xyp2_);
	        x.GetPlot(RZPar0)->Fill(CutFlowCand_TRG[i].rzp0_);
	        x.GetPlot(RZPar1)->Fill(CutFlowCand_TRG[i].rzp1_);
	        x.GetPlot(RZcurv)->Fill(CutFlowCand_TRG[i].rzp2_);
        	x.GetPlot(HcalIso)->Fill(CutFlowCand_TRG[i].hIso_);
        	x.GetPlot(E55)->Fill(CutFlowCand_TRG[i].e55_);
	        x.GetPlot(F51)->Fill(CutFlowCand_TRG[i].f51_);
        	x.GetPlot(ABCD)->Fill(CutFlowCand_TRG[i].f51_,CutFlowCand_TRG[i].dEdXSig_);
           	}
	   }
		count++;

	}
        sort(CutFlowCand_Dist.begin(),CutFlowCand_Dist.begin()+CutFlowCand_Dist.size());
	
	if(CutFlowCand_Dist.size()>0) 
	{
		Matching(0,CutFlowCand_Dist);
		Dist_count++;
	}

	//count cutflow TRG+Qual and get plots
        sort(CutFlowCand_HIso.begin(),CutFlowCand_HIso.begin()+CutFlowCand_HIso.size());
	if(CutFlowCand_HIso.size()>0) 
	{
		HIso_count++;	
		Matching(1,CutFlowCand_HIso);
	}
	
	//count cutflow TRG+Qual+Energy and get plots
        sort(CutFlowCand_XYPar0.begin(),CutFlowCand_XYPar0.begin()+CutFlowCand_XYPar0.size());
	if(CutFlowCand_XYPar0.size()>0)
	{
		XYPar0_count++;	
		Matching(2,CutFlowCand_XYPar0);
	}

	//TRG+Qual+Energy+F51
        sort(CutFlowCand_XYPar1.begin(),CutFlowCand_XYPar1.begin()+CutFlowCand_XYPar1.size());
	if(CutFlowCand_XYPar1.size()>0)
	{
		XYPar1_count++;	
		Matching(3,CutFlowCand_XYPar1);
	}
        sort(CutFlowCand_XYPar2.begin(),CutFlowCand_XYPar2.begin()+CutFlowCand_XYPar2.size());
	if(CutFlowCand_XYPar2.size()>0)
	{
		XYPar2_count++;	
		Matching(4,CutFlowCand_XYPar2);
	}
	
	sort(CutFlowCand_RZPar0.begin(),CutFlowCand_RZPar0.begin()+CutFlowCand_RZPar0.size());
	if(CutFlowCand_RZPar0.size()>0)
	{
		RZPar0_count++;	
		Matching(5,CutFlowCand_RZPar0);
	}
	sort(CutFlowCand_RZPar1.begin(),CutFlowCand_RZPar1.begin()+CutFlowCand_RZPar1.size());
	if(CutFlowCand_RZPar1.size()>0)
	{
		RZPar1_count++;	
		Matching(6,CutFlowCand_RZPar1);
	}
	sort(CutFlowCand_RZPar2.begin(),CutFlowCand_RZPar2.begin()+CutFlowCand_RZPar2.size());
	if(CutFlowCand_RZPar2.size()>0)
	{
		RZPar2_count++;	
		Matching(7,CutFlowCand_RZPar2);
	}
	sort(CutFlowCand_Energy.begin(),CutFlowCand_Energy.begin()+CutFlowCand_Energy.size());
	if(CutFlowCand_Energy.size()>0) {
		E_count++;
		Matching(8,CutFlowCand_Energy);
	}
	sort(CutFlowCand_F51.begin(),CutFlowCand_F51.begin()+CutFlowCand_F51.size());
	if(CutFlowCand_F51.size()>0) {
		f51_count++;
		Matching(9,CutFlowCand_F51);
	}
	sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());
	if(CutFlowCand_Dedx.size()>0)
	{
		Dedx_count++;	
		Matching(10,CutFlowCand_Dedx);
	}
	/////////(/////////////////////////////////////////
	/////////  N1 Cut Plots and  Count   //////////////
	///////////////////////////////////////////////////

        sort(N1CutCand_Dist.begin(),N1CutCand_Dist.begin()+N1CutCand_Dist.size());
	if(N1CutCand_Dist.size()>0) NoDist++;
	sort(N1CutCand_HIso.begin(),N1CutCand_HIso.begin()+N1CutCand_HIso.size());
	if(N1CutCand_HIso.size()>0) NoHIso++;	
        sort(N1CutCand_XYPar0.begin(),N1CutCand_XYPar0.begin()+N1CutCand_XYPar0.size());
	if(N1CutCand_XYPar0.size()>0)	NoXYPar0++;	
        sort(N1CutCand_XYPar1.begin(),N1CutCand_XYPar1.begin()+N1CutCand_XYPar1.size());
	if(N1CutCand_XYPar1.size()>0)	NoXYPar1++;	
	sort(N1CutCand_XYPar2.begin(),N1CutCand_XYPar2.begin()+N1CutCand_XYPar2.size());
	if(N1CutCand_XYPar2.size()>0)   NoXYPar2++;	
        sort(N1CutCand_RZPar0.begin(),N1CutCand_RZPar0.begin()+N1CutCand_RZPar0.size());
	if(N1CutCand_RZPar0.size()>0)   NoRZPar0++;
        sort(N1CutCand_RZPar1.begin(),N1CutCand_RZPar1.begin()+N1CutCand_RZPar1.size());
	if(N1CutCand_RZPar1.size()>0) 	NoRZPar1++;
        sort(N1CutCand_RZPar2.begin(),N1CutCand_RZPar2.begin()+N1CutCand_RZPar2.size());
	if(N1CutCand_RZPar2.size()>0)   NoRZPar2++;
        sort(N1CutCand_TRG.begin(),N1CutCand_TRG.begin()+N1CutCand_TRG.size());
	if(N1CutCand_TRG.size()>0) NoTRG++;

  }
  bool operator<(const MonoCandidate &mc)const{
   if(dEdXSig_>mc.dEdXSig_) return true;
   else if(dEdXSig_==mc.dEdXSig_){
        if(f51_>mc.f51_) return true;
        else return false;
        }
    else return false;
  }

  void FillN1Histogram(int n, vector<MonoCandidate> N1CutCand){
	PlotSet &z = n_1Plot[n];
	for(int i=0; i < N1CutCand.size() ;i++){
	  z.GetPlot(FracSatVNstrips)->Fill(N1CutCand[i].subHits_,N1CutCand[i].subSatHits_/N1CutCand[i].subHits_);
	  z.GetPlot(DedXSig)->Fill(N1CutCand[i].dEdXSig_);
	  z.GetPlot(RZcurv)->Fill(N1CutCand[i].rzp2_);
          z.GetPlot(E55)->Fill(N1CutCand[i].e55_);
	  z.GetPlot(F51)->Fill(N1CutCand[i].f51_);
          z.GetPlot(HcalIso)->Fill(N1CutCand[i].hIso_);
          z.GetPlot(ABCD)->Fill(N1CutCand[i].f51_,N1CutCand[i].dEdXSig_);
	}
  }
  void FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand){
	PlotSet &z =CutFlow[n];
	for(int i=0; i < CutFlowCand.size() ;i++){
	  z.GetPlot(FracSatVNstrips)->Fill(CutFlowCand[i].subHits_,CutFlowCand[i].subSatHits_/CutFlowCand[i].subHits_);
	  z.GetPlot(DedXSig)->Fill(CutFlowCand[i].dEdXSig_);
	  z.GetPlot(XYPar0)->Fill(CutFlowCand[i].xyp0_);
	  z.GetPlot(XYPar1)->Fill(CutFlowCand[i].xyp1_);
	  z.GetPlot(XYPar2)->Fill(CutFlowCand[i].xyp2_);
	  z.GetPlot(RZPar0)->Fill(CutFlowCand[i].rzp1_);
	  z.GetPlot(RZPar1)->Fill(CutFlowCand[i].rzp0_);
	  z.GetPlot(RZcurv)->Fill(CutFlowCand[i].rzp2_);
          z.GetPlot(E55)->Fill(CutFlowCand[i].e55_);
	  z.GetPlot(F51)->Fill(CutFlowCand[i].f51_);
          z.GetPlot(HcalIso)->Fill(CutFlowCand[i].hIso_);
          z.GetPlot(ABCD)->Fill(CutFlowCand[i].f51_,CutFlowCand[i].dEdXSig_);
	}
  }
  void  Matching(int n, vector<MonoCandidate> CutFlowCand){
	PlotSet &z =CutFlow[n];
	   for(int i=0; i<CutFlowCand.size();i++){

		double m_deltaR=0;
		double am_deltaR=0;
		m_deltaR = sqrt(pow(CutFlowCand[i].eta_-CutFlowCand[0].mono_eta_,2)+
				pow(CutFlowCand[i].phi_-CutFlowCand[0].mono_phi_,2));
		am_deltaR= sqrt(pow(CutFlowCand[i].eta_-CutFlowCand[0].amon_eta_,2)+
                                pow(CutFlowCand[i].phi_-CutFlowCand[0].amon_phi_,2));

		if(m_deltaR<0.15 || am_deltaR<0.15){
	           z.GetPlot(FracSatVNstrips)->Fill(CutFlowCand[i].subHits_,CutFlowCand[i].subSatHits_/CutFlowCand[i].subHits_);
	           z.GetPlot(DedXSig)->Fill(CutFlowCand[i].dEdXSig_);
	           z.GetPlot(XYPar0)->Fill(CutFlowCand[i].xyp0_);
	           z.GetPlot(XYPar1)->Fill(CutFlowCand[i].xyp1_);
	           z.GetPlot(XYPar2)->Fill(CutFlowCand[i].xyp2_);
	           z.GetPlot(RZPar0)->Fill(CutFlowCand[i].rzp1_);
	           z.GetPlot(RZPar1)->Fill(CutFlowCand[i].rzp0_);
	           z.GetPlot(RZcurv)->Fill(CutFlowCand[i].rzp2_);
                   z.GetPlot(E55)->Fill(CutFlowCand[i].e55_);
	           z.GetPlot(F51)->Fill(CutFlowCand[i].f51_);
                   z.GetPlot(HcalIso)->Fill(CutFlowCand[i].hIso_);
                   z.GetPlot(ABCD)->Fill(CutFlowCand[i].f51_,CutFlowCand[i].dEdXSig_);
		}	
	   }
	
		
  }
  void Clear(){

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
  }
  void WritePlots(TFile *oFile){
	oFile->cd(trigName_.c_str());
	Plot[0].WritePlot();
	for(unsigned c=0; c<nCut; c++) n_1Plot[c].WritePlot();
        cout<<"n-1cut pass writeplots func"<<endl;

	for(unsigned c=0; c<nCut; c++) CutFlow[c].WritePlot();
	cout<<"cutflow pass writeplots func"<<endl;
	
  }
  
  void PrintInfo(){	
//		cout<<"          candidate           monoGen         antimonoGen"<<endl;
//		cout<<"eta      "<<setprecision(5)<<CutFlowCand_Dedx[i].eta_<<setw(20)<<CutFlowCand_Dedx[i].mono_eta_
//				 <<setw(20)<<CutFlowCand_Dedx[i].amon_eta_<<endl;
//		cout<<"phi      "<<setprecision(5)<<CutFlowCand_Dedx[i].phi_<<setw(20)<<CutFlowCand_Dedx[i].mono_phi_
//				 <<setw(20)<<CutFlowCand_Dedx[i].amon_phi_<<endl;
//		cout<<"m deltaR "<<m_deltaR<<endl;
//		cout<<"a deltaR "<<am_deltaR<<endl;
//                cout<<"ev "<<ev<<" has monopole"<<endl;
//                cout<<"event tag"<<CutFlowCand_Dedx[0].event_<<endl;
//
		for(int i=0;i<CutFlowCand_Dedx.size();i++){
			cout<<i<<endl;
			cout<<"     eta "<<setprecision(5)<<CutFlowCand_Dedx[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<CutFlowCand_Dedx[i].phi_<<endl;
			cout<<"      Et "<<setprecision(5)<<CutFlowCand_Dedx[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<CutFlowCand_Dedx[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<CutFlowCand_Dedx[i].f51_<<endl;
			cout<<"   Swiss "<<setprecision(5)<<CutFlowCand_Dedx[i].Cross_<<endl;
			cout<<"----------------------------"<<endl;
		}
  }

  void PrintPhoton(){	
		cout<<"Print out photon information"<<endl;
		for(int j=0;j<HighPtPhoton.size();j++){
			cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                	cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                	cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
		}
		cout<<endl;
		cout<<"=================================="<<endl;
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
/*	cout<<"Relative effciency count"<<endl;
        cout<<"------------------------------------------"<<endl;
	cout<<"     No TRG "<<NoTRG<<endl;
	cout<<"       Dist "<<NoDist<<endl;
	cout<<"	      HIso "<<NoHIso<<endl;
	cout<<"     XYPar0 "<<NoXYPar0<<endl;
	cout<<"     XYPar0 "<<NoXYPar1<<endl;
	cout<<"     XYPar2 "<<NoXYPar2<<endl;
	cout<<"     RZPar0 "<<NoRZPar0<<endl;
	cout<<"     RZPar1 "<<NoRZPar1<<endl;
	cout<<"     RZPar2 "<<NoRZPar2<<endl;*/
//        cout<<" No Quality "<<NoQual<<endl;
//        cout<<"    No ECut "<< NoE <<endl;
//        cout<<"  No F51Cut "<<NoF51<<endl;
//        cout<<" No dEdXSig "<<NodEdXCut<<endl;
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
  static const unsigned nCut = 11U;
  string cutName_[nCut];
  string VectorName[20];


  //cuts analysis
  bool evalQuality(MonoCandidate &mc) { return TMath::Abs(mc.xyp0_ )< xyp0Cut_/* && TMath::Abs(mc.xyp1_ )< xyp1Cut_*/
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
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
  int TotalEvents = 15000;
 
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

void MonoAnalyzerPhoton_TrackParameter()
{
	TFile *oFile = new TFile("MonoPhotonAnalysis_2018_1000_TrackParameter_withMatching_XYPar2inverse.root","recreate");

	TFile *fin = new TFile("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MonoNtuple2018/1000/MonoNtuple2018_MC_1000.root");
        TTree *tree = (TTree*)fin->Get("monopoles");
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

	vector<double> * Cross = 0;
	unsigned event;
	unsigned NPV;
	
        double mono_eta;
        double mono_phi;
        double amon_eta;
        double amon_phi;

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
        tree->SetBranchAddress("cand_SwissCross",&Cross);

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
	
	
	MonoCuts noTrgAnalysis("NOTRG",oFile);
	MonoCuts TrgAnalysis("HLT_Photon200",oFile);

	vector<MonoCandidate> cand(10);	
	vector<Photon> photon(0);
	
        for(unsigned ev=0; ev<NEvents;ev++){
//             	cout<<"ev "<<ev<<endl;
		tree->GetEntry(ev);
		
		if(nCandidates>Et.size()) Et.resize(nCandidates);
		if(nCandidates>cand.size()) cand.resize(nCandidates);

		if(nPhoton>photon.size()) photon.resize(nPhoton);

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
		(*Cross)[i],
//		Et[i],
	        (*e55)[i],
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
			photon[j] = Photon(
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
	noTrgAnalysis.SignalEff("NOTRG");
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
