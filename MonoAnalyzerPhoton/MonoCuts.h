#ifndef _MONOCUTS_H_
#define _MONOCUTS_H_
#include "Candidate.h"
#include "PlotSet.h"
class MonoCuts:public MonoCandidate, public Photon
{
public:
  MonoCuts(){}

  MonoCuts(const string &trName,TFile *openFile):trigName_(trName){
	
	//create a new directory in the output file
	openFile->mkdir(trName.c_str());
	cout<<"create directory "<<trName.c_str()<<endl;	
	//create histrogram N-1 cutf for No trigger
	
	NoCutPlot.resize(1U);
	PlotSet &x = NoCutPlot[0];
        x.CreatPlot(FracSatVNstrips,new TH2D("FracSatVNstrips","",100,0,1000,100,0,1));
        x.CreatPlot(DedXSig,new TH1D("DedXSig","",100,0,30));
        x.CreatPlot(XYPar0,new TH1D("XYPar0","",100,-1,1));
        x.CreatPlot(XYPar2,new TH1D("XYPar2","",100,-2000,2000));
        x.CreatPlot(RZPar0,new TH1D("RZPar0","",100,-20,20));
        x.CreatPlot(RZPar1,new TH1D("RZPar1","",100,-1,100));
        x.CreatPlot(RZcurv,new TH1D("RZcurv","",100,-0.01,0.01));
	x.CreatPlot(E55,new TH1D("E55","",100,-1,1200));
        x.CreatPlot(F51,new TH1D("F51","",100,0.2,1.1));
        x.CreatPlot(HcalIso,new TH1D("HcalIso","",100,-1,10));
        x.CreatPlot(ABCD,new TH2D("ABCD","",100,0,1.1,100,0,30));
	 
	cutName_[0] = "Quality_";
	 cutName_[1] = "Energy_";
	 cutName_[2] = "F51_";
	 cutName_[3] = "dEdxSig_";
	 cutName_[4] = "HLT_";
	
        n_1Plot.resize(nCut);

	for( int c = 0; c <nCut ;c++){

	string cutn1name = "N1_"+ cutName_[c];	//N_1_HLT_Quality_FracSatVNstrips
	PlotSet &z = n_1Plot[c];//[1] = Quality,[2]energy, [3] f51, [4]dEdx
	
	z.CreatPlot(FracSatVNstrips,new TH2D((cutn1name+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
	z.CreatPlot(DedXSig,new TH1D((cutn1name+"DedXSig").c_str(),"",100,0,30));
        z.CreatPlot(XYPar0,new TH1D((cutn1name+"XYPar0").c_str(),"",100,-1,1));
        z.CreatPlot(XYPar2,new TH1D((cutn1name+"XYPar2").c_str(),"",100,-2000,2000));
        z.CreatPlot(RZPar0,new TH1D((cutn1name+"RZPar0").c_str(),"",100,-20,20));
        z.CreatPlot(RZPar1,new TH1D((cutn1name+"RZPar1").c_str(),"",100,-1,100));
	z.CreatPlot(RZcurv,new TH1D((cutn1name+"RZcurv").c_str(),"",100,-0.01,0.01));
        z.CreatPlot(E55,new TH1D((cutn1name+"E55").c_str(),"",100,-1,1200));
        z.CreatPlot(F51,new TH1D((cutn1name+"F51").c_str(),"",100,0.2,1.1));
        z.CreatPlot(HcalIso,new TH1D((cutn1name+"HcalIso").c_str(),"",100,-1,10));
        z.CreatPlot(ABCD,new TH2D((cutn1name+"ABCD").c_str(),"",100,0,1.1,100,0,30));
	}

	CutFlow.resize(nCut); 
	for( int c = 0;c<nCut-1;c++){
	// nCut-1 is for ignoring HLT_ ,omit extra plots	
	PlotSet &y = CutFlow[c];
	string cutflowName = "Flow_"+cutName_[c];
        y.CreatPlot(FracSatVNstrips,new TH2D((cutflowName+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
        y.CreatPlot(DedXSig,new TH1D((cutflowName+"DedXSig").c_str(),"",100,0,30));
        y.CreatPlot(XYPar0,new TH1D((cutflowName+"XYPar0").c_str(),"",100,-1,1));
        y.CreatPlot(XYPar2,new TH1D((cutflowName+"XYPar2").c_str(),"",100,-2000,2000));
        y.CreatPlot(RZPar0,new TH1D((cutflowName+"RZPar0").c_str(),"",100,-20,20));
        y.CreatPlot(RZPar1,new TH1D((cutflowName+"RZPar1").c_str(),"",100,-1,100));
        y.CreatPlot(RZcurv,new TH1D((cutflowName+"RZcurv").c_str(),"",100,-0.01,0.01));
        y.CreatPlot(E55,new TH1D((cutflowName+"E55").c_str(),"",100,-1,1200));
        y.CreatPlot(F51,new TH1D((cutflowName+"F51").c_str(),"",100,0.2,1.1));
        //y.CreatPlot(Spike,new TH1D((cutflowName+"Spike").c_str(),"",100,-5,5));
        y.CreatPlot(HcalIso,new TH1D((cutflowName+"HcalIso").c_str(),"",100,-1,10));
        y.CreatPlot(ABCD,new TH2D((cutflowName+"ABCD").c_str(),"",100,0,1.1,100,0,30));
        cout<<"good cutflow"<<cutflowName<<endl;
	}

  }

  ~MonoCuts(){}
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev);
  void FillNoCutHistogram(int n,vector<MonoCandidate> Cand);
  void FillN1Histogram(int n, vector<MonoCandidate> N1CutCand);
  void FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand);
  void Matching(vector<MonoCandidate> CutFlowCand);
  void Clear();
  void WritePlots(TFile *oFile);
  void PrintInfo(){	
        
	cout<<"event tag"<<CutFlowCand_TRG[0].event_<<endl;
	for(int i=0;i<CutFlowCand_TRG.size();i++){
    	        //cout<<"          candidate           monoGen         antimonoGen"<<endl;
    	   	//cout<<"eta      "<<setprecision(5)<<CutFlowCand_Dedx[i].eta_<<setw(20)<<CutFlowCand_Dedx[i].mono_eta_
    	   	//		 <<setw(20)<<CutFlowCand_Dedx[i].amon_eta_<<endl;
    	   	//cout<<"phi      "<<setprecision(5)<<CutFlowCand_Dedx[i].phi_<<setw(20)<<CutFlowCand_Dedx[i].mono_phi_
    	   	//		 <<setw(20)<<CutFlowCand_Dedx[i].amon_phi_<<endl;
    	   	//cout<<"m deltaR "<<m_deltaR<<endl;
    	   	//cout<<"a deltaR "<<am_deltaR<<endl;
		cout<<i<<endl;
		cout<<"     eta "<<setprecision(5)<<CutFlowCand_TRG[i].eta_<<endl;
		cout<<"     phi "<<setprecision(5)<<CutFlowCand_TRG[i].phi_<<endl;
		cout<<"     E55 "<<setprecision(5)<<CutFlowCand_TRG[i].e55_<<endl;
		cout<<" dEdxSig "<<setprecision(5)<<CutFlowCand_TRG[i].dEdXSig_<<endl;
		cout<<"     f51 "<<setprecision(5)<<CutFlowCand_TRG[i].f51_<<endl;
		cout<<"   Swiss "<<setprecision(5)<<CutFlowCand_TRG[i].Cross_<<endl;
		cout<<"----------------------------"<<endl;
	}
  }
  void PrintPhoton(){	
	//Just Print out all photon info in  one event
	for(int j=0;j<HighPtPhoton.size();j++){
		cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
        	cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
        	cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
	}

	//define: photon-monopole angle<0.15 is photon-like
//	for(int i=0;i<HighPtPhoton.size();i++){
//		double photonMono_deltaR = 0;
//		photonMono_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-HighPtPhoton[j].pho_eta_,2)+
//        	               	         pow(CutFlowCand_Dedx[i].phi_-HighPtPhoton[j].pho_phi_,2));
//		if(photonMono_deltaR<0.15){
//			cout<<"deltaR "<<photonMono_deltaR<<endl;
//			cout<<"photon-like monopole"<<endl;
//			photonLike++;// otherwise not photon-like monopole
//		}
//	}
	cout<<endl;
	cout<<"=================================="<<endl;
  }
  void SignalEff(const string trName, double TotalEvents);

  bool operator<(const MonoCandidate &mc)const{
   if(dEdXSig_>mc.dEdXSig_) return true;
   else if(dEdXSig_==mc.dEdXSig_){
        if(f51_>mc.f51_) return true;
        else return false;
        }
    else return false;
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
  // No Cut plot
  vector<PlotSet> NoCutPlot;
  //N-1 cut plot box
  vector<PlotSet> n_1Plot;
  //cutflow plot box
  vector<PlotSet> CutFlow;
  static const unsigned nCut = 5U;
  string cutName_[nCut];
  


  //cuts analysis
  bool evalQuality(MonoCandidate &mc) { return TMath::Abs(mc.xyp0_) < xyp0Cut_ && TMath::Abs(mc.xyp2_) > xyp2Cut_ 
			&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_  
			&& TMath::Abs( mc.rzp2_) < rzp2Cut_ && TMath::Abs(mc.rzp1_) < rzp1Cut_ && TMath::Abs(mc.rzp0_) < rzp0Cut_;  }
  bool evalE(MonoCandidate &mc) { return mc.e55_ > e55Cut_; }
  bool evalF51(MonoCandidate &mc) { return mc.f51_ > f51Cut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
 
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
  int SpikeLike=0;
  int EBarrel=0; 
  int Reco=0;
  int photonLike=0;
  //relative eff without HLT
  int NoTRG=0;
  int NoQual=0;
  int NoE=0;
  int NoF51=0;
  int NodEdXCut=0;
};

#endif
