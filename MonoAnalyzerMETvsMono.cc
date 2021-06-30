//////////////////////////////////
//
//	File name: MonoAnalyzerMETvsMono.cc
//	author: Shih Lin
//	Content: Matching
//		 Photon information(21/4/24 update) 
//		 MET / MET trigger path(21/5/11 update)
//	
//	purpose: Ecal spike-like / shower-like 
//		 trk(dedx quality)
//		 GenMET PFMET CaloMET
//	Output: Monopole_MET_analysis_*000GeV.root
//	(using MonoDrawPlots.cc to complete the plots)
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
  PFMETvsEcalTwoShower=0,//2
  PFMETvsEcalOneShowerOneSpike,//2
  PFMETvsEcalTwoSpike,//2
  PFMETvsTrkTwo,//2
  PFMETvsEcalOneShower_TrkOne,//2
  PFMETvsEcalOneSpike_TrkOne,//2
  PFMETvsEcalOneShower,//1
  PFMETvsEcalOneSpike,//1
  PFMETvsTrkOne,//1
  PFMETvsZeroMono,//0

 GenMETvsEcalTwoShower,//2
 GenMETvsEcalOneShowerOneSpike,//2
 GenMETvsEcalTwoSpike,//2
 GenMETvsTrkTwo,//2
 GenMETvsEcalOneShower_TrkOne,//2
 GenMETvsEcalOneSpike_TrkOne,//2
 GenMETvsEcalOneShower,//1
 GenMETvsEcalOneSpike,//1
 GenMETvsTrkOne,//1
 GenMETvsZeroMono,//0
  
  CaloMETvsEcalTwoShower,//2
  CaloMETvsEcalOneShowerOneSpike,//2
  CaloMETvsEcalTwoSpike,//2
  CaloMETvsTrkTwo,//2
  CaloMETvsEcalOneShower_TrkOne,//2
  CaloMETvsEcalOneSpike_TrkOne,//2
  CaloMETvsEcalOneShower,//1
  CaloMETvsEcalOneSpike,//1
  CaloMETvsTrkOne,//1
  CaloMETvsZeroMono,//0
  GenMET,
  GenMETminusTwoMonopole,

  MonoPtvsTracker,
  AmonPtvsTracker,
  MonoPtvsEcal,
  AmonPtvsEcal,
  MonoEtavsSpike,
  AmonEtavsSpike,
  MonoEtavsShower,
  AmonEtavsShower

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
    double e55, double SwissCross, double hiso, double eta,
    double phi, double mono_eta, double mono_phi, double mono_px, double mono_py, double mono_p,
    double amon_eta, double amon_phi, double amon_px, double amon_py, double amon_p,
    double event,double NPV, double MET_pt, double MET_phi, double GenMET_px,double GenMET_py,double GenMET_pt,double GenMET_phi,
    double CaloMET_px,double CaloMET_py,double CaloMET_pt,double CaloMET_phi):
  subHits_(sh),subSatHits_(satsh),dEdXSig_(dedxsig),tIso_(tiso),xyp0_(xyp0),
  xyp1_(xyp1),xyp2_(xyp2),rzp0_(rzp0),rzp1_(rzp1),rzp2_(rzp2),
  dist_(dist),f51_(f51),f15_(f15),e55_(e55),SwissCross_(SwissCross),hIso_(hiso),
  eta_(eta),phi_(phi),mono_eta_(mono_eta), mono_phi_(mono_phi),mono_px_(mono_px),mono_py_(mono_py),mono_p_(mono_p),
  amon_eta_(amon_eta), amon_phi_(amon_phi),amon_px_(amon_px),amon_py_(amon_py),amon_p_(amon_p),
  event_(event),NPV_(NPV),MET_pt_(MET_pt),MET_phi_(MET_phi), GenMET_px_(GenMET_px),GenMET_py_(GenMET_py),GenMET_pt_(GenMET_pt),GenMET_phi_(GenMET_phi),
  CaloMET_px_(CaloMET_px),CaloMET_py_(CaloMET_py),CaloMET_pt_(CaloMET_pt),CaloMET_phi_(CaloMET_phi) { }
  //This will be used in comparing with cut
  MonoCandidate(const MonoCandidate &mc) : 
    subHits_(mc.subHits_),subSatHits_(mc.subSatHits_),dEdXSig_(mc.dEdXSig_),tIso_(mc.tIso_),
    xyp0_(mc.xyp0_),xyp1_(mc.xyp1_),xyp2_(mc.xyp2_),
    rzp0_(mc.rzp0_),rzp1_(mc.rzp1_),rzp2_(mc.rzp2_),
    dist_(mc.dist_),f51_(mc.f51_),f15_(mc.f15_),e55_(mc.e55_),SwissCross_(mc.SwissCross_),
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
  double SwissCross_;
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
	cout<<"create directory "<<trName.c_str()<<endl;	
	//create histrogram N-1 cutf for No trigger
	Plot.resize(1U);
	PlotSet &x = Plot[0];
 	x.CreatPlot(PFMETvsEcalTwoShower,new TH1D("PFMETvsEcalTwoShower","2 Showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(PFMETvsEcalTwoSpike,new TH1D("PFMETvsEcalTwoSpike","2 Spikelike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneShowerOneSpike,new TH1D("PFMETvsEcalOneShowerOneSpike","1 Spikelike & 1 showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(PFMETvsTrkTwo,new TH1D("PFMETvsTrkTwo","2 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneShower_TrkOne,new TH1D("PFMETvsEcalOneShower_TrkOne","1 showerlike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneSpike_TrkOne,new TH1D("PFMETvsEcalOneSpike_TrkOne","1 spikelike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneShower,new TH1D("PFMETvsEcalOneShower","1 showerlike Monopole in Ecal",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneSpike,new TH1D("PFMETvsEcalOneSpike","1 spikelike Monopole in ECal",100,0,800));	
 	x.CreatPlot(PFMETvsTrkOne,new TH1D("PFMETvsTrkOne","1 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsZeroMono,new TH1D("PFMETvsZeroMono","No monopole is reconstructed in this event",100,0,800));	
	
 	x.CreatPlot(GenMETvsEcalTwoShower,new TH1D("GenMETvsEcalTwoShower","2 Showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(GenMETvsEcalTwoSpike,new TH1D("GenMETvsEcalTwoSpike","2 Spikelike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(GenMETvsEcalOneShowerOneSpike,new TH1D("GenMETvsEcalOneShowerOneSpike","1 Spikelike & 1 showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(GenMETvsTrkTwo,new TH1D("GenMETvsTrkTwo","2 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(GenMETvsEcalOneShower_TrkOne,new TH1D("GenMETvsEcalOneShower_TrkOne","1 showerlike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(GenMETvsEcalOneSpike_TrkOne,new TH1D("GenMETvsEcalOneSpike_TrkOne","1 spikelike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(GenMETvsEcalOneShower,new TH1D("GenMETvsEcalOneShower","1 showerlike Monopole in Ecal",100,0,800));	
 	x.CreatPlot(GenMETvsEcalOneSpike,new TH1D("GenMETvsEcalOneSpike","1 spikelike Monopole in ECal",100,0,800));	
 	x.CreatPlot(GenMETvsTrkOne,new TH1D("GenMETvsTrkOne","1 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(GenMETvsZeroMono,new TH1D("GenMETvsZeroMono","No monopole is reconstructed in this event",100,0,800));	

/* 	x.CreatPlot(CaloMETvsEcalTwoShower,new TH1D("CaloMETvsEcalTwoShower","2 Showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(CaloMETvsEcalTwoSpike,new TH1D("CaloMETvsEcalTwoSpike","2 Spikelike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(CaloMETvsEcalOneShowerOneSpike,new TH1D("CaloMETvsEcalOneShowerOneSpike","1 Spikelike & 1 showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(CaloMETvsTrkTwo,new TH1D("CaloMETvsTrkTwo","2 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(CaloMETvsEcalOneShower_TrkOne,new TH1D("CaloMETvsEcalOneShower_TrkOne","1 showerlike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(CaloMETvsEcalOneSpike_TrkOne,new TH1D("CaloMETvsEcalOneSpike_TrkOne","1 spikelike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(CaloMETvsEcalOneShower,new TH1D("CaloMETvsEcalOneShower","1 showerlike Monopole in Ecal",100,0,800));	
 	x.CreatPlot(CaloMETvsEcalOneSpike,new TH1D("CaloMETvsEcalOneSpike","1 spikelike Monopole in ECal",100,0,800));	
 	x.CreatPlot(CaloMETvsTrkOne,new TH1D("CaloMETvsTrkOne","1 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(CaloMETvsZeroMono,new TH1D("CaloMETvsZeroMono","No monopole is reconstructed in this event",100,0,800));	
 	
	x.CreatPlot(GenMET,new TH1D("GenMET","GenMET",100,0,800));	
	x.CreatPlot(GenMETminusTwoMonopole,new TH1D("GenMETminusTwoMonopole","GenMET - TwoMonopole",100,0,800));	
*/
	x.CreatPlot(MonoPtvsTracker,new TH1D("MonoPtvsTracker","",100,0,2500));	
	x.CreatPlot(AmonPtvsTracker,new TH1D("AmonPtvsTracker","",100,0,2500));	
	x.CreatPlot(MonoPtvsEcal,new TH1D("MonoPtvsEcal","",100,0,2500));	
	x.CreatPlot(AmonPtvsEcal,new TH1D("AmonPtvsEcal","",100,0,2500));	

        x.CreatPlot(MonoEtavsSpike,new TH1D("MonoEtavsSpike","",100,-4,4));
        x.CreatPlot(AmonEtavsSpike,new TH1D("AmonEtavsSpike","",100,-4,4));
        x.CreatPlot(MonoEtavsShower,new TH1D("MonoEtavsShower","",100,-4,4));
        x.CreatPlot(AmonEtavsShower,new TH1D("AmonEtavsShower","",100,-4,4));


	cout<<"the address of plots x "<<&x<<endl;

	//Name the cut name
	 cutName_[0] = "Quality_";
	 cutName_[1] = "Energy_";
	 cutName_[2] = "F51_";
	 cutName_[3] = "dEdxSig_";
	 cutName_[4] = "HLT_";
	
  }

  ~MonoCuts(){}
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev, double PFMET_pt, 
  double GenMET_pt, double CaloMET_pt)
  {
	
	cout<<"enter doAnalysis successfully"<<endl;	
	EcalShower.clear();
	EcalSpike.clear();
	Trk.clear();
	Quality.clear();	

	CutFlowCand_TRG.clear();
        HighPtPhoton.clear();


	EcalSpikeNum=0;
	EcalShowerNum=0;
	TrackerMono=0;

     for(unsigned c=0;c<nCandidates;c++){
	MonoCandidate &cands = cand[c];
	bool QualCut = evalQuality(cands);
	bool ECut = evalE(cands);
	bool F15Cut = evalF15(cands);
	bool F51Cut = evalF51(cands);
	bool SwissCrossCut = evalSwissCross(cands);
	bool dEdXCut = evaldEdX(cands);
	cout<<"hi"<<endl; 
	//no cut
	if(TRG ) CutFlowCand_TRG.push_back(cands);
	cout<<"TRG"<<endl;
	//cut all 
	if(TRG && QualCut && ECut && F51Cut && dEdXCut) CutFlowCand_Dedx.push_back(cands);

	//MET
	/*if(QualCut)*/ 			 			Quality.push_back(cands); 
//	if(QualCut && dEdXCut && ECut)   			EcalShower.push_back(cands);//need to add f51<0.85 in the following code)
//	if(QualCut && dEdXCut && ECut && SwissCrossCut)		EcalSpike.push_back(cands); 
//	if(QualCut && dEdXCut)	         			Trk.push_back(cands);
	if(dEdXCut)	         			Trk.push_back(cands);
	if(dEdXCut && SwissCrossCut)		EcalSpike.push_back(cands); 
	if(dEdXCut ) 				EcalShower.push_back(cands); 

     }//end of cand scan	
    if(nPhoton!=0){
    for(unsigned p=0;p<nPhoton;p++){
        Photon &photon = pho[p];
        bool PhotonPtCut = evalPhoton(photon);
        //HighPtPhoton.push_back(photon);
        if(PhotonPtCut) HighPtPhoton.push_back(photon);
    }
    }

	

	sort(Quality.begin(),Quality.begin()+Quality.size());

        sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	
	double mono_pt =0;
	mono_pt = sqrt(pow(CutFlowCand_TRG[0].mono_px_,2)+pow(CutFlowCand_TRG[0].mono_py_,2));
	double amon_pt =0;
	amon_pt = sqrt(pow(CutFlowCand_TRG[0].amon_px_,2)+pow(CutFlowCand_TRG[0].amon_py_,2));
	
	PlotSet &x = Plot[0];
	cout<<"ev "<<ev<<endl;

	
			
	sort(EcalSpike.begin(),EcalSpike.begin()+EcalSpike.size());
	if(EcalSpike.size()>0)
	{
	    
	    for(int i=0; i<EcalSpike.size();i++){
		double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(EcalSpike[i].eta_-EcalSpike[0].mono_eta_,2)+
                                pow(EcalSpike[i].phi_-EcalSpike[0].mono_phi_,2));
                am_deltaR= sqrt(pow(EcalSpike[i].eta_-EcalSpike[0].amon_eta_,2)+
                                pow(EcalSpike[i].phi_-EcalSpike[0].amon_phi_,2));
                if(m_deltaR<0.15) 
                {
			EcalSpikeNum++;
			cout<<"Spikelike m in Ecal "<<EcalSpikeNum<<endl;
                        x.GetPlot(MonoEtavsSpike)->Fill(mono_eta_);
			
		}
		else if(am_deltaR<0.15)
		{
			EcalSpikeNum++;
                        x.GetPlot(AmonEtavsSpike)->Fill(amon_eta_);
			cout<<"Spikelike am in Ecal "<<EcalSpikeNum<<endl;
		}
		x.GetPlot(MonoPtvsEcal)->Fill(mono_pt);
		x.GetPlot(AmonPtvsEcal)->Fill(amon_pt);
	      //}
	    }
	}


	//Ecal Shower 
	sort(EcalShower.begin(),EcalShower.begin()+EcalShower.size());
	if(EcalShower.size()>0)
	{
	    PlotSet &x = Plot[0];
	    for(int i=0; i<EcalShower.size();i++){
	
		if(EcalShower[i].SwissCross_<0.95){
	        double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(EcalShower[i].eta_-EcalShower[0].mono_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].mono_phi_,2));
                am_deltaR= sqrt(pow(EcalShower[i].eta_-EcalShower[0].amon_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].amon_phi_,2));
                if(m_deltaR<0.15) 
                {
			EcalShowerNum++;
			cout<<"Showerlike m in Ecal "<<EcalShowerNum<<endl;
                        x.GetPlot(MonoEtavsShower)->Fill(mono_eta_);
			
		}
		else if(am_deltaR<0.15)
		{
			EcalShowerNum++;
                        x.GetPlot(AmonEtavsShower)->Fill(amon_eta_);
			cout<<"Showerlike am in Ecal "<<EcalShowerNum<<endl;
		}
		
		x.GetPlot(MonoPtvsEcal)->Fill(mono_pt);
		x.GetPlot(AmonPtvsEcal)->Fill(amon_pt);
	    	}		
	   }
	}

	
	sort(Trk.begin(),Trk.begin()+Trk.size());
	if(Trk.size()>0)
	{
	    
	    PlotSet &x = Plot[0];
	    for(int i=0; i<Trk.size();i++){
           	
		
		if(Trk[i].e55_<0){
		
		double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(Trk[i].eta_-Trk[0].mono_eta_,2)+
                                pow(Trk[i].phi_-Trk[0].mono_phi_,2));
                am_deltaR= sqrt(pow(Trk[i].eta_-Trk[0].amon_eta_,2)+
                                pow(Trk[i].phi_-Trk[0].amon_phi_,2));
		
                if(m_deltaR<0.15) 
                {
			TrackerMono++;
			cout<<" m in Tracker "<<TrackerMono<<endl;
		}
		else if(am_deltaR<0.15)
		{
			TrackerMono++;
			cout<<"am in Tracker "<<TrackerMono<<endl;
		}
		x.GetPlot(MonoPtvsTracker)->Fill(mono_pt);
		x.GetPlot(AmonPtvsTracker)->Fill(amon_pt);
	    }
	  }
	}
	
		//2
		if( EcalShowerNum==0 &&  TrackerMono==0 && EcalSpikeNum==2 ){ 
			x.GetPlot(PFMETvsEcalTwoSpike)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsEcalTwoSpike)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsEcalTwoSpike)->Fill(CaloMET_pt);
                        cout<<"PFMET_pt   "<<PFMET_pt<<endl;
                        cout<<"GenMET_pt  "<<GenMET_pt<<endl;
                        cout<<"CaloMET_pt "<<CaloMET_pt<<endl;
			TwoSpikeLike++;
		}
		else if(EcalShowerNum == 2 && TrackerMono==0 && EcalSpikeNum ==0){
			x.GetPlot(PFMETvsEcalTwoShower)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsEcalTwoShower)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsEcalTwoShower)->Fill(CaloMET_pt);
			cout<<"PFMET   "<<PFMET_pt<<endl;
			cout<<"GenMET  "<<GenMET_pt<<endl;
			cout<<"CaloMET "<<CaloMET_pt<<endl;
			TwoShowerLike++;
		}
		else if(EcalShowerNum == 1 && TrackerMono ==0 && EcalSpikeNum == 1 ){
			x.GetPlot(PFMETvsEcalOneShowerOneSpike)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsEcalOneShowerOneSpike)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsEcalOneShowerOneSpike)->Fill(CaloMET_pt);
			cout<<"PFMET   "<<PFMET_pt<<endl;
			cout<<"GenMET  "<<GenMET_pt<<endl;
			cout<<"CaloMET "<<CaloMET_pt<<endl;
			OneShowerOneSpike++;
		}
		else if(EcalShowerNum==0 && TrackerMono==2 && EcalSpikeNum==0){
                        x.GetPlot(PFMETvsTrkTwo)->Fill(PFMET_pt);
                        x.GetPlot(GenMETvsTrkTwo)->Fill(GenMET_pt);
                //        x.GetPlot(CaloMETvsTrkTwo)->Fill(CaloMET_pt);
                        cout<<"PFMET    "<<PFMET_pt<<endl;
                        cout<<"GenMET   "<<GenMET_pt<<endl;
                        cout<<"CaloMET  "<<CaloMET_pt<<endl;
			TwoTrk++;
		}
		else if(EcalShowerNum ==1 && TrackerMono==1 && EcalSpikeNum==0){
			x.GetPlot(PFMETvsEcalOneShower_TrkOne)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsEcalOneShower_TrkOne)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsEcalOneShower_TrkOne)->Fill(CaloMET_pt);
                        cout<<"PFMET   "<<PFMET_pt<<endl;
                        cout<<"GenMET  "<<GenMET_pt<<endl;
                        cout<<"CaloMET "<<CaloMET_pt<<endl;
			OneShowerOneTrk++;
		}
		else if(EcalShowerNum == 0 && TrackerMono==1 && EcalSpikeNum ==1) {
			x.GetPlot(PFMETvsEcalOneSpike_TrkOne)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsEcalOneSpike_TrkOne)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsEcalOneSpike_TrkOne)->Fill(CaloMET_pt);
                        cout<<"PFMET    "<<PFMET_pt<<endl;
                        cout<<"GenMET   "<<GenMET_pt<<endl;
                        cout<<"CaloMET  "<<CaloMET_pt<<endl;
			OneSpikeOneTrk++;
		}

		//1
		else if(EcalShowerNum==0 && TrackerMono==0 && EcalSpikeNum ==1 ){ 
			x.GetPlot(PFMETvsEcalOneSpike)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsEcalOneSpike)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsEcalOneSpike)->Fill(CaloMET_pt);
                        cout<<"PFMET     "<<PFMET_pt<<endl;
			cout<<"GenMET    "<<GenMET_pt<<endl;
                        cout<<"CaloMET     "<<CaloMET_pt<<endl;
			OneSpikeLike++;
		}
		else if(EcalShowerNum==1 && TrackerMono==0 && EcalSpikeNum ==0){ 
			x.GetPlot(PFMETvsEcalOneShower)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsEcalOneShower)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsEcalOneShower)->Fill(CaloMET_pt);
                        cout<<"PFMET     "<<PFMET_pt<<endl;
                        cout<<"GenMET    "<<GenMET_pt<<endl;
                        cout<<"CaloMET   "<<CaloMET_pt<<endl;
			 OneShowerLike++;
		}
		else if( EcalShowerNum == 0 && TrackerMono==1 && EcalSpikeNum ==0 ){
			x.GetPlot(PFMETvsTrkOne)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsTrkOne)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsTrkOne)->Fill(CaloMET_pt);
                        cout<<"PFMET   "<<PFMET_pt<<endl;
                        cout<<"GenMET  "<<GenMET_pt<<endl;
                        cout<<"CaloMET "<<CaloMET_pt<<endl;
			OneTrk++;
		}
		else{ //if(TrackerMono==0 && EcalMono ==0){
			x.GetPlot(PFMETvsZeroMono)->Fill(PFMET_pt);
			x.GetPlot(GenMETvsZeroMono)->Fill(GenMET_pt);
		//	x.GetPlot(CaloMETvsZeroMono)->Fill(CaloMET_pt);
			cout<<"PFMET   "<<PFMET_pt<<endl;
			cout<<"GenMET  "<<GenMET_pt<<endl;
			cout<<"CaloMET "<<CaloMET_pt<<endl;
			ZeroMono++;
		}

	
		bool RealEvt=0;	
	cout<<"Quality size "<<Quality.size()<<endl;
	for(int i=0;i<Quality.size();i++){
                double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(Quality[i].eta_-Quality[0].mono_eta_,2)+
                                pow(Quality[i].phi_-Quality[0].mono_phi_,2));
                am_deltaR= sqrt(pow(Quality[i].eta_-Quality[0].amon_eta_,2)+
                                pow(Quality[i].phi_-Quality[0].amon_phi_,2));
                if(m_deltaR<0.15||am_deltaR<0.15)
                {
			cout<<i<<endl;
			cout<<"     eta "<<setprecision(5)<<Quality[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<Quality[i].phi_<<endl;
			cout<<"      Et "<<setprecision(5)<<Quality[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<Quality[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<Quality[i].f51_<<endl;
			cout<<"     f15 "<<setprecision(5)<<Quality[i].f15_<<endl;
			cout<<"   Cross "<<setprecision(5)<<Quality[i].SwissCross_<<endl;
			cout<<"     MET "<<setprecision(5)<<Quality[i].MET_pt_<<endl;
			RealNum++;
			RealEvt=1;
		
              		cout<<"          candidate           monoGen         antimonoGen"<<endl;
                	cout<<"eta      "<<setprecision(5)<<Quality[i].eta_<<setw(20)<<Quality[i].mono_eta_
                        	         <<setw(20)<<Quality[i].amon_eta_<<endl;
                	cout<<"phi      "<<setprecision(5)<<Quality[i].phi_<<setw(20)<<Quality[i].mono_phi_
                       		         <<setw(20)<<Quality[i].amon_phi_<<endl;
			cout<<"pT	"<<setw(27)<<mono_pt<<setw(20)<<amon_pt<<endl;
                	cout<<"m deltaR "<<m_deltaR<<endl;
                	cout<<"a deltaR "<<am_deltaR<<endl;
			cout<<"--------------------------------------------"<<endl;
		}
	}

	/*		double GenMETMinusTwoM = 0;		
			GenMETMinusTwoM = sqrt(pow(CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].mono_px_-CutFlowCand_TRG[0].amon_px_,2)+
					       pow(CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].mono_py_-CutFlowCand_TRG[0].amon_py_,2));
	
			cout<<"	      px	py      "<<endl;
			cout<<"GenMET      "<<CutFlowCand_TRG[0].GenMET_px_<<setw(10)<<CutFlowCand_TRG[0].GenMET_py_<<endl;
			cout<<"GenMono     "<<CutFlowCand_TRG[0].mono_px_<<setw(10)<<CutFlowCand_TRG[0].mono_py_<<endl;
			cout<<"GenAmon     "<<CutFlowCand_TRG[0].amon_px_<<setw(10)<<CutFlowCand_TRG[0].amon_py_<<endl;
			cout<<"GenMET+Mono "<<CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].mono_px_<<setw(10)
					    <<CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].mono_py_<<endl;
			cout<<"GenMET+Amon "<<CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].amon_px_<<setw(10)
					    <<CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].amon_py_<<endl;
			cout<<"GenMET+two  "<<CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].mono_px_-CutFlowCand_TRG[0].amon_px_<<setw(10)
					    <<CutFlowCand_TRG[0].GenMET_px_-CutFlowCand_TRG[0].mono_py_-CutFlowCand_TRG[0].amon_py_<<endl;

			x.GetPlot(GenMET)->Fill(GenMET_pt);
			x.GetPlot(GenMETminusTwoMonopole)->Fill(GenMETMinusTwoM);
	*/		
		cout<<"Print out Photon information----------"<<endl;        
                for(int j=0;j<HighPtPhoton.size();j++){
                        cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                        cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                        cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
                }

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


        	cout<<"event tag"<<CutFlowCand_TRG[0].event_<<endl;
		if(RealEvt==1) RealEvNum++;//If there is any monopole in event(where in tracker or ecal or number), the event is taken
		RealEvt=0;
			cout<<"=================================================="<<endl;		
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
  
  void SignalEff(const string trName,unsigned TotalEvents)
  {
        //signal efficiency = no. of events after all selection cuts/all events
	
	cout<<endl;
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

	cout<<"---------------2---------------------------"<<endl;
	cout<<"2 spike -like              "<<TwoSpikeLike<<endl;
	cout<<"2 shower-like              "<<TwoShowerLike<<endl;
	cout<<"1 shower-like 1 spike-like "<<OneShowerOneSpike<<endl;
	cout<<"1 shower-like 1 in Tracker "<<OneShowerOneTrk<<endl;
	cout<<"1 spike -like 1 in Tracker "<<OneSpikeOneTrk<<endl;
	cout<<"2 in Tracker               "<<TwoTrk<<endl;
	cout<<"---------------1---------------------------"<<endl;
	cout<<"1 Shower-like only         "<<OneShowerLike<<endl;
	cout<<"1 spike -like only         "<<OneSpikeLike<<endl;
	cout<<"1 in Tracker  only         "<<OneTrk<<endl;
	cout<<"---------------0---------------------------"<<endl;
	cout<<"0 monopole in the event    "<<ZeroMono<<endl;		
	
	
	cout<<"----------------------------------"<<endl;
	cout<<"NO monopole in this ev "<<ZeroMono<<endl;
	cout<<"Number of monopole ev  "<<RealEvNum<<endl;
        cout<<"Real monopole          "<<RealNum<<endl;
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
  static const double f15Cut_ ;
  static const double f51Cut_ ;
  static const double SwissCrossCut_ ;
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
  bool evalSwissCross(MonoCandidate &mc) { return mc.SwissCross_ > SwissCrossCut_ ; }
  bool evalF15(MonoCandidate &mc) { return mc.f15_ > f15Cut_ ; }
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

  int count_test=0;
  int TRGminitest=0;
  int Qualminitest=0;


//  int TwoMonoInEcal=0;
//  int TwoMonoInTracker=0;
//  int OneTrkOneEcal=0;
//  int OneMonoInEcal=0;
//  int OneMonoInTracker=0;
//  int ZeroMono=0;

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
  const double MonoCuts::e55Cut_ = 0;
//   const double MonoCuts::e55Cut_ = 200;
//  const double MonoCuts::f51Cut_ = 0.6;
  const double MonoCuts::f15Cut_ = 0.95;
  const double MonoCuts::f51Cut_ = 0.95;
  const double MonoCuts::SwissCrossCut_ = 0.95;
  const double MonoCuts::photonCut_ = 200;
void MonoAnalyzerMETvsMono()
{
	TFile *oFile = new TFile("Monopole_MET_analysis_TESTforSwissCross.root","recreate");
	TFile *fin = new TFile("MonoNtuple2018_MC_1000_TestSwisscross2.root");
        TTree *tree = (TTree*)fin->Get("monopoles");
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
        vector<double> * SwissCross =0;
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
        tree->SetBranchAddress("cand_SwissCross",&SwissCross);
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
	MonoCuts TrgAnalysis_Photon                                 ("HLT_Photon200",oFile);
/*	MonoCuts TrgAnalysis_PFMET300		                    ("HLT_PFMET300",oFile);
	MonoCuts TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned          ("HLT_PFMET170_HBHE_BeamHaloCleaned",oFile);
	MonoCuts TrgAnalysis_PFMET140_PFMHT140_IDTight              ("HLT_PFMET140_PFMHT140_IDTight",oFile);
	MonoCuts TrgAnalysis_PFMET250_HBHECleaned                   ("HLT_PFMET250_HBHECleaned",oFile);
	MonoCuts TrgAnalysis_PFMET300_HBHECleaned                   ("HLT_PFMET300_HBHECleaned",oFile);
	MonoCuts TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned          ("HLT_PFMET200_HBHE_BeamHaloCleaned",oFile);
	MonoCuts TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned   ("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",oFile);
	MonoCuts TrgAnalysis_CaloMET300_HBHECleaned                 ("HLT_CaloMET300_HBHECleaned",oFile);
*/
	vector<MonoCandidate> cand(10);	
	vector<Photon> photon(0);
	
        for(unsigned ev=0; ev<NEvents;ev++){
//		cout<<"Enter event loop"<<endl;
//		cout<<"event number "<<NEvents<<endl;
             tree->GetEntry(ev);
		
		if(nCandidates>Et.size()) Et.resize(nCandidates);
		if(nCandidates>cand.size()) cand.resize(nCandidates);

		if(nPhoton>photon.size()) photon.resize(nPhoton);
//		cout<<"nPhoton "<<nPhoton<<endl;	
//		cout<<"photon size "<<photon.size()<<endl;

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
		Et[i],
		(*SwissCross)[i],
//	        (*e55)[i],
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
			noTrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,true,ev,MET_pt,GenMET_pt,CaloMET_pt);
		cout<<"pass noTRG analysis"<<endl;		

/*         TrgAnalysis_PFMET300                            .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300,ev);
         TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET170_HBHE_BeamHaloCleaned,ev);
         TrgAnalysis_PFMET140_PFMHT140_IDTight           .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET140_PFMHT140_IDTight,ev);
         TrgAnalysis_PFMET250_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET250_HBHECleaned,ev);
         TrgAnalysis_PFMET300_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300_HBHECleaned,ev);
         TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET200_HBHE_BeamHaloCleaned,ev);
         TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned,ev);
         TrgAnalysis_CaloMET300_HBHECleaned              .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_CaloMET300_HBHECleaned,ev);
*/      //   TrgAnalysis_Photon.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon200,ev,MET_pt,GenMET_pt,CaloMET_pt);
		cout<<"pass TRG analysis"<<endl;		
	}//for every event loop
	
			TrgAnalysis_Photon         			  .SignalEff("HLT_Photon200",NEvents); 
//		        TrgAnalysis_PFMET300		                  .SignalEff("HLT_PFMET300");
//		        TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET170_HBHE_BeamHaloCleaned");
//		        TrgAnalysis_PFMET140_PFMHT140_IDTight             .SignalEff("HLT_PFMET140_PFMHT140_IDTight");
//		        TrgAnalysis_PFMET250_HBHECleaned                  .SignalEff("HLT_PFMET250_HBHECleaned");
//		        TrgAnalysis_PFMET300_HBHECleaned                  .SignalEff("HLT_PFMET300_HBHECleaned");
//		        TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET200_HBHE_BeamHaloCleaned");
//		        TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned  .SignalEff("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned");
//		        TrgAnalysis_CaloMET300_HBHECleaned                .SignalEff("HLT_CaloMET300_HBHECleaned");                  


         TrgAnalysis_Photon                              .WritePlots(oFile);
//         TrgAnalysis_PFMET300		                 .WritePlots(oFile);       
//         TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .WritePlots(oFile);     
//         TrgAnalysis_PFMET140_PFMHT140_IDTight           .WritePlots(oFile);
//         TrgAnalysis_PFMET250_HBHECleaned                .WritePlots(oFile);   
//         TrgAnalysis_PFMET300_HBHECleaned                .WritePlots(oFile);
//         TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .WritePlots(oFile);   
//         TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.WritePlots(oFile);
//         TrgAnalysis_CaloMET300_HBHECleaned              .WritePlots(oFile);   
	                                                     

	noTrgAnalysis.WritePlots(oFile);
	noTrgAnalysis.SignalEff("NOTRG",NEvents);
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
