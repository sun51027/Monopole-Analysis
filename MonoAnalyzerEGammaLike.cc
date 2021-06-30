//////////////////////////////////
//
//	File name: MonoAnalyzerEgammalikevsMET.cc
//	author: Shih Lin
//	Content: Matching
//		 Reco Photon information(21/4/24 update) 
//		 eta vs egamma like monopole (done in 21/6/11)
//		 pt  vs egamma like monopole (done in 21/6/11)
//	Purpose: eff of photon-like/electron-like
//		 compare them to the MET
//		 Should be modify to event 	  
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
  
  PhotonlikePt,
  ElectronlikePt,
  EGlikePt,
  PhotonlikeEta,
  ElectronlikeEta,
  EGlikeEta,


PFMETvsTwoPhotonLike,		
PFMETvsTwoElectronLike,		
PFMETvsOnePhotonLike,		
PFMETvsOneElectronLike,		
PFMETvsOnePhotonOneElectronlike,	
PFMETvsTwoNonSpikeLike,		
PFMETvsOneNonSpikeLike,		
PFMETvsOnePhotonOneNonSpikeLike,	
PFMETvsOneElectronOneNonSpikeLike,	
PFMETvsOneNonSpikeOneSpike,		
PFMETvsOneElectronOneSpike,		
PFMETvsOnePhotonOneSpike,		
PFMETvsOneSpike,			
PFMETvsTwoSpike,			
  MonoPtvsTracker,
  AmonPtvsTracker,
  MonoPtvsEcal,
  AmonPtvsEcal

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
class Electron
{
public:
  Electron(){}
  Electron(double ele_eta, double ele_phi, double ele_pt):ele_eta_(ele_eta),ele_phi_(ele_phi),ele_pt_(ele_pt){}
  Electron(const Electron &mc):ele_eta_(mc.ele_eta_),ele_phi_(mc.ele_phi_),ele_pt_(mc.ele_pt_){}
  double ele_eta_;
  double ele_phi_;
  double ele_pt_;

  ~Electron(){}

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
 
	

/*	x.CreatPlot(ElectronlikePt,new TH1D("ElectronlikePt","",100,0,800));	
	x.CreatPlot(PhotonlikePt,new TH1D("PhotonlikePt","",100,0,800));	
	x.CreatPlot(EGlikePt,new TH1D("EGlikePt","",100,0,800));	
	x.CreatPlot(TwoEGlikeMET,new TH1D("TwoEGlikeMET","",100,0,800));	
	x.CreatPlot(OneEGlikeMET,new TH1D("OneEGlikeMET","",100,0,800));	

	x.CreatPlot(ElectronlikeEta,new TH1D("ElectronlikeEta","",100,-4,4));	
	x.CreatPlot(PhotonlikeEta,new TH1D("PhotonlikeEta","",100,-4,4));	
	x.CreatPlot(EGlikeEta,new TH1D("EGlikeEta","",100,-4,4));	
*/	
/*	x.CreatPlot(PFMETvsEcalTwoShower,new TH1D("PFMETvsEcalTwoShower","2 Showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(PFMETvsEcalTwoSpike,new TH1D("PFMETvsEcalTwoSpike","2 Spikelike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneShowerOneSpike,new TH1D("PFMETvsEcalOneShowerOneSpike","1 Spikelike & 1 showerlike Monopole in Ecal ",100,0,800));	
 	x.CreatPlot(PFMETvsTrkTwo,new TH1D("PFMETvsTrkTwo","2 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneShower_TrkOne,new TH1D("PFMETvsEcalOneShower_TrkOne","1 showerlike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneSpike_TrkOne,new TH1D("PFMETvsEcalOneSpike_TrkOne","1 spikelike Monopole in Ecal & 1 in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneShower,new TH1D("PFMETvsEcalOneShower","1 showerlike Monopole in Ecal",100,0,800));	
 	x.CreatPlot(PFMETvsEcalOneSpike,new TH1D("PFMETvsEcalOneSpike","1 spikelike Monopole in ECal",100,0,800));	
 	x.CreatPlot(PFMETvsTrkOne,new TH1D("PFMETvsTrkOne","1 Monopole in Tracker",100,0,800));	
 	x.CreatPlot(PFMETvsZeroMono,new TH1D("PFMETvsZeroMono","No monopole is reconstructed in this event",100,0,800));	
	
//	x.CreatPlot(MonoPtvsTracker,new TH1D("MonoPtvsTracker","",100,0,2500));	
//	x.CreatPlot(AmonPtvsTracker,new TH1D("AmonPtvsTracker","",100,0,2500));	
//	x.CreatPlot(MonoPtvsEcal,new TH1D("MonoPtvsEcal","",100,0,2500));	
//	x.CreatPlot(AmonPtvsEcal,new TH1D("AmonPtvsEcal","",100,0,2500));	
*/

// 	x.CreatPlot(PFMETvsTwoPhotonLike,  new TH1D("PFMETvsTwoPhotonLike","",40,0,800));
// 	x.CreatPlot(PFMETvsOnePhotonLike,  new TH1D("PFMETvsOnePhotonLike","",40,0,800));
// 	x.CreatPlot(PFMETvsTwoElectronLike,new TH1D("PFMETvsTwoElectronLike","",40,0,800));
// 	x.CreatPlot(PFMETvsOneElectronLike,new TH1D("PFMETvsOneElectronLike","",40,0,800));
// 	x.CreatPlot(PFMETvsTwoNonSpikeLike,new TH1D("PFMETvsTwoNonSpikeLik","",40,0,800));
// 	x.CreatPlot(PFMETvsOneNonSpikeLike,new TH1D("PFMETvsOneNonSpikeLike","",40,0,800));
// 	x.CreatPlot(PFMETvsOnePhotonOneElectronlike,new TH1D("PFMETvsOnePhotonOneElectronlike","",40,0,800));
	cout<<"the address of plots x "<<&x<<endl;
        x.CreatPlot(PFMETvsTwoPhotonLike,		new TH1D("PFMETvsTwoPhotonLike","",40,0,800));	
        x.CreatPlot(PFMETvsTwoElectronLike,		new TH1D("PFMETvsTwoElectronLike","",40,0,800));		
        x.CreatPlot(PFMETvsOnePhotonLike,		new TH1D("PFMETvsOnePhotonLike","",40,0,800));		
        x.CreatPlot(PFMETvsOneElectronLike,		new TH1D("PFMETvsOneElectronLike","",40,0,800));		
        x.CreatPlot(PFMETvsOnePhotonOneElectronlike,	new TH1D("PFMETvsOnePhotonOneElectronlike","",40,0,800));	
        x.CreatPlot(PFMETvsTwoNonSpikeLike,		new TH1D("PFMETvsTwoNonSpikeLike","",40,0,800));		
        x.CreatPlot(PFMETvsOneNonSpikeLike,		new TH1D("PFMETvsOneNonSpikeLike","",40,0,800));		
        x.CreatPlot(PFMETvsOnePhotonOneNonSpikeLike,	new TH1D("PFMETvsOnePhotonOneNonSpikeLike","",40,0,800));	
        x.CreatPlot(PFMETvsOneElectronOneNonSpikeLike,	new TH1D("PFMETvsOneElectronOneNonSpikeLike","",40,0,800));	
        x.CreatPlot(PFMETvsOneNonSpikeOneSpike,		new TH1D("PFMETvsOneNonSpikeOneSpike","",40,0,800));		
        x.CreatPlot(PFMETvsOneElectronOneSpike,		new TH1D("PFMETvsOneElectronOneSpike","",40,0,800));		
        x.CreatPlot(PFMETvsOnePhotonOneSpike,		new TH1D("PFMETvsOnePhotonOneSpike","",40,0,800));		
        x.CreatPlot(PFMETvsOneSpike,			new TH1D("PFMETvsOneSpike","",40,0,800));			
        x.CreatPlot(PFMETvsTwoSpike,			new TH1D("PFMETvsTwoSpike","",40,0,800));			
	//Name the cut name
	 cutName_[0] = "Quality_";
	 cutName_[1] = "Energy_";
	 cutName_[2] = "F51_";
	 cutName_[3] = "dEdxSig_";
	 cutName_[4] = "HLT_";

  }

  ~MonoCuts(){}
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, vector<Electron> &ele, unsigned nCandidates,unsigned nPhoton, unsigned nElectron, bool TRG, unsigned ev, double PFMET_pt)
  {
	
	//cout<<"enter doAnalysis successfully"<<endl;	
	EcalShower.clear();
	EcalSpike.clear();
	Trk.clear();
	Quality.clear();
	Shower.clear();	

	CutFlowCand_TRG.clear();
        HighPtPhoton.clear();
        HighPtElectron.clear();

	EgammaLike=0;
	SpikeLike=0;
	//EcalShowerNum=0;
	TrackerMono=0;
	photonLike=0;
	electronLike=0;
	NonSpikeLike=0;

     for(unsigned c=0;c<nCandidates;c++){
	MonoCandidate &cands = cand[c];
	bool QualCut = evalQuality(cands);
	bool ECut = evalE(cands);
	bool F15Cut = evalF15(cands);
	bool F51Cut = evalF51(cands);
	bool SwissCrossCut = evalSwissCross(cands);
	bool dEdXCut = evaldEdX(cands);
 
	//no cut
	if(TRG) CutFlowCand_TRG.push_back(cands);

	//cut all 
	if(TRG && ECut && F51Cut && dEdXCut) CutFlowCand_Dedx.push_back(cands);

	if(QualCut) 			 		Quality.push_back(cands); 
	
	//MET
	if(QualCut && dEdXCut && ECut )	   		EcalShower.push_back(cands);//need to add f51<0.85 in the following code)
	if(QualCut && dEdXCut && ECut && SwissCrossCut)	EcalSpike.push_back(cands); 
	if(dEdXCut)	        			Trk.push_back(cands);
//	if(QualCut && dEdXCut)	         			Trk.push_back(cands);
//	if(QualCut && dEdXCut && ECut && F51Cut)		EcalSpike.push_back(cands); 
//	if(QualCut && dEdXCut && ECut ) 			EcalShower.push_back(cands); 

     }//end of cand scan	
    if(nPhoton!=0){
    for(unsigned p=0;p<nPhoton;p++){
        Photon &photon = pho[p];
        bool PhotonPtCut = evalPhoton(photon);
        HighPtPhoton.push_back(photon);
    }
    }

    if(nElectron!=0){
        for(unsigned p=0;p<nElectron;p++){
                Electron &electron = ele[p];
                HighPtElectron.push_back(electron);
        }
    }
	

	sort(Quality.begin(),Quality.begin()+Quality.size());

        sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	
	
	PlotSet &x = Plot[0];
	cout<<"ev "<<ev<<endl;
	
			
	sort(EcalShower.begin(),EcalShower.begin()+EcalShower.size());
	if(EcalShower.size()>0)
	{
	    cout<<"EcalShower size "<<EcalShower.size()<<endl; 
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
			Shower.push_back(EcalShower[i]);// put matched EcalShower Monopole into "Shower" to match with Egamma
		}
		else if(am_deltaR<0.15)
		{
			EcalShowerNum++;
			Shower.push_back(EcalShower[i]);
			cout<<"Showerlike am in Ecal "<<EcalShowerNum<<endl;
		}
//		x.GetPlot(MonoPtvsEcal)->Fill(EcalShower[0].mono_pt_);
//		x.GetPlot(AmonPtvsEcal)->Fill(amon_pt);
	    }
	    }
	}

	//Replace showe-like with Egamma like
		cout<<"Shower size "<<Shower.size()<<endl;


	//-------Photon-like-------//
        //to see monopole is photon-like or spike-like
        //define: photon-monopole angle<0.15 is photon-like
	cout<<"--------Print out Photon information----------"<<endl;        
	cout<<"Photon size "<<HighPtPhoton.size()<<endl;		 
        for(int j=0;j<HighPtPhoton.size();j++){
                cout<<j<<endl;
		cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
        }
	cout<<endl;
               
	for(int j=0;j<Shower.size();j++){
		if(HighPtPhoton.size()>0){
		for(int i=0;i<HighPtPhoton.size();i++){
                        double photonMono_deltaR = 0;
                        photonMono_deltaR = sqrt(pow(Shower[j].eta_-HighPtPhoton[i].pho_eta_,2)+
                                                 pow(Shower[j].phi_-HighPtPhoton[i].pho_phi_,2));
                        if(photonMono_deltaR < 0.15 ){
                                cout<<"deltaR "<<photonMono_deltaR<<"  photon-like monopole"<<endl;
				photonLike++;// otherwise not photon-like monopole
				EgammaLike++;
		//		x.GetPlot(PhotonlikePt)->Fill(HighPtPhoton[i].pho_pt_);
		//		x.GetPlot(PhotonlikeEta)->Fill(HighPtPhoton[i].pho_eta_);
		//		x.GetPlot(EGlikePt)->Fill(HighPtPhoton[i].pho_pt_);
		//		x.GetPlot(EGlikeEta)->Fill(HighPtPhoton[i].pho_eta_);
			}
			else{
		//		NonSpikeLike++;
                                cout<<"deltaR "<<photonMono_deltaR<<"  No photon-like"<<endl;	
			}	
                }
		}
	}
	//-------Electron-like-------//
        //to see monopole is electron-like or spike-like
        //define: electron-monopole angle<0.15 is electron-like
//	cout<<"-------Print out electron information----------"<<endl;        
	cout<<"Electron size "<<HighPtElectron.size()<<endl;		 
        for(int i=0;i<HighPtElectron.size();i++){
                cout<<i<<endl;
		cout<<"ele  eta "<<HighPtElectron[i].ele_eta_<<endl;
                cout<<"ele  phi "<<HighPtElectron[i].ele_phi_<<endl;
                cout<<"ele   pt "<<HighPtElectron[i].ele_pt_<<endl;
        }
	cout<<endl;
	for(int j=0;j<Shower.size();j++){
                if(HighPtElectron.size()!=0){ 
       		for(int i=0;i<HighPtElectron.size();i++){
        	        double electronMono_deltaR = 0;
        	        electronMono_deltaR = sqrt(pow(Shower[j].eta_-HighPtElectron[i].ele_eta_,2)+
        	                                   pow(Shower[j].phi_-HighPtElectron[i].ele_phi_,2));
        	        if(electronMono_deltaR < 0.15 )
			{
				electronLike++;// otherwise not eleton-like monopole
        	                cout<<"deltaR "<<electronMono_deltaR<<"  electron-like monopole"<<endl;
				EgammaLike++;
				//x.GetPlot(ElectronlikePt)->Fill(HighPtElectron[i].ele_pt_);
				//x.GetPlot(ElectronlikeEta)->Fill(HighPtElectron[i].ele_eta_);
				//x.GetPlot(EGlikePt)->Fill(HighPtElectron[i].ele_pt_);
				//x.GetPlot(EGlikeEta)->Fill(HighPtElectron[i].ele_eta_);
			}
			else{
			//	NonSpikeLike++;
        	                cout<<"deltaR "<<electronMono_deltaR<<"  No electron-like"<<endl;
			}
		}
		}
	}

	//Ecal Spike 
	sort(EcalSpike.begin(),EcalSpike.begin()+EcalSpike.size());
	if(EcalSpike.size()>0)	{
	    PlotSet &x = Plot[0];
	    for(int i=0; i<EcalSpike.size();i++){
	
	        double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(EcalSpike[i].eta_-EcalSpike[0].mono_eta_,2)+
                                pow(EcalSpike[i].phi_-EcalSpike[0].mono_phi_,2));
                am_deltaR= sqrt(pow(EcalSpike[i].eta_-EcalSpike[0].amon_eta_,2)+
                                pow(EcalSpike[i].phi_-EcalSpike[0].amon_phi_,2));
                if(m_deltaR<0.15) 
                {
			SpikeLike++;
			cout<<"Spikelike m in Ecal "<<SpikeLike<<endl;
			
		}
		else if(am_deltaR<0.15)
		{
			SpikeLike++;
			cout<<"Spikelike am in Ecal "<<SpikeLike<<endl;
		}
		
			//x.GetPlot(MonoPtvsEcal)->Fill(mono_pt);
			//x.GetPlot(AmonPtvsEcal)->Fill(amon_pt);
	    	}		
	}
	if(electronLike==0 && photonLike==0 && Shower.size()==2)	NonSpikeLike=2;
        else if(electronLike==0 && photonLike==0 && Shower.size()==1)        NonSpikeLike=1;
        else if(electronLike==1 && photonLike==0 && Shower.size()==1)        NonSpikeLike=1;
        else if(electronLike==0 && photonLike==1 && Shower.size()==1)        NonSpikeLike=1;
        else if(electronLike==1 && photonLike==1 && Shower.size()==1)        electronLike=0;

	
	cout<<"NonSpikeLike "<<NonSpikeLike<<endl;
	cout<<"electronLike "<<electronLike<<endl;
	cout<<"photonLike   "<<photonLike<<endl;
	cout<<"SpikeLike    "<<SpikeLike<<endl;
	//calculate the efficiency of egamma-like event
	if(     NonSpikeLike ==0 && electronLike ==0 && photonLike == 2 && SpikeLike == 0 ){
		 TwoPhotonLikeEvent++;
		x.GetPlot(PFMETvsTwoPhotonLike)->Fill(PFMET_pt);
	}
	else if(NonSpikeLike ==0 && electronLike ==0 && photonLike == 1 && SpikeLike == 0){
		 OnePhotonLikeEvent++;
		x.GetPlot(PFMETvsOnePhotonLike)->Fill(PFMET_pt);
	}
	else if(NonSpikeLike ==0 && electronLike ==2 && photonLike == 0 && SpikeLike == 0){
		 TwoElectronLikeEvent++;
		x.GetPlot(PFMETvsTwoElectronLike)->Fill(PFMET_pt);
	}
	else if(NonSpikeLike ==0 && electronLike ==1 && photonLike == 0 && SpikeLike == 0){
		 OneElectronLikeEvent++;
		x.GetPlot(PFMETvsOneElectronLike)->Fill(PFMET_pt);
	}
	else if(NonSpikeLike ==0 && electronLike ==1 && photonLike == 1 && SpikeLike == 0){
		OnePhotonOneElectronlikeEvent++;	
		x.GetPlot(PFMETvsOnePhotonOneElectronlike)->Fill(PFMET_pt);
	}
	//Also need to calculate the eff which shower-like/nonspikelike monopole are not identified as egamma.
	else if(NonSpikeLike ==2 && electronLike ==0 && photonLike == 0 && SpikeLike == 0){
		TwoNonSpikeLikeEvent++;
		x.GetPlot(PFMETvsTwoNonSpikeLike)->Fill(PFMET_pt);
	}
	else if(NonSpikeLike ==1 && electronLike ==0 && photonLike == 0 && SpikeLike == 0){
		OneNonSpikeLikeEvent++;
		x.GetPlot(PFMETvsOneNonSpikeLike)->Fill(PFMET_pt);
	}
        else if(NonSpikeLike ==1 && electronLike ==1 && photonLike == 0 && SpikeLike == 0){
                OneElectronOneNonSpikeLikeEvent++;
		x.GetPlot(PFMETvsOneElectronOneNonSpikeLike)->Fill(PFMET_pt);
        }
        else if(NonSpikeLike ==1 && electronLike ==0 && photonLike == 1 && SpikeLike == 0){
                OnePhotonOneNonSpikeLikeEvent++;
		x.GetPlot(PFMETvsOnePhotonOneNonSpikeLike)->Fill(PFMET_pt);
        }
        else if(NonSpikeLike ==1 && electronLike ==0 && photonLike == 0 && SpikeLike == 1){
		OneNonSpikeOneSpikeEvent++;
		x.GetPlot(PFMETvsOneNonSpikeOneSpike)->Fill(PFMET_pt);
	}
	// Need to add the category of spike vs showe like monopole with PFMET
        else if(NonSpikeLike ==0 && electronLike ==0 && photonLike == 1 && SpikeLike == 1){
		OnePhotonOneSpikeEvent++;
		x.GetPlot(PFMETvsOnePhotonOneSpike)->Fill(PFMET_pt);
	}
        else if(NonSpikeLike ==0 && electronLike ==1 && photonLike == 0 && SpikeLike == 1){
		OneElectronOneSpikeEvent++;
		x.GetPlot(PFMETvsOneElectronOneSpike)->Fill(PFMET_pt);
	}
        else if(NonSpikeLike ==0 && electronLike ==0 && photonLike == 0 && SpikeLike == 1){
		OneSpikeEvent++;
		x.GetPlot(PFMETvsOneSpike)->Fill(PFMET_pt);
	}
        else if(NonSpikeLike ==0 && electronLike ==0 && photonLike == 0 && SpikeLike == 2){
		TwoSpikeEvent++;
		x.GetPlot(PFMETvsTwoSpike)->Fill(PFMET_pt);
	}


	
/*	sort(Trk.begin(),Trk.begin()+Trk.size());
	if(Trk.size()>0)
	{
	    
	    PlotSet &x = Plot[0];
	    for(int i=0; i<Trk.size();i++){
           	
		
		if(Trk[i].e55_<=0){
		
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
	

		if(p==2||e==2||(p==1 && e==1)){
			x.GetPlot(TwoEGlikeMET)->Fill(PFMET_pt);
		}	
		else if((p==1 && e == 0)||(e==1 && p==0)){
			x.GetPlot(OneEGlikeMET)->Fill(PFMET_pt);
		}
		
			
*/
		//2
/*		if( EgammaLike==0 &&  TrackerMono==0 && EcalSpikeNum==2 ){ 
			x.GetPlot(PFMETvsEcalTwoSpike)->Fill(PFMET_pt);
                        cout<<"PFMET_pt   "<<PFMET_pt<<endl;
			TwoSpikeLike++;
		}
		else if(EgammaLike == 2 && TrackerMono==0 && EcalSpikeNum ==0){
			x.GetPlot(PFMETvsEcalTwoShower)->Fill(PFMET_pt);
			cout<<"PFMET   "<<PFMET_pt<<endl;
			TwoShowerLike++;
		}
		else if(EgammaLike == 1 && TrackerMono ==0 && EcalSpikeNum == 1 ){
			x.GetPlot(PFMETvsEcalOneShowerOneSpike)->Fill(PFMET_pt);
			cout<<"PFMET   "<<PFMET_pt<<endl;
			OneShowerOneSpike++;
		}
		else if(EgammaLike==0 && TrackerMono==2 && EcalSpikeNum==0){
                        x.GetPlot(PFMETvsTrkTwo)->Fill(PFMET_pt);
                        cout<<"PFMET    "<<PFMET_pt<<endl;
			TwoTrk++;
		}
		else if(EgammaLike ==1 && TrackerMono==1 && EcalSpikeNum==0){
			x.GetPlot(PFMETvsEcalOneShower_TrkOne)->Fill(PFMET_pt);
                        cout<<"PFMET   "<<PFMET_pt<<endl;
			OneShowerOneTrk++;
		}
		else if(EgammaLike == 0 && TrackerMono==1 && EcalSpikeNum ==1) {
			x.GetPlot(PFMETvsEcalOneSpike_TrkOne)->Fill(PFMET_pt);
                        cout<<"PFMET    "<<PFMET_pt<<endl;
			OneSpikeOneTrk++;
		}

		//1
		else if(EgammaLike==0 && TrackerMono==0 && EcalSpikeNum ==1 ){ 
			x.GetPlot(PFMETvsEcalOneSpike)->Fill(PFMET_pt);
                        cout<<"PFMET     "<<PFMET_pt<<endl;
			OneSpikeLike++;
		}
		else if(EgammaLike==1 && TrackerMono==0 && EcalSpikeNum ==0){ 
			x.GetPlot(PFMETvsEcalOneShower)->Fill(PFMET_pt);
                        cout<<"PFMET     "<<PFMET_pt<<endl;
			 OneShowerLike++;
		}
		else if( EgammaLike == 0 && TrackerMono==1 && EcalSpikeNum ==0 ){
			x.GetPlot(PFMETvsTrkOne)->Fill(PFMET_pt);
                        cout<<"PFMET   "<<PFMET_pt<<endl;
			OneTrk++;
		}
		else{ 
			x.GetPlot(PFMETvsZeroMono)->Fill(PFMET_pt);
			cout<<"PFMET   "<<PFMET_pt<<endl;
			ZeroMono++;
		}
*/
	
		bool RealEvt=0;	
	cout<<"---------------------------------------------"<<endl;
	cout<<"EcalShower size "<<EcalShower.size()<<endl;
	for(int i=0;i<EcalShower.size();i++){
                double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(EcalShower[i].eta_-EcalShower[0].mono_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].mono_phi_,2));
                am_deltaR= sqrt(pow(EcalShower[i].eta_-EcalShower[0].amon_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].amon_phi_,2));
        //        if(m_deltaR<0.15||am_deltaR<0.15)
          //      {
			cout<<i<<endl;
			cout<<"     eta "<<setprecision(5)<<EcalShower[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<EcalShower[i].phi_<<endl;
			cout<<"      Et "<<setprecision(5)<<EcalShower[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<EcalShower[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<EcalShower[i].f51_<<endl;
			cout<<"     f15 "<<setprecision(5)<<EcalShower[i].f15_<<endl;
			cout<<"   Cross "<<setprecision(5)<<EcalShower[i].SwissCross_<<endl;
			cout<<"     MET "<<setprecision(5)<<EcalShower[i].MET_pt_<<endl;
			RealNum++;
			RealEvt=1;
		
              		cout<<"          candidate           monoGen         antimonoGen"<<endl;
                	cout<<"eta      "<<setprecision(5)<<EcalShower[i].eta_<<setw(20)<<EcalShower[i].mono_eta_
                        	         <<setw(20)<<EcalShower[i].amon_eta_<<endl;
                	cout<<"phi      "<<setprecision(5)<<EcalShower[i].phi_<<setw(20)<<EcalShower[i].mono_phi_
                       		         <<setw(20)<<EcalShower[i].amon_phi_<<endl;
			//cout<<"pT	"<<setw(27)<<mono_pt<<setw(20)<<amon_pt<<endl;
                	cout<<"m deltaR "<<m_deltaR<<endl;
                	cout<<"a deltaR "<<am_deltaR<<endl;
			cout<<"--------------------------------------------"<<endl;
	//	}
	}


        	cout<<"event tag"<<CutFlowCand_TRG[0].event_<<endl;
		if(RealEvt==1) RealEvNum++;//If there is any monopole in event(where in tracker or ecal or number), the event is taken
		RealEvt=0;
			cout<<"=================================================="<<endl;		
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
/*        cout<<"        TRG "<<count<<endl;
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
	*/
	cout<<"# of EcalShower : pass dEdx>9, E55 > 0, NO Quality, f51..."<<endl;
	cout<<"EcalShowerNum "<<EcalShowerNum<<endl;
        cout<<"photonLike "<<photonLike<<endl;
        cout<<"electronLike "<<electronLike<<endl;
        cout<<"ratio of photonlike "<<(double)photonLike/(double)EcalShowerNum<<endl;
        cout<<"ratio of electronlike "<<(double)electronLike/(double)EcalShowerNum<<endl;
	cout<<"--------------------------------------------"<<endl;
        cout<<"TwoPhotonLikeEvent		"<<(double)TwoPhotonLikeEvent			<<endl;		
        cout<<"TwoElectronLikeEven		"<<(double)TwoElectronLikeEvent			<<endl;	
        cout<<"OnePhotonLikeEvent		"<<(double)OnePhotonLikeEvent			<<endl;	
        cout<<"OneElectronLikeEvent		"<<(double)OneElectronLikeEvent			<<endl;	
        cout<<"OnePhotonOneElectronlikeEvent	"<<(double)OnePhotonOneElectronlikeEvent	<<endl;	
        cout<<"TwoNonSpikeLikeEvent		"<<(double)TwoNonSpikeLikeEvent			<<endl;	
        cout<<"OneNonSpikeLikeEvent		"<<(double)OneNonSpikeLikeEvent			<<endl;	
        cout<<"OnePhotonOneNonSpikeLikeEvent	"<<(double)OnePhotonOneNonSpikeLikeEvent	<<endl;	
        cout<<"OneElectronOneNonSpikeLikeEvent	"<<(double)OneElectronOneNonSpikeLikeEvent	<<endl;
        cout<<"OneNonSpikeOneSpikeEvent		"<<(double)OneNonSpikeOneSpikeEvent		<<endl;	
        cout<<"OneElectronOneSpikeEvent		"<<(double)OneElectronOneSpikeEvent		<<endl;	
        cout<<"OnePhotonOneSpikeEvent		"<<(double)OnePhotonOneSpikeEvent		<<endl;	
        cout<<"OneSpikeEvent			"<<(double)OneSpikeEvent			<<endl;	
        cout<<"TwoSpikeEvent			"<<(double)TwoSpikeEvent			<<endl;	
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
  static const double SwissCrossCut_ ;
  static const double f51Cut_ ;
  static const double photonCut_ ;
  static const double electronCut_ ;
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
  bool evalF15(MonoCandidate &mc) { return mc.f15_ > f15Cut_ ; }
  bool evalSwissCross(MonoCandidate &mc) { return mc.SwissCross_ > SwissCrossCut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
  bool evalElectron(Electron &mc){ return mc.ele_pt_ > electronCut_; }
 
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
  vector<Electron> HighPtElectron;

  vector<MonoCandidate> Shower;
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
  int EgammaLike=0;
  int electronLike=0;
  int NonSpikeLike=0;
  int SpikeLike=0;

	int TwoPhotonLikeEvent=0;
	int TwoElectronLikeEvent=0;
	int OnePhotonLikeEvent=0;
	int OneElectronLikeEvent=0;
	int OnePhotonOneElectronlikeEvent=0;
	int TwoNonSpikeLikeEvent=0;
	int OneNonSpikeLikeEvent=0;
	int OnePhotonOneNonSpikeLikeEvent=0;
	int OneElectronOneNonSpikeLikeEvent=0;
	int OneNonSpikeOneSpikeEvent=0;	
	int OneElectronOneSpikeEvent=0;
	int OnePhotonOneSpikeEvent=0;
	int OneSpikeEvent=0;
	int TwoSpikeEvent=0;
//  int TwoMonoInEcal=0;
//  int TwoMonoInTracker=0;
//  int OneTrkOneEcal=0;
//  int OneMonoInEcal=0;
//  int OneMonoInTracker=0;
//  int ZeroMono=0;

  int TrackerMono=0;
  int EcalShowerNum=0;
 
/*  int TwoSpikeLike=0;
  int TwoShowerLike=0;
  int OneShowerOneSpike=0;
  int OneShowerOneTrk=0;
  int OneSpikeOneTrk=0;
  int TwoTrk=0;
  int OneShowerLike=0;
  int OneSpikeLike=0;
  int OneTrk=0;
  int ZeroMono=0;		
*/
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
  const double MonoCuts::SwissCrossCut_ = 0.95;
  const double MonoCuts::f51Cut_ = 0.95;
  const double MonoCuts::photonCut_ = 0;
  const double MonoCuts::electronCut_ = 0;
void MonoAnalyzerEGammaLike()
{
	cout<<"hello"<<endl;
	TFile *oFile = new TFile("Monopole_2018_EGLike_test.root","recreate");
	cout<<"new file line pass"<<endl;

	TFile *fin = new TFile("MonoNtuple2018_MC_1000_TestSwisscross2.root");
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
	vector<double> * SwissCross = 0;
        vector<double> * hIso = 0;
        vector<double> * eta = 0;
        vector<double> * phi = 0;
	vector<double> * pho_eta = 0;
	vector<double> * pho_phi = 0;
	vector<double> * pho_pt = 0;
        unsigned nPhoton;
        vector<double> * ele_eta = 0;
        vector<double> * ele_phi = 0;
        vector<double> * ele_pt = 0;
        unsigned nElectron;

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
	
        tree->SetBranchAddress("ele_N",&nElectron);
        tree->SetBranchAddress("ele_eta",&ele_eta);
        tree->SetBranchAddress("ele_phi",&ele_phi);
        tree->SetBranchAddress("ele_pt",&ele_pt);

	double subHits_count=0;
	double subSatHits_count=0;
        

	const unsigned NEvents = tree->GetEntries();
	vector<double> Et;
	
	
	MonoCuts noTrgAnalysis("NOTRG",oFile);
/*	MonoCuts TrgAnalysis_Photon                                 ("HLT_Photon200",oFile);
	MonoCuts TrgAnalysis_PFMET300		                    ("HLT_PFMET300",oFile);
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
        vector<Electron> electron(0);
	
        for(unsigned ev=0; ev<NEvents;ev++){
//		cout<<"Enter event loop"<<endl;
//		cout<<"event number "<<NEvents<<endl;
             tree->GetEntry(ev);
		
		if(nCandidates>Et.size()) Et.resize(nCandidates);
		if(nCandidates>cand.size()) cand.resize(nCandidates);
                if(nElectron>electron.size()) electron.resize(nElectron);
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
                if(nElectron!=0){
                for(unsigned j=0;j<nElectron;j++){
                        electron[j] = Electron(
                        (*ele_eta)[j],
                        (*ele_phi)[j],
                        (*ele_pt)[j]
                        );
                        }
                }
			noTrgAnalysis.doAnalysis(cand,photon,electron,nCandidates,nPhoton,nElectron,true,ev,MET_pt);
				
/*         TrgAnalysis_PFMET300                            .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300,ev);
         TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET170_HBHE_BeamHaloCleaned,ev);
         TrgAnalysis_PFMET140_PFMHT140_IDTight           .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET140_PFMHT140_IDTight,ev);
         TrgAnalysis_PFMET250_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET250_HBHECleaned,ev);
         TrgAnalysis_PFMET300_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300_HBHECleaned,ev);
         TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET200_HBHE_BeamHaloCleaned,ev);
         TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned,ev);
         TrgAnalysis_CaloMET300_HBHECleaned              .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_CaloMET300_HBHECleaned,ev);
         TrgAnalysis_Photon		 		 .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon200,ev);
*/	}//for every event loop
	
 /*			TrgAnalysis_Photon         			  .SignalEff("HLT_Photon200"); 
		        TrgAnalysis_PFMET300		                  .SignalEff("HLT_PFMET300");
		        TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET170_HBHE_BeamHaloCleaned");
		        TrgAnalysis_PFMET140_PFMHT140_IDTight             .SignalEff("HLT_PFMET140_PFMHT140_IDTight");
		        TrgAnalysis_PFMET250_HBHECleaned                  .SignalEff("HLT_PFMET250_HBHECleaned");
		        TrgAnalysis_PFMET300_HBHECleaned                  .SignalEff("HLT_PFMET300_HBHECleaned");
		        TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET200_HBHE_BeamHaloCleaned");
		        TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned  .SignalEff("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned");
		        TrgAnalysis_CaloMET300_HBHECleaned                .SignalEff("HLT_CaloMET300_HBHECleaned");                  


         TrgAnalysis_Photon                              .WritePlots(oFile);
         TrgAnalysis_PFMET300		                 .WritePlots(oFile);       
         TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .WritePlots(oFile);     
         TrgAnalysis_PFMET140_PFMHT140_IDTight           .WritePlots(oFile);
         TrgAnalysis_PFMET250_HBHECleaned                .WritePlots(oFile);   
         TrgAnalysis_PFMET300_HBHECleaned                .WritePlots(oFile);
         TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .WritePlots(oFile);   
         TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.WritePlots(oFile);
         TrgAnalysis_CaloMET300_HBHECleaned              .WritePlots(oFile);   
	                                                     
*/
	noTrgAnalysis.WritePlots(oFile);
	noTrgAnalysis.SignalEff("NOTRG",NEvents);
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
