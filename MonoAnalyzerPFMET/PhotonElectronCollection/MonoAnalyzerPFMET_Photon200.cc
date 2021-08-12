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
//	**The category include e/gamma-like spike-like nonspike-like
//	  too chaos**
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
  PFMETvsOnePhotonOneCand,	
  PFMETvsOneElectronOneCand,	
  PFMETvsTwoCand,
  PFMETvsOneCand,
  PFMETvsZeroMonopole	

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
  bool operator<(const Photon &mc)const{
   if(pho_pt_>mc.pho_pt_) return true;
    else return false;
  }

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
  bool operator<(const Electron &mc)const{
   if(ele_pt_>mc.ele_pt_) return true;
    else return false;
  }

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
 
	x.CreatPlot(ElectronlikePt,new TH1D("ElectronlikePt","",30,0,800));	
	x.CreatPlot(PhotonlikePt,new TH1D("PhotonlikePt","",30,0,800));	
	x.CreatPlot(ElectronlikeEta,new TH1D("ElectronlikeEta","",30,-4,4));	
	x.CreatPlot(PhotonlikeEta,new TH1D("PhotonlikeEta","",30,-4,4));	
	x.CreatPlot(EGlikePt,new TH1D("EGlikePt","",30,0,800));	
	x.CreatPlot(EGlikeEta,new TH1D("EGlikeEta","",30,-4,4));	
	
        x.CreatPlot(PFMETvsTwoPhotonLike,		new TH1D("PFMETvsTwoPhotonLike","",30,0,1100));	
        x.CreatPlot(PFMETvsTwoElectronLike,		new TH1D("PFMETvsTwoElectronLike","",30,0,1100));		
        x.CreatPlot(PFMETvsOnePhotonLike,		new TH1D("PFMETvsOnePhotonLike","",30,0,1100));		
        x.CreatPlot(PFMETvsOneElectronLike,		new TH1D("PFMETvsOneElectronLike","",30,0,1100));		
        x.CreatPlot(PFMETvsOnePhotonOneElectronlike,	new TH1D("PFMETvsOnePhotonOneElectronlike","",30,0,1100));	
        x.CreatPlot(PFMETvsOnePhotonOneCand,		new TH1D("PFMETvsOnePhotonOneCand","",30,0,1100));	
        x.CreatPlot(PFMETvsOneElectronOneCand,		new TH1D("PFMETvsOneElectronOneCand","",30,0,1100));	
        x.CreatPlot(PFMETvsTwoCand,			new TH1D("PFMETvsTwoCand","",30,0,1100));	
        x.CreatPlot(PFMETvsOneCand,			new TH1D("PFMETvsOneCand","",30,0,1100));	
        x.CreatPlot(PFMETvsZeroMonopole,		new TH1D("PFMETvsZeroMonopole","",30,0,1100));	

  }

  ~MonoCuts(){	}
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, vector<Electron> &ele, unsigned nCandidates,unsigned nPhoton, unsigned nElectron, bool TRG, unsigned ev, double PFMET_pt)
  {
	
	//cout<<"enter doAnalysis successfully"<<endl;	
	EcalShower.clear();
	Shower.clear();	

	CutFlowCand_TRG.clear();
        HighPtPhoton.clear();
        HighPtElectron.clear();

	recoMono=0;
	photonLike=0;
	electronLike=0;

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
//	if(TRG && ECut && F51Cut && dEdXCut) CutFlowCand_Dedx.push_back(cands);

	
	//MET
	if(QualCut && dEdXCut && ECut )	   		EcalShower.push_back(cands);//need to add f51<0.85 in the following code)
     
    }//end of cand scan	
   
    if(nPhoton!=0){
    for(unsigned p=0;p<nPhoton;p++){
        Photon &photon = pho[p];
	bool PhotonCut = evalPhoton(photon);
	if(PhotonCut)        HighPtPhoton.push_back(photon);
    }
    }

    if(nElectron!=0){
        for(unsigned p=0;p<nElectron;p++){
                
		Electron &electron = ele[p];
		bool ElectronCut = evalElectron(electron);
                if(ElectronCut)		HighPtElectron.push_back(electron);
        }
    }
	
	sort(HighPtPhoton.begin(),HighPtPhoton.begin()+HighPtPhoton.size());
	sort(HighPtElectron.begin(),HighPtElectron.begin()+HighPtElectron.size());

        sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	
	
	PlotSet &x = Plot[0];
	cout<<"ev "<<ev<<endl;
	
			
	sort(EcalShower.begin(),EcalShower.begin()+EcalShower.size());
	if(EcalShower.size()>0)
	{
//	    cout<<"EcalShower size "<<EcalShower.size()<<endl; 
	    for(int i=0; i<EcalShower.size();i++){
		double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(EcalShower[i].eta_-EcalShower[0].mono_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].mono_phi_,2));
                am_deltaR= sqrt(pow(EcalShower[i].eta_-EcalShower[0].amon_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].amon_phi_,2));
                if(m_deltaR<0.15) 
                {
			EcalShowerNum++;
//			cout<<"Showerlike m in Ecal "<<EcalShowerNum<<endl;
			Shower.push_back(EcalShower[i]);// put matched EcalShower Monopole into "Shower" to match with Egamma
		}
		else if(am_deltaR<0.15)
		{
			EcalShowerNum++;
			Shower.push_back(EcalShower[i]);
//			cout<<"Showerlike am in Ecal "<<EcalShowerNum<<endl;
		}
	    }
	}

	//Replace showe-like with Egamma like
	//	cout<<"Shower size "<<Shower.size()<<endl;


	//-------Photon-like-------//
        //to see monopole is photon-like or spike-like
        //define: photon-monopole angle<0.15 is photon-like
//	cout<<"--------Print out Photon information----------"<<endl;        
//	cout<<"Photon size "<<HighPtPhoton.size()<<endl;		 
       for(int j=0;j<HighPtPhoton.size();j++){
                cout<<j<<endl;
		cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
        }
   	double photonMono=0;
	for(int j=0;j<Shower.size();j++){
		if(HighPtPhoton.size()>0){
		for(int i=0;i<HighPtPhoton.size();i++){
        		double photonMono_deltaR = 0;
                        photonMono_deltaR = sqrt(pow(Shower[j].eta_-HighPtPhoton[i].pho_eta_,2)+
                                                 pow(Shower[j].phi_-HighPtPhoton[i].pho_phi_,2));
                        if(photonMono_deltaR < 0.15 ){
//                                cout<<"deltaR "<<photonMono_deltaR<<"  photon-like monopole"<<endl;
					if(Shower[j].SwissCross_ > 0.95){
					    if( Shower[j].eta_ <1.479 && Shower[j].eta_ > -1.479) 
						EMShowerButSpikelikeinBarrel++;
					    else
                	                        EMShowerButSpikelikeinEE++;
					}
	        	                    
					photonLike++;// otherwise not photon-like monopole
					EgammaLike++;
					TotalPhotonLike++;
					photonMono = photonMono_deltaR;
	//				x.GetPlot(PhotonlikePt)->Fill(HighPtPhoton[i].pho_pt_);
	//				x.GetPlot(PhotonlikeEta)->Fill(HighPtPhoton[i].pho_eta_);
					x.GetPlot(EGlikePt)->Fill(HighPtPhoton[i].pho_pt_);
					x.GetPlot(EGlikeEta)->Fill(HighPtPhoton[i].pho_eta_);
			}
			else{
//				cout<<"deltaR "<<photonMono_deltaR<<"  No photon-like"<<endl;	
			}	
                }
		}
	}
	//-------Electron-like-------//
        //to see monopole is electron-like or spike-like
        //define: electron-monopole angle<0.15 is electron-like
//	cout<<"-------Print out electron information----------"<<endl;        
//	cout<<"Electron size "<<HighPtElectron.size()<<endl;		 
        for(int i=0;i<HighPtElectron.size();i++){
                cout<<i<<endl;
		cout<<"ele  eta "<<HighPtElectron[i].ele_eta_<<endl;
                cout<<"ele  phi "<<HighPtElectron[i].ele_phi_<<endl;
                cout<<"ele   pt "<<HighPtElectron[i].ele_pt_<<endl;
        }
	cout<<endl;

	double electronMono=0;
	for(int j=0;j<Shower.size();j++){
                if(HighPtElectron.size()!=0){ 
       		for(int i=0;i<HighPtElectron.size();i++){
        		double electronMono_deltaR = 0;
        	        electronMono_deltaR = sqrt(pow(Shower[j].eta_-HighPtElectron[i].ele_eta_,2)+
        	                                   pow(Shower[j].phi_-HighPtElectron[i].ele_phi_,2));
        	        if(electronMono_deltaR < 0.15 )
			{
//        		                cout<<"deltaR "<<electronMono_deltaR<<"  electron-like monopole"<<endl;
					if(Shower[j].SwissCross_ > 0.95){
					    if( Shower[j].eta_ <1.479 && Shower[j].eta_ > -1.479) 
						EMShowerButSpikelikeinBarrel++;
					    else
                	                        EMShowerButSpikelikeinEE++;
					}
					electronLike++;
					EgammaLike++;
					electronMono = electronMono_deltaR;
			
					if(abs(electronMono-photonMono)>0.1){
						x.GetPlot(ElectronlikePt)->Fill(HighPtElectron[i].ele_pt_);
						x.GetPlot(ElectronlikeEta)->Fill(HighPtElectron[i].ele_eta_);
						x.GetPlot(EGlikePt)->Fill(HighPtElectron[i].ele_pt_);
						x.GetPlot(EGlikeEta)->Fill(HighPtElectron[i].ele_eta_);
						TotalElectronLike++;
					}
			}
	
        	        else{        
//				cout<<"deltaR "<<electronMono_deltaR<<"  No electron-like"<<endl;
			}
		}
		}
	}

	if(electronLike==0 && photonLike==0 && Shower.size()==2)	recoMono=2;
        else if(electronLike==0 && photonLike==0 && Shower.size()==1)        recoMono=1;
        else if(electronLike==1 && photonLike==0 && Shower.size()==2)        recoMono=1;
        else if(electronLike==0 && photonLike==1 && Shower.size()==2)        recoMono=1;
        else if(electronLike==1 && photonLike==1 && Shower.size()==1)        electronLike=0;
	
	if(electronLike==1 && photonLike==1 && Shower.size()==2){	
		if(abs(electronMono-photonMono)<0.1){
			electronLike=0;
			recoMono=1;
		}
		else cout<<"WWWWWTTTTTFFFFF"<<endl;	
	}
	else if( electronLike==1 && photonLike==2 && Shower.size()==2){
		if(abs(electronMono-photonMono)<0.1){
                	photonLike=1;	
		}
		else cout<<"WWWWTTTTFFFF2"<<endl;
	}
	else if( electronLike==2 && photonLike==2 && Shower.size()==2){
		if(abs(electronMono-photonMono)<0.1){
                	photonLike=1;
			electronLike=1;
		}
		else cout<<"WWWWTTTTFFFF3"<<endl;
	}
	else if( electronLike==2 && photonLike==1 && Shower.size()==2){
		if(abs(electronMono-photonMono)<0.1){
                	photonLike=1;
			electronLike=1;
		}
		else cout<<"WWWWTTTTFFFF4"<<endl;
	}


	cout<<"electronLike "<<electronLike<<endl;
	cout<<"photonLike   "<<photonLike<<endl;
	cout<<"reco Mono    "<<recoMono<<endl;





	//calculate the efficiency of egamma-like event
	if(PFMET_pt>200){
	if(     recoMono ==0 && electronLike ==0 && photonLike == 2  ){
		x.GetPlot(PhotonlikePt)->Fill(HighPtPhoton[0].pho_pt_);
		x.GetPlot(PhotonlikePt)->Fill(HighPtPhoton[0].pho_pt_);
		x.GetPlot(PhotonlikeEta)->Fill(HighPtPhoton[1].pho_eta_);
		x.GetPlot(PhotonlikeEta)->Fill(HighPtPhoton[1].pho_eta_);
	//	if(PFMET_pt > 200){ 
		TwoPhotonLikeEvent++;
		x.GetPlot(PFMETvsTwoPhotonLike)->Fill(PFMET_pt);
	//	}
	}
	else if(recoMono ==0 && electronLike ==0 && photonLike == 1 ){
		x.GetPlot(PhotonlikePt)->Fill(HighPtPhoton[0].pho_pt_);
		x.GetPlot(PhotonlikeEta)->Fill(HighPtPhoton[0].pho_eta_);
	//	if(PFMET_pt > 200){ 
		 OnePhotonLikeEvent++;
		x.GetPlot(PFMETvsOnePhotonLike)->Fill(PFMET_pt);
	//	}
	}
	else if(recoMono ==0 && electronLike ==2 && photonLike == 0 ){
	//	if(PFMET_pt > 200){ 
		 TwoElectronLikeEvent++;
		x.GetPlot(PFMETvsTwoElectronLike)->Fill(PFMET_pt);
	//	}
	}
	else if(recoMono ==0 && electronLike ==1 && photonLike == 0 ){
	//	if(PFMET_pt > 200){ 
		 OneElectronLikeEvent++;
		x.GetPlot(PFMETvsOneElectronLike)->Fill(PFMET_pt);
	//	}
	}
	else if(recoMono ==0 && electronLike ==1 && photonLike == 1 ){
		x.GetPlot(PhotonlikePt)->Fill(HighPtPhoton[0].pho_pt_);
		x.GetPlot(PhotonlikeEta)->Fill(HighPtPhoton[0].pho_eta_);
	//	if(PFMET_pt > 200){ 
		OnePhotonOneElectronlikeEvent++;		
		x.GetPlot(PFMETvsOnePhotonOneElectronlike)->Fill(PFMET_pt);
	//	}
	}
	else if(recoMono ==2 && electronLike ==0 && photonLike == 0 ){
	//	if(PFMET_pt > 200){ 
		TwoCandEvent++;
		x.GetPlot(PFMETvsTwoCand)->Fill(PFMET_pt);
	//	}
	}
	else if(recoMono ==1 && electronLike ==0 && photonLike == 0 ){
	//	if(PFMET_pt > 200){ 
		OneCandEvent++;
		x.GetPlot(PFMETvsOneCand)->Fill(PFMET_pt);
	//	}
	}
        else if(recoMono ==1 && electronLike ==1 && photonLike == 0 ){
	//	if(PFMET_pt > 200){ 
                OneElectronOneCandEvent++;
		x.GetPlot(PFMETvsOneElectronOneCand)->Fill(PFMET_pt);
	//	}
        }
        else if(recoMono ==1 && electronLike ==0 && photonLike == 1 ){
		x.GetPlot(PhotonlikePt)->Fill(HighPtPhoton[0].pho_pt_);
		x.GetPlot(PhotonlikeEta)->Fill(HighPtPhoton[0].pho_eta_);
	//	if(PFMET_pt > 200){ 
                OnePhotonOneCandEvent++;
		x.GetPlot(PFMETvsOnePhotonOneCand)->Fill(PFMET_pt);
	//	}
        }
	}
	//PFMET_pt < 200
        else if(recoMono ==0 && electronLike ==0 && photonLike == 0 ){
	//	if(PFMET_pt > 200){ 
		ZeroMonopoleEvent++;
		x.GetPlot(PFMETvsZeroMonopole)->Fill(PFMET_pt);
	//	}
        }
	
	else{
	//	if(PFMET_pt > 200){ 
		cout<<"there are some accident..."<<endl;
		      cout<<"electronLike "<<electronLike<<endl;
		      cout<<"photonLike   "<<photonLike<<endl;
		      cout<<"reco Mono    "<<recoMono<<endl;
		
		      
  	     for(int j=0;j<HighPtPhoton.size();j++){
                cout<<j<<endl;
		cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
        	}
	     for(int i=0;i<HighPtElectron.size();i++){
                cout<<i<<endl;
		cout<<"ele  eta "<<HighPtElectron[i].ele_eta_<<endl;
                cout<<"ele  phi "<<HighPtElectron[i].ele_phi_<<endl;
                cout<<"ele   pt "<<HighPtElectron[i].ele_pt_<<endl;
        	}

		for(int i=0;i<Shower.size();i++){
                double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(Shower[i].eta_-Shower[0].mono_eta_,2)+
                                pow(Shower[i].phi_-Shower[0].mono_phi_,2));
                am_deltaR= sqrt(pow(Shower[i].eta_-Shower[0].amon_eta_,2)+
                                pow(Shower[i].phi_-Shower[0].amon_phi_,2));
			cout<<i<<endl;
			cout<<"     eta "<<setprecision(5)<<Shower[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<Shower[i].phi_<<endl;
			cout<<"      Et "<<setprecision(5)<<Shower[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<Shower[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<Shower[i].f51_<<endl;
			cout<<"     f15 "<<setprecision(5)<<Shower[i].f15_<<endl;
			cout<<"   Cross "<<setprecision(5)<<Shower[i].SwissCross_<<endl;
			cout<<"     MET "<<setprecision(5)<<Shower[i].MET_pt_<<endl;
			RealNum++;
		
              		cout<<"          candidate           monoGen         antimonoGen"<<endl;
                	cout<<"eta      "<<setprecision(5)<<Shower[i].eta_<<setw(20)<<Shower[i].mono_eta_
                        	         <<setw(20)<<Shower[i].amon_eta_<<endl;
                	cout<<"phi      "<<setprecision(5)<<Shower[i].phi_<<setw(20)<<Shower[i].mono_phi_
                       		         <<setw(20)<<Shower[i].amon_phi_<<endl;
                	cout<<"m deltaR "<<m_deltaR<<endl;
                	cout<<"a deltaR "<<am_deltaR<<endl;
			cout<<"--------------------------------------------"<<endl;
		}	
	//	}

	}	
	
	
	for(int i=0;i<EcalShower.size();i++){
                double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(EcalShower[i].eta_-Shower[0].mono_eta_,2)+
                                pow(EcalShower[i].phi_-Shower[0].mono_phi_,2));
                am_deltaR= sqrt(pow(EcalShower[i].eta_-Shower[0].amon_eta_,2)+
                                pow(EcalShower[i].phi_-Shower[0].amon_phi_,2));
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
		
              		cout<<"          candidate           monoGen         antimonoGen"<<endl;
                	cout<<"eta      "<<setprecision(5)<<EcalShower[i].eta_<<setw(20)<<Shower[i].mono_eta_
                        	         <<setw(20)<<EcalShower[i].amon_eta_<<endl;
                	cout<<"phi      "<<setprecision(5)<<EcalShower[i].phi_<<setw(20)<<Shower[i].mono_phi_
                       		         <<setw(20)<<EcalShower[i].amon_phi_<<endl;
                	cout<<"m deltaR "<<m_deltaR<<endl;
                	cout<<"a deltaR "<<am_deltaR<<endl;
			cout<<"--------------------------------------------"<<endl;
	}

		
        	cout<<"event tag"<<CutFlowCand_TRG[0].event_<<endl;
	//		cout<<"=================================================="<<endl;		
			cout<<"============================================="<<endl;		

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
//	cout<<"# of EcalShower : pass dEdx>9, E55 > 0, NO Quality, f51..."<<endl;
//	cout<<"EcalShowerNum "<<EcalShowerNum<<endl;
//        cout<<"photonLike "<<photonLike<<endl;
//        cout<<"electronLike "<<electronLike<<endl;
//        cout<<"ratio of photonlike "<<(double)photonLike/(double)EcalShowerNum<<endl;
//        cout<<"ratio of electronlike "<<(double)electronLike/(double)EcalShowerNum<<endl;
	cout<<"--------------------------------------------"<<endl;
        cout<<"TwoPhotonLikeEvent		"<<(double)TwoPhotonLikeEvent			<<endl;		
        cout<<"TwoElectronLikeEven		"<<(double)TwoElectronLikeEvent			<<endl;	
        cout<<"OnePhotonLikeEvent		"<<(double)OnePhotonLikeEvent			<<endl;	
        cout<<"OneElectronLikeEvent		"<<(double)OneElectronLikeEvent			<<endl;	
        cout<<"OnePhotonOneElectronlikeEvent	"<<(double)OnePhotonOneElectronlikeEvent	<<endl;	
	cout<<"OnePhotonOneCandEvent    "<<(double)OnePhotonOneCandEvent<<endl;
	cout<<"OneElectronOneCandEvent  "<<(double)OneElectronOneCandEvent<<endl;	 
	cout<<"TwoCandEvent	        "<<(double)TwoCandEvent<<endl;
	cout<<"OneCandEvent		"<<(double)OneCandEvent<<endl;
	cout<<"ZeroMonopoleEvent 	"<<(double)ZeroMonopoleEvent<<endl;	
	cout<<endl;
	cout<<"Pass PFMET > 200 Total events  "<<TwoPhotonLikeEvent+TwoElectronLikeEvent+OnePhotonLikeEvent+OneElectronLikeEvent+OnePhotonOneElectronlikeEvent
						+OnePhotonOneCandEvent+OneElectronOneCandEvent+TwoCandEvent+OneCandEvent<<endl;

	cout<<endl<<"TotalPhotonLike event	"<<TwoPhotonLikeEvent+OnePhotonLikeEvent+OnePhotonOneElectronlikeEvent+OnePhotonOneCandEvent<<endl;
	cout<<"TotalElectronLike	"<<TwoElectronLikeEvent+OneElectronLikeEvent+OnePhotonOneElectronlikeEvent+OneElectronOneCandEvent<<endl;
//	cout<<endl<<"EM shower but spike like in EB "<<EMShowerButSpikelikeinBarrel;
//	cout<<endl<<"EM shower but spike like in EE "<<EMShowerButSpikelikeinEE;
//	cout<<endl<<"**Note that No spike def in EE, only in Ecal, so mainly see how many spike-like EMshower in EB"<<endl;
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
  int recoMono=0;
	int TwoPhotonLikeEvent = 0;	
	int TwoElectronLikeEvent = 0;		
	int OnePhotonLikeEvent = 0;	
	int OneElectronLikeEvent = 0;		
	int OnePhotonOneElectronlikeEvent = 0;	
	int OnePhotonOneCandEvent = 0;	
	int OneElectronOneCandEvent = 0;		
	int TwoCandEvent = 0;	
	int OneCandEvent = 0;
	int ZeroMonopoleEvent = 0;		

int TotalPhotonLike=0;
int TotalElectronLike=0;

int EMShowerButSpikelikeinBarrel=0;
int EMShowerButSpikelikeinEE=0;


  int TrackerMono=0;
  int EcalShowerNum=0;
 
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
  const double MonoCuts::f51Cut_ = 0.85;
  const double MonoCuts::photonCut_ = 0;
  const double MonoCuts::electronCut_ = 0;
void MonoAnalyzerPFMET_Photon200()
{
	cout<<"hello"<<endl;
	TFile *oFile = new TFile("Monopole_PFMET_Trigger_1000.root","recreate");
	cout<<"new file line pass"<<endl;

	TFile *fin = new TFile("MonoNtuple2018_MC_1000.root");
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
//	MonoCuts TrgAnalysis_Photon                                 ("HLT_Photon200",oFile);
//	MonoCuts TrgAnalysis_PFMET300		                    ("HLT_PFMET300",oFile);
//	MonoCuts TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned          ("HLT_PFMET170_HBHE_BeamHaloCleaned",oFile);
//	MonoCuts TrgAnalysis_PFMET140_PFMHT140_IDTight              ("HLT_PFMET140_PFMHT140_IDTight",oFile);
//	MonoCuts TrgAnalysis_PFMET300_HBHECleaned                   ("HLT_PFMET300_HBHECleaned",oFile);
	MonoCuts TrgAnalysis_PFMET250_HBHECleaned                   ("HLT_PFMET250_HBHECleaned",oFile);
	MonoCuts TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned          ("HLT_PFMET200_HBHE_BeamHaloCleaned",oFile);
	MonoCuts TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned   ("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",oFile);
	MonoCuts TrgAnalysis_CaloMET300_HBHECleaned                 ("HLT_CaloMET300_HBHECleaned",oFile);

	vector<MonoCandidate> cand(10);	
	vector<Photon> photon(0);
        vector<Electron> electron(0);
	
        for(unsigned ev=0; ev<NEvents;ev++){
	//	cout<<"Enter event loop"<<endl;
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
	//	cout<<"enter photon loop"<<endl;
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
				
//         TrgAnalysis_PFMET300                            .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300,ev);
//         TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET170_HBHE_BeamHaloCleaned,ev);
//         TrgAnalysis_PFMET140_PFMHT140_IDTight           .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET140_PFMHT140_IDTight,ev);
//         TrgAnalysis_PFMET300_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET300_HBHECleaned,ev);

//         TrgAnalysis_PFMET250_HBHECleaned                .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET250_HBHECleaned,ev);
//         TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMET200_HBHE_BeamHaloCleaned,ev);
//         TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned,ev);
//         TrgAnalysis_CaloMET300_HBHECleaned              .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_CaloMET300_HBHECleaned,ev);
//         TrgAnalysis_Photon		 		 .doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon200,ev);
	}//for every event loop
	
// 			TrgAnalysis_Photon         			  .SignalEff("HLT_Photon200"); 
//		        TrgAnalysis_PFMET300		                  .SignalEff("HLT_PFMET300");
//		        TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET170_HBHE_BeamHaloCleaned");
//		        TrgAnalysis_PFMET140_PFMHT140_IDTight             .SignalEff("HLT_PFMET140_PFMHT140_IDTight");
//		        TrgAnalysis_PFMET300_HBHECleaned                  .SignalEff("HLT_PFMET300_HBHECleaned");

//		        TrgAnalysis_PFMET250_HBHECleaned                  .SignalEff("HLT_PFMET250_HBHECleaned");
//		        TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned         .SignalEff("HLT_PFMET200_HBHE_BeamHaloCleaned");
//		        TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned  .SignalEff("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned");
//		        TrgAnalysis_CaloMET300_HBHECleaned                .SignalEff("HLT_CaloMET300_HBHECleaned");                  


//         TrgAnalysis_Photon                              .WritePlots(oFile);
//         TrgAnalysis_PFMET300		                 .WritePlots(oFile);       
//         TrgAnalysis_PFMET170_HBHE_BeamHaloCleaned       .WritePlots(oFile);     
//         TrgAnalysis_PFMET140_PFMHT140_IDTight           .WritePlots(oFile);
//         TrgAnalysis_PFMET300_HBHECleaned                .WritePlots(oFile);

//         TrgAnalysis_PFMET250_HBHECleaned                .WritePlots(oFile);   
//         TrgAnalysis_PFMET200_HBHE_BeamHaloCleaned       .WritePlots(oFile);   
//         TrgAnalysis_PFMETTypeOne200_HBHE_BeamHaloCleaned.WritePlots(oFile);
//         TrgAnalysis_CaloMET300_HBHECleaned              .WritePlots(oFile);   
	                                                     

	noTrgAnalysis.WritePlots(oFile);
	noTrgAnalysis.SignalEff("NOTRG",NEvents);
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
