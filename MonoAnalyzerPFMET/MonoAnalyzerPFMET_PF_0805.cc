//////////////////////////////////////////////////////
//
//	File name: MonoAnalyzerPFMET_PF.cc
//	author: Shih Lin
//	Purpose: 
//		 ** The previous analysis use Photon/Electron collection
//		 but PFMET use "PF" photon/electron, so change matching method.
//
//		 THe cause of large MET is one of monopoles 
//		 identify as e/gamma. This script is to understand
//		 how many monopole are identified as PF and cause MET.
//		 Also calcualte the efficiency of e/gamma like.
//		 (MET is not deal with SPIKE)
//	output:Monopole_PFMET_PF_*000.root
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
#include <fstream>
using namespace std;
enum PlotName{
	PFMETvsMonopoleCand,
	ParticleflowlikePt,
	ParticleflowlikeEta,
	pfElectronLikePt,
	pfPhotonLikePt,
	pfUnknownLikePt,
	
	
	PFMETvsOneCandOne1ChargedHadron,
	PFMETvsOneCandOne2Electron,
	PFMETvsOneCandOne3Muon,
	PFMETvsOneCandOne4Photon,
	PFMETvsOneCandOne5NeutralHadron,
	PFMETvsOneCandOne6HadronicEnergyinHF,
	PFMETvsOneCandOne7EMEnergyinHF,
	
	PFMETvsTwo1ChargedHadron,
	PFMETvsTwo2Electron,
	PFMETvsTwo3Muon,
	PFMETvsTwo4Photon,
	PFMETvsTwo5NeutralHadron,
	PFMETvsTwo6HadronicEnergyinHF,
	PFMETvsTwo7EMEnergyinHF,

	PFMETvsOnlyOne1ChargedHadron,
	PFMETvsOnlyOne2Electron,
	PFMETvsOnlyOne3Muon,
	PFMETvsOnlyOne4Photon,
	PFMETvsOnlyOne5NeutralHadron,
	PFMETvsOnlyOne6HadronicEnergyinHF,
	PFMETvsOnlyOne7EMEnergyinHF,

	PFMETvsZeroCand,
	PFMETvsOnlyOneCand,
	PFMETvsTwoCand
	
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
  void CreatProfile(const PlotName pn, TProfile* h){ 
	profile_[pn] = h;
  }
  TH1 * GetPlot(const PlotName pn){ return plots_[pn]; }
  TProfile *GetProfile(const PlotName f){ return profile_[f]; }
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
  vector<TProfile*> profile_;
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
class Particleflow
{
public:
  Particleflow(){}
  Particleflow(double pf_eta, double pf_phi, double pf_pt, double pf_E, double pf_pdgId):pf_eta_(pf_eta),pf_phi_(pf_phi),
										pf_pt_(pf_pt),pf_E_(pf_E),pf_pdgId_(pf_pdgId){}
  Particleflow(const Particleflow &mc):pf_eta_(mc.pf_eta_),pf_phi_(mc.pf_phi_),pf_pt_(mc.pf_pt_),pf_E_(mc.pf_E_),pf_pdgId_(mc.pf_pdgId_){}
  double pf_eta_;
  double pf_phi_;
  double pf_pt_;
  double pf_E_;
  double pf_pdgId_;

  ~Particleflow(){}

  bool operator<(const Particleflow &mc)const{
   if(pf_pt_>mc.pf_pt_) return true;
   else if(pf_pt_==mc.pf_pt_){
        if(pf_E_>mc.pf_E_) return true;
        else return false;
        }
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
 	
	//This is monopole candidate pass preselection, dEdxSig, f51, e55

	x.CreatPlot(PFMETvsMonopoleCand,new TH1D("PFvsMonopoleCand","",30,0,1100));
	
	
//	x.CreatPlot(ElectronlikePt,new TH1D("ElectronlikePt","",30,0,800));	
//	x.CreatPlot(PhotonlikePt,new TH1D("PhotonlikePt","",30,0,800));	
//	x.CreatPlot(ElectronlikeEta,new TH1D("ElectronlikeEta","",30,-4,4));	
//	x.CreatPlot(PhotonlikeEta,new TH1D("PhotonlikeEta","",30,-4,4));	
//	x.CreatPlot(EGlikePt,new TH1D("EGlikePt","",30,0,800));	
//	x.CreatPlot(EGlikeEta,new TH1D("EGlikeEta","",30,-4,4));

//	x.CreatPlot(ParticleflowlikePt,new TH1D("ParticleflowlikePt","",30,0,800));
//	x.CreatPlot(ParticleflowlikeEta,new TH1D("ParticleflowlikeEta","",30,-4,4));
//	x.CreatPlot(pfElectronLikePt,new TH1D("pfElectronLikePt","",30,0,800));
//	x.CreatPlot(pfPhotonLikePt,new TH1D("pfPhotonLikePt","",30,0,800));
//	x.CreatPlot(pfUnknownLikePt,new TH1D("pfUnknownLikePt","",30,0,800));
		
        x.CreatPlot(PFMETvsOnlyOne1ChargedHadron,	new TH1D("PFMETvsOnlyOne1ChargedHadron","",30,0,1100));	
        x.CreatPlot(PFMETvsOnlyOne2Electron,     	new TH1D("PFMETvsOnlyOne2Electron,","",30,0,1100));		
        x.CreatPlot(PFMETvsOnlyOne3Muon,        	new TH1D("PFMETvsOnlyOne3Muon","",30,0,1100));		
        x.CreatPlot(PFMETvsOnlyOne4Photon,       	new TH1D("PFMETvsOnlyOne4Photon","",30,0,1100));		
        x.CreatPlot(PFMETvsOnlyOne5NeutralHadron,	new TH1D("PFMETvsOnlyOne5NeutralHadron","",30,0,1100));	
        x.CreatPlot(PFMETvsOnlyOne6HadronicEnergyinHF,  new TH1D("PFMETvsOnlyOne6HadronicEnergyinHF","",30,0,1100));	
        x.CreatPlot(PFMETvsOnlyOne7EMEnergyinHF,	new TH1D("PFMETvsOnlyOne7EMEnergyinHF","",30,0,1100));	
	
        x.CreatPlot(PFMETvsOneCandOne1ChargedHadron,	new TH1D("PFMETvsOneCandOne1ChargedHadron","",30,0,1100));	
        x.CreatPlot(PFMETvsOneCandOne2Electron,     	new TH1D("PFMETvsOneCandOne2Electron,","",30,0,1100));		
        x.CreatPlot(PFMETvsOneCandOne3Muon,        	new TH1D("PFMETvsOneCandOne3Muon","",30,0,1100));		
        x.CreatPlot(PFMETvsOneCandOne4Photon,       	new TH1D("PFMETvsOneCandOne4Photon","",30,0,1100));		
        x.CreatPlot(PFMETvsOneCandOne5NeutralHadron,	new TH1D("PFMETvsOneCandOne5NeutralHadron","",30,0,1100));	
        x.CreatPlot(PFMETvsOneCandOne6HadronicEnergyinHF,new TH1D("PFMETvsOneCandOne6HadronicEnergyinHF","",30,0,1100));	
        x.CreatPlot(PFMETvsOneCandOne7EMEnergyinHF,	new TH1D("PFMETvsOneCandOne7EMEnergyinHF","",30,0,1100));	
        
	x.CreatPlot(PFMETvsTwo1ChargedHadron,	new TH1D("PFMETvsTwo1ChargedHadron","",30,0,1100));	
        x.CreatPlot(PFMETvsTwo2Electron,     	new TH1D("PFMETvsTwo2Electron,","",30,0,1100));		
        x.CreatPlot(PFMETvsTwo3Muon,        	new TH1D("PFMETvsTwo3Muon","",30,0,1100));		
        x.CreatPlot(PFMETvsTwo4Photon,       	new TH1D("PFMETvsTwo4Photon","",30,0,1100));		
        x.CreatPlot(PFMETvsTwo5NeutralHadron,	new TH1D("PFMETvsTwo5NeutralHadron","",30,0,1100));	
        x.CreatPlot(PFMETvsTwo6HadronicEnergyinHF, new TH1D("PFMETvsTwo6HadronicEnergyinHF","",30,0,1100));	
        x.CreatPlot(PFMETvsTwo7EMEnergyinHF,	new TH1D("PFMETvsTwo7EMEnergyinHF","",30,0,1100));	
        x.CreatPlot(PFMETvsTwoCand,			new TH1D("PFMETvsTwoCand","",30,0,1100));	
        x.CreatPlot(PFMETvsOnlyOneCand,			new TH1D("PFMETvsOnlyOneCand","",30,0,1100));	
        x.CreatPlot(PFMETvsZeroCand,		new TH1D("PFMETvsZeroCand","",30,0,1100));	

  }

  ~MonoCuts(){	}
  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, vector<Electron> &ele,vector<Particleflow> &pf, unsigned nCandidates,unsigned nPhoton, unsigned nElectron, unsigned nParticleflow, bool TRG, unsigned ev, double PFMET_pt, double PFMET_phi)
  {
	
//	cout<<"enter doAnalysis successfully"<<endl;	
	EcalShower.clear();
	Shower.clear();	
	AntiShower.clear();	

	CutFlowCand_TRG.clear();
	CutFlowCand_Dedx.clear();
        HighPtPhoton.clear();
        HighPtElectron.clear();
	HighPtParticleflow.clear();
	pfPhotonlike.clear();
	pfElectronlike.clear();
	pfUnknown.clear();
	pfAntiUnknown.clear();

	pfParticle_Mono_deltaR.clear() ;
	pfParticle_Amon_deltaR.clear() ;
	recoMono=0;
	pfMono=0;
	pfAmon=0;
     for(unsigned c=0;c<nCandidates;c++){
	MonoCandidate &cands = cand[c];
	bool QualCut = evalQuality(cands);
	bool ECut = evalE(cands);
	bool F15Cut = evalF15(cands);
	bool F51Cut = evalF51(cands);
	bool SwissCrossCut = evalSwissCross(cands);
	bool dEdXCut = evaldEdX(cands);
	bool PFMETCut = evalPFMET(cands); 
	//no cut
	if(TRG) CutFlowCand_TRG.push_back(cands);

	//cut all 
	if(TRG && QualCut && ECut && F51Cut && dEdXCut) CutFlowCand_Dedx.push_back(cands);

	
	//MET
	if(dEdXCut &&  ECut  )	 EcalShower.push_back(cands);
     
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
	
    if(nParticleflow!=0){
        for(unsigned p=0;p<nParticleflow;p++){
                
		Particleflow &particleflow = pf[p];
		bool ParticleflowCut = evalParticleflow(particleflow);
                if(ParticleflowCut)	HighPtParticleflow.push_back(particleflow);
        }
    }


        sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	sort(HighPtParticleflow.begin(),HighPtParticleflow.begin()+HighPtParticleflow.size());

        sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());
	PlotSet &x = Plot[0];
	cout<<"ev "<<ev<<endl;
	//Just want to see monopole candidate(pass all selection) vs PFMET	
	for(int i=0; i<CutFlowCand_Dedx.size();i++){
		double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].mono_eta_,2)+
                                pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].mono_phi_,2));
                am_deltaR= sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].amon_eta_,2)+
                                pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].amon_phi_,2));
                if(m_deltaR<0.15||am_deltaR<0.15) 
                {
			x.GetPlot(PFMETvsMonopoleCand)->Fill(PFMET_pt);	
		}
	}		

	//First, find the candidates that pass the selection
	//If the cands also match to GenMonopole, they defines as "Shower" and "AntiShower"
	sort(EcalShower.begin(),EcalShower.begin()+EcalShower.size());
	if(EcalShower.size()>0)
	{
	    for(int i=0; i<EcalShower.size();i++){
		double m_deltaR=0;
                double am_deltaR=0;
                m_deltaR = sqrt(pow(EcalShower[i].eta_-EcalShower[0].mono_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].mono_phi_,2));
                am_deltaR= sqrt(pow(EcalShower[i].eta_-EcalShower[0].amon_eta_,2)+
                                pow(EcalShower[i].phi_-EcalShower[0].amon_phi_,2));
                if(m_deltaR<0.15) 
                {
			Shower.push_back(EcalShower[i]);// put matched EcalShower Monopole into "Shower" to match with Egamma
		}
		else if(am_deltaR<0.15)
		{
			AntiShower.push_back(EcalShower[i]);
		}
	    }
	}
	if(Shower.size()>1) Shower.pop_back();
	if(AntiShower.size()>1) AntiShower.pop_back();
	//to look monopole correspond to what particle in the PF
	//Don't mind whether monopole is photon or not, it could be anything in the PF

	/*----------Match PF with monopole----------*/
	for(int j=0;j<Shower.size();j++){
       		for(int i=0;i<HighPtParticleflow.size();i++){
				
			double R = 0;
        	        R = sqrt(pow(Shower[j].eta_-HighPtParticleflow[i].pf_eta_,2)+
  	      	                 pow(Shower[j].phi_-HighPtParticleflow[i].pf_phi_,2));
			if( R  < 0.15)
			{
				pfParticle_Mono_deltaR.push_back(R);
				pfUnknown.push_back(HighPtParticleflow[i]);
			}
		}
	}
	//I only care about if there is any PFCand MATCH to recoMono
	//Since one monopole could be break into several PFCands.
	//Thus, if there at least one PFcand match to reco Mono,
	//I just mark it as pfMono=1 !!!!
	if(pfUnknown.size()!=0) pfMono=1;
	
	for(int j=0;j<AntiShower.size();j++){
       		for(int i=0;i<HighPtParticleflow.size();i++){
				
			double R = 0;
        	        R = sqrt(pow(AntiShower[j].eta_-HighPtParticleflow[i].pf_eta_,2)+
  	      	                 pow(AntiShower[j].phi_-HighPtParticleflow[i].pf_phi_,2));
			if( R  < 0.15)
			{
				pfParticle_Amon_deltaR.push_back(R);
				pfAntiUnknown.push_back(HighPtParticleflow[i]);
			}
		}
	}
	if(pfAntiUnknown.size()!=0) pfAmon = 1;
	
	/*--------Print-----------*/
	for(int i=0;i<pfUnknown.size();i++){	
		cout<<i<<" Mono -------------------------------------"<<endl;
	        cout<<setw(12)<<"eta"<<setw(9)<<"phi"<<setw(9)<<"pT"<<setw(9)<<"E"<<setw(9)<<"ID"<<endl;
	        cout<<"PFCand  "<<setprecision(4)<<pfUnknown[i].pf_eta_<<setw(9)
		        <<pfUnknown[i].pf_phi_<<setw(9)<<pfUnknown[i].pf_pt_<<setw(9)
			<<pfUnknown[i].pf_E_<<setw(6)
			<<pfUnknown[i].pf_pdgId_<<endl;
		cout<<"deltaR "<<pfParticle_Mono_deltaR[i]<<endl;
	}
	for(int i=0;i<pfAntiUnknown.size();i++){	
		cout<<i<<" Amon --------------------------------------"<<endl;
	        cout<<setw(12)<<"eta"<<setw(9)<<"phi"<<setw(9)<<"pT"<<setw(9)<<"E"<<setw(9)<<"ID"<<endl;
	        cout<<"PFCand  "<<setprecision(4)<<pfAntiUnknown[i].pf_eta_<<setw(9)
		        <<pfAntiUnknown[i].pf_phi_<<setw(9)<<pfAntiUnknown[i].pf_pt_<<setw(9)
			<<pfAntiUnknown[i].pf_E_<<setw(6)
			<<pfAntiUnknown[i].pf_pdgId_<<endl;
		cout<<"deltaR "<<pfParticle_Amon_deltaR[i]<<endl;
	}

	//pfUnknown should have only 2 particles
	//If pfUnknown have only 1 particle, and it also counted into PF, it should have large MET.
	//MET large case 1
//if(PFMET_pt>200){	

	cout<<setw(5)<<" ---------------------"<<endl;		
	cout<<setw(5)<<" | Shower          "<<Shower.size()<<" |"<<endl;
	cout<<setw(5)<<" | AntiShower 	   "<<AntiShower.size()<<" |"<<endl;
	cout<<setw(5)<<" | pfMono          "<<pfMono<<" |"<<endl;
	cout<<setw(5)<<" | pfAmon          "<<pfAmon<<" |"<<endl;
	cout<<setw(5)<<" ---------------------"<<endl;		


	if(Shower.size()==1 && AntiShower.size() == 1 && pfMono==1 && pfAmon==1){
		TwoPfEvent++;		
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | MET should be small"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;

			if(pfUnknown[0].pf_pdgId_ == 1) x.GetPlot(PFMETvsTwo1ChargedHadron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 2) x.GetPlot(PFMETvsTwo2Electron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 3) x.GetPlot(PFMETvsTwo3Muon)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 4) x.GetPlot(PFMETvsTwo4Photon)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 5) x.GetPlot(PFMETvsTwo5NeutralHadron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 6) x.GetPlot(PFMETvsTwo6HadronicEnergyinHF)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 7) x.GetPlot(PFMETvsTwo7EMEnergyinHF)->Fill(PFMET_pt);
			
			if(pfAntiUnknown[0].pf_pdgId_ == 1) x.GetPlot(PFMETvsTwo1ChargedHadron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 2) x.GetPlot(PFMETvsTwo2Electron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 3) x.GetPlot(PFMETvsTwo3Muon)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 4) x.GetPlot(PFMETvsTwo4Photon)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 5) x.GetPlot(PFMETvsTwo5NeutralHadron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 6) x.GetPlot(PFMETvsTwo6HadronicEnergyinHF)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 7) x.GetPlot(PFMETvsTwo7EMEnergyinHF)->Fill(PFMET_pt);


	  if(PFMET_pt > 300){
		Print();
		WeirdTwoPfEvent++;
		if(HighPtPhoton.size()!=0)	PrintPhoton();
  	   }
	}
	else if((Shower.size()==1 && AntiShower.size() == 1 && pfMono==1 && pfAmon==0) 
	     ||( Shower.size()==1 && AntiShower.size() == 1 && pfMono==0 && pfAmon==1)){
		OneCandOnePfEvent++;
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | MET should be large"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;
		recoMono=1;
		if(pfMono!=0){
	
			if(pfUnknown[0].pf_pdgId_ == 1) x.GetPlot(PFMETvsOneCandOne1ChargedHadron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 2) x.GetPlot(PFMETvsOneCandOne2Electron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 3) x.GetPlot(PFMETvsOneCandOne3Muon)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 4) x.GetPlot(PFMETvsOneCandOne4Photon)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 5) x.GetPlot(PFMETvsOneCandOne5NeutralHadron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 6) x.GetPlot(PFMETvsOneCandOne6HadronicEnergyinHF)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 7) x.GetPlot(PFMETvsOneCandOne7EMEnergyinHF)->Fill(PFMET_pt);

		}
		if(pfAmon!=0){

			if(pfAntiUnknown[0].pf_pdgId_ == 1) x.GetPlot(PFMETvsOneCandOne1ChargedHadron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 2) x.GetPlot(PFMETvsOneCandOne2Electron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 3) x.GetPlot(PFMETvsOneCandOne3Muon)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 4) x.GetPlot(PFMETvsOneCandOne4Photon)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 5) x.GetPlot(PFMETvsOneCandOne5NeutralHadron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 6) x.GetPlot(PFMETvsOneCandOne6HadronicEnergyinHF)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 7) x.GetPlot(PFMETvsOneCandOne7EMEnergyinHF)->Fill(PFMET_pt);
		}
	} 
	else if(Shower.size()==1 && AntiShower.size() == 0 && pfMono==1 && pfAmon==0){
		OnePfEvent++;
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | MET should be large"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;
			if(pfUnknown[0].pf_pdgId_ == 1) x.GetPlot(PFMETvsOnlyOne1ChargedHadron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 2) x.GetPlot(PFMETvsOnlyOne2Electron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 3) x.GetPlot(PFMETvsOnlyOne3Muon)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 4) x.GetPlot(PFMETvsOnlyOne4Photon)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 5) x.GetPlot(PFMETvsOnlyOne5NeutralHadron)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 6) x.GetPlot(PFMETvsOnlyOne6HadronicEnergyinHF)->Fill(PFMET_pt);
			if(pfUnknown[0].pf_pdgId_ == 7) x.GetPlot(PFMETvsOnlyOne7EMEnergyinHF)->Fill(PFMET_pt);
	}
	else if(Shower.size()==0 && AntiShower.size() == 1 && pfMono==0 && pfAmon==1){
		
		OnePfEvent++;
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | MET should be large"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;
			if(pfAntiUnknown[0].pf_pdgId_ == 1) x.GetPlot(PFMETvsOnlyOne1ChargedHadron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 2) x.GetPlot(PFMETvsOnlyOne2Electron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 3) x.GetPlot(PFMETvsOnlyOne3Muon)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 4) x.GetPlot(PFMETvsOnlyOne4Photon)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 5) x.GetPlot(PFMETvsOnlyOne5NeutralHadron)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 6) x.GetPlot(PFMETvsOnlyOne6HadronicEnergyinHF)->Fill(PFMET_pt);
			if(pfAntiUnknown[0].pf_pdgId_ == 7) x.GetPlot(PFMETvsOnlyOne7EMEnergyinHF)->Fill(PFMET_pt);
	} 
	//MET should be small case
	else if(Shower.size()==1 && AntiShower.size() == 1 && pfMono==0 && pfAmon==0){

		TwoCandEvent++;
		recoMono=2;
		x.GetPlot(PFMETvsTwoCand)->Fill(PFMET_pt);
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | MET should be small"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;

	  if(PFMET_pt > 300){
		Print();
		WeirdTwoCandEvent++;
		if(HighPtPhoton.size()!=0)	PrintPhoton();
  	   }
	}
	else if(Shower.size()==1 && AntiShower.size() == 0 && pfMono==0 && pfAmon==0){

		OneCandEvent++;
		recoMono=1;
		x.GetPlot(PFMETvsOnlyOneCand)->Fill(PFMET_pt);
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | MET should be small"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;
	  if(PFMET_pt > 300){
		Print();
		WeirdOneCandEvent++;
		if(HighPtPhoton.size()!=0)	PrintPhoton();
  	   }
	}
	else if(Shower.size()==0 && AntiShower.size() == 1 && pfMono==0 && pfAmon==0){

		OneCandEvent++;
		recoMono=1;
		x.GetPlot(PFMETvsOnlyOneCand)->Fill(PFMET_pt);
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | MET should be small"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;
	  if(PFMET_pt > 300){
		Print();
		WeirdOneCandEvent++;
		if(HighPtPhoton.size()!=0)	PrintPhoton();
  	   }
	}
	else if(Shower.size()==0 && AntiShower.size()==0){
		
		ZeroCandEvent++;
		x.GetPlot(PFMETvsZeroCand)->Fill(PFMET_pt);
		cout<<setw(5)<<" ---------------------"<<endl;		
		cout<<setw(5)<<" | No monopole reconstructed"<<endl;
		cout<<setw(5)<<" | MET should be small"<<endl;
		cout<<setw(5)<<" | MET "<<PFMET_pt<<endl;
		cout<<setw(5)<<" | phi "<<PFMET_phi<<endl;
	  if(PFMET_pt > 300){
		WeirdZeroCandEvent++;
		Print();
		if(HighPtPhoton.size()!=0)	PrintPhoton();
  	   }
	}
	else{
		Bug++;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"THERE IS BUG*************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<"Shower.size     "<<Shower.size()<<endl;
		cout<<"AntiShower.size "<<AntiShower.size()<<endl;
		cout<<"pfMono          "<<pfMono<<endl;
		cout<<"pfAmon          "<<pfAmon<<endl;
	}
	
	//debug
	//To see weird case MET phi
	
//}//PFMET_pt>200		
//        	cout<<"event tag"<<CutFlowCand_TRG[0].event_<<endl;
			cout<<"============================================="<<endl;		
}
  void WritePlots(TFile *oFile ){
	oFile->cd(trigName_.c_str());
	
	Plot[0].WritePlot();
	
  }
 void Print(){
	
	for(int i=0;i<Shower.size();i++){
                double m_deltaR=0;
                m_deltaR = sqrt(pow(Shower[i].eta_-Shower[0].mono_eta_,2)+
                                pow(Shower[i].phi_-Shower[0].mono_phi_,2));
		if(m_deltaR<0.15){
			cout<<i+1<<" Mono +++++++++++++++++++++++++++++++++++"<<endl;
			cout<<"     eta "<<setprecision(5)<<Shower[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<Shower[i].phi_<<endl;
			cout<<"     e55 "<<setprecision(5)<<Shower[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<Shower[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<Shower[i].f51_<<endl;
			cout<<"     f15 "<<setprecision(5)<<Shower[i].f15_<<endl;
			cout<<"   Cross "<<setprecision(5)<<Shower[i].SwissCross_<<endl;
			cout<<"     MET "<<setprecision(5)<<Shower[i].MET_pt_<<endl;
			cout<<" MET phi "<<setprecision(5)<<Shower[i].MET_phi_<<endl;
		
              		cout<<"          candidate           monoGen         antimonoGen"<<endl;
                	cout<<"eta      "<<setprecision(5)<<Shower[i].eta_<<setw(20)<<Shower[i].mono_eta_
                        	         <<setw(20)<<Shower[i].amon_eta_<<endl;
                	cout<<"phi      "<<setprecision(5)<<Shower[i].phi_<<setw(20)<<Shower[i].mono_phi_
                       		         <<setw(20)<<Shower[i].amon_phi_<<endl;
                	cout<<"m deltaR "<<m_deltaR<<endl;
		}
	}
	for(int i=0;i<AntiShower.size();i++){
                double am_deltaR=0;
                am_deltaR= sqrt(pow(AntiShower[i].eta_-AntiShower[0].amon_eta_,2)+
                                pow(AntiShower[i].phi_-AntiShower[0].amon_phi_,2));
		if(am_deltaR<0.15){
			cout<<i+1<<" Amon ++++++++++++++++++++++++++++++++++++"<<endl;
			cout<<"     eta "<<setprecision(5)<<AntiShower[i].eta_<<endl;
			cout<<"     phi "<<setprecision(5)<<AntiShower[i].phi_<<endl;
			cout<<"     e55 "<<setprecision(5)<<AntiShower[i].e55_<<endl;
			cout<<" dEdxSig "<<setprecision(5)<<AntiShower[i].dEdXSig_<<endl;
			cout<<"     f51 "<<setprecision(5)<<AntiShower[i].f51_<<endl;
			cout<<"     f15 "<<setprecision(5)<<AntiShower[i].f15_<<endl;
			cout<<"   Cross "<<setprecision(5)<<AntiShower[i].SwissCross_<<endl;
			cout<<"     MET "<<setprecision(5)<<AntiShower[i].MET_pt_<<endl;
			cout<<" MET phi "<<setprecision(5)<<AntiShower[i].MET_phi_<<endl;
		
              		cout<<"          candidate           monoGen         antimonoGen"<<endl;
                	cout<<"eta      "<<setprecision(5)<<AntiShower[i].eta_<<setw(20)<<AntiShower[i].mono_eta_
                        	         <<setw(20)<<AntiShower[i].amon_eta_<<endl;
                	cout<<"phi      "<<setprecision(5)<<AntiShower[i].phi_<<setw(20)<<AntiShower[i].mono_phi_
                       		         <<setw(20)<<AntiShower[i].amon_phi_<<endl;
                	cout<<"a deltaR "<<am_deltaR<<endl;
		}
	}
  }
  void PrintPhoton(){
       for(int j=0;j<HighPtPhoton.size();j++){
                cout<<j<<endl;
                cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
                cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
                cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
        }	
  }
  void SignalEff(const string trName,unsigned TotalEvents)
  {
        //signal efficiency = no. of events after all selection cuts/all events
	
	cout<<endl;
	cout<<trName<<" ============================="<<endl;
	cout<<"Statistical Result for Mass 1000GeV "<<endl;
	cout<<" _____________________________"<<endl;
	cout<<"Large MET case --------------"<<endl;
	cout<<" | 1 PFlike 1 Reco Mono "<<OneCandOnePfEvent<<setprecision(4)<<setw(10)<<(double)OneCandOnePfEvent/TotalEvents<<endl;
	cout<<" | 1 PFlike             "<<OnePfEvent<<setprecision(4)<<setw(10)<<(double)OnePfEvent/TotalEvents<<endl;
	cout<<" ____________________________"<<endl;
	cout<<"Small MET case --------------"<<endl;
	cout<<" | 2 PFlike 	       "<<TwoPfEvent<<setprecision(4)<<setw(10)<<(double)TwoPfEvent/TotalEvents<<endl;
	cout<<" | 	   2 Reco Mono "<<TwoCandEvent<<setprecision(4)<<setw(10)<<(double)TwoCandEvent/TotalEvents<<endl;
	cout<<" | 	   1 Reco Mono "<<OneCandEvent<<setprecision(4)<<setw(10)<<(double)OneCandEvent/TotalEvents<< endl;
	cout<<" | 0 Monopole in event "<<ZeroCandEvent<<setprecision(4)<<setw(10)<<(double)ZeroCandEvent/TotalEvents<<endl;
	cout<<" ____________________________"<<endl;
	cout<<endl;
	cout<<"Weird Small MET case --------------"<<endl;
	cout<<" | 2 PFlike 	       "<<WeirdTwoPfEvent<<setprecision(4)<<setw(10)<<(double)WeirdTwoPfEvent/TotalEvents<<endl;
	cout<<" | 	   2 Reco Mono "<<WeirdTwoCandEvent<<setprecision(4)<<setw(10)<<(double)WeirdTwoCandEvent/TotalEvents<<endl;
	cout<<" | 	   1 Reco Mono "<<WeirdOneCandEvent<<setprecision(4)<<setw(10)<<(double)WeirdOneCandEvent/TotalEvents<< endl;
	cout<<" | 0 Monopole in event "<<WeirdZeroCandEvent<<setprecision(4)<<setw(10)<<(double)WeirdZeroCandEvent/TotalEvents<<endl;
	cout<<" ____________________________"<<endl;
	cout<<"Bug "<<Bug<<endl;
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
  static const double PFMETCut_ ;
  static const double photonCut_ ;
  static const double electronCut_ ;
  static const double particleflowCut_ ;
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
  bool evalPFMET(MonoCandidate &mc) { return mc.MET_pt_ > PFMETCut_ ; }
  bool evalF15(MonoCandidate &mc) { return mc.f15_ > f15Cut_ ; }
  bool evalSwissCross(MonoCandidate &mc) { return mc.SwissCross_ > SwissCrossCut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
  bool evalElectron(Electron &mc){ return mc.ele_pt_ > electronCut_; }
  bool evalParticleflow(Particleflow &mc){ return mc.pf_E_ > particleflowCut_; }
	
  ofstream BugInfo;//("BugInfo.txt");

  vector<MonoCandidate> CutFlowCand_TRG;
  vector<MonoCandidate> CutFlowCand_Qual;
  vector<MonoCandidate> CutFlowCand_Energy;
  vector<MonoCandidate> CutFlowCand_F51;
  vector<MonoCandidate> CutFlowCand_Dedx;
	
  vector<Photon> HighPtPhoton;
  vector<Electron> HighPtElectron;
  vector<Particleflow> HighPtParticleflow;
  vector<Particleflow> pfPhotonlike;
  vector<Particleflow> pfElectronlike;
  vector<Particleflow> pfUnknown;
  vector<Particleflow> pfAntiUnknown;

  vector<MonoCandidate> Shower;
  vector<MonoCandidate> AntiShower;
  vector<MonoCandidate> EcalShower;

  vector<double> pfParticle_Mono_deltaR;
  vector<double> pfParticle_Amon_deltaR ;
 
  int Bug=0;
  int pfphotonlike=0; 
  int pfelectronlike=0; 
  int pfMono=0;
  int pfAmon=0;
  int photonLike=0;
  int EgammaLike=0;
  int electronLike=0;
  int recoMono=0;
  //large
	int OneCandOnePfEvent = 0;
	int OnePfEvent = 0;
  //small 
	int TwoPfEvent =0; 
	int TwoCandEvent = 0;	
	int OneCandEvent = 0;
	int ZeroCandEvent = 0;		

	int WeirdTwoPfEvent =0; 
	int WeirdTwoCandEvent = 0;	
	int WeirdOneCandEvent = 0;
	int WeirdZeroCandEvent = 0;		

  int TotalParticleflowLike=0;
  int EMShowerButSpikelikeinBarrel=0;
  int EMShowerButSpikelikeinEE=0;


  int TrackerMono=0;
  int EcalShowerNum=0;
 

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
  const double MonoCuts::e55Cut_ = 100;
//  const double MonoCuts::f51Cut_ = 0.6;
  const double MonoCuts::f15Cut_ = 0.95;
  const double MonoCuts::SwissCrossCut_ = 0.95;
  const double MonoCuts::f51Cut_ = 0.85;
  const double MonoCuts::PFMETCut_ = 200;//we only look to high MET case
  const double MonoCuts::photonCut_ = 100;
  const double MonoCuts::electronCut_ = 0;
  const double MonoCuts::particleflowCut_ = 100;//the lowest require for PF cut
void MonoAnalyzerPFMET_PF_0805()
{
	cout<<"hello"<<endl;
	TFile *oFile = new TFile("Monopole_PFMET_PF_1000.root","recreate");
	cout<<"new file line pass"<<endl;

	TFile *fin = new TFile("MonoNtuple2018_MC_1000_0804.root");
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
        vector<double> * pf_eta = 0;
        vector<double> * pf_phi = 0;
        vector<double> * pf_pt = 0;
        vector<double> * pf_E = 0;
        vector<int> * pf_pdgId = 0;
        unsigned nParticleflow;
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

        tree->SetBranchAddress("pf_N",&nParticleflow);
        tree->SetBranchAddress("pf_eta",&pf_eta);
        tree->SetBranchAddress("pf_phi",&pf_phi);
        tree->SetBranchAddress("pf_pt",&pf_pt);
        tree->SetBranchAddress("pf_E",&pf_E);
        tree->SetBranchAddress("pf_pdgId",&pf_pdgId);

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
        vector<Particleflow> particleflow(0);
	
	ofstream BugInfo ("BugInformation.txt");
	if(BugInfo.is_open())
	{
		BugInfo<<"test";
	}
        for(unsigned ev=0; ev<NEvents;ev++){
	//	cout<<"Enter event loop"<<endl;
//		cout<<"event number "<<NEvents<<endl;
             tree->GetEntry(ev);
		
		if(nCandidates>Et.size()) Et.resize(nCandidates);
		if(nCandidates>cand.size()) cand.resize(nCandidates);
                if(nElectron>electron.size()) electron.resize(nElectron);
                if(nParticleflow>particleflow.size()) particleflow.resize(nParticleflow);
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
//		Et[i],
	        (*e55)[i],
		(*SwissCross)[i],
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
 
               if(nParticleflow!=0){
                for(unsigned j=0;j<nParticleflow;j++){
                        particleflow[j] = Particleflow(
                        (*pf_eta)[j],
                        (*pf_phi)[j],
                        (*pf_pt)[j],
			(*pf_E)[j],
                        (*pf_pdgId)[j]
                        );
                        }
                }
			noTrgAnalysis.doAnalysis(cand,photon,electron,particleflow,nCandidates,nPhoton,
						nElectron,nParticleflow,true,ev,MET_pt,MET_phi);
				
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
	//BugInfo.close();
	cout<<"end of the code"<<endl;
}
