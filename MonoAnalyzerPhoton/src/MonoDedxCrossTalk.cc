//////////////////////////////////
//
//	File name: MonoDedxCrosstalk.cc
//	author: Shih Lin
//
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
#include "/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MonoAnalyzerPhoton/interface/Candidate.h"
#include "/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MonoAnalyzerPhoton/interface/PlotSet.h"
#include "/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MonoAnalyzerPhoton/interface/MonoCuts.h"
  void MonoCuts::doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev)
  {
	Clear();	
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
	if(TRG && dEdXCut) CrossTalk_QualDedx.push_back(cands);
	
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

	//Only use for CrossTalk 
	sort(CrossTalk_QualDedx.begin(),CrossTalk_QualDedx.begin()+CrossTalk_QualDedx.size());
	if(CrossTalk_QualDedx.size()>0){
		crosstalk_count++;
	}
	
	sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());

		
	//Cut all
	if(CutFlowCand_Dedx.size()>0)
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
		//PrintInfo();

	   for(int i=0; i<CutFlowCand_Dedx.size();i++){

		double m_deltaR=0;
		double am_deltaR=0;
		m_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].mono_eta_,2)+
				pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].mono_phi_,2));
		am_deltaR= sqrt(pow(CutFlowCand_Dedx[i].eta_-CutFlowCand_Dedx[0].amon_eta_,2)+
                                pow(CutFlowCand_Dedx[i].phi_-CutFlowCand_Dedx[0].amon_phi_,2));

		if(m_deltaR<0.15) 	Reco++;
		else if(am_deltaR<0.15)	Reco++;
	
		if(abs(CutFlowCand_Dedx[i].eta_ ) < 1.48){
		   
		   EBarrel++;

  		   if(CutFlowCand_Dedx[i].Cross_ > 0.95){
			SpikeLike++;
			x.GetPlot(Spike)->Fill(CutFlowCand_Dedx[i].eta_);
		   }
		}
	   }
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
  
		//to see monopole is photon-like or spike-like
		//define: photon-monopole angle<0.5 is photon-like
		//for(int j=0;j<HighPtPhoton.size();j++){
		//	double photonMono_deltaR = 0;
		//	photonMono_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-HighPtPhoton[j].pho_eta_,2)+
                //	               	         pow(CutFlowCand_Dedx[i].phi_-HighPtPhoton[j].pho_phi_,2));
		//	if(photonMono_deltaR<0.15){
		//		cout<<"deltaR "<<photonMono_deltaR<<endl;
		//		cout<<"photon-like monopole"<<endl;
		//		photonLike++;// otherwise not photon-like monopole
		//	}
		//}
		cout<<endl;
		cout<<"=================================="<<endl;
  }
  void SignalEff(const string trName)
  {
        //signal efficiency = no. of events after all selection cuts/all events
	cout<<trName<<" ================================="<<endl;
        cout<<"        TRG "<<count<<endl;
        cout<<"QualityCuts "<<Qual_count<<endl;
        cout<<"       ECut "<<E_count<<endl;
        cout<<"     F51Cut "<<f51_count<<endl;
        cout<<" dEdXSigCut "<<dEdX_count<<endl;
	cout<<"  CrossTalk "<<crosstalk_count<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"Signal efficiency = "<<(double)dEdX_count/(double)TotalEvents<<endl;
        cout<<endl;
	cout<<"Relative effciency count"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"     No TRG "<<NoTRG<<endl;
        cout<<" No Quality "<<NoQual<<endl;
        cout<<"    No ECut "<< NoE <<endl;
        cout<<"  No F51Cut "<<NoF51<<endl;
        cout<<" No dEdXSig "<<NodEdXCut<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"Relative efficiency"<<endl;
        cout<<"------------------------------------------"<<endl;
        cout<<"        TRG "<<(double)dEdX_count/(double)NoTRG<<endl;
        cout<<"QualityCuts "<<(double)dEdX_count/(double)NoQual<<endl;
        cout<<"       ECut "<<(double)dEdX_count/(double)NoE<<endl;
        cout<<"     F51Cut "<<(double)dEdX_count/(double)NoF51<<endl;
        cout<<" dEdXSigCut "<<(double)dEdX_count/(double)NodEdXCut<<endl;

        cout<<endl;
	cout<<"Total reconstructed event "<<Reco<<endl;
	cout<<"EBarrel "<<EBarrel<<endl;
	cout<<"SpikeLike in EB "<<SpikeLike<<endl;
	  
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
  bool evalQuality(MonoCandidate &mc) { return TMath::Abs(mc.xyp0_) < xyp0Cut_&& TMath::Abs(mc.xyp2_) > xyp2Cut_ 
			&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_ 
			&& TMath::Abs(mc.rzp2_) < rzp2Cut_ && TMath::Abs(mc.rzp1_)<rzp1Cut_ && TMath::Abs(mc.rzp0_) < rzp0Cut_;  }
  bool evalE(MonoCandidate &mc) { return mc.e55_ > e55Cut_; }
  bool evalF51(MonoCandidate &mc) { return mc.f51_ > f51Cut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
  int TotalEvents = 16389;
 
  vector<MonoCandidate> CutFlowCand_TRG;
  vector<MonoCandidate> CutFlowCand_Qual;
  vector<MonoCandidate> CutFlowCand_Energy;
  vector<MonoCandidate> CutFlowCand_F51;
  vector<MonoCandidate> CutFlowCand_Dedx;
  vector<MonoCandidate> CrossTalk_QualDedx;

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
  int crosstalk_count=0; 
  int SpikeLike=0;
  int EBarrel=0; 
 
  int Reco=0;
  int photonLike=0;

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

void MonoDedxCrossTalk()
{
	TFile *oFile = new TFile("MonoPhotonAnalysis_2018_1000.root","recreate");
//	TFile *oFile = new TFile("DedxCrossTalk_X0down_Analysis_2018_1000.root","recreate");
//	TFile *fin = new TFile("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/Systematic/DedxCrossTalk/X0/Down/1000/DedxCrossTalk_X0down_2018_1000.root");
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
	
	
	MonoCuts noTrgAnalysis("NO TRG",oFile);
	MonoCuts TrgAnalysis("HLT_Photon200",oFile);

	vector<MonoCandidate> cand(10);	
	vector<Photon> photon(0);
	cout<<"total evt  "<<NEvents<<endl;	
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
	//TrgAnalysis.SignalEff("HLT_Photon200");

	noTrgAnalysis.WritePlots(oFile);
	noTrgAnalysis.SignalEff("NO TRG");
	oFile->Close();	
	cout<<"end of the code"<<endl;
}
