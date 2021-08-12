////////////////////////////////////
//
//
//      File name: MonoDrawKinematics.cc
//      author: Shih Lin
//      Purpose: Draw Kinematics variables to understand gen Mono more
//      Date: 2021.7.28
//
///////////////////////////////////////////////////

#include "iostream"
#include "TAttMarker.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "math.h"
#include <algorithm>

using namespace std;
void MonoDrawKinematics()
{
        TFile *f1 = new TFile("MonoNtuple2018_MC_1000_0722.root");
        TTree *tree1 = (TTree*)f1->Get("monopoles");
        TFile *f2 = new TFile("MonoNtuple2018_MC_2000_0722.root");
        TTree *tree2 = (TTree*)f2->Get("monopoles");
	TFile *f3 = new TFile("MonoNtuple2018_MC_3000_0722.root");
        TTree *tree3 = (TTree*)f3->Get("monopoles");
	TFile *f4 = new TFile("MonoNtuple2018_MC_4000_0722.root");
        TTree *tree4 = (TTree*)f4->Get("monopoles");

	double Mass1;
	double Mass2;
	double Mass3;
	double Mass4;

	double px1;
	double py1;
	double px2;
        double py2;
	double px3;
	double py3;
	double px4;
	double py4;

	double eta1;
        double eta2;
        double eta3;
        double eta4;

	double Momentum1;
	double Momentum2;
	double Momentum3;
	double Momentum4;

	double MomentumTransverse1;
	double MomentumTransverse2;
	double MomentumTransverse3;
	double MomentumTransverse4;
	
	double EnergyTransverse1;
	double EnergyTransverse2;
	double EnergyTransverse3;
	double EnergyTransverse4;
	

	double Energy1;
	double Energy2;
	double Energy3;
	double Energy4;

	double MassTransverse1;
	double MassTransverse2;
	double MassTransverse3;
	double MassTransverse4;

	tree1->SetBranchAddress("mono_px",&px1);
        tree1->SetBranchAddress("mono_py",&py1);
        tree1->SetBranchAddress("mono_eta",&eta1);
        tree1->SetBranchAddress("mono_p",&Momentum1);
        tree1->SetBranchAddress("mono_pt",&MomentumTransverse1);
        tree1->SetBranchAddress("mono_Et",&EnergyTransverse1);
        tree1->SetBranchAddress("mono_E",&Energy1);
        tree1->SetBranchAddress("mono_mt",&MassTransverse1);
        tree1->SetBranchAddress("mono_m",&Mass1);

	tree2->SetBranchAddress("mono_px",&px2);
        tree2->SetBranchAddress("mono_py",&py2);
        tree2->SetBranchAddress("mono_eta",&eta2);
        tree2->SetBranchAddress("mono_p",&Momentum2);
        tree2->SetBranchAddress("mono_pt",&MomentumTransverse2);
        tree2->SetBranchAddress("mono_Et",&EnergyTransverse2);
        tree2->SetBranchAddress("mono_E",&Energy2);
        tree2->SetBranchAddress("mono_mt",&MassTransverse2);
        tree2->SetBranchAddress("mono_m",&Mass2);
	
	tree3->SetBranchAddress("mono_px",&px3);
        tree3->SetBranchAddress("mono_py",&py3);
        tree3->SetBranchAddress("mono_eta",&eta3);
        tree3->SetBranchAddress("mono_p",&Momentum3);
        tree3->SetBranchAddress("mono_pt",&MomentumTransverse3);
        tree3->SetBranchAddress("mono_Et",&EnergyTransverse3);
        tree3->SetBranchAddress("mono_E",&Energy3);
        tree3->SetBranchAddress("mono_mt",&MassTransverse3);
        tree3->SetBranchAddress("mono_m",&Mass3);
	
	tree4->SetBranchAddress("mono_px",&px4);
        tree4->SetBranchAddress("mono_py",&py4);
        tree4->SetBranchAddress("mono_eta",&eta4);
        tree4->SetBranchAddress("mono_p",&Momentum4);
        tree4->SetBranchAddress("mono_pt",&MomentumTransverse4);
        tree4->SetBranchAddress("mono_Et",&EnergyTransverse4);
        tree4->SetBranchAddress("mono_E",&Energy4);
        tree4->SetBranchAddress("mono_mt",&MassTransverse4);
        tree4->SetBranchAddress("mono_m",&Mass4);

	TCanvas *canvas = new TCanvas("canvas","canvas",800,600);

	TH1D *pT1 = new TH1D("pT1","",100,-1,4000);
        TH1D *pT2 = new TH1D("pT2","",100,-1,4000);
	TH1D *pT3 = new TH1D("pT3","",100,-1,4000);
	TH1D *pT4 = new TH1D("pT4","",100,-1,4000);

	TH1D *Eta1 = new TH1D("Eta1","",100,-4,4);
        TH1D *Eta2 = new TH1D("Eta2","",100,-4,4);
        TH1D *Eta3 = new TH1D("Eta3","",100,-4,4);
        TH1D *Eta4 = new TH1D("Eta4","",100,-4,4);

	TH1D *P1 = new TH1D("P1","",100,-1,4000);
        TH1D *P2 = new TH1D("P2","",100,-1,4000);
        TH1D *P3 = new TH1D("P3","",100,-1,4000);
        TH1D *P4 = new TH1D("P4","",100,-1,4000);
	
	TH1D *Et1 = new TH1D("Et1","",100,-1,5000);
        TH1D *Et2 = new TH1D("Et2","",100,-1,5000);
        TH1D *Et3 = new TH1D("Et3","",100,-1,5000);
        TH1D *Et4 = new TH1D("Et4","",100,-1,5000);

	TH1D *E1 = new TH1D("E1","",100,990,5600);
        TH1D *E2 = new TH1D("E2","",100,990,5600);
        TH1D *E3 = new TH1D("E3","",100,990,5600);
        TH1D *E4 = new TH1D("E4","",100,990,5600);
	
	TH1D *Mt1 = new TH1D("Mt1","",100,990,4600);
        TH1D *Mt2 = new TH1D("Mt2","",100,990,4600);
        TH1D *Mt3 = new TH1D("Mt3","",100,990,4600);
        TH1D *Mt4 = new TH1D("Mt4","",100,990,4600);





	cout<<"total event "<<tree1->GetEntries()<<endl;
	for(unsigned ev=0;ev<15000;ev++){
		cout<<"ev "<<ev<<endl;
		tree1->GetEntry(ev);
		pT1->Fill(MomentumTransverse1);
  		Eta1->Fill(eta1);
  		P1->Fill(Momentum1);
  		Et1->Fill(EnergyTransverse1);
  		E1->Fill(Energy1);
  		Mt1->Fill(MassTransverse1);
		
				

                tree2->GetEntry(ev);
		pT2->Fill(MomentumTransverse2);
  		Eta2->Fill(eta2);
  		P2->Fill(Momentum2);
  		Et2->Fill(EnergyTransverse2);
  		E2->Fill(Energy2);
  		Mt2->Fill(MassTransverse2);
		
		tree3->GetEntry(ev);
		pT3->Fill(MomentumTransverse3);
  		Eta3->Fill(eta3);
  		P3->Fill(Momentum3);
  		Et3->Fill(EnergyTransverse3);
  		E3->Fill(Energy3);
  		Mt3->Fill(MassTransverse3);
		
		tree4->GetEntry(ev);
		pT4->Fill(MomentumTransverse4);
  		Eta4->Fill(eta4);
  		P4->Fill(Momentum4);
  		Et4->Fill(EnergyTransverse4);
  		E4->Fill(Energy4);
  		Mt4->Fill(MassTransverse4);
	}
        gStyle->SetOptStat(0);

	cout<<"pT----------------------"<<endl;
	cout<<"Mean 1 : "<<pT1->GetMean()<<endl;
	cout<<"Mean 2 : "<<pT2->GetMean()<<endl;
	cout<<"Mean 3 : "<<pT3->GetMean()<<endl;
	cout<<"Mean 4 : "<<pT4->GetMean()<<endl;
	pT1->Draw();
	pT1->SetXTitle("Transverse momentum pT (GeV)");
	pT1->SetYTitle("Entries");
	pT1->SetLineColorAlpha(46,1.5);
	pT2->SetLineColorAlpha(9,1.5);
	pT3->SetLineColorAlpha(8,1.5);
	pT4->SetLineColorAlpha(41,1.5);
	pT2->Draw("same");
	pT3->Draw("same");
	pT4->Draw("same");
	mgr::SetSinglePad(canvas);
        mgr::SetAxis(pT1);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("15000 Events 2018 13TeV");
        TLegend *leg2 = new TLegend(0.55,0.65,0.8,0.85);
        leg2->SetBorderSize( 0 );
        leg2->SetFillColorAlpha( 0, 0 );
        leg2->AddEntry(pT1,"1000 GeV ","f");
        leg2->AddEntry(pT2,"2000 GeV ","f");
	leg2->AddEntry(pT3,"3000 GeV ","f");
	leg2->AddEntry(pT4,"4000 GeV ","f");
        leg2->Draw();
        canvas->SaveAs("GenMonoPt.pdf");

	cout<<"P----------------------"<<endl;
	cout<<"Mean 1 : "<<P1->GetMean()<<endl;
	cout<<"Mean 2 : "<<P2->GetMean()<<endl;
	cout<<"Mean 3 : "<<P3->GetMean()<<endl;
	cout<<"Mean 4 : "<<P4->GetMean()<<endl;
	P1->Draw();
	P1->SetXTitle("Momentum (GeV)");
	P1->SetYTitle("Entries");
	P1->SetLineColorAlpha(46,1.5);
	P2->SetLineColorAlpha(9,1.5);
	P3->SetLineColorAlpha(8,1.5);
	P4->SetLineColorAlpha(41,1.5);
	P2->Draw("same");
	P3->Draw("same");
	P4->Draw("same");
	mgr::SetSinglePad(canvas);
        mgr::SetAxis(P1);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("15000 Events 2018 13TeV");
        TLegend *leg1 = new TLegend(0.55,0.65,0.8,0.85);
        leg1->SetBorderSize( 0 );
        leg1->SetFillColorAlpha( 0, 0 );
        leg1->AddEntry(P1,"1000 GeV ","f");
        leg1->AddEntry(P2,"2000 GeV ","f");
	leg1->AddEntry(P3,"3000 GeV ","f");
	leg1->AddEntry(P4,"4000 GeV ","f");
        leg1->Draw();
        canvas->SaveAs("GenMonoMomentum.pdf");

	cout<<"Et----------------------"<<endl;
	cout<<"Mean 1 : "<<Et1->GetMean()<<endl;
	cout<<"Mean 2 : "<<Et2->GetMean()<<endl;
	cout<<"Mean 3 : "<<Et3->GetMean()<<endl;
	cout<<"Mean 4 : "<<Et4->GetMean()<<endl;
	Et1->Draw();
	Et1->SetXTitle("Transverse Energy (GeV)");
	Et1->SetYTitle("Entries");
	Et1->SetLineColorAlpha(46,1.5);
	Et2->SetLineColorAlpha(9,1.5);
	Et3->SetLineColorAlpha(8,1.5);
	Et4->SetLineColorAlpha(41,1.5);
	Et2->Draw("same");
	Et3->Draw("same");
	Et4->Draw("same");
	mgr::SetSinglePad(canvas);
        mgr::SetAxis(Et1);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("15000 Events 2018 13TeV");
        TLegend *leg4 = new TLegend(0.55,0.65,0.8,0.85);
        leg4->SetBorderSize( 0 );
        leg4->SetFillColorAlpha( 0, 0 );
        leg4->AddEntry(Et1,"1000 GeV ","f");
        leg4->AddEntry(Et2,"2000 GeV ","f");
	leg4->AddEntry(Et3,"3000 GeV ","f");
	leg4->AddEntry(Et4,"4000 GeV ","f");
        leg4->Draw();
        canvas->SaveAs("GenMonoEnergyTransverse.pdf");

	cout<<"E----------------------"<<endl;
	cout<<"Mean 1 : "<<E1->GetMean()<<endl;
	cout<<"Mean 2 : "<<E2->GetMean()<<endl;
	cout<<"Mean 3 : "<<E3->GetMean()<<endl;
	cout<<"Mean 4 : "<<E4->GetMean()<<endl;
	E1->Draw();
	E1->SetXTitle("Energy (GeV)");
	E1->SetYTitle("Entries");
	E1->SetLineColorAlpha(46,1.5);
	E2->SetLineColorAlpha(9,1.5);
	E3->SetLineColorAlpha(8,1.5);
	E4->SetLineColorAlpha(41,1.5);
	E2->Draw("same");
	E3->Draw("same");
	E4->Draw("same");
	mgr::SetSinglePad(canvas);
        mgr::SetAxis(E1);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("15000 Events 2018 13TeV");
        TLegend *leg5 = new TLegend(0.55,0.65,0.8,0.85);
        leg5->SetBorderSize( 0 );
        leg5->SetFillColorAlpha( 0, 0 );
        leg5->AddEntry(E1,"1000 GeV ","f");
        leg5->AddEntry(E2,"2000 GeV ","f");
	leg5->AddEntry(E3,"3000 GeV ","f");
	leg5->AddEntry(E4,"4000 GeV ","f");
        leg5->Draw();
        canvas->SaveAs("GenMonoEnergy.pdf");

	cout<<"Mt----------------------"<<endl;
	cout<<"Mean 1 : "<<Mt1->GetMean()<<endl;
	cout<<"Mean 2 : "<<Mt2->GetMean()<<endl;
	cout<<"Mean 3 : "<<Mt3->GetMean()<<endl;
	cout<<"Mean 4 : "<<Mt4->GetMean()<<endl;
	Mt1->Draw();
	Mt1->SetXTitle("Transverse Mass (GeV)");
	Mt1->SetYTitle("Entries");
	Mt1->SetLineColorAlpha(46,1.5);
	Mt2->SetLineColorAlpha(9,1.5);
	Mt3->SetLineColorAlpha(8,1.5);
	Mt4->SetLineColorAlpha(41,1.5);
	Mt2->Draw("same");
	Mt3->Draw("same");
	Mt4->Draw("same");
	mgr::SetSinglePad(canvas);
        mgr::SetAxis(Mt1);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("15000 Events 2018 13TeV");
        TLegend *leg6 = new TLegend(0.55,0.65,0.8,0.85);
        leg6->SetBorderSize( 0 );
        leg6->SetFillColorAlpha( 0, 0 );
        leg6->AddEntry(Mt1,"1000 GeV ","f");
        leg6->AddEntry(Mt2,"2000 GeV ","f");
	leg6->AddEntry(Mt3,"3000 GeV ","f");
	leg6->AddEntry(Mt4,"4000 GeV ","f");
        leg6->Draw();
        canvas->SaveAs("GenMonoTransverseMass.pdf");

	Eta4->Draw();
        Eta4->SetXTitle("eta");
        Eta4->SetYTitle("Entries");
        Eta1->SetLineColorAlpha(46,1.5);
        Eta2->SetLineColorAlpha(9,1.5);
        Eta3->SetLineColorAlpha(8,1.5);
        Eta4->SetLineColorAlpha(41,1.5);
        Eta1->Draw("same");
        Eta2->Draw("same");
        Eta3->Draw("same");
        mgr::SetSinglePad(canvas);
	mgr::SetAxis(Eta4);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("15000 Events 2018 13TeV");
	TLegend *leg3 = new TLegend(0.70,0.60,0.90,0.80);
        leg3->SetBorderSize( 0 );
        leg3->SetFillColorAlpha( 0, 0 );
        leg3->AddEntry(Eta1,"1000 GeV ","f");
        leg3->AddEntry(Eta2,"2000 GeV ","f");
        leg3->AddEntry(Eta3,"3000 GeV ","f");
        leg3->AddEntry(Eta4,"4000 GeV ","f");
        leg3->Draw();
        canvas->SaveAs("GenMonoEta.pdf");



}	
	


