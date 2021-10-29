////////////////////////////////////
//
//
//      File name: MonoDrawEtaPt.cc
//      author: Shih Lin
//      Content: Compare the eta between different mass
//      	(Gen Monopole)
//               
//      Output: /Plots/2016MC1TeV/*.pdf
//      Date: 2021.6.4
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
void MonoDrawEtaPt()
{
        TFile *f1 = new TFile("MonoNtuple2018_MC_1000.root");
        TTree *tree1 = (TTree*)f1->Get("monopoles");
        TFile *f2 = new TFile("MonoNtuple2018_MC_2000.root");
        TTree *tree2 = (TTree*)f2->Get("monopoles");
	TFile *f3 = new TFile("MonoNtuple2018_MC_3000.root");
        TTree *tree3 = (TTree*)f3->Get("monopoles");
	TFile *f4 = new TFile("MonoNtuple2018_MC_4000.root");
        TTree *tree4 = (TTree*)f4->Get("monopoles");

	

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

	tree1->SetBranchAddress("mono_px",&px1);
        tree1->SetBranchAddress("mono_py",&py1);
        tree1->SetBranchAddress("mono_eta",&eta1);

        tree2->SetBranchAddress("mono_px",&px2);
        tree2->SetBranchAddress("mono_py",&py2);
        tree2->SetBranchAddress("mono_eta",&eta2);

	tree3->SetBranchAddress("mono_px",&px3);
        tree3->SetBranchAddress("mono_py",&py3);
        tree3->SetBranchAddress("mono_eta",&eta3);

	tree4->SetBranchAddress("mono_px",&px4);
        tree4->SetBranchAddress("mono_py",&py4);
        tree4->SetBranchAddress("mono_eta",&eta4);
	
	TCanvas *canvas = new TCanvas("canvas","canvas",800,600);

	TH1D *h_energy = new TH1D("h_energy","",100,-1,4000);
        TH1D *h_energy2 = new TH1D("h_energy2","",100,-1,4000);
	TH1D *h_energy3 = new TH1D("h_energy3","",100,-1,4000);
	TH1D *h_energy4 = new TH1D("h_energy4","",100,-1,4000);
	TH1D *h_eta1 = new TH1D("h_eta1","",100,-4,4);
        TH1D *h_eta2 = new TH1D("h_eta2","",100,-4,4);
        TH1D *h_eta3 = new TH1D("h_eta3","",100,-4,4);
        TH1D *h_eta4 = new TH1D("h_eta4","",100,-4,4);
	TProfile *test = new TProfile("test","",4,0,4000);

	cout<<"total event "<<tree1->GetEntries()<<endl;
	for(unsigned ev=0;ev<10000;ev++){
		cout<<"ev "<<ev<<endl;
		tree1->GetEntry(ev);
		h_energy->Fill(sqrt(px1*px1+py1*py1));
  		h_eta1->Fill(eta1);

                tree2->GetEntry(ev);
                h_energy2->Fill(sqrt(px2*px2+py2*py2));
		h_eta2->Fill(eta2);
		
		tree3->GetEntry(ev);
                h_energy3->Fill(sqrt(px3*px3+py3*py3));
		h_eta3->Fill(eta3);
		
		tree4->GetEntry(ev);
                h_energy4->Fill(sqrt(px4*px4+py4*py4));
		h_eta4->Fill(eta4);
	}
        gStyle->SetOptStat(0);

	h_energy->Draw();
	h_energy->SetXTitle("Transverse momentum(pt)");
	h_energy->SetYTitle("Entries");
	h_energy->SetLineColorAlpha(46,1.5);
	h_energy2->SetLineColorAlpha(9,1.5);
	h_energy3->SetLineColorAlpha(8,1.5);
	h_energy4->SetLineColorAlpha(41,1.5);
	h_energy2->Draw("same");
	h_energy3->Draw("same");
	h_energy4->Draw("same");
	mgr::SetSinglePad(canvas);
        mgr::SetAxis(h_energy);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("10000 Events 2018 13TeV");
        TLegend *leg2 = new TLegend(0.55,0.65,0.8,0.85);
        leg2->SetBorderSize( 0 );
        leg2->SetFillColorAlpha( 0, 0 );
        leg2->AddEntry(h_energy, "1000 GeV ","f");
        leg2->AddEntry(h_energy2,"2000 GeV ","f");
	leg2->AddEntry(h_energy3,"3000 GeV ","f");
	leg2->AddEntry(h_energy4,"4000 GeV ","f");
        leg2->Draw();
        canvas->SaveAs("MonoPt.pdf");


	h_eta4->Draw();
        h_eta4->SetXTitle("eta");
        h_eta4->SetYTitle("Entries");
        h_eta1->SetLineColorAlpha(46,1.5);
        h_eta2->SetLineColorAlpha(9,1.5);
        h_eta3->SetLineColorAlpha(8,1.5);
        h_eta4->SetLineColorAlpha(41,1.5);
        h_eta1->Draw("same");
        h_eta2->Draw("same");
        h_eta3->Draw("same");
        mgr::SetSinglePad(canvas);
	mgr::SetAxis(h_eta4);
        mgr::DrawCMSLabelOuter(SIMULATION);
        mgr::DrawEntryRight("10000 Events 2018 13TeV");
	TLegend *leg3 = new TLegend(0.70,0.60,0.90,0.80);
        leg3->SetBorderSize( 0 );
        leg3->SetFillColorAlpha( 0, 0 );
        leg3->AddEntry(h_eta1,"1000 GeV ","f");
        leg3->AddEntry(h_eta2,"2000 GeV ","f");
        leg3->AddEntry(h_eta3,"3000 GeV ","f");
        leg3->AddEntry(h_eta4,"4000 GeV ","f");
        leg3->Draw();
        canvas->SaveAs("MonoEta.pdf");

}	
	


