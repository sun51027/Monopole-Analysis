////////////////////////////////////
//
//
//      File name: MonoDrawEtaPt.cc
//      author: Shih Lin
//      Content: Compare the eta between different mass
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
void MonoDrawEgPt()
{
        TFile *f1 = new TFile("MonoNtuple2018_MC_1000.root");
        TTree *tree1 = (TTree*)f1->Get("monopoles");
        TFile *f2 = new TFile("MonoNtuple2018_MC_2000.root");
        TTree *tree2 = (TTree*)f2->Get("monopoles");
	TFile *f3 = new TFile("MonoNtuple2018_MC_3000.root");
        TTree *tree3 = (TTree*)f3->Get("monopoles");
	TFile *f4 = new TFile("MonoNtuple2018_MC_4000.root");
        TTree *tree4 = (TTree*)f4->Get("monopoles");

	

	vector<double> photon_pt1;
	vector<double> electron_pt1;
	vector<double> photon_pt2;
        vector<double> electron_pt2;
	vector<double> photon_pt3;
	vector<double> electron_pt3;
	vector<double> photon_pt4;
	vector<double> electron_pt4;


	tree1->SetBranchAddress("pho_pt",&photon_pt1);
        tree1->SetBranchAddress("ele_pt",&electron_pt1);

        tree2->SetBranchAddress("pho_pt",&photon_pt2);
        tree2->SetBranchAddress("ele_pt",&electron_pt2);

	tree3->SetBranchAddress("pho_pt",&photon_pt3);
        tree3->SetBranchAddress("ele_pt",&electron_pt3);

	tree4->SetBranchAddress("pho_pt",&photon_pt4);
        tree4->SetBranchAddress("ele_pt",&electron_pt4);
	
	TCanvas *canvas = new TCanvas("canvas","canvas",800,600);

	TH1D *Photon = new TH1D("Photon",  "",40,0,800);
        TH1D *Photon2 = new TH1D("Photon2","",40,0,800);
	TH1D *Photon3 = new TH1D("Photon3","",40,0,800);
	TH1D *Photon4 = new TH1D("Photon4","",40,0,800);
	TH1D *Electron1 = new TH1D("Electron1","",40,0,800);
        TH1D *Electron2 = new TH1D("Electron2","",40,0,800);
        TH1D *Electron3 = new TH1D("Electron3","",40,0,800);
        TH1D *Electron4 = new TH1D("Electron4","",40,0,800);

	cout<<"total event "<<tree1->GetEntries()<<endl;
	for(unsigned ev=0;ev<tree1->GetEntries();ev++){
		cout<<"ev "<<ev<<endl;
		tree1->GetEntry(ev);
  		Photon->Fill(photon_pt1);
		Electron1->Fill(electron_pt1);

                tree2->GetEntry(ev);
  		Photon2->Fill(photon_pt2);
		Electron2->Fill(electron_pt2);
		
		tree3->GetEntry(ev);
  		Photon3->Fill(photon_pt3);
		Electron3->Fill(electron_pt3);
		
		tree4->GetEntry(ev);
  		Photon4->Fill(photon_pt4);
		Electron4->Fill(electron_pt3);
	}
        gStyle->SetOptStat(0);

	Photon->Draw();
	Photon->SetXTitle("Transverse momentum(pt)");
	Photon->SetYTitle("Entries");
	Photon->SetLineColorAlpha(46,1.5);
	Photon2->SetLineColorAlpha(9,1.5);
	Photon3->SetLineColorAlpha(8,1.5);
	Photon4->SetLineColorAlpha(41,1.5);
	Photon2->Draw("same");
	Photon3->Draw("same");
	Photon4->Draw("same");
	mgr::SetSinglePad(canvas);
        mgr::SetAxis(Photon);
        mgr::DrawCMSLabelOuter(SIMULATION);
        TLegend *leg2 = new TLegend(0.55,0.65,0.8,0.85);
        leg2->SetBorderSize( 0 );
        leg2->SetFillColorAlpha( 0, 0 );
        leg2->AddEntry(Photon, "1000 GeV ","f");
        leg2->AddEntry(Photon2,"2000 GeV ","f");
	leg2->AddEntry(Photon3,"3000 GeV ","f");
	leg2->AddEntry(Photon4,"4000 GeV ","f");
        leg2->Draw();
        canvas->SaveAs("MonoPt.pdf");


	Electron4->Draw();
        Electron4->SetXTitle("pT(GeV)");
        Electron4->SetYTitle("Entries");
        Electron1->SetLineColorAlpha(46,1.5);
        Electron2->SetLineColorAlpha(9,1.5);
        Electron3->SetLineColorAlpha(8,1.5);
        Electron4->SetLineColorAlpha(41,1.5);
        Electron1->Draw("same");
        Electron2->Draw("same");
        Electron3->Draw("same");
        mgr::SetSinglePad(canvas);
	mgr::SetAxis(Electron4);
        mgr::DrawCMSLabelOuter(SIMULATION);
        TLegend *leg3 = new TLegend(0.70,0.60,0.90,0.80);
        leg3->SetBorderSize( 0 );
        leg3->SetFillColorAlpha( 0, 0 );
        leg3->AddEntry(Electron1,"1000 GeV ","f");
        leg3->AddEntry(Electron2,"2000 GeV ","f");
        leg3->AddEntry(Electron3,"3000 GeV ","f");
        leg3->AddEntry(Electron4,"4000 GeV ","f");
        leg3->Draw();
        canvas->SaveAs("Monoelectron_pt.pdf");

}	
	


