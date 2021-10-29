////////////////////////////////////
//
//
//      File name: MonoDrawE55.cc
//      author: Shih Lin
//      Content: Compare the rec mono e55
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
#include "TH1D.h"
#include "math.h"
#include <algorithm>

using namespace std;
void MonoDrawE55()
{
        TFile *f1 = new TFile("MonoNtuple2018_MC_1000_0722.root");
        TTree *tree1 = (TTree*)f1->Get("monopoles");
        TFile *f2 = new TFile("MonoNtuple2018_MC_2000_0722.root");
        TTree *tree2 = (TTree*)f2->Get("monopoles");
	TFile *f3 = new TFile("MonoNtuple2018_MC_3000_0722.root");
        TTree *tree3 = (TTree*)f3->Get("monopoles");
	TFile *f4 = new TFile("MonoNtuple2018_MC_4000_0722.root");
        TTree *tree4 = (TTree*)f4->Get("monopoles");


	vector<double> *e55_1 = 0;	
	vector<double> *e55_2 = 0;	
	vector<double> *e55_3 = 0;	
	vector<double> *e55_4 = 0;	
        unsigned nCandidates1;
        unsigned nCandidates2;
        unsigned nCandidates3;
        unsigned nCandidates4;

        tree1->SetBranchAddress("cand_N",&nCandidates1);
        tree2->SetBranchAddress("cand_N",&nCandidates2);
        tree3->SetBranchAddress("cand_N",&nCandidates3);
        tree4->SetBranchAddress("cand_N",&nCandidates4);
	tree1->SetBranchAddress("cand_e55",&e55_1);
        tree2->SetBranchAddress("cand_e55",&e55_2);
	tree3->SetBranchAddress("cand_e55",&e55_3);
	tree4->SetBranchAddress("cand_e55",&e55_4);
	
	TCanvas *canvas = new TCanvas("canvas","canvas",800,600);

	TH1D *h_energy = new TH1D("h_energy","",30,-1,1000);
        TH1D *h_energy2 = new TH1D("h_energy2","",30,-1,1000);
	TH1D *h_energy3 = new TH1D("h_energy3","",30,-1,1000);
	TH1D *h_energy4 = new TH1D("h_energy4","",30,-1,1000);
	
	cout<<"total event "<<tree1->GetEntries()<<endl;
	for(unsigned ev=0;tree1->GetEntries();ev++){
		cout<<"ev "<<ev<<endl;

		tree1->GetEntry(ev);
                for(unsigned i=0;i<nCandidates1;i++){
			h_energy->Fill((*e55_1)[i]);
		}
	
                tree2->GetEntry(ev);
                for(unsigned i=0;i<nCandidates2;i++){
                	h_energy2->Fill((*e55_2)[i]);
		}	
	
		tree3->GetEntry(ev);
                for(unsigned i=0;i<nCandidates3;i++){
                	h_energy3->Fill((*e55_3)[i]);
		}		

		tree4->GetEntry(ev);
                for(unsigned i=0;i<nCandidates4;i++){
                	h_energy4->Fill((*e55_4)[i]);
		}
	}
        gStyle->SetOptStat(0);
	gPad->SetLogy();

	h_energy->Draw();
	h_energy->SetXTitle("e55(GeV)");
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
        mgr::DrawEntryRight("15000 Events 2018 13TeV");
        TLegend *leg2 = new TLegend(0.55,0.65,0.8,0.85);
        leg2->SetBorderSize( 0 );
        leg2->SetFillColorAlpha( 0, 0 );
        leg2->AddEntry(h_energy, "1000 GeV ","f");
        leg2->AddEntry(h_energy2,"2000 GeV ","f");
	leg2->AddEntry(h_energy3,"3000 GeV ","f");
	leg2->AddEntry(h_energy4,"4000 GeV ","f");
        leg2->Draw();
        canvas->SaveAs("MonoE55.pdf");



}	
	


