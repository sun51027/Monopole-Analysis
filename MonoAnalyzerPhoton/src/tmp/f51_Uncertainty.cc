
using namespace std;
void f51_Uncertainty(string year, string mass)
{

        TFile *fin1 = new TFile (("./output/MonoPhotonAnalysis_"+year+"_"+mass+".root").c_str());
        TDirectory *NOTRG1 = (TDirectory*)fin1->Get("HLT_Photon200");
        NOTRG1->cd();
        TFile *fin2 = new TFile (("./output/DeltaRayOff_Analysis_"+year+"_"+mass+".root").c_str());
        TDirectory *NOTRG2 = (TDirectory*)fin2->Get("HLT_Photon200");
        NOTRG2->cd();
        TFile *fin3 = new TFile (("./output/EcalUncertainty_Analysis_"+year+"_"+mass+".root").c_str());
        TDirectory *NOTRG3 = (TDirectory*)fin3->Get("HLT_Photon200");
        NOTRG3->cd();

        TH1D *f51; 
        TH1D *f51_deltaRayOff; 
        TH1D *f51_spikeOff; 
        f51  = (TH1D*)NOTRG1->Get("Flow_dEdxSig_F51");	 
        f51_deltaRayOff  = (TH1D*)NOTRG2->Get("Flow_dEdxSig_F51");	 
        f51_spikeOff  = (TH1D*)NOTRG3->Get("Flow_dEdxSig_F51");	 
	f51->Scale(1/f51->Integral());
	f51_deltaRayOff->Scale(1/f51_deltaRayOff->Integral());
	f51_spikeOff->Scale(1/f51_spikeOff->Integral());
	
	ofstream fout("f51_uncert"+mass+".csv");

	fout<<year<<""<<mass+",f51 deltaRay Off,f51 spike Off,f51 default"<<endl;
	fout<<"Mean (from f51 distribution),"<<f51_deltaRayOff->GetMean()<<","<<f51_spikeOff->GetMean()<<","<<f51->GetMean()<<endl;
	fout<<"Uncertainty,"<<(f51_deltaRayOff->GetMean()-f51->GetMean())/f51->GetMean()<<","<<(f51_spikeOff->GetMean()-f51->GetMean())/f51->GetMean()<<endl;

	cout<<"f51 default "<<f51->GetMean()<<endl;
	cout<<"f51 deltaRay Off "<<f51_deltaRayOff->GetMean()<<endl;
	cout<<"f51 spike Off "<<f51_spikeOff->GetMean()<<endl;
}
