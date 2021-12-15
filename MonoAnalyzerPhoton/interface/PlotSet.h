#ifndef _PLOTSET_H_
#define _PLOTSET_H_
static const unsigned nPlot = 100U;
enum PlotName{
  FracSatVNstrips = 0, // fraction of saturated strips vs. number of strips
  DedXSig,             // dE/dX significance
  XYPar0,//d0
  XYPar1,//phi0
  XYPar2,//radius
  RZPar0,
  RZPar1,
  RZcurv,              // RZ curvature/uncertainty
  E55,
  F51,                 // frac 51
  HcalIso,             // Hcal Iso
  Dist,
  ABCD,
  Spike,

// for data
  F51_slice1,
  F51_slice2,
  F51_slice3,
  F51_slice4,
  F51_slice5,
  F51_slice6,
  dEdXSig_slice1,
  dEdXSig_slice2,
  dEdXSig_slice3,
  dEdXSig_slice4,
  dEdXSig_slice5,
  dEdXSig_slice6,
  dEdXSig_slice7,

  EcalBarrel,
  EcalEndCup,
  EcalAll
};
class PlotSet
{
public:
  PlotSet(){
	plots_.resize(nPlot);
	profile_.resize(nPlot);}
  ~PlotSet(){}
  void CreatPlot(const PlotName pn, TH1* h){ 
	plots_[pn] = h;
  }
  void CreatProfile(const PlotName pn, TProfile *f){ 
	profile_[pn] = f;
  }
  TH1 * GetPlot(const PlotName pn){ return plots_[pn]; }
  TH1 * GetProfile(const PlotName pn){ return profile_[pn]; }
  void WritePlot(){
	for(int pn=0;pn<nPlot;pn++){
	TH1 *h = plots_[pn];
	if(h){ 
	  h->Write();
	  }
	}
  }
  void WriteProfile(){
	for(int pn=0;pn<nPlot;pn++){
 	   TProfile *f = profile_[pn];
	   if(f)  f->Write();
	}
  }
private:
  vector<TH1*> plots_;
  vector<TProfile*> profile_;
};
/*class DrawPlot
{
public:
  DrawPlot()
  ~DrawPlot()
  void DrawPlot(const PlotName pn, TH1 *h){
	plots_[pn] = h;
	Scale_(plots_[pn];	
	plots_[pn]->Draw();
  }
  void Scale_(TH1 *h){
	Double_t factor = 1.;
	h->Scale(factor/h->Integral());
  }
}*/
#endif

