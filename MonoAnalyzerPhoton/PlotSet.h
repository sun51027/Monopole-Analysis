#ifndef _PLOTSET_H_
#define _PLOTSET_H_
static const unsigned nPlot = 20U;
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
  Spike
};
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
#endif

