/*******************************************************************************
*
*  Filename    : Common.hpp
*  Description : Typical helper setups for plotting objects
*  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
*
*  Making it easier to comply with standards in
*   https://ghm.web.cern.ch/ghm/plots/
*
*******************************************************************************/
#ifndef MANAGERUTILS_PLOTUTILS_COMMON_HPP
#define MANAGERUTILS_PLOTUTILS_COMMON_HPP

/*******************************************************************************
*   Including commong ROOT plotting objects
*******************************************************************************/
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TPaveText.h"

/*******************************************************************************
*   Including custom defined functions
*******************************************************************************/
#include "./ManagerUtils/PlotUtils/interface/Constants.hpp"
#include "./ManagerUtils/PlotUtils/interface/DrawCMS.hpp"
#include "./ManagerUtils/PlotUtils/interface/GetRange.hpp"
#include "./ManagerUtils/PlotUtils/interface/LatexMgr.hpp"
#include "./ManagerUtils/PlotUtils/interface/PlotObj.hpp"
#include "./ManagerUtils/PlotUtils/interface/SaveUtils.hpp"
#include "./ManagerUtils/PlotUtils/interface/SetAxis.hpp"
#include "./ManagerUtils/PlotUtils/src/SetAxis.ipp"
#include "./ManagerUtils/SysUtils/interface/PathUtils.hpp"

#endif/* end of include guard: MANAGERUTILS_PLOTUTILS_COMMON_HPP */
