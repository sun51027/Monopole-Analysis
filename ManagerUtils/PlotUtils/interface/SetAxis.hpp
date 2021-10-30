/*******************************************************************************
*
*  Filename    : SetAxis.hpp
*  Description : Template functions for styling Axis, using template functiosn
*                for ROOT, ROOFIT object compatibility.
*  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
*
*******************************************************************************/
#ifndef MANAGERUTILS_PLOTUTILS_SETAXIS_HPP
#define MANAGERUTILS_PLOTUTILS_SETAXIS_HPP

#include "TGaxis.h"

namespace mgr{
    template<class T>
    void SetAxis( T* plot );

    template<class T>
    void SetTopPlotAxis( T* plot );

    template<class T>
    void SetBottomPlotAxis( T* plot );

    template<class T>
    void SetNormToUnity( T* plot );
};

#include "./ManagerUtils/PlotUtils/src/SetAxis.ipp"

#endif/* end of include guard: MANAGERUTILS_PLOTUTILS_SETAXIS_HPP */
