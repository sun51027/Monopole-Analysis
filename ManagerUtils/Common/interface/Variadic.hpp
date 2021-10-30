/*******************************************************************************
*
*  Filename    : Variadic.hpp
*  Description : Template functions for aiding with Variadic functions of same type
*                by contracting all arguments into a std::vector<T>
*  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
*
*******************************************************************************/
#ifndef MANAGERUTILS_COMMON_VARIADIC_HPP
#define MANAGERUTILS_COMMON_VARIADIC_HPP

#include <vector>

namespace mgr{
    /*******************************************************************************
    *   Single Entry function
    *******************************************************************************/
    template<typename OBJTYPE, typename ARGTYPE>
    std::vector<OBJTYPE> MakeVector( ARGTYPE first );

    /*******************************************************************************
    *   Multi entry argument input
    *******************************************************************************/
    template<typename OBJTYPE, typename ARGTYPE, typename ... VARTYPE>
    std::vector<OBJTYPE> MakeVector( ARGTYPE first, VARTYPE ... arg );
};/* mgr */

/*******************************************************************************
*   Including implementation src code
*******************************************************************************/
#include "/wk_cms2/shihlin0314/CMSSW_8_0_29/src/ManagerUtils/Common/src/Variadic.ipp"


#endif/* end of include guard: MANAGERUTILS_COMMON_VARIADIC_HPP */
