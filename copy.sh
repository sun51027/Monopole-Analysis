#!/bin/bash
set -e  #if there is error , stop running
#scp -rp lshih@lxplus.cern.ch:~/CMSSW_10_6_8_patch1/src/Monopoles . 
#scp -rp lshih@lxplus.cern.ch:~/test1/CMSSW_8_0_33_UL/src/linshih_crab/crab_Monopole_SpinHalf_M-1000_DY_13TeV_v1/ ./MonoNtuple2016/
#scp -rp lshih@lxplus.cern.ch:~/OldestTest/CMSSW_10_6_8_patch1/src/MonopoleNtuple1000_2017_1000evt.root .
#scp -rp lshih@lxplus.cern.ch:/eos/user/l/lshih/MonoNtuple2018_MC_1000.root ./SmallerData/
#scp -rp lshih@lxplus.cern.ch:/eos/user/l/lshih/MonoNtuple2018_MC_3000.root .
#scp -rp lshih@lxplus.cern.ch:/eos/user/l/lshih/MonoNtuple2018_MC_2000.root .
scp -rp lshih@lxplus.cern.ch:/eos/user/l/lshih/MonoNtuple2018_MC_1000_TestSwisscross2.root .
