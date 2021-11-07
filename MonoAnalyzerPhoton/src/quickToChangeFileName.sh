#!/bin/bash
sed -n '701p' MonoDedxCrossTalk.cc
read -p "Enter old X: " oldX
read -p "Enter new X: " newX
read -p "Enter old up/down : " oldUpDown
read -p "Enter old up/down again with capital: " oldCapitalUpDown
read -p "Enter new up/down : " newUpDown
read -p "Enter new up/down again with capital: " newCapitalUpDown
#read -p "please enter y or n: " yn
#[ "${yn}" == "Y" -o "${yn}" == "y" ] && echo "OK, continue" && exit 0
#[ "${yn}" == "N" -o "${yn}" == "n" ] && echo "Oh, interrupt!" && exit 0
sed -i.bak "701s|${oldX}${oldUpDown}|${newX}${newUpDown}|" MonoDedxCrossTalk.cc #change output file name
sed -i.bak "702s|${oldX}|${newX}|1" MonoDedxCrossTalk.cc #change directory
sed -i.bak "702s|${oldCapitalUpDown}|${newCapitalUpDown}|1" MonoDedxCrossTalk.cc #change directory
sed -i.bak "702s|${oldX}${oldUpDown}|${newX}${newUpDown}|" MonoDedxCrossTalk.cc #change input file name
echo "#####Finish the change!#####"
sed -n '701,702p' MonoDedxCrossTalk.cc
echo ""
echo "--Start to run root --"
root -l MonoDedxCrossTalk.cc -q

