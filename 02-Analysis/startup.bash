#!/bin/bash

stm=$PWD
#echo $stm

dir="$stm"
echo "Set Analysis path:"
AnalysisPATH="$dir"
echo "    " $AnalysisPATH
export AnalysisPATH

add="/01-Code"
#echo $add
dir="$stm$add"
#echo $dir
echo "Set AnalysisCPATH:"
AnalysisCPATH="$dir"
echo "    " $AnalysisCPATH
export AnalysisCPATH
