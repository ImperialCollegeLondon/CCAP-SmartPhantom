#!/bin/bash

g++ $AnalysisCPATH/RunControl.cpp $AnalysisCPATH/Analysis.cpp  $AnalysisPATH/03-Skeleton/Analysis-skeleton.cpp \
    -I$AnalysisCPATH `root-config --cflags --libs` -o $AnalysisPATH/12-Bin/Analysis-skeleton.exe
