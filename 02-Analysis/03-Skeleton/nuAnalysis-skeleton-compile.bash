#!/bin/bash

g++ $nuAnalysisCPATH/RunControl.cpp $nuAnalysisCPATH/nuAnalysis.cpp ./nuAnalysis-skeleton.cpp -I$nuAnalysisCPATH `root-config --cflags --libs` -o nuAnalysis-skeleton.exe
