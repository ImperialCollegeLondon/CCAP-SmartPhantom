#!/bin/bash

g++ $AnalysisCPATH/RunControl.cpp $AnalysisCPATH/Analysis.cpp ./Analysis-skeleton.cpp -I$AnalysisCPATH `root-config --cflags --libs` -o Analysis-skeleton.exe
