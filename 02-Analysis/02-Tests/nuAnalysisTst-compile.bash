#!/bin/bash

g++ 01-Code/RunControl.cpp 01-Code/nuAnalysis.cpp 02-Tests/nuAnalysisTst.cpp -I01-Code `root-config --cflags --libs` -o 12-Bin/nuAnalysisTst.exe
