#!/bin/bash

g++ 01-Code/RunControl.cpp 01-Code/Analysis.cpp 02-Tests/AnalysisTst.cpp -I01-Code `root-config --cflags --libs` -o 12-Bin/AnalysisTst.exe
