#!/bin/bash

stm=$PWD
add="/12-Bin"
dir="$stm$add"

if [[ ! -d $dir ]]
then
  echo "$dir does exists on your filesystem."
  echo "Making directory..."
  mkdir $dir
fi

g++ 01-Code/RunControl.cpp 02-Tests/RunControlTst.cpp -I01-Code `root-config --cflags --libs` -o 12-Bin/RunControlTst.exe
