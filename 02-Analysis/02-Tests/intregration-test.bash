#!/bin/bash
#
#..  Integration test script.
#
#    Runs each test and creates integrated log file.  Check of log
#    file constitutes the test.
#
echo " CppAnal: integration tests starts"
echo " ========"
echo
#--> RunControl:
#.. Basic:
12-Bin/RunControlTst.exe
#.. No file and Bad file:
12-Bin/RunControlTst.exe -f || echo " Successfully caught no file error"
12-Bin/RunControlTst.exe -f crp || echo " Successfully caught bad file error"
#.. No directory and Bad Directory:
12-Bin/RunControlTst.exe -c
12-Bin/RunControlTst.exe -c crp || echo " Successfully caught bad directory error"
#.. No o/p directory and Bad o/p Directory:
12-Bin/RunControlTst.exe -o
12-Bin/RunControlTst.exe -o crp || echo " Successfully caught bad o/p directory error"
#.. Good commands:
touch TestFile.tmp
12-Bin/RunControlTst.exe -f TestFile.tmp && echo " Successfully loaded file"
rm TestFile.tmp
mkdir TestDir
12-Bin/RunControlTst.exe -c TestDir && echo " Successfully loaded directory"
rmdir TestDir
#.. Everything!
touch TestFile.tmp
mkdir TestDir
12-Bin/RunControlTst.exe -d -h -f TestFile.tmp -c TestDir -o 99-Scratch
rm TestFile.tmp
rmdir TestDir
#
#--> Done:
echo
echo " CppAnal: integration tests done"
echo " ========"
