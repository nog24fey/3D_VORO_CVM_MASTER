#!/bin/bash
icpc -O3 -xHOST -ipo -no-prec-div -Wall -std=c++11 -I../../VOROPPINSTALL40/voro++-0.4.6/src/ -L../../VOROPPINSTALL40/voro++-0.4.6/src/ -L/usr/local/lib64/ langevin1.cc -o langevin1 -lm -lvoro++
