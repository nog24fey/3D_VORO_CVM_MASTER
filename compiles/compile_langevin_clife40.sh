#!/bin/bash
icpc -O3 -xHOST -ipo -no-prec-div -Wall -std=c++11 -I../../../VOROPPINSTALL40/voro++-0.4.6icc/src/ -L../../../VOROPPINSTALL40/voro++-0.4.6icc/src/ ../ccfiles/langevin.cc -o langevin -lm -lvoro++
