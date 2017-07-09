#!/bin/bash
icpc -O3 -xHOST -ipo -no-prec-div -Wall -std=c++11 -I../../../VOROPPINSTALL90/voro++-0.4.6icc/src/ -L../../../VOROPPINSTALL90/voro++-0.4.6icc/src/ ../forceshorten.cc -o forceshorten -lm -lvoro++
