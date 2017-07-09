#!/bin/bash

g++ -O3 -Wall -std=c++11 ../ccfiles/langevin.cc -o langevin -lm -lvoro++
