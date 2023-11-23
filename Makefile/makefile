#!/bin/bash
all:
	g++ ansatz.cpp main.cpp per_analyzer.cpp PIC_general.cpp PIC_LVP.cpp PIC_VP.cpp rad_fct.cpp steadystate.cpp XV.cpp -O3 -pthread -static-libstdc++
run: 
	./a.out

.PHONY: all run
