#!/bin/bash

rm -rf results

mkdir -p results/bacteria results/cytokine results/apc_n results/apc_a results/antibody results/tissue
mkdir -p results/b_memory results/b_active results/b_naive results/apc_a_linf results/antibody_linf  

gcc -o execute main.c source/*.c -Iinclude -fopenmp -lm -g
