# rmc-irrep-analysis
Representational Analysis of Extended Disorder in Atomistic Ensembles

Journal Reference: J. Appl. Cryst. 48, 1560-72 (2015)

This page serves as a repository for all of the necessary computer code (written in C) to extract quantitative information regarding structural disorder in big boxes (atomistic ensembles) through the use of a locally symmetry adapted tight binding basis.


Usage

To use the code, download the source. Modify Makefile to match the local compiler settings.

Build the executables using "make", and then run rmc-symmetry on each of the datasets. More complex problems may require you to use the individual transform/rebuild routines.


Important Notes

The included example demonstrates the format of the necessary input files. 
