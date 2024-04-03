====== Description ======
The C-code DFMSPH22 is designed to obtain the nucleus-nucleus potential using the double folding model (DFM) 
and in particular to find the Coulomb barrier.

The code consists of the following files (subfolder \DFMSPH22\0_DFMSPH22_code):
1) Header file <dfmsph_def.h> contains the instructions for the preprocessor, the definitions of some constants, 
	the descriptions of the global variables, and declarations of all functions.
2) File <dfmsph_mai.c> is the main file containing several functions including the function main(). 
	These functions calculate the potential and find the Coulomb barrier.
3) File <dfmsph_fun.c> contains a set of auxiliary functions.
4) File <dfmsph_inp.c> is responsible for reading the input information from all input files and 
	checking the values of the input parameters.
5) File <dfmsph_pot.c> is responsible for calculation of the DFM potential.
6) File <dfmsph_out.c> is responsible for printing the output information into the output files.
7) File <dfmsph_exc.c> contains a number of functions which are designed to calculate the finite-range exchange part 
	of the nuclear term and the exchange part of the Coulomb term of the DFM potential.
Input and output files (folders \DFMSPH22\Test_run_1 and \DFMSPH22\Test_run_2):
1) Input file <INP_DFMSPH22.c> is the principal input file. The parameters of the target and projectile nuclei, 
	the kind of the NN forces, etc. are defined here. There are also several key parameters. 
	An instruction for using <INP_DFMSPH22.c> is located in file <Program_changes.txt>.
2) Input files <inp_rhoP.c>, <inp_rhoT.c> contain the input densities for the projectile (<inp_rhoP.c>) and 
	target (<inp_rhoT.c>) nuclei as the functions of distance from the center of the nucleus. 
	More detailed description is located in file <Program_changes.txt>.
3) Input file <INP_NN_forces.c> contains the coefficients of the M3Y, RMF, and Migdal NN forces. 
	This file is not expected to be modified by the user although one can alternate the sets of the RMF forces 
        keeping the values of key_vNN.
4) Input file <NUMBER.c> provides the ordering number of the actual run of the program.
5) Output file <out_dfmsph22.c> is the principal (cumulative) output file. After each run of the code, 
	a new line appears at the end of the file. This line contains about 50 output quantities and some text information. 
	A more detailed description is located in file <Program_changes.txt>.
6) Output file <out_inp_dfm22.c> contains the input information used in each run. After each run of the code, 
	a new line appears at the end of this file. This line contains about 50 input quantities and some text information. 
	A more detailed description is located in file <Program_changes.txt>.
7) Output file <out_vNN_.c> contains the r-dependence of the NN-interaction used for calculating the
          nucleus-nucleus potential. Some input and output parameters are printed at the end of this file.
8) Output file <OUT_U_R.c> contains the components of the nucleus-nucleus potential as the functions 
	of the center of mass distance and the fits of the nuclear term. Some input and output parameters are printed 
	at the end of this file.

====== How to run ======
1) replace the values of the input parameters in the file <INP_DFMSPH22.c> by those corresponding to your problem;
2) prepare the density profiles for projectile and target nuclei (files <inp_rhoP.c> and <inp_rhoT.c>, respectively);
3) in the subfolder \DFMSPH22\0_DFMSPH22_code, run the file <compile.bat> under Windows or <makefile> under 
	MacOS/Linux/Unix to produce executable file;
4) move the resulting executable file <dfmsph22.exe> or <dfmsph22> into the external folder with the input files, 
	e.g. into \DFMSPH22\Test_run_1;
5) run the executable file.
