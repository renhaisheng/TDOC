# TDOC
Automated thermdyamic data generation

TDOC: Thermodynamic Data Offline Calculator. 
Date: July 18, 2022.
Author: Huajie Xu, Zihan Xu, Lu Liu, Zerong Li, Quan Zhu, Haisheng Ren.
Email: renhs@scu.edu.cn.
Copyright: Center for Combustion Dynamics, Sichuan University.

## Characteristic:
This program of TDOC is developed to obtain accurate thermodynamic parameters of mechanisms within chemical accuracy at CCSD(T)/CBS level with corrections of BAC. 
TDOC can automatically build input files and deal with occurred errors, such as spin contaminations, imaginary frequencies, multireference effects, and so on. After all errors are eliminated, the accurate thermodynamic data will be generated in the format of 14 parameters for Chemkin use.


## Directories:

1. chemkin_dat: Final thermodynamic paramters file for CHEMKIN format.
2. gaussian_out: Input and output files of calculations for Linux system.
3. preexisted-dat: Already existed thermodyamic paramters file which will not be calculated any more.
4. program_out: Some viewable results of program output.
5. submitted_inp: Submitted files to be calculated in Linux system.
6. SRC: Source code files.


## Preparation

Install software and add channel to environment variables on Windows platforms.
1. Obtain conda by Anaconda or Miniconda (https://docs.conda.io/en/latest/) where the Python version is greater than 3.7.

The required modules are rdkit, utilspie, ase.
1. Use "conda" or "pip" command to install any vacant module.

Install essential packages in Linux platforms for high performance computing.
1. For M062X calculations, the Gaussian of g09 or g16 is required. You can change the  default word of "g09" to "g16" in submit.py file if you have Gaussian copyright of g16 version.  
2. For CCSD(T) computation, the Molpro is indispensable.
3. For GFNFF conformers sampling, the CREST is available from GitHub (https://github.com/grimme-lab/crest). The required xtb binary may also be taken from GitHub (https://github.com/grimme-lab/xtb).

Note: The operation on Windows system may also be replaced in Linux system by making minor modifications to source code. This program is open source and can be modified as required.

## Usage

Before running  program, the smiles.txt file which contains SMILES of all species in the mechanism is required in main directory. Moreover, the input parameters can be modified in the input.ini file as needed. 

To run TDOC.py in main directory, you can type command of 'python  TDOC.py' in shell or operate directly in python interpreter. Or you can double click the TDOC.exe to complete the above operation.

