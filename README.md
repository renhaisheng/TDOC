## TDOC:
TDOC: Thermodynamic Data Off-line Calculator 

Date: July 25, 2022.

Author: Huajie Xu, Zihan Xu, Lu Liu, Zerong Li, Quan Zhu, Haisheng Ren.

Email: renhs@scu.edu.cn.

Copyright: Center for Combustion Dynamics, Sichuan University.


## Characteristic:
TDOC is developed to obtain accurate thermodynamic parameters of mechanisms within chemical accuracy for CHO carbon hydrogen.
It can automatically build input files and deal with occurred errors, such as spin contaminations, imaginary frequencies, 
multireference effects, and so on. After all errors are eliminated, the accurate thermodynamic data will be generated in the
format of 14 parameters for Chemkin use.

For flexible larger molecules, the enthalpies of formation of larger molecules are derived by CBH-3 extrapolation with bond 
additivity corrections, dramatically decreasing computational costs. For aromatic molecules, the enthalpies of formation 
are processed by the reservation of aromatic rings in CBH-3 rung, while the smaller aromatic molecules are directly calculated
by CCSD(T)/CBS method by symmetry acceleration. The extremely large polycyclic aromatic molecules are limited for its too 
expensive calculations of aromatic rings. Combined the corrections of hindered rotor and conformational sampling, the calculated
accuracy generally meets the requirements of chemical accuracy.


## Preparation:
1. Platform: Windows or Linux. Note that the Linux system is required to get QM calculations of Gaussian and Molpro.
2. Environment: Conda with Python (>=3.7). Conda can be obtained by Anaconda or Miniconda (https://docs.conda.io/en/latest/).
3. Modules: requirements.txt. Enter "conda install --yes -c conda-forge --file requirements.txt" or "python setup.py install".
4. Extensions: Gaussian, Mopro, xtb, CREST packages. Gaussian and Mopro packages are indispensable to get optimized structure 
   and single point energies. For conformational sampling, open-source CREST with xtb are taken from https://github.com/grimme-lab/.


## setup
Enter "python setup.py install" in conda shell for windows or bash shell for Linux.


## Usage:
1. Enter "python tdoc.py smiles.txt -p parameters.py" in conda shell for windows or bash shell for Linux to automatically generate
   input scripts and files. 
2. Copy these files to Linux platform for Gaussian, Molpro, and CREST calculations.
3. Input "chmod +x job-nohup && nohup ./job-nohup &" to run jobs. In particular, enter "chmod +x job-qsub && nohup ./job-qsub &" 
   to start jobs for a cluster of PBS.
4. Copy all output files to the directory of gaussian_out.
5. Repeat previous steps until the all tasks are completed.


## Directories:
1. preexisted-dat: Thermodynamic data files for bath gas, such as Ar, N2 and so on.
2. chemkin_dat: Generated thermodynamic data files with CHEMKIN format by TDOC.
3. gaussian_out: Input and output files of QM calculations for Linux system.
4. program_out: Some viewable and summary analysis for all output results of TDOC.
5. submitted_inp: Submitted scripts and files for QM calculations on Linux system.

