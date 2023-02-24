## TDOC-2.0:
TDOC: Thermodynamic Data Off-line Calculator 

Date: February 23, 2023.

Author: Huajie Xu, Bo wang, Quan Zhu, Haisheng Ren.

Email: renhs@scu.edu.cn.

Copyright: Center for Combustion Dynamics, Sichuan University.

Cite: 10.1016/j.fuel.2022.125431


## Characteristic:
TDOC-2.0 is developed to obtain accurate thermodynamic parameters of mechanisms within chemical accuracy for CHO carbon hydrogen.
It can automatically build input files and deal with occurred errors, such as spin contaminations, imaginary frequencies and so on.
After all errors are eliminated, the accurate thermodynamic data will be generated in the format of 14 parameters for Chemkin use.

The enthalpies of formation of larger molecules are derived by CBH-BDC method with BAC-EP correction, while those of small molecules
are evaluated by CCSD(T)/CBS method. Combined the corrections of hindered rotor and conformational sampling, the calculated accuracy
of all thermodynamic data generally meets the requirements of chemical accuracy.


## Preparation:
1. Platform: Windows or Linux. Note that the Linux system is required to get QM calculations of Gaussian and Molpro.
2. Environment: Conda with Python (>=3.7). Conda can be obtained by Anaconda or Miniconda (https://docs.conda.io/en/latest/).
3. Modules: requirements.txt. Enter "conda install --yes -c conda-forge --file requirements.txt".
4. Extensions: Gaussian, Mopro, xtb, CREST packages. Gaussian and Mopro packages are indispensable to get optimized structure 
   and single point energies. For conformational sampling, open-source CREST with xtb are taken from https://github.com/grimme-lab/.


## Setup
1. Unpack "tdoc-1.0.tar.gz" or "tdoc-1.0.zip" file and enter the main directory of "tdoc-1.0".
2. Enter "conda install --yes -c conda-forge --file requirements.txt" to install required modules.
3. Add the main directory to system environment. For Linux system, the TDOC.py file needs key words of ':set ff=unix' in vim editor 
   to convert dos to unix. Additionally, enter "chmod +x tdoc.py" to make an excutable file in Linux shell. 


## Usage:
1. Enter "tdoc.py smiles.txt -p parameters.py" in another directory to automatically generate input scripts and files. 
2. Copy submitted files in the newly generated "submitted_inp" directory to Linux platform for Gaussian, Molpro, and CREST calculations.
3. Input "chmod +x job-nohup && nohup ./job-nohup &" to run jobs. In particular, enter "chmod +x job-qsub && nohup ./job-qsub &" 
   to start jobs for a cluster of PBS.
4. Copy all output files to the directory of "gaussian_out".
5. Repeat previous steps until the all tasks are completed.
6. If BAC or BDC parameters are missing, one can utilize the corresponding scripts to train parameters based on available data.


## Examples:
1. smiles.txt: The basic file for input species with their SMILES. It is based on manually input by users.
2. parameters.py: The control parameters should be set according to system configuration and user demand.
3. BAC_parameters.txt: The BAC paramters for bond additivity correction which are continuously updated. 
4. BDC_parameters.txt: The BDC paramters for bond difference correction which are continuously updated. 
5. chemkin_dat: The data directory to store generated thermodynamic parameters files by TDOC for CHEMKIN use.
6. gaussian_out: The output directory to deposit input and output files of QM calculations for Linux system.
7. program_out: The result directory to contain some summary analysis and overall data files for species.
8. submitted_inp: The submitted directory for QM calculations on Linux system. It will be removed if all data are completed.
