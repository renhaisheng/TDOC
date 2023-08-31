"""
#######################################  Basic parameters  ######################################
#                                                                                               #
#  The input_parameters should be modified according to computer configuration or user's demand.#
#  The default_parameters should not be changed except for H_atoms of 'M062X' which should keep #
#  pace with the used method of key_words in input_parameters.                                  #
#                                                                                               #
#################################################################################################

***************************************  input_parameters  **************************************
*                                                                                               *
*   number_task: The number of tasks to run in parallel on a singel computer.                   *
*   number_process: The number of CPU cores to run in parallel for each job.                    *
*   submitted_type: The submitted type of job scripts in submitted_scripts.                     *
*   conformer_method: The constructed method of structure, 1 by openbabel, 2 by rdkit.          *
*   gaussian_version: The Gaussian version on the Linux system for M062X calculations.          *
*   key_words: The key words to control low-level calculations of optimization and frequency.   *
*                                                                                               *
*************************************************************************************************

***************************************  input_parameters  **************************************
*                                                                                               *
*   multiplicities: The spin multiplicities in smiles to classify different species.            *
*   submitted_scripts: The submitted scripts to run jobs on the Linux system.                   *
*   electrons_atoms: The electrons of atoms. This parameter should not be changed.              *
*   bond_marks: The specific marks of bonds for linear bond compensation method.                *
*   Hf_atoms: The experimental enthalpies of formation for atoms in atomization method.         *
*   H_atoms: The enthalpies of atoms in atomization method. H_atoms = SPE + H_corr(0-->298.15K).*
*                                                                                               *
*************************************************************************************************
"""

input_parameters(
    number_task = 4,
    number_process = 8,
    submitted_type = 1,
    conformer_method = 1,
    gaussian_version= 'g16',
    key_words = '# opt freq=hinderedrotor B3LYP/6-31g(d,p) symm=veryloose em=gd3bj',
    zpe_scale = 0.9838
)

default_parameters(
    multiplicities = {'S': 1, 'D': 2, 'T': 3, 'Q': 4, 'P': 5},
    submitted_scripts = {1: 'job-nohup', 2: 'job-qsub'},
    electrons_atoms = {'C': 6, 'H': 1, 'O': 8},
    bond_marks = {1.0: '-', 2.0: '=', 3.0: '#', 1.5: '$'},
    Hf_atoms = {'C': 0.2729652, 'H': 0.0830318, 'O': 0.0949038},
    H_atoms = {'MP2' : {'C': -37.7288488, 'H': -0.4958729, 'O': -74.8756186},
                'CCSDT' : {'C': -37.7898614, 'H': -0.4976322, 'O': -75.004833},
                },
)
