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
*   calculated_method: The calculated method of high level in calcualted_levels.                *
*   conformer_method: The constructed method of structure, 1 by openbabel, 2 by rdkit.          *
*   gaussian_version: The Gaussian version on the Linux system for M062X calculations.          *
*   aro_rings: The reserved aromatic rings of larger aromatic molecules in CBH extrapolation.   *
*   key_words: The key words to control low-level calculations of optimization and frequency.   *
*                                                                                               *
*************************************************************************************************

***************************************  input_parameters  **************************************
*                                                                                               *
*   multiplicities: The spin multiplicities in smiles to classify different species.            *
*   submitted_scripts: The submitted scripts to run jobs on the Linux system.                   *
*   calcualted_levels: The calculated levels to obtain accurate enthalpies of formation.        *
*   electrons_atoms: The electrons of atoms. This parameter should not be changed.              *
*   Hf_atoms: The experimental enthalpies of formation for atoms in atomization method.         *
*   H_atoms: The enthalpies of atoms in atomization method. H_atoms = SPE + H_corr(0-->298.15K).*
*                                                                                               *
*************************************************************************************************
"""

input_parameters(
    number_task = 4,
    number_process = 8,
    submitted_type = 1,
    calculated_method = 1,
    conformer_method = 1,
    gaussian_version= 'g09',
    aro_rings = ['C12=CC=CC=C1C=CC=C2','C1=CC=CC=C1','C1=COC=C1'],
    key_words = '# opt freq=hinderedrotor M062X/6-311++g(d,p) symm=veryloose',
)

default_parameters(
    multiplicities = {'S': 1, 'D': 2, 'T': 3, 'Q': 4, 'P': 5},
    submitted_scripts = {1: 'job-nohup', 2: 'job-qsub'},
    calcualted_levels = {1: 'CCSDT', 2: 'G4', 3: 'CBSQB3'},
    electrons_atoms = {'C': 6, 'H': 1, 'O': 8},
    Hf_atoms = {'C': 0.2729652, 'H': 0.0830318, 'O': 0.0949038},
    H_atoms = {'M062X' : {'C': -37.8380477, 'H': -0.4958344, 'O': -75.0586831},
                'CCSDT' : {'C': -37.7900384, 'H': -0.4976322, 'O': -75.0050381},
                'G4' : {'C': -37.8318079, 'H': -0.4990600, 'O': -75.0431406},
                'CBSQB3' : {'C': -37.7830150, 'H': -0.4974574, 'O': -74.9852586},
                },
)
