# -*- coding: utf-8 -*-
"""
TDOC: Thermodynamic Data Offline Calculator. 
Date: March 3, 2022.
Author: Huajie Xu, Zihan Xu, Lu Liu, Quan Zhu, Haisheng Ren.
Email: renhs@scu.edu.cn.
Copyright: Center for Combustion Dynamics, Sichuan University.
"""

from src._thermo import THERMO
from shutil import which
import logging as log
log.getLogger().setLevel(log.INFO)


def main():
    log.info('\n\nTDOC: Thermodynamic Data offline Calculator.\nDate: March 24, 2022.\nEmail: renhs@scu.edu.cn.\n'
            'Author: Huajie Xu, Zihan Xu, Lu Liu, Zerong Li, Quan Zhu, Haisheng Ren.\nCopyright: Center for Combustion '
            'Dynamics, Sichuan University.\n')
    
    log.info('\n\nThis program of TDOC is developed to obtain accurate thermodynamic parameters of mechanisms within '
            'chemical accuracy at CCSD(T)/CBS level with corrections of BAC. TDOC can automatically build input files'
            ' and deal with occurred errors, such as spin contaminations, imaginary frequencies, multireference effects,'
            'and so on. After all errors are eliminated, the accurate thermodynamic data will be generated '
            'in the format of 14 parameters for Chemkin use.\n')

    log.info('\n============ START ============\n')
    data = THERMO()
    data.Get_submit_out()
    data.Get_thermo_out()
    data.Write_thermo_dat()
    log.info('\n============  END  ============\n')

if __name__ == '__main__':
    main()




    