# -*- coding: utf-8 -*-
from os.path import abspath,dirname

from tdoc._submit import Submit
from tdoc._thermo import Thermo
from tdoc._parameters import Parameters



class TDOC(object):
    
    """ To initiate parameters for TDOC. """
    def __init__(self, input_file, para_file):
        super(TDOC, self).__init__()
        self.input_file = input_file
        self.para_file = para_file
        self.work_path = dirname(abspath(self.input_file))


 
    """ To execute programm. """
    def execute(self):
        parameters = Parameters(self.input_file,self.para_file,self.work_path)
        self.para = parameters.get_all_parameters()
        submit = Submit(self.para)
        submit.get_submitted_out()
        thermo = Thermo(self.para)
        thermo.get_all_thermodat()
      