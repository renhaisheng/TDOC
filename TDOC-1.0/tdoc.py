#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

from tdoc._tdoc import TDOC



def main():
    
    # Provide arguments.
    parser = argparse.ArgumentParser(usage='%(prog)s smiles.txt -p parameters.py')
    parser.add_argument('smiles_file',  default='smiles.txt', help='Input SMILES file')
    parser.add_argument('-p', '--para', default='parameters.py', help='Input parameters file')
    
    # Obtain arguments.
    args = parser.parse_args()
    input_file=args.smiles_file
    para_file=args.para
    
    # run procedure.
    work = TDOC(input_file, para_file)
    work.execute()



if __name__ == '__main__':
    main()
