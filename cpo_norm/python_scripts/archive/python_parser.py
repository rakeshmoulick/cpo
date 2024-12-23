"""
Plots the dispersion graph of from the electric field data
Use %reset -f to clear the python workspace
Data File Invoked: processed_results_all.npz
Run as: 
"""
# TO RUN: python3 dispersion_npz.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser


config_path = sys.argv[1]

# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)


NUM_TS = config.getint('time','NUM_TS')
DT_coeff = config.getfloat('diagnostics','DT_coeff')
write_interval = config.get('diagnostics','write_interval')
n0 = config.get('population','density')
mi = config.get('population','massI')
NC = config.get('grid','NC')

print(NUM_TS, type(NUM_TS))
print(DT_coeff, type(DT_coeff))
print(write_interval, type(write_interval))
print(n0, type(n0))
print(mi, type(mi))
print(NC, type(NC))

