# %%
from numpy import *
import pandas as pd
import tkinter as tk

# Dye dictionaries

af_350 = {
    'Name':'AF 350',
    'MW': 410,
    'Lambda max': 346,
    'Emission max': 442,
    'Extinction coefficient': 19000,
    'Correction factor': 0.19
}
af_488 = {
    'Name':'AF 488',
    'MW': 885,
    'Lambda max': 494,
    'Emission max': 519,
    'Extinction coefficient': 71000,
    'Correction factor': 0.11
}
af_532 = {
    'Name':'AF 532',
    'MW': 724,
    'Lambda max': 530,
    'Emission max': 554,
    'Extinction coefficient': 81000,
    'Correction factor': 0.09
}

all_dyes = [af_350, af_488, af_532]

test_dict = {'AF 350':af_350, 'AF 488':af_488}


dye_options = [dye['Name'] for dye in all_dyes]

print(dye_options)

print(test_dict['AF 350']['MW'])




