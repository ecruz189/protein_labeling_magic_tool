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
af_546 = {
    'Name':'AF 546',
    'MW': 1160,
    'Lambda max': 554,
    'Emission max': 570,
    'Extinction coefficient': 104000,
    'Correction factor': 0.12
}
af_555 = {
    'Name':'AF 555',
    'MW': 1260,
    'Lambda max': 555,
    'Emission max': 565,
    'Extinction coefficient': 150000,
    'Correction factor': 0.08
}
af_568 = {
    'Name':'AF 568',
    'MW': 792,
    'Lambda max': 577,
    'Emission max': 603,
    'Extinction coefficient': 91300,
    'Correction factor': 0.46
}
af_594 = {
    'Name':'AF 594',
    'MW': 820,
    'Lambda max': 590,
    'Emission max': 617,
    'Extinction coefficient': 73000,
    'Correction factor': 0.56
}
af_633 = {
    'Name':'AF 633',
    'MW': 1200,
    'Lambda max': 632,
    'Emission max': 647,
    'Extinction coefficient': 100000,
    'Correction factor': 0.55
}
af_647 = {
    'Name':'AF 647',
    'MW': 1300,
    'Lambda max': 650,
    'Emission max': 668,
    'Extinction coefficient': 239000,
    'Correction factor': 0.03
}
af_660 = {
    'Name':'AF 660',
    'MW': 1100,
    'Lambda max': 663,
    'Emission max': 690,
    'Extinction coefficient': 132000,
    'Correction factor': 0.1
}
af_680 = {
    'Name':'AF 680',
    'MW': 1150,
    'Lambda max': 679,
    'Emission max': 702,
    'Extinction coefficient': 184000,
    'Correction factor': 0.05
}
fluorescein_EX = {
    'Name':'Fluorescein-EX',
    'MW': 591,
    'Lambda max': 494,
    'Emission max': 518,
    'Extinction coefficient': 68000,
    'Correction factor': 0.2
}
oregon_green_488 = {
    'Name':'Oregon Green 488',
    'MW': 509,
    'Lambda max': 496,
    'Emission max': 524,
    'Extinction coefficient': 70000,
    'Correction factor': 0.12
}
pacific_blue = {
    'Name':'Pacific Blue',
    'MW': 339,
    'Lambda max': 410,
    'Emission max': 455,
    'Extinction coefficient': 30000,
    'Correction factor': 0.2
}

all_dyes = [af_350, af_488, af_532, af_546, af_555, af_568, af_594, af_633, af_647, af_660, af_680, fluorescein_EX, oregon_green_488, pacific_blue]

test_dict = {'AF 350':af_350, 'AF 488':af_488, 'AF 532':af_532, 'AF 546':af_546, 
            'AF 555':af_555, 'AF 568':af_568, 'AF 594':af_594, 'AF 633':af_633, 
            'AF 647':af_647, 'AF 660':af_660, 'AF 680':af_680, 'Fluorescein-EX':fluorescein_EX,
            'Oregon Green 488':oregon_green_488, 'Pacific Blue':pacific_blue}

##################################  Function to calculate Degree of Labeling (DOL) and protein concentration  ##################################

def calculate_dol():
    default_antibody = antibody_boolvar.get()
    if default_antibody:
        prot_mw = 150000
        prot_ext_coeff = 203000
    else:
        prot_mw = int(protein_mw.get(1.0, "end-1c"))
        prot_ext_coeff = int(protein_ext_coeff.get(1.0, "end-1c"))

    dye = clicked.get()
    dye_mw = test_dict[dye]['MW']
    dye_ext_coeff = test_dict[dye]['Extinction coefficient']
    dye_cf = test_dict[dye]['Correction factor']
    prot_abs_280 = float(absorbance_280.get(1.0, "end-1c"))
    dye_abs = float(absorbance_dye.get(1.0, "end-1c"))

    prot_conc = (prot_abs_280 - (dye_abs * dye_cf)) / prot_ext_coeff
    prot_conc_uM = prot_conc * 1000000
    prot_conc_mgmL = prot_conc * prot_mw

    dol = dye_abs / (dye_ext_coeff * prot_conc)

    dol_result = tk.Label(
            master=frame_six,
            text='DOL:  {:.3f}\n\nProtein concentration:  {:.3f} µM  |  {:.3f} mg/mL\n\nProt MW:  {} Da   |   Prot ext coeff:  {} M-1 cm-1   |   Dye:   {}\nDye MW:  {} Da   |   Dye ext coeff:  {} M-1 cm-1   |   Dye CF (280 nm):  {}\nProt Abs:  {}   |   Dye Abs:  {}'.format(dol, prot_conc_uM, prot_conc_mgmL, prot_mw, prot_ext_coeff, dye, dye_mw, dye_ext_coeff, dye_cf, prot_abs_280, dye_abs),
            width=65,
            borderwidth=2,
            height=7,
            bg='white',
            pady=10
    )

    dol_result.grid(row=1, column=0)

##################################  GUI  ##################################

root = tk.Tk()

root.columnconfigure(0, weight=1, minsize=150)

##################################  First frame  ##################################

logo_frame = tk.Frame(master=root, borderwidth=1, padx=5, pady=5)
logo_frame.grid(row=0, column=0, sticky='w')

logo = tk.Label(
        text=" Esteban's protein labeling magic tool ",
        fg='white',
        bg='#505050',
        width=67,
        height=2,
        master=logo_frame,
        relief=tk.RIDGE,
        borderwidth=2
)

logo.pack()

##################################  First frame (dye selection and antibody default parameters)  ##################################

frame_one = tk.Frame(master=root, borderwidth=1, padx=5, pady=10)
frame_one.grid(row=1, column=0)

clicked = tk.StringVar()
clicked.set("Dye")

dye_options = [dye['Name'] for dye in all_dyes]

initial_settings = tk.OptionMenu(frame_one, clicked, *dye_options)

initial_settings.grid(row=0, column=0, padx=5)

antibody_boolvar = tk.BooleanVar()
antibody_boolvar.set(False)

antibody_option = tk.Checkbutton(
        master=frame_one,
        text='Use antibody default parameters                 \n(MW 150.000 da | molar ext. coeff. 203.000)',
        variable=antibody_boolvar,
        onvalue=True,
)

antibody_option.grid(row=0, column=1, sticky='w', padx=20)

##################################  Second frame (Absorbance values) ##################################

frame_two = tk.Frame(master=root, borderwidth=1, padx=5, pady=5, relief=tk.RIDGE)
frame_two.grid(row=2, column=0, sticky='w')

absorbance_280_lbl = tk.Label(
        master=frame_two,
        text='Absorbance 280 nm',
        bg='#4682B4',
        fg='white',
        borderwidth=2,
        width=24,
        height=1
)
absorbance_280_lbl.grid(row=0, column=0)

absorbance_280 = tk.Text(
        master=frame_two,
        width=6,
        height=1
)
absorbance_280.grid(row=0, column=1)

spacer_frame_two = tk.Label(
        master=frame_two,
        text='   ',
        borderwidth=2,
        width=2,
        height=1
)
spacer_frame_two.grid(row=0, column=2)

absorbance_dye_lbl = tk.Label(
        master=frame_two,
        text='Absorbance dye lambda max',
        bg='#4682B4',
        fg='white',
        borderwidth=2,
        width=24,
        height=1
)
absorbance_dye_lbl.grid(row=0, column=3)

absorbance_dye = tk.Text(
        master=frame_two,
        width=6,
        height=1
)
absorbance_dye.grid(row=0, column=4)

##################################  Spacer frame  ##################################

frame_three = tk.Frame(master=root, borderwidth=1, padx=5)
frame_three.grid(row=3, column=0, sticky='w')

optional_parameters_lbl = tk.Label(master=frame_three, text='', height=1)
optional_parameters_lbl.grid(row=0, column=0)

##################################  Third frame (Optional protein parameters) ##################################

frame_three = tk.Frame(master=root, borderwidth=1, padx=5)
frame_three.grid(row=4, column=0, sticky='w')

optional_parameters_lbl = tk.Label(master=frame_three, text='Optional protein parameters (if antibody default parameters not selected):', height=1, pady=5)
optional_parameters_lbl.grid(row=0, column=0)

##################################  Fourth frame (Optional protein parameters input) ##################################

frame_four = tk.Frame(master=root, borderwidth=1, padx=5, pady=5, relief=tk.RIDGE)
frame_four.grid(row=5, column=0)

protein_mw_lbl = tk.Label(
        master=frame_four,
        text='Protein MW (da)',
        bg='#4682B4',
        fg='white',
        borderwidth=2,
        width=24,
        height=1
)
protein_mw_lbl.grid(row=0, column=0)

protein_mw = tk.Text(
        master=frame_four,
        width=6,
        height=1
)
protein_mw.grid(row=0, column=1)

spacer_frame_four = tk.Label(
        master=frame_four,
        text='   ',
        borderwidth=2,
        width=2,
        height=1
)
spacer_frame_four.grid(row=0, column=2)

protein_ext_coeff_lbl = tk.Label(
        master=frame_four,
        text='Protein molar ext. coeff.',
        bg='#4682B4',
        fg='white',
        borderwidth=2,
        width=24,
        height=1
)
protein_ext_coeff_lbl.grid(row=0, column=3)

protein_ext_coeff = tk.Text(
        master=frame_four,
        width=6,
        height=1
)
protein_ext_coeff.grid(row=0, column=4)

##################################  Spacer frame  ##################################

frame_five = tk.Frame(master=root, borderwidth=1, padx=5)
frame_five.grid(row=6, column=0, sticky='w')

optional_parameters_lbl = tk.Label(master=frame_five, text='', height=1)
optional_parameters_lbl.grid(row=0, column=0)

##################################  Sixth frame (Calculate button) ##################################

frame_six = tk.Frame(master=root, borderwidth=1, padx=5, pady=5)
frame_six.grid(row=7, column=0)

calculate_btn = tk.Button(
        text='Calculate DOL and protein concentration',
        width=35,
        height=1,
        bg='#505050',
        fg='white',
        master=frame_six,
        relief=tk.RAISED,
        borderwidth=2,
        command=calculate_dol
)
calculate_btn.grid(row=0, column=0)

##################################  Main loop  ##################################

root.mainloop()
