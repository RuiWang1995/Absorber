import tkinter as tk

from Model.Gases import Gases
from Model.Properties import Properties
from Control.SolDepart import solDepart
from Control.SolvePrEos import SolvePrEos
from Control.ThermoIdeal import ThermoIdeal

fields = ('Temperature_1', 'Temperature_2', 'Pressure_1', 'Pressure_2', 'Name(UpperCase)', 'Molecular Weight', 'Critical Temperature',
          'Critical Pressure', 'Acentric Factor', 'Dissociation Energy', 'Is linear (True/False)',
          'Phase_1 (gas/liquid)', 'Phase_2 (gas/liquid)', 'Molar Volume', 'Molar Volume_2', 'H_real - H_ideal',
          'H_real - H_ideal_2', 'S_real - S_ideal', 'S_real - S_ideal_2', 'G_real - G_ideal', 'G_real - G_ideal_2',
          'A_ideal gas', 'A_ideal gas_2', 'U_ideal gas', 'U_ideal gas_2', 'Cv', 'Cv_2', 'S_ideal gas', 'S_ideal gas_2',
          'H_ideal gas', 'H_ideal gas_2', 'G_ideal gas', 'G_ideal gas_2',
          'Delta_H', 'Delta_S', 'Delta_G')
unit = ['K', 'K', 'Pa', 'Pa', '-', 'g/mol', 'K',
        'Pa', '-', 'J/mol', '-',
        '-', '-', 'm^3/mol', 'm^3/mol', 'J/mol', 'J/mol', 'J/mol/K', 'J/mol/K', 'J/mol', 'J/mol',
        'J/mol', 'J/mol', 'J/mol', 'J/mol', 'J/mol/K', 'J/mol/K', 'J/mol/K', 'J/mol/K', 'J/mol', 'J/mol', 'J/mol', 'J/mol',
        'J/mol', 'J/mol/K', 'J/mol']


def calculated_molar_volume(entries):
    # get the parameters
    temperature = float(entries['Temperature_1'].get())
    temperature2 = float(entries['Temperature_2'].get())
    pressure = float(entries['Pressure_1'].get())
    pressure2 = float(entries['Pressure_2'].get())
    name = entries['Name(UpperCase)'].get()
    mw = float(entries['Molecular Weight'].get())
    tc = float(entries['Critical Temperature'].get())
    pc = float(entries['Critical Pressure'].get())
    af = float(entries['Acentric Factor'].get())
    do = float(entries['Dissociation Energy'].get())
    if entries['Is linear (True/False)'].get() == 'True':
        is_linear = True
    else:
        is_linear = False
    gas_or_liquid = entries['Phase_1 (gas/liquid)'].get()
    gas_or_liquid2 = entries['Phase_2 (gas/liquid)'].get()
    # period rate:
    gas_chosen = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid)
    gas_chosen2 = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid2)
    PT = Properties(temperature, pressure)
    PT2 = Properties(temperature2, pressure2)
    # molar volume
    pr_eos = SolvePrEos(gas_chosen, PT)
    pr_eos2 = SolvePrEos(gas_chosen2, PT2)
    vm = pr_eos.solve()
    vm2 = pr_eos2.solve()
    # vm = ("%8.2f" % monthly).strip()
    entries['Molar Volume'].delete(0, tk.END)
    entries['Molar Volume'].insert(0, vm)
    entries['Molar Volume_2'].delete(0, tk.END)
    entries['Molar Volume_2'].insert(0, vm2)


def departure_function(entries):
    temperature = float(entries['Temperature_1'].get())
    temperature2 = float(entries['Temperature_2'].get())
    pressure = float(entries['Pressure_1'].get())
    pressure2 = float(entries['Pressure_2'].get())
    name = entries['Name(UpperCase)'].get()
    mw = float(entries['Molecular Weight'].get())
    tc = float(entries['Critical Temperature'].get())
    pc = float(entries['Critical Pressure'].get())
    af = float(entries['Acentric Factor'].get())
    do = float(entries['Dissociation Energy'].get())
    if entries['Is linear (True/False)'].get() == 'True':
        is_linear = True
    else:
        is_linear = False
    gas_or_liquid = entries['Phase_1 (gas/liquid)'].get()
    gas_or_liquid2 = entries['Phase_2 (gas/liquid)'].get()
    gas_chosen = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid)
    gas_chosen2 = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid2)
    PT = Properties(temperature, pressure)
    PT2 = Properties(temperature2, pressure2)
    # Departure Function
    sol_depart = solDepart(gas_chosen, PT)
    sol_depart2 = solDepart(gas_chosen2, PT2)
    deltaH = sol_depart.deltaH()
    deltaS = sol_depart.deltaS()
    deltaG = sol_depart.deltaG()
    deltaH2 = sol_depart2.deltaH()
    deltaS2 = sol_depart2.deltaS()
    deltaG2 = sol_depart2.deltaG()
    entries['H_real - H_ideal'].delete(0, tk.END)
    entries['H_real - H_ideal'].insert(0, deltaH)
    entries['S_real - S_ideal'].delete(0, tk.END)
    entries['S_real - S_ideal'].insert(0, deltaS)
    entries['G_real - G_ideal'].delete(0, tk.END)
    entries['G_real - G_ideal'].insert(0, deltaG)
    entries['H_real - H_ideal_2'].delete(0, tk.END)
    entries['H_real - H_ideal_2'].insert(0, deltaH2)
    entries['S_real - S_ideal_2'].delete(0, tk.END)
    entries['S_real - S_ideal_2'].insert(0, deltaS2)
    entries['G_real - G_ideal_2'].delete(0, tk.END)
    entries['G_real - G_ideal_2'].insert(0, deltaG2)


def ideal_thermodynamic_properties(entries):
    temperature = float(entries['Temperature_1'].get())
    temperature2 = float(entries['Temperature_2'].get())
    pressure = float(entries['Pressure_1'].get())
    pressure2 = float(entries['Pressure_2'].get())
    name = entries['Name(UpperCase)'].get()
    mw = float(entries['Molecular Weight'].get())
    tc = float(entries['Critical Temperature'].get())
    pc = float(entries['Critical Pressure'].get())
    af = float(entries['Acentric Factor'].get())
    do = float(entries['Dissociation Energy'].get())
    if entries['Is linear (True/False)'].get() == 'True':
        is_linear = True
    else:
        is_linear = False
    gas_or_liquid = entries['Phase_1 (gas/liquid)'].get()
    gas_or_liquid2 = entries['Phase_2 (gas/liquid)'].get()
    gas_chosen = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid)
    gas_chosen2 = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid2)
    PT = Properties(temperature, pressure)
    PT2 = Properties(temperature2, pressure2)
    # Departure Function
    thermodynamic_ideal = ThermoIdeal(gas_chosen, PT)
    thermodynamic_ideal2 = ThermoIdeal(gas_chosen2, PT2)
    aig = thermodynamic_ideal.Aig()
    uig = thermodynamic_ideal.Uig()
    cv = thermodynamic_ideal.Cv()
    sig = thermodynamic_ideal.Sig()
    hig = thermodynamic_ideal.Hig()
    gig = thermodynamic_ideal.Gig()
    aig2 = thermodynamic_ideal2.Aig()
    uig2 = thermodynamic_ideal2.Uig()
    cv2 = thermodynamic_ideal2.Cv()
    sig2 = thermodynamic_ideal2.Sig()
    hig2 = thermodynamic_ideal2.Hig()
    gig2 = thermodynamic_ideal2.Gig()
    entries['A_ideal gas'].delete(0, tk.END)
    entries['A_ideal gas'].insert(0, aig)
    entries['U_ideal gas'].delete(0, tk.END)
    entries['U_ideal gas'].insert(0, uig)
    entries['Cv'].delete(0, tk.END)
    entries['Cv'].insert(0, cv)
    entries['S_ideal gas'].delete(0, tk.END)
    entries['S_ideal gas'].insert(0, sig)
    entries['H_ideal gas'].delete(0, tk.END)
    entries['H_ideal gas'].insert(0, hig)
    entries['G_ideal gas'].delete(0, tk.END)
    entries['G_ideal gas'].insert(0, gig)
    entries['A_ideal gas_2'].delete(0, tk.END)
    entries['A_ideal gas_2'].insert(0, aig2)
    entries['U_ideal gas_2'].delete(0, tk.END)
    entries['U_ideal gas_2'].insert(0, uig2)
    entries['Cv_2'].delete(0, tk.END)
    entries['Cv_2'].insert(0, cv2)
    entries['S_ideal gas_2'].delete(0, tk.END)
    entries['S_ideal gas_2'].insert(0, sig2)
    entries['H_ideal gas_2'].delete(0, tk.END)
    entries['H_ideal gas_2'].insert(0, hig2)
    entries['G_ideal gas_2'].delete(0, tk.END)
    entries['G_ideal gas_2'].insert(0, gig2)


def delta_H(entries):
    temperature = float(entries['Temperature_1'].get())
    temperature2 = float(entries['Temperature_2'].get())
    pressure = float(entries['Pressure_1'].get())
    pressure2 = float(entries['Pressure_2'].get())
    name = entries['Name(UpperCase)'].get()
    mw = float(entries['Molecular Weight'].get())
    tc = float(entries['Critical Temperature'].get())
    pc = float(entries['Critical Pressure'].get())
    af = float(entries['Acentric Factor'].get())
    do = float(entries['Dissociation Energy'].get())
    if entries['Is linear (True/False)'].get() == 'True':
        is_linear = True
    else:
        is_linear = False
    gas_or_liquid = entries['Phase_1 (gas/liquid)'].get()
    gas_or_liquid2 = entries['Phase_2 (gas/liquid)'].get()
    gas_chosen = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid)
    gas_chosen2 = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid2)
    PT = Properties(temperature, pressure)
    PT2 = Properties(temperature2, pressure2)
    # Departure Function
    thermodynamic_ideal1 = ThermoIdeal(gas_chosen, PT)
    thermodynamic_ideal2 = ThermoIdeal(gas_chosen2, PT2)
    prEos1 = solDepart(gas_chosen, PT)
    prEos2 = solDepart(gas_chosen2, PT2)
    deltaH = round(thermodynamic_ideal2.Hig() + prEos2.deltaH() - thermodynamic_ideal1.Hig() - prEos1.deltaH(), 5)
    entries['Delta_H'].delete(0, tk.END)
    entries['Delta_H'].insert(0, deltaH)


def delta_S(entries):
    temperature = float(entries['Temperature_1'].get())
    temperature2 = float(entries['Temperature_2'].get())
    pressure = float(entries['Pressure_1'].get())
    pressure2 = float(entries['Pressure_2'].get())
    name = entries['Name(UpperCase)'].get()
    mw = float(entries['Molecular Weight'].get())
    tc = float(entries['Critical Temperature'].get())
    pc = float(entries['Critical Pressure'].get())
    af = float(entries['Acentric Factor'].get())
    do = float(entries['Dissociation Energy'].get())
    if entries['Is linear (True/False)'].get() == 'True':
        is_linear = True
    else:
        is_linear = False
    gas_or_liquid = entries['Phase_1 (gas/liquid)'].get()
    gas_or_liquid2 = entries['Phase_2 (gas/liquid)'].get()
    gas_chosen = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid)
    gas_chosen2 = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid2)
    PT = Properties(temperature, pressure)
    PT2 = Properties(temperature2, pressure2)
    # Departure Function
    thermodynamic_ideal1 = ThermoIdeal(gas_chosen, PT)
    thermodynamic_ideal2 = ThermoIdeal(gas_chosen2, PT2)
    prEos1 = solDepart(gas_chosen, PT)
    prEos2 = solDepart(gas_chosen2, PT2)
    deltaS = round(thermodynamic_ideal2.Sig() + prEos2.deltaS() - thermodynamic_ideal1.Sig() - prEos1.deltaS(), 5)
    entries['Delta_S'].delete(0, tk.END)
    entries['Delta_S'].insert(0, deltaS)


def delta_G(entries):
    temperature = float(entries['Temperature_1'].get())
    temperature2 = float(entries['Temperature_2'].get())
    pressure = float(entries['Pressure_1'].get())
    pressure2 = float(entries['Pressure_2'].get())
    name = entries['Name(UpperCase)'].get()
    mw = float(entries['Molecular Weight'].get())
    tc = float(entries['Critical Temperature'].get())
    pc = float(entries['Critical Pressure'].get())
    af = float(entries['Acentric Factor'].get())
    do = float(entries['Dissociation Energy'].get())
    if entries['Is linear (True/False)'].get() == 'True':
        is_linear = True
    else:
        is_linear = False
    gas_or_liquid = entries['Phase_1 (gas/liquid)'].get()
    gas_or_liquid2 = entries['Phase_2 (gas/liquid)'].get()
    gas_chosen = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid)
    gas_chosen2 = Gases(name, mw, tc, pc, af, do, is_linear, gas_or_liquid2)
    PT = Properties(temperature, pressure)
    PT2 = Properties(temperature2, pressure2)
    # Departure Function
    thermodynamic_ideal1 = ThermoIdeal(gas_chosen, PT)
    thermodynamic_ideal2 = ThermoIdeal(gas_chosen2, PT2)
    prEos1 = solDepart(gas_chosen, PT)
    prEos2 = solDepart(gas_chosen2, PT2)
    deltaG = round(thermodynamic_ideal2.Gig() + prEos2.deltaG() - thermodynamic_ideal1.Gig() - prEos1.deltaG(), 5)
    entries['Delta_G'].delete(0, tk.END)
    entries['Delta_G'].insert(0, deltaG)


def makeform(root, fields):
    entries = {}
    photo = tk.PhotoImage(file="Methane-3D-balls.png")
    imgLabel = tk.Label(root, image=photo)
    i = 0
    for i in range(len(fields)):
        row = tk.Frame(root)
        lab = tk.Label(row, width=16, text=fields[i] + ": ", anchor='w')
        lab_indentation = tk.Label(row, width=20)
        lab_for_unit = tk.Label(row, width=7, text=unit[i], anchor='w')
        lab_for_unit2 = tk.Label(row, width=4, text=unit[i], anchor='w')
        ent = tk.Entry(row)
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        if i in [0, 2, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]:
            lab.pack(side=tk.LEFT)
            lab_for_unit.pack(side=tk.LEFT)
            ent.pack(side=tk.LEFT, expand=tk.YES)
            lab2 = tk.Label(row, width=16, text=fields[i+1] + ": ", anchor='w')
            lab2.pack(side=tk.LEFT)
            lab_for_unit2.pack(side=tk.LEFT)
            ent2 = tk.Entry(row)
            ent2.pack(side=tk.LEFT, expand=tk.YES, fill=tk.X)
            entries[fields[i]] = ent
            entries[fields[i+1]] = ent2
            if i == 0:
                ent.insert(0, "298.15")
                ent2.insert(0, "298.15")
            if i == 2:
                ent.insert(0, "101325")
                ent2.insert(0, "101325")
            if i == 11:
                ent.insert(0, "gas")
                ent2.insert(0, "gas")
            continue
        if not (i in [1, 3, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32]):
            if i == 4:
                ent.insert(0, "CH4")
            if i == 5:
                ent.insert(0, "16.043")
            if i == 6:
                ent.insert(0, "190.6")
            if i == 7:
                ent.insert(0, "4610000")
            if i == 8:
                ent.insert(0, "0.0115")
            if i == 9:
                ent.insert(0, "1640546.1")
            if i == 10:
                ent.insert(0, "False")
            if 4 <= i <= 12 or i >= 23:
                lab_indentation.pack(side=tk.LEFT)
                lab.pack(side=tk.LEFT)
                lab_for_unit.pack(side=tk.LEFT)
                ent.pack(side=tk.LEFT)
                if i == 12:
                    imgLabel.pack(side=tk.LEFT)
            else:
                lab.pack(side=tk.LEFT)
                lab_for_unit.pack(side=tk.LEFT)
                ent.pack(side=tk.LEFT)
            entries[fields[i]] = ent
        i += 1
    return entries


if __name__ == '__main__':
    root = tk.Tk()
    root.title('Peng-Robinson Equation of State')
    ents = makeform(root, fields)
    b2 = tk.Button(root, text='Molar Volume',
                   command=(lambda e=ents: calculated_molar_volume(e)))
    b2.pack(side=tk.LEFT, padx=5, pady=5)
    b1 = tk.Button(root, text='Departure Function',
                   command=(lambda e=ents: departure_function(e)))
    b1.pack(side=tk.LEFT, padx=5, pady=5)
    b4 = tk.Button(root, text='Ideal Property',
                   command=(lambda e=ents: ideal_thermodynamic_properties(e)))
    b4.pack(side=tk.LEFT, padx=5, pady=5)

    b5 = tk.Button(root, text='Delta_H',
                   command=(lambda e=ents: delta_H(e)))
    b5.pack(side=tk.LEFT, padx=5, pady=5)
    b6 = tk.Button(root, text='Delta_S',
                   command=(lambda e=ents: delta_S(e)))
    b6.pack(side=tk.LEFT, padx=5, pady=5)
    b7 = tk.Button(root, text='Delta_G',
                   command=(lambda e=ents: delta_G(e)))
    b7.pack(side=tk.LEFT, padx=5, pady=5)
    b3 = tk.Button(root, text='Quit', command=root.quit)
    b3.pack(side=tk.LEFT, padx=5, pady=5)

    root.mainloop()
