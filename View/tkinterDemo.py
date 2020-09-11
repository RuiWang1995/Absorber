import tkinter as tk

from Model.Gases import Gases
from Model.Properties import Properties
from Control.SolvePrEos import SolvePrEos

on_hit = False


def calculated_molar_volume(gas_chosen, temperature, pressure):
    CH3CL = Gases('CH3Cl', 50.49, 412.6, 6680000, 0.151, 1552046.088274, False, 'gas')
    CH4 = Gases('CH4', 16.043, 190.6, 4610000, 0.0115, 1552046.088274, True, 'gas')
    PT = Properties(temperature, pressure)
    if gas_chosen == 'CH3CL':
        solPrEos1 = SolvePrEos(CH3CL, PT)
        return solPrEos1.solve()
    if gas_chosen == 'CH4':
        solPrEos1 = SolvePrEos(CH4, PT)
        return solPrEos1.solve()


window = tk.Tk()
window.title('Peng-Robinson Equation of State')
window.geometry('340x500+30+30')

temperature = tk.StringVar()
temperature_label = tk.Label(window, text='Temperature', bg='white', fg='black', font=('Arial', 12), width=10)
temperature_label.place(x=20, y=30, width=120, height=25)
temperature_text = tk.Entry(window, show=None, font=('Arial', 14), textvariable=temperature)
temperature_text.focus_get()

pressure = tk.StringVar()
pressure_label = tk.Label(window, text='Pressure', bg='white', fg='black', font=('Arial', 12), width=10)
pressure_label.place(x=20, y=70, width=120, height=25)
pressure_text = tk.Entry(window, show=None, font=('Arial', 14), textvariable=pressure)
pressure_text.place(x=150, y=70, width=120, height=25)
pressure_text.focus_get()


var2 = tk.StringVar()
var2.set(('CH4', 'CH3CL'))
gases = tk.Listbox(window, listvariable=var2)
gases.place(x=20, y=110, width=120, height=50)

gas_chosen = 'CH4'

vm = calculated_molar_volume(gas_chosen, temperature.get(), pressure.get())
molar_volume = tk.DoubleVar()


def calculate():
    global on_hit
    if not on_hit:
        on_hit = True
        molar_volume.set(vm)
    else:
        on_hit = False
        molar_volume.set('')


calculate_button = tk.Button(window, text='calculate', width=15, height=2,
                             command=calculate)
calculate_button.place(x=85, y=170, width=80, height=25)

molar_volume_label = tk.Label(window, textvariable=molar_volume,
                              bg='white', fg='black', font=('Arial', 12), width=10)
molar_volume_label.place(x=20, y=205, width=120, height=25)
window.mainloop()
