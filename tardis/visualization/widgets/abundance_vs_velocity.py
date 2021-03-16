import pandas as pd
import matplotlib.pyplot as plt

def plot_abundance_vs_velocity(sim):
    #loading abundance sata in 'df'
    df = (sim.model.abundance).T
    
    #saving the velocity array in 'velocity'
    velocity = sim.model.velocity
    
    #obtaing a map of atomic numbers to element symbol
    #the index of a will be the atomic number and corresponding value will be element symbol
    a = atomic_number_to_element_symbol()
    
    #plotting the graph of abundance of element vs velocity
    
    #for loop iterate and plot abundance of each element at a time
    for index,row in df.T.iterrows():
        plt.plot(velocity[:-1] , df[index] , 'o-' , label=a[index])
    plt.legend()
    plt.title('Abundance of elements v/s Velocity')
    plt.xlabel('Velocity')
    plt.ylabel('Abundance (in Fraction)')
    plt.show()


def atomic_number_to_element_symbol():
    a = [None]*119
    a[1:16] = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P']
    a[16:31] = ['S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
    a[31:45] = ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru']
    a[45:59] = ['Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce']
    a[59:73] = ['Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf']
    a[73:87] = ['Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']
    a[87:101] = ['Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm']
    a[101:114] = ['Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut']
    a[114:119] = ['Fl', 'Uup', 'Lv', 'Uus', 'Uuo']
    return a



