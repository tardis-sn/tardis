import numpy as np
import pandas as pd

class Parser:
    """
    Class to get data from 'si2_osc_kurucz' in a tabular form.

    Example:
            temp = Parser
            print(temp.energy_levels)
            print(temp.oscillator_strengths)

    """

    def find_header_and_row( data, text1, text2):
        #Returns values from header.
        with open(data) as d:
            for l2 in d:            
                if text1 in l2:
                    break

        n = int(l2.split()[0])
    
        with open(data) as d:
            rows = 0
            for l1 in d:        

                rows += 1
                if text2 in l1:
                    break
        
        return(n),(rows-1)

    file = '/home/atul/Tardis/atomic_data_15nov16/atomic/SIL/II/16sep15/si2_osc_kurucz'
    args = {}
    
    #Table 1
    args['header'] = None
    args['delim_whitespace'] = True
    
    # Search the total number of energy levels in header and the first energy level = 0.000
    args['nrows'], args['skiprows'] = find_header_and_row(file, "Number of energy levels", "0.000")

    energy_levels = pd.read_csv(file, **args)
    energy_levels.columns = ['Energy Level', 'g', 'E(cm^-1)', '10^15 Hz', 'eV', 'Lam(A)', 'ID', 'ARAD', 'C4', 'C6']

    
    #Table 2
    #  Search number of transitions in the second header and the line where oscillator strengths table starts
    args['nrows'],args['skiprows']  = find_header_and_row(file, "Number of transitions", "Transition")
    args['skiprows'] = args['skiprows'] +1
    
    oscillator_strengths = pd.read_csv(file, **args)

    args['delim_whitespace'] = False
    args['sep'] = '(?<=[^E])-(?:[ ]{1,})?|(?<!-)[ ]{2,}[-,\|]?'

    oscillator_strengths = pd.read_csv(file, **args)
    oscillator_strengths.columns = ['State A', 'State B', 'f', 'A', 'Lam(A)', 'i', 'j', 'Lam(obs)', '% Acc', '?']    


temp = Parser
print(temp.energy_levels)
print(temp.oscillator_strengths)






