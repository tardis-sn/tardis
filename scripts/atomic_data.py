"""
=====================
TARDIS Atomic_data module
=====================
created on Feb 29, 2020
"""

import numpy as np
import pandas as pd

class Atomic_Data: 
    """ Will serve as base class for all atomic data contained in 
    'si2_osc_kurucz' in order to provide consistent Tabular data.
    
    Class Members
        ----------
        data : dataframe
            raw unprocessed data
            
        first_table : dataframe
            Contains Energy Levels and statistical weights for Si II
            
        second_table : dataframe
            Oscillator strengths for LLL 
            
    Example to use
         ----------
         
         x = Atomic_Data()
         
         print ( x.first_table )
         
         """
    
    data = pd.read_csv( 'si2_osc_kurucz' , header = None )
    data = data[0] .str. split( expand = True)
    data. replace( "|" , np.nan , inplace = True )
    data ["NaN"] = data .isnull(). sum( axis = 1 )
    
    """
    data conatins two tables.
    Need to collect them separately

    """
    t1 = 0
    t2 = 0
    t = []
    
    for i in range( data.shape[0]-1 ):
        if( abs(data["NaN"][i+1] - data["NaN"][i]) > 2 ):
            if(t2 - t1 > 5):
                t.append([t1,t2])
            t1 = i+1
            t2 = i+1
        else:
            t2 = t2 + 1
    t.append([t1,data.shape[0]])
    
    first_table = data[ t[0][0] : t[0][1]+1 ][ : ]
    first_table. reset_index( inplace = True)
    first_table = first_table. drop( ["NaN"] , axis = 1 )
    
    second_table = data [ t[1][0]+1 : t[1][1]+1] [ : ]
    second_table. reset_index ( inplace = True )
    
    for i in range(second_table.shape[0]):
        if ( second_table["NaN"][i] > 3):
            first_atr , second_atr = second_table[0][i]. split( "]-" )
            first_atr += "]"
            c = "-" + second_atr
            j = 7
            while ( j > 0 ):
                second_table .at[ i , j+1 ] = second_table[j][i]
                j -= 1
            second_table. at[ i , 0 ] = first_atr
            second_table. at[ i , 1 ] = c
    
    for i in list(second_table.columns):
        if(second_table .isnull(). sum( axis = 0 )[i] > second_table.shape[0]//2):
            second_table = second_table. drop( [i] , axis = 1 )
    
    second_table = second_table. drop( ["NaN"] , axis = 1 )
    