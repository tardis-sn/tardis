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
         import atomic
         
         x = Atomic_Data()
         
         print ( x.first_table )
         
         """
    
    data = pd.read_csv( 'si2_osc_kurucz' , header = None )

    """
    data conatins two tables.
    Need to collect them separately

    """

    first_table = data[ 12 : 169 ][ : ]
    first_table. reset_index( inplace = True)
    first_table = first_table[0]. str. split( expand = True )

    second_table = data [ 176 : ] [ : ]
    second_table. reset_index ( inplace = True )
    second_table = second_table[0].str.split( expand = True )
    second_table. replace( "|" , np.nan , inplace = True )
    second_table [10] = second_table .isnull(). sum( axis = 1 )
    
    for i in range(4196):
        if ( second_table[10][i] > 3):
            first_atr , second_atr = second_table[0][i]. split( "]-" )
            first_atr += "]"
            c = "-" + second_atr
            j = 7
            while ( j > 0 ):
                second_table .at[ i , j+1 ] = second_table[j][i]
                j -= 1
            second_table. at[ i , 0 ] = first_atr
            second_table. at[ i , 1 ] = c
            
    second_table = second_table. drop( [7,8,9,10] , axis = 1 )
