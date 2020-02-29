"""
=====================
TARDIS data_table module
=====================
created on Mar 1, 2020

"""

import numpy as np
import pandas as pd

class DataTable:

    """ Class that can be used to keep track of units associated with the
    column attributes. Also it can be used to preserve this unit's data 
    while copying in another object of same class.
    
    Class Members
        ----------
        data : dataframe
            It is dataframe that actually stores the information
            in tabular form
            
        units : array
            It contains the units where index of units coincides with
            index of column attribute
            
        attribute_unit : dictionary
            Dictionary having key attribute name and unit as value 
    
    Example to use
         ----------
         df = pd.DataFrame( 
                data = pd.np.random.randint(0, 100, (10, 5)) , 
                columns = list('ABCED') )
         dt1 = DataTable(df,['m', 's', 't', 'm/s', 'kg'])
         
         
    Example to copy
         ----------
         dt2 = DataTable (dt1.data , dt1.units)
         
         """
    
     
    def __init__(self, data, units):
        if( len(list( data.columns )) == len(units) ):
            self.data = data
            self.units = units
            columns = list( self.data.columns )
            self.unit = {}
            for i in range( len(units) ):
                a = columns[i]
                self.unit[a] = units[i]
        else:
            print("Invalid Input")
            print(" The dimension of units and number of attributes didn't match.")
    
    
    def get_unit(self,attr):
        
        """
        Return unit of a particular attribute.

        Parameters
        ----------

        attr: string or int
            name of the attribute
            
        Returns: string
            unit of the attribute
        """
        return self.unit[attr]    
  