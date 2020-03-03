"""
=====================
TARDIS datatable module
=====================
created on Mar 2, 2020
"""

import pandas as pd
import numpy as np

class DataTable(pd.DataFrame):

    """ Subclass of Dataframe. This class extends pandas can be used to 
    keep track of units associated with the column attributes. Also it 
    can be used to preserve this unit's data while copying in another 
    object of same class.
    
    Class Members
        ----------
        data : dataframe
            It is dataframe that actually stores the information
            in tabular form. It can be simply accessed by the name of the
            object.
            
        units : Series
            It contains the units where index of units coincides with
            index of column attribute. It is the metadata and hence remains 
            conserved during manipulation. 
            
    Example to use
         ----------
         df = pd.DataFrame( 
                data = pd.np.random.randint(0, 100, (10, 5)) , 
                columns = list('ABCED') )
         
         datatable = DataTable (df, units=['m', 's', 't', 'm/s', 'kg']) 
         
         
    Example to set units afterwards
         ----------
         datatable.units = pd.Series(['m', 's', 't', 'm/s', 'kg'])
         
         """

    _metadata = ['units']


    @property
    def _constructor(self):
        return DataTable


    def __init__(self, *args, **kwargs):
        self.units = pd.Series(kwargs.pop("units", None))
        super().__init__(*args, **kwargs)


    def copy(self):
        """
        Creates a copy of the given class
        Returns
        --------
        Object of DataTable class
        
        """
        x = DataTable(self)
        x.units = self.units
        return x


    def get_unit(self, attr):

        """
        Return unit of a particular attribute.
        Parameters
        ----------
        attr: string or int
            name of the attribute
            
        Returns
        ----------
        unit of the attribute
        
        """
        ind = list (self.columns).index(attr)
        if len(self.units) > max (0 ,ind) :
            if len (self.units[ind]) > 0 :
                return self.units[ind]
            else :
                print( "The unit for the following attribute is not set." )
        else :
            print( "The unit for the following attribute is not set." )
            
    
    def scalar_mul(self, const, unit, attr):
        
        """
        Multiplies a constant having unit to entire datatable
        
        Parameters
        ----------
        const: int
            value of the scalar
            
        unit: string
            unit of the scalar
        
        attr: array of string or int
            list of columns to which scalar is multiplied
        
        """        
        new_units = []
        for i in range(len(self.units)):
            if( list(self.columns) [i] in attr):
                new_units.append(self.units[i]+unit)
            else:
                new_units.append(self.units[i])
        
        for i in attr:
            if(i not in self.columns):
                print("Wrong Attributes entered")
                break
            self[i] = const*self[i]
            
        self.units = pd.Series(new_units)   
        
        
    def series_dot_product(self, series, unit, attr):
        
        """
        Multiplies a constant having unit to some columns of datatable
        
        Parameters
        ----------
        series: Series or array
            values of the series
            
        unit: string
            unit of the scalar
            
        attr: array of string or int
            list of columns to which series is multiplied
            
        Returns
        ----------
        object of class DataTable
        
        """
        new_units = []
        for i in range(len(self.units)):
            if(list(self.columns)[i] in attr):
                new_units.append(self.units[i]+unit)
            else:
                new_units.append(self.units[i])
        
        if self.shape[0] != len(list(series)):
            print("The series dimension is not correct")
        
        else:
            data2 = []
            for i in attr:
                if(i not in self.columns):
                    print("Wrong Attributes entered")
                    break
                else:
                    res = 0
                    for j in range(self.shape[0]):
                        res += self[i][j]*series[j]
                    data2.append(res)
        
            dataframe = pd.DataFrame(data = [data2] ,columns=self.columns)
            result = DataTable(dataframe)
            result.units = pd.Series(new_units)
            return result
        
