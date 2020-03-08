"""
=====================
TARDIS datatable module
=====================
created on Mar 2, 2020
"""

import pandas as pd
import numpy as np
import astropy.units as u
import uuid
from pandas import HDFStore
import hashlib

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
            
        units : Series of astropy.units
            It contains the units where index of units coincides with
            index of column attribute. It is the metadata and hence remains 
            conserved during manipulation. 
            
    Example to use
         ----------
         df = pd.DataFrame( 
                data = pd.np.random.randint(0, 100, (10, 5)) , 
                columns = list('ABCED') )
         
         datatable = DataTable (df, units=[u.meter, u.second, u.kg, u.meter/u.second, u.Ohm]) 
         
         
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
        if len((self.units)) > max (0 ,ind) :
            return self.units[ind]
        else :
            print( "The unit for the following attribute is not set." )
            
    
    def scalar_mul(self, const, unit, attr):
        
        """
        Multiplies a constant having unit to entire datatable
        
        Parameters
        ----------
        const: int
            value of the scalar
            
        unit: astropy
            unit of the scalar
        
        attr: array of string or int
            list of columns to which scalar is multiplied
        
        """        
        new_units = []
        for i in range(len(self.units)):
            if( list(self.columns) [i] in attr):
                new_units.append(self.units[i]*unit)
            else:
                new_units.append(self.units[i])
        
        for i in attr:
            if(i not in self.columns):
                print("Wrong Attributes entered")
                break
            self[i] = const*self[i]
            
        self.units = pd.Series(new_units)   
        
        
    def __private_dot(self,series):
        """
        It is a private method that computes Matrix Multiplication.
        
        Parameters
        ----------
        series : Series, DataFrame or array-like
            The other object to compute the matrix product with.
        
        Returns
        -------
        Series or DataFrame
        Matrix product    
        
        """
        if isinstance(series, (pd.Series, pd.DataFrame)):
            common = self.columns.union(series.index)
            if len(common) > len(self.columns) or len(common) > len(series.index):
                raise ValueError("matrices are not aligned")

            left = self.reindex(columns=common, copy=False)
            right = series.reindex(index=common, copy=False)
            lvals = left.values
            rvals = right.values
        else:
            left = self
            lvals = self.values
            rvals = np.asarray(series)
            if lvals.shape[1] != rvals.shape[0]:
                raise ValueError(
                    f"Dot product shape mismatch, {lvals.shape} vs {rvals.shape}"
                )

        if isinstance(series, pd.DataFrame):
            return self._constructor(
                np.dot(lvals, rvals), index=left.index, columns=series.columns
            )
        elif isinstance(series, pd.Series):
            return pd.Series(np.dot(lvals, rvals), index=left.index)
        elif isinstance(rvals, (np.ndarray, pd.Index)):
            result = np.dot(lvals, rvals)
            if result.ndim == 2:
                return self._constructor(result, index=left.index)
            else:
                return pd.Series(result, index=left.index)
        else:  # pragma: no cover
            raise TypeError(f"unsupported type: {type(series)}")
            
      
    def dot(self, series, unit, attr):
        
        """
        Multiplies a constant having unit to some columns of datatable
        
        Parameters
        ----------
        series: Series or array
            values of the series
            
        unit: astropy.units
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
                new_units.append(self.units[i]*unit)
            else:
                new_units.append(self.units[i])
        
        flag=1
        for i in attr:
            if(i not in self.columns):
                print("Wrong Attributes entered")
                flag=0
                break
        if(flag): 
            dataframe = self.__private_dot(series)
            result = DataTable(dataframe)
            result.units = pd.Series(new_units)
            return result
    
    def add(self, series, unit, attr):
        
        """
        Add a column in datatable 
        
        Parameters
        ----------
        series: Series or array
            values of the series
            
        unit: astropy.units
            unit of the scalar
            
        attr: string or int 
            name of the new attribute
            
        """
        flag = 1
        for i in self.columns:
            if i == attr:
                flag = 0
        
        if flag == 0 :
            print("The attribute name already exists")
            print("Kindly change the attribute name")
            print("Or use replace function if you want to change the value")
        else:
            self [attr] = series
            self . units = list(self.units)
            self. units. append( unit ) 
            self.units = pd.Series(self.units)
            
        
    def delete(self, attr):
        
        """
        Deletes a column
        
        Parameters
        ----------
        attr: string or int 
            name of the new attribute
            
        """
        flag = 1
        index = 0
        for i in self.columns:
            if i == attr:
                flag = 0
            else:
               index+=1
        
        if flag == 1 :
            print("The attribute does not exists")
            print("Kindly change the attribute name")
        
        else:
            self.units = list(self.units)
            del self.units[list(self.columns).index(attr)]
            self.drop(columns = attr) 
            self.units = pd.Series(self.units)
        
    def replace(self, series, unit, attr):
        
        """
        Replace a column in datatable 
        
        Parameters
        ----------
        series: Series or array
            values of the series
            
        unit: astropy.units
            unit of the scalar
            
        attr: string or int 
            name of the new attribute that needs to be replaced
            
        """
        flag = 1
        for i in self.columns:
            if i == attr:
                flag = 0
    
                
        if flag == 1 :
            print("The attribute does not exists")
            print("Kindly change the attribute name")
            print("Or use add function if you want to change the value")
        else:
            ind = list(self.columns).index(attr)
            self [attr] = series
            self . units = list(self.units)
            self. units[ind] = unit  
            self.units = pd.Series(self.units)
            
    def join(self, othertable):
        
        """
        Performs a join of self with othertable and returns the same.
        Parameters
        ----------
        othertable: object of class DataTable
            The table with which the join is to be taken
            
        Returns
        ----------
        object of class DataTable
        DataTable that is the result of join of self with othertable
            
        """
        columns1 = self.columns
        columns2 = othertable.columns
        common_columns = []
        for i in columns1:
            for j in columns2:
                if i==j :
                    common_columns.append(i)
                    
        resDf = pd.merge(self, othertable, on=common_columns)
        print(resDf)
        result = DataTable(resDf)
        res_units = []
        for i in list(resDf.columns):
            if(i in self.columns):
                ind = list(self.columns).index(i)
                res_units.append(list(self.units)[ind])
            elif(i in othertable.columns):
                ind = list(othertable.columns).index(i)
                res_units.append(list(othertable.units)[ind])
        result.units = pd.Series(res_units)
        return result
        
    def writeToHDF5(self, path):
        """
        Write the content of current datatable into a HDF% file.
        
        Parameters
        ----------
        path: String
            path of the HDF5 file
            
        """
        with HDFStore(path) as store:
            store.put("units", self.units)
            md5_hash = hashlib.md5()
            for key in store.keys():
                tmp = np.ascontiguousarray(store[key].values.data)
                md5_hash.update(tmp)
            
            uuid1 = uuid.uuid1().hex

            print("Signing AtomData: \nMD5: {}\nUUID1: {}".format(
                md5_hash.hexdigest(), uuid1))

            store.root._v_attrs['md5'] = md5_hash.hexdigest().encode('ascii')
            store.root._v_attrs['uuid1'] = uuid1.encode('ascii')
            store.close()
