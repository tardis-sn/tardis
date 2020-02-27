#importing neccessary library
import pandas as pd

#reading data
data = pd.read_csv('si2_osc_kurucz',header=None)

"""
data conatins two tables.
Need to collect them separately

"""

#first_table
#selecting relevant rows
first_table = data[12:169][:]
#resuming the indexing 
first_table.reset_index(inplace=True)
#spliting the columns
first_table=first_table[0].str.split(expand=True)


#second_table
#selecting relevant rows
second_table = data[176:][:]
#resuming the indexing 
second_table.reset_index(inplace=True)
#spliting the columns
second_table=second_table[0].str.split(expand=True)
#selecting relevant columns
second_table = second_table.drop([7,8,9],axis=1)