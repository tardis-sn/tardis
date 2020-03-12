from DataUnits import DataUnit
import pytest as pt



def test_ans():
	o1=DataUnit({"col1":[1,2,2,4],"col2":[2,3,4,5],"col3":[5,2,3,4]},units=["kg","m","A"])
	o2=DataUnit({"col1":[1,2,2,4],"col2":[2,3,4,5],"col3":[5,2,3,4]},units=["m","s","N"])
	o3=DataUnit({"col1":[1,2,2,4],"col2":[2,3,4,5],"col3":[5,2,3,4]},units=["P","N","A"])
	obj=[o1,o2,o3]
	t=[]
	for o in obj:
		t.append(o.get_data_unit(['col1']))

	assert t==[['kg'],['m'],['P']]





