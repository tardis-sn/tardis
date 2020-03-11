from DataUnits import DataUnit
import pytest as pt

def test(O):
	t=[]
	for o in O:
		t.append(o.get_data_unit(['col1']))
	return t

def test_ans(O):
	assert test(O)==[['kg'],['m'],['P']]

if __name__=="__main__":
	o1=DataUnit({"col1":[1,2,2,4],"col2":[2,3,4,5],"col3":[5,2,3,4]},units=["kg","m","A"])
	o2=DataUnit({"col1":[1,2,2,4],"col2":[2,3,4,5],"col3":[5,2,3,4]},units=["m","s","N"])
	o3=DataUnit({"col1":[1,2,2,4],"col2":[2,3,4,5],"col3":[5,2,3,4]},units=["P","N","A"])
	obj=[o1,o2,o3]
	test_ans(obj)
	print("Passed")



