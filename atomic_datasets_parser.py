from bs4 import BeautifulSoup
import urllib2
import sys
import re
reload(sys)
sys.setdefaultencoding('utf-8')

url="http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=html"

content=urllib2.urlopen(url).read()

soup=BeautifulSoup(content, "lxml")

#f=open("atom.txt","w")
#f.write(soup.prettify())

tags=soup.find_all("tr")

def atomic_number(tag):
	if (tag.has_attr('rowspan') and tag.has_attr('valign') and tag.has_attr('align') and tag["align"]=="right"):
		return True
	else:
		return False 
def atom_name(tag):
	if ((tag.has_attr('rowspan') and tag.has_attr('valign') and tag.has_attr('align') and tag["align"]=="center") or (tag.has_attr('align') and tag["align"]=="center")):
		return True
	else:
		return False 
def mass_num(tag):
	if(tag.has_attr('align') and tag["align"]=="right"):
		return True
	else:
		return False
def has_colspan(tag):
	if tag.has_attr('colspan'):
		return True
	else:
		return False

f1=open("data3.txt","w")
flag=0
flag1=0
beg=0
for line in tags:
	#tag=line.findChildren("td")
	list1=""
	for child in line.find_all("td"):
		alphabets=child.string
		if (alphabets==None):
			continue
			
		if(atomic_number(child)==True):
			alphabets=child.string
			flag+=1
			beg+=1
			if(flag>=2):
				list1=list1[:-5]+"\n"+"Atomic_Number:"+alphabets
				beg=0
			else:
				list1=list1[:-3]+"Atomic_Number:"+alphabets
				beg=0
		elif(atom_name(child)==True):
			alphabets=child.string
			if(alphabets[1].isupper()==True):
				list1+="\n"+"Atom_Name:"+alphabets
				
		elif(mass_num(child)==True):
			alphabets=child.string
			list1+="\n"+alphabets		
		else:
			alphabets=child.string
			alphabets.strip()
			
			if(alphabets==None):
				continue
			else:
				for sib in child.next_siblings:
					alph=sib.string
					if(alph !=None):
						try:
							if(sib.has_attr("td")==True):
								list1+="\n\t"
						except AttributeError:
							pass
				
				alphabets=alphabets.replace(' ','')
				alphabets.strip()
				listnum=re.findall(r'\d+\.\d+|\d+', alphabets)
				a1="".join(listnum)
				flag1+=1
				if(flag1>=2):
					list1+=a1+" "
				else:
					list1+=a1
		
	f1.write(list1)
	#f1.write("\n"+list1)	
		

f1.close()
import pandas as pd
with open('data3.txt') as f:
    table = pd.read_table(f, sep=' ', index_col=0, header=None, names=['A','B','C','D','E','F','G','H','I','J','K'],
                          lineterminator='\n')
    table=table.fillna(" ")
print table
f.close()
 
