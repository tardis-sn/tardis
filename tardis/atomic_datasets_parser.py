from bs4 import BeautifulSoup
import urllib2
import sys
reload(sys)
sys.setdefaultencoding('utf-8')

url="http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=html"

content=urllib2.urlopen(url).read()

soup=BeautifulSoup(content, "lxml")

#f=open("atom.txt","w")
#f.write(soup.prettify())

tag=soup.find_all("td")

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
def has_colspan(tag):
	if tag.has_attr('colspan'):
		return True
	else:
		return False
list1=[]
f1=open("data.txt","w")

for child in tag:
	if (has_colspan(child)==False):
		alphabets=child.string
		if (alphabets==None):
			continue
		elif(alphabets.find(']')):
			f1.write("\n")
		if(atomic_number(child)==True):
			alphabets=child.string
			alphabets.strip()
			f1.write("Atomic number:"+alphabets)
		elif(atom_name(child)==True):
			alphabets=child.string
			if(alphabets[1].isupper()==True):
				f1.write("Atom name:"+alphabets+"\n")
			"""child=child.find_next("td")
			alphabets=child.string
			if (alphabets==None):
				continue
			else:
				f1.write("Mass Number:"+alphabets+"\n")
			child=child.find_next("td")
			alphabets=child.string
			if (alphabets==None):
				continue
			else:
				f1.write("Relative Atomic Mass:"+alphabets+"\n")
			child=child.find_next("td")
			alphabets=child.string
			if (alphabets==None):
				continue
			else:
				f1.write("Isotopic Composition:"+alphabets+"\n")
			child=child.find_next("td")
			alphabets=child.string
			if (alphabets==None):
				continue
			else:
				f1.write("Standard Atomic Weight:"+alphabets+"\n")
			child1=child	"""		
		else:
			alphabets=child.string
			alphabets.strip()
			if (alphabets==None) and c:
				continue
			else:
				alphabets=alphabets.replace(' ','')
				alphabets.strip()
				f1.write(alphabets)
	else:
		f1.write('\n\n')
			

#f.close()
f1.close()
