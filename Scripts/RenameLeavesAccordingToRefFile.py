#! /usr/bin/env python
#renames leaves in a tree (or any other type of file)

import math
import re
import sys
print (sys.argv)

argv=sys.argv[1:]
file=argv[0]
out=argv[1]
ref=argv[2]

try:
    f=open(file, 'r')
except IOError:
    print ("Unknown file: ",file)
    sys.exit()

try:
    fout=open(out, 'w')
except IOError:
    print ("Unknown file: ",out)
    sys.exit()

try:
    fref=open(ref, 'r')
except IOError:
    print ("Unknown file: ",ref)
    sys.exit()


corresp = dict()
for l in fref:
    li=l.split("\t")
    if l.startswith("accession"):
        pass
    elif len(li)>=2:
        sp_name =""
        # sp_name is going to be accession_numer%taxid%taxonomy
        if "d__Bacteria" in li[-6]:
            sp_name = li[0] + "%" + li[51] + "%"+li[-6].strip().replace("d__Bacteria;","").replace(";","%").replace(":","_").replace("(","%").replace(")","%") + "%" + li[-5].strip().replace("(","%").replace(")","%")
        elif "d__Bacteria" in li[-5]:
            sp_name = li[0] + "%" + li[51] + "%"+ li[-5].strip().replace("d__Bacteria;","").replace(";","%").replace(":","_").replace("(","%").replace(")","%") + "%" + li[-4].strip().replace("(","%").replace(")","%")
        elif "d__Archaea" in li[-5]:
            sp_name = li[0] + "%" + li[51] + "%"+ li[-5].strip().replace("d__Archaea;","").replace(";","%").replace(":","_").replace("(","%").replace(")","%") + "%" + li[-4].strip().replace("(","%").replace(")","%")
        else:
            print("\n\n\t\tPROBLEM LINE: " + l)
            exit(-1)
        corresp[ li[0] ] = sp_name
        print (li[0] +" Corresponds to "+ sp_name)

print("ref file read.")

#regex = re.compile("'\w -:\.'")
regex = re.compile("'*'")
for l in f:
#    ma=regex.match( l)
#    print(ma)
    l=re.sub("'[a-zA-Z0-9:;\-_\. ]*'", "", l)
    print("AFTER regex:")
    print(l)
    i = 0
    for k in corresp.keys():
#        print k +" indeed corresponds to "+ corresp[k]
        l = l.replace(k, corresp[k])
        i = i+1
        if i % 1000 == 0:
            print (i)

    fout.write(l)
f.close()
fout.close()
fref.close()
