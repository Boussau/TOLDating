#! /usr/bin/env python
#renames leaves in a tree (or any other type of file)

import math
import re
import sys
print (sys.argv)

argv=sys.argv[1:]
file=argv[0]
out=argv[1]

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

#regex = re.compile("'\w -:\.'")
regex = re.compile("'*'")
for l in f:
#    ma=regex.match( l)
#    print(ma)
    l=re.sub("'[a-zA-Z0-9:;\-_\. ]*'", "", l)
    print("AFTER regex:")
    print(l)
    i = 0
    fout.write(l)
f.close()
fout.close()
