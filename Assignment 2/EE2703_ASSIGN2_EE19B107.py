import sys
import numpy as np
from math import *

count1 = 0                                  #Creating empty variables
count2 = 0
freq = 0
list = np.array(['None']*100)
print(pi)

try:                                        #Extracting data from sys
    if(len(sys.argv)!=2):
       sys.exit('Invalid number of arguements')
    name = sys.argv[1]
    with open(name) as f:
        lines = f.readlines()
          
    class comp:                             #Creating object of the class that
         cname = ' '                        #stores the component details 
         comp_n1 = ' '
         comp_n2 = ' '
         stype = ' '
         k = 0
         value1 = 0
         value2 = 0
         
         def __init__(self,i):
                 self.i = i
                 
         def assign(self,str):
              m = len(str.split())
              if len(str.split())==3:
                  self.stype = 'ac'
                  self.cname = str.split()[1]
                  freq = float(str.split()[2])
              else:
                  self.cname = str.split()[0]
                  self.comp_n1 = str.split()[1]     
                  self.comp_n2 = str.split()[2]
                  if len(str.split()) > 4:
                    self.stype = str.split()[3]
                    if self.stype=='ac':
                       self.k = int(str.strip()[1])
                       self.value1 = int(str.split()[m-2])
                       self.value2 = float(str.split()[m-1])
                    elif self.stype=='dc':
                       self.value1 = float(str.split()[m-1])                    
                  else:
                      self.value1 = float(str.split()[3])

    
                 
    for i,str in enumerate(lines):       #Removing unwanted lines from input  
      try:                               #that doesn't contain circuit details
          if str.split()[0]=='.circuit':
            start=i
          elif str.split()[0]=='.end':
            end=i
      except IndexError:
          continue

    
    ckt=lines[start+1:end]              #Creating alist of circuit components
    try:
      if lines[end+1].split()[0]=='.ac':
        ckt = ckt + lines[end+1]
    except IndexError:
        pass
except NameError:
    print('Invalid format\n')
    exit()

    for i,str in enumerate(ckt):       #Assigning objects to the list ckt and
        list[i] = comp(i)              #finding the values of n and k
        list[i].assign(str)
        if list[i].cname.strip()[0]=='V':
          count2 = count2 + 1
        for j in range(i):
          if list[i].comp_n1==list[j].comp_n1:
              continue
        count1 = count1 + 1
               

    hash = {'GND':0}                   #Creating dictionary

    c = count1 + count2

    for j in range(len(count1)-1):
        hash['j+1'] = j+1
    for k in range(len(count2)):
        hash['V(k+1)'] = j+k+1
        j = j+1


    M = np.array([zeros(c) for i in range(c)], dtype=complex) #Creating the arrays

    b = np.array([zeros(c)], dtype=complex)

    for i in range(c):                  #Assigning values to the arrays
        if list[i].cname.strip()[0]=='V':
            r = cmath.rect(list[i].value1/2,list[i].value2)
            b[hash['V(list[i].k)']] = r
        elif list[i].cname.strip()[0]=='I':
            r = cmath.rect(list[i].value1/2,list[i].value2)
            b[hash['list[i].comp_n1']] = -r
            b[hash['list[i].comp_n2']] = r
        else:
          for k in range(c):
            p = hash['list[i].comp_n1']; q = hash['list[i].comp_n2']; a = hash['list[i].k']
            if list[i].cname.strip()[0]=='R':
              M[p][q] = -1/list[i].value1 + M[p][q]
              M[q][p] = -1/list[i].value1 + M[q][p]
              M[p][p] = 1/list[i].value1 + M[p][p]
              M[q][q] = 1/list[i].value1 + M[q][q]
            elif list[i].cname.strip()[0]=='C':
              M[p][q] = complex(0,(list[i].value1)*freq*2*pi) + M[p][q]
              M[q][p] = complex(0,(list[i].value1)*freq*2*pi) + M[q][p]
              M[p][p] = complex(0,-(list[i].value1)*freq*2*pi) + M[p][p]
              M[q][q] = complex(0,-(list[i].value1)*freq*2*pi) + M[q][q]
            elif list[i].cname.strip()[0]=='L':
              M[p][q] = complex(0,-(1/list[i].value1)*freq*2*pi) + M[p][q]
              M[q][p] = complex(0,-(1/list[i].value1)*freq*2*pi) + M[q][p]
              M[p][p] = complex(0,(1/list[i].value1)*freq*2*pi) + M[p][p]
              M[q][q] = complex(0,(1/list[i].value1)*freq*2*pi) + M[q][q]
            elif list[i].cname.strip()[0]=='V':
              M[p][a] = 1 + M[p][a]
              M[q][a] = -1 + M[q][a]
              M[a][p] = 1 + M[a][p]
              M[a][q] = -1 + M[a][q]

    x = np.linalg.solve(M,b)           #Solving the matrix equation

    for i in range(count1):            #Printing the output
      print('V{} = x[i]' .format(i+1))
    for i in range(count2):
      print('I{} = x[count1+i]' .format( i+1))


except IOError:
  print('Invalid file\n')
      
             