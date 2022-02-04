import sys

CIRCUIT = '.circuit'
END = '.end'

try:
    if(len(sys.argv)<2):
    sys.exit('Invalid number of arguements')
    name = sys.argv[1]
    with open('name') as f:
      lines = f.readlines()
          
        
    for i,str in enumerate(lines):
      try:
        if str.split()[0]==CIRCUIT:
          start=i
        elif str.split()[0]==END:
          end=i
      except IndexError:
          continue
               
    try:
      ckt=lines[start+1:end]
    except NameError:
          print('Invalid format\n')
          exit()
          
    for l in reversed(ckt):
          tokens=l.split()
          i=len(tokens)
          reverse = ' '
          
    for j in range(len(tokens)):
        if tokens[j][0]=='#':
          index = j
          break
            
    for j in range(index):
          reverse1 = tokens[index-(j+1)]
          reverse = reverse + reverse1 + ' '
               
    print(reverse)
      
except IOError:
   print('Invalid file\n')
      
      
                
            
     
