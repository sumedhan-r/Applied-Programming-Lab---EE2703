from statistics import *               #The required plots in this program
from array import *                    #come one after the other by closing
from numpy import *                    #the previous plot
from pylab import *
import scipy.special as sp
from scipy.linalg import *
from io import StringIO

t = array(zeros(101))                  #Defining empty arrays 
Y = array(zeros(11), dtype=object)
n = array(zeros(10), dtype=object)
sigma = array(zeros(9))


Arr = loadtxt("./fitting.dat")         #Loading the data
for i in range(101):
    t[i] =  Arr[i][0]


def g(t,A,B):                          #Defining the function g(t,A,B)
    return (A*sp.jn(2,t) + B*t)


y = g(t,1.05,-0.105)                   #Writing the data into variables for
for i in range(9):                     #easy plotting and finding noise
  	Y[i+1] = (Arr[:, [i+1]])
  	n[i+1] = [x1 - x2 for (x1,x2) in zip(Y[i+1],y)]
  	sigma[i] = std(n[i+1])

for i in range(9):                     #Plotting the first graph
    k = i+1
    plot(t,Y[k],label = r'$\sigma_{%d}$ = %.2f' %( k, sigma[i]))
plot(t,y,label = 'True value')
xlabel(r'$t\ \longrightarrow$' ,size=20)
ylabel(r'$f(t)\ +\ n\ \longrightarrow$' ,size=20)
legend(loc = 'upper right')
grid(True)
savefig('Bessel.jpeg', bbox_inches = 'tight')

                                       #Plotting the second graph
figure()
errorbar(t[::5],Y[1][::5],sigma[0],fmt='ro',label = r'$Errorbar$')
plot(t,y,label = r'$f(t)$')
xlabel(r'$t\ \longrightarrow$' ,size=20)
legend(loc = 'upper right')
grid(True)
savefig('Bessel_errorbar.jpeg', bbox_inches = 'tight')

                                       #Proving the two vectors to be equal
A0 = float(input('Enter the value of A0 : '))
B0 = float(input('Enter the value of B0 : '))
x = sp.jn(2,t) 
M = c_[x,t]
v = array([A0,B0])
v = v.astype('float64')
ans = dot(M,v)
if array_equal(ans, g(t,A0,B0)) is True:
	print('The two vectors are similar')
else:
	print('The two vectors are not similar')

Alist = arange(0.0,2.1,0.1)            #Creating values of A and B and
Blist = arange(-0.2,0.0,0.01)          #also the meshgrid
Blist = append(Blist,0)
A1, B1 = meshgrid(Alist,Blist)

error = array(zeros((21,21)))          #Finding the mean squared error

for i in range(21):
	for j in range(21):
		for k in range(101):
			y0 = g(t[k],Alist[i],Blist[j])
			error[i][j] = error[i][j] + pow((Arr[k][1] - y0),2) 
	error[i][j] = error[i][j]/101
figure()
co = contour(A1,B1,error)   
clabel(co, inline = True)                #Plotting the mean squared error
plot(1.05, -0.105, 'ro')
annotate(r'$Exact\ Location$', xy=(1.05,-0.105))
xlabel(r'$A\ \longrightarrow$' ,size=20)
ylabel(r'$B\ \longrightarrow$' ,size=20)
savefig('Error.jpeg', bbox_inches = 'tight')

A2 = array(zeros(9))                   #Finding and assigning least square 
B2 = array(zeros(9))                   #values of A and B

for i in range(9):
	A2[i], B2[i] = lstsq(M,Y[i+1])[0]
	A2[i] = A2[i] - 1.05
	B2[i] = B2[i] + 0.105
figure()
plot(sigma, A2, 'ro--', label = r'$A_{err}$') #Plotting the lstsq errors
plot(sigma, B2, 'go--', label = r'$B_{err}$')
xlabel(r'$Noise\ Standard\ Deviation \longrightarrow$' ,size=20)
ylabel(r'$MS\ Error \longrightarrow$' ,size=20)
legend(loc = 'upper left')
grid(True)
savefig('Least_square.jpeg', bbox_inches  = 'tight')

figure()
grid(linestyle = ':')
loglog(sigma, A2, 'ro--', label = r'$A_{err}$') #Plotting th lstsq errors
loglog(sigma, B2, 'go--', label = r'$B_{err}$') #on the log scale

errorbar(sigma, A2, A2, fmt = 'ro', uplims = True, lolims = True, linestyle = ':')
errorbar(sigma, B2, B2, fmt = 'ro', uplims = True, lolims = True, linestyle = ':')

xscale('log')
yscale('log')
xlabel(r'$\sigma_{n}\ \longrightarrow$' ,size=20)
ylabel(r'$MS\ Error \longrightarrow$' ,size=20)
legend(loc = 'upper left')
savefig('Least_square_log.jpeg', bbox_inches = 'tight')













