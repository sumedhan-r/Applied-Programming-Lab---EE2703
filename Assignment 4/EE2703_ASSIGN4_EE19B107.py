import sys
import re
from numpy import *
from  pylab import *
from scipy.integrate import *
#from math import *

t = []
inp = linspace(-2*pi, 4*pi, num=400, endpoint=False)

def f(x):
	t.clear()
	if size(x) == 1:
		return(exp(x))
	else:	
		for i in range(size(x)):
			t.append(exp(x[i]))
		return(t)

def g(x):
	t.clear()
	if size(x) == 1:
		return(cos(cos(x)))
	else:
		for i in range(size(x)):
			t.append(cos(cos(x[i])))
		return(t)

figure(num = 1, clear = True)
semilogy(inp[0: int(400/3)], f(inp[int(400/3): int(800/3)]), 'g', label = r'$e^{x}$')
semilogy(inp[int(400/3): int(800/3)], f(inp[int(400/3): int(800/3)]), 'g')
semilogy(inp[int(800/3)+1: ], f(inp[int(400/3): int(800/3)]), 'g')
xlim(-2*pi, 4*pi)
xlabel(r'$x\ \longrightarrow$', size=20)
ylabel(r'$f(x)\ \longrightarrow$', size=20)
legend( loc = 'upper right')
grid(True)
savefig('exp_actual.jpeg', bbox_inches = 'tight')

figure(num = 2, clear = True)
plot(inp, g(inp), 'orange', label = r'$cos(cos(x))$')
xlim(-2*pi, 4*pi)
xlabel(r'$x\ \longrightarrow$', size=20)
ylabel(r'$g(x)\ \longrightarrow$', size=20)
legend( loc = 'upper right')
grid(True)
savefig('cos_cos_actual.jpeg', bbox_inches = 'tight')

def u_1(x,k):
	return(f(x)*cos(k*x))
def v_1(x,k):
	return(f(x)*sin(k*x))
def u_2(x,k):
	return(g(x)*cos(k*x))
def v_2(x,k):
	return(g(x)*sin(k*x))

n1 = linspace(0, 25, num=26, endpoint=True)
n2 = linspace(1, 25, num=25, endpoint=True)
vf = vg = array(zeros(51))
v1 = array(zeros(26))
v2 = array(zeros(25))

v1[0] = vf[0] = quad(f, 0, 2*pi)[0]/(2*pi)
for i in range(25):
	v1[i+1] = vf[2*i+1] = quad(u_1, 0, 2*pi, args=(i+1))[0]/(2*pi)
	v2[i] = vf[2*(i+1)] = quad(v_1, 0, 2*pi, args=(i+1))[0]/(2*pi)
	if vf[2*i+1] < 0:
		v1[i] = -1*vf[2*(i+1)]
	elif vf[2*(i+1)] < 0:
		v2[i] = -1*vf[2*(i+1)]

figure(num = 3)
semilogy(n1, v1, 'ro', label = r'$A$' )
semilogy(n2, v2, 'ro', label = r'$B$')
xlabel(r'$N\ \longrightarrow$', size = 20)
ylabel(r'$Coefficients\ Magnitude\ \longrightarrow$', size = 20)
legend(loc = 'upper right')
grid(True)
savefig('exp_semi.jpeg', bbox_inches = 'tight')

figure(num = 4)
loglog(n1, v1, 'ro', label = r'$A$' )
loglog(n2, v2, 'ro', label = r'$B$')
xlabel(r'$N\ \longrightarrow$', size = 20)
ylabel(r'$Coefficients\ Magnitude\ \longrightarrow$', size = 20)
legend(loc = 'upper right')
grid(True)
savefig('exp_log.jpeg', bbox_inches = 'tight')

v1[0] = vg[0] = quad(g, 0, 2*pi)[0]/(2*pi)
for i in range(25):
	v1[i+1] = vg[2*i+1] = quad(u_2, 0, 2*pi, args=(i+1))[0]/(2*pi)
	v2[i] = vg[2*(i+1)] = quad(v_2, 0, 2*pi, args=(i+1))[0]/(2*pi)
	if vg[2*i+1] < 0:
		v1[i] = -1*vg[2*(i+1)]
	elif vg[2*(i+1)] < 0:
		v2[i] = -1*vg[2*(i+1)]

figure(num = 5)
semilogy(n1, v1, 'ro', label = r'$A$' )
semilogy(n2, v2, 'ro', label = r'$B$')
xlabel(r'$N\ \longrightarrow$', size = 20)
ylabel(r'$Coefficients\ Magnitude\ \longrightarrow$', size = 20)
legend(loc = 'upper right')
grid(True)
savefig('cos_cos_semi.jpeg', bbox_inches = 'tight')

figure(num = 6)
loglog(n1, v1, 'ro', label = r'$A$' )
loglog(n2, v2, 'ro', label = r'$B$')
xlabel(r'$N\ \longrightarrow$', size = 20)
ylabel(r'$Coefficients\ Magnitude\ \longrightarrow$', size = 20)
legend(loc = 'upper right')
grid(True)
savefig('cos_cos_log.jpeg', bbox_inches = 'tight')

x = linspace(0, 2*pi, num=401)
x = x[:-1] 
b = f(x)   
A = zeros((400,51))     
A[:,0] = 1              
for k in range(1,26):
	A[:,2*k-1] = cos(k*x) 
	A[:,2*k] = sin(k*x)   
cf = lstsq(A, b, rcond = None)[0] 
for i in range(25):
	v1[i+1] = cf[2*i+1]
	v2[i] = cf[2*(i+1)]
v1[0] = cf[0]
figure(num = 7)
plot(n1, v1, 'go', label = r'$A$')
plot(n2, v2, 'go', label = r'$B$')
xlabel(r'$N\ \longrightarrow$', size=20)
ylabel(r'$Coefficients\ \longrightarrow$', size=20)
grid(True)
savefig('exp_lstsq.jpeg', bbox_inches = 'tight')

b = g(x)            
cg = lstsq(A, b, rcond = None)[0] 
for i in range(25):
	v1[i+1] = cg[2*i+1]
	v2[i] = cg[2*(i+1)]
v1[0] = cg[0]
figure(num = 1, clear = True)
plot(n1, v1, 'go', label = r'$A$')
plot(n2, v2, 'go', label = r'$B$')
xlabel(r'$N\ \longrightarrow$', size=20)
ylabel(r'$Coefficients\ \longrightarrow$', size=20)
grid(True)
savefig('cos_cos_lstsq.jpeg', bbox_inches = 'tight')

diff_f = [abs(a-b) for a, b in zip(cf, vf)]
diff_g = [abs(a-b) for a, b in zip(cg, vg)]
n = linspace(0, 50, num=51, endpoint=True)
s1 = std(diff_f)
s2 = std(diff_g)
diff_f.sort()
diff_g.sort()
print("The largest deviations are %f and %f " %(diff_f[-1], diff_g[-1]))
print("The standard deviations are %f and %f" %(s1, s2))

pf = dot(A,cf)
pg = dot(A,cg)

figure(num = 1, clear = True)
semilogy(inp[0: int(400/3)], f(inp[int(400/3): int(800/3)]), 'g', label = r'$True\ Values$')
semilogy(inp[int(400/3): int(800/3)], f(inp[int(400/3): int(800/3)]), 'g')
semilogy(inp[int(800/3)+1: ], f(inp[int(400/3): int(800/3)]), 'g')
semilogy(inp[0: int(400/3)], pf[int(400/3): int(800/3)], 'orange', label = r'$Estimated\ Values$')
semilogy(inp[int(400/3): int(800/3)], pf[int(400/3): int(800/3)], 'orange')
semilogy(inp[int(800/3)+1: ], pf[int(400/3): int(800/3)], 'orange')
xlim(-2*pi, 4*pi)
xlabel(r'$x\ \longrightarrow$', size = 20)
ylabel(r'$f(x)\ \longrightarrow$', size = 20)
legend(loc = 'upper right')
grid(True)
savefig('exp_estimated.jpeg', bbox_inches = 'tight')

figure(num = 2, clear = True)
plot(inp, g(inp), 'g', label = r'$True\ Values$')
plot(inp, pg, 'orange', label = r'$Estimated\ Values$')
xlim(-2*pi, 4*pi)
xlabel(r'$x\ \longrightarrow$', size = 20)
ylabel(r'$g(x)\ \longrightarrow$', size = 20)
legend(loc = 'upper right')
grid(True)
savefig('cos_cos_estimated.jpeg', bbox_inches = 'tight')









	
