from pylab import *
from numpy import *
import scipy.signal as sp
import sympy as sy

def solve_x_t(d):
	x_s = sp.lti([1,d],np.polymul([1,0,2.25], np.polyadd(np.polymul([1, d],[1, d]),[2.25])))
	t_s = np.linspace(0,50,1000)
	t,x = sp.impulse(x_s,None,t_s )
	return t,x

t_1, x_1 = solve_x_t(0.5)
figure()
plot(t_1,x_1)
xlabel(r'$t\ \longrightarrow$', size=20)
ylabel(r'$Time\ Response\ \longrightarrow$', size=20)
grid(True)
savefig('spring1.jpeg')

t_2, x_2 = solve_x_t(0.05)
figure()
plot(t_2,x_2)
xlabel(r'$t\ \longrightarrow$', size=20)
ylabel(r'$Time\ Response\ \longrightarrow$', size=20)
grid(True)
savefig('spring2.jpeg')


H_f = sp.lti([1], [1,0,2.25])
t = linspace(0,50,10000)
figure()
for  i in arange(0, 3, 1):
	f = cos((1.4 + i*0.05)*t)*exp(-0.05*t)
	t0, y, sv = sp.lsim2(H_f, U = f, T = t)
	plot(t0 ,y ,label = r'$frequency(w) = {}$'.format(1.4 + i*0.05))
xlabel(r'$t\ \longrightarrow$', size=20)
ylabel(r'$Plots\ for\ varying\ \omega \longrightarrow$', size=20)
legend(loc = 'upper left')
grid(True)
savefig('spring_v.jpeg')


t = linspace(0,20,2000)
X_0 = sp.lti([1,0,2],[1,0,3,0])
Y_0 = sp.lti([2],[1,0,3,0])
t_1, x_0 = sp.impulse(X_0, None, t)
t_2, y_0 = sp.impulse(Y_0, None, t)
figure()
plot(t_1, x_0, 'g')
plot(t_2, y_0, 'y')
xlabel(r'$t\ \longrightarrow$', size=20)
ylabel(r'$Coupled\ Spring\ Response\ \longrightarrow$', size=20)
grid(True)
savefig('spring_co.jpeg')


L = 1e-6
C = 1e-6
R = 100
H_rlc = sp.lti([1], [L*C, R*C, 1])
w, H_rlc_mod, H_rlc_p = H_rlc.bode()

figure()
semilogx(w, H_rlc_mod)
xlabel(r'$\omega \ \longrightarrow$', size=20)
ylabel(r'$Magnitude\ \longrightarrow$', size=20)
grid(True)
savefig('RLC_mag.jpeg')

figure()
semilogx(w, H_rlc_p)
xlabel(r'$\omega \ \longrightarrow$', size=20)
ylabel(r'$Phase\ \longrightarrow$', size=20)
grid(True)
savefig('RLC_ph.jpeg')


t_rlc = arange(0,30e-6,1e-8)
v_i_st = cos(1e3*t_rlc) - cos(1e6*t_rlc)
t_rlc_2 = arange(0,0.01,1e-7)
v_i_lt = cos(1e3*t_rlc_2) - cos(1e6*t_rlc_2)

t, v_o, sv = sp.lsim(H_rlc, v_i_st, t_rlc)
figure()
plot(t, v_o)
xlabel(r'$t \ \longrightarrow$', size=20)
ylabel(r'$Output\ to\ the\ circuit\ \longrightarrow$', size=20)
grid(True)
savefig('RLC_s.jpeg')

t, v_o, sv = sp.lsim(H_rlc, v_i_lt, t_rlc_2)
figure()
plot(t, v_o)
xlabel(r'$t \ \longrightarrow$', size=20)
ylabel(r'$Output\ to\ the\ circuit\ \longrightarrow$', size=20)
grid(True)
savefig('RLC_l.jpeg')






