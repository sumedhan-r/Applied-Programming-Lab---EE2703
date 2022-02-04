from sympy import *
import pylab as p
import numpy as np
import scipy.signal as sp

Pi = np.pi
t = np.linspace(0,0.01,int(1e6))
s = symbols("s")

def lowpass(R1,R2,C1,C2,G,Vi):
	A = Matrix([[0,0,1,-1/G],
				[-1/(1+s*R2*C2),1,0,0],
				[0,-G,G,1],
				[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
	b = Matrix([0,0,0,-Vi/R1])
	V = A.inv()*b
	return V[3]


def sympy_to_LTI(H):
    num, den = H.as_numer_denom()
    numer = Poly(num, s)
    denom = Poly(den, s)
    numeratorCoeffs = numer.all_coeffs()
    denominatorCoeffs = denom.all_coeffs()
    for i in range(len(numeratorCoeffs)):
        x = float(numeratorCoeffs[i])
        numeratorCoeffs[i] = x
    for j in range(len(denominatorCoeffs)):
        x = float(denominatorCoeffs[j])
        denominatorCoeffs[j] = x
    return sp.lti(numeratorCoeffs, denominatorCoeffs)


Vi = np.heaviside(t,1)*(np.sin(2e3*Pi*t) + np.cos(2e6*Pi*t))

Vo_1 = lowpass(10000,10000,1e-9,1e-9,1.586,1)
print(Vo_1)

f_1 = sympy_to_LTI(Vo_1)
w_1, mag_1, phi_1 = sp.bode(f_1, np.linspace(1, 1e6, int(1e6)))

p.figure()
p.semilogx(w_1, mag_1)
p.title('Magnitude plot for lowpass')
p.xlabel(r'$\omega \ \to$')
p.ylabel(r'$20log(\|H(j\omega)\|)$')
p.grid(True)
p.savefig('lowpass_mag.jpeg')
#p.show()

p.figure()
p.semilogx(w_1, phi_1)
p.title('Phase plot for lowpass')
p.xlabel(r'$\omega \ \to$')
p.ylabel(r'$\angle H(j\omega)$')
p.grid(True)
p.savefig('lowpass_ph.jpeg')
#p.show()

t, vo_1, svec = sp.lsim(f_1, Vi, t)
p.figure()
p.plot(t, vo_1)
p.title('Output of Vi in Lowpass')
p.xlabel(r'$t\ \to$')
p.ylabel(r'$V_o(t)\ \to$')
p.xlim(0,1e-3)
p.grid(True)
p.savefig('lowpass_o.jpeg')
#p.show()

t, x = sp.step(f_1, None, t)
p.figure()
p.plot(t,x)
p.title('Step response of Lowpass filter')
p.xlabel(r'$t\ \to$')
p.ylabel(r'$V_o(t)\ \to$')
p.xlim(0, 1e-3)
p.grid(True)
p.savefig('lowpass_s.jpeg')
#p.show()





def highpass(R1,R3,C1,C2,G,Vi):
	A = Matrix([[0, -1, 0, 1/G],
				[-(s*R3*C2)/(1+s*R3*C2), 0, -1, 0],
				[0, G, -G, 1],
				[-1/R1-s*C1-s*C2, 0, s*C2, 1/R1]])
	b = Matrix([0,0,0,-Vi*s*C1])
	V = A.inv()*b
	return V[3]


Vo_2 = highpass(10000,10000,1e-9,1e-9,1.586,1)
print(Vo_2)

f_2 = sympy_to_LTI(Vo_2)
w_2, mag_2, phi_2 = sp.bode(f_2, np.linspace(1, 1e6, int(1e6)))

p.figure()
p.semilogx(w_2, mag_2)
p.title('Magnitude plot for highpass')
p.xlabel(r'$\omega \ \to$')
p.ylabel(r'$20log(\|H(j\omega)\|)$')
p.grid(True)
p.savefig('highpass_mag.jpeg')
#p.show()

p.figure()
p.semilogx(w_2, phi_2)
p.title('Phase plot for highpass')
p.xlabel(r'$\omega \ \to$')
p.ylabel(r'$\angle H(j\omega)$')
p.grid(True)
p.savefig('highpass_ph.jpeg')
#p.show()

t, vo_2, svec = sp.lsim(f_2, Vi, t)
p.figure()
p.plot(t, vo_2)
p.title('Output of Vi in Highpass')
p.xlabel(r'$t\ \to$')
p.ylabel(r'$V_o(t)\ \to$')
p.xlim(0,1e-4)
p.grid(True)
p.savefig('highpass_o.jpeg')
#p.show()


Vi_lf = np.exp(-t)*np.sin(2*Pi*t) # decay coefficient = 1

Vi_hf = np.exp(-100*t)*np.sin(2e6*Pi*t) # decay coefficient = 100

t, dam_1,svec = sp.lsim(f_2, Vi_lf, np.linspace(0,1e-4,int(1e6)))

p.figure()
p.plot(t, dam_1)
p.title('Circuit_2_damped_input_lowfreq')
p.xlabel(r'$t\ \to$')
p.ylabel(r'$V_o(t)\ \to$')
p.grid(True)
p.savefig('highpass_dl.jpeg')
#p.show()

t, dam_2,svec = sp.lsim(f_2, Vi_hf, np.linspace(0,1e-4,int(1e6)))

p.figure()
p.plot(t, dam_2)
p.title('Circuit_2_damped_input_highfreq')
p.xlabel(r'$t\ \to$')
p.ylabel(r'$V_o(t)\ \to$')
p.grid(True)
p.savefig('highpass_dh.jpeg')
#p.show()


t, x = sp.step(f_2, None, t)
p.figure()
p.plot(t,x)
p.title('Step response of Highpass filter')
p.xlabel(r'$t\ \to$')
p.ylabel(r'$V_o(t)\ \to$')
p.xlim(0, 1e-4)
p.grid(True)
p.savefig('highpass_s.jpeg')
#p.show()




