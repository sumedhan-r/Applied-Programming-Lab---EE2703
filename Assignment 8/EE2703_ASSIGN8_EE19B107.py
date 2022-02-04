from pylab import*
from numpy.fft import *

t=linspace(-4*pi,4*pi,513)
t=t[:-1]
y=(1+0.1*cos(t))*cos(10*t)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-15,15])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig('eg.jpeg')
#show()

y=sin(t)**3
Y=fftshift(fft(y))/512.0
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\sin^{3}\left(t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig('sin.jpeg')
#show()

y=cos(t)**3
Y=fftshift(fft(y))/512.0
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\cos^{3}\left(t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig('cos.jpeg')
#show()

y=cos(20*t + 5*cos(t))
Y=fftshift(fft(y))/512.0
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-30,30])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\cos\left(20t+5\cos\left(t\right)\right)$")
grid(True)
subplot(2,1,2)
ii = where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'ro',lw=2)
xlim([-30,30])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig('cos_d.jpeg')
#show()



t=linspace(-65536*pi,65536*pi,8388609)
t=t[:-1]
y=exp(-0.5*t**2)
Y=fftshift(fft(y))/8388608.0
w=linspace(-128,128,8388609)
w=w[:-1]
figure()
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\exp\left(-t^{2}/2\right)$")
xlabel(r"$\omega$",size=16)
grid(True)
savefig('gauss.jpeg')
#show()