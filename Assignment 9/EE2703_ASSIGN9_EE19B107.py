from pylab import*
from numpy.fft import *
import sys
import mpl_toolkits.mplot3d.axes3d as p3


def wnd(n):
	N = n.size
	window = zeros(N)
	window = 0.54 + 0.46*cos((2*pi*n)/(N-1))
	return fftshift(window)

n=arange(256)

t=linspace(-4*pi,4*pi,257)
t=t[:-1]
dt=t[1]-t[0]
fmax=1/dt
y=sin(sqrt(2)*t)
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/256.0
w=linspace(-pi*fmax,pi*fmax,257)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b-',lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega\ \to$",size=16)
grid(True)
savefig('sin_root.jpeg')
#show()

y=sin(sqrt(2)*t)
y=y*wnd(n)
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/256.0
w=linspace(-pi*fmax,pi*fmax,257)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b-',lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega\ \to$",size=16)
grid(True)
savefig('sin_root_wnd.jpeg')
#show()

y=cos(0.86*t)**3
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/256.0
w=linspace(-pi*fmax,pi*fmax,257)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b-',lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\cos^{3}\left(\omega_{o}t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega\ \to$",size=16)
grid(True)
savefig('cos_3.jpeg')
#show()

y=cos(0.86*t)**3
y=y*wnd(n)
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/256.0
w=linspace(-pi*fmax,pi*fmax,257)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b-',lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\cos^{3}\left(\omega_{o}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega\ \to$",size=16)
grid(True)
savefig('cos_3_wnd.jpeg')
#show()




def estimateWandD(w, wo, Y, do, pow=2):
    wEstimate = sum(abs(Y)**pow * abs(w))/sum(abs(Y)**pow) # weighted average
    print("wo = {:.03f}\t\two (Estimated) = {:.03f}".format(wo, wEstimate))

    t = linspace(-pi,pi, 129)[:-1]
    y = cos(wo*t + do)

    c1 = cos(wEstimate*t)
    c2 = sin(wEstimate*t)
    A = c_[c1, c2]
    vals = lstsq(A, y, rcond = None)[0]
    dEstimate = arctan2(-vals[1], vals[0])
    print("do = {:.03f}\t\tdo (Estimated) = {:.03f}".format(do, dEstimate))

wo = 1.35
d = pi/2

t = linspace(-pi, pi, 129)[:-1]
trueCos = cos(wo*t + d)
fmax = 1.0/(t[1]-t[0])
n = arange(128)
y = trueCos.copy()*wnd(n)
y = fftshift(y)
Y = fftshift(fft(y))/128.0
w = linspace(-pi*fmax, pi*fmax, 129)[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b-',lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $cos(\omega_o t + \delta) \cdot w(t))$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega\ \to$",size=16)
grid(True)
savefig('cos_wo.jpeg')
#show()

estimateWandD(w, wo, Y, d, pow=1.75)


trueCos = cos(wo*t + d)
noise = 0.1*randn(128)   ####
n = arange(128)
y = (trueCos + noise)*wnd(n)
fmax = 1.0/(t[1]-t[0])
y = fftshift(y)
Y = fftshift(fft(y))/128.0
w = linspace(-pi*fmax, pi*fmax, 129)[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b-',lw=2)
xlim([-4,4])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $(cos(\omega_o t + \delta) + noise) \cdot w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega\ \to$",size=16)
grid(True)
savefig('cos_wo_noise.jpeg')
#show()

estimateWandD(w, wo, Y, d, pow=2.5)





def chirp(t):
	return cos(16*(1.5 + t/(2*pi))*t)


t=linspace(-pi,pi,1025)
t=t[:-1]
dt=t[1]-t[0]
fmax=1/dt
y=chirp(t)
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/1024.0
w=linspace(-pi*fmax,pi*fmax,1025)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b-',lw=2)
xlim([-75,75])
ylabel(r"$\|Y\|$",size=16)
title(r"Spectrum of $\cos\left(16\left(1.5+\dfrac{t}{2\pi}\right)t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-75,75])
ylabel(r"$\angle Y$",size=16)
xlabel(r"$\omega\ \to$",size=16)
grid(True)
savefig('chirp.jpeg')
#show()


def STFT(x, t, batchSize=64):
    t_batch = split(t, 1024//batchSize)
    x_batch = split(x, 1024//batchSize)
    X = np.zeros((1024//batchSize, batchSize), dtype=complex)
    for i in range(1024//batchSize):
        X[i] = fftshift(fft(x_batch[i]))/batchSize
    return X

t_1 = t[::64]
w = linspace(-fmax*pi,fmax*pi,65)[:-1]
T, W = meshgrid(t_1, w)

x = chirp(t)
X = STFT(x, t)
fig1 = figure(figsize = (8,7))
ax = fig1.add_subplot(211, projection='3d')
#title('The 3-D surface plot of the potential')
surf = ax.plot_surface(W, T, abs(X).T, cmap=cm.viridis)
xlabel(r"Frequency $\to$")
ylabel(r"Time $\to$")
title(r"Magnitude $\|Y\|$")
grid(True)

ax = fig1.add_subplot(212, projection='3d')
surf = ax.plot_surface(W, T, angle(X).T, cmap=cm.viridis)
xlabel(r"Frequency $\to$")
ylabel(r"Time $\to$")
title(r"Angle $\angle Y$")
grid(True)
savefig('chirp_3d.jpeg')
#show()



n = arange(1024)
x = chirp(t)*wnd(n)
X = STFT(x, t)
fig2 = figure(figsize = (8,7))
ax = fig2.add_subplot(211, projection='3d')
#title('The 3-D surface plot of the potential')
surf = ax.plot_surface(W, T, abs(X).T, cmap=cm.plasma)
xlabel(r"Frequency $\to$")
ylabel(r"Time $\to$")
title(r"Magnitude $\|Y\|$")
grid(True)

ax = fig2.add_subplot(212, projection='3d')
surf = ax.plot_surface(W, T, angle(X).T, cmap=cm.plasma)
xlabel(r"Frequency $\to$")
ylabel(r"Time $\to$")
title(r"Angle $\angle Y$")
grid(True)
savefig('chirp_3d_wnd.jpeg')
#show()
