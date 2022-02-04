import sys
import mpl_toolkits.mplot3d.axes3d as p3
from pylab import *
from numpy import *
from numpy.linalg import *

try:
	if __name__ == '__main__':
		Nx = 25
		Ny = 25
		radius = 8
		Niter = 1500
	if(len(sys.argv)) == 2:
		Nx = int(sys.argv[1])
	elif(len(sys.argv)) == 3:
		Nx = int(sys.argv[1])
		Ny = int(sys.argv[2])
	elif(len(sys.argv)) == 4:
		Nx = int(sys.argv[1])
		Ny = int(sys.argv[2])
		radius = int(sys.argv[3])
	elif(len(sys.argv)) == 5:
		Nx = int(sys.argv[1])
		Ny = int(sys.argv[2])
		radius = int(sys.argv[3])
		Niter = int(sys.argv[4])

except IOError:
	print('Invalid file\n')

phi = zeros((Ny, Nx))
x = linspace(-(Nx-1)/2, (Nx-1)/2, num=Nx)
y = linspace(-(Ny-1)/2, (Ny-1)/2, num=Ny)
(X,Y) = meshgrid(x, y)
ii = where(X*X + Y*Y <= (radius)**2)
phi[ii] = 1.0
figure(1)
plot(ii[0] - (Nx-1)/2, ii[1] - (Ny-1)/2, 'ro')
contour(x,y,phi)
grid(True)
savefig('contour1.jpeg', bbox_inches = 'tight')

errors = zeros(Niter)
for k in range(Niter):
	oldphi = phi.copy()
	phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2] + phi[1:-1,2:] + phi[0:-2,1:-1] + phi[2:,1:-1])
	phi[:,0] = phi[:,1]
	phi[:,-1] = phi[:,-2]
	phi[-1,:] = phi[-2,:]
	phi[ii] = 1.0
	errors[k] = (abs(phi - oldphi)).max()


N = linspace(0, Niter-1, num=Niter)
N1 = array([50*k for k in range(0, int(Niter/50))])
N2 = array([50*k for k in range(10, int(Niter/50))])
errors_1 = array([errors[50*k] for k in range(0, int(Niter/50))])
figure(2)
semilogy(N, errors)
xlabel(r'$No.\ of\ Iterations$', size=20)
ylabel(r'$Error\ Values$', size=20)
grid(True)
savefig('error_semilogy.jpeg', bbox_inches = 'tight')

figure(3)
loglog(N, errors)
xlabel(r'$No.\ of\ Iterations$', size=20)
ylabel(r'$Error\ Values$', size=20)
grid(True)
savefig('error_loglog.jpeg', bbox_inches = 'tight')


M1 = zeros([Niter,2])
M1[:,0] = 1
M1[:,1] = N
M2 = zeros([Niter-500,2])
M2[:,0] = 1
M2[:,1] = N[500:]
o1 = lstsq(M1, log(errors), rcond=-1)[0]
o2 = lstsq(M2, log(errors[500:]), rcond=-1)[0]
error_1 = exp(matmul(M1, o1))
error_2 = exp(matmul(M2, o2))
errors_n1 = array([error_1[50*k] for k in range(0, int(Niter/50))])
errors_n2 = array([error_2[50*k] for k in range(0, int((Niter-500)/50))])

figure(4)
semilogy(N1, errors_1, 'ro', label=r'$errors$')
semilogy(N1, errors_n1, 'g', label=r'$fit1$')
semilogy(N2, errors_n2, 'b', label=r'$fit2$')
xlabel(r'$No.\ of\ Iterations$', size=20)
ylabel(r'$Error\ Values$', size=20)
legend(loc = 'upper right')
grid(True)
savefig('error_n_semilogy.jpeg', bbox_inches = 'tight')

figure(5)
loglog(N1, errors_1, 'ro', label=r'$errors$')
loglog(N1, errors_n1, 'g', label=r'$fit1$')
loglog(N2, errors_n2, 'b', label=r'$fit2$')
xlabel(r'$No.\ of\ Iterations$', size=20)
ylabel(r'$Error\ Values$', size=20)
legend(loc = 'upper right')
grid(True)
savefig('error_n_loglog.jpeg', bbox_inches = 'tight')

figure(6)
plot(ii[0] - (Nx-1)/2, ii[1] - (Ny-1)/2, 'ro')
contour(x,y,phi)
grid(True)
savefig('contour2.jpeg', bbox_inches = 'tight')


fig1=figure(7)
ax = p3.Axes3D(fig1)
title('The 3-D surface plot of the potential')
surf = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1, cmap=cm.jet)
savefig('surface.jpeg', bbox_inches = 'tight')

figure(8)
Jx = zeros((Ny, Nx))
Jy = zeros((Ny, Nx))
Jx[:,1:-1] = 0.5*(phi[:,0:-2] - phi[:,2:])
Jy[1:-1,:] = 0.5*(phi[0:-2,:] - phi[2:,:])
Jy[-1,:] = 0
Jx[:,0] = 0
Jx[:,-1] = 0
Jy[0,:] = 0
quiver(X, Y, Jx, Jy, scale = 6)
plot(ii[0] - (Nx-1)/2, ii[1] - (Ny-1)/2, 'ro')
grid(True)
savefig('vector.jpeg', bbox_inches = 'tight')




