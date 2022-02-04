import sys
from numpy import *
from pylab import *


if len(sys.argv)==6:
    n=sys.argv[0] 
    M=sys.argv[1]    
    nk=sys.argv[2] 
    u0=sys.argv[3]  
    p=sys.argv[4] 
    Msig=sys.argv[5] 
    params=[n,M,nk,u0,p,Msig]
else:
    params=[[100,5,500,5,0.25,2],[100,5,500,5,0.50,2],
        [100,5,500,5,1,2],[100,5,500,10,0.25,2],
        [100,5,500,10,0.50,2],[100,5,500,10,1,2]]



def loop(params, q):

	n = params[0]
	M = params[1]
	nk = params[2]
	u0 = params[3]
	p = params[4]
	Msig = params[5]

	xx = zeros((n*M))
	u = zeros((n*M))
	dx = zeros((n*M))

	I = []
	X = []
	V = []

	for i in range(1,nk):
		ii = where(xx>0)
		dx[ii] = u[ii] + 0.5
		xx[ii] = xx[ii] + dx[ii]
		u[ii] = u[ii] + 1
		j = where(xx>n)
		dx[ii[0][j]] = 0
		xx[ii[0][j]] = 0
		u[ii[0][j]] = 0
		kk = where(u>u0)
		ll = where(rand(len(kk[0]))<=p)
		kl = kk[0][ll]
		u[kl] = 0
		if q == False:
			rho = rand(len(kl))
		if q == True:
			rho = power(0.5,len(kl))
		xx[kl] = xx[kl] - dx[kl]*rho
		I.extend(xx[kl].tolist())
		m = int(rand()*Msig + M)
		n0 = where(xx==0)
		nv = min((n*M - len(n0), m))
		xx[n0[:nv]] = 1
		u[n0[:nv]] = 0
		dx[n0[:nv]] = 0
		X.extend(xx[ii].tolist())
		V.extend(u[ii].tolist())
	return X,V,I

def plot_no_of_elec(X,u0,p):
    figure(1)
    hist(X,bins = arange(0,101,0.5),rwidth = 0.8,color = 'g')
    title('Number of Electrons vs $x$ with $u_0=$%f and p=%f'%(u0,p))
    xlabel('$x$')
    ylabel('Number of electrons')
    show()

def plot_intensity_map(I,u0,p):
    histogram_= hist(I,bins = arange(0,101,1),rwidth = 0.8,color = 'r')
    x = histogram_[1][1:]
    y = histogram_[0]
    fig, (ax) = subplots(nrows = 1, sharex = True)
    extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
    intensity = ax.imshow(y[np.newaxis,:], cmap = "gray", aspect = "auto", extent = extent)
    ax.set_yticks([])
    ax.set_xlim(extent[0], extent[1])
    title('Intensity map with $u_0=$%f and p=%f'%(u0,p))
    xlabel('$x$')
    colorbar(intensity)
    tight_layout()
    show()
    
def plot_intensity(X,V,u0,p):
    figure(0)
    histogram = hist(I,bins = arange(0,101,0.5),rwidth = 0.8,color = 'r')
    title('Intensity histogram with $u_0=$%f and p=%f'%(u0,p))
    xlabel('$x$')
    ylabel('Intensity')
    show()
    return histogram

def plot_phase(X,V,u0,p):
    figure(2)
    plot(X,V,'bo')
    title('Electron Phase Space with $u_0=$%f and p=%f'%(u0,p))
    xlabel('$x$')
    ylabel('Velocity-$v$')
    show()

    
param=[100,5,500,5,0.25,2]
X,V,I=loop(param,False)
histo=plot_intensity(X,V,param[3],param[4])
plot_no_of_elec(X,param[3],param[4])
plot_phase(X,V,param[3],param[4])
plot_intensity_map(I,param[3],param[4])

param=[100,5,500,5,0.25,2]
X,V,I=loop(param,True)
plot_intensity(X,V,param[3],param[4])
plot_no_of_elec(X,param[3],param[4])
plot_phase(X,V,param[3],param[4])
plot_intensity_map(I,param[3],param[4])

for param in params:
    X,V,I=loop(param,False)
    plot_intensity(X,V,param[3],param[4])
    plot_no_of_elec(X,param[3],param[4])
    plot_phase(X,V,param[3],param[4])
    plot_intensity_map(I,param[3],param[4])



