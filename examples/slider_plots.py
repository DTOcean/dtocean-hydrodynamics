#!/usr/local/epd/canopy/User/bin/python
from numpy import transpose
from pylab import (plot, imshow, colorbar, show, subplot,
                   xlabel, ylabel, ylim, autoscale)
#from matplotlib.pyplot import plot, imshow, colorbar, show, subplot, title, figure, xlabel, ylabel

from matplotlib.widgets import Slider
from matplotlib.pyplot import axes, subplots

from dtocean_tidal.submodel.ParametricWake import (initiate_reader,
                                                   delete_reader)

##
##  Set Ct and TI values here
##

ct=0.9
ti=0.051

read_db = initiate_reader()
#read_db.nx=800
#read_db.ny=320
#read_db.alloc_arrays() 

fig, ax = subplots()

ct0 = ct
ti0 = ti

#axcolor = 'lightgoldenrodyellow'
#axti = axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
#axct  = axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

axti = axes([0.25, 0.1, 0.65, 0.03])
axct  = axes([0.25, 0.15, 0.65, 0.03])

sct = Slider(axct, 'ct', 0.1, 1.0, valinit=ct0)
sti = Slider(axti, 'ti', 0.0, 0.3, valinit=ti0)

# call fortran routine to average the nearby solutions
read_db.read_surf(ct,ti)
#print 'success'
u=read_db.u
x=read_db.x
y=read_db.y

cline = [ 1.0 - u[i,160] for i in range(len(u[:,1]))]

#figure()
#figure(facecolor="white")
fig.set_facecolor("white")
#print 'U =',u[0:4],x[0:4],y[0:4]
subplot(311)
bounds=[0,1,100]
img=imshow(transpose(u),origin='lower',cmap='jet',
           extent = [-2.0, 18.0, -4.0, 4.0] )
#imshow(u,cmap='jet')
#colorbar(img, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0, 5, 10])
#colorbar(boundaries=bounds)
#clim(0,1)
autoscale()
colorbar()

subplot(312)
xlabel('x/D')
ylabel('1 - U/Uinf')
ylim([0,1])
plt, = plot(x,cline)

#fig.tight_layout()

def update(val):
    ct = sct.val
    ti = sti.val

    print 'ct = ', ct
    print 'ti = ', ti

    # call fortran routine to average the nearby solutions
    read_db.read_surf(ct,ti)
    #print 'success'

    cline = [ 1.0 - u[i,160] for i in range(len(u[:,1]))]

    img.set_data(transpose(u))
    img.autoscale()

    plt.set_ydata(cline)


sti.on_changed(update)
sct.on_changed(update)

show()

delete_reader()
