#!/usr/local/epd/canopy/User/bin/python
from numpy import transpose
from pylab import (plot, imshow, colorbar, show, subplot, title, figure,
                   xlabel, ylabel, savefig, rc)

from dtocean_tidal.submodel.ParametricWake import (initiate_reader,
                                                   delete_reader)

# set the default math text to be the same as the normal
rc('mathtext', fontset='stixsans', default='regular')

##
##  Set Ct and TI values here
##

ct=0.82
ti=0.09

PlotTitle='Ct = '+str(ct)+', TI = '+str(ti)
fpref = 'Ct_'+str(ct)+'_TI_'+str(ti)

#nx=200
#ny=80
read_db = initiate_reader()
#read_db.nx=800
#read_db.ny=320
#read_db.alloc_arrays() 

# call fortran routine to average the nearby solutions
read_db.read_surf(ct,ti)

print 'success'
u=read_db.u
x=read_db.x
y=read_db.y

#cline = [ 1.0 - u[i,160] for i in range(len(u[:,1]))]
cline = [ u[i,160] for i in range(len(u[:,1]))]

f1 = open(fpref+'_def.dat', 'w')
f2 = open(fpref+'_rel.dat', 'w')
for i in range(len(cline)):
  string=str(x[i]) + '   ' + str(cline[i]) + '\n'
  f1.write(string)
  string=str(x[i]) + '   ' + str(1.0-cline[i]) + '\n'
  f2.write(string)

figure(facecolor="white")
#print 'U =',u[0:4],x[0:4],y[0:4]
subplot(211)
title(PlotTitle)
imshow(transpose(u),cmap='jet',extent = [-2.0, 18.0, -4.0, 4.0] )
#imshow(u,cmap='jet')
colorbar()

subplot(212)
xlabel('x/D')
#ylabel('1 - U/Uinf')
ylabel(r'$U_{rel}$')

plot(x,cline)

#savefig('Interp.png')
savefig(fpref+'.png')
show()

delete_reader()
