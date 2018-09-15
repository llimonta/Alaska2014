import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import trapz

def QECurveFit(x,a,b,c,d,e,f,g,h,i,l,m,n):
    return n+a*x+b*x**2+c*x**3+d*x**4+e*x**5+f*x**6+g*x**7+h*x**8+i*x**9+l*x**10+m*x**11

def TrapzIntegral(x,y):
    delta = np.subtract(x[1:-1],x[0:-2])
    Il=np.multiply(delta,y[0:-2])
    Ir=np.multiply(delta,y[1:-1])
    return sum((Il+Ir))/2

# Read Data
qeEXdata = pd.read_csv('/home/limo/Documents/Thesis/CodeForPlots/CameraQEExExf.csv')
qeEXx_points = qeEXdata['x']*10
qeEXy_points = qeEXdata['Curve1']
qeBVdata = pd.read_csv('/home/limo/Documents/Thesis/CodeForPlots/CameraQEBvBvf.csv')
qeBVx_points = qeBVdata['x']*10
qeBVy_points = qeBVdata['Curve2']
MeteorSpectrographyData = pd.read_csv('/home/limo/Documents/Thesis/CodeForPlots/Vojacek2015MeteorSx336.csv')
MeteorX_points = MeteorSpectrographyData['x']
MeteorY_points = MeteorSpectrographyData['Curve1']
# Plot Camera Quantum Efficiency
fig = plt.figure(figsize=(8, 6))
plt.plot(qeEXx_points,qeEXy_points,label='Ex/Exf',linewidth=3)
plt.plot(qeBVx_points,qeBVy_points,label='Bv/Bvf',linewidth=3)
plt.xlabel('Wavelength [$\AA$]',fontsize=20)
plt.ylabel('Quantum Efficiency',fontsize=20)
plt.title('Quantum Efficiency Andor iXon Ultra 897',fontsize=24)
plt.legend(loc='best',fontsize=15)

# Plot SuperImposition Of Spectrum and Camera QE
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(qeBVx_points,qeBVy_points,color='orange',label='Bv/Bvf',linewidth=3)
ax2.plot(MeteorX_points,MeteorY_points,color = 'g',label='Sx336 Spectrum',linewidth = 3)
ax1.set_xlabel('Wavelength [$\AA$]',fontsize=24)
ax1.set_ylabel('Quantum Efficiency Andor Camera',fontsize=24)
ax2.set_ylabel('Device Raw Units', color='g',fontsize=24)
plt.title('Meteor spectrum and Camera QE',fontsize=24)
plt.legend(loc='best')

# Plot expected read count given previous spectrum

# Fit QE in BVF
plt.figure()
params1,cov1=curve_fit(QECurveFit,qeBVx_points,qeBVy_points,)
# print params1, cov1
# print len(params1)
r2_1=(sum((QECurveFit(qeBVx_points,*params1)-np.mean(qeBVy_points))**2)/sum((qeBVy_points-np.mean(qeBVy_points))**2))
print 'r2_1',r2_1
print 1-(1-r2_1)*(len(qeBVy_points))/(len(qeBVy_points)-len(params1)-1)
plt.plot(qeBVx_points,QECurveFit(qeBVx_points,*params1),'g')
plt.plot(qeBVx_points,qeBVy_points,color='orange')

# Convert in Raw count
fig = plt.figure(figsize=(8, 6))
plt.plot(MeteorX_points,MeteorY_points,color = 'g',label='Sx336 Original Spectrum',linewidth = 3)
plt.plot(MeteorX_points,MeteorY_points*QECurveFit(MeteorX_points,*params1)/100,'r',label='Sx336 Device Spectrum',linewidth = 3)
Meteor_Camera_Response=MeteorY_points*QECurveFit(MeteorX_points,*params1)/100
print QECurveFit(MeteorX_points,*params1)/100
plt.xlabel('Wavelength [$\AA$]',fontsize=20)
plt.ylabel('Device Raw Units',fontsize=20)
plt.title('Sx336 Meteor Spectum',fontsize=20)
plt.legend(loc='best',fontsize=15)

# Compute difference between real spectrum and detected spectrum
# print MeteorY_points
# print MeteorX_points
plt.figure(),plt.plot(MeteorX_points)
MeteorTotal = trapz(np.array(MeteorY_points),np.array(MeteorX_points))

print 'MeteorTotal',MeteorTotal
CameraMeteor = trapz(Meteor_Camera_Response,MeteorX_points)

print 'CameraMeteor',CameraMeteor
print 'Loss',CameraMeteor/MeteorTotal










plt.show()
