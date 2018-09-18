# -*- coding: utf-8 -*-
import numpy as np
import pylab as plt
import radar5
from build_c import build_c_py

code = """//#include <math.h>
void equation(double *F, const double *Y, double t, double (*lag)(int,int,void*), void* ctx, void *params) 
{
	double *rpar = params; 
	double Er=Y[0], Ei=Y[1], n=Y[2];
	double delay_Er=lag(0,0,ctx), delay_Ei=lag(0,1,ctx);
	double T = rpar[0], K=rpar[1], p=rpar[2], alpha=rpar[3], Phi=rpar[4];
	F[0] = (Er-alpha*Ei)*n/2 + K*cos(Phi)*delay_Er - K*sin(Phi)*delay_Ei;
	F[1] = (Ei+alpha*Er)*n/2 + K*cos(Phi)*delay_Ei + K*sin(Phi)*delay_Er;
	F[2] = (p - n - (1+n)*(Er*Er + Ei*Ei))/T;
}"""

def py_equation(Y,t,lag,params):
	Er, Ei, n = Y
	T,K,p,alpha,Phi=params
	delay_Er, delay_Ei = lag(0,0), lag(0,1)
	return [ (Er-alpha*Ei)*n/2 + K*np.cos(Phi)*delay_Er - K*np.sin(Phi)*delay_Ei,
	  (Ei+alpha*Er)*n/2 + K*np.cos(Phi)*delay_Ei + K*np.sin(Phi)*delay_Er,
	  (p - n - (1+n)*(Er*Er + Ei*Ei))/T ]
	

c_equation = build_c_py(code, 'equation', cflags = ['-g'])


#t0 = np.linspace(0,100,1000)
#y0 = np.zeros((len(t0),3))
#y0[:,0] = np.sin(t0/20)

#print y0.shape

res = radar5.radar5(code, [1,0,0], 1000, lags=[100], lagvars=[0,1], rpar=[200, 0.1, 1, 1, 0], ngrid=10000,verbose=1)

t0 = res[:,0]
y0 = res[:,1:]

res = radar5.radar5(c_equation, res, 1000, lags=[100], lagvars=[0,1], rpar=[200, 0.1, 1, 1, 0], ngrid=10000,verbose=1)

res2 = radar5.radar5(c_equation, [1,0,0], 2000, lags=[100], lagvars=[0,1], rpar=[200, 0.1, 1, 1, 0], ngrid=10000,verbose=1)

#plt.plot(res[:,0], res[:,1]**2 + res[:,2]**2,'b')
#plt.plot(res2[:,0]-1000, res2[:,1]**2 + res2[:,2]**2,'b:')
#plt.plot(t0-t0[-1], y0[:,0]**2 + y0[:,1]**2, 'b')
#plt.twinx()
#plt.plot(res[:,0], np.angle(res[:,1] + 1j*res[:,2]), 'g')
#plt.plot(res2[:,0]-1000, np.angle(res2[:,1] + 1j*res2[:,2]), 'g:')
#plt.plot(t0-t0[-1], np.angle(y0[:,0] + 1j*y0[:,1]), 'g')

plt.plot(res[:,0], res[:,1],'b')
plt.plot(res2[:,0]-1000, res2[:,1],'b:')
plt.plot(t0-t0[-1], y0[:,0], 'b')

plt.plot(res[:,0], res[:,2],'g')
plt.plot(res2[:,0]-1000, res2[:,2],'g:')
plt.plot(t0-t0[-1], y0[:,1], 'g')

plt.show()
