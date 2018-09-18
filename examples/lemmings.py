# -*- coding: utf-8 -*-
import numpy as np
import pylab as plt
import radar5
#from build_c import build_c_py

r = 3.5
m = 19

def equation(x,t,lag):
    return r*x*(1.-lag(0,0)/m)

code = """void equation(double *F, const double *Y, double t, double (*lag)(int,int,void*), void* ctx, void *params) 
{
    F[0] = Y[1];
    F[1] = -lag(0,0,ctx);
}"""

init=[19.001]
res = radar5.radar5(equation,init,np.linspace(0,40,1000),lags=[.74],lagvars=[0],verbose=1)

plt.plot(res[:,0], res[:,1],'x-')
#plt.plot(res[:,0], res[:,2])

plt.show()
