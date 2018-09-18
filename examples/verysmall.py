# -*- coding: utf-8 -*-
import numpy as np
import pylab as plt
import radar5

def equation(x,t,lag):
    return [x[1],-lag(0,0)]

code = """void equation(double *F, const double *Y, double t, double (*lag)(int,int,void*), void* ctx, void *params) 
{
    F[0] = Y[1];
    F[1] = -lag(0,0,ctx);
}"""

init=[0,1]
res = radar5.radar5(code,init,np.linspace(0,100,100),lags=[.1],lagvars=[0],verbose=1,full_output=True)

plt.plot(res[:,0], res[:,1])
plt.plot(res[:,0], res[:,2])

plt.show()

