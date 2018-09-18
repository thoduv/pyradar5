# -*- coding: utf-8 -*-
import numpy as np
import pylab as plt
import radar5
#from build_c import build_c_py

def equation(x,t,lag):
    return 0.2*lag(0,0)/(1.+lag(0,0)**10)-0.1*x

code = """void equation(double *F, const double *Y, double t, double (*lag)(int,int,void*), void* ctx, void *params) 
{
    F[0] = Y[1];
    F[1] = -lag(0,0,ctx);
}"""

init=[.5]
res = radar5.radar5(equation,init,np.linspace(0,100,10000),lags=[14],lagvars=[0],verbose=1)

plt.plot(res[:,0], res[:,1],'x-')
#plt.plot(res[:,0], res[:,2])

plt.show()
