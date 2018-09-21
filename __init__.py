# -*- coding: utf-8 -*-
try:
	import radar5.__config__
except:
	import __config__
from .radar5 import radar5
radar5 = radar5

def test():
    import numpy as np
    import pylab as plt

    code = """void equation(double *F, const double *Y, double t, double (*lag)(int,int,void*), void* ctx, void *params) 
    {
        F[0] = Y[1];
        F[1] = -lag(0,0,ctx);
    }"""

    res = radar5(code,[0,1],np.linspace(0,100,100),lags=[.1],lagvars=[0],verbose=1,full_output=True)

    plt.plot(res[:,0], res[:,1])
    plt.plot(res[:,0], res[:,2])

    plt.show()
    
