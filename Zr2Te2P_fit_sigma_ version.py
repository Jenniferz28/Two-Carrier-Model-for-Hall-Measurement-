# -*- coding: utf-8 -*-
"""
@author: qzhou
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import lmfit

   
def f_sigma_xx (B, n_e, mu_e, n_h, mu_h):
    return 1.602*10**5*n_e*mu_e/(1.+(mu_e*B)**2)+1.602*10**5*n_h*mu_h/(1.+(mu_h*B)**2)
    
def f_sigma_xy (B, n_e, mu_e, n_h, mu_h):
    return 1.602*10**5*n_e*mu_e**2*B/(1+(mu_e*B)**2)-1.602*10**5*n_h*mu_h**2*B/(1+(mu_h*B)**2)

    
os.chdir("Z:\\CMS\\Qiong Zhou\\For Kuan-wen\\Zr2Te2P")

#2K data
data=pd.read_csv("2K_Zr2Te2P.dat",sep="\t", header=0)
B=np.array(data.iloc[:,12])
p_xx=np.array(data['p_xx Ohm*cm'])
p_xy=np.array(data['p_xy Ohm*cm'])
sigma_xx = p_xx/(p_xx**2 + p_xy**2)
sigma_xy =-p_xy/(p_xx**2 + p_xy**2)
#rho = p_xx +1j* p_xy
#sigma=rho**(-1)
#sigma_xx = sigma.real
#sigma_xy = sigma.imag

def fit_function(params, x=B, dat1=sigma_xx, dat2=sigma_xy):
    n_e = params['n_e']
    mu_e=params['mu_e']
    n_h=params['n_h']
    mu_h=params['mu_h']
    model1 = 1.602*10**5*n_e*mu_e/(1.+(mu_e*B)**2)+1.602*10**5*n_h*mu_h/(1.+(mu_h*B)**2)
    model2 = 1.602*10**5*n_e*mu_e**2*B/(1+(mu_e*B)**2)-1.602*10**5*n_h*mu_h**2*B/(1+(mu_h*B)**2)
    resid1 = dat1 - model1
    resid2 = (dat2 - model2)*10
    return np.concatenate((resid1, resid2))
#n_e=params[0]
#mu_e=params[1]
#n_h=params[2]
#mu_h=params[3]
# model1 "sigma_xx"; model2"sigma_xy"

p0 = [1,1,1,1] #initial fitting parameters
params = Parameters()
params.add('n_e', value=1)
params.add('mu_e', value=0.1)
params.add('n_h', value=1)
params.add('mu_h', value=0.1)
fit = minimize(fit_function, params, args=(B,sigma_xx,sigma_xy))
lmfit.printfuncs.report_fit(fit.params, min_correl=0.5)

#define plot function
def myplot(B,sigma_xx,sigma_xy,params):  
    ax1=plt.subplot(2,1,1)
    plt.ylabel('Sigma_xy', fontsize = 16)
    plt.xlabel('B', fontsize = 16)
    plt.xlim(0,10)
    plt.plot(B, f_sigma_xy(B, fit.params['n_e'], fit.params['mu_e'], fit.params['n_h'], fit.params['mu_h']),'r',\
B, sigma_xy,'b',linewidth=2.0)
    plt.legend(["fitting data","real data"])

    ax2=plt.subplot(2,1,2)
    plt.ylabel('Sigma_xx', fontsize = 16)
    plt.xlabel('B', fontsize = 16)
    plt.xlim(0,10)
    plt.plot(B, f_sigma_xx(B,fit.params['n_e'],fit.params['mu_e'], fit.params['n_h'], fit.params['mu_h']),'r',
B, sigma_xx,'b',linewidth=2.0)
    plt.legend(["fitting data","real data"])
    plt.plot()
#plot 2K fits and curves
f2_K=myplot(B,sigma_xx,sigma_xy,params)



#loop through all tempeartures and plot the curves
for i in [2,10,30,50,70,100,150,200,250,300]:
    data=pd.read_csv("%sK_Zr2Te2P.dat" %i,sep="\t", header=0)
    B=np.array(data.iloc[:,12])
    p_xx=np.array(data['p_xx Ohm*cm'])
    p_xy=np.array(data['p_xy Ohm*cm'])
    sigma_xx = p_xx/(p_xx**2 + p_xy**2)
    sigma_xy =-p_xy/(p_xx**2 + p_xy**2)
    # fit on data and print output
    fit = minimize(fit_function, params, args=(B,sigma_xx,sigma_xy))
    print str(i)+"K",lmfit.printfuncs.report_fit(fit.params, min_correl=0.5)
    #plot the fits and curves
    myplot(B,sigma_xx,sigma_xy,params)
    #save fitted data to txt file
    d1=np.column_stack((B,sigma_xx,sigma_xy, f_sigma_xx(B, fit.params['n_e'], fit.params['mu_e'], fit.params['n_h'], fit.params['mu_h']), f_sigma_xy(B, fit.params['n_e'], fit.params['mu_e'], fit.params['n_h'], fit.params['mu_h'])))
    df=pd.DataFrame(d1,columns=['B','sigma_xx','sigma_xy','f_sigma_xx','f_sigma_xy'])
    df.to_csv('fitted_sigma_{}.txt'.format(i),sep='\t',header=True)
    
    





