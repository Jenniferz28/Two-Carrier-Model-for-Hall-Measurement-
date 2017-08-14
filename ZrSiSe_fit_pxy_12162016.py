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

   
def f_pxy (B, n_e, mu_e, n_h, mu_h):
    return (10.**(-5.)/1.60217662)* B*((n_h*(mu_h**2)-n_e*(mu_e**2)) + (n_h-n_e)*(mu_h**2)*(mu_e**2)*(B**2))/((n_h*mu_h+n_e*mu_e)**2 + ((n_h-n_e)**2) *(mu_e**2)*(mu_h**2))


os.chdir("Z:\\CMS\\Qiong Zhou\\For Yu-che\\10262016\\fit_pxy_12162016")


#2K data
data=pd.read_csv("MC006_S7_cutS7-2K_(Rhoxx,Rhoxy).dat",sep="\t", header=0)
B=np.array(data["H(Tesla)"])
#pxx=np.array(data["S7-2K_Rhoxx(Ohm-cm)"])
pxy=np.array(data["S7-2K_Rhoxy(Ohm-cm)"])
    #define fitting and plot
def fit_function(params, x=B, dat2=pxy):
    n_e = params['n_e']
    mu_e=params['mu_e']
    n_h=params['n_h']
    mu_h=params['mu_h']
    model2 = (10.**(-5.)/1.60217662)* B*((n_h*(mu_h**2)-n_e*(mu_e**2)) + (n_h-n_e)*(mu_h**2)*(mu_e**2)*(B**2))/((n_h*mu_h+n_e*mu_e)**2 + ((n_h-n_e)**2) *(mu_e**2)*(mu_h**2))
    resid2 = dat2 - model2  
    return np.array(resid2)
    
params = Parameters()
params.add('n_e', value=0.1)
params.add('mu_e', value=0.01)
params.add('n_h', value=0.1)
params.add('mu_h', value=0.01)
fit = minimize(fit_function, params, args=(B,pxy))
lmfit.printfuncs.report_fit(fit.params, min_correl=0.5)     
   
def myplot(B,pxy,params):  
    plt.ylabel('p_xy', fontsize = 16)
    plt.xlabel('B', fontsize = 16)
    plt.xlim(0,10)
    plt.plot(B, f_pxy(B, fit.params['n_e'], fit.params['mu_e'], fit.params['n_h'], fit.params['mu_h']),'r',\
B, pxy,'b',linewidth=2.0)
    plt.legend(["fitting data","real data"],loc=2)
    plt.plot()
    
#plot 2K fits and curves
f2_K=myplot(B,pxy, params)    

    
#loop through all tempeartures
for i in [2,3,5,10,15,20,25,30,40,50]:
    data=pd.read_csv("MC006_S7_cutS7-%sK_(Rhoxx,Rhoxy).dat" %i,sep="\t", header=0)
    B=np.array(data["H(Tesla)"])
    pxx=np.array(data["S7-%sK_Rhoxx(Ohm-cm)" %i])
    pxy=np.array(data["S7-%sK_Rhoxy(Ohm-cm)" %i])
       
    params = Parameters()
    params.add('n_e', value=0.1)
    params.add('mu_e', value=0.1)
    params.add('n_h', value=0.1)
    params.add('mu_h', value=0.1)    
    # fit on data and print output
    fit = minimize(fit_function, params, args=(B,pxy))
    print str(i)+"K",lmfit.printfuncs.report_fit(fit.params, min_correl=0.5)
    myplot(B, pxy, params)
#    #save fitted data to txt file
    d1=np.column_stack((B,pxy,f_pxy(B, fit.params['n_e'], fit.params['mu_e'], fit.params['n_h'], fit.params['mu_h'])))
    df=pd.DataFrame(d1,columns=['B','pxy','f_pxy'])
    df.to_csv('fitted_pxy_{}.txt'.format(i),sep='\t',header=True)
    




for i in [75,100,125,150,175,200,225,250]:
    data=pd.read_csv("MC006_S7_cutS7-%sK_(Rhoxx,Rhoxy).dat" %i,sep="\t", header=0)
    B=np.array(data["H(Tesla)"])
    pxx=np.array(data["S7-%sK_Rhoxx(Ohm-cm)" %i])
    pxy=np.array(data["S7-%sK_Rhoxy(Ohm-cm)" %i])
       
    params = Parameters()
    params.add('n_e', value=0.1)
    params.add('mu_e', value=0.1)
    params.add('n_h', value=0.1)
    params.add('mu_h', value=0.1)    
    # fit on data and print output
    fit = minimize(fit_function, params, args=(B,pxy))
    print str(i)+"K",lmfit.printfuncs.report_fit(fit.params, min_correl=0.5)
    myplot(B, pxy, params)
#    #save fitted data to txt file
    d1=np.column_stack((B,pxy,f_pxy(B, fit.params['n_e'], fit.params['mu_e'], fit.params['n_h'], fit.params['mu_h'])))
    df=pd.DataFrame(d1,columns=['B','pxy','f_pxy'])
    df.to_csv('fitted_pxy_{}.txt'.format(i),sep='\t',header=True)
    









