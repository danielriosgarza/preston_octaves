# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 21:28:53 2016

@author: user
"""
import random as rd
import math as m
import numpy as np
from pylab import *

def polar_normal_variable(mu, sigma):
    #two_pi = 2.0*3.14159265358979323846;
    s = 2
    while s>1:
        a,b = rd.uniform(0,1), rd.uniform(0,1)
        v1 = (2*a)-1
        v2 = (2*b)-1
        s = (v1**2) + (v2**2)
        if s<1:
            x = m.sqrt((-2*m.log(s))/s)*v1
    return mu + m.sqrt(sigma)*x
    
def interval_counter(x,interval):
    count =0
    for i in x:
        if (i<interval[1] and i>interval[0]):
            count+=1
        elif i==interval[0] or i==interval[1]:
            count+=0.5
    return count

def interval_counter2(x,interval):
    count =0
    for i in x:
        if (i<interval[1] and i>interval[0]):
            count+=i
        elif i==interval[0] or i==interval[1]:
            count+=0.5*i
    return count
    

def preston_octaves(species_counts):
    ma = int(m.log((max(species_counts)+1), 2))
    species_counts = np.log2(np.array(species_counts))
    x,y,z=[],[],[]
    for i in xrange(ma):
        b = interval_counter(species_counts, [i,i+1])
        c = interval_counter2(species_counts, [i,i+1])
        x.append(i)
        y.append(b)
        z.append(c)
    return np.array(x),np.array(y),np.array(z)

def fitted_curve(y_0, l, a):
    y =[]
    ls = np.linspace(0,2*l, 20000)
    for i in ls:
        r = abs(l-i)
        y.append(y_0*m.exp(-((a*r)**2)))
    ls = ls-l
    return ls,y


def p_curve():
   ls = np.linspace(-5, 5, 10000)
   ls2 = ls**2
   ls2 = ls2/(-2.)
   return ls, 1./sqrt(2*pi)*np.exp(ls2)        

def get_x(p_curve, sig_l):
    x = 0
    y = 0
    while y<sig_l:
        x+=0.001
        ls = np.linspace(-x, x, 1000)
        ls2 = ls**2
        ls2 = ls2/(-2.)
        y = np.trapz(1./sqrt(2*pi)*np.exp(ls2), ls)
    return x


def canonical_curve(n):
    u = (n-1)/float(n)
    x = get_x(p_curve, u)
    rmax = 1.442695*(x**2)
    a = 0.490129/x
    N = (0.276525716866*n)/x
    xp,yp = fitted_curve(N, rmax, a)
    sigma = 1.42695*x
    r = 2**(rmax)
    pred_commonest = rmax*r
    pred_rarest = rmax/float(r)
    r_2 = pred_commonest/pred_rarest
    predicted_n_ind = .5*sqrt(2*pi)*pred_rarest*r_2*sigma
    I_over_m = .5*sqrt(2*pi)*r_2*sigma
    pseudo_sample = np.array([2**(polar_normal_variable(rmax, sigma**2)) for i in xrange(int(n))])
    preston_x, preston_y, preston_z = preston_octaves(pseudo_sample)
    return {'species_x':xp, 'species_y': yp, 'x':x, 'a':a, 'rmax':rmax, 'y0':N, 'x':x, 'sigma':sigma, \
    'preston_x':preston_x, 'preston_y':preston_y, 'preston_z':preston_z, 'pred_common':pred_commonest, \
    'pred_rarest' : pred_rarest, 'r':r, 'r2' : r_2, 'pred_individuals': predicted_n_ind, 'I/m':I_over_m}

