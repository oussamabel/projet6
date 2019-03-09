#-*- coding: utf-8 -*-
from scipy.integrate import odeint as edo
import numpy as np

def solution(dgamma, ddgamma, mu, m, t1, n, c):
    """
    Méthode afin de calculer la solution approchée d'une EDO.
    Cette méthode retourne une liste de n couples (S(t), dS(t)) avec t variant
    de 0 à t1. Elle dépend des fonctions dgamma (dérivée de gamma) et ddgamma
    (derivée seconde de gamma) ainsi que d'un coefficient de frottement mu,
    d'une masse m et du paramètre c.
    """
    y0=[c,0]
    t = np.linspace(0, t1, n)
    return edo(f, y0, t, args=(m, mu, dgamma, ddgamma), mxstep=10000)

def f(y, t, m, mu, dgamma, ddgamma):
    """
    Résout l'EDO selon les parametres y (point de départ), m (la masse), mu
    (le coefficient de frottement), dgamma et ddgamma (voir plus haut).
    """
    S,dS=y
    c=np.dot(dgamma(S),ddgamma(S))
    normcarre=np.dot(dgamma(S),dgamma(S))
    k=c/normcarre
    l=(9.81*dgamma(S)[1])/normcarre
    z=(mu/m)
    dydt=[dS,-(k*dS*dS)-l-(z*dS)]
    return dydt
