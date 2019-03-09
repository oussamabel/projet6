from scipy.integrate import odeint as edo
import numpy as np

def solution(dgamma, ddgamma, mu, m, t1, n, y0):
    """
    Methode afin de calculer la solution d'une EDO en un point.
    Cette methode retourne une liste de n couples (S(t), dS(t)) avec t variant
    de 0 a t1. Elle depend des fonctions dgamma (derivee de gamma) et ddgamma (derivee
    seconde de gamma) ainsi que d'un coefficient de frottements mu, une masse m
    et d'un point de depart y0.
    """
    t = np.linspace(0, float(t1),float(n))
    return edo(f, y0, t, args=(m, mu, dgamma, ddgamma), mxstep=10000)

def f(y, t, m, mu, dgamma, ddgamma):
    """
    Methode fct(t) afin de resoudre l'EDO depandant des parametres y0 (point de
    depart), m (la masse), mu (le coefficient de frottement), df (derivée de
    f) et ddf (derivee seconde de f).
    """
    s, ds = y
    dgamma_s = dgamma(s)
    ddgamma_s = ddgamma(s)
    a = np.dot(dgamma_s, ddgamma_s) * ds**2
    b = 9.81 * dgamma_s[1]
    c = (mu / m) * ds
    d = np.dot(dgamma_s, dgamma_s)
    res = ((-a - b)/d) - c
    dydt = [ds, res]
    return dydt
