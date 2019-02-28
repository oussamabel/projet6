from scipy.integrate import odeint as edo
import numpy as np

def solution(dgamma, ddgamma, mu, m, t1, n, y0):
    """
    Methode afin de calculer la solution d'une EDO en un point.
    Cette methode retourne une liste de n couples (S(t), dS(t)) avec t variant
    de 0 a t1. Elle depend des fonctions dgamma (derivee de gamma) et ddgamma
    (derivee seconde de gamma) ainsi que d'un coefficient de frottements mu,
    une masse m et d'un point de depart y0.
    """
    t = np.linspace(0, t1, n)
    return edo(f, y0, t, args=(m, mu, dgamma, ddgamma), mxstep=10000)
edo

def f(y, t, m, mu, dgamma, ddgamma):
    """
    Methode f(t) afin de resoudre l'EDO depandant des parametres y0 (point de
    depart), m (la masse), mu (le coefficient de frottement), dgamma (derivée
    de gamma) et ddgamma (derivee seconde de gamma).
    """
    s, ds = y
    dgamma_s = dgamma(s)
    ddgamma_s = ddgamma(s)
    a = prod(dgamma_s, ddgamma_s) * ds**2
    b = 9.81 * dgamma_s[1]
    c = (mu / m) * ds
    d = prod(dgamma_s, dgamma_s)
    res = ((-a - b)/d) - c
    dydt = [ds, res]
    return dydt

def prod(u, v):
    """
    Methode permetant de retourner le produit scalaire entre deux vecteurs de
    dimension 2.
    """
    return u[0] * v[0] + u[1] * v[1]
