from scipy.integrate import odeint as edo
import numpy as np

def solution(df, ddf, mu, m, t1, n, y0):
    """
    Methode afin de calculer la solution d'une EDO en un point.
    Cette methode retourne une liste de n couples (S(t), dS(t)) avec t variant
    de 0 a t1. Elle depend des fonctions dxi (derivee de xi) et ddxi (derivee
    seconde de xi) ainsi que d'un coefficient de frottements mu, une masse m
    et d'un point de depart y0.
    """
    t = np.linspace(0, float(t1),float(n))
    return edo(fct, y0, t, args=(m, mu, df, ddf), mxstep=10000)

def fct(y, t, m, mu, df, ddf):
    """
    Methode fct(t) afin de resoudre l'EDO depandant des parametres y0 (point de
    depart), m (la masse), mu (le coefficient de frottement), df (derivée de
    f) et ddf (derivee seconde de f).
    """
    s, ds = y
    df_s = df(s)
    ddf_s = ddf(s)
    a = prod(df_s, ddf_s) * ds**2
    b = 9.81 * df_s[1]
    c = (mu / m) * ds
    d = prod(df_s, df_s)
    res = ((-a - b)/d) - c
    dydt = [ds, res]
    return dydt

def prod(u, v):
    """
    Methode permetant de retourner le produit scalaire entre deux vecteurs de
    dimension 2.
    """
    return u[0] * v[0] + u[1] * v[1]
