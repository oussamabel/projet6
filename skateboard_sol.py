from math import sin, cos, pi
import sys, skateboard
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

def xi(s):
    """
    Fonction xi(s).
    """
    return (s, (1 / 2) * (s + 2)**2)

def dxi(s):
    """
    Derivee de la fonction xi(s).
    """
    return (1, s+2)

def ddxi(s):
    """
    Derivee seconde de la fonction xi(s).
    """
    return (0, 1)

def eta(s):
    return (s, (1 + cos(s))/(1 + 0.1*s))

def deta(s, num = False):
    tmp = -sin(s) * (1 + 0.1*s)- 0.1 - 0.1*cos(s)
    if (num):
        return tmp
    else:
        return (1, tmp / (1 + 0.1*s)**2)

def ddeta(s):
    sin_s = sin(s)
    cos_s = cos(s)
    return (0, ((-cos_s - (cos_s * 0.1 * s + sin_s * 0.1) + 0.1 * sin_s)
                *(1 + 0.1*s) - deta(s, True) * 0.2) / (1 + 0.1 * s)**3)

def s(m, t, mu, y0):
    """
    Fonction dont on va rechercher les racines (en fonction du temps) dependant
    de plusieurs parametres : m (la masse), mu (le coefficient de frottements)
    et y0 (le point de depart).
    La racine correspondra au moment ou le mobile passe par le point p1.
    """
    res = skateboard.solution(lambda x: dxi(x), lambda x: ddxi(x),
                               mu, m, t, 2, y0)
    return res[1][0] #retourne s(t), donc la premiere composante du deuxieme
    # couple du tableau

def ds(m, t, mu, y0):
    """
    Fonction dont on va rechercher les racines (en fonction du temps) dependant
    de plusieurs parametres : m (la masse), mu (le coefficient de frottements)
    et y0 (le point de depart).
    La racine correspondra au moment ou le mobile passe par le point p1.
    """
    res = skateboard.solution(lambda x: dxi(x), lambda x: ddxi(x),
                               mu, m, t, 2, y0)
    return res[1][1] #retourne ds(t), donc la deuxieme composante du deuxieme
    # couple du tableau

def afficher_s(tab):
    res = ""
    for i in range(0, 1000):
        res += str(i) + " : " + str(tab[i][0])
        res += "\n"
    print(res)
def solve(mu,m,c,is_list):
    if is_list:
        return lambda l:skateboard.solution(dxi,ddxi,mu,m,l,42,c)[:,0]
    else:
        return lambda t:skateboard.solution(dxi,ddxi,mu,m,[0,t],42,c)[1][0]
    
def time_p0(m, c, mu):
    """
    Fonction principale afin de rechercher le moment ou le mobile passe par
    le point p0 dependant de plusieurs parametres : m (la masse) et a
    (parametre de la fonction a).
    """
    y0 = [c, 0]
    if (c > -4):
        return None
    S = solve(mu,m,c,False)
    """
    t=np.linspace(0,100,200)
    plt.plot(t,[S(t) for t in t])
    plt.plot([t[0],t[-1]],[0,0])
    plt.xlabel("t")
    plt.savefig("/home/obi/public_html/plot.png")
    print("plotted")
    """

    b=0
    while S(b)<0:
        b=b-(c*mu*$) #FIXME robuste mais lent...
    root=brentq(S,0,b)#S(0)=c<0 ok ,donc il suffit de chercher l'autre borne
    return root #bon en réalité on est p-ê passé au-dessus de qqch...

def q4(m, c, mu):
    t = time_p0(m, c, mu)
    ds_t= skateboard.solution(dxi,ddxi,mu,m,[0,t],42,c)[1][1]
    if t != None:
        return t,(ds_t,2*ds_t)

def q5_time_and_dist(m, c, mu):
    v = q4(m, c, mu)
    if v != None:
        v_x, v_y = v[1]
        t = 2 * v_y / 9.81
        return t, t * v_x

def q5(m, c, mu):
    return q5_time_and_dist(m, c, mu)[1]

def q6(m, c, mu):
    """
    Fonction principale afin de rechercher le moment ou le mobile passe par
    le point p0 dependant de plusieurs parametres : m (la masse) et a
    (parametre de la fonction a).
    """
    res=q4(m,c,mu)
    if res != None:
        t1=res[0]
        t2= q5_time_and_dist(m, c, mu)[0]
        y0 = [0, 0] #FIXME pas zéro la vitesse
        #np.dot(deta(0), v) * deta(0)?
        sol = solve(mu,m,c,False)
        S=lambda t:sol(t)-pi
        b=1
        while S(b)<0:
            b=b+1
        t3=brentq(S,0,b)
        return t1+t2+t3


if __name__ == "__main__":
    """
    Methode "main", elle lancera la methode principale de time_p0 avec les
    valeurs donnees en ligne de commande.
    """
    c = float(sys.argv[1])
    mu = float(sys.argv[2])
    t = q4(1, c, mu)
    print(t)
