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
    b=1
    while S(b)<0:
        b=b+1 #FIXME ou 0.001 si ça marche pas
    root=brentq(S,0,b)#S(0)=c<0 ok ,donc il suffit de chercher l'autre borne
    return root #bon en réalité on est p-ê passé au-dessus de qqch...

def q4(m, c, mu):
    t = time_p0(m, c, mu)
    ds_t= skateboard.solution(dxi,ddxi,mu,m,[0,t],42,c)[1][1]
    if t != None:
        return t,(ds_t,2*ds_t)

def q5(m, c, mu,  time = False):
    v = q4(m, c, mu)[1]
    v_x = v[0]
    v_y = v[1]
    t = 2 * v_y / 9.81
    d = v_x * t
    if (time):
        return (d, t)
    else:
        return d

def q6(m, c, mu): #a adapter
    (d, t) = q5(m, c, mu, True)
    """
        Fonction principale afin de rechercher le moment ou le mobile passe par
        le point p0 dependant de plusieurs parametres : m (la masse) et a
        (parametre de la fonction a).
        """
    y0 = [0, 0]
    if (c > -4):
        return "Pas de solution"
    if (mu == 0):
        t = 10
    else:
        t = 15 * (1 / mu)
    res = skateboard.solution(lambda x: deta(x), lambda x: ddeta(x),
                              mu, m, t, 10001, y0)
    # (S(t), dS(t))
    n = 1
    while (n < 10001 and res[n][0] < 0):
        n += 1
    t0 = (t / 10000) * (n - 1)
    t = (t / 10000) * n
    # afficher_s(res)
    # print("indice du tableau", n, "sur", 10000)
    # print("bissection avec t0 = ", t0, " et t1 = ", t,"\net ", "S(t0) = ", res[n-1][0], " S(t1) = ", res[n][0])
    if (s(m, t, mu, y0) < 0):
        return "Pas de solution"  # attention chercher la bonne borne sur le temps
    else:
        x = (t0 + t) / 2
        while (t - t0 > 1e-8 and abs(s(m, x, mu, y0)) > 1e-8):
            if (s(m, x, mu, y0) > 0):
                t = x
            else:
                t0 = x
            x = (t0 + t) / 2
        return x


if __name__ == "__main__":
    """
    Methode "main", elle lancera la methode principale de time_p0 avec les
    valeurs donnees en ligne de commande.
    """
    c = float(sys.argv[1])
    mu = float(sys.argv[2])
    t = q4(1, c, mu)
    print(t)
