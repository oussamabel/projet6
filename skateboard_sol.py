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
def solve(mu,m,y0,is_list):
    if is_list:
        return lambda l:skateboard.solution(dxi,ddxi,mu,m,l,42,y0)[:,0]
    else:
        return lambda t:skateboard.solution(dxi,ddxi,mu,m,[0,t],42,y0)[1][0]
    
def time_p0(m, c, mu):
    """
    Fonction principale afin de rechercher le moment ou le mobile passe par
    le point p0 dependant de plusieurs parametres : m (la masse) et a
    (parametre de la fonction a).
    """
    y0 = [c, 0]
    if (c > -4):
        return None
    S = solve(mu,m,[c,0],False)
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
        b=b-((mu+0.01)*c)
    root=brentq(S,0,b)#S(0)=c<0 ok ,donc il suffit de chercher l'autre borne
    return root #bon en réalité on est p-ê passé au-dessus de qqch...

def q4(m, c, mu):
    t = time_p0(m, c, mu)
    ds_t= skateboard.solution(dxi,ddxi,mu,m,[0,t],42,[c,0])[1][1]
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

def q6(m, c, mu1,mu2):
    """
    Fonction principale afin de rechercher le moment ou le mobile passe par
    le point p1 dependant de plusieurs parametres : m (la masse) et a
    (parametre de la fonction a).
    """
    v = q4(m, c, mu1)
    if v != None:
        v_x, v_y = v[1]
        t1=v[0]
        t2= q5_time_and_dist(m, c, mu1)[0]
        ds_0=(v_x-(v_y*0.2))/1.04
        y0 = [0, ds_0] 
        sol = solve(mu2,m,y0,False)
        S=lambda t:sol(t)-pi
        b=0
        while S(b)<0:
            b=b-((mu2+0.01)*c)
        t3=brentq(S,0,b)
        return t1+t2+t3

def further_than_p2(c, m, mu1, mu2, t, incr):
    try:
         v = q4(m, c, mu)
         v_x, v_y = v[1]
         t2= q5_time_and_dist(m, c, mu)[0]
         t1=v[0]
    except:
        return False
    ds_0 = (v_x - 0.2 * v_y) / 1.04
    y0 = [0, ds_0]
    S = solve(mu2,m,[0,ds_0],False)
    while True:
        s1 = S(t)
        s2 = S(t+incr)
        print(s1, s2)
        if s1 <= s2 and s2 < 2*pi: #si on avance sur la courbe "vers" p2
            t += incr
        elif s1 < 2*pi and s1 > s2: #si on a pas atteint p2
            return False
        else: #si on a atteint p2
            return True

def q7(m, mu1, mu2):
    borne_sup =-4
    pas=-5
    borne_inf=-4
    while not(further_than_p2(borne_inf, m, mu1, mu2, 500, 500)):
        borne_inf += pas
        print(borne_inf)
    f = lambda x: 1 if further_than_p2(x, m, mu1, mu2, 100, 100) else -1
    return brentq(f, borne_inf, borne_sup)


if __name__ == "__main__":
    """
    Methode "main", elle lancera la methode principale de time_p0 avec les
    valeurs donnees en ligne de commande.
    """
    c = float(sys.argv[1])
    mu1 = float(sys.argv[2])
    mu2=float(sys.argv[3])
    #t = q4(1, c, mu1)
    t=q6(1,c,mu1,mu2)
    print(t)
