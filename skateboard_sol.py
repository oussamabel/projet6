from math import sin, exp, cos, pi
import sys, trajectoire

def xi(s):
    """
    Fonction xi(s).
    """
    return (s, (1/2)*(s+2)*(s+2))

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

def s(m, t, mu, y0):
    """
    Fonction dont on va rechercher les racines (en fonction du temps) dependant
    de plusieurs parametres : m (la masse), mu (le coefficient de frottements)
    et y0 (le point de depart).
    La racine correspondra au moment ou le mobile passe par le point p1.
    """
    res = trajectoire.solution(lambda x: dxi(x), lambda x: ddxi(x),
                               mu, m, t, 2, y0)
    return res[1][0]

def time_p0(m, c, mu):
    """
    Fonction principale afin de rechercher le moment ou le mobile passe par
    le point p1 dependant de plusieurs parametres : m (la masse) et a
    (parametre de la fonction a).
    """
    y0 = xi(c)
    if(mu == 0):
        if (c > -4):
            return "Pas de solution"
        else :
            t = 10 * (-c)
    else:
        t = 100 * (-c)/mu
    res = trajectoire.solution(lambda x: dxi(x), lambda x: ddxi(x),
                               mu, m, t, 10001, y0)
    n = 1
    while (n < 10001 and res[n][0] < 0):
        n += 1
    t0 = (t / 10000) * (n - 1)
    t = (t / 10000) * n
    print(n)
    if (s(m, t, mu, y0) < 0):
        return "Pas de solution"
    else:
        x = (t0 + t) / 2
        while (t - t0 > 1e-8 and abs(s(m, x, mu, y0)) > 1e-100):
            if (s(m, x, mu, y0) > 0):
                t = x
            else:
                t0 = x
            x = (t0 + t) / 2
        return x

if __name__ == "__main__":
    """
    Methode "main", elle lancera la methode principale de time_p1 avec les
    valeurs donnees en ligne de commande.
    """
    t1 = time_p0(1, -3.745, 0)
    print(t1)