from math import sin, exp, cos, pi
import sys, skateboard

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

def time_p0(m, c, mu):
    """
    Fonction principale afin de rechercher le moment ou le mobile passe par
    le point p0 dependant de plusieurs parametres : m (la masse) et a
    (parametre de la fonction a).
    """
    y0 = xi(c)
    if (c > -4):
        return "Pas de solution"
    if (mu == 0):
        t = 10
    else:
        t = 10* (1/mu) * (1/abs(c))
    res = skateboard.solution(lambda x: dxi(x), lambda x: ddxi(x),
                               mu, m, t, 10001, y0)
    #(S(t), dS(t))
    n = 1
    while (n < 10001 and res[n][0] < 0):
        n += 1
    t0 = (t / 10000) * (n - 1)
    t = (t / 10000) * n
    #afficher_s(res)
    print("indice du tableau", n, "sur", 10000)
    print("bissection avec t0 = ", t0, " et t1 = ", t,"\net ", "S(t0) = ", res[n-1][0], " S(t1) = ", res[n][0])
    if (s(m, t, mu, y0) < 0):
        return "Pas de solution" #attention chercher la bonne borne sur le temps
    else:
        x = (t0 + t) / 2
        while (t - t0 > 1e-8 and abs(s(m, x, mu, y0)) > 1e-100):
            if (s(m, x, mu, y0) > 0):
                t = x
            else:
                t0 = x
            x = (t0 + t) / 2
        return x

def q4(m, c, mu):
    t = time_p0(m, c, mu)
    y0 = [c,0]
    s_t = s(m, t, mu, y0)
    ds_t = ds(m, t, mu, y0)
    return (t, (ds_t, (s_t +2)*ds_t)) #t t.q. S(t)=0 et le vecteur vitesse en ce temps

if __name__ == "__main__":
    """
    Methode "main", elle lancera la methode principale de time_p0 avec les
    valeurs donnees en ligne de commande.
    """
    c = float(sys.argv[1])
    mu = float(sys.argv[2])
    t = q4(1, c, mu)
    print(t)
