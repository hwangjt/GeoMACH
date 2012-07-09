from __future__ import division
import numpy, pylab
        

def cone(Lc, rz, ry, d, e, yy, zz):
    if d==0:
        d = 1
    if e==0:
        e = 1
    Q = (1 + 1/d**2 + 1/e**2)**0.5
    a = abs(Lc)/(1-1/Q)
    b = ry*(d**2+1)
    c = rz*(e**2+1)
    if Lc<0:
        yy *= -1
    s = -yy*b/a/d
    t = -zz*c/a/e
    den = (1/a**2+s**2/b**2+t**2/c**2)**0.5
    if Lc>0:
        x = a - 1/den
        y = -s/den
        z = -t/den
    else:
        x = abs(Lc) - a + 1/den
        y = s/den
        z = t/den
    return x, y, z

def cone2(Lc, rz, ry, s, t, i, j):
    q = s**2/ry**2 + t**2/rz**2
    a = 1/2/q*(abs(Lc)*q + (q*(Lc**2*q+4*(Lc**2*q+1)**0.5+4))**0.5)
    b = s/(s**2/ry**2-1/a**2)**0.5
    c = t/(t**2/rz**2-1/a**2)**0.5
    if Lc>0:
        s = s*(2*i-1)
        t = -t*j
    else:
        s = s*(1-2*i)
        t = t*(1-j)
    den = (1/a**2+s**2/b**2+t**2/c**2)**0.5
    if Lc>0:
        x = a - 1/den
        y = -s/den
        z = -t/den
    else:
        x = abs(Lc) - a + 1/den
        y = s/den
        z = t/den
    return x, y, z
