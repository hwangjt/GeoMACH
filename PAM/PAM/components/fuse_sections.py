from __future__ import division
import numpy, pylab



def circular(rz, ry, part, i):
    cos = numpy.cos
    sin = numpy.sin
    pi = numpy.pi
    if part==0:
        t = pi/2.0 - pi/4.0*i
    elif part==1:
        t = pi/4.0 - pi/2.0*i
    elif part==2:
        t = -pi/4.0 - pi/4.0*i
    z = rz*cos(t)
    y = ry*sin(t)
    return z, y

def rounded2(rz, ry, part, i):
    if part==0:
        z,y = circular(rz, ry, part, i)
    elif part==1 and i < 0.5:
        z,y = circular(rz, ry, part, i)
    else:
        z,y = rounded(rz, ry, part, i)
    return z,y

def rounded3(rz, ry, part, i):
    if part==2:
        z,y = circular(rz, ry, part, i)
    elif part==1 and i > 0.5:
        z,y = circular(rz, ry, part, i)
    else:
        z,y = rounded(rz, ry, part, i, fz=0.5, fy=0.2)
    return z,y        

def rounded(rz, ry, part, i, fz=0.5, fy=0.5):
    cos = numpy.cos
    sin = numpy.sin
    pi = numpy.pi
    if part==0:
        if i > 1-fz:
            t = pi/2.0 - pi/4.0*(i - (1-fz))/fz
            z = rz*(1-fz) + rz*fz*cos(t)
            y = ry*(1-fy) + ry*fy*sin(t)
        else:
            z = rz*i
            y = ry
    elif part==1:
        if i < fy/2.0:
            t = pi/4.0 - pi/4.0*i/(fy/2.0)
            z = rz*(1-fz) + rz*fz*cos(t)
            y = ry*(1-fy) + ry*fy*sin(t)
        elif i > 1-fy/2.0:
            t = -pi/4.0*(i - (1-fy/2.0))/(fy/2.0)
            z = rz*(1-fz) + rz*fz*cos(t)
            y = -ry*(1-fy) + ry*fy*sin(t)
        else:
            z = rz
            y = ry - 2*ry*i
    elif part==2:
        if i < fz:
            t = -pi/4.0 - pi/4.0*i/fz
            z = rz*(1-fz) + rz*fz*cos(t)
            y = -ry*(1-fy) + ry*fy*sin(t)
        else:
            z = rz*(1-i)
            y = -ry
    return z, y

def cone(Lc, rz, ry, d, e, i, j):
    Q = (1 + 1/d**2 + 1/e**2)**0.5
    a = abs(Lc)/(1-1/Q)
    b = ry*(d**2+1)
    c = rz*(e**2+1)
    s = b/a/d
    t = c/a/e
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
