from __future__ import division
import numpy, pylab



def isContained(t, t1, t2):
    return t1 <= t and t <= t2

def nMap(t, t1, t2):
    return (t-t1)/(t2-t1)

def rounded2(rz, ry, t):
    if isContained(t,0,1):
        z,y = circular(rz, ry, t)
    else:
        z,y = rounded(rz, ry, t)
    return z,y

def rounded3(rz, ry, t):
    if isContained(t,1,2):
        z,y = circular(rz, ry, t)
    else:
        z,y = rounded(rz, ry, t)
    return z,y

def rounded4(rz, ry, t):
    if isContained(t,0,1):
        z,y = rounded(rz, ry, t, t1=0.5, t2=0.5)
    else:
        z,y = rounded(rz, ry, t)
    return z,y 

def rounded5(rz, ry, t):
    z,y = rounded(rz, ry, t, t1=0, t2=0.6)
    return z,y 

def rounded6(rz, ry, t):
    z,y = rounded(rz, ry, t, t1=0, t2=0.55)
    return z,y 

def circular(rz, ry, t):
    z = rz*numpy.cos(t*numpy.pi)
    y = ry*numpy.sin(t*numpy.pi)
    return z, y

def rounded(rz, ry, t, t1, t2):
    t1 /= 2.0
    t2 /= 2.0
    cos = numpy.cos
    sin = numpy.sin
    tan = numpy.tan
    pi = numpy.pi
    z0 = ry*tan((0.5-t2)*pi)
    y0 = rz*tan(t1*pi)
    sz = rz - z0
    sy = ry - y0
    if isContained(t,0,t1):
        z = rz
        y = rz*tan(t*pi)
    elif isContained(t,t1,t2):
        t = nMap(t,t1,t2)/2.0
        z = z0 + sz*cos(t*pi)
        y = y0 + sy*sin(t*pi)
    elif isContained(t,t2,1-t2):
        z = ry*tan((0.5-t)*pi)
        y = ry
    elif isContained(t,1-t2,1-t1):
        t = nMap(t,1-t2,1-t1)/2.0 + 0.5
        z = -z0 + sz*cos(t*pi)
        y = y0 + sy*sin(t*pi)
    elif isContained(t,1-t1,1+t1):
        z = -rz
        y = rz*tan((1-t)*pi)
    elif isContained(t,1+t1,1+t2):
        t = nMap(t,1+t1,1+t2)/2.0 + 1
        z = -z0 + sz*cos(t*pi)
        y = -y0 + sy*sin(t*pi)
    elif isContained(t,1+t2,2-t2):
        z = -ry*tan((1.5-t)*pi)
        y = -ry
    elif isContained(t,2-t2,2-t1):
        t = nMap(t,2-t2,2-t1)/2.0 + 1.5
        z = z0 + sz*cos(t*pi)
        y = -y0 + sy*sin(t*pi)
    elif isContained(t,2-t1,2):
        z = rz
        y = rz*tan(t*pi)  
    return z,y
        

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
