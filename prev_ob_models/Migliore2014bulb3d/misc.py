# packages of various functions.
from copy import copy
from math import *

def distance(p, q):
    x = p[0] - q[0]
    y = p[1] - q[1]
    z = p[2] - q[2]
    return sqrt(x ** 2 + y ** 2 + z ** 2)

class Spherical:
    @staticmethod
    def to(p, center):
        rho = distance(p, center)
        
        p = copy(p); p[0] -= center[0]; p[1] -= center[1]; p[2] -= center[2]
        
        phi = atan2(p[1], p[0])
        try:
            theta = acos(p[2] / rho)
        except ZeroDivisionError:
            theta = acos(p[2] / 1e-8)
            
        return rho, phi, theta

    @staticmethod
    def xyz(rho, phi, theta, center):
        x = rho * cos(phi) * sin(theta) + center[0]
        y = rho * sin(phi) * sin(theta) + center[1]
        z = rho * cos(theta) + center[2]
        return [ x, y, z ]

def centroid(pts):
    x = 0.
    y = 0.
    z = 0.
    for p in pts:
        x += p[0]
        y += p[1]
        z += p[2]
    x /= len(pts)
    y /= len(pts)
    z /= len(pts)
    return [ x, y, z ]

# for elliptical coords
class Ellipse:
    
    def __init__(self, pos, axis):
        self.__pos = copy(pos)
        halfAxis = copy(axis)
        for i in range(3): halfAxis[i] /= 2.
        self.__inverse = halfAxis[0] < halfAxis[1]
        self.__halfAxis = halfAxis
        
        # eccentricity
        a = 0; b = 1;
        if halfAxis[a] < halfAxis[b]: b = 0; a = 1
        
        self.__eccen = sqrt(halfAxis[a] ** 2 - halfAxis[b] ** 2) / halfAxis[a]

    def R(self, phi): return self.__halfAxis[0] / sqrt(1 - (self.__eccen * sin(phi)) ** 2)

    # from elliptical to cartesian
    def toXYZ(self, h, lamb, phi):
        N = self.R(phi)
        XYProj = (N + h) * cos(phi)
        p = [ XYProj * cos(lamb),
              XYProj * sin(lamb),
              ((1 - self.__eccen ** 2) * N + h) * sin(phi) ]
        if self.__inverse: aux = p[0]; p[0] = p[1]; p[1] = aux
        for i in range(3): p[i] += self.__pos[i]
        return p

    # from cartesian to elliptical
    def toElliptical(self, pt):
        x = pt[0] - self.__pos[0]
        y = pt[1] - self.__pos[1]
        z = pt[2] - self.__pos[2]
        
        if self.__inverse: aux = y; y = x; x = aux
            
        lamb = atan2(y, x)

        p = sqrt(x ** 2 + y ** 2)
        try:
            phi = atan(z / ((1 - self.__eccen ** 2) * p))
        except ZeroDivisionError:
            phi = atan(z / 1e-8)

        MAXIT = int(1e+4)
        for i in range(MAXIT):
            phi1 = phi
            N = self.R(phi)
            h = p / cos(phi) - N
            try:
                phi = atan(z /  ((1 - self.__eccen ** 2 * N / (N + h)) * p))
            except ZeroDivisionError:
                phi = atan(z / 1e-8)
            if abs(phi - phi1) < 1e-8: break
        return h, lamb, phi

    def getZ(self, pt):            
        x = pt[0]; y = pt[1]
        try:
            return self.__pos[2] + self.__halfAxis[2] * sqrt(1 - ((x - self.__pos[0]) / self.__halfAxis[0]) ** 2 - ((y - self.__pos[1]) / self.__halfAxis[1]) ** 2)
        except ValueError:
            return None

    def normalRadius(self, pt):
        r = 0.
        for i in range(3): r += ((pt[i] - self.__pos[i]) / self.__halfAxis[i]) ** 2
        return r


# return versor between two points
def versor(p, q):
        d = distance(p, q)
        v = [ 0 ] * 3
        for i in range(3): v[i] = (p[i] - q[i]) / d
        return v
        
      
def ellipseLineIntersec(u, p, center, axis): 
    p = copy(p)
    A = 0.
    B = 0.
    C = -1
    v = []

    for i in range(3):
        _ax = axis[i] / 2
        A += u[i] ** 2 / _ax ** 2
        B += 2 * u[i] * (p[i] - center[i]) / _ax ** 2
        C += (p[i] - center[i]) ** 2 / _ax ** 2
          
    delta = B ** 2 - 4 * A * C
    t0 = (-B + sqrt(delta)) / (2 * A)
    t1 = (-B - sqrt(delta)) / (2 * A)
    t = t1
    if abs(t0) < abs(t1): t = t0
    for i in range(3): p[i] += t * u[i]
    return p

def ellipseIntersec(p, center, axis):
    e = Ellipse(center, axis)
    h, lamb, phi = e.toElliptical(p)
    u = [ sin(lamb) * cos(phi), cos(lamb) * cos(phi), sin(phi) ]
    p = ellipseLineIntersec(u, p, center, axis)
    return p, lamb, phi      


    
    
