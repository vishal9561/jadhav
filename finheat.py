# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:49:17 2016

@author: Sagar
"""
from __future__ import division
from math import sqrt
#Variable declaration:
'''h1 = 13.0                       #Air-side heat transfer coefficient (W/m^2.K)
A = 1.0                         #Base wall area (m^2)
L = 2.5/100                     #Length of steel fins (m)
L2 = 1.5/10**3                  #Length of steel wall (m)
k = 13.0                        #Thermal conductivity of fin (W/m.K)
k1 = 38.0                       #Thermal conductivity of steel wall (W/m.K)
h2 = 260.0                      #Water side heat transfer coefficient (W/m^2.K)
T4 = 19.0                       #Air temperature (°C)
T1 = 83.0                       #Water temperature (°C)
t = 1.3/10**3                   #Thickness of steel fins (m)
w = 1.0                         #Width of wall (m)
S = 1.3/100                     #Fin pitch(m)'''

h1=raw_input("give value of h1:") 
A=raw_input("give value of A:")
L=raw_input("give value of L:")
L2=raw_input("give value of L1:")
k=raw_input("give value of K:")
k1=raw_input("give value of K1:")
h2=raw_input("give value of h2:")
T4=raw_input("give value of T4:")
T1=raw_input("give value of T1:")
t=raw_input("give value of t:")
w=raw_input("give value of w:")
S=raw_input("give value of s:")
h1=float(h1)
A=float(A)
L=float(L)
L2=float(L2)
k=float(k)
k1=float(k1)
h2=float(h2)
T4=float(T4)
T1=float(T1)
t=float(t)
w=float(w)
S=float(S)

#Calculation:
R1 = 1/(h1*A)                   #Air resistance (°C/W) (part 1)
R2 = L2/(k1*A)                  #Conduction resistance (°C/W)
R3 = 1/(h2*A)                   #Water resistance (°C/W)
Rt = (R1+R3)                    #Total resistance (°C/W) (part 2)
Q = (T1-T4)/Rt                  #Total heat transfer (W)
Nf = 1/S                        #Number of fins (part 3)
Lbe = w - Nf*t                  #Unfinned exposed base surface
Abe = w*Lbe                     #Exposed base surface area (m^2)
Lc = L+t/2                      #Corrected length (m)
Ap = Lc*t                       #Profile area (m^2)
Af = 2*w*Lc                     #Fin surface area (m^2)
Bi = h1*(t/2)/k1                #Biot number
a = sqrt(Lc**3*h1/(k*Ap))       #Abscissa of the fin efficiency
#
nf = 0.88                       #Fin efficiency-from correlation
Rb = 1/(h1*Abe)                 #Air thermal resistance of base wall (°C/W)
Rf = 1/(h1*Nf*Af*nf)            #Air thermal resistance of fins (°C/W)
RT1 = 1/(1/Rb+1/Rf)             #Total outside resistance of the fin array (°C/W)
Rt3 = RT1+R3                    #Total resistance on air side fins (°C/W)
Qt = (T1-T4)/round(Rt3,5)       #Heat transfer rate on air side fins (W)
I = (Qt/Q - 1.0)*100.0              #Percent increase in heat transfer rate to air side fins (W)
A = sqrt(Lc**3*h2/(k1*Ap))      #Abscissa of the new fin efficiency (part 4)
#
nf2 = 38.0                      #New fin efficiency-from correlation
Rb2 = 1/(h2*Abe)                #Thermal resistance of base wall (°C/W)
Rf2 = 1/(h2*Nf*Af*nf2)          #Thermal resistance of fins (°C/W)
Rt4 = 1/(1/Rb2+1/Rf2)           #Total resistance of the finned surface (°C/W)
Rt5 = R1+Rt4                    #Total resistance on water side fins (°C/W)
QT1 = (T1-T4)/Rt5               #Heat transfer rate on water side fins (W)
I2 = (QT1/Q - 1)*100            #Percent increase in heat transfer rate to water side fins (W)

#Result:
if (R2<R1 or R2<R3):
    print "1. The conduction resistance may be neglected."
else:
    print "1. The conduction resistance can not be neglected."
print "2. The rate of heat transfer from water to air is :",round(Q,1)," W ."
print "3. The percent increase in steady-state heat transfer rate by adding fins to the air side of the plane wall is :",round(I,1)," % ."
print "4. The percent increase in steady-state heat transfer rate by adding fins to the water side of the plane wall is :",round(I2,1)," % ."


