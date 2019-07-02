from __future__ import division, print_function
from visual import *
from visual.graph import *
scene.width = 800

##Constants
massAu = (79+118)*1.7e-27
massAlpha = (2+2)*1.7e-27
qAu = 2*1.6e-19
qAlpha = 79*1.6e-19
oofpez = 9e9     ## one-over-four-pi-epsilon-zero
deltat = 1e-23

##Objects
b = 5e-15
Au = sphere(pos=vector(0,0,0), radius=8e-15, color=color.yellow, make_trail=True)
Alpha = sphere(pos=vector(-1e-13,b,0), radius=2e-15, color=color.red, make_trail=True)

#Initial Values
pAu = massAu*vector(0,0,0)           
pAlpha = vector(1.043e-19,0,0)
t = 0

X_momentum = gdisplay(x=0,y=400,height=200)
px_Alpha = gcurve(color=color.red)
px_Au = gcurve(color=color.blue)
px_total = gcurve(color=color.cyan)
Y_momentum = gdisplay(x=0,y=600, height=200)
py_Alpha = gcurve(color=color.red)
py_Au = gcurve(color=color.blue)
py_total = gcurve(color=color.cyan)

##Calculation Loop
while t < 1e-20:
    rate(1000)
    vecAuAlpha = Au.pos - Alpha.pos
    fAuAlpha = (oofpez * qAu * qAlpha / mag2(vecAuAlpha)) * norm(vecAuAlpha)
    fAlphaAu = -fAuAlpha
    pAu = pAu + fAuAlpha * deltat
    pAlpha = pAlpha + fAlphaAu * deltat
    ptotal = pAu + pAlpha
    px_Alpha.plot(pos = (t, pAlpha.x))
    py_Alpha.plot(pos = (t, pAlpha.y))
    px_Au.plot(pos = (t, pAu.x))
    py_Au.plot(pos = (t, pAu.y))
    px_total.plot(pos = (t, ptotal.x))
    py_total.plot(pos = (t, ptotal.y))
    Au.pos = Au.pos + (pAu / massAu) * deltat
    Alpha.pos = Alpha.pos + (pAlpha/massAlpha)*deltat
    t = t + deltat

phat = norm(pAlpha)
theta = acos(phat.x)
print("theta1 = ",theta)

##Objects
b = 1e-14
posAu = vector(0,0,0)
posAlpha = vector(-1e-13,b,0)

##Values
pAu = massAu*vector(0,0,0)           
pAlpha = vector(1.043e-19,0,0)
t = 0

##Calculation Loop
while t < 1e-20:
    rate(1000)
    vecAuAlpha = posAu - posAlpha
    fAuAlpha = (oofpez * qAu * qAlpha / mag2(vecAuAlpha)) * norm(vecAuAlpha)
    fAlphaAu = -fAuAlpha
    pAu = pAu + fAuAlpha * deltat
    pAlpha = pAlpha + fAlphaAu * deltat
    ptotal = pAu + pAlpha
    posAu = posAu + (pAu / massAu) * deltat
    posAlpha = posAlpha + (pAlpha/massAlpha)*deltat
    t = t + deltat

phat = norm(pAlpha)
theta = acos(phat.x)
print("theta2 = ",theta)

##Objects
b = 1e-13
posAu = vector(0,0,0)
posAlpha = vector(-1e-13,b,0)

##Values
pAu = massAu*vector(0,0,0)           
pAlpha = vector(1.043e-19,0,0)
t = 0

##Calculation Loop
while t < 1e-20:
    rate(1000)
    vecAuAlpha = posAu - posAlpha
    fAuAlpha = (oofpez * qAu * qAlpha / mag2(vecAuAlpha)) * norm(vecAuAlpha)
    fAlphaAu = -fAuAlpha
    pAu = pAu + fAuAlpha * deltat
    pAlpha = pAlpha + fAlphaAu * deltat
    ptotal = pAu + pAlpha
    posAu = posAu + (pAu / massAu) * deltat
    posAlpha = posAlpha + (pAlpha/massAlpha)*deltat
    t = t + deltat

phat = norm(pAlpha)
theta = acos(phat.x)
print("theta3 = ",theta)

##from:http://hyperphysics.phy-astr.gsu.edu/hbase/Nuclear/impar.html
##b = (z1 * z2 * k * e^2 / m * v0^2) * cot(theta/2)
def fb(theta):
    return (qAu * qAlpha * oofpez / (mag(pAlpha) * (mag(pAlpha) / massAlpha) * tan(theta / 2)))

def fb2(theta):
    m = (massAlpha * massAu)/(massAlpha + massAu)
    v0 = pAlpha / massAlpha
    return (qAu * qAlpha * oofpez / (m * mag2(v0) * tan(theta / 2)))

b1 = fb(pi/4)
b2 = fb(pi/2)
b3 = fb(3 * pi/4)

##b4 = fb(2.59489129099)
##b5 = fb(2.11621090035)
##b6 = fb(0.251486441186)

##rb4 = fb2(2.59489129099)
##rb5 = fb2(2.11621090035)
##rb6 = fb2(0.251486441186)
print("-----------------------")
print("b1 = ",b1)
print("b2 = ",b2)
print("b3 = ",b3)
print("-----------------------")
##print(b4)
##print(b5)
##print(b6)
##print("-----------------------")
##print("reduced mass",rb4)
##print("reduced mass",rb5)
##print("reduced mass",rb6)
