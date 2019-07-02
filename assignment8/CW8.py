from __future__ import division, print_function
from visual import *
from visual.graph import *
scene.width = scene.height = 800

## constants
oofpez = 9e9 # stands for One Over Four Pi Epsilon-Zero
qe = 1.6e-19  
s = 4e-11   
scalefactor = 3e-20 # you may need to change this

p1 = sphere(pos=vector(0,s/2,0), radius=1e-11, color=color.red)
q1 = qe
p2 = sphere(pos=vector(0,-s/2,0), radius=1e-11, color=color.blue)
q2 = -qe

proton = sphere(pos = vector(-3 * s, 0, 0), radius = 2e-11, color = color.cyan, make_trail = true)
mproton = 1.6726219e-27
vproton = vector(0, 0, 0)
qproton = q1
## for plotting

egraphs=gdisplay(x=600, width=600, height = 600)   ## move graph so it's not on top of scene
Ug = gcurve(color=color.yellow)
Kg = gcurve(color=color.cyan)
KUg = gcurve(color=color.magenta)

R = 3e-10
theta = 0
while(theta < 2 * pi):
    #rate(300)
    position = R * vector(cos(theta), cos(pi/2 - theta), 0)
    E1 = oofpez * q1 * norm(position - p1.pos) / mag2(position - p1.pos)
    E2 = oofpez * q2 * norm(position - p2.pos) / mag2(position - p2.pos)
    arrow(pos = position, color = color.orange, axis = scalefactor * (E1 + E2))
    position = R * vector(0, cos(theta), cos(pi / 2 - theta))
    E1 = oofpez * q1 * norm(position - p1.pos) / mag2(position - p1.pos)
    E2 = oofpez * q2 * norm(position - p2.pos) / mag2(position - p2.pos)
    arrow(pos = position, color = vector(0.5, 0.8, 0.2), axis = scalefactor * (E1 + E2))
    theta += pi / 6

t = 0
dt = 1e-17
counter = 0
last = 0
while(counter < 6):
    rate(300)
    t += dt
    F1 = oofpez * q1 * qproton * norm(proton.pos - p1.pos) / mag2(proton.pos - p1.pos)
    F2 = oofpez * q2 * qproton * norm(proton.pos - p2.pos) / mag2(proton.pos - p2.pos)
    vproton += (F1 + F2) * dt / mproton
    if(last > 0 and vproton.y < 0):
        counter += 1
    last = vproton.y
    proton.pos += vproton * dt
    K = mproton * mag2(vproton) / 2
    Kg.plot(pos = (t, K))
    U = (oofpez * q1 * qproton / mag(proton.pos - p1.pos)) + (oofpez * q2 * qproton / mag(proton.pos - p2.pos)) + (oofpez * q1 * q2 / mag(p1.pos - p2.pos))
    Ug.plot(pos = (t, U))
    KUg.plot(pos = (t, U + K))
    
