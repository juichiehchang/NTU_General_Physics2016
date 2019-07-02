from __future__ import division, print_function
from visual import *
scene.width = 1000
scene.height = 700
R = 0.08
Q = 50e-9
Nq = 25   
dtheta = 2*pi/Nq
dQ = Q/Nq
k = 9e9
scalefactor = 5e-5
## source charges
ring(pos=vector(0,0,0), radius=R, color=color.red, thickness=0.005)
sources = []
obs = []
offaxis = []
angles=arange(0,2*pi,dtheta)
for theta in angles:
    a = sphere(pos=vector(0, R*cos(theta),R*sin(theta)), radius=0.01, color=color.cyan, q=dQ)
    sources.append(a)

for i in range(-8, 9):
    position = vector(R * i / 4, 0, 0)
    Enet = vector(0, 0, 0)
    for j in range(0, Nq):
        Enet += (9 * 1e9 * sources[j].q / mag(position - sources[j].pos)) * (position - sources[j].pos)
    spot = arrow(pos = position, axis = Enet * scalefactor, opacity = 0.5)
    obs.append(spot)

for i in range(0, 7):
    position = vector(R / 2, 0, -2 * R + (R * 2 / 3) * i)
    Enet = vector(0, 0, 0)
    for j in range(0, Nq):
        Enet += (9 * 1e9 * sources[j].q / mag(position - sources[j].pos)) * (position - sources[j].pos)
    spot = arrow(pos = position, axis = Enet * scalefactor, color = color.cyan, opacity = 0.5)
    offaxis.append(spot)

mouseposition = scene.waitfor('click').pos     ## wait for mouse click, save the mouse position
electron = sphere(pos = mouseposition, radius = R / 20, q = -1.6e-19, color = vector(0.7, 0.7, 0), make_trail = True, v = vector(0, 0, 0), m = 9.10938356e-31)
t = 0
dt = 1e-10
while True:
    rate(500)
    force = vector(0, 0, 0)
    for i in range(0, Nq):
        vec = electron.pos - sources[i].pos
        force += (k * sources[i].q * electron.q / mag(vec)) * (vec)
    electron.v += (force / electron.m) * dt
    electron.pos += electron.v * dt
