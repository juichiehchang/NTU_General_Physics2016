from visual import *
from random import random
from visual.graph import *
import numpy as np

N = 50                              # number of the He atoms
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 # side length of the cubic container box

m, size = 4E-3/6E23, 310E-12        # He atoms mass, radius are made 10 times bigger for easiear collision
L_size=L-size                       
k, T = 1.38E-23, 298.0              # k = Boltzmann constant, and T = temperature in unit K

t, dt = 0, 0.5E-13                  # dt = 0.5E-13 second
vrms = (3*k*T/m)**0.5               # root mean square velocity at T
atoms = []                          # list to contain the 50 atoms

# histogram initialization, more on this see http://vpython.org/contents/docs/graph.html
dv = 10.
deltav = 100.
vdist = gdisplay(x=800, y=0, ymax = N*deltav/1000.,width=500, height=300, xtitle='v', ytitle='dN')
theory = gcurve(color=color.cyan) # plot the theoretical speed distribution dv = 10.
for v in arange(0.,3001.+dv,dv):
    theory.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*v**2*dv))
observation = ghistogram(bins=arange(0.,3000.,deltav),accumulate=1, average=1, color=color.red)


# initialization of display, setting up for the random position distribution and random velocity direction of atoms
scene = display(width=800, height=800,background=(0.2,0.2,0))
container = box(length = 2*L, height = 2*L, width = 2*L, opacity=0.2, color = color.yellow )

poslist = []
vlist = []

# Initialize atom position and velocity
for i in range(N):
    position = vector(-L_size+2*L_size*random(),-L_size+2*L_size*random(),-L_size+2*L_size*random())
    poslist.append(position)
    if i== N-1:
        atom = sphere(pos=position, radius = size, color=color.yellow, make_trail = True, retain = 600)
    else:
        atom = sphere(pos=position, radius = size, color=(random(), random(), random()))
        theta, phi = pi*random(), 2*pi*random()
    atom.m, atom.v = m, vector(vrms*sin(theta)*cos(phi), vrms*sin(theta)*sin(phi), vrms*cos(theta))
    vlist.append(atom.v)
    atoms.append(atom)
    
posnp = np.array(poslist)
rnp = np.full((50), size)
vnp = np.array(vlist)
lastcollision = np.zeros(50)
lasthit = np.zeros(50)
hitcount = 0
totalpath = 0
def vcollision(a1,a2):
    '''
    Function to find the velocities of atoms after each collision
    '''
    v1prime = a1.v - 2 * a2.m/(a1.m+a2.m) *(a1.pos-a2.pos) \
    * dot (a1.v-a2.v, a1.pos-a2.pos) / abs(a1.pos-a2.pos)**2
    v2prime = a2.v - 2 * a1.m/(a1.m+a2.m) *(a2.pos-a1.pos)  \
    * dot (a2.v-a1.v, a2.pos-a1.pos) / abs(a2.pos-a1.pos)**2
    return v1prime, v2prime
    
pressure = 0
counter = 0 
while True:
    t += dt
    counter += 1
    rate(1000)
    #calculate new positions for all atoms and plot histogram
    v=[]
    #### calculate new positions for atoms
    posnp = np.add(posnp, dt * vnp)
    lastcollision = np.add(lastcollision, dt * sqrt(sum(square(vnp), -1)))
    for i in range (N):
        atoms[i].pos = posnp[i]
        atoms[i].v = vnp[i]
        v.append(mag(atoms[i].v))
    observation.plot(data=v)

    #### find collisions between atoms, and then handle their collisions
    r = posnp-posnp[:,newaxis] # all pairs of atom-to-atom vectors
    rmag = sqrt(sum(square(r),-1)) # atom-to-atom scalar distances
    hit = less_equal(rmag,rnp+rnp[:,newaxis])-identity(N)
    hitlist = sort(nonzero(hit.flat)[0]).tolist() # i,j encoded as i*N+j
    for i in range (len(hitlist)):
        a1 = hitlist[i] / 50
        a2 = hitlist[i] % 50
        if(mag(posnp[a1] - posnp[a2]) < size * 2):
            if(lasthit[a1] == 1):
                totalpath += lastcollision[a1]
                hitcount += 1
            if(lasthit[a2] == 1):
                totalpath += lastcollision[a2]
                hitcount += 1
            lasthit[a1] = 1
            lasthit[a2] = 1
            lastcollision[a1] = 0
            lastcollision[a2] = 0
            colt = 0
            while(mag(posnp[a1] - posnp[a2]) < size * 2):
                posnp[a1] -= vnp[a1] * dt * 0.1
                posnp[a2] -= vnp[a2] * dt * 0.1
                colt += 1
            atoms[a1].pos = posnp[a1]
            atoms[a2].pos = posnp[a2]
            vnp[a1], vnp[a2] = vcollision(atoms[a1], atoms[a2])
            ##posnp[a1] += vnp[a1] * dt * (1 - 0.002 * colt)
            ##posnp[a2] += vnp[a2] * dt * (1 - 0.002 * colt)
            atoms[a1].pos = posnp[a1]
            atoms[a2].pos = posnp[a2]
            atoms[a1].v = vnp[a1]
            atoms[a2].v = vnp[a2]
    #### find collisions between atoms and walls, and then handle their collision
    outside = less_equal(posnp,size - L) # walls closest to origin
    vnp1 = vnp*outside
    vnp = vnp-vnp1+abs(vnp1) # force p component inward
    pressure += 2 * abs(sum(vnp1))
    outside = greater_equal(posnp,L-size) # walls farther from origin
    vnp1 = vnp*outside
    pressure += 2 * abs(sum(vnp1))
    vnp = vnp-vnp1-abs(vnp1) # force p component inward
    ##for i in range (N):
       ## pressure += mag(atoms[i].v - vnp[i]) * m / dt
    #### calculate the momentum transferred to the walls and obtain the pressure
    #### print the averaged pressure on the walls every 1000*dt
    if(counter % 1000 == 0):
        print("pressure")
        ##print(pressure)
        print((pressure * m) / t / 24 / L / L)
        print("free path")
        if(hitcount != 0):
            print(totalpath / hitcount)
        else:
            print("no hit")
