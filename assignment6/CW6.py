from __future__ import division, print_function
from visual.graph import *
from visual.factorial import *
from math import*

kb = 1.38e-23
hbar = 1.05e-34
def microstates(q,N):
   """ returns the number of microstates Omega for a system with 
       q energy quanta and N oscillators 
       INPUT: 
          q: number of energy quanta
          N: number of oscilators
   """
   return (combin(q + N - 1, q))

def entropy(q,N):
   """ returns the entropy S for a system with 
       q energy quanta and N oscillators 
       INPUT: 
          q: number of energy quanta
          N: number of oscilators
   """
   return (kb * log(microstates(q, N)))

def temperature(q,N,w):
   """ returns the temperature T for a system with 
       q energy quanta and N oscillators of energy quantum \Delta E=\hbar w
       INPUT: 
          q: number of energy quanta
          N: number of oscilators
          w: oscillator angular frequency
   """
   #print("entropy")
   #print(entropy(q+1, N))
   #print(entropy(q-1, N))
   #print("chances")
   #print(log(microstates(q + 1, N)))
   #print(log(microstates(q - 1, N)))
   return (hbar * 2 * w / (entropy(q + 1, N) - entropy(q - 1, N)))


waygraph = gvbars(delta=0.7, color=color.red) # to make vertical bar graph
Ntotal = 6
N1 = 3
N2 = Ntotal - N1
qtotal = 4
q1 = 0

while q1 <= qtotal:
    q2 = qtotal - q1
    ways1 = microstates(q1, N1)
    ways2 = microstates(q2, N2)
    chance = ways1 * ways2
    waygraph.plot( pos=(q1, chance)  )
    q1 = q1+1

g1 = gdisplay(x = 0, y = 400, height = 300, width = 400)
waygraph1 = gvbars(delta=0.7, color=color.red) # to make vertical bar graph
Ntotal = 700
N1 = 400
N2 = Ntotal - N1
qtotal = 100
q1 = 0
mostchance = 0
mostq1 = 0
mostq2 = 0
halfq1 = 0
halfq2 = 0
A = []

while q1 <= qtotal:
    q2 = qtotal - q1
    ways1 = microstates(q1, N1)
    ways2 = microstates(q2, N2)
    chance = ways1 * ways2
    A.append(chance)
    waygraph1.plot( pos=(q1, chance)  )
    if chance >= mostchance:
       mostchance = chance
       mostq1 = q1
       mostq2 = q2
    q1 = q1+1

delta = 1e150
for i in range (0, 101):
   #print(A[i])
   if(abs(mostchance/2 - A[i]) <= delta):
      delta = abs(mostchance/2 - A[i])
      halfq1 = i
      halfq2 = qtotal - i
print("most probable: q1 = %d, q2 = %d" %(mostq1, mostq2))
print("half most probable: q1 = %d, q2 = %d" %(halfq1, halfq2))

g = gdisplay(x = 400, y = 400, height = 300, width = 400)
gc1= gcurve(color = color.cyan)
gc2 = gcurve(color = color.red)
gc3 = gcurve(color = color.blue)

Ntotal = 700
N1 = 300
N2 = Ntotal - N1
qtotal = 100
q1 = 0

mostchance = 0
mostq1 = 0
mostq2 = 0

while q1 <= qtotal:
    q2 = qtotal - q1
    S1 = log(microstates(q1, N1))
    S2 = log(microstates(q2, N2))
    S3 = S1 + S2

    if S3 >= mostchance:
       mostchance = S3
       mostq1 = q1
       mostq2 = q2

    gc1.plot(pos = (q1, S1))
    gc2.plot(pos = (q1, S2))
    gc3.plot(pos = (q1, S3))
    q1 = q1+1

print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("maximunm value of ln = ", mostchance)
print("most probable q1 = %d, q2 = %d" %(mostq1, mostq2))
q = 150
N = 105
ksi = 16
mAl = 26.981
w0 = sqrt(4*ksi/mAl/(1.7e-27))
#print("hbar * w0 = ", hbar * w0)
print("the temperature of Al = ", temperature(q, N, w0))
