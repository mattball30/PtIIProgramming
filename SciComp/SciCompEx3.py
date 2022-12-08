 # defining the euler method in order to time develop our system: 
import numpy as np
import matplotlib.pyplot as plt
import math


# defining our rate constants:
def K_fR15(ur): 
  return 26000 * math.exp(-1.68*ur)
def K_fR16(ur): 
  return 730 * math.exp(-1.72*ur)
def K_uR15(ur): 
  return 6e-2 * math.exp(0.95*ur)
def K_uR16(ur):  
  return 7.5e-4 * math.exp(1.2*ur)

#defining our new class: 
class Cell:
  def __init__(self, reaction, k):
    # reaction is passed as an array of 2 elements - the reactants and the products
    # k = int 
      self.reaction = reaction
      self.k = k
      self.reactants = reaction.split('->')[0].split('+')
      self.products = reaction.split('->')[1].split('+')

def euler(f_user, U0, n, T): 
  f = lambda u , t: np.asarray(f_user(u,t))
  # this equation will solve u' = f(u,t), with n steps up to T

  t = np.zeros(n+1) 
  t[0] = 0 
  dt = T / n
  # now we generate our solution matrix - U, with the dimensions determined by the number of initial conditions (ie number of reagents) 
  neq = len(U0) 
  u = np.zeros((n+1, neq)) 
  u[0, :] = U0

  for i in range(n): 
    u[i+1, :] = u[i, :] + f(u[i, :],t[i]) * dt
    t[i+1] = t[i] + dt
  return u, t



# installising the reaction (and setting the initial value for the concentration of urea):

a = Cell('D->I', K_fR15)
b = Cell('I->D',K_uR15)
c = Cell('I->N', K_fR16)
d = Cell('N->I', K_uR16)

# defining the different rates of change of concentrations:
def dD(u): 
  return b.k(ur) * u[1] - a.k(ur) * u[0]
def dI(u): 
  return a.k(ur) * u[0] - b.k(ur) * u[1] - c.k(ur) * u[1] + d.k(ur) * u[2]
def dN(u): 
  return c.k(ur) * u[1] - d.k(ur) * u[2]

# now we define our input functions: 
def f(u,t): 
  return [dD(u),dI(u),dN(u)] 

# we have returned this as an array, with each index corresponding to the derivative of each species
# u[0] = D, u[1] = I, u[2] = N

# defining the parameters for the set of reactions: 
T = 1
n = 10000
U0 = [1,0,0]
jump = 0.25
urea = np.arange(0,10,jump)
u_eq = np.zeros((len(urea), len(U0)))

# generating the different equilibrium values:
for i in urea:
  ur = urea[int(i/jump)]
  u, t = euler(f, U0, n, T)
  u_eq[int(i/jump),0] = u[n,0]
  u_eq[int(i/jump),1] = u[n,1]
  u_eq[int(i/jump),2] = u[n,2]



# mapping the reactants to the appropriate columns in our u matrix:
# this is D:
vars()[a.reactants[0]] = u[:,0]
# this is I: 
vars()[a.products[0]] = u[:,1]
# this is N:
vars()[c.products[0]] = u[:,2]

# plotting the resulting function
plt.plot(urea, u_eq[:,0], 'r', marker = 'o')
plt.plot(urea, u_eq[:,1], 'b', marker = 'o')
plt.plot(urea, u_eq[:,2], 'g', marker = 'o')
plt.xlabel('Concentration of urea / M')
plt.ylabel('Equilibrium fraction of species')
plt.title('Equilbrium fractions of protein states against concentration of urea')
plt.legend(['D','I','N'])
plt.show()

# now moving onto the second task: 
# we now have a new set of reactions - with 5 variables in total: 
# we have now set up the appropriate code, so we can repeat with the oreganator: 
Ore = [Cell('A+Y->X+P',1.34), Cell('X+Y->P',1.6*10**9), Cell('B+X->2X+Z',8*10**3), Cell('2X->Q',4*10**7), Cell('Z->Y', 1)]

# defining our new rates of change: 
def dA(u): 
  return - Ore[0].k * u[0] * u[5]
def dB(u):
  return - Ore[2].k * u[1] * u[4]
def dP(u): 
  return Ore[0].k * u[0] * u[5] + Ore[1].k * u[4] * u[5]
def dQ(u):
  return Ore[3].k * u[4]**2
def dX(u): 
  return Ore[0].k * u[0] * u[5] - Ore[1].k * u[4] * u[5] + Ore[2].k * u[1] * u[4] - 2 * Ore[3].k * u[4]**2
def dY(u): 
  return - Ore[0].k * u[0] * u[5] - Ore[1].k * u[4] * u[5] + Ore[4].k * u[6]
def dZ(u): 
  return - Ore[4].k * u[6] + Ore[2].k * u[1] * u[4]

# defining our differential equations for the euler solver: 
def f(u,t): 
  return [dA(u), dB(u), dP(u), dQ(u), dX(u), dY(u), dZ(u)]


# defining the parameters for the set of reactions: 
T = 80
n = 80_000_000
U0 = [0.06,0.06,0,0,10**-9.8,10**-6.52, 10**-7.32]

u, t = euler(f, U0, n, T)


# once again assigning our values: 
vars()[Ore[0].products[0]] = u[:,4]
vars()[Ore[0].reactants[1]] = u[:,5]
vars()[Ore[4].reactants[0]] = u[:,6]

plt.plot(t, u[:,4], 'b', marker='x',markevery=8_000_000)
plt.plot(t, u[:,5], 'g', marker='x',markevery=8_000_000)
plt.plot(t, u[:,6], 'r', marker='x',markevery=8_000_000)
plt.xlabel('t / s')
plt.ylabel('Concentration / M')
plt.title('Concentrations of the relevant species in the Oreganator over time')
plt.legend(['X','Y','Z'])
plt.yscale("log")
plt.show()

