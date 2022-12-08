# In this exercise we are creating a program to simulate the interactions between particles in a system based on 2 different potentials.
# importing the relevant packages: 
from numpy.core.memmap import ndarray
import numpy as np 
import random

# defining the distance between 2 vectors and the length of a vector:   
def r_ab(a, b): 
    r = a - b
    d = 0
    for i in range(len(r)): 
          d += r[i]**2
    return np.sqrt(d)

def length(x):
#Returns the length of the vector."""
  a = 0
  for i in range(len(x)): 
      a += x[i]**2
  return np.sqrt(a)

# defining the 2 different potentials for the vectors : 
def LJ(a,b):  
  return 4*((1/r_ab(a,b))**12-(1/r_ab(a,b))**6)

def morse(a,b,sigma): 
  return (1-np.exp(r_ab(a,b) - sigma))**2

# defining a new class for the system, which is an n by 3 array containing the coordinates of the relevant particles
class System(np.ndarray):      
  def __new__(cls, n):
    pos = np.array(np.zeros((n,3)), dtype=float)
    x = np.ndarray.__new__(cls, shape=(n,3), dtype=float, buffer=pos)
    return x

  #generating our starting positions for the optimisation
  def start(self): 
    n = len(self)
    for i in range(n): 
        self[i] = vec((random.uniform(-1.3,1.3), random.uniform(-1.3,1.3), random.uniform(-1.3,1.3)))
    return self
  
  # generating the energies of the system based on the ij interactions. 
  def U(self, pot):
    U = 0 
    n = len(self)
    for i in range(0,n): 
        for j in range(i+1,n): 
            U += pot(self[i],self[j])

    return U 

  # defining the infinitesimal change in the energy
  def dU(self, pot, delta, atom, coord): #here coord, x = 0, y = 1, z = 2 and atom refers to the atom in question 

    # copying the left matrix here so we can act on it without changing the matrix elsewhere in the system.
    dummy_array = self.copy()
    dummy_array[atom][coord] = dummy_array[atom][coord] - delta
    U0 = dummy_array.U(pot) 

    # applying a small change to this matrix
    dummy_array[atom][coord] = dummy_array[atom][coord] + 2*delta

    # calculating the new potential energy and returning the difference: 
    U1 = dummy_array.U(LJ)
    
    return (U1 - U0)/(2*delta)

  # defining the corresponding infinitesimal change in the gradient 
  def dR(self, pot, delta, atom): 
        dR = [self.dU(pot, delta, atom, 0), self.dU(pot, delta, atom, 1), self.dU(pot, delta, atom, 2)] 
        return dR 

  # defining the optimization method: 
  def opt(self, pot, lam ,delta): 
    # these are the 2 parameters selected which give the best convergence / performance
    # lam = 0.001 
    # delta = 1e-8
    
    for i in range(len(self)): 
        print('inital particle ' + str(i) + ' position = ', self[i])

    # this creates a matrix to store the gradient vector for each atom
    x = np.zeros((len(self),3))
    for i in range(len(self)):
        x[i] = self.dR(LJ, delta, i) 

    # this modulates the length of our gradient vectors to give convergence around the equilibrium bond length
        if length(x[i]) > 1: 
            x[i] = x[i] / length(x[i])

    # this is creating the initial condition for the while loop: 
    sum = 0
    for i in range(len(self)):
        sum += length(x[i])
    
    # i.e. this is the check for convergence by looking at the length of the gradient vectors and ensuring it is small: 
    while sum > 3e-7: 
        x = np.zeros((len(self),3))
        sum = 0
        for i in range(len(self)):
            x[i] = self.dR(LJ, delta, i)
            if length(x[i]) > 1: 
                x[i] = x[i] / length(x[i])

    # this adjusts the matrix with the small change in the gradient 
        self = self - lam * x

    # this is the check for convergence 
        for i in range(len(self)):
            sum += length(x[i]) 

    # finally shifting everything to the origin:
    dummy = self.copy() 
    for i in range(len(self)): 
        self[i] = self[i] - dummy[0]
    
    # outputting our data as an XYZ file: 
    with open("output.txt", "w") as txt_file:
        txt_file.write(str(len(self))+'\n')
        txt_file.write('The lowest energy conformation for '+str(len(self))+' particles with the Lennard-Jones potential:'+'\n')
        for i in range(len(self)):
            txt_file.write('Particle' + str(i) + '  ' + str(self[i][0]) + ' ' + str(self[i][1]) + ' ' + str(self[i][2])  + "\n") # works with any number of elements in a line

    return self

# defining the program to run this code: 

def prog(): 
    while True:
        try:
            x = int(input('Enter the number of particles in the system = ')) 
        
        except ValueError:
            print("Please enter a valid integer")
            continue 

        X = System(x).start()
        delta , lam = float(input('Delta = ')), float(input('Lambda = '))
        pot = input('Enter the potential for the system (either LJ or morse) = ') 
        break
        
    return X.opt(pot, lam, delta)


prog() 
