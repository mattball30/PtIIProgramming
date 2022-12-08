# We will attempt to create a program which finds Huckel pi energies and degeneracies for different scenarios.
# We know the energy of overlap of orbitals on the same carbon is a, and on adjacent carbons is b
# We also know that the eigenvalues of the huckel matrix will give the energies, whilst the degeneracies of the eigenvalues will give the degeneracies of the energies.

import numpy as np

# Hard coding a specific matrix and finding its eigenvalues (using the same example as the notes):
Mz = np.zeros((4,4))
lc = Mz.shape[0]
for i in range(lc):
   Mz[i, lc-i-1] = 1 
MI = np.identity(4)
M2 = Mz + MI
evals, evecs = np.linalg.eig(M2)
print(sorted(evals))

# Creating a function that will give the eigenvalues of a matrix:

def get_evals(x_): 
  evals, evecs = np.linalg.eig(x_)
# converting our function output into a list (previously a NoneType object):
  evals_list = list(sorted(evals))
#also noting due to the non exact nature of finding the eigenvalues the energies all appear to be non degenerate - we can round these eigenvalues
  evround = []
  for i in range(len(evals_list)):
   evround = np.append(evround,round(evals_list[i],3))
  return list(evround)

#defining a counter, using the dictionary in python, if the element is already in k, add one, otherwise add a new element and give it a value of 1:
def degeneracy(a):
  k = {}
  for m in a:
    if m in k:
      k[m] +=1
    else:
      k[m] =1
  return k

# We can check this function works by checking the eigenvalues of a set matrix:
get_evals(M2)

# This matches the eigenvalues calculated so we can assume the function works (note the values match without rounding - that has been added for the degeneracies step later)!

def huck_2d():
  chain = input("Please input polyene type:")
  if chain == "linear":
     print("How many carbons in the", chain, "polyene chain?:")
     n = int(input())
     Hz = np.zeros((n,n))
     print("alpha value:")
     a = int(input())
     print("beta value:")
     b = int(input())
     # Now writing a function to create a general Huckel matrix of n dimensions. 
     for i in range(n):
        for j in range(n): 
          if j == i:
            Hz[i,j] = a
          elif abs(i-j) == 1:
            Hz[i,j] = b
          else: Hz[i,j] = 0
     # checking the huckel matrix is  correct (checking with small n)
     if n <=5:
       print(Hz)
     # solving our generalised huckel matrix:
     print("Eigenvalues for n =", n, ":")
     print(get_evals(Hz))
     # adding in our degeneracies:
     print("Degeneracy =")
     print(degeneracy(get_evals(Hz)))

  elif chain == "cyclic":
      print("How many carbons in the", chain, "polyene chain?:")
      n = int(input())
      Hz = np.zeros((n,n))
      print("alpha value:")
      a = int(input())
      print("beta value:")
      b = int(input())
      # the code for the cylic polyene is very similar, however we now have an interaction between our 'ends' of the equivalent linear polyene: 
      for i in range(n):
        for j in range(n): 
          if j == i:
            Hz[i,j] = a
          elif abs(i-j) == 1:
            Hz[i,j] = b
      # note here the extra condition to get the end contribution
          elif abs(i-j) == n-1:  
            Hz[i,j] = b
          else: Hz[i,j] = 0

      # checking the huckel matrix is  correct (this can be deleted before the n = 100 test)
      if n <=6:
       print(Hz)
    
      # solving our generalised huckel matrix:
      print("Eigenvalues for n =", n, ":")
      print(get_evals(Hz))
      # adding in our degeneracies:
      print("Degeneracy =")
      print(degeneracy(get_evals(Hz)))

  else: 
      print("not a valid input, please choose cyclic or linear")
      
# now defining the huckel matrix for the 3D platonic solids: 
def huck_3d(): 
 x = input("Please input 3D platonic solid type: ")
 if x == "tetrahedron":
     Hz = np.zeros((4,4))
     print("alpha value:")
     a = int(input())
     print("beta value:")
     b = int(input())
     # Now writing a function to create a huckel matrix for a tetrahedron. We note all atoms are directly connected to all other atoms.  
     for i in range(4):
        for j in range(4): 
          if j == i:
            Hz[i,j] = a
          else: Hz[i,j] = b
     # checking the huckel matrix is  correct (checking with small n)
     print(Hz)
     # solving our generalised huckel matrix:
     print("Eigenvalues for", x,":")
     print(get_evals(Hz))
     # adding in our degeneracies:
     print("Degeneracy =")
     print(degeneracy(get_evals(Hz)))

 elif x == "cube":
      print("alpha value:")
      a = int(input())
      print("beta value:")
      b = int(input())
      # the matrix for the cube is much harder to program - so we manually hard code it: 
      Hz = np.array(([a,b,0,b,b,0,0,0],[b,a,b,0,0,b,0,0],[0,b,a,b,0,0,b,0],[b,0,b,a,0,0,0,b],[b,0,0,0,a,b,0,b],[0,b,0,0,b,a,b,0],[0,0,b,0,0,b,a,b],[0,0,0,b,b,0,b,a]))
      # checking the huckel matrix is  correct
      print(Hz)
      # solving our generalised huckel matrix:
      print("Eigenvalues for", x, ":")
      print(get_evals(Hz))
      # adding in our degeneracies:
      print("Degeneracy =")
      print(degeneracy(get_evals(Hz)))

 elif x == "dodecahedron":
      print("alpha value:")
      a = int(input())
      print("beta value:")
      b = int(input())
      # the matrix for the dodecahedron can be programmed as: 
      Hz = np.zeros((20,20))
      for i in range(20):
       for j in range(20): 
        if i == j:
         Hz[i,j] = a
        elif abs(i-j) == 1: 
         Hz[i,j] = b 
        else: 
         Hz[i,j] = 0 
      for [i,j] in {(0,5),(0,8),(1,9),(2,11),(3,13),(5,14),(6,16),(8,17),(10,18),(12,19),(15,19)}:
       Hz[i,j] = Hz[j,i] = b
      # solving our generalised huckel matrix:
      print("Eigenvalues for", x, ":")
      print(get_evals(Hz))
      # adding in our degeneracies:
      print("Degeneracy =")
      print(degeneracy(get_evals(Hz)))

 else: 
      print("not a valid input, please choose tetrahedron, cube or dodecahedron.")

# We can now combine the two into a single program:
def prog():
 y = input("Is the pi system 2d or 3d? : ")
 if y == "2d":
  huck_2d()
 elif y == "3d": 
  huck_3d()
 else: 
   print("Input must be 2d or 3d, please try again.")
   
prog()
