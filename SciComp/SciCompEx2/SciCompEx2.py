import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate
from matplotlib import cm
import re
import numpy as np

# firstly we need to import and read our files
f = open("Python/SciCompEx2/H2O/H2O.r0.70theta70.0.out","r")
for line in f:
   if "Done" in line: 
      l = line.split()
      # this is a check of the energy of the conformation.
      energy = l[-5]
      print(energy)
# this prints the energy value from the file therefore this function works

def prog(): 
   # we first need to define our dictionary: 
   dic = {"Energy":[],"Angle":[],"Length":[]}
   # the bond lengths jump in 0.05 angstroms, the bond angles jump in 1 degree increments, we can set these 2 different ranges to iterate over opening the files and recording the energies:   
   for i in np.arange(l_min,l_max,0.05):
      for j in np.arange(70,161,1):
         if np.round(i*10, 2) == int(np.round(i*10, 2)):
            f = open("/Users/mattball/Desktop/Python/SciCompEx2/"+str(mol)+"/"+str(mol)+".r"+str(np.round(i,2))+"0theta"+str(np.round(j,2))+".0.out", "r")
            for line in f:
               if "Done" in line: 
                  l = line.split()
                  e = l[-5]
                  dic["Energy"].append(e)
                  dic["Angle"].append(j)
                  dic["Length"].append(i)
               else: 
                  pass
         else: 
            f = open("/Users/mattball/Desktop/Python/SciCompEx2/"+str(mol)+"/"+str(mol)+".r"+str(np.round(i,2))+"theta"+str(np.round(j,2))+".0.out", "r")
            for line in f:
               if "Done" in line: 
                  l = line.split()
                  e = l[-5]
                  dic["Energy"].append(e)
                  dic["Angle"].append(j)
                  dic["Length"].append(i)
               else: 
                  pass

   # we have extracted the minimum energy, we now need to find the values of the bond length and bond angle that this corresponds to. 
   min_energy = max(dic["Energy"])
   print("minimum energy =", min_energy) 

   # finding and printing the optimum geometry
   prog.bond_length = dic["Length"][dic["Energy"].index(max(dic["Energy"]))]
   prog.bond_angle = dic["Angle"][dic["Energy"].index(max(dic["Energy"]))]
   optimum_geom = (str(float(prog.bond_length)) + str(" Angstrom"), str(prog.bond_angle) + str(" Degrees"))
   print("The optimum geometry is : ", optimum_geom)
               

   # we can now print our contour plot:
   # converting our lists into array so they can pass through the function
   x = np.array(dic["Angle"])
   y = np.array(dic["Length"])
   z = np.array(dic["Energy"])

   #gridding the data (and relating the energies to each individual point)
   xgrid = np.array(dic["Angle"])
   ygrid = np.array(dic["Length"])
   xgrid, ygrid = np.meshgrid(xgrid, ygrid)
   zgrid = np.array(interpolate.griddata((x,y),z, (xgrid, ygrid)))

   # plotting the contour graph
   fig, ax = plt.subplots(subplot_kw={"projection": "3d"}) 
   surf = ax.plot_surface(xgrid, ygrid, zgrid, cmap=cm.coolwarm,linewidth=0.5, edgecolor ='black', ccount = 6)
   plt.xlabel('Bond Angle / Degrees')
   plt.ylabel('Bond Length / Angstrom')
   ax.set_zlabel('Energy / Hartree', rotation=180)
   plt.title('Potential Energy Surface of '+str(mol))
   fig.colorbar(surf, shrink=0.5, aspect=10)
   plt.show()

# finally we can determine the vibrational frequencies for the symmetric stretch and the bending mode by making some assumptions: 
# we can fit a 2d polynomial to the data around the minimum point to extract the force constants for a change in bond length k_r and bond angle k_theta.
# we can then extract the vibrational frequencies from the parameters we find. 
def vib_freq():
   prog() 
   E_eq = []
   BA_eq = []
   BL_eq = []
   dic_eq = {"Energy":[],"Angle":[],"Length":[]}

   # the bond lengths jump in 0.05 angstroms, the bond angles jump in 1 degree increments, we can set these 2 different ranges to iterate over opening the files and recording the energies:   
   for i in np.arange(prog.bond_length - 0.05 ,prog.bond_length +0.10,0.05):
      for j in np.arange(int(prog.bond_angle) -5,int(prog.bond_angle) +6,1):
         if np.round(i*10, 2) == int(np.round(i*10, 2)):
            f = open("/Users/mattball/Desktop/Python/SciCompEx2/"+str(mol)+"/"+str(mol)+".r"+str(np.round(i,2))+"0theta"+str(np.round(j,2))+".0.out", "r")
            for line in f:
               if "Done" in line: 
                  l = line.split()
                  e = float(l[-5])
                  dic_eq["Energy"].append(e)
                  dic_eq["Angle"].append(j*np.pi /180)
                  dic_eq["Length"].append(i)
               else: 
                  pass
         else: 
            f = open("/Users/mattball/Desktop/Python/SciCompEx2/"+str(mol)+"/"+str(mol)+".r"+str(np.round(i,2))+"theta"+str(np.round(j,2))+".0.out", "r")
            for line in f:
               if "Done" in line: 
                  l = line.split()
                  e = float(l[-5])
                  dic_eq["Energy"].append(e)
                  dic_eq["Angle"].append(j*np.pi /180)
                  dic_eq["Length"].append(i)
               else: 
                  pass

   x = np.array(dic_eq["Angle"]) - prog.bond_angle * np.pi /180
   y = np.array(dic_eq["Length"]) - prog.bond_length
   z = np.array(dic_eq["Energy"]) 

   #gridding the data (and relating the energies to each individual point)
   xgrid = np.array(dic_eq["Angle"]) - prog.bond_angle * np.pi /180
   ygrid = np.array(dic_eq["Length"]) - prog.bond_length
   xgrid, ygrid = np.meshgrid(xgrid, ygrid)
   zgrid = np.array(interpolate.griddata((x,y),z, (xgrid, ygrid)))

   X = xgrid.flatten()
   Y = ygrid.flatten()
   A = np.array([X*0+1, 1/2 * X**2, 1/2 *Y**2]).T
   B = zgrid.flatten()

   # fitting the function to the quadratic
   coeff, r, rank, s = np.linalg.lstsq(A, B, rcond=-1)

   # converting the coefficients to SI units : 
   coeff_ = [[0],[0],[0]]
   coeff_[0]  = coeff[0] 
   coeff_[1] = coeff[1] * 4.35974e-18
   coeff_[2] = coeff[2] * 4.35974e2 
   print(coeff_)

   nu_1 = (1/(2*np.pi) * np.sqrt(coeff_[2]/(2*1.66e-27)))/2.998e10
   nu_2 = (1/(2*np.pi) * np.sqrt(coeff_[1]/((0.5*1.66e-27)*((prog.bond_length*1e-10)**2)))) / 2.998e10
   print(nu_1, nu_2)

# Creating the program:
mol = input('Please type the molecule you are analysing: ') 
if mol == 'H2S': 
   l_max = 1.85
   l_min =  0.60
   vib_freq()

elif mol == 'H2O': 
   l_max = 1.95
   l_min = 0.70
   vib_freq()

else: 
   print('Not a valid input, type H2S or H2O ')
