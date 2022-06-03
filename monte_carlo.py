# This code runs a simulated annealing version of Monte-Carlo simulation for 13 particles.
# The simulation finds the minimum of a cluster of the 13 particles 
# interacting via isotropic pair potential V(x), using Lennard Jones parameters.

# IMRORTANT NOTE: running the code may result in overflow issues (it should still work).
# Specifically the exponent part for the simulated annealing: the exponent might by too big.
# Also, after about 50-60 iterations, the guesses stop changing. not sure why that happens.

import random
import numpy as np

N_particles = 13
T = 298.15 # Temperature in Kelvin

# define bounds of r as well as the "r" array
# coord_max defines the boundary for this problem. no coordinate can exceed 5 units.
coord_min, coord_max = 0, 5

coord = [] # initialize the random X, Y, and Z coordinates. Should be 3N coordinates.
for i in range(0,N_particles*3):
   n = coord_max*random.random()
   coord.append(n)

# Creating the r array (78 values)
# r = r1_2 r1_3 r1_4 ... r12_13

r = []
K = 1
for J in range(0, N_particles):
   if J < K:
      continue
   for i in range(K, N_particles):
      g = ((coord[3*K] - coord[3*J])**2 + (coord[(3*K)+1] - coord[(3*J)+1])**2 + (coord[(3*K)+2] - coord[(3*J)+2])**2)**(1/2) 
      r.append(g)
      K += 1
   K = J + 1

# ^this part^ doesn't work, not sure why. it generates the right amount of values
# (78) but it consistently generates 12 zeroes, which is not statistically probable


# Setting up the isotropic pair potential equation
# with Lennard Jones parameters:

def V(x):
   epsilon = 1
   sigma = 1
   return 4*epsilon*((sigma/x)**12 - (sigma/x)**6)

def sim(N, step_size, constraint):
# N = number of iterations 
# step_size is how far we want to search away from the guess
# constraint is the limit that the simulated annealing can guess without turning the current guess into an extremely large number
 
   best_eval = 0
   best = []
   counter = 0
   for i in range(0, len(r)):
      if r[i] == 0.0: # this part accounts for the zeroes by skipping them. not ideal, but i'll fix later
         continue
      best_eval += V(r[i]) 
      best.append(r[i])
      counter += 1
   print("Initial Guess = ", abs(best_eval))
# first, we define the "best" guess by calculating the isotropic pair potential using the r array
 
   curr, curr_eval = best, best_eval 
# This sets the working current guess. 

   for i in range(1, N+1):
      W = []
      for J in range(0, counter):
         e = step_size*random.uniform(-1, 1)
         W.append(e) 
# W array is the guess array, with values between the negative step size and the positive step size   

      candidate = curr + W
      candidate_eval = 0
      for j in range(0, counter):
         candidate_eval += V(candidate[i])
# W is added to the current guess array and the new guess is evaluated using the isotropic pair potential equation
 
      if abs(candidate_eval) < abs(best_eval):
         best, best_eval = candidate, candidate_eval 
# if the new guess array is smaller than the best guess, then the "best" guess is replaced by the current guess

      diff = abs(candidate_eval) - abs(curr_eval)
      if diff < 0: 
         curr, curr_eval = candidate, candidate_eval
# This part is Monte-Carlo. if the new guess is better than the "current" guess, the current guess is updated to the new guess

      if diff > 0:
         x = random.random()
         p = np.exp((abs(candidate_eval) - abs(curr_eval))/T)
         if p > x and candidate_eval <= constraint:
            curr, curr_eval = candidate, candidate_eval
# This part is simulated annealing. If the new guess is worse than the "current" guess, then by chance, the current guess might
# still by updated to the new guess anyway.

      print("Iteration", i, ", Current guess = ", abs(curr_eval))
         
   return [best,best_eval]

[x,y] = sim(100,0.2,1.0) # runs the simulation for a 100 iterations, using a 0.2 step size  
print("After 100 iterations, the best guess is ", abs(y))
