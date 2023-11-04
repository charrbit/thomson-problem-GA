from argparse import ArgumentParser

from Population import *

# Parse arguments
parser = ArgumentParser(description='Genetic Algorithm Approach to the Tompson Problem')
parser.add_argument('-np', '--nPoints', type=int, default=2, help='Number of points to place on the sphere')
parser.add_argument('-ni', '--nIndividuals', type=int, default=50, help='Number of individuals in the population')
parser.add_argument('-ng', '--nGenerations', type=int, default=25, help='Number of generations of evolution to simulate')
parser.add_argument('-mp', '--mutationProb', type=float, default=0.1, help='Probability of an individual mutating')
# parser.add_argument('-t', '--tol', type=float, default=1e-12, help='Tolerence between successive generation average scores before accepting convergence')
args=parser.parse_args()

# Set arguments
numPoints = args.nPoints
numIndividuals = args.nIndividuals
numGenerations = args.nGenerations
mutationProb = args.mutationProb
# tol = args.tol

# Display parsed arguments
print('Genetic Algorithm Approach to the Thomson Problem\n')
print('Number of points to place on the sphere: ', numPoints)
print('Number of individuals in the population: ', numIndividuals)
print('Number of generations of evolution to simulate: ', numGenerations)
print('Probability of an individual mutating: ', mutationProb)

# Start the algorithm
print('\nStarting . . . \n')
population = Population(numPoints, numIndividuals, mutationProb)
print('Initial population generated!\n')
while (population.currGeneration < numGenerations):
    print(f'*** Starting Generation: {population.currGeneration} ***\n')    
    population.makeSelection()
    population.crossParents()
    mutCount = population.applyMutations()
    print(f'There were {mutCount} mutations this generation!\n')
    population.currGeneration += 1