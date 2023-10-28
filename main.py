from argparse import ArgumentParser

from Population import *

# Parse arguments
parser = ArgumentParser(description='Genetic Algorithm Approach to the Tompson Problem')
parser.add_argument('-np', '--nPoints', type=int, default=16, help='Number of points to place on the sphere')
parser.add_argument('-ni', '--nIndividuals', type=int, default=50, help='Number of individuals in the population')
parser.add_argument('-ng', '--nGenerations', type=int, default=100, help='Number of generations of evolution to simulate')
parser.add_argument('-mp', '--mutationProb', type=float, default=0.1, help='Probability of an individual mutating')
parser.add_argument('-t', '--tol', type=float, default=1e-12, help='Tolerence between successive generation average scores before accepting convergence')
args=parser.parse_args()

# Set arguments
numPoints = args.nPoints
numIndividuals = args.nIndividuals
numGenerations = args.nGenerations
mutationProb = args.mutationProb
tol = args.tol

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
population.printBestIndividual()
print('\n')
prevAvgScore = 0
while (population.currGeneration < numGenerations) and (abs(population.getAverageScore()-prevAvgScore) > tol):
    print('*** Starting Generation: ', population.currGeneration + 1, '***\n')    
    prevAvgScore = population.getAverageScore()
    population.makeSelection()
    population.crossParents()
    mutCount = population.applyMutations()
    print('There were ', mutCount, 'mutations this generation!\n')
    population.printBestIndividual()
    population.currGeneration += 1
population.plotBestIndividual()
    