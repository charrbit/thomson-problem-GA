from argparse import ArgumentParser

from Population import *

# Parse arguments
parser = ArgumentParser(description='Genetic Algorithm Approach to the Thomson Problem')
parser.add_argument('-np', '--nPoints', type=int, default=2, help='Number of points to place on the sphere')
parser.add_argument('-ni', '--nIndividuals', type=int, default=50, help='Number of individuals in the population')
parser.add_argument('-ng', '--nGenerations', type=int, default=25, help='Number of generations of evolution to simulate')
parser.add_argument('-mp', '--mutationProb', type=float, default=0.1, help='Probability of an individual mutating')
args=parser.parse_args()

# Set arguments
numPoints = args.nPoints
numIndividuals = args.nIndividuals
numGenerations = args.nGenerations
mutationProb = args.mutationProb

# Display parsed arguments
print('Genetic Algorithm Approach to the Thomson Problem\n')
print('Number of points to place on the sphere: ', numPoints)
print('Number of individuals in the population: ', numIndividuals)
print('Number of generations of evolution to simulate: ', numGenerations)
print('Probability of an individual mutating: ', mutationProb)

# Start the algorithm
print('\nStarting . . . \n')
population = Population(numPoints, numIndividuals, mutationProb)
print('Initial population generated!')
print(f'Initial best fitness: {population.getBest().fitnessScore}\n')
while (population.currGeneration < numGenerations):
    print(f'*** Starting Generation: {population.currGeneration} ***')
    # Generate new children from the best individuals
    population.crossParents()
    # Mutate the population
    print(f'\tThere were {population.applyMutations()} mutations this generation!')
    # Print the best individual
    print(f'\tBest fitness: {population.getBest().fitnessScore}\n')
    population.currGeneration += 1
print(f'Final best fitness: {population.getBest().fitnessScore}')
population.getBest().plot()