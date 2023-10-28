import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from Individual import *

class Population:
    def __init__(self, numPoints, numIndividuals, mutationProb):
        self.individals = []
        self.numPoints = numPoints
        self.numIndividuals = numIndividuals
        self.currGeneration = 1
        self.mutationRate = mutationProb

        # Generate the population
        for i in range(numIndividuals):
            self.individals.append(Individual(numPoints))

    def sortPopulation(self):
        # Sort individuals in ascending order of fitness score
        self.individals.sort(key=lambda x: x.fitnessScore)

    def makeSelection(self):
        self.sortPopulation()
        # Keep only the top 50% to be the parents
        self.individals = self.individals[:int(self.numIndividuals/2)]

    def crossParents(self):
        ''' Randomly pairs up two individuals within the population to cross. Whether
            the child inherits the beginning or end of a particular parent is randomly chosen '''
        # While the population is not back to its original size
        while len(self.individals) < self.numIndividuals:
            # Randomly choose two parents to cross
            parent1index = rng.integers(0, len(self.individals) - 1)
            parent2index = rng.integers(0, len(self.individals) - 1)
            # Ensure the same parent is not chosen twice
            while parent1index != parent2index:
                parent2index = rng.integers(0, len(self.individals))

            parent1chromosome = self.individals[parent1index].chromosome
            parent2chromosome = self.individals[parent2index].chromosome
            newChild = Individual(self.numPoints)
            # Randomly choose which parent the child gets the first half of its chromosome from
            if rng.random() <= 0.5:
                newChild.chromosome = parent1chromosome[:int(len(parent1chromosome)/2)] + parent2chromosome[int(len(parent2chromosome)/2):]
            else:
                newChild.chromosome = parent2chromosome[:int(len(parent2chromosome)/2)] + parent1chromosome[int(len(parent1chromosome)/2):]

            if newChild.isBadChromosome():
                # Generate a new valid chromosome
                newChild.chromosome = Individual(self.numPoints).chromosome
            newChild.fitnessScore = newChild.getFitness()
            self.individals.append(newChild)

    def applyMutations(self):
        count = 0
        for i in range(self.numIndividuals):
            if rng.random() <= self.mutationRate:
                count += 1
                self.individals[i].mutate()
        return count

    def getAverageScore(self):
        total = 0
        for i in range(self.numIndividuals):
            total += self.individals[i].fitnessScore
        return total / self.numIndividuals