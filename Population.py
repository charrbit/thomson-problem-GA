from Individual import *

class Population:
    ''' A population of the genetic algorithm consists of a set of possible solutions to the Thomson problem
        where point coordinates are "genetically" crossed and mutated towards a fitness minimum '''
    
    def __init__(self, numPoints, numIndividuals, mutationProb):
        self.individals = []
        self.numPoints = numPoints
        self.numIndividuals = numIndividuals
        self.currGeneration = 1
        self.mutationRate = mutationProb

        # Generate the population of individuals
        for i in range(numIndividuals):
            self.individals.append(Individual(self.numPoints))

    def sortPopulation(self):
        ''' Sorts individuals in ascending order of fitness score (the
            Thomson problem seeks the minimum electrostatic potential) '''
        self.individals.sort(key=lambda x: x.fitnessScore)

    def makeSelection(self):
        ''' Splits individuals in half keeping only those with the lowest fitness scores '''
        self.sortPopulation()
        self.individals = self.individals[:int(self.numIndividuals/2)]
    
    def applyMutations(self):
        ''' Applies mutations to individuals based on the mutationRate
            and returns the number of mutations that occured '''
        count = 0
        for i in range(self.numIndividuals):
            if rng.random() <= self.mutationRate:
                count += 1
                self.individals[i].mutate()
        return count

    def crossParents(self):
        ''' Randomly pairs up two parent individuals within the population to cross and produce a
            child individual until the population size has returned to numIndividuals. A cross consists
            of randomly choosing points from both parents until the child has reached numPoints '''
        self.makeSelection()
        numParents = len(self.individals)
        # While the population is not back to its original size
        while len(self.individals) < self.numIndividuals:
            # Randomly choose two parents to cross
            parent1index, parent2index = rng.choice(numParents, size=2, replace=False)
            parent1 = self.individals[parent1index]
            parent2 = self.individals[parent2index]
            newChild = Individual(self.numPoints)
            # Get a list of all points contained in both parents
            parentPoints = parent1.chromosome + parent2.chromosome
            # Randomly select distinct points for the child to inherit
            for i in range(self.numPoints):
                pointIndex = rng.integers(len(parentPoints))
                newChild.chromosome[i].setPoint(parentPoints.pop().getPoint())
            # Add the child to the population
            newChild.updateFitness()
            self.individals.append(newChild)