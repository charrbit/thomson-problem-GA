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
        numParents = len(self.individals)
        # While the population is not back to its original size
        while len(self.individals) < self.numIndividuals:
            # Randomly choose two parents to cross
            parent1index = rng.integers(0, numParents - 1)
            parent2index = rng.integers(0, numParents - 1)
            # Ensure the same parent is not chosen twice
            while parent1index == parent2index:
                parent2index = rng.integers(0, numParents - 1)

            # Split each parent chromosome into two halves for easier crossing
            parent1chromosome = splitBitArray(self.individals[parent1index].chromosome, 2)
            parent2chromosome = splitBitArray(self.individals[parent2index].chromosome, 2)

            # Generate a new individual in the population
            newChild = Individual(self.numPoints)
            # Randomly choose which parent the child gets the first half of its chromosome from
            # Assume it's parent1 at first
            newChild.chromosome = joinBitArrays([parent1chromosome[0], parent2chromosome[1]])
            if rng.random() <= 0.5:
                newChild.chromosome = joinBitArrays([parent2chromosome[0], parent1chromosome[1]])

            # Ensure the chromosome cross didn't create any invalid points
            points = splitBitArray(newChild.chromosome, self.numPoints)
            indicesToReplace = []
            for i, point in enumerate(points):
                if not isValidPoint(point):
                    indicesToReplace.append(i)
            # If an invalid point is found
            if len(indicesToReplace) > 0:
                for index in indicesToReplace:
                    # Replace it with a new random valid point
                    points[index] = pointToBitArray(generateRandomPoint(rng))
            # Recombine the chromosome and update the fitness
            newChild.chromosome = joinBitArrays(points)
            newChild.updateFitness()

            self.individals.append(newChild)

    def applyMutations(self):
        count = 0
        for i in range(self.numIndividuals):
            if rng.random() <= self.mutationRate:
                count += 1
                self.individals[i].mutate()
        return count