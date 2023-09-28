from bitstring import BitArray
import random as rng
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Individual:
    # r = ( theta, phi ) = ( k * pi, k' * 2pi ) = ( k, k' ) where k and k' in [0,1]
    def __init__(self, numPoints, SCBL):
        self.chromosome = ''
        self.numPoints = numPoints
        self.SCBL = SCBL # Single Coordinate Bit Length
        self.SPBL = self.SCBL * 2 # Single Point Bit Length

        # Generate a random chromosome for this individual
        for i in range(self.numPoints):
            k = rng.random()
            kprime = rng.random()
            # Convert these values to a bitstring and append it to the chromosome
            self.chromosome += BitArray(float=k, length=self.SCBL).bin + BitArray(float=kprime, length=self.SCBL).bin
        # Set the initial fitness value
        self.fitnessScore = self.getFitness()

    def fitness(self, point1, point2):
        ''' Calculates the electrostatic potential between two points where
            the points are given as bitstrings containing a k and k' value'''
        # Extract k and kprime for each point i, j in bitstring form
        k_i = np.array([point1[0 : self.SCBL], point1[self.SCBL:]])
        k_j = np.array([point2[0 : self.SCBL], point2[self.SCBL:]])
        # Convert the k and kprime values from bitstring to a float represented as a string
        for i in range(len(k_i)):
            k_i[i] = BitArray(bin=k_i[i]).float
            k_j[i] = BitArray(bin=k_j[i]).float
        # Convert the points to floats then find the difference between the k and kprime values
        diff = k_i.astype(np.float) - k_j.astype(np.float)
        # Convert k and kprime values to spherical coordinates theta and psi)
        diff *= np.array([np.pi, 2*np.pi])
        return 1 / np.linalg.norm(diff)

    def getFitness(self):
        ''' Calculates the total electrostatic potential for the configuration '''
        # Seperate each point within the chromosome bitstring
        points = []
        for i in range(self.numPoints):
            points.append(self.chromosome[i*self.SPBL : (i+1)*self.SPBL])
        
        # Calculate the fitness score for each combination of points
        fitnessScore = 0
        for i in range(self.numPoints - 1):
            for j in range(i+1, self.numPoints):
                fitnessScore += self.fitness(points[i], points[j])
        return fitnessScore

    def mutate(self):
        mutIndex = rng.randint(0, len(self.chromosome) - 1)
        # Get chromosome as a list so it can be modified
        # (Python strings are immutable)
        chromList = list(self.chromosome)
        if self.chromosome[mutIndex] == '1':
            chromList[mutIndex] = '0'
        else:
            chromList[mutIndex] = '1'
        self.chromosome = ''.join(chromList)
        if self.isBadChromosome():
            # Generate a new valid chromosome
            self.chromosome = Individual(self.numPoints, self.SCBL).chromosome
        self.fitnessScore = self.getFitness()

    def isBadChromosome(self):
        ''' When mutating or crossing chromosomes, it is possible to get a value that is
            not within the appropriate [0,1] range for k or kprime. If this happens, the 
            chromosome is considered bad and is replaced with a random valid one '''
        for i in range(self.numPoints):
            point = self.chromosome[i*(self.SCBL*2) : (i+1)*(self.SCBL*2)]
            k = np.array([point[0 : self.SCBL], point[self.SCBL:]])
            k[0] = BitArray(bin=k[0]).float
            k[1] = BitArray(bin=k[1]).float
            k = k.astype(np.float)
            if (k[0] < 0) or (k[0] > 1) or (k[1] < 0) or (k[1] > 1):
                return True
        return False

class Population:
    def __init__(self, numPoints, numIndividuals, mutationProb, SCBL):
        self.individals = []
        self.numPoints = numPoints
        self.numIndividuals = numIndividuals
        self.currGeneration = 0
        self.mutationRate = mutationProb
        self.SCBL = SCBL

        # Generate the population
        for i in range(numIndividuals):
            self.individals.append(Individual(numPoints, SCBL))

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
            parent1index = rng.randint(0, len(self.individals) - 1)
            parent2index = rng.randint(0, len(self.individals) - 1)
            # Ensure the same parent is not chosen twice
            while parent1index != parent2index:
                parent2index = rng.randint(0, len(self.individals))

            parent1chromosome = self.individals[parent1index].chromosome
            parent2chromosome = self.individals[parent2index].chromosome
            newChild = Individual(self.numPoints, self.SCBL)
            # Randomly choose which parent the child gets the first half of its chromosome from
            if rng.random() <= 0.5:
                newChild.chromosome = parent1chromosome[:int(len(parent1chromosome)/2)] + parent2chromosome[int(len(parent2chromosome)/2):]
            else:
                newChild.chromosome = parent2chromosome[:int(len(parent2chromosome)/2)] + parent1chromosome[int(len(parent1chromosome)/2):]

            if newChild.isBadChromosome():
                # Generate a new valid chromosome
                newChild.chromosome = Individual(self.numPoints, self.SCBL).chromosome
            newChild.fitnessScore = newChild.getFitness()
            self.individals.append(newChild)

    def applyMutations(self):
        count = 0
        for i in range(self.numIndividuals):
            if rng.random() <= self.mutationRate:
                count += 1
                self.individals[i].mutate()
        return count

    def printBestIndividual(self):
        self.sortPopulation()
        bestIndividual = self.individals[0]
        print('The best individual has a fitness score of: ', bestIndividual.fitnessScore)
        # Seperate each point within the chromosome bitstring
        for i in range(self.numPoints):
            point = bestIndividual.chromosome[i*(self.SCBL*2) : (i+1)*(self.SCBL*2)]
            k = np.array([point[0 : self.SCBL], point[self.SCBL:]])
            k[0] = BitArray(bin=k[0]).float
            k[1] = BitArray(bin=k[1]).float
            r = k.astype(np.float) * np.array([np.pi, 2*np.pi])
            print('Point', i, '| Theta: ', r[0], 'rad\n\t  Psi: ', r[1], 'rad\n')
        print('-------------------------------------------------------')

    def getAverageScore(self):
        total = 0
        for i in range(self.numIndividuals):
            total += self.individals[i].fitnessScore
        return total / self.numIndividuals

    def plotBestIndividual(self):
        # Generate the unit sphere
        r = 1
        phi, theta = np.mgrid[0.0:np.pi:200j, 0.0:2.0*np.pi:200j]
        x = r*np.sin(phi)*np.cos(theta)
        y = r*np.sin(phi)*np.sin(theta)
        z = r*np.cos(phi)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x,y,z, color='white', alpha=0.6, linewidth=0)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
        ax.set_title('Number of Points = ' + str(self.numPoints))

        # Plot the points associated with the best individual on the sphere
        self.sortPopulation()
        bestIndividual = self.individals[0]
        for i in range(self.numPoints):
            point = bestIndividual.chromosome[i*(self.SCBL*2) : (i+1)*(self.SCBL*2)]
            k = np.array([point[0 : self.SCBL], point[self.SCBL:]])
            k[0] = BitArray(bin=k[0]).float
            k[1] = BitArray(bin=k[1]).float
            pos = k.astype(np.float) * np.array([np.pi, 2*np.pi])
            pointX = np.sin(pos[1])*np.cos(pos[0])
            pointY = np.sin(pos[1])*np.sin(pos[0])
            pointZ = np.cos(pos[1])
            ax.scatter(pointX, pointY, pointZ, s=15, color="r")
        plt.show()