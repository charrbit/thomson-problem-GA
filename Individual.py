from bitstring import BitArray
import random as rng
import numpy as np

class Individual:
    ''' An individual of the genetic algorithm consits of an approximate solution to the Thomson
        problem for a given number of electrons (i.e. a set of points on the sphere) '''
    # A point on a sphere, in a spherical coordinate system, can be described by three values
        # P = (R, theta, phi)
    # For the unit sphere, R = 1 => P = (1, theta, phi) => P(theta, phi)
        # Note: theta in [0, pi] and phi in [0, 2pi]
    # Then, a point on the unit sphere can be described with two real numbers, k and k'
        # P(theta, phi) = P(k * pi, k' * 2pi) where k and k' in [0, 1] => P(k, k')
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
        diff = k_i.astype(float) - k_j.astype(float)
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
            k = k.astype(float)
            if (k[0] < 0) or (k[0] > 1) or (k[1] < 0) or (k[1] > 1):
                return True
        return False