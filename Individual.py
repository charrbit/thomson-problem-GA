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
    # Then, a point on the unit sphere can be described with two real numbers, thetaReal and phiReal
        # P(theta, phi) = P(thetaReal * pi, phiReal * 2pi) => P(thetaReal, phiReal) where thetaReal and phiReal in [0, 1]
    # Let each real coordinate be represented by a 64-bit modifiable "chromosome" (BitArray)
    def __init__(self, numPoints):
        self.chromosome = ''
        self.numPoints = numPoints
        self.chromLength = 64

        # Generate a random chromosome for this individual
        for i in range(self.numPoints):
            thetaReal = rng.random()
            phiReal = rng.random()
            # Convert these values to a bitstring and append it to the chromosome
            self.chromosome += BitArray(float=thetaReal, length=self.chromLength).bin + BitArray(float=phiReal, length=self.chromLength).bin
        # Set the initial fitness value
        self.fitnessScore = self.getFitness()

    def fitness(self, point1, point2):
        ''' Calculates the electrostatic potential between two points where
            the points are given as bitstrings containing a thetaReal and phiReal value'''
        # Extract thetaReal and phiReal for each point in bitstring form
        P1 = np.array([point1[0 : self.chromLength], point1[self.chromLength:]])
        P2 = np.array([point2[0 : self.chromLength], point2[self.chromLength:]])
        # Convert the thetaReal and phiReal values from bitstring to a float represented as a string
        for i in range(len(P1)):
            P1[i] = BitArray(bin=P1[i]).float
            P2[i] = BitArray(bin=P2[i]).float
        # Convert the points to floats then find the difference between the thetaReal and phiReal values
        diff = P1.astype(float) - P2.astype(float)
        # Convert thetaReal and phiReal values to spherical coordinates theta and phi
        diff *= np.array([np.pi, 2*np.pi])
        return 1 / np.linalg.norm(diff)

    def getFitness(self):
        ''' Calculates the total electrostatic potential for the configuration '''
        # Seperate each point within the chromosome bitstring
        points = []
        for i in range(self.numPoints):
            points.append(self.chromosome[i*(self.chromLength * 2) : (i+1)*(self.chromLength * 2)])
        
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
            self.chromosome = Individual(self.numPoints).chromosome
        self.fitnessScore = self.getFitness()

    def isBadChromosome(self):
        ''' When mutating or crossing chromosomes, it is possible to get a value that is
            not within the appropriate [0,1] range for thetaReal or phiReal. If this happens, the 
            chromosome is considered bad and is replaced with a random valid one '''
        for i in range(self.numPoints):
            point = self.chromosome[i*(self.chromLength*2) : (i+1)*(self.chromLength*2)]
            P = np.array([point[0 : self.chromLength], point[self.chromLength:]])
            P[0] = BitArray(bin=P[0]).float
            P[1] = BitArray(bin=P[1]).float
            P = P.astype(float)
            if (P[0] < 0) or (P[0] > 1) or (P[1] < 0) or (P[1] > 1):
                return True
        return False