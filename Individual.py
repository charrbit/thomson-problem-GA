import numpy as np

from helperMethods import *

# Create random number generator
rng = np.random.default_rng()

class Individual:
    ''' An individual of the genetic algorithm consits of an approximate solution to the Thomson
        problem for a given number of electrons (i.e. a set of points on the sphere) '''

    # A point on a sphere, in a spherical coordinate system, can be described by three values
        # P = (R, theta, phi)
    # For the unit sphere, R = 1 => P = (1, theta, phi) => P(theta, phi)
        # Note: theta in [0, pi] and phi in [0, 2pi]
    # Then, a point on the unit sphere can be described with two real numbers, thetaReal and phiReal
        # P(theta, phi) = P(thetaReal * pi, phiReal * 2pi) => P(thetaReal, phiReal) where thetaReal and phiReal in [0, 1]

    # Let each real coordinate be represented by a 64-bit BitArray
    # Let each point be represented by a 128-bit BitArray
    # Let the chromosome be the combination of all points in the configuration be represented by a (128*numPoints)-bit BitArray

    def __init__(self, numPoints):
        self.chromosome = BitArray(bin='0b') # Empty
        self.numPoints = numPoints
        self.chromLength = 64

        # Generate a random chromosome for this individual
        for i in range(self.numPoints):
            # Generate a new random point on the unit sphere
            newPoint = generateRandomPoint(rng)
            # Convert the point into its binary representation (BitArray)
            newPointBitArray = pointToBitArray(newPoint)
            # Add it to the chromosome
            self.chromosome.append(newPointBitArray)
        # Set the initial fitness value
        self.fitnessScore = self.getFitness()

    def fitness(self, point1, point2):
        ''' Calculates the electrostatic potential between two points where
            the points are given as BitArrays each containing a thetaReal and phiReal value'''
        # Extract thetaReal and phiReal values for both points
        P1 = bitArrayToPoint(point1)
        P2 = bitArrayToPoint(point2)
        # Convert the points to numpy arrays to find the vector difference
        diff = np.subtract(np.array(P1), np.array(P2))
        # Convert thetaReal and phiReal values to spherical coordinates theta and phi
        diff *= np.array([np.pi, 2*np.pi])
        return 1 / np.linalg.norm(diff)

    def getFitness(self):
        ''' Calculates the total electrostatic potential for the configuration '''
        # Seperate each point within the chromosome BitArray
        points = splitBitArray(self.chromosome, self.numPoints)
        
        # Calculate the fitness score for each combination of points
        fitnessScore = 0
        for i in range(self.numPoints - 1):
            for j in range(i+1, self.numPoints):
                fitnessScore += self.fitness(points[i], points[j])
        return fitnessScore

    def mutate(self):
        # Pick a random index in the chromosome BitArray to mutate
        mutIndex = rng.integers(0, len(self.chromosome) - 1)
        # Flip the bit
        self.chromosome.invert(mutIndex)
        # Check if this mutation lead to a coordinate value outside the acceptable range of [0, 1]
        if self.isBadChromosome():
            # Generate a new valid chromosome
            self.chromosome = Individual(self.numPoints).chromosome
        self.fitnessScore = self.getFitness()

    def isBadChromosome(self):
        ''' When mutating or crossing chromosomes, it is possible to get a point value that is
            not within the appropriate [0, 1] range for thetaReal or phiReal. If this happens, the 
            chromosome is considered bad and should be replaced with a valid one '''
        points = splitBitArray(self.chromosome, self.numPoints)
        for point in points:
            thetaReal, phiReal = bitArrayToPoint(point) 
            if (0 > thetaReal > 1) or (0 > phiReal > 1):
                return True
        return False