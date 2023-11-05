import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

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
        self.chromosome = []
        self.numPoints = numPoints
        self.fitnessScore = 0

        # Generate a random chromosome for this individual
        self.generateChromosome()
        # Set the initial fitness value
        self.updateFitness()

    def generateChromosome(self):
        ''' Generates a random list of points on the unit sphere '''
        # Reset the current chromosome
        self.chromosome.clear()
        for i in range(self.numPoints):
            # Generate a random point on the unit sphere and add it to the chromosome
            chromosome.append(generateRandomPoint())

    def updateFitness(self):
        ''' Calculates the total electrostatic potential for the configuration of points'''
        # Reset the current fitness score
        self.fitnessScore = 0
        # Calculate the fitness score for each pair of points
        for i in range(self.numPoints):
            for j in range(i + 1, self.numPoints):
                self.fitnessScore += electroStaticPotential(self.chromosome[i], self.chromosome[j])

    def mutate(self):
        ''' Randomly changes the value of one of the chromosome point coordinates '''
        # Pick a random point from the chromosome list
        pointIndex = rng.integers(self.numPoints)
        point = self.chromosome[pointIndex]
        # Pick a random coordinate from the point
        coordIndex = rng.integers(2)
        coord = point.getPointBitArray()[coordIndex]
        # Pick a random bit in the coordinate to mutate
        coord.invert(rng.integers(len(coord)))
        # Get the new mutated coordinate value
        mutVal = coord.float64

        # Update the point with the new coordinate
        newPointCoords = point.getPoint()
        newPointCoords[coordIndex] = mutVal
        point.setPoint(newPointCoords)
        # Check if this mutation lead to a coordinate value outside the acceptable range of [0, 1]
        if not point.isValid():
            point = generateRandomPoint()

        # Replace the old point in the list
        self.chromosome[pointIndex] = point
        # Update the fitness score for the new configuration
        self.updateFitness()