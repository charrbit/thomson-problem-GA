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

        # Generate a random chromosome for this individual
        self.generateNewChromosome()
        # Set the initial fitness value
        self.updateFitness()

    def generateNewChromosome(self):
        # Clear the current chromosome (if it exists)
        self.chromosome.clear()
        for i in range(self.numPoints):
            # Generate a new random point on the unit sphere
            newPoint = generateRandomPoint(rng)
            # Convert the point into its binary representation (BitArray)
            newPointBitArray = pointToBitArray(newPoint)
            # Add it to the chromosome
            self.chromosome.append(newPointBitArray)

    def updateFitness(self):
        ''' Calculates the total electrostatic potential for the configuration '''
        # Seperate each point within the chromosome BitArray
        points = splitBitArray(self.chromosome, self.numPoints)
        
        # Calculate the fitness score for each pair of points
        fitnessScore = 0
        for i in range(self.numPoints - 1):
            for j in range(i+1, self.numPoints):
                fitnessScore += electroStaticPotential(points[i], points[j])
        self.fitnessScore = fitnessScore

    def mutate(self):
        # Seperate the chromosome into a list of point BitArrays
        points = splitBitArray(self.chromosome, self.numPoints)
        # Pick a random point in the chromosome BitArray
        pointIndex = rng.integers(0, self.numPoints - 1)
        point = points[pointIndex]
        # Pick a random index in the point BitArray to mutate
        mutIndex = rng.integers(0, len(point.bin) - 1)
        point.invert(mutIndex)
        # Check if this mutation lead to a coordinate value outside the acceptable range of [0, 1]
        if not isValidPoint(point):
            # Generate a new valid point
            point = pointToBitArray(generateRandomPoint(rng))
        # Replace the old point in the list
        points[pointIndex] = point
        # Rebuild the chromosome
        self.chromosome = joinBitArray(points)
        # Update the fitness score for the new configuration
        self.updateFitness()