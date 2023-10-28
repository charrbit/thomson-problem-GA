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
        self.chromosome = joinBitArrays(points)
        # Update the fitness score for the new configuration
        self.updateFitness()

    def plot(self):
        # Generate grid points around the unit sphere
        theta = np.linspace(0, np.pi, 16)
        phi = np.linspace(0, 2*np.pi, 32)
        # Convert the spherical grid points to Cartesian
        x = np.outer(np.sin(theta), np.cos(phi))
        y = np.outer(np.sin(theta), np.sin(phi))
        z = np.outer(np.cos(theta), np.ones_like(phi))

        # Generate the 3D figure
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.set(xticklabels=[],
               yticklabels=[],
               zticklabels=[],
               title=f'Number of Points: {self.numPoints}')
        # Add the unit sphere grid to the figure
        ax.plot_wireframe(x, y, z, color='k', alpha=0.3)

        # Add the points of the individual onto the surface of the sphere
        points = splitBitArray(self.chromosome, self.numPoints)
        points = [bitArrayToPoint(point) for point in points]
        points = np.array(points) * np.array([np.pi, 2*np.pi])
        theta, phi = points[:, 0], points[:, 1]
        x, y, z = sphericalToCartesian(1, theta, phi)
        ax.scatter(x, y, z, s=20, color='b')

        # Show the figure
        plt.show()

    def print(self):
        print(f'Number of points: {self.numPoints}')
        # Get each point as the spherical theta and phi values
        points = splitBitArray(self.chromosome, self.numPoints)
        points = [bitArrayToPoint(point) for point in points]
        points = np.array(points) * np.array([np.pi, 2*np.pi])
        for i, point in enumerate(points):
            print(f'Point {i}:')
            print(f'\tTheta:{point[0]} rad')
            print(f'\tPhi:{point[1]} rad')