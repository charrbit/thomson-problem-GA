import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from helperFunctions import *

class Individual:
    ''' An individual of the genetic algorithm consists of a possible solution 
        to the Thomson problem for a given number of points on the unit sphere '''

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
            self.chromosome.append(generateRandomPoint())

    def updateFitness(self):
        ''' Calculates the total electrostatic potential for the configuration of points '''
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
    
    def print(self):
        print(f'Number of points: {self.numPoints}')
        print(f'Fitness score: {self.fitnessScore}')
        # Get each point as the real-valued thetaReal and phiReal coordinates
        points = [point.getPoint() for point in self.chromosome]
        for i, point in enumerate(points):
            print(f'Point {i}:')
            print(f'\tTheta: {point[0]} rad')
            print(f'\tPhi: {point[1]} rad')

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
                title=f'Number of Points: {self.numPoints}\nFitness Score: {self.fitnessScore}')
        # Add the unit sphere grid to the figure
        ax.plot_wireframe(x, y, z, color='k', alpha=0.3)

        # Add the points of the individual onto the surface of the sphere
        points = [point.getPointSpherical() for point in self.chromosome]
        theta = [point[0] for point in points]
        phi = [point[1] for point in points]
        x, y, z = sphericalToCartesian(1, theta, phi)
        ax.scatter(x, y, z, s=26, color='b')

        # Show the figure
        plt.show()