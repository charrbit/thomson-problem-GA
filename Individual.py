from Point import *

# Create random number generator
rng = np.random.default_rng()

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

# Helper Methods
def generateRandomPoint():
    ''' Generates a random point on the unit sphere '''
    thetaReal, phiReal = rng.random(2)
    return Point(thetaReal, phiReal)

def electroStaticPotential(point1, point2):
    ''' Calculates the electrostatic potential between two points  '''
    # Convert the points to spherical coordinates then find the distance between them
    dist = np.subtract(np.array(point1.getPointSpherical()), np.array(point2.getPointSpherical()))
    return 1 / np.linalg.norm(dist)