import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from Individual import *

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
            r = k.astype(float) * np.array([np.pi, 2*np.pi])
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
            pos = k.astype(float) * np.array([np.pi, 2*np.pi])
            pointX = np.sin(pos[1])*np.cos(pos[0])
            pointY = np.sin(pos[1])*np.sin(pos[0])
            pointZ = np.cos(pos[1])
            ax.scatter(pointX, pointY, pointZ, s=15, color="r")
        plt.show()