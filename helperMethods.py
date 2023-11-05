import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

# Convert a float into a BitArray (Should be 64-bit)
def floatToBitArray(floatToConvert):
    return BitArray(float=floatToConvert, length=64)

# Convert a BitArray into a float
def bitArrayToFloat(bitArrayToConvert):
    return bitArrayToConvert.float64

# Convert a point containing two float coordinates into a single BitArray
def pointToBitArray(pointToConvert):
    # Convert the thetaReal and phiReal coordinates into BitArrays
    point = [floatToBitArray(x) for x in pointToConvert]
    # Combine the coordinates into a single BitArray
    point = joinBitArrays(point)
    return point

# Convert a single BitArray point into two float coordinates
def bitArrayToPoint(bitArrayToConvert):
    # Split the BitArray point into its two component BitArrays
    point = splitBitArray(bitArrayToConvert, 2)
    # Convert each component into a float
    point = [bitArrayToFloat(x) for x in point]
    return point

# Split a BitArray into a list of numSegment BitArrays (Should have len(BitArray) % numSegments = 0)
def splitBitArray(bitArrayToSplit, numSegments):
    bitArray = bitArrayToSplit.bin
    segmentLength = int(len(bitArray) / numSegments)
    bitArrayList = [bitArray[segmentLength*i : segmentLength*(i+1)] for i in range(numSegments)]
    bitArrayList = [BitArray(bin=f'0b{item}') for item in bitArrayList]
    return bitArrayList

# Combine a list of BitArrays into a single BitArray
def joinBitArrays(bitArrayList):
    bitarray = BitArray(bin='0b')
    for item in bitArrayList:
        bitarray.append(item)
    return bitarray

def isValidPoint(pointToCheck):
    thetaReal, phiReal = bitArrayToPoint(pointToCheck)
    if thetaReal < 1e-16 or thetaReal > 1:
        return False
    elif phiReal < 1e-16 or phiReal > 1:
        return False
    return True
    
def sphericalToCartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return [x, y, z]

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
    print(f'Fitness score: {self.fitnessScore}')
    print(f'Number of points: {self.numPoints}')
    # Get each point as the spherical theta and phi values
    points = splitBitArray(self.chromosome, self.numPoints)
    points = [bitArrayToPoint(point) for point in points]
    points = np.array(points) * np.array([np.pi, 2*np.pi])
    for i, point in enumerate(points):
        print(f'Point {i}:')
        print(f'\tTheta:{point[0]} rad')
        print(f'\tPhi:{point[1]} rad')