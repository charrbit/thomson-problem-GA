from Point import *

# Create random number generator
rng = np.random.default_rng()

def generateRandomPoint():
    ''' Generates a random point on the unit sphere '''
    thetaReal, phiReal = rng.random(2)
    return Point(thetaReal, phiReal)

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

def electroStaticPotential(point1, point2):
    ''' Calculates the electrostatic potential between two points  '''
    # Convert the points to spherical coordinates then find the distance between them
    dist = np.subtract(np.array(point1.getPointSpherical()), np.array(point2.getPointSpherical()))
    return 1 / np.linalg.norm(dist)

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