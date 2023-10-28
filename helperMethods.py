import numpy as np
from bitstring import BitArray

# Generate two random numbers in [0, 1] that represent the coordinates of a point on the unit sphere
def generateRandomPoint(rng):
    thetaReal, phiReal = rng.random(2)
    return [thetaReal, phiReal]

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
    ''' Calculates the electrostatic potential between two points where
        the points are given as BitArrays each containing a thetaReal and phiReal value '''
    # Extract thetaReal and phiReal values for both points
    P1 = bitArrayToPoint(point1)
    P2 = bitArrayToPoint(point2)
    # Convert the points to numpy arrays to find the vector difference
    diff = np.subtract(np.array(P1), np.array(P2))
    # Convert thetaReal and phiReal values to spherical coordinates theta and phi
    diff *= np.array([np.pi, 2*np.pi])
    return 1 / np.linalg.norm(diff)

def isValidPoint(pointToCheck):
    thetaReal, phiReal = bitArrayToPoint(pointToCheck) 
    if (0 > thetaReal > 1) or (0 > phiReal > 1):
        return False
    return True