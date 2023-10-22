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
    # Create a new BitArray containing the thetaReal coordinate
    point = floatToBitArray(pointToConvert[0])
    # Add the phiReal coordinate to the point
    point.append(floatToBitArray(pointToConvert[1]))
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