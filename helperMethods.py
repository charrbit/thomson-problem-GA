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
def bitArrayToPoint(bitArrayToConvert, coordBitLength):
    # Split the BitArray point into its two component BitArrays
    point = [bitArrayToConvert.bin[:coordBitLength], bitArrayToConvert.bin[coordBitLength:]]
    point = [BitArray(bin=f'0b'{x}) for x in point]
    # Convert each component into a float
    point = [bitArrayToFloat(x) for x in point]
    return point