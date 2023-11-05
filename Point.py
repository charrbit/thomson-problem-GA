import numpy as np
from bitstring import BitArray

# A point on a sphere, in a spherical coordinate system, can be described by three values
    # P = (R, theta, phi)
# For the unit sphere, R = 1 => P = (1, theta, phi) => P(theta, phi)
    # Note: theta in [0, pi] and phi in [0, 2pi)
# Then, a point on the unit sphere can be described with two real numbers, thetaReal and phiReal
    # P(theta, phi) = P(thetaReal * pi, phiReal * 2pi) => P(thetaReal, phiReal) where thetaReal and phiReal in [0, 1]

class Point:
    ''' A point on the unit sphere represented with two real-valued coordinates in the Unit Interval [0, 1] '''

    def __init__(self, thetaReal, phiReal):
        self.thetaReal = thetaReal
        self.phiReal = phiReal

    # Get the point as the real-valued thetaReal and phiReal coordinates
    def getPoint(self):
        return [self.thetaReal, self.phiReal]

    # Set the point with real-valued thetaReal and phiReal coordinates
    def setPoint(self, point):
        self.thetaReal, self.phiReal = point

    # Get the point as the spherical theta and phi coordinates
    def getPointSpherical(self):
        return np.array(self.getPoint()) * np.array([np.pi, 2*np.pi])

    # Get the point as a 128-bit BitArray (binary string)
    def getPointBitArray(self):
        # Convert the real valued-coordinates into 64-bit BitArrays
        thetaReal, phiReal = [BitArray(float=coord, length=64) for coord in self.getPoint()]
        # Join them together
        point = BitArray(bin='0b') # Empty
        point.append(thetaReal)
        point.append(phiReal)
        return point

    # Check if the coordinates of the point are within the Unit Interval
    def isValid(self):
        if self.thetaReal < 0 or self.thetaReal > 1:
            return False
        elif self.phiReal < 0 or self.phiReal > 1:
            return False
        return True