from Point import *

# Create singleton random number generator
rng = np.random.default_rng()

# Helper Functions
def generateRandomPoint():
    ''' Generates a random point on the unit sphere '''
    thetaReal, phiReal = rng.random(2)
    return Point(thetaReal, phiReal)

def electroStaticPotential(point1, point2):
    ''' Calculates the electrostatic potential between two points '''
    # Convert the points to spherical coordinates then find the distance between them
    dist = np.subtract(np.array(point1.getPointSpherical()), np.array(point2.getPointSpherical()))
    return 1 / np.linalg.norm(dist)

def sphericalToCartesian(r, theta, phi):
    ''' Converts spherical coordinates into cartesian coordinates '''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return [x, y, z]