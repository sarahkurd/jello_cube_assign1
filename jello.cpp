/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

  Your name:
  Sarah Kurdoghlian

*/

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 6;

//PerformanceCounter performanceCounter = PerformanceCounter();

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int shear=0, bend=0, structural=1, pause=0, viewingMode=1, saveScreenToFile=0;

struct world jello;

// use this list to access the spring
struct spring *springs = NULL;

double *intervals;

struct point pointsOnInclinePlane[3];
struct point inclinePlaneNormal;

int windowWidth, windowHeight;

/**
 * Get the rest length of the input spring type
 *
 * @param s The type of SPRING, an enum
 * @return A double representing the given spring rest length
 */
double getSpringLength(SPRING s){
  double length;
  switch (s) {
    case STRUCTURAL:
        length = 1.0 / 7;
        break;
    case SHEAR_FACE:
        length = sqrt(2) / 7;
        break;
    case SHEAR_DIAGONAL:
        length = sqrt(3) / 7;
        break;
    case BEND:
        length = 2.0 / 7;
        break;
    case COLLISION:
        length = 0;
        break;
  }
  return length;
}

/**
 * Connect two points with the given spring
 * Stores the spring in a linked list called [spring]
 *
 * @param currPoint The origin point
 * @param i_connect the x coordinate of the point we are connecting to
 * @param j_connect the y coordinate of the point we are connecting to
 * @param k_connect the z coordinate of the point we are connecting to
 * @param springType The [SPRING] enum to be used to get coefficient values during force calculations
 */
void connectSprings(int i, int j, int k, int i_connect, int j_connect, int k_connect, SPRING springType) {
    if (!((i_connect>7) || (i_connect<0) 
          || (j_connect>7) || (j_connect<0) 
          || (k_connect>7) || (k_connect<0))) 
    {
        // These are the indices of the points in the jello.p[][][] list
        struct index p1 = {i, j, k};
        struct index p2 = {i_connect, j_connect, k_connect};

        struct spring *springNew;
        springNew = (struct spring*) malloc(sizeof(struct spring));
        springNew->p1 = p1;
        springNew->p2 = p2;
        springNew->springType = springType;
        springNew->restLength = getSpringLength(springType);
        springNew->next = NULL;

        // add new spring at the end of the linked list
        if (springs == NULL) {
            springs = springNew;
        } else {
            struct spring *temp = springs;
            while (temp->next != NULL) {
                temp = temp->next;
            }
            temp->next = springNew;
        }
    }
}

/**
 * Connect the jello cube with structural, shear, and bend springs.
 * The list of springs starts at [springList].
 */
void buildSpringNetwork() {
  for (int i=0; i<=7; i++)
    for (int j=0; j<=7; j++)
      for (int k=0; k<=7; k++)
      {
          // STRUCTURAL
          connectSprings(i, j, k, i + 1, j, k, STRUCTURAL);
          connectSprings(i, j, k, i, j + 1, k, STRUCTURAL);
          connectSprings(i, j, k,  i, j, k + 1, STRUCTURAL);

          // SHEAR
          connectSprings(i, j, k, i + 1, j + 1, k, SHEAR_FACE);
          connectSprings(i, j, k, i + 1, j - 1, k, SHEAR_FACE);
          connectSprings(i, j, k, i, j + 1, k + 1, SHEAR_FACE);
          connectSprings(i, j, k, i, j + 1, k - 1, SHEAR_FACE);
          connectSprings(i, j, k, i + 1, j, k + 1, SHEAR_FACE);
          connectSprings(i, j, k, i + 1, j, k - 1, SHEAR_FACE);

          // SHEAR DIAGONAL
          connectSprings(i, j, k, i + 1, j + 1, k + 1, SHEAR_DIAGONAL);
          connectSprings(i, j, k, i + 1, j - 1, k + 1, SHEAR_DIAGONAL);
          connectSprings(i, j, k, i + 1, j + 1, k - 1, SHEAR_DIAGONAL);
          connectSprings(i, j, k, i + 1, j - 1, k - 1, SHEAR_DIAGONAL);

          // BEND
          connectSprings(i, j, k, i + 2, j, k, BEND);
          connectSprings(i, j, k, i, j + 2, k, BEND);
          connectSprings(i, j, k, i, j, k + 2, BEND);
      }
}

/**
 * Given the resolution of the force field,
 * create a list of the intervals in the 3D force grid.
 */
void createIntervalsFromResolution() {
    intervals = (double *)malloc((jello.resolution) * sizeof(double));
    for (int i = 0; i < jello.resolution; i++) {
        double coordinate = -2.0 + i * 4.0 / (jello.resolution - 1);
        intervals[i] = coordinate;
    }
}

/**
 * Given a coordinate, determine what interval it falls in.
 * This will be the min and max coordinates of the point in either the x,y,z dimensions.
 * @param xyz A double coordinate
 * @param xInterval Used to set the minimum and maximum values of the interval for this given coordinate
 */
void getForceFieldInterval(double xyz, double xInterval[2]) {
    for (int i = 0; i < jello.resolution; i++) {
        if (xyz >= intervals[i] && xyz <= intervals[i + 1]) {
            xInterval[0] = intervals[i];
            xInterval[1] = intervals[i + 1];
            break;
        }
    }
}

/**
 * Given a coordinate, use the given formula to determine which index it is in the jello.forceField array
 * @param coordinate A double coordinate value
 * @return An integer representing an index into the jello.forceField array
 */
int getIndexFromCoordinate(double coordinate) {
    return ((jello.resolution - 1) * (coordinate + 2)) / 4;
}

/**
 * Iterate through each point in the cube and interpolate the forces around it.
 * Assign the interpolated value as the external force field on that point,
 * and add to the force in jello.particleForces[i][j][k] .
 */
void getExternalForceAtPoints() {
    for (int i=0; i<=7; i++)
        for (int j=0; j<=7; j++)
            for (int k=0; k<=7; k++)
            {
                struct point currPoint = jello.p[i][j][k];
                // if point is not outside of the bounding cube, then we can apply external force on it
                if (!(currPoint.x < -2.0 || currPoint.x > 2.0 || currPoint.y < -2.0
                        || currPoint.y > 2.0 || currPoint.z < -2.0 || currPoint.z > 2.0)) {
                    double xInterval[2];
                    double yInterval[2];
                    double zInterval[2];
                    getForceFieldInterval(currPoint.x, xInterval);
                    getForceFieldInterval(currPoint.y, yInterval);
                    getForceFieldInterval(currPoint.z, zInterval);

                    int index_x_min = getIndexFromCoordinate(xInterval[0]);
                    int index_x_max = getIndexFromCoordinate(xInterval[1]);
                    int index_y_min = getIndexFromCoordinate(yInterval[0]);
                    int index_y_max = getIndexFromCoordinate(yInterval[1]);
                    int index_z_min = getIndexFromCoordinate(zInterval[0]);
                    int index_z_max = getIndexFromCoordinate(zInterval[1]);

                    struct point force1, force2, force3, force4, force5, force6, force7, force8;

                    // force field 1: x_min, y_min, z_min coordinates
                    force1.x = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_min].x;
                    force1.y = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_min].y;
                    force1.z = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_min].z;

                    // force field 2: x_min, y_min, z_max coordinates
                    force2.x = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_max].x;
                    force2.y = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_max].y;
                    force2.z = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_max].z;

                    // force field 3: x_min, y_max, z_min coordinates
                    force3.x = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_min].x;
                    force3.y = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_min].y;
                    force3.z = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_min].z;

                    // force field 4: x_min, y_max, z_max coordinates
                    force4.x = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_max].x;
                    force4.y = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_max].y;
                    force4.z = jello.forceField[index_x_min * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_max].z;

                    // force field 5: x_max, y_min, z_min coordinates
                    force5.x = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_min].x;
                    force5.y = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_min].y;
                    force5.z = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_min].z;

                    // force field 6: x_max, y_min, z_max coordinates
                    force6.x = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_max].x;
                    force6.y = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_max].y;
                    force6.z = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_min * jello.resolution + index_z_max].z;

                    // force field 7: x_max, y_max, z_min coordinates
                    force7.x = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_min].x;
                    force7.y = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_min].y;
                    force7.z = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_min].z;

                    // force field 8: x_max, y_max, z_max coordinates
                    force8.x = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_max].x;
                    force8.y = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_max].y;
                    force8.z = jello.forceField[index_x_max * jello.resolution * jello.resolution
                                                + index_y_max * jello.resolution + index_z_max].z;

                    // Now, we do trilinear interpolation
                    double alpha = (currPoint.x - xInterval[0]) / (xInterval[1] - xInterval[0]);
                    double beta = (currPoint.y - yInterval[0]) / (yInterval[1] - yInterval[0]);
                    double gamma = (currPoint.z - zInterval[0]) / (zInterval[1] - zInterval[0]);

                    struct point finalTrilinearForce;
                    finalTrilinearForce.x =
                            ((1 - alpha) * (1 - beta) * (1 - gamma) * force1.x)
                            + ((alpha) * (1 - beta) * (1 - gamma) * force5.x)
                            + ((alpha) * (beta) * (1 - gamma) * force7.x)
                            + ((1 - alpha) * (beta) * (1 - gamma) * force3.x)
                            + ((1 - alpha) * (1 - beta) * (gamma) * force2.x)
                            + ((1 - alpha) * (beta) * (gamma) * force4.x)
                            + ((alpha) * (1 - beta) * (gamma) * force6.x)
                            + ((alpha) * (beta) * (gamma) * force8.x);

                    finalTrilinearForce.y =
                            ((1 - alpha) * (1 - beta) * (1 - gamma) * force1.y)
                            + ((alpha) * (1 - beta) * (1 - gamma) * force5.y)
                            + ((alpha) * (beta) * (1 - gamma) * force7.y)
                            + ((1 - alpha) * (beta) * (1 - gamma) * force3.y)
                            + ((1 - alpha) * (1 - beta) * (gamma) * force2.y)
                            + ((1 - alpha) * (beta) * (gamma) * force4.y)
                            + ((alpha) * (1 - beta) * (gamma) * force6.y)
                            + ((alpha) * (beta) * (gamma) * force8.y);

                    finalTrilinearForce.z =
                            ((1 - alpha) * (1 - beta) * (1 - gamma) * force1.z)
                            + ((alpha) * (1 - beta) * (1 - gamma) * force5.z)
                            + ((alpha) * (beta) * (1 - gamma) * force7.z)
                            + ((1 - alpha) * (beta) * (1 - gamma) * force3.z)
                            + ((1 - alpha) * (1 - beta) * (gamma) * force2.z)
                            + ((1 - alpha) * (beta) * (gamma) * force4.z)
                            + ((alpha) * (1 - beta) * (gamma) * force6.z)
                            + ((alpha) * (beta) * (gamma) * force8.z);

                    // Add this force to the particle forces
                    jello.particleForces[i][j][k].x += finalTrilinearForce.x;
                    jello.particleForces[i][j][k].y += finalTrilinearForce.y;
                    jello.particleForces[i][j][k].z += finalTrilinearForce.z;
                }
            }
}

/**
 * Given the position of a point in the jello cube, depending on where that point collided,
 * create the collision point by clamping to -2 or 2 and setting it as the collision point.
 */
void findCollisionPoint(int i, int j, int k, struct point *collisionPoint) {
    struct point p = jello.p[i][j][k];
    if (p.x >= 2.0) {
        collisionPoint->x = 2.0;
        collisionPoint->y = p.y;
        collisionPoint->z = p.z;
    } else if (p.x <= -2.0) {
        collisionPoint->x = -2.0;
        collisionPoint->y = p.y;
        collisionPoint->z = p.z;
    } else if (p.y >= 2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = 2.0;
        collisionPoint->z = p.z;
    } else if (p.y <= -2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = -2.0;
        collisionPoint->z = p.z;
    } else if (p.z >= 2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = p.y;
        collisionPoint->z = 2.0;
    } else if (p.z <= -2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = p.y;
        collisionPoint->z = -2.0;
    }
}

/**
 * If incline plane exists, call this function once.
 * Determines 3 points on the incline plane.
 */
void getPointsOnInclinePlane() {
    // set 2 coordinates to 0, and solve for the other coordinate
    double z = (-1.0 * jello.d) / jello.c;
    double y = (-1.0 * jello.d) / jello.b;
    double x = (-1.0 * jello.d) / jello.a;
    struct point coord1 = {0.0, 0.0, z};
    struct point coord2 = {0.0, y, 0.0};
    struct point coord3 = {x, 0.0, 0.0};
    pointsOnInclinePlane[0] = coord1;
    pointsOnInclinePlane[1] = coord2;
    pointsOnInclinePlane[2] = coord3;
}

/**
 * If incline plane exists, run this function once.
 * Returns a vector that is normal to the incline plane.
 * @return  A struct point representing the normal vector.
 */
struct point getPlaneNormal() {
    struct point normal;
    struct point vec1;
    struct point vec2;
    pDIFFERENCE(pointsOnInclinePlane[1], pointsOnInclinePlane[0], vec1)
    pDIFFERENCE(pointsOnInclinePlane[2], pointsOnInclinePlane[0], vec2)
    CROSSPRODUCTp(vec1, vec2, normal);
    return normal;
}

/**
 * Given a point on the cube, first find the vector between that point
 * and a point on the plane.
 * Do a dot product with the normal vector of the plane to determine
 * if the point is on the normal side or non normal side of the plane.
 */
bool isPointCollidedWithPlane(struct point cubePoint) {
    struct point vec;
    pDIFFERENCE(cubePoint, pointsOnInclinePlane[0], vec)
    double dot = dotProduct(vec, inclinePlaneNormal);
    if (dot <= 0) {
        return true;
    } else {
        return false;
    }
}

/**
 * Find the closest point on the incline plane to the given cube point
 * We get the equation of a line that goes through  the cube point and is perpendicular to the plane.
 * @param cubePoint
 * @return
 */
struct point getClosestPointOnPlaneToConnectSpring(struct point cubePoint) {
    double x = cubePoint.x * jello.a;
    double xt = inclinePlaneNormal.x * jello.a;
    double y = cubePoint.y * jello.b;
    double yt = inclinePlaneNormal.y * jello.b;
    double z = cubePoint.z * jello.c;
    double zt = inclinePlaneNormal.z * jello.c;
    double finalT = ((-1.0 * jello.d) - (x + y + z)) / (xt + yt + zt);

    struct point closestPoint;
    closestPoint.x = inclinePlaneNormal.x * finalT + cubePoint.x;
    closestPoint.y = inclinePlaneNormal.y * finalT + cubePoint.y;
    closestPoint.z = inclinePlaneNormal.z * finalT + cubePoint.z;
    return closestPoint;
}

/**
 * Check for any collisions between the cube and the bounding box.
 * Create a new spring if there is a collisions,
 * and the spring to our list of springs.
 */
void checkCollision() {
    for (int i=0; i<=7; i++)
        for (int j=0; j<=7; j++)
            for (int k=0; k<=7; k++)
            {
                struct point p = jello.p[i][j][k];
                bool planeCollision = false;
                if (jello.incPlanePresent == 1) {
                    planeCollision = isPointCollidedWithPlane(p);
                }
                if (p.x >= 2.0 || p.x <= -2.0 || p.y >= 2.0
                        || p.y <= -2.0 || p.z >= 2.0 || p.z <= -2.0 ||
                        planeCollision
                   ) {
                    struct index cubePoint = {i, j, k};
                    struct point collisionPoint = {0};
                    if (planeCollision) {
                        collisionPoint = getClosestPointOnPlaneToConnectSpring(p);
                    } else {
                        findCollisionPoint(i, j, k, &collisionPoint);
                    }

                    struct spring *collisionSpring;
                    collisionSpring = (struct spring*) malloc(sizeof(struct spring));
                    collisionSpring->p1 = cubePoint;
                    collisionSpring->collisionPoint = collisionPoint;
                    collisionSpring->restLength = getSpringLength(COLLISION);
                    collisionSpring->springType = COLLISION;
                    collisionSpring->next = NULL;

                    // add new spring at the end of the linked list
                    if (springs == NULL) {
                        springs = collisionSpring;
                    } else {
                        struct spring *temp = springs;
                        while (temp->next != NULL) {
                            temp = temp->next;
                        }
                        temp->next = collisionSpring;
                    }
                }
            }
}

/**
 * Clears the force on each particle from the previous iteration
 * We do not want to accumulate forces over time.
 * Each time step is a new set of calculations.
 */
void clearParticleForces() {
    for (int i=0; i<=7; i++)
        for (int j=0; j<=7; j++)
            for (int k=0; k<=7; k++)
            {
                jello.particleForces[i][j][k].x = 0.0;
                jello.particleForces[i][j][k].y = 0.0;
                jello.particleForces[i][j][k].z = 0.0;
            }
}

/**
 * Iterate through the list of springs and calculate hooks force and damping force.
 * Store the total force on each particle in a new 3D array of forces.
 * Lastly, if external forces exist, call the function to calculate external forces at the end.
 */
void calculateForcesOnParticles() {
    struct spring *currSpring = springs;
    struct spring *previous = currSpring;
    struct point hooksForce = {0.0};
    struct point dampingForce = {0.0};
    struct point vectorBetweenPoints = {0.0};
    while (currSpring != NULL) {
        // using my indices stored in the spring to index into the jello.p[][][] array
        struct point p1 = jello.p[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k];
        struct point v1 = jello.v[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k];
        if (currSpring->springType == COLLISION) {
            struct point v2 = {0.0}; // collision point has no velocity
            struct point p2 = currSpring->collisionPoint;
            vectorBetweenPoints = difference(p2, p1);
            hooksForce = computeHooks(jello.kCollision, currSpring->restLength, vectorBetweenPoints);
            dampingForce = computeDamping(jello.dCollision, vectorBetweenPoints, v2, v1);
            jello.particleForces[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k].x +=
                    (-hooksForce.x) + (-dampingForce.x);
            jello.particleForces[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k].y +=
                    (-hooksForce.y) + (-dampingForce.y);
            jello.particleForces[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k].z +=
                    (-hooksForce.z) + (-dampingForce.z);
        } else {
            struct point p2 = jello.p[currSpring->p2.i][currSpring->p2.j][currSpring->p2.k];
            struct point v2 = jello.v[currSpring->p2.i][currSpring->p2.j][currSpring->p2.k];
            vectorBetweenPoints = difference(p2, p1);
            hooksForce = computeHooks(jello.kElastic, currSpring->restLength, vectorBetweenPoints);
            dampingForce = computeDamping(jello.dElastic, vectorBetweenPoints, v2, v1);
            jello.particleForces[currSpring->p2.i][currSpring->p2.j][currSpring->p2.k].x +=
                    hooksForce.x + dampingForce.x;
            jello.particleForces[currSpring->p2.i][currSpring->p2.j][currSpring->p2.k].y +=
                    hooksForce.y + dampingForce.y;
            jello.particleForces[currSpring->p2.i][currSpring->p2.j][currSpring->p2.k].z +=
                    hooksForce.z + dampingForce.z;

            jello.particleForces[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k].x +=
                    (-hooksForce.x) + (-dampingForce.x);
            jello.particleForces[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k].y +=
                    (-hooksForce.y) + (-dampingForce.y);
            jello.particleForces[currSpring->p1.i][currSpring->p1.j][currSpring->p1.k].z +=
                    (-hooksForce.z) + (-dampingForce.z);
        }

        // Unlink the collision node from the linked list.
        // We only want it in our list of springs once, at the time of collision.
        if (currSpring->springType == COLLISION) {
            previous->next = currSpring->next;
        } else {
            previous = currSpring;
        }
        currSpring = currSpring->next;
    }
    // Last step to get total force on particles
    // Add the external force field force on each particle
    if (jello.resolution != 0) {
        getExternalForceAtPoints();
    }
}

void myinit()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90.0,1.0,0.01,1000.0);

    // set background color to grey
    glClearColor(0.0, 0.5, 0.5, 0.0);

    //glCullFace(GL_BACK);
    //glEnable(GL_CULL_FACE);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    return;
}

void reshape(int w, int h) 
{
    // Prevent a divide by zero, when h is zero.
    // You can't make a window of zero height.
    if(h == 0)
    h = 1;

    glViewport(0, 0, w, h);

    // Reset the coordinate system before modifying
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Set the perspective
    double aspectRatio = 1.0 * w / h;
    gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    windowWidth = w;
    windowHeight = h;

    glutPostRedisplay();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // camera parameters are Phi, Theta, R
    gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
            0.0,0.0,0.0, 0.0,0.0,1.0);


    /* Lighting */
    /* You are encouraged to change lighting parameters or make improvements/modifications
     to the lighting model .
     This way, you will personalize your assignment and your assignment will stick out.
    */

    // global ambient light
    GLfloat aGa[] = { 1.0, 1.0, 1.0, 1.0 };

    // light 's ambient, diffuse, specular
    GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd0[] = { 0.0, 0.5, 0.7, 1.0 };
    GLfloat lKs0[] = { 0.0, 0.5, 0.7, 1.0 };

    GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd1[] = { 1.0, 0.0, 1.0, 1.0 };
    GLfloat lKs1[] = { 1.0, 0.0, 1.0, 1.0 };

    GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd2[] = { 0.0, 1.0, 1.0, 1.0 };
    GLfloat lKs2[] = { 0.0, 1.0, 1.0, 1.0 };

    GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd3[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat lKs3[] = { 1.0, 1.0, 1.0, 1.0 };

    GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd4[] = { 1.0, 0.0, 0.0, 1.0 };
    GLfloat lKs4[] = { 1.0, 0.0, 0.0, 1.0 };

    GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 };
    GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 };

    GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd6[] = { 1.0, 1.0, 0.0, 1.0 };
    GLfloat lKs6[] = { 1.0, 1.0, 0.0, 1.0 };

    GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat lKd7[] = { 1.0, 0.0, 0.0, 1.0 };
    GLfloat lKs7[] = { 1.0, 0.0, 0.0, 1.0 };

    GLfloat lP0[] = { -3.999, -3.999, -3.999, 1.0 };
    GLfloat lP1[] = { 2.999, -2.999, -2.999, 1.0 };
    GLfloat lP2[] = { 2.999, 2.999, -2.999, 1.0 };
    GLfloat lP3[] = { -2.999, 2.999, -2.999, 1.0 };
    GLfloat lP4[] = { -2.999, -2.999, 2.999, 1.0 };
    GLfloat lP5[] = { 2.999, -2.999, 2.999, 1.0 };
    GLfloat lP6[] = { 3.999, 3.999, 3.999, 1.0 };
    GLfloat lP7[] = { 0.0, 0.0, 2.999, 1.0 };
  
    // jelly material color
    GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mKd[] = { 0.4, 0.4, 0.4, 0.7 };
    GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 };

    /* set up lighting */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
    glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

    // set up cube color
    glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
    glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
    glMaterialf(GL_FRONT, GL_SHININESS, 120);
    
    // macro to set up light i
    #define LIGHTSETUP(i)\
    glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
    glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
    glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
    glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
    glEnable(GL_LIGHT##i)

    LIGHTSETUP (0);
    LIGHTSETUP (1);
    LIGHTSETUP (2);
    LIGHTSETUP (3);
    LIGHTSETUP (4);
    LIGHTSETUP (5);
    LIGHTSETUP (6);
    LIGHTSETUP (7);

    // enable lighting
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    // show the cube
    showCube(&jello);

    glDisable(GL_LIGHTING);

    // show the bounding box
    showBoundingBox(&jello);

    glutSwapBuffers();
}

void doIdle()
{
//    performanceCounter.StopCounter();
//    printf("Frames per second = %f\n", 1.0/performanceCounter.GetElapsedTime());
//    performanceCounter.StartCounter();

    char s[20]="picxxxx.ppm";
    int i;

    // save screen to file
    s[3] = 48 + (sprite / 1000);
    s[4] = 48 + (sprite % 1000) / 100;
    s[5] = 48 + (sprite % 100 ) / 10;
    s[6] = 48 + sprite % 10;

    if (saveScreenToFile==1)
    {
    saveScreenshot(windowWidth, windowHeight, s);
    saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
    sprite++;
    }

    if (sprite >= 300) // allow only 300 snapshots
    {
    exit(0);
    }

    if (pause == 0)
    {
        // insert code which appropriately performs one step of the cube simulation:
        for (i = 1; i <= jello.n; i++) {
            clearParticleForces(); // clear all previous forces from previous iterations
            checkCollision();
            calculateForcesOnParticles();
            if (jello.integrator[0]=='E') {// Euler
                Euler(&jello);
            } else if (jello.integrator[0]=='R') { // RK4
                RK4(&jello);
            }
        }
    }
    glutPostRedisplay();
}

/**
 * If the normal vector of the plane is facing in the opposite
 * direction of the jello cube to start the simulation,
 * then flip the direction of the plane's normal vector.
 */
void checkIfPlaneNormalCorrectToStart() {
    struct point cubePoint = jello.p[0][0][0];
    struct point vec;
    vec.x = cubePoint.x - pointsOnInclinePlane[0].x;
    vec.y = cubePoint.y - pointsOnInclinePlane[0].y;
    vec.z = cubePoint.z - pointsOnInclinePlane[0].z;
    double dot = dotProduct(vec, inclinePlaneNormal);
    if (dot < 0) { //this means th cube is on the opposite side of the normal to start
        // flip the direction of the normal
        inclinePlaneNormal.x = -1.0 * inclinePlaneNormal.x;
        inclinePlaneNormal.y = -1.0 * inclinePlaneNormal.y;
        inclinePlaneNormal.z = -1.0 * inclinePlaneNormal.z;
    }
}

int main (int argc, char ** argv)
{
    if (argc<2)
    {
    printf ("Oops! You didn't say the jello world file!\n");
    printf ("Usage: %s [worldfile]\n", argv[0]);
    exit(0);
    }

    readWorld(argv[1],&jello);

    // build spring network once we have loaded the world
    buildSpringNetwork();
    createIntervalsFromResolution();
    if (jello.incPlanePresent == 1) {
      getPointsOnInclinePlane();
      inclinePlaneNormal = getPlaneNormal();
      // now, check that the cube is on the side of the normal to start.
      // If it is not, then flip the direction of the plane normal.
      checkIfPlaneNormalCorrectToStart();
    }

    glutInit(&argc,argv);

    /* double buffered window, use depth testing, 640x480 */
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    windowWidth = 640;
    windowHeight = 480;
    glutInitWindowSize (windowWidth, windowHeight);
    glutInitWindowPosition (0,0);
    glutCreateWindow ("Jello cube");

    /* tells glut to use a particular display function to redraw */
    glutDisplayFunc(display);

    /* replace with any animate code */
    glutIdleFunc(doIdle);

    /* callback for mouse drags */
    glutMotionFunc(mouseMotionDrag);

    /* callback for window size changes */
    glutReshapeFunc(reshape);

    /* callback for mouse movement */
    glutPassiveMotionFunc(mouseMotion);

    /* callback for mouse button changes */
    glutMouseFunc(mouseButton);

    /* register for keyboard events */
    glutKeyboardFunc(keyboardFunc);

    /* do initialization */
    myinit();

    /* forever sink in the black hole */
    glutMainLoop();

    return(0);
}

