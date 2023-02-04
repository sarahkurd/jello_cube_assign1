/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);
struct point computeHooks(double k, double restLength, struct point *vectorBetweenPoints);
struct point computeDamping(double k, struct point *vectorBetweenPoints, struct point *velocityPointA, struct point *velocityPointB);
double dotProduct(struct point *vector1, struct point *vector2);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif

