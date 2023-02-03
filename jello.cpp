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

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int shear=0, bend=0, structural=1, pause=0, viewingMode=0, saveScreenToFile=0;

struct world jello;

// use this list to access the spring
struct spring *springs;

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
        length = 1.0 / 7.0;
        break;
    case SHEAR_FACE:
        length = sqrt(2.0) / 7.0;
        break;
    case SHEAR_DIAGONAL:
        length = sqrt(3.0) / 7.0;
        break;
    case BEND:
        length = 2.0 / 7.0;
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
        struct point p1;
        struct point p2;
        p1.x = i;
        p1.y = j;
        p1.z = k;
        p2.x = i_connect;
        p2.y = j_connect;
        p2.z = k_connect;

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
 * Connect the jello cube with structural, shear, and bend spring
 * The list of spring starts at [springList]
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
          connectSprings(i, j, k,  i + 1, j + 1, k, SHEAR_FACE);
          connectSprings(i, j, k, i + 1, j, k + 1, SHEAR_FACE);
          connectSprings(i, j, k, i, j + 1, k + 1, SHEAR_FACE);
          connectSprings(i, j, k, i - 1, j + 1, k, SHEAR_FACE);
          connectSprings(i, j, k, i - 1, j, k + 1, SHEAR_FACE);
          connectSprings(i, j, k, i, j - 1, k + 1, SHEAR_FACE);
          connectSprings(i, j, k, i - 1, j + 1, k + 1, SHEAR_DIAGONAL);
          connectSprings(i, j, k, i - 1, j - 1, k + 1, SHEAR_DIAGONAL);
          connectSprings(i, j, k, i + 1, j + 1, k + 1, SHEAR_DIAGONAL);
          connectSprings(i, j, k, i + 1, j - 1, k + 1, SHEAR_DIAGONAL);
          // BEND
          connectSprings(i, j, k, i + 2, j, k, BEND);
          connectSprings(i, j, k, i, j + 2, k, BEND);
          connectSprings(i, j, k, i, j, k + 2, BEND);
      }
}

void clearParticleForces() {
    for (int i=0; i<=7; i++)
        for (int j=0; j<=7; j++)
            for (int k=0; k<=7; k++)
            {
                jello.particleForces[i][j][k].x = 0;
                jello.particleForces[i][j][k].y = 0;
                jello.particleForces[i][j][k].z = 0;
            }
}

/**
 * Iterate through the list of springs and calculate hooks force and damping force.
 * Store the total force on each particle in a new 3D array of forces
 */
void calculateForcesOnParticles() {
    // clear all previous forces from previous iterations
    clearParticleForces();
    struct spring *currSpring = springs;
    struct spring *previous = currSpring;
    while (currSpring != NULL) {
        struct point vectorBetweenPoints;
        struct point hooksForce;
        struct point dampingForce;

        // using my indices stored in the spring to index into the jello.p[][][] array
        struct point p1 = jello.p[(int) currSpring->p1.x][(int) currSpring->p1.y][(int) currSpring->p1.z];
        struct point v1 = jello.v[(int) currSpring->p1.x][(int) currSpring->p1.y][(int) currSpring->p1.z];
        if (currSpring->springType == COLLISION) {
            struct point v2; // collision point has no velocity
            v2.x = 0;
            v2.y = 0;
            v2.z = 0;
            struct point p2 = currSpring->p2;
            vectorBetweenPoints.x = p2.x - p1.x;
            vectorBetweenPoints.y = p2.y - p1.y;
            vectorBetweenPoints.z = p2.z - p1.z;
            hooksForce = computeHooks(jello.kCollision, currSpring->restLength, &vectorBetweenPoints);

            vectorBetweenPoints.x = p2.x - p1.x;
            vectorBetweenPoints.y = p2.y - p1.y;
            vectorBetweenPoints.z = p2.z - p1.z;
            dampingForce = computeDamping(jello.dCollision, &vectorBetweenPoints, &v1, &v2);
        } else {
            struct point p2 = jello.p[(int) currSpring->p2.x][(int) currSpring->p2.y][(int) currSpring->p2.z];
            struct point v2 = jello.v[(int) currSpring->p2.x][(int) currSpring->p2.y][(int) currSpring->p2.z];
            vectorBetweenPoints.x = p2.x - p1.x;
            vectorBetweenPoints.y = p2.y - p1.y;
            vectorBetweenPoints.z = p2.z - p1.z;
            hooksForce = computeHooks(jello.kElastic, currSpring->restLength, &vectorBetweenPoints);

            // reset the vector between points because it was modified in the hook's calculations
            vectorBetweenPoints.x = p2.x - p1.x;
            vectorBetweenPoints.y = p2.y - p1.y;
            vectorBetweenPoints.z = p2.z - p1.z;
            dampingForce = computeDamping(jello.dElastic, &vectorBetweenPoints, &v1, &v2);
            jello.particleForces[(int) currSpring->p2.x][(int) currSpring->p2.y][(int) currSpring->p2.z].x +=
                    hooksForce.x + dampingForce.x;
            jello.particleForces[(int) currSpring->p2.x][(int) currSpring->p2.y][(int) currSpring->p2.z].y +=
                    hooksForce.y + dampingForce.y;
            jello.particleForces[(int) currSpring->p2.x][(int) currSpring->p2.y][(int) currSpring->p2.z].z +=
                    hooksForce.z + dampingForce.z;
        }
        jello.particleForces[(int) currSpring->p1.x][(int) currSpring->p1.y][(int) currSpring->p1.z].x +=
                -1.0 * hooksForce.x + -1.0 * dampingForce.x;
        jello.particleForces[(int) currSpring->p1.x][(int) currSpring->p1.y][(int) currSpring->p1.z].y +=
                -1.0 * hooksForce.y + -1.0 * dampingForce.y;
        jello.particleForces[(int) currSpring->p1.x][(int) currSpring->p1.y][(int) currSpring->p1.z].z +=
                -1.0 * hooksForce.z + -1.0 * dampingForce.z;

        // Unlink the collision node from the linked list.
        // We only want it in our list of springs once, at the time of collision.
        if (currSpring->springType == COLLISION) {
            previous->next = currSpring->next;
        } else {
            previous = currSpring;
        }
        currSpring = currSpring->next;
    }
}

void findCollisionPoint(int i, int j, int k, struct point *collisionPoint) {
    struct point p = jello.p[i][j][k];
    if (p.x > 2.0) {
        collisionPoint->x = 2.0;
        collisionPoint->y = p.y;
        collisionPoint->z = p.z;
    } else if (p.x < -2.0) {
        collisionPoint->x = -2.0;
        collisionPoint->y = p.y;
        collisionPoint->z = p.z;
    } else if (p.y > 2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = 2.0;
        collisionPoint->z = p.z;
    } else if (p.y < -2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = -2.0;
        collisionPoint->z = p.z;
    } else if (p.z > 2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = p.y;
        collisionPoint->z = 2.0;
    } else if (p.z < -2.0) {
        collisionPoint->x = p.x;
        collisionPoint->y = p.y;
        collisionPoint->z = -2.0;
    }
}

void checkCollision() {
    for (int i=0; i<=7; i++)
        for (int j=0; j<=7; j++)
            for (int k=0; k<=7; k++)
            {
                struct point p = jello.p[i][j][k];
                if (p.x > 2.0 || p.x < -2.0 || p.y > 2.0 || p.y < -2.0 || p.z > 2.0 || p.z < -2.0) {
                    // create a point, the collision point, that will have the
                    struct point cubePoint;
                    struct point collisionPoint;
                    cubePoint.x = i;
                    cubePoint.y = j;
                    cubePoint.z = k;
                    findCollisionPoint(i, j, k, &collisionPoint);

                    struct spring *collisionSpring;
                    collisionSpring = (struct spring*) malloc(sizeof(struct spring));
                    collisionSpring->p1 = cubePoint;
                    collisionSpring->p2 = collisionPoint;
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

void myinit()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(90.0,1.0,0.01,1000.0);

  // set background color to grey
  glClearColor(0.0, 0.5, 0.5, 0.0);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

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
  //GLfloat aGa[] = { 0.0, 0.0, 0.0, 1.0 };

   // light 's ambient, diffuse, specular
  // GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd0[] = { 1.0, 1.0, 1.0, 1.0 };
  // GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 };

  // GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 };
  // GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 };

  // GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd2[] = { 1.0, 1.0, 0.0, 1.0 };
  // GLfloat lKs2[] = { 1.0, 1.0, 0.0, 1.0 };

  // GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd3[] = { 0.0, 1.0, 1.0, 1.0 };
  // GLfloat lKs3[] = { 0.0, 1.0, 1.0, 1.0 };

  // GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd4[] = { 0.0, 0.0, 1.0, 1.0 };
  // GLfloat lKs4[] = { 0.0, 0.0, 1.0, 1.0 };

  // GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 };
  // GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 };

  // GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd6[] = { 1.0, 1.0, 1.0, 1.0 };
  // GLfloat lKs6[] = { 1.0, 1.0, 1.0, 1.0 };

  // GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
  // GLfloat lKd7[] = { 0.0, 1.0, 1.0, 1.0 };
  // GLfloat lKs7[] = { 0.0, 1.0, 1.0, 1.0 };
  
  // light 's ambient, diffuse, specular
  GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd0[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs0[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd2[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs2[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd3[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs3[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd4[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs4[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd5[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs5[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd6[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs6[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd7[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs7[] = { 1.0, 0.0, 0.0, 1.0 };

  // light positions and directions
  GLfloat lP0[] = { -1.999, -1.999, -1.999, 1.0 };
  GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
  GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
  GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
  GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
  GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
  GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
  GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };
  
  // jelly material color
  GLfloat mKa[] = { 0.0, 0.7, 0.0, 1.0 };
  //GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat mKd[] = { 0.5, 0.5, 0.5, 1.0 };
  //GLfloat mKd[] = { 0.3, 0.3, 0.3, 1.0 };
  GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mKe[] = { 0.1, 0.1, 0.1, 1.0 };
  //GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 };

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
  showBoundingBox();
 
  glutSwapBuffers();
}

void doIdle()
{
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
    checkCollision();
    calculateForcesOnParticles();
    if (strcmp(jello.integrator, "RK4")) {
      RK4(&jello);
    } else {
      Euler(&jello);
    }
  }

  glutPostRedisplay();
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

