CSCI 520, Assignment 1

Name: Sarah Kurdoghlian
Date: February 12, 2023

================

# ---------- Accomplishments and Features ------------- #
- A mass spring network connecting the 512 vertices with structural, shear, and bend springs
- A collision detection system with the 2x2x2 bounding cube
- Math for computing hooks and damping forces
- A way of mapping a given point in the cube to the 3D force grid
- Trilinear interpolation of eigth external forces in the 3D force grid to determine the external force on a given point in the cube
- A collision detection system with a plane of the format ax + by + cz + d = 0

# ------------- Extra Credit - Incline Plane Collision -------------- #

1) -------
I drew an inclined plane using the given a,b,c,d constants in the world files. To decide how to draw it,
I solved for 4 points on the plane. I did this by setting x=2, y=2 and then solving for z, for example. Another point
was found by setting x=2, y=-2 and again solving for z. This made sure the plane was drawing at reasonable points inside the bounding cube.
Then, I connected my 4 points with a GL_TRIANGLE_STRIP.

To program collisions with the plane, the first step was calculating the plane normal.
If the plane normal was facing the opposite direction of the cube starting position, then I flipped the normal's direction.

At each timestep, I checked every point in the cube and calculated a vector between that point and a point
on the incline plane. For example if my cube point is "p" and my plane point is "a", I did p - a
to get the vector pointing from a to p. The important step was the doing the dot product between this vector and the plane normal.
If the result of the dot product was greater than 0, then the cube point was in the normal direction of the plane.
But, if the dot product was less than or equal to 0, then the cube point falls on the opposite side of the plane,
away from the plane normal.

Now that I know the cube collided with the incline plane, the next step was deciding what point on the plane to attach my collision spring to.
I used a formula to find the closest point on a plane to a given point. I attached the cube's collision point to the closest
point on the plane, created a spring, and added it to my list of springs. These incline plane collisions used the same
hooks and damping forces as regular collision springs.

2) --------
I made some updates to the openGL code.
I update the positions of the cameras to be further away from the cube.
I added transparency so I can add the alpha values to the RGBA colors.
