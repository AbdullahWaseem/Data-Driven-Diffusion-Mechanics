//
L = 1; 			//length of the unit-cell, the units are arbitrary, multiply accordingly		
r = 0.25;		//the radius of the circle
ne = 1; 		//the number of elements in each direction
mesh_size = L/ne;	//this the general size of the mesh. change this value if more elements are requried arround a point

//points are defined first
Point(1) = {0,0,0,mesh_size}; 		//points for square
Point(2) = {L*10,0,0,mesh_size};		//
Point(3) = {L*10,L,0,mesh_size};		//
Point(4) = {0,L,0,mesh_size};		//
//Point(5) = {L,L,0,mesh_size};	

//lines are defined second
Line(1) = {1,2};	//lines for square
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

//line loops are defined third
Line Loop(1) = {1,2,3,4};

//plane surface OR surfaces are defined out of the line loops fourth
Plane Surface(1) = {1};

//points to refine/define mesh on a point/line is defined fifth
//Point{5} In Surface{1};

//apply physical tags for boundary conditions and material distributions

// for boundary conditions
base  = 11;
right = 12;
top   = 13;
left  = 14;

// physical tags for boundary conditions
Physical Line(base)  = {1};
Physical Line(right) = {2};
Physical Line(top)   = {3};
Physical Line(left)  = {4};

// for material properties
matrix    = 1;

// physical tags for material properties
Physical Surface(matrix) = {1};

