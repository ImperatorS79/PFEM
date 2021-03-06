//+
L = 1;
//+
H = 0.75;
//+
R = 0.02;
//+
D = 0.3;
//+
E = 0.5;
//+
d = 0.01;
//+
Point(1) = {0, -H/2, 0, d};
//+
Point(2) = {0, H/2, 0, d};
//+
Point(3) = {L, -H/2, 0, d};
//+
Point(4) = {L, H/2, 0, d};
//+
Point(5) = {D, 0, 0, d};
//+
Point(6) = {D - R, 0, 0, d};
//+
Point(7) = {D + R, 0, 0, d};
//+
Point(8) = {E, -H/2, 0, d};
//+
Point(9) = {E, H/2, 0, d};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 9};
//+
Line(3) = {9, 4};
//+
Line(4) = {1, 8};
//+
Line(5) = {8, 3};
//+
Line(6) = {8, 9};
//+
Circle(7) = {6, 5, 7};
//+
Circle(8) = {7, 5, 6};
//+
Physical Curve("FluidInput") = {1};
//+
Physical Curve("FreeSurface") = {6};
//+
Physical Curve("Walls") = {2, 3, 4, 5};
//+
Physical Curve("CylinderWall") = {8, 7};
//+
Curve Loop(1) = {1, 2, -6, -4};
//+
Curve Loop(2) = {7, 8};
//+
Surface(1) = {1, 2};
//+
Physical Surface("Fluid", 9) = {1};
