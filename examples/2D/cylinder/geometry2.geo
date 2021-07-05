//+
L = 1;
//+
H = 0.5;
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
Circle(7) = {6, 5, 7};
//+
Circle(8) = {7, 5, 6};
//+
Line(9) = {2, 1};
//+
Line(10) = {1, 8};
//+
Line(11) = {8, 9};
//+
Line(12) = {9, 2};
//+
Curve Loop(1) = {9, 10, 11, 12};
//+
Curve Loop(2) = {8, 7};
//+
Surface(1) = {1, 2};
//+
Physical Surface("Fluid") = {1};
//+
Physical Curve("FreeSurface") = {12, 11, 10};
//+
Physical Curve("FluidInput") = {9};
//+
Physical Curve("CylinderWall") = {8, 7};
