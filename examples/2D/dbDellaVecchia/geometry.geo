//+
d = 0.01;
//+
LBox = 0.98;
//+
HBox = 0.35;
//+
L0 = 0.32;
//+
H0 = 0.267;
//+
Point(0) = {0, 0, 0, d};
//+
Point(1) = {0, HBox, 0, d};
//+
Point(2) = {LBox, HBox, 0, d};
//+
Point(3) = {LBox, 0, 0, d};
//+
Point(4) = {0, H0, 0, d};
//+
Point(5) = {L0, H0, 0, d};
//+
Point(6) = {L0, 0, 0, d};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 0};
//+
Line(3) = {0, 6};
//+
Line(4) = {6, 3};
//+
Line(5) = {3, 2};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 6};
//+
Curve Loop(1) = {2, 3, -7, -6};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Walls", 8) = {1, 2, 3, 4, 5};
//+
Physical Curve("FreeSurface", 9) = {6, 7};
//+
Physical Surface("Fluid", 10) = {1};
