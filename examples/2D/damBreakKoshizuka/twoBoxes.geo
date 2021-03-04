//+
L = 0.01;
//+
d = L/40;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, L, 0, d};
//+
Point(5) = {L, -L, 0, d};
//+
Point(6) = {2*L, -L, 0, d};
//+
Point(7) = {2*L, 0, 0, d};
//+
Point(8) = {2*L, L, 0, d};
//+
Point(9) = {L, L, 0, d};
//+
Point(10) = {L, 0, 0, d};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 10};
//+
Line(3) = {10, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {2, 9};
//+
Line(8) = {9, 10};
//+
Curve Loop(1) = {1, 2, -8, -7};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Boundary") = {1, 2, 3, 4, 5, 6};
//+
Physical Curve("FreeSurface") = {7, 8};
//+
Physical Surface("Fluid") = {1};
