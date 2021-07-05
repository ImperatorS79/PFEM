//+
H = 0.005;
//+
W = 0.02;
//+
d = 0.0001;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, H, 0, d};
//+
Point(3) = {W, 0, 0, d};
//+
Point(4) = {W, H, 0, d};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Left") = {1};
//+
Physical Curve("FreeSurface") = {2};
//+
Physical Curve("Right") = {3};
//+
Physical Curve("Bottom") = {4};
//+
Physical Surface("Fluid") = {1};
