//+
d = 0.02;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, 1, 0, d};
//+
Point(3) = {4, 0, 0, d};
//+
Point(4) = {4, 1, 0, d};
//+
Point(5) = {0, 1.1, 0, d};
//+
Point(6) = {4, 1.1, 0, d};//+
Line(1) = {5, 2};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 6};
//+
Line(6) = {4, 2};
//+
Curve Loop(1) = {2, 3, 4, 6};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Left", 7) = {1, 2};
//+
Physical Curve("Right", 8) = {5, 4};
//+
Physical Curve("Bottom", 9) = {3};
//+
Physical Curve("FreeSurface", 10) = {6};
//+
Physical Surface("Fluid", 11) = {1};
