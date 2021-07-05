//+
d = 0.05;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {-0.62035, 0, 0, d};
//+
Point(3) = {0.62035, 0, 0, d};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 2};
//+
Curve Loop(1) = {2, 1};
//+
Surface(1) = {1};
//+
Surface(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("FreeSurface", 3) = {2, 1};
//+
Physical Surface("Fluid", 4) = {1};
