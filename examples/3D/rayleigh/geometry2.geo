//+
d = 0.1;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, 4, 0, d};
//+
Point(3) = {4, 0, 0, d};
//+
Point(4) = {4, 4, 0, d};
//+
Point(5) = {0, 0, 1, d};
//+
Point(6) = {0, 4, 1, d};
//+
Point(7) = {4, 0, 1, d};
//+
Point(8) = {4, 4, 1, d};
//+
Point(9) = {0, 0, 1.5, d};
//+
Point(10) = {0, 4, 1.5, d};
//+
Point(11) = {4, 0, 1.5, d};
//+
Point(12) = {4, 4, 1.5, d};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 2};
//+
Line(5) = {6, 5};
//+
Line(6) = {5, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 6};
//+
Line(9) = {6, 2};
//+
Line(10) = {5, 1};
//+
Line(11) = {7, 3};
//+
Line(12) = {8, 4};
//+
Curve Loop(1) = {5, 10, -1, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -11, -6, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 3, -12, -7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 9, -4, -12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 6, 7, 8};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 5, 2, 6, 3, 4};
//+
Volume(1) = {1};

//+
Line(13) = {5, 9};
//+
Line(14) = {6, 10};
//+
Line(15) = {7, 11};
//+
Line(16) = {8, 12};
//+
Line(17) = {11, 9};
//+
Line(18) = {9, 10};
//+
Line(19) = {10, 12};
//+
Line(20) = {12, 11};
//+
Curve Loop(7) = {17, -13, 6, 15};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {7, 16, 20, -15};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {16, -19, -14, -8};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {14, -18, -13, -5};
//+
Plane Surface(10) = {10};
//+
Physical Surface("Bottom", 21) = {6};
//+
Physical Surface("Lateral", 22) = {1, 10, 2, 7, 3, 8, 4, 9};
//+
Physical Surface("FreeSurface", 23) = {5};
//+
Physical Volume("Fluid", 24) = {1};
