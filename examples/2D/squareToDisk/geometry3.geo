//+
d = 0.05;
//+
Point(1) = {0.1, 0, 0, d};
//+
Point(2) = {0, 0.1, 0, d};
//+
Point(3) = {0, 0.9, 0, d};
//+
Point(4) = {0.1, 1, 0, d};
//+
Point(5) = {0.9, 0, 0, d};
//+
Point(6) = {1, 0.1, 0, d};
//+
Point(7) = {0.9, 1, 0, d};
//+
Point(8) = {1, 0.9, 0, d};
//+
Point(9) = {0.9, 0.9, 0, d};
//+
Point(10) = {0.1, 0.1, 0, d};
//+
Point(11) = {0.9, 0.1, 0, d};
//+
Point(12) = {0.1, 0.9, 0, d};
//+
Line(1) = {3, 2};
//+
Line(2) = {1, 5};
//+
Line(3) = {6, 8};
//+
Line(4) = {7, 4};
//+
Circle(5) = {3, 12, 4};
//+
Circle(6) = {8, 9, 7};
//+
Circle(7) = {6, 11, 5};
//+
Circle(8) = {1, 10, 2};
//+
Line(9) = {3, 12};
//+
Line(10) = {12, 4};
//+
Line(11) = {12, 9};
//+
Line(12) = {9, 7};
//+
Line(13) = {9, 8};
//+
Line(14) = {9, 11};
//+
Line(15) = {11, 10};
//+
Line(16) = {10, 12};
//+
Line(17) = {2, 10};
//+
Line(18) = {1, 10};
//+
Line(19) = {5, 11};
//+
Line(20) = {11, 6};
//+
Curve Loop(2) = {16, 11, 14, 15};
//+
Plane Surface(1) = {2};
//+
Curve Loop(3) = {18, -15, -19, -2};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {20, 7, 19};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {3, -13, 14, 20};
//+
Plane Surface(4) = {5};
//+
Curve Loop(6) = {12, 4, -10, 11};
//+
Plane Surface(5) = {6};
//+
Curve Loop(7) = {5, -10, -9};
//+
Plane Surface(6) = {7};
//+
Curve Loop(8) = {12, -6, -13};
//+
Plane Surface(7) = {8};
//+
Curve Loop(9) = {1, 17, 16, -9};
//+
Plane Surface(8) = {9};
//+
Curve Loop(10) = {18, -17, -8};
//+
Plane Surface(9) = {10};
//+
Physical Curve("FreeSurface", 21) = {1, 5, 4, 6, 3, 2, 7, 8};
//+
Physical Surface("Fluid", 22) = {8, 1, 5, 6, 4, 7, 2, 9, 3};
