//+
H = 0.02;
//+
W = 0.02;
//+
d = 0.00025;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, H, 0, d};
//+
Point(3) = {W, 0, 0, d};
//+
Point(4) = {W, H, 0, d};
//+
Point(5) = {0, H+0.3*H, 0, d};
//+
Point(6) = {W, H+0.3*H, 0, d};
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
//+
Line(5) = {2, 5};
//+
Line(6) = {4, 6};
//+
Physical Curve("Left", 1) += {5};
//+
Physical Curve("Right", 3) += {6};
//+
Point(7) = {-0.1*W, 1.3*H, 0, d};
//+
Point(8) = {1.1*W, 1.3*H, 0, d};
//+
Line(7) = {5, 7};
//+
Line(8) = {6, 8};
//+
Physical Curve("Left", 1) += {7};
//+
Physical Curve("Right", 3) += {8};
