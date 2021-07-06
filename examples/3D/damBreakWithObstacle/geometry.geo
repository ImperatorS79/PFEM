//+
a = 1.228;
//+
b = 1.248;
//+
c = 0.744;
//+
d = 0.4;
//+
e = 1;
//+
s = 0.16;
//+
h = 0.55;
//+
h2 = 0.75;
//+
de = 0.1;
//+
Point(1) = {0, -e/2, 0, de};
//+
Point(2) = {0, e/2, 0, de};
//+
Point(3) = {a, -e/2, 0, de};
//+
Point(4) = {a, e/2, 0, de};
//+
Point(5) = {a + b + c, -e/2, 0, de};
//+
Point(6) = {a + b + c, e/2, 0, de};
//+
Point(7) = {0, -e/2, h, de};
//+
Point(8) = {0, e/2, h, de};
//+
Point(9) = {a, -e/2, h, de};
//+
Point(10) = {a, e/2, h, de};
//+
Point(11) = {0, -e/2, h2, de};
//+
Point(12) = {0, e/2, h2, de};
//+
Point(13) = {a + b + c, -e/2, h2, de};
//+
Point(14) = {a + b + c, e/2, h2, de};
//+
Point(15) = {a + b - s/2, d/2, 0, de};
//+
Point(16) = {a + b - s/2, -d/2, 0, de};
//+
Point(17) = {a + b + s/2, d/2, 0, de};
//+
Point(18) = {a + b + s/2, -d/2, 0, de};
//+
Point(19) = {a + b - s/2, d/2, s, de};
//+
Point(20) = {a + b - s/2, -d/2, s, de};
//+
Point(21) = {a + b + s/2, d/2, s, de};
//+
Point(22) = {a + b + s/2, -d/2, s, de};
//+
Point(23) = {a, -e/2, h2, de};
//+
Point(24) = {a, e/2, h2, de};
//+
Point(25) = {a + b - s/2, e/2, 0, de};
//+
Point(26) = {a + b - s/2, -e/2, 0, de};
//+
Point(27) = {a + b + s/2, e/2, 0, de};
//+
Point(28) = {a + b + s/2, -e/2, 0, de};
//+
Point(29) = {a + b + c, -e/2, h, de};
//+
Point(30) = {a + b + c, e/2, h, de};
//+
Point(31) = {a + b + c, d/2, 0, de};
//+
Point(32) = {a + b + c, -d/2, 0, de};
//+
Point(33) = {a, d/2, 0, de};
//+
Point(34) = {a, -d/2, 0, de};
//+
Point(35) = {0, d/2, 0, de};
//+
Point(36) = {0, -d/2, 0, de};
//+
Point(37) = {a, d/2, h, de};
//+
Point(38) = {a, -d/2, h, de};
//+
Point(39) = {0, d/2, h, de};
//+
Point(40) = {0, -d/2, h, de};
//+
Point(41) = {a + b - s/2, e/2, h2, de};
//+
Point(42) = {a + b - s/2, -e/2, h2, de};
//+
Point(43) = {a + b + s/2, e/2, h2, de};
//+
Point(44) = {a + b + s/2, -e/2, h2, de};
//+
Point(45) = {a + b - s/2, e/2, h, de};
//+
Point(46) = {a + b - s/2, -e/2, h, de};
//+
Point(47) = {a + b + s/2, e/2, h, de};
//+
Point(48) = {a + b + s/2, -e/2, h, de};
//+
Point(49) = {a + b + c, d/2, h, de};
//+
Point(50) = {a + b + c, -d/2, h, de};
//+
Point(51) = {a + b + c, d/2, h2, de};
//+
Point(52) = {a + b + c, -d/2, h2, de};
//+
Point(53) = {0, d/2, h2, de};
//+
Point(54) = {0, -d/2, h2, de};
//+
Line(1) = {11, 7};
//+
Line(2) = {7, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 9};
//+
Line(5) = {9, 7};
//+
Line(6) = {23, 9};
//+
Line(7) = {3, 34};
//+
Line(8) = {34, 33};
//+
Line(9) = {33, 4};
//+
Line(10) = {4, 10};
//+
Line(11) = {10, 24};
//+
Line(12) = {10, 37};
//+
Line(13) = {37, 38};
//+
Line(14) = {38, 9};
//+
Line(15) = {1, 36};
//+
Line(16) = {36, 35};
//+
Line(17) = {35, 2};
//+
Line(18) = {2, 4};
//+
Line(19) = {12, 8};
//+
Line(20) = {2, 8};
//+
Line(21) = {8, 10};
//+
Line(22) = {8, 39};
//+
Line(23) = {39, 40};
//+
Line(24) = {40, 7};
//+
Line(25) = {36, 34};
//+
Line(26) = {33, 35};
//+
Line(27) = {40, 38};
//+
Line(28) = {37, 39};
//+
Line(30) = {11, 23};
//+
Line(31) = {23, 24};
//+
Line(32) = {24, 12};
//+
Line(33) = {3, 26};
//+
Line(34) = {26, 28};
//+
Line(35) = {28, 5};
//+
Line(36) = {5, 32};
//+
Line(37) = {32, 31};
//+
Line(38) = {31, 6};
//+
Line(39) = {6, 27};
//+
Line(40) = {27, 25};
//+
Line(41) = {25, 4};
//+
Line(42) = {34, 16};
//+
Line(43) = {15, 33};
//+
Line(44) = {18, 32};
//+
Line(45) = {31, 17};
//+
Line(46) = {26, 16};
//+
Line(47) = {15, 16};
//+
Line(48) = {15, 25};
//+
Line(49) = {27, 17};
//+
Line(50) = {17, 18};
//+
Line(51) = {18, 28};
//+
Line(52) = {18, 16};
//+
Line(53) = {15, 17};
//+
Line(54) = {17, 21};
//+
Line(55) = {19, 21};
//+
Line(56) = {21, 22};
//+
Line(57) = {22, 20};
//+
Line(58) = {20, 19};
//+
Line(59) = {20, 16};
//+
Line(60) = {18, 22};
//+
Line(61) = {19, 15};
//+
Line(62) = {9, 46};
//+
Line(63) = {46, 48};
//+
Line(64) = {48, 29};
//+
Line(65) = {29, 5};
//+
Line(66) = {46, 26};
//+
Line(67) = {48, 28};
//+
Line(68) = {23, 42};
//+
Line(69) = {42, 44};
//+
Line(70) = {44, 13};
//+
Line(71) = {13, 29};
//+
Line(72) = {42, 46};
//+
Line(73) = {44, 48};
//+
Line(74) = {24, 41};
//+
Line(75) = {41, 43};
//+
Line(76) = {43, 14};
//+
Line(77) = {14, 30};
//+
Line(78) = {30, 6};
//+
Line(79) = {10, 45};
//+
Line(80) = {45, 47};
//+
Line(81) = {30, 47};
//+
Line(82) = {41, 45};
//+
Line(83) = {45, 25};
//+
Line(84) = {27, 47};
//+
Line(85) = {47, 43};
//+
Line(86) = {29, 50};
//+
Line(87) = {50, 49};
//+
Line(88) = {49, 30};
//+
Line(89) = {13, 52};
//+
Line(90) = {52, 51};
//+
Line(91) = {51, 14};
//+
Line(92) = {52, 50};
//+
Line(93) = {50, 32};
//+
Line(94) = {31, 49};
//+
Line(95) = {49, 51};
//+
Line(96) = {34, 38};
//+
Line(97) = {33, 37};
//+
Line(98) = {36, 40};
//+
Line(99) = {39, 35};
//+
Curve Loop(1) = {2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14, -4, 7, 96};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {13, -96, 8, 97};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, -97, 9, 10};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {15, 25, -7, -3};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {8, 26, -16, 25};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {26, 17, 18, -9};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {24, 2, 15, 98};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {16, -99, 23, -98};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {17, 20, 22, 99};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {20, 21, -10, -18};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {5, -24, 27, 14};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {23, 27, -13, 28};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {12, 28, -22, 21};
//+
Plane Surface(14) = {14};
//+
Line(100) = {12, 53};
//+
Line(101) = {53, 54};
//+
Line(102) = {54, 11};
//+
Line(103) = {40, 54};
//+
Line(104) = {39, 53};
//+
Curve Loop(15) = {102, 1, -24, 103};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {23, 103, -101, -104};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {100, -104, -22, -19};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {1, -5, -6, -30};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {32, 19, 21, 11};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {6, 62, -72, -68};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {72, 63, -73, -69};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {73, 64, -71, -70};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {4, 62, 66, -33};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {66, 34, -67, -63};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {35, -65, -64, 67};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {65, 36, -93, -86};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {86, -92, -89, 71};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {87, -94, -37, -93};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {38, -78, -88, -94};
//+
Plane Surface(29) = {29};
//+
Curve Loop(30) = {81, -84, -39, -78};
//+
Plane Surface(30) = {30};
//+
Curve Loop(31) = {84, -80, 83, -40};
//+
Plane Surface(31) = {31};
//+
Curve Loop(32) = {83, 41, 10, 79};
//+
Plane Surface(32) = {32};
//+
Curve Loop(33) = {74, 82, -79, 11};
//+
Plane Surface(33) = {33};
//+
Curve Loop(34) = {82, 80, 85, -75};
//+
Plane Surface(34) = {34};
//+
Curve Loop(35) = {81, 85, 76, 77};
//+
Plane Surface(35) = {35};
//+
Curve Loop(36) = {33, 46, -42, -7};
//+
Plane Surface(36) = {36};
//+
Curve Loop(37) = {44, -36, -35, -51};
//+
Plane Surface(37) = {37};
//+
Curve Loop(38) = {51, -34, 46, -52};
//+
Plane Surface(38) = {38};
//+
Curve Loop(39) = {44, 37, 45, 50};
//+
Plane Surface(39) = {39};
//+
Curve Loop(40) = {49, -45, 38, 39};
//+
Plane Surface(40) = {40};
//+
Curve Loop(41) = {43, 9, -41, -48};
//+
Plane Surface(41) = {41};
//+
Curve Loop(42) = {49, -53, 48, -40};
//+
Plane Surface(42) = {42};
//+
Curve Loop(43) = {42, -47, 43, -8};
//+
Plane Surface(43) = {43};
//+
Curve Loop(44) = {61, 47, -59, 58};
//+
Plane Surface(44) = {44};
//+
Curve Loop(45) = {52, -59, -57, -60};
//+
Plane Surface(45) = {45};
//+
Curve Loop(46) = {56, -60, -50, 54};
//+
Plane Surface(46) = {46};
//+
Curve Loop(47) = {54, -55, 61, 53};
//+
Plane Surface(47) = {47};
//+
Curve Loop(48) = {55, 56, 57, 58};
//+
Plane Surface(48) = {48};
//+
Surface Loop(1) = {1, 8, 5, 6, 7, 10, 11, 9, 13, 12, 2, 3, 4, 14};
//+
Volume(1) = {1};
//+
Curve Loop(49) = {90, -95, -87, -92};
//+
Plane Surface(49) = {49};
//+
Curve Loop(50) = {91, 77, -88, 95};
//+
Plane Surface(50) = {50};
//+
Physical Volume("Fluid", 105) = {1};
//+
Physical Surface("FreeSurface", 106) = {12, 13, 14, 4, 3, 2};
//+
Physical Surface("Walls", 107) = {1, 18, 8, 9, 10, 17, 16, 15, 11, 19, 7, 6, 5, 20, 32, 23, 33, 34, 35, 30, 21, 22, 24, 31, 25, 39, 40, 37, 38, 42, 28, 26, 27, 49, 29, 50, 48, 45, 47};
//+
Physical Surface("Walls", 107) += {41, 43, 36};
//+
Physical Surface("Walls", 107) += {46, 44};
