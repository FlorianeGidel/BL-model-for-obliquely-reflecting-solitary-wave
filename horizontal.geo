ls = 1.0;

Point(1) = {0.0, 0.0, 0.0, ls};
Point(2) = {5.0, 0.0, 0.0, ls};
Point(3) = {427.620397047, 244.0, 0.0, ls};
Point(4) = {427.620397047, 249.4, 0.0, ls};
Point(5) = {427.620397047, 250.0, 0.0, ls};
Point(6) = {5.0, 250.0, 0.0, ls};

Point(7) = {0.0, 250.0, 0.0, ls};

Point(8) = {0.0, 5.4, 0.0, ls};

Point(9) = {5.0, 5.4, 0.0, ls};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {6,9};

Line(10) = {9,2};

Line(11) = {8,9};

Line(12) = {9,4};

Line Loop(13) = {6, 7, 11, -9};
Plane Surface(14) = {13};

Line Loop(15) = {9, 12, 4, 5};
Plane Surface(16) = {15};

Line Loop(17) = {8, 1, -10, -11};
Plane Surface(18) = {17};

Line Loop(19) = {12, -3, -2, -10};
Plane Surface(20) = {19};

Transfinite Line{1, 6, 11} = 12;
Transfinite Line{2, 12, 5} = 1952;
Transfinite Line{3, 8, 10} = 21;
Transfinite Line{4, 7, 9} = 163;
Transfinite Surface{14, 16, 18, 20};
Recombine Surface{14, 16, 18, 20};

