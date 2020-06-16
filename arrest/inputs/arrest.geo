//all dimensions are in mm

//major dimensions
a = 140;//length of specimen
b = 70.0;//height of specimen
d1 = 30; x1 = 45; y1 = 35;//hole-1
d2 = 23; x2 = 90; y2 = 50;//hole-2
l = 15; w = 0.5; xl = 75;//initial crack
r1 = d1/2; dy = w/2; dx =(r1*r1-dy^2)^0.5;

e = 40; f = (b-e)/2; g = 1.5;

//---- mesh size, smaller finer meshes ------------
l0 = 0.0375;
h = 0.75; hh = l0/4;


//points
Point(1)={0,0,0,h};
Point(2)={xl,0,0,h};
Point(3)={a,0,0,h};
Point(4)={a,f,0,h};
Point(5)={a,b-f,0,h};
Point(6)={a,b,0,h};
Point(7)={xl,b,0,h};
Point(8)={0,b,0,h};

Point(30)={x1,y1,0,h};
Point(31)={x1+d1/2,y1,0,h};
Point(9)={x1+dx,y1+dy,0,h};
Point(10)={x1+dx,y1-dy,0,h};
Point(11) = {xl, y1+dy,0,hh};
Point(12)= {xl, y1-dy,0,hh};
Point(32)={x1,y1+d1/2,0,h};
Point(33)={x1-d1/2,y1,0,h};
Point(34)={x1,y1-d1/2,0,h};

Point(40)={x2,y2,0,h};
Point(41)={x2+d2/2,y2,0,h};
Point(42)={x2,y2+d2/2,0,h};
Point(43)={x2-d2/2,y2,0,hh*5};
Point(44)={x2,y2-d2/2,0,hh*2};

Point(52)={x1,b,0,h};
Point(53)={0,y1,0,h};
Point(54)={x1,0,0,h};

Point(62)={x2,b,0,h};
Point(63)={xl,y2,0,hh*5};
Point(64)={x2,0,0,h};
Point(641)= {x2, y1+dy,0,hh};
Point(642)= {x2, y1-g,0,hh};
Point(643)= {xl, y1-g,0,hh};

Point(13)={0,f,0,h};
Point(14)={0,b-f,0,h};

//lines
Line(1)={1,54}; Line(101)={54,2}; 
Line(2)={2,64}; Line(201)={64,3};
Line(3)={3,4}; Line(4)={4,5}; Line(5)={5,6}; 
Line(6)={6,62}; Line(601)={62,7};
Line(7)={7,52}; Line(701)={52,8};
Line(801)={8,14}; Line(8)={14,53}; Line(9)={53,13}; Line(901)={13,1};

Line(10)={10,12}; Line(11)={9,11};
Line(12)={2,643}; Line(121)={643,12}; 
Line(13)={12,11}; 
Line(14)={11,63}; Line(15)={63,7};



Line(52)={32,52}; Line(53)={33,53}; Line(54)={34,54};
Line(62)={42,62}; Line(63)={43,63}; 
Line(64)={44,641}; Line(643)={641,642}; Line(641)={642,64}; Line(642)={11,641}; 
Line(644)={643,642};

Circle(31) = {9,30,32};
Circle(32) = {32,30,33};
Circle(33) = {33,30,34};
Circle(34) = {34,30,10};

Circle(41) = {41,40,42};
Circle(42) = {42,40,43};
Circle(43) = {43,40,44};
Circle(44) = {44,40,41};

//line loop
Line Loop(1)={53,9,901,1,-54,-33};
Line Loop(2)={54,101,12,121,-10,-34};
Line Loop(3)={11,14,15,7,-52,-31};
Line Loop(4)={52,701,801,8,-53,-32};
Line Loop(5)={2,-641,-644,-12};
Line Loop(9)={644,-643,-642,-13,-121};
Line Loop(6)={201,3,4,5,6,-62,-41,-44,64,643,641};
Line Loop(7)={601,-15,-63,-42,62};
Line Loop(8)={63,-14,642,-64,-43};

//plane surfaces
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};Plane Surface(4)={4};
Plane Surface(5)={5};Plane Surface(6)={6};Plane Surface(7)={7};Plane Surface(8)={8};
Plane Surface(9)={9};

//for quadrlateral elements
//Recombine Surface{1,2,3,4,5,6,7,8,9};


//physical items
Physical Surface("bulk")={1,2,3,4,5,6,7,8,9};

Physical Line("RightNodes")={4}; //right-fixed
Physical Line("LeftNodes")={8,9}; //left-load

//Physical Line(444)={31,32,33,34}; //CircleEdge1
//Physical Line(555)={41,42,43,44}; //CircleEdge2

//Physical Point(666)={9}; //cracknode-1
//Physical Point(777)={10}; //cracknode-2

