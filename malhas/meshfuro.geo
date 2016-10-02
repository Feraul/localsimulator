//--------------------------------------------------------------------
//---------------- UNIVERSIDADE FEDERAL DE PERNAMBUCO ----------------
//---------------- CENTRO DE TECNOLOGIA E GEOCIENCIAS ----------------
//---------- PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL -----------
//--------------------------------------------------------------------

//Work developed by: Marcio Souza and Luiz E. Queiroz
//Adviser Professors: Paulo Lyra & Darlan Carvalho
//Create date: 2015/8/31;	hour: 21:15h

//--------------------------------------------------------------------
//This file has CAD parameters. It is related to building of domain

//"cl1" corresponds to element size attributed in "Start.dat";
cl1 = 0.125000;

Point(1) = {0.000000, 0.000000, 0.000000, cl1};
Point(2) = {1.000000, 0.000000, 0.000000, cl1};
Point(3) = {1.000000, 1.000000, 0.000000, cl1};
Point(4) = {0.000000, 1.000000, 0.000000, cl1};
Point(5) = {0.444400, 0.444400, 0.000000, cl1};
Point(6) = {0.555500, 0.444400, 0.000000, cl1};
Point(7) = {0.555500, 0.555500, 0.000000, cl1};
Point(8) = {0.444400, 0.555500, 0.000000, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line Loop(11) = {1, 2, 3, 4, -8, -7, -6, -5};
Plane Surface(11) = {11};
Physical Point(101) = {1, 2, 3, 4};
Physical Point(102) = {5, 6, 7, 8};
Physical Line(101) = {1, 2, 3, 4};
Physical Line(102) = {5, 6, 7, 8};
Physical Surface(1) = {11};

