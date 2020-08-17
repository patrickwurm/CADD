xmin = 0;

x0 = xmin;

Geometry.CopyMeshingMethod = 1;

//numy: number of atomic cells in y-direction
// some choices for numx and numy lead to "less structured meshes"
numy = 23; //must be an even number
divy = numy;

a = 2.6491000*1e-10;

dy = 1/2*Sqrt(3)*a;
dx = 1*a;

//dx = 1.03*dx; //250 K PRESTRAIN
//dx = 1.0165*dx; //150 K PRESTRAIN
//dx = 1.0085*dx; //50 K PRESTRAIN
//dy = (1-(dx-a)*1/3)*dy;

dist = 110*dx;

ymin = -(3*Sqrt(3)*a) ;

dlx = 0;
dly = 0;

lc1 = 10e5;

ymax = ymin+numy*Sqrt(3)*a;

//Lx: length of the body in x-direction
//Ly: length of the body in y-direction
//Lx = 160*10^-10;
//Ly = 130*10^-10;
Lx = -610*10^-10;

Lx = -75*dx;

P1 = 110000;
P2 = 120000;
P1R = 130000;
P2R = 140000;

L1 = 210000;
L2 = 220000;
L3 = 230000;
L1R = 240000;
L2R = 250000;
L3R = 260000;

actual_rings=-1;

For ring In {0:10000:1}
	actual_rings = actual_rings+1;

	If( xmin > -50*dx)
	    //dlx = 1.15*dlx;
	    //divy = 3;
        ymax = ymax + dlx;
        ymin = ymin - dlx;
	EndIf



	
	If (xmin-dlx-Lx<dlx)
		Dx = xmin-Lx;
		xmin = xmin-Dx;
		ring = 10000;
	Else
		xmin = xmin-dlx;
	EndIf
	
	
	Point(P1) = {xmin, ymax, 0, lc1};
	Point(P2) = {xmin, ymin, 0, lc1};	
	Point(P1R) = {dist + Fabs(xmin), ymax, 0, lc1};
	Point(P2R) = {dist + Fabs(xmin), ymin, 0, lc1};	
	
	If (ring==1)
        	dlx = 0.999999995*dlx;
        	//numx = 2*numx;
        	//numy = 2*numy;
	EndIf
	
	If (ring==0)
		Line(L1) = {P1, P2};
		Transfinite Line(L1) = numy+1;
		Line(L1R) = {P1R, P2R};
		Transfinite Line(L1R) = numy+1;
		dlx=dx;
	Else
		Line(L1) = {P1, P2};
		P1m = P1-1;
		P2m = P2-1;
		Line(L2) = {P1m, P1};
		Line(L3) = {P2, P2m};
		
		Line(L1R) = {P1R, P2R};
		P1mR = P1R-1;
		P2mR = P2R-1;
		Line(L2R) = {P1mR, P1R};
		Line(L3R) = {P2R, P2mR};
		
		modu = Fmod(ring-1,3);
		If (modu==0)
			If (ring==1)
			Else
				//divx = divx-3;
				//divy = divy-2;
				//divx = Floor((divx-1)/2)+2;
				divx = 1;
			EndIf
		Else
			divx = divx + 2;
		EndIf
		
		If (ring==1)
			//divx = Floor(numx/2)+2;
			divx = 1;
		EndIf
		
		If (ring==15)
			divy = 4;
		EndIf

		If (ring==10)
			divy = 8;
		EndIf

		If (ring==6)
			divy = 12;
		EndIf

		If (ring>4)
			dlx = 1.35*dlx;
		EndIf
	
		//If (ring<6 && ring > 0)
		//	divy = divy+2;
		//EndIf
		
		Transfinite Line(L1) = divy+1;
		Transfinite Line(L2) = 1;
		Transfinite Line(L3) = 1;
		Transfinite Line(L1R) = divy+1;
		Transfinite Line(L2R) = 1;
		Transfinite Line(L3R) = 1;
		
		modu = Fmod(ring-1,3);
		Printf("ring = %f", ring);
		Printf("fmod = %f", modu);
		If (modu==0)
			If (ring==1)
			Else
			    //If(ring>5)
			    If(xmin < (-1)*15*a && xmin > -50*a)
				    dlx=1.15*dlx;
				    If(dlx > 10*dx)
				        dlx=10*dx;
				    EndIf
				EndIf
			EndIf
		EndIf

	EndIf

	If (ring!=10000)
		P1 = P1+1;
		P2 = P2+1;
	
		L1 = L1+1;
		L2 = L2+1;
		L3 = L3+1;
		
		P1R = P1R+1;
		P2R = P2R+1;
	
		L1R = L1R+1;
		L2R = L2R+1;
		L3R = L3R+1;
	EndIf
EndFor

Mesh.Algorithm = 1; //1 = MeshAdapt Algorithm

lineloop = {};
bottom = {};
top = {};
left = {};
left = L1;

lineloopR = {};
bottomR = {};
topR = {};
leftR = {};
leftR = L1R;

lineloop[0] = 210000;
lineloopR[0] = 240000;
Printf("lineloop_test = %f", lineloop[0]);

L3=230001;
L3R=260001;
For ring In {1:actual_rings:1}
	lineloop[ring] = -L3;
	bottom[ring-1] = L3;
	L3++; 
	lineloopR[ring] = -L3R;
	bottomR[ring-1] = L3R;
	L3R++; 
EndFor

ring = actual_rings;
lineloop[ring+1] = -L1;
lineloopR[ring+1] = -L1R;
Printf("lineloop_test = %f", lineloop[ring+1]);

index1 = -1;
index = ring+1;
For ring In {actual_rings:1:-1}
	index=index+1;
	index1 = index1+1;
	lineloop[index] = -L2;
	top[index1] = L2;
	L2--; 
	lineloopR[index] = -L2R;
	topR[index1] = L2R;
	L2R--; 
EndFor


Line Loop(1) = lineloop[];
Line Loop(2) = lineloopR[];
Plane Surface(1) = {1};
Plane Surface(2) = {2};

L1 = 210001;
L1R = 240001;
For ring In {1:actual_rings:1}
	Line{L1} In Surface{1};
	L1++;
	Line{L1R} In Surface{2};
	L1R++;
EndFor

Printf("top = %f", top[0]);
Printf("top = %f", top[1]);
Printf("top = %f", top[2]);
Printf("top = %f", top[3]);
Printf("top = %f", top[4]);
Printf("top = %f", top[5]);
Printf("top = %f", top[6]);
Printf("top = %f", top[7]);
Printf("top = %f", top[8]);
Printf("top = %f", top[9]);
Printf("top = %f", top[10]);
Printf("top = %f", top[11]);

Mesh.Algorithm = 6; //1 = MeshAdapt Algorithm, 6 = Frontal Algorithm
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthMin = 1e20;

Physical Line ("L: Atoms right") = lineloop[0];
Physical Line ("L: Left") = left[];
//Physical Line ("L: Bottom") = bottom[];
//Physical Line ("L: Top") = top[];
Physical Line ("L: Bottom") = bottom[0];
Physical Line ("L: Top") = top[0];


Physical Line ("R: Atoms right") = lineloopR[0];
Physical Line ("R: Left") = leftR[];
//Physical Line ("R: Bottom") = bottomR[];
//Physical Line ("R: Top") = topR[];
Physical Line ("R: Bottom") = bottomR[0];
Physical Line ("R: Top") = topR[0];

Physical Surface("Mesh") = {1, 2}; 




