#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"

#define AMBIENT 0.2
#define MAX_DIFFUSE 0.5
#define SPECPOW 50

#define hither 2
#define yon 100
#define REF_DEPTH 4 // reflection depth
#define H 0.57735 //tan(30) half angle

double obmat[100][4][4] ;
double obinv[100][4][4] ;
double color[100][3] ;
int    num_objects ;
int    obtype[100] ; // 0 line segment, 1 circle, 2 hyperbola

double reflectivity[100];
double light_in_eye_space[3];
double light_in_world_space[3];

double ARGB[3];

double T;
int OB;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int Light_Model(double irgb[3], double s[3], double p[3], double n[3], double argb[3]) {

	double len ;
	double N[3] ; 
	len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
	if (len == 0) { return 0 ; }
	N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

	double E[3] ;
	E[0] = s[0] - p[0] ; 
	E[1] = s[1] - p[1] ; 
	E[2] = s[2] - p[2] ; 
	len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
	if (len == 0) { return 0 ; }
	E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
	double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

	double L[3] ;
	L[0] = light_in_eye_space[0] - p[0] ; 
	L[1] = light_in_eye_space[1] - p[1] ; 
	L[2] = light_in_eye_space[2] - p[2] ; 
	len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
	if (len == 0) { return 0 ; }
	L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
	double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;


	double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
	// this needs to occur BEFORE you possibly jump to LLL below



	double intensity ;
	if (NdotL*NdotE < 0) {
		// eye and light are on opposite sides of polygon
		intensity = AMBIENT ; 
		goto LLL ;
	} else if ((NdotL < 0) && (NdotE < 0)) {
		// eye and light on same side but normal pointing "wrong" way
		N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
		NdotL *= (-1.0) ;
	}


	// ignore Blinn's variant
	double R[3] ; // Reflection vector of incoming light
	R[0] = 2*NdotL*N[0] - L[0] ;
	R[1] = 2*NdotL*N[1] - L[1] ;
	R[2] = 2*NdotL*N[2] - L[2] ;

	double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

	double diffuse ;
	if (NdotL <= 0.0) { 
		diffuse = 0.0 ; 
	} else { 
		diffuse = MAX_DIFFUSE*NdotL ; 
	}

	double specular ;
	if (EdotR <= 0.0) { 
		specular = 0.0 ; 
	} else { 
		specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;
	}

	// printf("%lf %lf\n",diffuse,specular) ;
	intensity = AMBIENT + diffuse + specular ;



	LLL : ;

	double f,g ;
	if (intensity <= max_ambient_and_diffuse) {
		f = intensity / max_ambient_and_diffuse ;
		argb[0] = f * irgb[0] ;
		argb[1] = f * irgb[1] ;
		argb[2] = f * irgb[2] ;
	} else {
		f = (intensity - max_ambient_and_diffuse) / 
				   (1.0 - max_ambient_and_diffuse) ;
		g = 1.0 - f ;
		argb[0] = g * irgb[0] + f ;
		argb[1] = g * irgb[1] + f ;
		argb[2] = g * irgb[2] + f ;
	}

	return 1 ;
}

// irgb == inherent color of object (input to this function)
// p0, p1, p2 points in the object
// argb == actual color of object (output of this function)
void light_model (double irgb[3], double p0[3], double p1[3], double p2[3], double argb[3]) {
	double Eye[3] ;
	Eye[0] = 0 ; Eye[1] = 0 ; Eye[2] = 0 ; 

	double P[3]  ;
	P[0] = p0[0] ;  P[1] = p0[1] ;  P[2] = p0[2] ;

	double a[3] ;
	a[0] = p1[0] - p0[0] ;  a[1] = p1[1] - p0[1] ;  a[2] = p1[2] - p0[2] ;

	double b[3] ;
	b[0] = p2[0] - p0[0] ;  b[1] = p2[1] - p0[1] ;  b[2] = p2[2] - p0[2] ;

	double N[3] ;
	M3d_x_product (N, a,b) ;

	Light_Model (irgb, Eye, P, N, argb) ;
}



double magnitude(double v[3]) {

	double mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) ;

	return mag ;

}


int find_reflection(double intersection[3], double source[3], double N[3], double R[3]) {
	// L ray vector in obj space, N normal in obj space, R reflection in obj space
	// L flipped direction

	double L[3] ;

	L[0] = -intersection[0] + source[0] ;
	L[1] = -intersection[1] + source[1] ;
	L[2] = -intersection[2] + source[2] ;

	// unit vector
	double mag = magnitude(L) ;
	L[0] /= mag; L[1] /= mag; L[2] /= mag ;

	double NdotL = L[0] * N[0] + L[1] * N[1] + L[2] * N[2] ;

	// R unit vector

	R[0] = 2*NdotL*N[0] - L[0] ;
	R[1] = 2*NdotL*N[1] - L[1] ;
	R[2] = 2*NdotL*N[2] - L[2] ;

	return 1 ;

}



int find_normal(double intersect[3], double normal[3], double M[4][4], int obnum) {

	// pass through inverse obj matrix; swap along diagonal /

	////////////////////////////////////////////////////////////

	if(obtype[obnum] == 1) {
		// derivative = 2x, 2y
		normal[0] = 2*intersect[0] ;
		normal[1] = 2*intersect[1] ;
		normal[2] = 2*intersect[2] ;


		double nx = M[0][0]*normal[0] + M[1][0]*normal[1] + M[2][0]*normal[2] ;
		double ny = M[0][1]*normal[0] + M[1][1]*normal[1] + M[2][1]*normal[2] ;
		double nz = M[0][2]*normal[0] + M[1][2]*normal[1] + M[2][2]*normal[2] ;

		normal[0] = nx ;
		normal[1] = ny ;
		normal[2] = nz ;

		// unit normal
		double mag = magnitude(normal) ;
		normal[0] /= mag ;
		normal[1] /= mag ;
		normal[2] /= mag ;

		return 1 ;

	}

	///////////////////////////////////////////////////////////////////////

	else if(obtype[obnum] == 0) {
		// derivative = 2x, 2y
		normal[0] = 2*intersect[0] ;
		normal[1] = 0 ;
		normal[2] = 2*intersect[2] ;


		double nx = M[0][0]*normal[0] + M[1][0]*normal[1] + M[2][0]*normal[2] ;
		double ny = M[0][1]*normal[0] + M[1][1]*normal[1] + M[2][1]*normal[2] ;
		double nz = M[0][2]*normal[0] + M[1][2]*normal[1] + M[2][2]*normal[2] ;

		normal[0] = nx ;
		normal[1] = ny ;
		normal[2] = nz ;

		// unit normal
		double mag = magnitude(normal) ;
		normal[0] /= mag ;
		normal[1] /= mag ;
		normal[2] /= mag ;

		return 1 ;

	}

}



int find_intersection(double source[3], double tip[3], double intersect[3], int obnum) {
	// returns 1 if intersection is found; returns 0 if none

	double x0, y0, z0;
	double a, b, c;
	double t1, t2, t;

	double point_y1, point_y2, point_x;


	x0 = tip[0] - source[0];
	y0 = tip[1] - source[1];
	z0 = tip[2] - source[2];

	////////////////////////////////////////////////////////////////////////////////////

	if(obtype[obnum] == 1) {
		// quadratic equation for t
		a = x0*x0 + y0*y0 + z0*z0 ;
		b = 2*x0*source[0] + 2*y0*source[1] + 2*z0*source[2] ;
		c = source[0]*source[0] + source[1]*source[1] + source[2]*source[2] - 1 ;


		// if inside sqrt is negative return 0 (no intersections)

		if(b*b - 4*a*c > 0) {

			//    printf("checking\n");

			//printf("point found\n");

			t1 = (-b + sqrt(b*b - 4*a*c))/(2*a) ;
			t2 = (-b - sqrt(b*b - 4*a*c))/(2*a) ;

			// smallest t = closest intersection pt

			if(t1 <= t2) {
				t = t1 ;
			} else {
				t = t2 ;
			}

			//    printf("t = %lf\n", t);


			if(t < T && t > 0) {

				T = t ;
				OB = obnum ;

				//    printf("intersection found\n");

				intersect[0] = source[0] + x0*T ;
				intersect[1] = source[1] + y0*T ;
				intersect[2] = source[2] + z0*t ;
			}

		} // end if/else

	}

	/////////////////////////////////////////////////////////////////////////////////////////


	else if(obtype[obnum] == 0) {

		double point_y1, point_y2 ;
		// quadratic equation for t
		a = x0*x0 + z0*z0 ;
		b = 2*x0*source[0] + 2*z0*source[2] ;
		c = source[0]*source[0]+ source[2]*source[2] - 1 ;


		// if inside sqrt is negative return 0 (no intersections)

		if(b*b - 4*a*c > 0) {

			//    printf("checking\n");

			//printf("point found\n");

			t1 = (-b + sqrt(b*b - 4*a*c))/(2*a) ;
			t2 = (-b - sqrt(b*b - 4*a*c))/(2*a) ;

			// smallest t = closest intersection pt

			point_y1 = source[1] + t1*y0 ;
			point_y2 = source[1] + t2*y0 ;


			/* if(point_y1 <= 1 && point_y1 >= -1) {
			if(point_y2 <= 1 && point_y2 >= -1) {
			if(t1 <= t2 && t1 > 0)
			{t = t1;}
			else if (t2 > 0) {t = t2;}
			else if (t1 > 0) {t = t1;}
			}
			else if (t1 > 0) {
			t = t1;
			}
			}

			else if(point_y2 <= 1 && point_y2 >= -1 && t2 > 0) {
			t = t2;
			} */

			if(t1 <= t2) {
				t = t1 ;
			} else {
				t = t2 ;
			}

			//    printf("t = %lf\n", t);




			if(t < T && t > 0) {

				T = t ;
				OB = obnum ;

				//    printf("intersection found\n");

				intersect[0] = source[0] + x0*T ;
				intersect[1] = source[1] + y0*T ;
				intersect[2] = source[2] + z0*t ;
			} 
		} // end if b

	} // end if/else


	if (T== 100000000) {
		return 0 ;
	} else {
		return 1 ;
	}
}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int test01() {
	double vm[4][4], vi[4][4] ;
	double Tvlist[100] ;
	int Tn, Ttypelist[100] ;
	double m[4][4], mi[4][4] ;
	double Rsource[3] ;
	double Rtip[3] ;

	light_in_eye_space[0] = 5 ;
	light_in_eye_space[1] = 10 ;
	light_in_eye_space[2] = -4 ;

	double irgb[100][3] ;
	double argb[100][3] ;
	double lrgb[100][3] ;


	//////////////////////////////////////////////////////////////////////
	M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
	//////////////////////////////////////////////////////////////////////

	num_objects = 0 ;

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this into a floor plane like the ocean that is very reflective

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0.9 ;
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   4   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.1   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   4   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   5   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =   0   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  -1   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Unreflective moon

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0.0 ;
	irgb[num_objects][0] = 1;
	irgb[num_objects][1] = 0;
	irgb[num_objects][2] = 0; // red

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.5   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.5   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.5   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =   0.0   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =   0.0   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////


	G_rgb(0,0,0) ;
	G_clear() ;

	Rsource[0] =  0 ;  Rsource[1] =  0 ;  Rsource[2] = 0 ;
	Rtip[0] = -H ; Rtip[1] = -H ; Rtip[2] = 1 ;

	int obnum, insct, x, y, i ;
	double Rsource_new[3], Rtip_new[3], intersect[3] ;
	double ytip_save, xtip_save  ;

	ytip_save = Rtip[0] ;
	xtip_save = Rtip[0] ;

	double normal[3] ;
	double reflection[3] ;

	int objs_reflected[100] ;
	int num_intersections ;

	double proj_u, proj_v, u, v ;
	int id ;
	int dim[2] ;
	double width, height ;


	id = init_xwd_map_from_file("moonJ.xwd") ;
	if(id == -1) {
		printf("failure: can't init map\n") ;
	}
	i = get_xwd_map_dimensions(id, dim) ;
	if(i == -1) {
		printf("failure: can't dimensions\n") ;
	}
	width = dim[0]; height = dim[1] ;


	// double for loop moves through each pt of "film"


	for(x = -400; x < 400; x++) {
		for(y = -400; y < 400; y++) {
			Rtip[0] = x*(H/400);
			Rtip[1] = y*(H/400);
			Rtip[2] = 1;

			Rsource[0] =  0 ;  Rsource[1] =  0 ;  Rsource[2] = 0 ;

			insct = 1;
			num_intersections = 0;

			while (insct == 1 && num_intersections < REF_DEPTH) {

				T = 100000000;
				OB = -1;

				for (obnum = 0; obnum < num_objects; obnum++) {
					M3d_mat_mult_pt(Rsource_new, obinv[obnum], Rsource);
					M3d_mat_mult_pt(Rtip_new, obinv[obnum], Rtip);

					insct = find_intersection(Rsource_new, Rtip_new, intersect, obnum);

				}

				if(OB > -1) {

					// find unit normal in obj space
					i = find_normal(intersect, normal, obinv[OB], OB);

					// texture mapping on sphere 1


					if(OB == 1) { // was at 0, which made the floor plan textured
						u = 0.5 + (atan2(intersect[2], intersect[0]))/(2*M_PI);
						v = 0.5 + asin(intersect[1])/M_PI;

						proj_u = u*width;
						proj_v = v*height;

						i = get_xwd_map_color(id, proj_u, proj_v, irgb[OB]);
					}


					// send intersection to obj space and save it for later
					M3d_mat_mult_pt(intersect, obmat[OB], intersect);

					objs_reflected[num_intersections] = OB;
					num_intersections++;


					// find unit reflection in obj space
					i = find_reflection(intersect, Rsource, normal, reflection);


					// light model
					i = Light_Model(irgb[OB], Rsource, intersect, normal, lrgb[OB]);


					// make intersection pt the new source; reflection pt the new tip
					Rsource[0] = intersect[0] + 0.0001 * reflection[0];
					Rsource[1] = intersect[1] + 0.0001 * reflection[1];
					Rsource[2] = intersect[2] + 0.0001 * reflection[2];
					Rtip[0]    = intersect[0] + reflection[0];
					Rtip[1]    = intersect[1] + reflection[1];
					Rtip[2]    = intersect[2] + reflection[2];


				}

			} // end while

			if(num_intersections > 0) {
				obnum = objs_reflected[num_intersections - 1];
				ARGB[0] = lrgb[obnum][0];
				ARGB[1] = lrgb[obnum][1];
				ARGB[2] = lrgb[obnum][2];

				for(i = num_intersections - 1; i >= 0; i--){
					obnum = objs_reflected[i];

					ARGB[0] *= reflectivity[obnum];
					ARGB[0] += (1 - reflectivity[obnum])*lrgb[obnum][0];

					ARGB[1] *= reflectivity[obnum];
					ARGB[1] += (1 - reflectivity[obnum])*lrgb[obnum][1];

					ARGB[2] *= reflectivity[obnum];
					ARGB[2] += (1 - reflectivity[obnum])*lrgb[obnum][2];

				}


			} else {
				ARGB[0] = 0; ARGB[1] = 0; ARGB[2] = 0;
			}

			G_rgb(ARGB[0], ARGB[1], ARGB[2]);
			G_point(x + 400, y + 400) ;



		} // end for y

		//      printf("\n") ;
	} // end for x


	G_wait_key();
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


int main() {
	G_init_graphics(800,800);
	test01() ;
}
