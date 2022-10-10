/*********************************************************************************************
This is a rough C++ port of KS_circle_pack.c, the Knotscape hyperbolic circle packing algorithm.
Little attempt has been made to tidy or optimize the code, so that comparison with the original 
is easier.  However, only the necessary functions have been ported and the Knotscape data 
structures have been replaced with the same 'flowers', 'vertex_centre', 'radius' and so on 
used in circlepack.cpp.

The code for the read and write phases has been taken from circle_pack.cpp, as has the algorithm 
for selecting the next circle to place into the packing (it produced the same order of placing 
circles as in KS_circle_pack.c but uses much less opaque logic!)  Thus, this port really only 
retains the Newton-Raphson technique of the function h_riffle, together with all the supporting 
hyperbolic geometry, allowing a more direct comparison with circle_pack.cpp.

The functions, mob_norm, mob_norm_inv, and mob_trans duplicate the functionality of the earlier
mobius_transformation, inv_mobius_transformation and mobius_transform_to_origin (in drawfuns.cpp)
but have been retained, for now, as the are very small so represent only a little "code bloat".
**********************************************************************************************/

using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <iomanip>
#include <complex>


/********************* External variables ***********************/
extern ofstream     debug;

#include <util.h>
#include <matrix.h>
#include <draw.h>

#define MAX_SIZE 5001		/* max number of vertices of packing */
#define ITERATIONS 1000000	/* max passes through the repack routine */
#define NUM_PETALS 60		/* max number of neighbors of any circle */
#define MAX_COMPONENTS 1	/* number of components of complex should be simply connected here */
#define TOLER .000001		/* threshold for repack computations */
#define OKERR .00000001		/* Roundoff for zero */

extern double two_pi;

extern char const* circlepack_output_file;
extern char const* triangulation_output_file;

/* geometry control */
extern bool DRAW_IN_HYPERBOLIC_DISC;
extern int SHIFT_CIRCLE_PACKING_ORIGIN;

int MAX_H_RADCALC_ITERATIONS=2;			/* parameter of current packing routines */

/********************* Function prototypes ***********************/
void h_riffle(matrix<int>& flowers, vector<double>& radius, int num_type_1234_vertices, int passes);
double h_radcalc(matrix<int>& flowers, vector<double>& radius, int vertex);
void h_curvcalc(matrix<int>& flowers, vector<double>& radius, int i, double s, double* c, double* d);
double h_anglesum (matrix<int>& flowers, vector<double>& radius, int vertex);
void h_anglesum_deltas (matrix<int>& flowers, vector<double>& radius, vector<double>& delta, int num_type_1234_vertices);
double h_comp_cos(double s1,double s2,double s3);
int h_compcenter(complex<double> z1,complex<double> z2,complex<double>* z3, double s1,double s2,double* s3);
double h_cos_overlap(double s1,double s2,double s3,double t1,double t2,double t3);
double e_cos_overlap(double e1,double e2,double e3,double t1,double t2,double t3);
complex<double> mob_norm(complex<double> z,complex<double> a,complex<double> b); 
complex<double> mob_norm_inv(complex<double> w,complex<double> a,complex<double> b);
complex<double> mob_trans(complex<double>z,complex<double>a);
int circle_3(complex<double> z1,complex<double> z2,complex<double> z3,complex<double>* c,double* rad);
void compute_circle_pack_centres(matrix<int>& flowers, vector<double>& radius, int nodecount, int num_type_1234_vertices, matrix<double>& vertex_centre);
void h_to_e_data(complex<double> h_centre,double s_rad,complex<double>* e_centre,double* e_rad);
void h_norm_pack(matrix<double>& vertex_centre,int alpha,int gamma);


int KS_circle_pack(char const* infile, char const* outfile)
{
	int nodecount;
	int alpha;
	int beta;
	int gamma;

	/* Read triangulation data from triangulation_output_file */
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "AB_circle_pack: read triangulation data from" << triangulation_output_file << endl;

	ifstream input(triangulation_output_file);
	if(!input)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "AB_circle_pack: could not open " << triangulation_output_file << endl;
		exit(0);
	}		
	
	string next_line;
	getline(input,next_line);
	if (next_line.find("NODECOUNT") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> nodecount;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "AB_circle_pack: nodecount = " << nodecount << endl;
	}

	getline(input,next_line);
	if (next_line.find("ALPHA") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> alpha; //always 1, the first vertex
		iss >> beta; //the first boundary (type 4) vertex
		iss >> gamma; // IS THIS USED?
			
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "AB_circle_pack: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
	}
	
	/* store the flowers for each node in a matrix with the first column set to the number of neighbours
	   of that vertex.  No vertex is adjacent to itself, nor is any vertex adjacent to all other vertices
	   so the maximum number of neighbours is <= nodecount -2.  Recall that for interior vertices the first
	   petal appears at the end of the list as well, to simplify the calculation of angle sums.
	*/
	matrix<int> flowers(nodecount,nodecount);

	while (getline(input,next_line))
	{
		if (next_line.find("FLOWERS") != string::npos)
		{
			getline(input,next_line); // read the empty line
			
			/* now read the flowers */
			while(getline(input,next_line))
			{
				if (next_line.find("END") != string::npos)
					break;
				
				istringstream iss(next_line);
				int node;
				int count;
				int first_petal;
				int petal;
				
				iss >> node;
				iss >> count;
				iss >> first_petal;
				node--;
				
				flowers[node][0] = count;
				flowers[node][1] = first_petal;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "AB_circle_pack: node " << node << ", petal count = " << count << ": " << first_petal;
				
				for (int i=0; i< count; i++)
				{
					iss >> petal;
					flowers[node][i+2] = petal;
					
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << " " << petal;					
				}

				if (node < beta-1)
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "; interior vertex" << endl;					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "; boundary vertex" << endl;					
				}										
			}		
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "flowers" << endl;					
	print(flowers,debug,4,"");
}
	
	int num_type_1234_vertices = beta-1;
	
	
	/* set the intial circle radii; we are working with hyperbolic geometry and 
	   for a hyperbolic radius r, the value exp(-r) is stored for finite radii; 
	   for infinite hyp radius, we store the negative of the euclidean radius 
	   (which can then use euclidean data for plotting). */	
	
	vector<double> radius(nodecount);
	for (int i=0; i< num_type_1234_vertices; i++)
		radius[i] = 0.1; //set interior vertex initial radius
	for (int i = num_type_1234_vertices; i< nodecount; i++)
		radius[i] = -0.1; // set boundary vertex radius

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "AB_circle_pack: initial radii ";
	for (int i = 0; i< nodecount; i++)
		debug << radius[i] << ' ';
	debug << endl;
}

	h_riffle(flowers, radius, num_type_1234_vertices,ITERATIONS);


if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "AB_circle_pack: final vertex radii " << endl;
	for (int i = 0; i< nodecount; i++)
		debug << "AB_circle_pack: vertex " << i << " " << radius[i] << endl;
}

	/* compute the circle centres based on the final radii */
	matrix<double> vertex_centre(nodecount,2);
	compute_circle_pack_centres(flowers, radius, nodecount, num_type_1234_vertices,vertex_centre);

	/* normalize */
	h_norm_pack(vertex_centre,alpha-1,gamma-1);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "AB_circle_pack: h_norm_pack produced vertex centres:" << endl;
	print (vertex_centre, debug, 10,"AB_circle_pack: ");

	debug << "AB_circle_pack: final hyperbolic radii " << endl;
	for (int i = 0; i< nodecount; i++)
		debug << "AB_circle_pack: vertex " << i << " " << radius[i] << endl;
}

	
/* ---------- convert data to euclidean ------ */

    if (!DRAW_IN_HYPERBOLIC_DISC)
    {
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "AB_circle_pack: converting to euclidean centres with h_to_e_data" << endl;
		
		for (int i=0;i<nodecount;i++)
		{
			complex<double>h_centre;
			h_centre.real(vertex_centre[i][0]);
			h_centre.imag(vertex_centre[i][1]);
			
			double s_rad = radius[i];
			
			complex<double> e_centre;
			double e_rad;;
			
			h_to_e_data(h_centre,s_rad,&e_centre,&e_rad);

			vertex_centre[i][0] = e_centre.real();
			vertex_centre[i][1] = e_centre.imag();
			radius[i] = e_rad;
		}
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "AB_circle_pack: euclidean vertex centres:" << endl;
	print (vertex_centre, debug, 10,"AB_circle_pack: ");
	
	debug << "AB_circle_pack: euclidean radii " << endl;
	for (int i = 0; i< nodecount; i++)
		debug << "AB_circle_pack: vertex " << i << " " << radius[i] << endl;
}	

	}

 /*----------- write output phase ------------ */

	ofstream output(circlepack_output_file);
	
	if (!output)
	{
		cout << "\nError opening output file\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "AB_circle_pack: could not open " << circlepack_output_file << endl;
		exit(0);
	}		
	
	/* write the centres and radii of the type 4 vertices to the output file */
	output << nodecount << endl;

	output << fixed;

	for (int i= 0; i< nodecount; i++)
	{
		output << setprecision(6) << setw(10) << vertex_centre[i][0] 		
		       << setprecision(6) << setw(10) << vertex_centre[i][1] << endl;		
	}
	for (int i= 0; i< nodecount; i++)
		output << setprecision(6) << setw(10) << radius[i] << endl;		

	output.close();
	return 1;

} /* end of main */


 /* adjust s_radii to meet curvature targets in 'aim'. aim<0 means that radius is fixed - no adjustments.
    Bdry radii adjusted only if their aim >= 0. 
 */
void h_riffle(matrix<int>& flowers, vector<double>& radius, int num_type_1234_vertices, int passes)
{
	int count=0;
	
	/* establish a vector to record the delta between a vertex's angle sum and 2\pi */
	vector<double> delta(num_type_1234_vertices);
	vector<double> angle_sum(num_type_1234_vertices); // recorded for debugging purposes only
	
	h_anglesum_deltas (flowers, radius, delta, num_type_1234_vertices);

	double accum = 0;
	for (int i=0; i< num_type_1234_vertices; i++)
		accum += abs(delta[i]);
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "h_riffle: accum = " << accum << endl;

	/* we have cut = accnum/(3*num_type_1234_vertices) where accnum is the aggregate absolute delta
	   between the initial angle sum and the target angle sum.  The initial angle sum was calculated
	   by readpack from the initial radii.
	*/
	double cut=accum/3/num_type_1234_vertices;

	while (cut > TOLER && count < passes)
	{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "h_riffle: starting iteration with angle sum error threshold = " << cut << endl;

		for (int i=0;i<num_type_1234_vertices;i++)
	    {
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "h_riffle: node " << i << endl;
			double current_angle_sum = h_anglesum (flowers, radius, i);
			angle_sum[i] = current_angle_sum; // record for debugging
			delta[i] = current_angle_sum - two_pi;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "h_riffle: current angle sum = " << current_angle_sum << ", sum-2\\pi error = " << delta[i] << endl;

			if (abs(delta[i]) > cut)
			{
				radius[i] = h_radcalc(flowers, radius, i);
				count++;
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_riffle: radius adjusted to " << radius[i] << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_riffle: radius not adjusted" << endl;
			}
		}
		
//		h_anglesum_deltas (flowers, radius, delta, num_type_1234_vertices);

		accum=0;
		for (int i=0; i< num_type_1234_vertices; i++)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_riffle: node " << i << ": angle sum " << angle_sum[i] << ", difference " << delta[i] << endl;

			accum += abs(delta[i]);
		}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "h_riffle: revised aggregate absolute angle sum " << accum << endl;
		
		cut=accum/3/num_type_1234_vertices;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "h_riffle: angle sum error threshold adjusted to " << cut << endl;
		
	} /* end of while */
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "h_riffle: completed after " << count<< " iterations, cut threhold reached = " << (cut<=TOLER?"true":"false") << endl;
}	

} /* h_riffle */



/* returns best radius at vertex i to
achieve anglesum aim. Use neighboring s_radii from pack and starting 
s_radius s (first guess). This routine is a mixture of binary search and
then Newton's method. `MAX_H_RADCALC_ITERATIONS' gives the max number of iterations to
use in trying for solution. */
double h_radcalc(matrix<int>& flowers, vector<double>& radius, int vertex)  
{
	int num_hradcalc_iterations=0;
	double current_radius = radius[vertex];
	
	double angle_sum;
	double derivative;
	double test_radius;
	double test_angle_sum;

	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_radcalc: called with current current_radius = " << current_radius << endl;
   
	/* evaluate the  current angle sum in angle_sum and it'current_radius derivative in derivative */
	h_curvcalc(flowers,radius,vertex,current_radius,&angle_sum,&derivative); 
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_radcalc: current angle sum = " << angle_sum << ", initial f(r) = " << angle_sum-two_pi << endl;

    /* The current_radius calculation is based on the function f(r) = \theta(r;{r_i}) - 2\pi, the delta between
       the angle sum and 2\pi
       
       the routine is looking for a current_radius that minimizes f (similar to step 3 of section 3.3 but using 
       \theta not \hat\theta)
       
       while |f(current_radius)| > TOLER update the current current_radius by doing Newton-Raphson then 
       some binary searching.
    */
//	while ( (angle_sum > (two_pi+TOLER) || angle_sum < (two_pi-TOLER) ) && num_hradcalc_iterations < MAX_H_RADCALC_ITERATIONS )
	while ( abs(angle_sum-two_pi) > TOLER && num_hradcalc_iterations < MAX_H_RADCALC_ITERATIONS )
	{
		 
		/* do Newton-Raphson on f(r) = \theta(r;{r_i}) - 2\pi starting from r = current current_radius */
		test_radius=current_radius-(angle_sum-two_pi)/derivative;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_radcalc: Newton current_radius = " <<  test_radius << endl;

		if (test_radius<current_radius*0.5) 
		    test_radius=current_radius*0.5;
		    
		if (test_radius>(1+current_radius)*0.5) 
		    test_radius=(1+current_radius)*0.5;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_radcalc: initial trial current_radius  = adjusted Newton current_radius = " << test_radius << endl;

		/* evaluate angle sum with Newton Raphson current_radius in test_angle_sum */
		h_curvcalc(flowers,radius,vertex,test_radius,&test_angle_sum,&derivative);
		
		/* set err to the absolute delta between the current angle sum and 2\pi 
		   that is err = |f(current_radius)|
		*/
		double err=angle_sum-two_pi;
		
		if (err<0) 
		    err=(-err);
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_radcalc: err threshold = f(current_radius) = f(" << current_radius << ") = " << angle_sum-two_pi << endl;

		/* while |f(test_angle_sum)| > err, set the trial current_radius to the value halfway between
		   the trial value and the current current_radius
		   
		   DOESN'T THIS ALWAYS GO THE WRONG WAY?  In any case, the Newton current_radius is closer to the
		   root of f than the initial value, so surely we are never going to get |f(test_angle_sum)| > err,
		   so this code has no effect.
		*/
		while (abs(test_angle_sum-two_pi)>err && num_hradcalc_iterations < MAX_H_RADCALC_ITERATIONS)
		{
			num_hradcalc_iterations++;
			test_radius=(test_radius+current_radius)*0.5;
			h_curvcalc(flowers,radius,vertex,test_radius,&test_angle_sum,&derivative);
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_radcalc: setting trial current_radius to " << test_radius << " gives f(test_radius) = " << test_angle_sum-two_pi << endl;
		}
		num_hradcalc_iterations++;
		current_radius=test_radius;
		h_curvcalc(flowers,radius,vertex,current_radius, &angle_sum,&derivative);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_radcalc: set current_radius to " << current_radius << " with angle sum " << angle_sum << " and f(current_radius) = " << angle_sum-two_pi << endl;
	}
	return current_radius;
} /* h_radcalc */

 /* comp angle sum c and its derivative d (with
respect to the s_radius, for interior nodes in hyp setting. Deriv
depends on whether some of radii are infinite */
void h_curvcalc(matrix<int>& flowers, vector<double>& radius, int vertex, double s, double* c, double* d) 
{
	*c=0;*d=0;
	if (s<=0) return; /* infinite radius at vertex of interest (i.e. it's type 4 - shouldn't happen now */

	double s2=radius[flowers[vertex][1]-1];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_curvcalc: s2 initialized for node " << flowers[vertex][1] << ", radius " << s2 << endl;

	for (int k=2;k<=flowers[vertex][0]+1;k++)
	{
		double s1=s2;
		s2=radius[flowers[vertex][k]-1];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_curvcalc: next s2 is node " << flowers[vertex][k]-1 << ", radius " << s2 << endl;

		double cs=h_comp_cos(s,s1,s2);
		*c+=acos(cs);
		if ((s1>0) && (s2>0))
		 {
			double cc=(-s1*s1-s2*s2);
			double b=s1*s1*s2*s2;
			double a=(-2-cc-2*b);
			double e=cc-a;
			double ss=s*s;
			*d+=2*e*(1-b*ss*ss)/((1+cc*ss+b*ss*ss)*
				sqrt( e*(2+(cc+a)*ss+2*b*ss*ss) ));
		 }
		else if ((s1<=0) && (s2<=0)) *d+=2/sqrt(1-s*s);
		else if (s1<=0) *d+=sqrt(s2/s)/(1+s2*s);
		else *d+=sqrt(s1/s)/(1+s1*s);
	 }
} /*h_curvcalc */

 /* given ordered triple of s_radii, compute
the cosine of the angle at first circle in triangle formed by
mutually tangent circles. */
double h_comp_cos(double s1,double s2,double s3)
{
	double s1s,s2s,s3s,ans;

	if (s1<=0) return (1.0);
	s1s=s1*s1; 
	if ((s2<=0) && (s3<=0)) return (1-2*s1s);
	s3s=s3*s3;
	if (s2<=0) return ((1+s1s*(s3s-2))/(1-s1s*s3s));
	s2s=s2*s2;
 	if (s3<=0) return ((1+s1s*(s2s-2))/(1-s1s*s2s));
 	
	ans=((1+s1s*s2s)*(1+s1s*s3s)-2*s1s*(1+s2s*s3s))/
		((1-s1s*s2s)*(1-s1s*s3s));
	if (ans>1) return 1;
	if (ans<-1) return -1;
	return (ans);
} /* h_comp_cos */


/* converts circle data from hyp to eucl. */
void h_to_e_data(complex<double> h_centre,double s_rad,complex<double>* e_centre,double* e_rad) 
{
	double a,b,d,k,aec,ahc,g;

	if (s_rad<=0) /* infinite hyp radius */
	 {
		a=1+s_rad;
		*e_centre = h_centre*a;
//		e_centre->real()=h_centre.real()*a;
//		e_centre->imag()=h_centre.imag()*a;
		*e_rad=(-s_rad); /* assumes -s_rad is meaningful */
		return;
	 }
//	ahc=cAbs(h_centre);
	ahc=abs(h_centre);
	g=((1+ahc)/(1-ahc));
	a=g/s_rad;
	d=(a-1)/(a+1);
	b=g*s_rad;
	k=(b-1)/(b+1);
	*e_rad = (d-k)*0.5;
	aec=(d+k)*0.5;
	if (ahc<=OKERR)
	 {
		*e_centre = {0,0};
//		e_centre->real()=0;
//		e_centre->imag()=0;
	 }
	else
	 {
		b=aec/ahc;
		*e_centre = h_centre*b;
//		e_centre->real() = b*h_centre.real();
//		e_centre->imag()=b*h_centre.imag();
	 }
} /* h_to_e_data */

/* compute ang sum. */
// i interior vertex, numbered from 1
// s current radius
// c returned angle sum
double h_anglesum (matrix<int>& flowers, vector<double>& radius, int vertex)
{
	double angle_sum = 0;
	double x = radius[vertex];
	
	for (int i=1; i<=flowers[vertex][0]; i++)
	{
		double y = radius[flowers[vertex][i]-1];
		double z = radius[flowers[vertex][i+1]-1]; // requires first petal to be at the end of the list
		
		angle_sum += acos(h_comp_cos(x,y,z));
	}
	return angle_sum;
}

 /* given 2 hyp centers/s_radii and
s_radius of third in ordered triple, return hyp center and s_radius
for third circle. oj is cos of overlap angle opposite to circle j. 
There is no consistency check on first two circles' data. 

oj removed, as KS sets these values to 1.0*/
int h_compcenter(complex<double> z1,complex<double> z2,complex<double>* z3, double s1,double s2,double* s3)
{
	double flag=1,cc,rad,ahc,sc,s,ac,x0,x1,r,sidelength,ss1,ss3;
	double o1=1.0, o2=1.0, o3=1.0;
	complex<double> a,b,c,c1,w1,w2,w3,cent2,newcent3,par;

	a=z1;
	b=z2;
	if ((s1<=0) && (s2>0)) /* second circle finite, first not */
	 {
		a=z2; b=z1; s=s1; s1=s2; s2=s; s=o1; o1=o2; o2=s; flag=(-1); 
			/* interchange order */
	 }
	if (s1>0) /* first is now finite */
	 {
		cc=h_cos_overlap(s1,s2,*s3,o1,o2,o3);
		z3->real(cc);
		z3->imag(flag*sqrt(1-cc*cc));
		if (*s3>0)
		 {
			ss1=s1*s1;
			ss3=*s3*(*s3);
			sidelength=exp( acosh(
			  (1/(4.0*s1*(*s3)))*
			  ((1+ss3)*(1+ss1)+(1-ss3)*(1-ss1)*o2) ) );
			ahc=(sidelength-1)/(sidelength+1); /* abs value of hyp 
							center */
			*z3 *= ahc;
//			z3->real() *= ahc;
//			z3->imag() *= ahc; /* center as if z1 at origin */
			*z3=mob_norm_inv(*z3,a,b);   /* move to right place */
			return (0);
		 }
		r=(1-s1)/(1+s1);
		sc=(r*r+1+2*r*o2)/(2*(1+r*o2));
			/* abs value of eucl center c of third circle */
		c = (*z3)*sc;
//		c.real()=sc*z3->real();
//		c.imag()=sc*z3->imag();
		rad=1-sc; /* Now have c and its radius */
		w1.real(c.real()-rad);
		w1.imag(c.imag());
		w2.real(c.real()+rad);
		w2.imag(c.imag());
		w3.real(c.real());
		w3.imag(c.imag()+rad); /* three points on the circle */
		w1=mob_norm_inv(w1,a,b);
		w2=mob_norm_inv(w2,a,b);
		w3=mob_norm_inv(w3,a,b);       /* 3 pts on new circle */
		circle_3(w1,w2,w3,&c,&rad);
	 	*s3=(-rad); /* store new s_radius */
//		ac=cAbs(c);
		ac=abs(c);
		*z3 = c/ac;
//		z3->real()=c.real()/ac;
//		z3->imag()=c.imag()/ac; /* get hyp center on unit circle */
		return (0);
	 }
	/* remaining: first 2 have infinite rad.*/
	if (*s3<=0) /* third also infinite. As temp shortcut, we make 
		radius large but finite. */
		*s3=OKERR;
	/* Now, third finite.
		Pretend 3 is at origin, locate first/second on unit circle.
		Apply Mobius to make first correct eucl size (keeping real
		axis fixed), move back by rotation to get first to 
		original location. */
	cent2.real(h_cos_overlap(*s3,s1,s2,o3,o1,o2));
	cent2.imag(sqrt(1-cent2.real()*cent2.real())); /* circle 2 center in 
		normalize position */
	r=(1-*s3)/(1+*s3); /* eucl radius of circle 3 when at origin */
	x0=1-(r*r+1+2*r*o2)/(2*(1+r*o2)); /* eucl center of circle 1 
		in normalized position - tangent at 1. */
	x1=1-fabs(s1); 	 /* desired eucl center of cir 1 */
	if (x1<OKERR || x1>(1-OKERR)) x1=.5; /* in case s1 
		hadn't been set */
	par.real((x1-x0)/(x1*x0-1));
	par.imag(0.0); /* set parameter for mobius 
			(z-par)/(1-conj(par)*z) */
	newcent3.real((-1)*par.real());
	newcent3.imag(0.0); /* this is just mob_trans(ORIGIN,par) */
//	*z3=cmult(z1,newcent3); /* rotate to put circle 1 back in right location */
	*z3=z1*newcent3; /* rotate to put circle 1 back in right location */
	return (0);	
} /* h_compcenter */

 /* given ordered triple of s_radii and
cosines of overlap angles, compute cosine of angle at the first
circle of triangle they form. Note: tj is cos of overlap angle opposite 
circle j. */
double h_cos_overlap(double s1,double s2,double s3,double t1,double t2,double t3)
{
	double h1,h2,h3,e1,e2,e3,L,len,ans;

	if (s1<=0) return (1.0);
	e1=(1.0-s1)/(1.0+s1);
	h1=-log(s1);
	if (s2<=0) L=1.0;
	else 
	 {
		h2=-log(s2);
		len=exp(-h2-acosh(cosh(h1)*cosh(h2)+sinh(h1)*sinh(h2)*t3));
		L=(1.0-len)/(1.0+len);
	 }
	e2=(L*L-e1*e1)/(2.0*e1*t3+2.0*L);
	if (s3<=0) L=1.0;
	else 
	 {
		h3=-log(s3);
		len=exp(-h3-acosh(cosh(h1)*cosh(h3)+sinh(h1)*sinh(h3)*t2));
		L=(1.0-len)/(1.0+len);
	 }
	e3=(L*L-e1*e1)/(2.0*e1*t2+2.0*L);
		/* have euclidean radii now */
	ans=e_cos_overlap(e1,e2,e3,t1,t2,t3); /* find cos from eucl routine */
	return (ans);
} /* h_cos_overlap */

/* given three eucl radii and cosines
of opposite overlap angles, return cos of angle at e1 */
double e_cos_overlap(double e1,double e2,double e3,double t1,double t2,double t3) 
{
	double l1,l2,l3,ang;

	l3=e1*e1+e2*e2+2*e1*e2*t3;
	l2=e1*e1+e3*e3+2*e1*e3*t2;
	l1=e2*e2+e3*e3+2*e2*e3*t1;
	ang=(l2+l3-l1)/(2*sqrt(l2*l3));
	return (ang);
} /* e_cos_overlap */

complex<double> mob_norm_inv(complex<double> w,complex<double> a,complex<double> b) 
/* returns preimage of w under mobius of disc
which maps a to zero and b to positive x-axis */
{
	complex<double> z;
	double c;

	z=mob_trans(b,a);
//	c=cAbs(z);
	c=abs(z);
	z /= c;
//	z.real()/=c;
//	z.imag()/=c;
//	z=mob_trans(cmult(w,z),a);
	z=mob_trans(w*z,a);
	return (z);
} /* mob_norm_inv */

complex<double> mob_trans(complex<double>z,complex<double>a) 
/* return value for (a-z)/(1-z*conj(a)) */
{
	complex<double> w,COMP_UNIT;

	COMP_UNIT.real(1.0);
	COMP_UNIT.imag(0.0);
	
//	w=cdiv(csub(a,z),csub(COMP_UNIT,cmult(z,cconj(a)) ) );
	w= (a-z)/(COMP_UNIT-z*conj(a));
	return (w);
} /* mob_trans */

/* find eucl center and rad for circle through 3 given points. Returns 1 in case of success. */
int circle_3(complex<double> z1,complex<double> z2,complex<double> z3,complex<double>* c,double* rad) 
{
	double det,a1,a2,b1,b2,c1,c2,dum;

	a1=z2.real()-z1.real();
	a2=z3.real()-z2.real();
	b1=z2.imag()-z1.imag();
	b2=z3.imag()-z2.imag();
	det=2*(a1*b2-b1*a2);
	if (fabs(det)<OKERR) return 0;
//	dum=cAbs(z2);
	dum=abs(z2);
	dum=dum*dum;
//	c1=cAbs(z1);
	c1=abs(z1);
	c1=(-c1*c1);
	c1+=dum;
//	c2=cAbs(z3);
	c2=abs(z3);
	c2=c2*c2;
	c2-=dum;
	c->real((b2*c1-b1*c2)/det);
	c->imag((a1*c2-a2*c1)/det);
//	*rad=cAbs(csub(*c,z1));
	*rad=abs(*c-z1);
	return 1;
} /* circle_3 */


	/* we place vertex zero at the origin and the vertex in the flower of vertex zero along the x-axis, 
	   then work around the vertex flowers in triples: given two centres and radii, and a third radii, 
	   we calculate the third centre.	   
	*/
void compute_circle_pack_centres(matrix<int>& flowers, vector<double>& radius, int nodecount, int num_type_1234_vertices, matrix<double>& vertex_centre)
{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "compute_circle_pack_centres:  compute vertex centres from radii" << endl;
	for (int i = 0; i< nodecount; i++)
		debug << "compute_circle_pack_centres: vertex " << i << " " << radius[i] << endl;
}

	vertex_centre[0][0] = 0;
	vertex_centre[0][1] = 0;
	
	int first_petal = flowers[0][1]-1;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "compute_circle_pack_centres:   vertex 0 centre = (" << vertex_centre[0][0] << "," << vertex_centre[0][1] << ")" << endl;
	debug << "compute_circle_pack_centres:   first_petal = " << first_petal << endl;
}		
	/* we don't want to put a type 4 vertex on the x-axis, so make sure that the first petal is an interior vertex.
	   to accommodate infinite turning cycles with just two edges, we may have to look past the second petal in 
	   flowers as well as the first.
	*/
	if (first_petal >= num_type_1234_vertices)
	{
		first_petal = flowers[0][2]-1;

		if (first_petal >= num_type_1234_vertices)
			first_petal = flowers[0][3]-1;
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:   first_petal is a type 4 vertex, changing first petal to " << first_petal << endl;

	}

//	vertex_centre[first_petal][0] = radius[0]+radius[first_petal];
	{
		double s1 = radius[0];
		double s2 = radius[first_petal];
		double ss1=s1*s1;
		double ss2=s2*s2;
		double ovlp=1.0;
		double s=exp( acosh((1/(4.0*s1*s2))*((1+ss1)*(1+ss2)+(1-ss1)*(1-ss2)*ovlp) ) );
		vertex_centre[first_petal][0]=(s-1)/(s+1);
	 }
	vertex_centre[first_petal][1] = 0;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:   first_petal " << first_petal << " centre = (" << vertex_centre[first_petal][0] << "," << vertex_centre[first_petal][1] << ")" << endl;

	vector<bool> computed_centre(nodecount); // default boolean constructed as false
	computed_centre[0] = true;
	computed_centre[first_petal] = true;
	
	list<int> vertex_list;
	vertex_list.push_back(0); // start with the first vertex.
	vertex_list.push_back(first_petal); // we've just set the centre of the first petal.

	while (vertex_list.size() != 0)
	{
		int vertex = *vertex_list.begin();
		vertex_list.pop_front();

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "compute_circle_pack_centres:    pop " << vertex << " from front of vertex_list" << endl;
	debug << "compute_circle_pack_centres:    petals: ";
	for (int i=1; i<= flowers[vertex][0]; i++)
			debug << flowers[vertex][i] -1 << ' ';
	debug << endl;
}
		
		/* we only consider internal vertices because all type 4 vertices lie in the flower 
		   of an internal vertex, so we will compute the type 4 centres from an internal vertex
		*/
		if (vertex >= num_type_1234_vertices)
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:      vertex is a type 4 vertex, do nothing" << endl;
			continue;
		}
		
		/* we have already computed the centre of at least one of the neighbours of
		   this vertex but we need to identify where it is in the corresponding row
		   of flowers.  We store the index of the first neighbour for which the centre
		   is known in first_neighbour.
		*/
		int first_neighbour;  
		for (int i=1; i<flowers[vertex][0]; i++)
		{
			if (computed_centre[flowers[vertex][i]-1] == true)
			{
				first_neighbour = i;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "compute_circle_pack_centres:      first neighbour for which we already have computed the centre is " << flowers[vertex][i]-1
	      << ", index " << first_neighbour << endl;
}	      	      
				break;
			}
		}
		
		/* We may now work around the vertex from first neighbour computing any other centres 
		   we've not yet calculated.  If there are n neighbours of the vertex, we need to cycle
		   the columns of flowers indexed 1 to n starting from first_neighbour+1, which we do by 
		   cycling the indices modulo n allowing for the fact that column 0 stores the value n.
		   Note that we do not want to return to first_neighbour, so we cycle through vertices
		   1,...,n-1 mod n.
		*/
		int n = flowers[vertex][0];
		for (int i=1; i<n; i++)
		{
			int c2_index = ((first_neighbour-1)+i-1)%n +1;
			int c3_index = ((first_neighbour-1)+i)%n +1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:      c2_index = " << c2_index << ", c3_index = " << c3_index << endl;
			
			if (computed_centre[flowers[vertex][c3_index]-1] == false)
			{
				/* we use c1, c2, and c3 to identify the vertices with the corresponding centres */
				int c1 = vertex;
				int c2 = flowers[vertex][c2_index]-1;
				int c3 = flowers[vertex][c3_index]-1;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:      computing centre for vertex c3 = " << c3 << ", c2 = " << c2 << ", c1 = " << c1 << endl;

				complex<double>z1;
				z1.real(vertex_centre[c1][0]);
				z1.imag(vertex_centre[c1][1]);
				
				complex<double>z2;
				z2.real(vertex_centre[c2][0]);
				z2.imag(vertex_centre[c2][1]);

				complex<double> z3;
				double s3 = radius[c3];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "compute_circle_pack_centres:      s3 = " << s3 << " &s3 = " << &s3 << endl;
				h_compcenter(z1,z2,&z3,radius[c1],radius[c2],&s3);
				
				vertex_centre[c3][0] = z3.real();
				vertex_centre[c3][1] = z3.imag();
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "compute_circle_pack_centres:      s3 = " << s3 << " radius[c3] = " << radius[c3] << endl;
				radius[c3]=s3;
				
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:      c3 centre = (" << vertex_centre[c3][0] << "," << vertex_centre[c3][1] << ")" << endl;
				
				/* push vertex c3 onto the vertex_list */
				vertex_list.push_back(flowers[vertex][c3_index]-1);
				computed_centre[flowers[vertex][c3_index]-1] = true;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:    push " << flowers[vertex][c3_index]-1 << " onto back of vertex_list" << endl;

			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "compute_circle_pack_centres:      already computed centre for vertex c3 = " << flowers[vertex][c3_index]-1 << endl;				
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "compute_circle_pack_centres: vertex centres: " << endl;
	print (vertex_centre, debug, 10,"compute_circle_pack_centres: ");
}
}

void h_anglesum_deltas (matrix<int>& flowers, vector<double>& radius, vector<double>& delta, int num_type_1234_vertices)
{
	for (int i=0; i< num_type_1234_vertices; i++)
	{

		double angle_sum = h_anglesum (flowers, radius, i);
		double difference = angle_sum - two_pi;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "h_anglesum_deltas: vertex " << i << ", angle sum = " << angle_sum << ", difference = " << difference << endl;					

		delta[i]  = difference;
	}
}


/* normalizes hyperbolic data of pack p by 
putting point c at origin and d on pos y-axis */
void h_norm_pack(matrix<double>& vertex_centre,int alpha,int gamma) 
{
	int nodecount = vertex_centre.numrows();
	
	complex<double> c(vertex_centre[alpha][0],vertex_centre[alpha][1]);
	complex<double> d(vertex_centre[gamma][0],vertex_centre[gamma][1]);

	double abd;

//	if (cAbs(c)>=OKERR) /* j vertex not origin */
	if (abs(c)>=OKERR) 
	{
		/* j vertex not origin */
//		for (int i=1;i<=p->nodecount;i++)f
		for (int i=0;i< nodecount;i++)
		{
			complex<double> i_centre(vertex_centre[i][0],vertex_centre[i][1]);
			complex<double> shifted_i_centre=mob_norm(i_centre,c,d);
			vertex_centre[i][0] = shifted_i_centre.real();
			vertex_centre[i][1] = shifted_i_centre.imag();
		}
	}
	else if (d.imag()>=OKERR || d.imag()<=(-OKERR) || d.real()<=(-OKERR))		
	{
		/* just a rotation is needed */
//		abd=1.0/cAbs(d);
		abd=1.0/abs(d);
		d *= abd;
//		d.real() *= abd;
//		d.imag() *= abd;
		
//		for (int i=1;i<=p->nodecount;i++)
		for (int i=0;i< nodecount;i++)
		{
//			pR_ptr[i].center=cdiv(pR_ptr[i].center,d);
			complex<double> i_centre(vertex_centre[i][0],vertex_centre[i][1]);
			complex<double> shifted_i_centre= i_centre/d;
			vertex_centre[i][0] = shifted_i_centre.real();
			vertex_centre[i][1] = shifted_i_centre.imag();
		}
	 }

	/* now rotate another 90 degrees --- need to clean this up later */
//	for (int i=1;i<=p->nodecount;i++)
	for (int i=0;i< nodecount;i++)
	{
//		temp=pR_ptr[i].center;
//		pR_ptr[i].center.re()=(-temp.im());
//		pR_ptr[i].center.im()=temp.re();
		complex<double> temp(vertex_centre[i][0],vertex_centre[i][1]);
		complex<double> shifted_i_centre(-temp.imag(),temp.real());
		vertex_centre[i][0] = shifted_i_centre.real();
		vertex_centre[i][1] = shifted_i_centre.imag();
	}
} /* h_norm_pack */

/* returns value at z of mobius transform of unit
disc which maps a to zero and b to positive x-axis. */
complex<double> mob_norm(complex<double> z,complex<double> a,complex<double> b) 
{
	complex<double> w,v;
	double c;

	v=mob_trans(b,a);
//	c=cAbs(v);
	c=abs(v);
	if (c<OKERR) return (mob_trans(z,a));
//	w=cdiv(mob_trans(z,a),v);
	w=mob_trans(z,a)/v;
	w *=c;
//	w.real() *= c;
//	w.imag() *= c;
	return (w);
} /* mob_norm */
