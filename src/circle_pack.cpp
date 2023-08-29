/************************************************************************
                  Roger Fenn Circle packing algorithm


**************************************************************************/
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

extern bool ONE_METAPOST_PATH;
extern double badness_threshold; 
extern int metapost_coordinate_scale_factor;
extern int metapost_hyperbolic_scale_factor;
extern int DRAW_IN_HYPERBOLIC_DISC;


#include <util.h>
#include <matrix.h>
#include <draw.h>

double two_pi = 6.28318530718;
double epsilon = .000001;
extern int max_circle_packing_iterations;

extern char const* circlepack_output_file;
extern char const* triangulation_output_file;

/********************* Function prototypes ***********************/
void hyperbolic_representation(matrix<double>& centre, vector<double>& radius);
void mobius_transform_to_origin(complex<double>& a, complex<double>& b, complex<double> z);
complex<double> mobius_transformation(complex<double> a, complex<double> b, complex<double> z);
complex<double> inv_mobius_transformation(complex<double> a, complex<double> b, complex<double> z);
double max_delta(vector<double>& delta,double reference);
double abs_delta_sum(vector<double>& delta);
double calculate_angle_sum (matrix<int>& flowers, vector<double>& radius, int vertex);
void calculate_angle_sum_deltas (matrix<int>& flowers, vector<double>& radius, vector<double>& delta, int num_type_1234_vertices);

int circle_pack(char const* infile, char const* outfile)
{
	int nodecount;
	int alpha;
	int beta;
	int gamma;

	/* Read triangulation data from triangulation_output_file */
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack: read triangulation data from" << triangulation_output_file << endl;

	ifstream input(triangulation_output_file);
	if(!input)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "circle_pack: could not open" << triangulation_output_file << endl;
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
	debug << "circle_pack: nodecount = " << nodecount << endl;
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
	debug << "circle_pack: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
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
	debug << "circle_pack: node " << node << ", petal count = " << count << ": " << first_petal;
				
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
	
	/* set the intial circle radii */
	int num_type_1234_vertices = beta-1;
	vector<double> radius(nodecount);
	for (int i=0; i< num_type_1234_vertices; i++)
		radius[i] = 10; //set interior vertex initial radius
	for (int i = num_type_1234_vertices; i< nodecount; i++)
		radius[i] = 75; // set boundary vertex radius

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "circle_pack: initial radii ";
	for (int i = 0; i< nodecount; i++)
		debug << radius[i] << ' ';
	debug << endl;
}
	
	/* establish a vector to record the delta between a vertex's angle sum and 2\pi */
	vector<double> delta(num_type_1234_vertices);
	calculate_angle_sum_deltas (flowers, radius, delta, num_type_1234_vertices);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "circle_pack: initial angle sum deltas ";
	for (int i = 0; i< num_type_1234_vertices; i++)
		debug << setprecision(2) << delta[i] << ' ';
	debug << endl;
}
	
	int num_iterations = 0;

	/* Iterate until the maximum angle sum error vector is less than epsilon 
	
	original version:
	*/
	
	while (max_delta(delta,two_pi) >= epsilon && num_iterations < max_circle_packing_iterations)
	{
		bool positive_radii;
		int delta_factor = 2;
		vector<double> test_radius(num_type_1234_vertices);
		do 
		{
			positive_radii = true;

			for (int i=0; i< num_type_1234_vertices; i++)
			{
				test_radius[i] = radius[i]+delta[i]/delta_factor;
				
				if (test_radius[i] < 0)
				{
					positive_radii = false;
					delta_factor*=2;
					break;
				}
			}
			
		} while (positive_radii == false);
		
		for (int i=0; i< num_type_1234_vertices; i++)
			radius[i] = test_radius[i];
			
		calculate_angle_sum_deltas (flowers, radius, delta, num_type_1234_vertices);

		num_iterations++;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "circle_pack:  " << num_iterations << ": adjusted radii ";
	for (int i = 0; i< num_type_1234_vertices; i++)
		debug << radius[i] << ' ';
	debug << endl;
	debug << "circle_pack:  " << num_iterations << ": angle sum deltas ";
	for (int i = 0; i< num_type_1234_vertices; i++)
		debug << setprecision(2)<< delta[i] << ' ';
	debug << endl;
}
	}



	/* Knotscape approach 
	double absolute_aggregate_delta_sum = abs_delta_sum(delta);
	double test_threshold = absolute_aggregate_delta_sum/3/num_type_1234_vertices;
	
//	while (test_threshold >= epsilon && num_iterations < max_circle_packing_iterations)
	while (test_threshold >= 0.02 && num_iterations < max_circle_packing_iterations)
	{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "circle_pack:  " << num_iterations << " absolute aggregate delta = " << absolute_aggregate_delta_sum 
	      << ", test_threshold = " << test_threshold << endl;
}

		for (int i=0; i< num_type_1234_vertices; i++)
		{
			delta[i] = calculate_angle_sum (flowers, radius, i) - two_pi;
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "circle_pack:  " << num_iterations << ": vertex " << i << ", delta = " << delta[i];
			
			if (abs(delta[i]) > epsilon)
			{
				double test_radius;
				int delta_factor = 2;
				
				do 
				{
					test_radius = radius[i]+delta[i]/delta_factor;
						
					if (test_radius < 0)
						delta_factor*=2;
					
				} while (test_radius < 0);
				
				radius[i] = test_radius;
				
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << ", radius adjusted to " << radius[i] << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << ", not adjusted" << endl;
			}
		}

		absolute_aggregate_delta_sum = abs_delta_sum(delta);
		test_threshold = absolute_aggregate_delta_sum/3/num_type_1234_vertices;
		num_iterations++;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "circle_pack:  " << num_iterations << ": adjusted radii ";
	for (int i = 0; i< num_type_1234_vertices; i++)
		debug << radius[i] << ' ';
	debug << endl;
	debug << "circle_pack:  " << num_iterations << ": angle sum deltas ";
	for (int i = 0; i< num_type_1234_vertices; i++)
		debug << setprecision(2)<< delta[i] << ' ';
	debug << endl;
}
	}
*/

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "circle_pack: final angle sum deltas ";
	for (int i = 0; i< num_type_1234_vertices; i++)
		debug << delta[i] << ' ';
	debug << endl;
	debug << "circle_pack: num_iterations = " << num_iterations << endl;
}

	/* compute the circle centres based on the final radii */
	matrix<double> vertex_centre(nodecount,2);
	/* we place vertex zero at the origin and the vertex in the flower of vertex zero along the x-axis, 
	   then work around the vertex flowers in triples: given two centres and radii, and a third radii, 
	   we calculate the third centre.
	   
	   Let the known centres be c1 and c2, with radii r1 and r2 and suppose we wish to place an osculating
	   circle at c3 of radius r3, where we move anticlockwise around c1 to move from c2 to c3.
	   
	   The cosine rule gives us the angle c2c1c3, \alpha say, so we calculate c3 as follows:
	   
	    1. translate so that c1' lies at the origin
	    2. rotate so that c2' lies on the x-axis
	    3. use \alpha to determine c3' from r1 and r3
	    4. rotate and translate back
	*/
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:  compute vertex centres from radii" << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "circle_pack: final radii ";
	for (int i = 0; i< num_type_1234_vertices; i++)
		debug << radius[i] << ' ';
	debug << endl;
}


	vertex_centre[0][0] = 0;
	vertex_centre[0][1] = 0;
	
	int first_petal = flowers[0][1]-1;
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "circle_pack:   vertex 0 centre = (" << vertex_centre[0][0] << "," << vertex_centre[0][1] << ")" << endl;
	debug << "circle_pack:   first_petal = " << first_petal << endl;
}		
	/* we don't want to put a type 4 vertex on the x-axis, so make sure that the
	   first petal is an interior vertex.
	*/
	if (first_petal >= num_type_1234_vertices)
	{
		first_petal = flowers[0][2]-1;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:   first_petal is a type 4 vertex, changing first petal to " << first_petal << endl;

	}

	vertex_centre[first_petal][0] = radius[0]+radius[first_petal];
	vertex_centre[first_petal][1] = 0;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:   first_petal " << first_petal << " centre = (" << vertex_centre[first_petal][0] << "," << vertex_centre[first_petal][1] << ")" << endl;

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

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "circle_pack:    pop " << vertex << " from front of vertex_list" << endl;
	debug << "circle_pack:    petals: ";
	for (int i=1; i<= flowers[vertex][0]; i++)
			debug << flowers[vertex][i] -1 << ' ';
	debug << endl;
}
		
		/* we only consider internal vertices because all type 4 vertices lie in the flower 
		   of an internal vertex, so we will compute the type 4 centres from an internal vertex
		*/
		if (vertex >= num_type_1234_vertices)
		{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      vertex is a type 4 vertex, do nothing" << endl;
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
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "circle_pack:      first neighbour for which we already have computed the centre is " << flowers[vertex][i]-1
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

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      c2_index = " << c2_index << ", c3_index = " << c3_index << endl;
			
			if (computed_centre[flowers[vertex][c3_index]-1] == false)
			{
				/* we use c1, c2, and c3 to identify the vertices with the corresponding centres */
				int c1 = vertex;
				int c2 = flowers[vertex][c2_index]-1;
				int c3 = flowers[vertex][c3_index]-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      computing centre for vertex c3 = " << c3 << ", c2 = " << c2 << ", c1 = " << c1 << endl;
				
				/* evaluate c2', the image of c2 after translating c1 to the origin */
				double c2_prime_x = vertex_centre[c2][0] - vertex_centre[c1][0];
				double c2_prime_y = vertex_centre[c2][1] - vertex_centre[c1][1];

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "circle_pack:      c1 centre = (" << vertex_centre[c1][0] << "," << vertex_centre[c1][1] << ")" << endl;
	debug << "circle_pack:      c2 centre = (" << vertex_centre[c2][0] << "," << vertex_centre[c2][1] << ")" << endl;
	debug << "circle_pack:      c2_prime = (" << c2_prime_x << "," << c2_prime_y << ")" << endl;
}	

				
				/* the angle beta rotates c2' to the x-axis, the quadrant containing
				   c2' determines whether this is added or subtracted from 0, 90, 180
				   or 270 degrees.
				*/
				double beta = atan( abs(c2_prime_y)/abs(c2_prime_x));

				/* alpha is the angle c2c1c3 */
				double a = radius[c1]+radius[c3];
				double b = radius[c1]+radius[c2];
				double c = radius[c2]+radius[c3];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "circle_pack:      radius[c1] = " << radius[c1] << ", radius[C2] = " << radius[c2] << ", radius[c3] = " << radius[c3] << endl;
	debug << "circle_pack:      a = " << a << ", b = " << b << ", c = " << c << endl;
}

				double cosine = (a*a+b*b-c*c)/2/a/b;

				if (cosine > 1)
				{
					cosine = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to 1 " << endl;	
				}
				else if (cosine < -1)
				{
					cosine = -1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to -1 " << endl;
				}

				double alpha = acos(cosine);
					
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "circle_pack:      cosine = " << cosine << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      beta = " << beta << "(" << beta*360/two_pi << ")" << ", alpha = " << alpha << "(" << alpha*360/two_pi << ")" << endl;

				/* adjust alpha by beta before calculating c3' (thus our c3' is really the image of c3'
				   as described above after the rotation by beta)
				   
				   We have to consider each quadrant due to the ambiguity of the acos function.
				*/
				if (c2_prime_x > 0 && c2_prime_y > 0)
					alpha += beta;
				else if (c2_prime_x < 0 && c2_prime_y > 0)
					alpha = alpha + two_pi/2 - beta; //alpha + \pi - beta
				else if (c2_prime_x < 0 && c2_prime_y < 0)
					alpha = alpha + two_pi/2 + beta; //alpha + \pi + beta
				else // (c2_prime_x > 0 && c2_prime_y < 0)
					alpha -= beta;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      alpha after beta adjustment = " << alpha << "(" << alpha*360/two_pi << ")" << endl;

				
				/* now work out where c3' is with the new alpha */
				double c3_prime_x = a*cos(alpha);
				double c3_prime_y = a*sin(alpha);

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      c3_prime = (" << c3_prime_x << "," << c3_prime_x << ")" << endl;
				
				/* finally, translate back to c3 */
				vertex_centre[c3][0] = c3_prime_x + vertex_centre[c1][0];
				vertex_centre[c3][1] = c3_prime_y + vertex_centre[c1][1];

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      c3 centre = (" << vertex_centre[c3][0] << "," << vertex_centre[c3][1] << ")" << endl;
				
				/* push vertex c3 onto the vertex_list */
				vertex_list.push_back(flowers[vertex][c3_index]-1);
				computed_centre[flowers[vertex][c3_index]-1] = true;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:    push " << flowers[vertex][c3_index]-1 << " onto back of vertex_list" << endl;

			}
			else
			{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "circle_pack:      already computed centre for vertex c3 = " << flowers[vertex][c3_index]-1 << endl;				
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "circle_pack: vertex centres: " << endl;
	print (vertex_centre, debug, 10,"circle_pack: ");
}

	/* write the centres and radii of the interior vertices to the output file */
	ofstream output(circlepack_output_file);
	
	if (!output)
	{
		cout << "\nError opening output file\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "circle_pack: could not open " << circlepack_output_file << endl;
		exit(0);
	}		
	
/*
	output << num_type_1234_vertices << endl;
	for (int i=0; i< num_type_1234_vertices; i++)
	{
		output << setprecision(6) << setw(10) << vertex_centre[i][0] 		
		       << setprecision(6) << setw(10) << vertex_centre[i][1] << endl;		
	}

	for (int i=0; i< num_type_1234_vertices; i++)
		output << setprecision(6) << setw(10) << radius[i] << endl;		
*/

	/* write the centres and radii of the type 4 vertices to the output file */
	output << nodecount << endl;
//	for (int i= num_type_1234_vertices; i< nodecount; i++)
	output << std::fixed;
	
	for (int i= 0; i< nodecount; i++)
	{
		output << setprecision(6) << setw(10) << vertex_centre[i][0] 		
		       << setprecision(6) << setw(10) << vertex_centre[i][1] << endl;		
	}
	for (int i= 0; i< nodecount; i++)
		output << fixed << setprecision(6) << setw(10) << radius[i] << endl;		

	output.close();
	return 1;
}


double max_delta(vector<double>& delta,double reference)
{
	double max = abs(delta[0]-reference);
	
	for (unsigned int i=1; i< delta.size(); i++)
	{
		double loc = abs(delta[i]-reference);
		if (loc > max)
			max = loc;
	}
	
	return max;
}

double abs_delta_sum(vector<double>& delta)
{
	double aggregate_delta = 0;
	
	for (unsigned int i=0; i< delta.size(); i++)
		aggregate_delta += abs(delta[i]);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "abs_delta_sum: aggregate absolute delta = " << aggregate_delta << endl;
	
	return aggregate_delta;
}

void calculate_angle_sum_deltas (matrix<int>& flowers, vector<double>& radius, vector<double>& delta, int num_type_1234_vertices)
{
	for (int i=0; i< num_type_1234_vertices; i++)
	{

		double angle_sum = calculate_angle_sum (flowers, radius, i);
		double difference = angle_sum - two_pi;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum_deltas: vertex " << i << ", angle sum = " << angle_sum << ", difference = " << difference << endl;					

		delta[i]  = difference;
	}
}

double calculate_angle_sum (matrix<int>& flowers, vector<double>& radius, int vertex)
{
	double angle_sum = 0;
	double x = radius[vertex];
	
	for (int i=1; i<=flowers[vertex][0]; i++)
	{
		double y = radius[flowers[vertex][i]-1];
		double z = radius[flowers[vertex][i+1]-1]; // requires first petal to be at the end of the list
		
		double a = x+z;
		double b = x+y;
		double c = y+z;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: x = " << y << ", y = " << y << ", z = " << z << ", a = " << a << ", b = " << b << ", c = " << c << endl;					
		
		double cosine = (a*a+b*b-c*c)/2/a/b;

		if (cosine > 1)
		{
			cosine = 1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to 1 " << endl;	
		}
		else if (cosine < -1)
		{
			cosine = -1;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to -1 " << endl;
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", increment = " << acos(cosine) << endl;


		angle_sum += acos(cosine);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: angle_sum incremented to " << angle_sum << endl;
	}
	return angle_sum;
}
