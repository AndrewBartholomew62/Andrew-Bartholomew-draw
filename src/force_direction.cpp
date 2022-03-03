/************************************************************************
                  Force directed vertex placement


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
extern ofstream 	output;

#include <util.h>
#include <matrix.h>
#include <draw.h>

extern char const* circlepack_output_file;
extern char const* triangulation_output_file;
extern float boundary_vertex_retraction_factor;
extern float average_triangulation_length_factor;
extern bool INCLUDE_BOUNDARY_VERTICES;
extern bool DRAW_IN_HYPERBOLIC_DISC;
extern bool USE_CENTRE_OF_GRAVITY_PLACEMENT;
extern bool RETRACT_BOUNDARY_VERTICES;
extern bool TRACK_PLACEMENT_ITERATION;
extern bool USE_KEN_STEPHENSON_CIRCLE_PACKING;
extern bool PLESTENJAK_FORCE_DIRECTION;
extern bool APPLY_FORCES_TO_TYPE_12_ONLY;

extern int metapost_coordinate_scale_factor;
extern int placement_iteration_tracking_step;

double EPSILON = 0.001;
double ZERO_THRESHOLD = 1; // 0.001;
double pi = 3.14159265358979;

double repulsive_force(double distance);
double attractive_force (double distance);
double k_constant;


/********************* Function prototypes ***********************/
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, matrix<double> vcoords, vector<double> vertex_radius, vector<int>* _state=0);
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, 
                    char const* circlepack_output_file, char const* circlepack_output_savefile, vector<int>* _state=0);
void set_frame_corners(string coordinate_file, metapost_control& mp_control);
int circle_pack (char const* inputfile, char const* outputfile);
int KS_circle_pack (char const* inputfile, char const* outputfile);

void force_directed_placement(metapost_control& mp_control, generic_code_data& code_data, int num_iterations, string title)
{

	int num_crossings = code_data.num_crossings;
	int num_edges = 2 * num_crossings;
	int num_type12_vertices = num_crossings+num_edges;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement: code_data: ";
	write_code_data(debug,code_data);
	debug << endl;
	print_code_data(code_data, debug, "force_directed_placement: ");
	debug << "force_directed_placement: num_iterations = " << num_iterations << endl;
	if (TRACK_PLACEMENT_ITERATION)
		debug << "force_directed_placement: placement_iteration_tracking_step = " << placement_iteration_tracking_step << endl;
}

	/* record the frame corners for the force directed placement case */
	mp_control.frame_minx = -metapost_coordinate_scale_factor;
	mp_control.frame_maxx = metapost_coordinate_scale_factor;
	mp_control.frame_miny = -metapost_coordinate_scale_factor;
	mp_control.frame_maxy = metapost_coordinate_scale_factor;

	/* turn off unsupported options */
	DRAW_IN_HYPERBOLIC_DISC = false;
	mp_control.circle_packing = false;
		
	/* use circle packing to provide an initial placement of the vertices, where necessary */
	if (!PLESTENJAK_FORCE_DIRECTION)
	{
		if (USE_KEN_STEPHENSON_CIRCLE_PACKING)
			KS_circle_pack (triangulation_output_file, circlepack_output_file);
		else
			circle_pack (triangulation_output_file, circlepack_output_file);
	}

	/* Read triangulation data from triangulation_output_file */
	int nodecount;
	int alpha;
	int beta;
	int gamma;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "force_directed_placement: read triangulation data from " << triangulation_output_file << endl;

	ifstream triangulation(triangulation_output_file);
	if(!triangulation)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: could not open " << triangulation_output_file << endl;
		exit(0);
	}		
	
	string next_line;
	getline(triangulation,next_line);
	if (next_line.find("NODECOUNT") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> nodecount;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "force_directed_placement: nodecount = " << nodecount << endl;
	}

	getline(triangulation,next_line);
	if (next_line.find("ALPHA") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> alpha; //always 1, the first vertex
		iss >> beta; //the first boundary (type 5) vertex
		iss >> gamma; // IS THIS USED?
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "force_directed_placement: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
	}
	
	int num_type123_vertices = beta-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "force_directed_placement: number of type 1,2 & 3 vertices num_type123_vertices = " << num_type123_vertices << endl;

	/* Read the boundary orientation (i.e whether the infinite region was bounded by a left turning cycle or 
	   a right turning cycle) in order to place the boundary vertices correctly.
	   
	   Store the flowers for each node in a matrix with the first column set to the number of neighbours
	   of that vertex.  No vertex is adjacent to itself, nor is any vertex adjacent to all other vertices
	   so the maximum number of neighbours is <= nodecount -2.  Recall that for interior vertices the first
	   petal appears at the end of the list as well, to simplify the calculation of angle sums in circle packing.
	*/
	bool clockwise_boundary = true;
	
	matrix<int> flowers(nodecount,nodecount);

	while (getline(triangulation,next_line))
	{
		
		if (next_line.find("BOUNDARY") != string::npos)
		{
			if (next_line.find("anticlockwise") != string::npos)
				clockwise_boundary = false;
		}		
		else if (next_line.find("FLOWERS") != string::npos)
		{
			getline(triangulation,next_line); // read the empty line
			
			/* now read the flowers */
			while(getline(triangulation,next_line))
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

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: node " << node << ", petal count = " << count << ": " << first_petal;
				
				for (int i=0; i< count; i++)
				{
					iss >> petal;
					flowers[node][i+2] = petal;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " " << petal;					
				}

				if (node < beta-1)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "; interior vertex" << endl;					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "; boundary vertex" << endl;					
				}										
			}		
		}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "flowers" << endl;					
	print(flowers,debug,4,"");
}


    matrix<double> vcoords(nodecount,2);

	if (PLESTENJAK_FORCE_DIRECTION)
	{
		/* set the boundary vertices evenly around the unit circle and initialize
		   all other vertices to sit at the origin.
		*/
	
		for (int i=0; i< num_type123_vertices; i++)
		{
			vcoords[i][0] = 0;
			vcoords[i][1] = 0;
		}
	
	
		/* Standard type 5 vertices have five neighbours.  If the turning cycle bunding the infinite
		   region has length 2, two additonal type 5 vertices have been created, listed at the end of
		   flowers, and boundary vertices have just four neighbours, so we can identify them easily and
		   place them correctly around the boundary circle.  Note that their orientation is the same 
		   regardless of the value of clockwise_boundary.
		*/
		int num_standard_type4_vertices = nodecount-num_type123_vertices;
		bool infinite_cycle_length_2 = false;
		
		if (flowers[nodecount-1][0] == 4)
		{
			infinite_cycle_length_2 = true;
			num_standard_type4_vertices -= 2;
		}
		
		for (int i=0; i< num_standard_type4_vertices; i++)
		{
			double theta = (clockwise_boundary? num_standard_type4_vertices - i -1 : i)*2*pi/(num_standard_type4_vertices);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: boundary vertex " << i << ", theta = " << theta << endl;
	
			vcoords[i+num_type123_vertices][0] = cos(theta);
			vcoords[i+num_type123_vertices][1] = sin(theta);
		}
		
		if (infinite_cycle_length_2)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: setting additional boundary vertices" << endl;
	
			vcoords[nodecount-2][0] = 0;
			vcoords[nodecount-2][1] = 1.0;
			vcoords[nodecount-1][0] = 0;
			vcoords[nodecount-1][1] = -1.0;
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement:   vcoords[" << nodecount-2 << "] = (" << vcoords[nodecount-2][0] << "," << vcoords[nodecount-2][1] << ")" << endl;
	debug << "force_directed_placement:   vcoords[" << nodecount-1 << "] = (" << vcoords[nodecount-1][0] << "," << vcoords[nodecount-1][1] << ")" << endl;
}

		}
   
	}
	else
	{
		/* read the initial vertex placements from __circlepack.out */
		ifstream vertices;
		
		vertices.open("__circlepack.out");
		if (!vertices)
		{
			cout << "\nError opening output file __circlepack.out\n";
			exit(0);
	    }
	    
		int loc_num_vertices;
		vertices >> loc_num_vertices;
		
		for (int i=0; i< nodecount; i++)
		{
			for (int j=0; j< 2; j++)
				vertices >> vcoords[i][j];
		}

//	for (int i=0; i< nodecount; i++)
//	{
//		for (int j=0; j< 2; j++)
//			vcoords[i][j] *= force_directed_coordinate_scale_factor;
//	}
	
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement: vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "force_directed_placement:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
}


		/* if, required, reduce the scale of the type 5 vertices */
		if (RETRACT_BOUNDARY_VERTICES)
		{
			double boundary_cog_x=0;
			double boundary_cog_y=0;		
			for (int v=num_type123_vertices; v< nodecount; v++)
			{
				boundary_cog_x += vcoords[v][0];
				boundary_cog_y += vcoords[v][1];
			}
			boundary_cog_x /= (nodecount - num_type123_vertices);
			boundary_cog_y /= (nodecount - num_type123_vertices);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: boundary vertex centre of gravity = (" << boundary_cog_x << "," << boundary_cog_y << ")" << endl;
	
			for (int i=num_type123_vertices; i< nodecount; i++)
			{
				vcoords[i][0] = (vcoords[i][0]-boundary_cog_x)*boundary_vertex_retraction_factor + boundary_cog_x;
				vcoords[i][1] = (vcoords[i][1]-boundary_cog_y)*boundary_vertex_retraction_factor + boundary_cog_y;
			}
		}
	}

    matrix<double> initial_vcoords = vcoords;
	matrix<double> two_back_vcoords = vcoords;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "force_directed_placement: initial_vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "force_directed_placement:   v" << i << " = ("<< initial_vcoords[i][0] << ", " << initial_vcoords[i][1] << ")" << endl;
}

	int num_vertices;
	
	/* identify the number of vertices we're actually interested in applying forces to */
	if (INCLUDE_BOUNDARY_VERTICES)
	{
		num_vertices = nodecount;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: including boundary vertices in force directed placement, num_vertices = " << num_vertices << endl;
	}
	else
	{
		num_vertices = num_type123_vertices;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: no boundary vertices in force directed placement, num_vertices = " << num_vertices << endl;
	}
	
	/* set the k_constant, which is used in the default attractive and repulsive force functions.
	   as before we consider only the type 1, 2, and 3 vertices 
	*/
	
/* 7/6/18 NODECOUNT SHOULD PROBALBY BE NUM_VERTICES */

	k_constant = sqrt(nodecount/pi);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: k_constant  = " << k_constant << endl;


	/* Main iterative loop */	   
	for (int i=0; i< num_iterations; i++)
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement: iteration " << i << endl;

		/* determine the minimum and average length of an edge in the triangulation	*/
		double average_triangulation_edge_length = 0;
		double minimum_triangulation_edge_length = 0;
		
		int num_edges = 0;
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "force_directed_placement: evaluating the average triangulation edge length" << endl;

		for (int u=0; u < (!PLESTENJAK_FORCE_DIRECTION && APPLY_FORCES_TO_TYPE_12_ONLY?num_type12_vertices:num_vertices); u++)
		{			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "force_directed_placement:   considering the flower of vertex " << u  << endl;
	
			/* there is always one extra vertex recorded in the rows of flower, if it's a type 1, 2, 
			   or 3 vertex, there is a repeated petal, if it is a type 5 vertex, the vertex itself
			   appears at the end of the row.
			*/
			for (int j=1; j <= flowers[u][0]; j++)
			{				
				int v = flowers[u][j] - 1; // flowers records vertices numbered from 1

				/* ignore the edge if it connects to the wrong type of vertex, we only want to consider edges joining 
				   type 1 and type 2 vertices.  Note that u is known to be type 1 or 2 and that type 1 vertices are
				   not joined by an edge, so we only need check for when u is joined to a type 3 or 5 vertex, and when
				   u is type 2 and joined to another type 2.
				*/
				if (!PLESTENJAK_FORCE_DIRECTION && APPLY_FORCES_TO_TYPE_12_ONLY)
				{
					if (v >= num_type12_vertices)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "force_directed_placement:     vertex " << v << " is a type 3 or 5 vertex, ignoring" << endl;
						continue; // type 3 or 5 vertex						
					}
					else if (u >= num_crossings && v >= num_crossings) 
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "force_directed_placement:     both vertices " << u << " and " << v << " are type 2, ignoring " << v << endl;
						continue; // u is type 2 so v must be type 1
					}					
				}
				else if (v >= num_type123_vertices && !INCLUDE_BOUNDARY_VERTICES)
					continue; // type 5 vertex

				double delta_x = vcoords[v][0] - vcoords[u][0];
				double delta_y = vcoords[v][1] - vcoords[u][1];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "force_directed_placement:     vertices " << u  << ", " << v << endl;
	debug << "force_directed_placement:       v = (" << vcoords[v][0] << "," << vcoords[v][1] << ")" << endl;
	debug << "force_directed_placement:       u = (" << vcoords[u][0] << "," << vcoords[u][1] << ")" << endl;
	debug << "force_directed_placement:     delta = (" << delta_x << "," << delta_y << ")" << endl;
}
			
				double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);
				
				if (minimum_triangulation_edge_length)
				{
					if (minimum_triangulation_edge_length > mod_delta)
						minimum_triangulation_edge_length = mod_delta;
				}
				else
					minimum_triangulation_edge_length = mod_delta;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "force_directed_placement:     mod_delta =  " << mod_delta << endl;
	
				average_triangulation_edge_length += mod_delta;
				num_edges++;

			}
		}

		/* The above loop has counted every edge twice and accumulated
		   every edge length twice, so the factor of two cancels out 
		*/
		average_triangulation_edge_length /= num_edges;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement: average_triangulation_edge_length =  " << average_triangulation_edge_length << endl;
	debug << "force_directed_placement: average_triangulation_length_factor =  " << average_triangulation_length_factor << endl;
	debug << "force_directed_placement: average_triangulation_edge_length* average_triangulation_length_factor =  " << average_triangulation_edge_length*average_triangulation_length_factor << endl;
	debug << "force_directed_placement: minimum_triangulation_edge_length =  " << minimum_triangulation_edge_length << endl;
	debug << "force_directed_placement: applying attractive and repulsive forces to vertices:" << endl;
}

		matrix<double> displacement(nodecount,2);

		/* temperature decreases with i */
		double temperature = sqrt(pi/nodecount)/(1 + pow(i,1.5)*pi/nodecount);
		if (debug_control::DEBUG >= debug_control::SUMMARY)
			debug << "force_directed_placement:  operating temperature  = " << temperature << endl;

		/* Calculate the forces between vertices joined by an edge.
		   We calculate the forces on a vertex u from the vertices in its flower.
		*/
		
		for (int u=0; u < (!PLESTENJAK_FORCE_DIRECTION && APPLY_FORCES_TO_TYPE_12_ONLY?num_type12_vertices:num_type123_vertices); u++)
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "force_directed_placement:   considering the flower of vertex " << u  << endl;

			for (int j=1; j <= flowers[u][0]; j++)
			{				
				int v = flowers[u][j] - 1; // flowers records vertices numbered from 1

				/* ignore the edge if it connects to the wrong type of vertex, we only want to consider edges joining 
				   type 1 and type 2 vertices.  Note that u is known to be type 1 or 2 and that type 1 vertices are
				   not joined by an edge, so we only need check for when u is joined to a type 3 or 5 vertex, and when
				   u is type 2 and joined to another type 2.
				*/
				if (!PLESTENJAK_FORCE_DIRECTION && APPLY_FORCES_TO_TYPE_12_ONLY)
				{
					if (v >= num_type12_vertices)
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "force_directed_placement:     vertex " << v << " is a type 3 or 5 vertex, ignoring" << endl;
						continue; // type 3 or 5 vertex						
					}
					else if (u >= num_crossings && v >= num_crossings) 
					{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "force_directed_placement:     both vertices " << u << " and " << v << " are type 2, ignoring " << v << endl;
						continue; // u is type 2 so v must be type 1
					}
				}
				else if (!PLESTENJAK_FORCE_DIRECTION && !INCLUDE_BOUNDARY_VERTICES && v >= num_type123_vertices )
					continue; // type 5 vertex



				
				double delta_x = vcoords[v][0] - vcoords[u][0];
				double delta_y = vcoords[v][1] - vcoords[u][1];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "force_directed_placement:     vertices " << u  << ", " << v << endl;
	debug << "force_directed_placement:       v = (" << vcoords[v][0] << "," << vcoords[v][1] << ")" << endl;
	debug << "force_directed_placement:       u = (" << vcoords[u][0] << "," << vcoords[u][1] << ")" << endl;
	debug << "force_directed_placement:     delta = (" << delta_x << "," << delta_y << ")" << endl;
}
			
				double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);


				double force;
				
				if (PLESTENJAK_FORCE_DIRECTION)
				{
					force = attractive_force(mod_delta);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement:     mod_delta =  " << mod_delta << endl;
	debug << "force_directed_placement:     vertices " << u  << ", " << v << " attract with force " << force << endl;
}

				}
				else
				{
					if ( mod_delta < average_triangulation_edge_length * average_triangulation_length_factor)
//						force = repulsive_force(mod_delta);
						force = repulsive_force(average_triangulation_edge_length * average_triangulation_length_factor - mod_delta);
					else
						force = 0;

if (force && debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement:     mod_delta =  " << mod_delta << endl;
	debug << "force_directed_placement:     force based on length =  " << average_triangulation_edge_length * average_triangulation_length_factor - mod_delta << endl;
	debug << "force_directed_placement:     vertex " << u  << " repelled with force " << force << " by vertex " << v << endl;
}
				}

						
				double disp_x;
				double disp_y;
				
				disp_x = delta_x*force;
				disp_y = delta_y*force;

if (force && debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "force_directed_placement:   (" << u << " -> " << v << ") disp_x = " << disp_x << endl;
	debug << "force_directed_placement:   (" << u << " -> " << v << ") disp_y = " << disp_y << endl;
}

				if (PLESTENJAK_FORCE_DIRECTION)
				{
					/* attract vertices towards each other */
					displacement[u][0] += disp_x;
					displacement[u][1] += disp_y;
				}
				else
				{
					/* repulse vertices that are too close */
					displacement[u][0] -= disp_x;
					displacement[u][1] -= disp_y;
				}

if (force && debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_directed_placement:   displacement[" << u << "] = (" << displacement[u][0] << "," << displacement[u][1] << ")" << endl;

			}
		}		

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement: final displacements: " << endl;
	for (int j=0; j< nodecount; j++)
		debug << "force_directed_placement:   vertex " << j << ": " << displacement[j][0] << " " << displacement[j][1] << endl;
}
			
		/* Apply the displacements to the vertex positions, limiting the maximum
		   displacement to that determined by the temperature.
		*/
		for (int j=0; j< num_type123_vertices; j++)
		{			
			double displacement_length = sqrt(displacement[j][0]*displacement[j][0]+displacement[j][1]*displacement[j][1]);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement:   adjusting vertex " << j << endl; 
	debug << "force_directed_placement:     displacement length = " << displacement_length 
	      << (temperature < displacement_length ? " bounded" : " not bounded") << " by temperature = " << temperature << endl;
	debug << "force_directed_placement:     displacement factor = " << min(displacement_length, temperature)/displacement_length << endl;
}

			if (displacement_length != 0)
			{
				displacement[j][0] *= min(displacement_length, temperature)/displacement_length;
				displacement[j][1] *= min(displacement_length, temperature)/displacement_length;
				vcoords[j][0] += displacement[j][0]; 
				vcoords[j][1] += displacement[j][1];
			}
		}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "force_directed_placement: temperature adjusted final displacements: " << endl;
	for (int j=0; j< nodecount; j++)
		debug << "force_directed_placement:   vertex " << j << ": " << displacement[j][0] << " " << displacement[j][1] << endl;
		
	debug << "force_directed_placement: updated vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "force_directed_placement:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
}

		if (TRACK_PLACEMENT_ITERATION && i % placement_iteration_tracking_step == placement_iteration_tracking_step-1)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_direction: tracking force directed placement at iterative step i = " << i << endl;

			matrix<double> write_metapost_vcoords(2*nodecount,2);
			for (int j=0; j< nodecount; j++)
			{
				for (int k=0; k< 2; k++)
					write_metapost_vcoords[j][k] = vcoords[j][k];
			}
			for (int j=0; j< nodecount; j++)
			{
				for (int k=0; k< 2; k++)
					write_metapost_vcoords[j+nodecount][k] = initial_vcoords[j][k];
			}
			
			vector<double> vertex_radius(nodecount); // this will be initialized to zero, it will be unused by write_metapost.
			write_metapost(output, code_data, title, mp_control, write_metapost_vcoords, vertex_radius);	
		}

		
		/* check for oscillatory convergence.  Look for the maximum delta between the vertex coordinates
		   now and two steps back in the iteration.  If this is less than EISPILON, stop the iteration.
		*/
		if (i > 2)
		{
			double max_delta = abs(vcoords[0][0] - two_back_vcoords[0][0]);
			
			if (abs(vcoords[0][1] - two_back_vcoords[0][1]) > max_delta)
				max_delta = abs(vcoords[0][1] - two_back_vcoords[0][1]);
			
			for (int j=1; j< nodecount; j++)
			{
				for (int k=0; k< 2; k++)
				{
					if (abs(vcoords[j][k] - two_back_vcoords[j][k]) > max_delta)
						max_delta = abs(vcoords[j][k] - two_back_vcoords[j][k]);
				}
			}
			
			if (max_delta < EPSILON)
			{
				cout << "Force directed placement stopping after " << i << " iterations, oscillation convergence detected" << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "Force directed placement stopping after " << i << " iterations, oscillation convergence detected" << endl;
				
				break;
			}
		}

		two_back_vcoords = initial_vcoords;
		initial_vcoords = vcoords;
		
	}

	/* write the new location of the vertices to the output file */
	ofstream output(circlepack_output_file);
	
	if (!output)
	{
		cout << "\nError opening output file\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "force_direction: could not open " << circlepack_output_file << endl;
		exit(0);
	}		
	
	output << nodecount << endl;
	for (int i=0; i< nodecount; i++)
	{
		output << setprecision(6) << setw(15) << vcoords[i][0] << " "		
		       << setprecision(6) << setw(15) << vcoords[i][1] << endl;		
	}

	/* write the radii as zero, since these are not used 
	   with force direction, since mp_control.circle_packing == false. 
	*/
	for (int i=0; i< nodecount; i++)
		output << setprecision(6) << setw(10) << 0 << endl;		
	
	output.close();
}

void centre_of_gravity_placement(metapost_control& mp_control, generic_code_data& code_data, int num_iterations, string title)
{

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement: code_data: ";
	write_code_data(debug,code_data);
	debug << endl;
	print_code_data(code_data, debug, "centre_of_gravity_placement: ");
	debug << "centre_of_gravity_placement: num_iterations = " << num_iterations << endl;
	if (TRACK_PLACEMENT_ITERATION)
		debug << "centre_of_gravity_placement: placement_iteration_tracking_step = " << placement_iteration_tracking_step << endl;
}
	/* turn off unsupported options */
	DRAW_IN_HYPERBOLIC_DISC = false;
	mp_control.circle_packing = false;

	/* use circle packing to provide an initial placement of the vertices */
	if (USE_KEN_STEPHENSON_CIRCLE_PACKING)
		KS_circle_pack (triangulation_output_file, circlepack_output_file);
	else
		circle_pack (triangulation_output_file, circlepack_output_file);

	/* record the frame corners based on the final choice of infinite cycle */
	set_frame_corners(circlepack_output_file,mp_control); 

	/* if we're tracking the placement iteration, draw the starting point but don't attempt to
	   draw any comparative triangulation, because there isn't one yet!
	*/
	if (TRACK_PLACEMENT_ITERATION)
	{
		TRACK_PLACEMENT_ITERATION = false;
		write_metapost(output, code_data, title, mp_control, circlepack_output_file, circlepack_output_file);		
		TRACK_PLACEMENT_ITERATION = true;
	}

	/* Read triangulation data from triangulation_output_file */
	int nodecount;
	int alpha;
	int beta;
	int gamma;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "centre_of_gravity_placement: read triangulation data from " << triangulation_output_file << endl;

	ifstream triangulation(triangulation_output_file);
	if(!triangulation)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: could not open " << triangulation_output_file << endl;
		exit(0);
	}		
	
	string next_line;
	getline(triangulation,next_line);
	if (next_line.find("NODECOUNT") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> nodecount;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "centre_of_gravity_placement: nodecount = " << nodecount << endl;
	}

	getline(triangulation,next_line);
	if (next_line.find("ALPHA") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> alpha; //always 1, the first vertex
		iss >> beta; //the first boundary (type 5) vertex
		iss >> gamma; // IS THIS USED?
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "centre_of_gravity_placement: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
	}
	
	int num_type123_vertices = beta-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "centre_of_gravity_placement: number of type 1,2 & 3 vertices = " << num_type123_vertices << endl;

	/* Store the flowers for each node in a matrix with the first column set to the number of neighbours
	   of that vertex.  No vertex is adjacent to itself, nor is any vertex adjacent to all other vertices
	   so the maximum number of neighbours is <= nodecount -2.  Recall that for interior vertices the first
	   petal appears at the end of the list as well, to simplify the calculation of angle sums.
	*/
	matrix<int> flowers(nodecount,nodecount);

	while (getline(triangulation,next_line))
	{
		if (next_line.find("FLOWERS") != string::npos)
		{
			getline(triangulation,next_line); // read the empty line
			
			/* now read the flowers */
			while(getline(triangulation,next_line))
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

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: node " << node << ", petal count = " << count << ": " << first_petal;
				
				for (int i=0; i< count; i++)
				{
					iss >> petal;
					flowers[node][i+2] = petal;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " " << petal;					
				}

				if (node < beta-1)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "; interior vertex" << endl;					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "; boundary vertex" << endl;					
				}										
			}		
		}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "flowers" << endl;					
	print(flowers,debug,4,"");
}

	/* read the initial vertex placements from circlepack_output_file */
	ifstream vertices;
	
	vertices.open(circlepack_output_file);
	if (!vertices)
	{
		cout << "\nError opening output file " << circlepack_output_file << endl;
		exit(0);
    }
    
	int num_crossings = code_data.num_crossings;
	int num_edges = 2 * num_crossings;
	int num_vertices;
	vertices >> num_vertices;
	
    matrix<double> vcoords(nodecount,2);

	for (int i=0; i< nodecount; i++)
	{
		for (int j=0; j< 2; j++)
			vertices >> vcoords[i][j];
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement: vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "centre_of_gravity_placement:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
}


	/* if, required, reduce the scale of the type 5 vertices */
	if (RETRACT_BOUNDARY_VERTICES)
	{
		double boundary_cog_x=0;
		double boundary_cog_y=0;		
		for (int v=num_type123_vertices; v< nodecount; v++)
		{
			boundary_cog_x += vcoords[v][0];
			boundary_cog_y += vcoords[v][1];
		}
		boundary_cog_x /= (nodecount - num_type123_vertices);
		boundary_cog_y /= (nodecount - num_type123_vertices);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: boundary vertex centre of gravity = (" << boundary_cog_x << "," << boundary_cog_y << ")" << endl;
	
		for (int i=num_type123_vertices; i< nodecount; i++)
		{	
			vcoords[i][0] = (vcoords[i][0]-boundary_cog_x)*boundary_vertex_retraction_factor + boundary_cog_x;
			vcoords[i][1] = (vcoords[i][1]-boundary_cog_y)*boundary_vertex_retraction_factor + boundary_cog_y;
		}
	}
  
    matrix<double> initial_vcoords = vcoords;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "centre_of_gravity_placement: initial_vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "centre_of_gravity_placement:   v" << i << " = ("<< initial_vcoords[i][0] << ", " << initial_vcoords[i][1] << ")" << endl;
}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement: num_crossings = " << num_crossings << endl;
	debug << "centre_of_gravity_placement: num_edges = " << num_edges << endl;
}

	/* identify the number of vertices we're actually interested in applying forces to */
	if (INCLUDE_BOUNDARY_VERTICES)
	{
		num_vertices = nodecount;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: including boundary vertices in force directed placement, num_vertices = " << num_vertices << endl;
	}
	else
	{
		num_vertices = num_type123_vertices;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: no boundary vertices in force directed placement, num_vertices = " << num_vertices << endl;
	}
	   
	
	/* We will write inner hull vertices to a file for testing purposes*/
	ofstream IHoutput;
	IHoutput.open("__inner_hull.out");
	if (!IHoutput)
	{
		cout << "\nError opening output file __inner_hull.out\n";
		exit(0);
    }

	/*  Determing the radius of the boundary vertex circle. This is used if 
		we are considering boundary circles.
	*/

	double boundary_cog_x=0;
	double boundary_cog_y=0;
	double boundary_vertex_circle_radius = 0;
	
	if (INCLUDE_BOUNDARY_VERTICES)
	{			
		for (int v=num_type123_vertices; v< num_vertices; v++)
		{
			boundary_cog_x += vcoords[v][0];
			boundary_cog_y += vcoords[v][1];
		}
		boundary_cog_x /= (num_vertices - num_type123_vertices);
		boundary_cog_y /= (num_vertices - num_type123_vertices);
		
		for (int v=num_type123_vertices; v< num_vertices; v++)
		{
			double rv_x = vcoords[v][0] - boundary_cog_x;
			double rv_y = vcoords[v][1] - boundary_cog_y;
		    boundary_vertex_circle_radius += sqrt (rv_x*rv_x+rv_y*rv_y);
		}
		boundary_vertex_circle_radius /= (num_vertices - num_type123_vertices);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement: boundary centre of gravity  = (" << boundary_cog_x << "," << boundary_cog_y << 
	      "), boundary radius = "<<	boundary_vertex_circle_radius << endl;
}
	}
	else
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: not considering boundary vertices, boundary centre of gravity not evaluated" << endl;
	}

	/* write the boundary circle centre and radius to the inner hull output file */
	IHoutput << boundary_cog_x << " " << boundary_cog_y << " " << boundary_vertex_circle_radius << endl;

	
	/* Main iterative loop */	   
	for (int i=0; i< num_iterations; i++)
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: iteration " << i << endl;

		matrix<double> displacement(nodecount,2);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: applying cente of gravity placement" << endl;


		/* In centre of gravity placement we evaluate the inner hull of the closed star of each vertex.  For
		   interior vertices we move the vertex to the centre of gravity of the inner hull.  Since the inner hull
		   of boundary vertices have as a vertex the boundary vertex itself, moving that vertex to the centre of
		   gravity of the inner hull would compress the diagram.  We therefore evaluate the centre of gravity of
		   the boundary vertices and the average distance of the boundary vertices from that centre.  We then
		   project the two inner hull vertices adjacent to the boundary vertex onto a circle, centred at the 
		   boundary vertex centroid with radius this average distance, and move the boundary vertex to the midpoint
		   on the circle between these two projections.  In the first iteration, the circle approximates a circle 
		   on which the boundary vertices lie, in subsequent iterations it is exact, of course.
		*/

		matrix<double> cog_coords(nodecount,2);
		for (int j=0; j< num_vertices; j++)
		{

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement: vertex " << j << endl;
	debug << "centre_of_gravity_placement:   " << ( j < num_type123_vertices ? "interior" : "boundary") << " vertex" << endl;
}


			/* evaluate the inner hull of the closed star of vertex j, first calculate the min and max
			   x and y coordinates of the vertices in the link of vertex j;
			*/
			double minx = vcoords[flowers[j][1]-1][0];
			double miny = vcoords[flowers[j][1]-1][1];
			double maxx = vcoords[flowers[j][1]-1][0];
			double maxy = vcoords[flowers[j][1]-1][1];
			
			for (int k=1; k < (j < num_type123_vertices ? flowers[j][0]: flowers[j][0]+1); k++)
			{
				minx = min(minx, vcoords[flowers[j][k+1]-1][0]);
				maxx = max(maxx, vcoords[flowers[j][k+1]-1][0]);
				miny = min(miny, vcoords[flowers[j][k+1]-1][1]);
				maxy = max(maxy, vcoords[flowers[j][k+1]-1][1]);
			}
							
			/* Ukl will denote the intersection between the line joining vetices k and k+1 and
			   the line joining vertices l and l+1
			*/
			double Ukl_x;
			double Ukl_y;
			
			/* The maximum number of inner hull vertices is the number of vertices in the link of 
			   vertex j plus the number of pairwise intersections between the edges in the link.
			   If there are n vertices in the link, there are n edges, so the maximum number of
			   vertices in the inner hull is n + n(n-1)/2 = n(n+1)/2
			*/
			matrix<double> inner_hull_vertex(flowers[j][0]*(flowers[j][0]+1)/2, 2);
			int num_IH_vertices = 0;
			int boundary_vertex_IH_index = 0;  // used to identify the boundary vertex within its own innner hull
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:   supporting line intersections:" << endl;

				
			for (int k=0; k < (j < num_type123_vertices ? flowers[j][0]-1: flowers[j][0]); k++)
			{
				for (int l= k+1; l < (j < num_type123_vertices? flowers[j][0]: flowers[j][0]+1); l++)
				{
					if (j < num_type123_vertices && k == 0 && l == flowers[j][0]-1)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:   k = " << k << ", l = " << l << endl;
	debug << "centre_of_gravity_placement:   from " << flowers[j][1]-1 << " and from " << flowers[j][l+1]-1 << endl;
	debug << "centre_of_gravity_placement:   intersection at " << flowers[j][1]-1 << endl;
}

						Ukl_x = vcoords[flowers[j][1]-1][0];
						Ukl_y = vcoords[flowers[j][1]-1][1];
					}
					else if (j >= num_type123_vertices && k == 0 && l == flowers[j][0])
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:   k = " << k << ", l = " << l << endl;
	debug << "centre_of_gravity_placement:   from " << flowers[j][1]-1 << " and from " << flowers[j][l+1]-1 << endl;
	debug << "centre_of_gravity_placement:   intersection at " << flowers[j][1]-1 << endl;
}

						Ukl_x = vcoords[flowers[j][1]-1][0];
						Ukl_y = vcoords[flowers[j][1]-1][1];
					}
					else if (l == k+1)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:   k = " << k << ", l = " << l << endl;
	debug << "centre_of_gravity_placement:   from " << flowers[j][k+1]-1 << " and from " << flowers[j][l+1]-1 << endl;
	debug << "centre_of_gravity_placement:   intersection at " << flowers[j][l+1]-1 << endl;
}
						Ukl_x = vcoords[flowers[j][l+1]-1][0];
						Ukl_y = vcoords[flowers[j][l+1]-1][1];
						
						if (j >= num_type123_vertices && flowers[j][l+1]-1 == j)
						{
							/* we know this is going to be an inner hull vertex */
							boundary_vertex_IH_index = num_IH_vertices;
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:   boundary_vertex_IH_index = " << boundary_vertex_IH_index << endl;
						}
					}
					else
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:   k = " << k << ", l = " << l << endl;
	debug << "centre_of_gravity_placement:   from " << flowers[j][k+1]-1 << " and from " << flowers[j][l+1]-1 << endl;
}

						double Vk_x = vcoords[flowers[j][k+1]-1][0];
						double Vk_y = vcoords[flowers[j][k+1]-1][1];
						double Vkk_x = vcoords[flowers[j][k+2]-1][0];
						double Vkk_y = vcoords[flowers[j][k+2]-1][1];
						double Vl_x = vcoords[flowers[j][l+1]-1][0];
						double Vl_y = vcoords[flowers[j][l+1]-1][1];
						
//							double Vll_x = vcoords[flowers[j][(j < num_type123_vertices && l < 5? l+2: 1)]-1][0];
//							double Vll_y = vcoords[flowers[j][(j < num_type123_vertices && l < 5? l+2: 1)]-1][1];

						double Vll_x;
						if (j >= num_type123_vertices && l == 5)
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:     boundary vertex, terminal petal" << endl;							
							Vll_x = vcoords[flowers[j][1]-1][0];
						}
						else
						{
							Vll_x = vcoords[flowers[j][l+2]-1][0];
						}
						
						double Vll_y;
						if (j >= num_type123_vertices && l == 5)
							Vll_y = vcoords[flowers[j][1]-1][1];
						else
							Vll_y = vcoords[flowers[j][l+2]-1][1];

						
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:     Vk_x = " << Vk_x << ", Vk_y = " << Vk_y << ", Vkk_x = " << Vkk_x << ", Vkk_y = " << Vkk_y << endl;							
	debug << "centre_of_gravity_placement:     Vl_x = " << Vl_x << ", Vl_y = " << Vl_y << ", Vll_x = " << Vll_x << ", Vll_y = " << Vll_y << endl;							
}	
						
			            /* check whether the lines joining k, k+1 and l,l+1 are parallel: compare gradients */
			            double gk = (Vkk_y - Vk_y)/(Vkk_x - Vk_x); 
			            double gl = (Vll_y - Vl_y)/(Vll_x - Vl_x); 
		
			            double delta = gk - gl;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     gradient gk = " << gk << ", gl = " << gl << ", delta = " << delta << endl;
				
			            if (abs(delta) <= ZERO_THRESHOLD)
			            {
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     considered parallel, continue to consider next line" << endl;
							continue; // to next l
						}
				            
			            /* lines not parallel, evaluate intersection 
			            
			               The line joining vertices k and k+1 (Vk and Vkk) is given by y = gk (x - Vk_x) + Vk_y
			               and the line joining l and l+1 is y = gl (x - Vl_x) + Vl_y.  Writing these equations as
			               
			               gk x - y = gk*Vk_x - Vk_y
			               gl x - y = gl*Vl_x - Vl_y
			               
			               we have two simultaneous equations
			               
			               (a  b)   (x)   =  (u)
			               (c  d)   (y)   =  (v)
			               
			               whose solution is 
			               
			               (x)        1      ( d   -b) (u)          1    (v - u)
			                     =   ----                     =    ----          , since b=d=-1.
			               (y)      (ad-bc)  (-c    a) (v)         (c-a) (av-cu)					     					               
			            */
						double u = (gk*Vk_x - Vk_y);
						double v = (gl*Vl_x - Vl_y);
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     u= " << u << ", v = " << v << endl;
							Ukl_x = (v-u)/(gl-gk);
						Ukl_y = (gk*v-gl*u)/(gl-gk);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:   intersection at (" << Ukl_x << "," << Ukl_y << ")" << endl;
					}

			        /* Check whether Uij is a defining point of the inner hull */
			
			        /* If Ukl is a point of the inner hull it must at least lie within the 
			           min and max bounds for x and y.  We needed this check to stay within
			           Metapost's arithmetic bounds, so it seems reasonable to keep it here
			        */
			        if (Ukl_x < minx || Ukl_x > maxx || Ukl_y < miny || Ukl_y > maxy)
			        {
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:   intersection outside min-max bounds, continue to consider next line" << endl;
			            continue;  // to next l
			        }

			        /* Test Ukl against all support lines.  If the j-th vertex is a type 1, 2, or 3
			           vertex we know it lies in the interior of the inner hull and therefore reject 
			           Ukl if it lies on the opposite side of a support line to the j-th vertex.
			            
			           If the j-th vertex is a type 5 vertex, it lies within two of the support lines, 
			           joining the vertex to two other type 5 vertices.  For these lines we test against 
			           the intersection against the type 5 vertex NOT contained in the line.
			           
			        */
			        bool accept_Ukl = true;					
			        
		            for (int v = 0; v < (j < num_type123_vertices ? flowers[j][0]: flowers[j][0]+1); v++)
		            {
						int a = flowers[j][v+1]-1;
						int b;
						if (j >= num_type123_vertices && v == 5)
							b = flowers[j][1]-1;
						else
							b = flowers[j][v+2]-1;
							
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "centre_of_gravity_placement:   test intersection against line through " << a << " and " << b << endl;
							
						double alpha;
		                double g = (vcoords[b][1] - vcoords[a][1])/(vcoords[b][0] - vcoords[a][0]);
						
						if (j < num_type123_vertices)
						{								
							/* determine on which side of the line joining the a-th and b-th vertex 
							   the j-th vertex lies 
							*/									
							alpha = vcoords[j][1] - g*(vcoords[j][0] - vcoords[a][0]) - vcoords[a][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     reference vertex is " << j << endl;
						}
						else
						{
							if (v == 4)
							{
								/* reference vertex is flowers[j][1] */
								alpha = vcoords[flowers[j][1]-1][1] - g*(vcoords[flowers[j][1]-1][0] - vcoords[a][0]) - vcoords[a][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     reference vertex is " << flowers[j][1]-1 << endl;
							}
							else if (v == 5)
							{
								/* reference vertex is flowers[j][5] */
								alpha = vcoords[flowers[j][5]-1][1] - g*(vcoords[flowers[j][5]-1][0] - vcoords[a][0]) - vcoords[a][1];
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     reference vertex is " << flowers[j][5]-1 << endl;
							}
							else
							{
								/* reference vertex is the jth vertex */
								alpha = vcoords[j][1] - g*(vcoords[j][0] - vcoords[a][0]) - vcoords[a][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     reference vertex is " << j << endl;
							}
						}
							
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     reference vertex alpha value = " << alpha << endl;

		                /* we only need the sign and keep the numbers small to avoid arithmetic overflow */
		                if (abs(alpha) < ZERO_THRESHOLD) 
		                    alpha = 0; 
		                			
		                if (alpha !=0) 
		                    alpha /= abs(alpha); 
		                
		
		                /* determine on which side of the line joining the a-th and b-th vertex Ukl lies */
		                double beta = Ukl_y - g*(Ukl_x - vcoords[a][0]) - vcoords[a][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "centre_of_gravity_placement:     intersection's beta value = " << beta << endl;
			                
		                if (abs(beta) < ZERO_THRESHOLD) 
		                    beta = 0; 
		                			
		                if (beta !=0) 
		                    beta /= abs(beta); 
		
		                if (alpha * beta < 0)
		                {
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:     intersection lies on wrong side of the line to be an inner hull vertex" << endl;
		                    accept_Ukl = false;
		                    break;
		                }
		            }
			            
		            if (accept_Ukl)
		            {
			            inner_hull_vertex[num_IH_vertices][0] = Ukl_x;
			            inner_hull_vertex[num_IH_vertices][1] = Ukl_y;
			            num_IH_vertices++;
				            
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:   intersection is an inner hull vertex" << endl;
					}
				}
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:   number of inner hull vertices = " << num_IH_vertices << endl;
	debug << "centre_of_gravity_placement:   inner hull vertices:" << endl;
	for (int i=0; i< num_IH_vertices; i++)
		debug << "centre_of_gravity_placement:   (" << inner_hull_vertex[i][0] << "," << inner_hull_vertex[i][1] << ")" << endl;
}

			double cog_x = 0;
			double cog_y = 0;

			for (int v=0; v< num_IH_vertices; v++)
			{
				cog_x += inner_hull_vertex[v][0];
				cog_y += inner_hull_vertex[v][1];
			}
			
			cog_x /= num_IH_vertices;				
			cog_y /= num_IH_vertices;					

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:   centre of gravity = (" << cog_x << "," << cog_y << ")" << endl;

				
			/* We order the inner hull vertices for by looking at their argument around the centre of 
			   gravity of the inner hull (this approach works for both interior and boundary vertices).  				   
			*/				
			vector<double> theta(num_IH_vertices);
			
			for (int v = 0; v < num_IH_vertices; v++)
			{
				theta[v] = atan2(inner_hull_vertex[v][1]-cog_y,inner_hull_vertex[v][0]-cog_x);

				if (theta[v] < 0)
				    theta[v] += 2*3.14159265358979;

				theta[v] += 1;
			}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "centre_of_gravity_placement:   inner hull vertex theta values: ";
	for (int i=0; i< num_IH_vertices; i++)
		debug << theta[i] << " ";
	debug << endl;
}
				
			matrix<double> ordered_inner_hull_vertex(num_IH_vertices, 2);
			
			/* we only want to re-evaluate the boundary_vertex_IH_index once in the case of boundary vertices. */
			bool reordered_boundary_vertex_IH_index = false;
			
		    /* We select the maximum argument first, setting to zero those arguments we've already selected.
		       Due to this approach we must handle the case where one of the inner hull vertices
		       has argument zero.  To do this we simply increment all arguments by one, 
			   since we are only interested in their order, not their value.  Note that the angle
			   function returns a value in the range (-180,180] so we must adjust this to [0,360) 
			   before incrementing to avoid introducing a value of zero in the case that angle returns -1.
			*/
			for (int u=num_IH_vertices-1; u >= 0; u--)
			{
			    double maxt = theta[0];
			    int index = 0;
			    
			    for (int v=1; v < num_IH_vertices; v++)
			    {
			        if (theta[v] > maxt)
			        {
			            maxt = theta[v];
			            index=v;
			        }
			    }

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "centre_of_gravity_placement:   largest theta value found at " << index << endl;
				
			    ordered_inner_hull_vertex[u][0] = inner_hull_vertex[index][0];
			    ordered_inner_hull_vertex[u][1] = inner_hull_vertex[index][1];
			    theta[index] =0;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "centre_of_gravity_placement:   inner hull vertex theta values: ";
	for (int i=0; i< num_IH_vertices; i++)
		debug << theta[i] << " ";
	debug << endl;
}

				if (j >= num_type123_vertices && !reordered_boundary_vertex_IH_index && index == boundary_vertex_IH_index)
				{
					boundary_vertex_IH_index = u;
					reordered_boundary_vertex_IH_index = true;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "centre_of_gravity_placement:   boundary_vertex_IH_index re-ordered to " << boundary_vertex_IH_index << endl;
				}	
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:   ordered_inner hull vertices:" << endl;
	for (int i=0; i< num_IH_vertices; i++)
		debug << "centre_of_gravity_placement:   (" << ordered_inner_hull_vertex[i][0] << "," << ordered_inner_hull_vertex[i][1] << ")" << endl;
}
				
			/* For testing and demonstration purposes we write the location of the inner hull vertices
			   from the last iteration to an output file, so it can be read when we write the metapost
			   code.  
			*/
//				if (i==0)
			if (i == num_iterations-1)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:    write inner hull vertices: i = " << i << endl;
				IHoutput << j << " " << num_IH_vertices;
				for (int v=0; v< num_IH_vertices; v++)
				{					
					IHoutput << " " << ordered_inner_hull_vertex[(boundary_vertex_IH_index+v)%num_IH_vertices][0] 
					         << " " << ordered_inner_hull_vertex[(boundary_vertex_IH_index+v)%num_IH_vertices][1];						         
				}
				IHoutput << endl;
			}
			
			/* For interior vertices store coordinates of the centre of gravity of vertex's inner hull in cog_coords.
			   For boundary vertices we project the two inner hull vertices adjacent to the vertex onto a circle 
			   centred at the centroid of the type 5 vertices and of radius the average distance between a type 5 
			   vertex and this centroid. Then store the midpoint on this circle between these two projections in cog_coords.
			*/
			if (j < num_type123_vertices)
			{
				cog_coords[j][0] = cog_x;
				cog_coords[j][1] = cog_y;
			}
			else
			{
				int adj_vertex_1 = (boundary_vertex_IH_index - 1 + num_IH_vertices)%num_IH_vertices;
				int adj_vertex_2 = (boundary_vertex_IH_index + 1)%num_IH_vertices;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement:   inner hull vertices adjacent to boundary vertex " << j << ": (" 
		  << ordered_inner_hull_vertex[adj_vertex_1][0] << "," << ordered_inner_hull_vertex[adj_vertex_1][1] << ") and (" 
		  << ordered_inner_hull_vertex[adj_vertex_2][0] << "," << ordered_inner_hull_vertex[adj_vertex_2][1] << ")" << endl;
}
				/* evaluate the argument of vertex j, adj_vertex_1 and adj_vertex_2 about the boundary circle centre */					
				double theta0 = atan2(ordered_inner_hull_vertex[boundary_vertex_IH_index][1]-boundary_cog_y,ordered_inner_hull_vertex[boundary_vertex_IH_index][0]-boundary_cog_x);
				double theta1 = atan2(ordered_inner_hull_vertex[adj_vertex_1][1]-boundary_cog_y,ordered_inner_hull_vertex[adj_vertex_1][0]-boundary_cog_x);
				double theta2 = atan2(ordered_inner_hull_vertex[adj_vertex_2][1]-boundary_cog_y,ordered_inner_hull_vertex[adj_vertex_2][0]-boundary_cog_x);

				/* Move vertex j to the midpoint of the projections of the adjacent vertices on the boundary vertex circle.
				   If theta1 and theta2 have opposite signs, taking their average
				*/
				if (theta0 < 0)
					theta0 += 2*pi;
				if (theta1 < 0)
					theta1 += 2*pi;
				if (theta2 < 0)
					theta2 += 2*pi;
				double theta = (theta1+theta2)/2;

				if ((theta0 < theta1 && theta0 < theta2) ||  (theta0 > theta1 && theta0 > theta2))
					theta += pi;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:   theta0 = " << theta0 << ", theta1 = " << theta1 << ", theta2 = " << theta2 << ", theta = " << theta << ", boundary radius = " <<	boundary_vertex_circle_radius << endl;

				cog_coords[j][0] = boundary_vertex_circle_radius*cos(theta)+boundary_cog_x;
				cog_coords[j][1] = boundary_vertex_circle_radius*sin(theta)+boundary_cog_y;					

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement:   boundary vertex displaced to (" << cog_coords[j][0] << "," << cog_coords[j][1] << ")" << endl;
			}
		}
		
		/* if we're not considering the bundary vertices, they remain where they are */
		if (!INCLUDE_BOUNDARY_VERTICES)
		{
			for (int v=num_type123_vertices; v< nodecount; v++)
			{
				cog_coords[v][0] = vcoords[v][0];
				cog_coords[v][1] = vcoords[v][1];
			}				
		}
		
		/* now copy cog_coords to vcoords */
		 vcoords = cog_coords;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "centre_of_gravity_placement: updated vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "centre_of_gravity_placement:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
}
	
		if (TRACK_PLACEMENT_ITERATION && i % placement_iteration_tracking_step == placement_iteration_tracking_step-1)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: tracking force directed placement at iterative step i = " << i << endl;

			matrix<double> write_metapost_vcoords(2*nodecount,2);
			for (int j=0; j< nodecount; j++)
			{
				for (int k=0; k< 2; k++)
					write_metapost_vcoords[j][k] = vcoords[j][k];
			}
			for (int j=0; j< nodecount; j++)
			{
				for (int k=0; k< 2; k++)
					write_metapost_vcoords[j+nodecount][k] = initial_vcoords[j][k];
			}
			
			vector<double> vertex_radius(nodecount); // this will be initialized to zero, it will be unused by write_metapost.
			write_metapost(output, code_data, title, mp_control, write_metapost_vcoords, vertex_radius);	
		}
		
		initial_vcoords = vcoords;	
	}
	IHoutput.close();

	/* write the new location of the vertices to the output file */
	ofstream output(circlepack_output_file);
	
	if (!output)
	{
		cout << "\nError opening output file\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "centre_of_gravity_placement: could not open " << circlepack_output_file << endl;
		exit(0);
	}		
	
	output << nodecount << endl;
	for (int i=0; i< nodecount; i++)
	{
		output << setprecision(6) << setw(10) << vcoords[i][0] 		
		       << setprecision(6) << setw(10) << vcoords[i][1] << endl;		
	}

	/* write the radii as zero, since these are not used 
	   with force direction, since mp_control.circle_packing == false. 
	*/
	for (int i=0; i< nodecount; i++)
		output << setprecision(6) << setw(10) << 0 << endl;		
	
	output.close();
}

double repulsive_force(double distance)
{
	/* this is the value in the paper */
//	return k_constant*k_constant/distance;

	/* this is a test value */
	return k_constant*distance;
//return 0;

}

double attractive_force (double distance)
{
	return k_constant*distance*distance*distance;
}

