/**********************************************************************************
                            Region shrinking placement
                                   October 2018

This approach considers the area of the regions determined by turning cycles after 
initial vertex placement by circle packing and aims to achieve an aesthetic layout by 
adjusting the diagram so that the regions are all roughly of the same size.  It does
this by iteratively contracting a large area of the diagram, shrinking the 2-simplexes
in the sub-triangulation of the region towards their barycentre and reassembling the
shrunken triangulation by basing it at the midpoint of the region's edge "furthest 
away" from the congested area of the diagram.  

The algorithm terminates when all the compact regions of the diagram's complemente 
lie within region_shrinking_area_factor of their average area.  The amount by which
the 2-simplexes of "large" regions are shrunk is governed by region_shrinking_factor.

When shrinking a region, if that region includes a boundary edge in the infinite turning
cycle, we reassemble the triangulation from the type 2 vertex at that edge's midpoint.  
Note that there may be more than one such edge, as in the case of the Kishino knots, in which case 
the choice is the first such edge encountered as we traverse the region's boundary. If
there is no such edge but there is a type 1 vertex in the region boundary that lies in the
infinite turning cycle, we reassemble from that type 1 vertex.  Otherwise, we use a 
weighting of the region edges to determine the edge "furthest away" from 
the most congested part of the diagram, identified by compact regions of small area.  
We first assign a weight to each edge, which is the area of the adjacent region in 
which that edge also resides.  Then we assign to each edge an "opposite weight" by 
regarding the region as a polygon with n edges.  If n is even, the opposite weight is 
the weight of the edge directly opposite in the polygon and if n is odd it is the sum of 
the weights of the two edges directly opposite.  We consider the edge having the smallest 
opposite weight as being furthest away from the congested part of the diagram.

When reconstructing the triangulation, we reposition the vertices in the region's 
sub-triangulation, with the exception of any type 1 vertex lying in the infinite turning 
cycle.  This avoids squeezing the diagram by pulling in the boundary of the triangulated 
disc.

**********************************************************************************/
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
extern char const* triangulation_shrink_file;
extern float boundary_vertex_retraction_factor;
extern float region_shrinking_factor;
extern float region_shrinking_area_factor;
extern float average_triangulation_length_factor;
extern bool INCLUDE_BOUNDARY_VERTICES;
extern bool DRAW_IN_HYPERBOLIC_DISC;
extern bool USE_CENTRE_OF_GRAVITY_PLACEMENT;
extern bool RETRACT_BOUNDARY_VERTICES;
extern bool TRACK_PLACEMENT_ITERATION;
extern bool USE_KEN_STEPHENSON_CIRCLE_PACKING;
extern bool APPLY_FORCES_TO_TYPE_12_ONLY;
extern bool USE_SMALL_REGION_SHRINKING_PLACEMENT;

extern int metapost_coordinate_scale_factor;
extern int placement_iteration_tracking_step;

extern double EPSILON;
extern double ZERO_THRESHOLD;

/********************* Function prototypes ***********************/
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, matrix<double> vcoords, vector<double> vertex_radius, vector<int>* _state=0);
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, 
                    char const* circlepack_output_file, char const* circlepack_output_savefile, vector<int>* _state=0);
void set_frame_corners(string coordinate_file, metapost_control& mp_control);
int circle_pack (char const* inputfile, char const* outputfile);
int KS_circle_pack (char const* inputfile, char const* outputfile);
double area_of_triangle (double u_x, double u_y, double v_x, double v_y, double w_x, double w_y);
void evaluate_region_areas (vector<double>& region_area, matrix<double>& region_coords, matrix<int>& cycle, 
                            int num_cycles, int infinite_region, int num_crossings, vector<int>& type_3_vertex, matrix<double>& vcoords);
void get_region_coords (int region, matrix<double>& region_coords, matrix<int>& cycle, int num_cycles, int infinite_region, 
                        int num_crossings, vector<int>& type_3_vertex, matrix<double>& vcoords);
void shrink_triangulation (matrix<double>& shrunken_triangulation, int num_triangles, matrix<double>& region_coords);
void shrink_triangle (matrix<double>& shrunken_triangulation, int index, double u_x, double u_y, double v_x, double v_y, double w_x, double w_y);
int adjacent_region(matrix<int>& cycle, int num_cycles, int region, int edge, int& index);
double edge_weight(matrix<int>& cycle, int num_cycles, int region, int edge, vector<double>& region_area);


void region_shrinking_placement(metapost_control& mp_control, generic_code_data& code_data, matrix<int>& cycle, int num_cycles, 
                                int num_left_cycles, int infinite_region, int num_iterations, string title)
{

	int num_crossings = code_data.num_crossings;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	if (USE_SMALL_REGION_SHRINKING_PLACEMENT)
			debug << "region_shrinking_placement: using small region placement variant" << endl;
	debug << "region_shrinking_placement: code_data: ";
	write_code_data(debug,code_data);
	debug << endl;
	print_code_data(debug, code_data, "region_shrinking_placement: ");
	debug << "region_shrinking_placement: num_iterations = " << num_iterations << endl;
	if (TRACK_PLACEMENT_ITERATION)
		debug << "region_shrinking_placement: placement_iteration_tracking_step = " << placement_iteration_tracking_step << endl;
}

	/* record the frame corners for the region shrinking placement case
	mp_control.frame_minx = -metapost_coordinate_scale_factor;
	mp_control.frame_maxx = metapost_coordinate_scale_factor;
	mp_control.frame_miny = -metapost_coordinate_scale_factor;
	mp_control.frame_maxy = metapost_coordinate_scale_factor;
	*/

	/* turn off unsupported options */
	DRAW_IN_HYPERBOLIC_DISC = false;
	mp_control.circle_packing = false;
		
	if (USE_KEN_STEPHENSON_CIRCLE_PACKING)
		KS_circle_pack (triangulation_output_file, circlepack_output_file);
	else
		circle_pack (triangulation_output_file, circlepack_output_file);

	/* Read triangulation data from triangulation_output_file */
	int nodecount;
	int alpha;
	int beta;
	int gamma;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement: read triangulation data from " << triangulation_output_file << endl;

	ifstream triangulation(triangulation_output_file);
	if(!triangulation)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: could not open " << triangulation_output_file << endl;
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
	debug << "region_shrinking_placement: nodecount = " << nodecount << endl;
	}

	getline(triangulation,next_line);
	if (next_line.find("ALPHA") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> alpha; //always 1, the first vertex
		iss >> beta; //the first boundary (type 4) vertex
		iss >> gamma; // IS THIS USED?
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
	}
	
	int num_type123_vertices = beta-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement: number of type 1,2 & 3 vertices num_type123_vertices = " << num_type123_vertices << endl;

	/* Read the boundary orientation (i.e whether the infinite region was bounded by a left turning cycle or 
	   a right turning cycle) in order to place the boundary vertices correctly.
	   
	   Store the flowers for each node in a matrix with the first column set to the number of neighbours
	   of that vertex.  No vertex is adjacent to itself, nor is any vertex adjacent to all other vertices
	   so the maximum number of neighbours is <= nodecount -2.  Recall that for interior vertices the first
	   petal appears at the end of the list as well, to simplify the calculation of angle sums in circle packing.
	*/
//	bool clockwise_boundary = true;
	
	matrix<int> flowers(nodecount,nodecount);

	while (getline(triangulation,next_line))
	{
		
//		if (next_line.find("BOUNDARY") != string::npos)
//		{
//			if (next_line.find("anticlockwise") != string::npos)
//				clockwise_boundary = false;
//		}		
//		else if (next_line.find("FLOWERS") != string::npos)
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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "region_shrinking_placement: node " << node << ", petal count = " << count << ": " << first_petal;
				
				for (int i=0; i< count; i++)
				{
					iss >> petal;
					flowers[node][i+2] = petal;
					
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " " << petal;					
				}

				if (node < beta-1)
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "; interior vertex" << endl;					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "; boundary vertex" << endl;					
				}										
			}		
		}
	}
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "flowers" << endl;					
	print(flowers,debug,4,"");
}

	/* note the type 3 vertices corresponding to the barycentres
	   of regions having more than two edges.  Regions with two edges 
	   will be recorded as vertex zero but there is no ambiguity since
	   vertex zero always corresponds to crossing zero.
	*/
	vector<int> type_3_vertex(num_cycles);
	int temp = 3*num_crossings; // this is the next vertex number to allocate
	for (int i=0; i< num_cycles; i++)
	{
		if (cycle[i][0] > 2 && i != infinite_region )
			type_3_vertex[i] = temp++;
		else
			type_3_vertex[i] = 0;
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "region_shrinking_placement: type_3_vertex: ";
    for (int i=0; i< num_cycles; i++)
		debug << type_3_vertex[i] << ' ';
	debug << endl;
}
	
    matrix<double> vcoords(nodecount,2);

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
	
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "region_shrinking_placement: vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "region_shrinking_placement:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
}


	/* if, required, reduce the scale of the type 4 vertices */
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
	debug << "region_shrinking_placement: boundary vertex centre of gravity = (" << boundary_cog_x << "," << boundary_cog_y << ")" << endl;
	
		for (int i=num_type123_vertices; i< nodecount; i++)
		{
			vcoords[i][0] = (vcoords[i][0]-boundary_cog_x)*boundary_vertex_retraction_factor + boundary_cog_x;
			vcoords[i][1] = (vcoords[i][1]-boundary_cog_y)*boundary_vertex_retraction_factor + boundary_cog_y;
		}
	}

    matrix<double> initial_vcoords = vcoords;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "region_shrinking_placement: initial_vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "region_shrinking_placement:   v" << i << " = ("<< initial_vcoords[i][0] << ", " << initial_vcoords[i][1] << ")" << endl;
}
	
	/* Main iterative loop */	   
	bool algorithm_converged = false;
	
//	for (int i=0; i< num_iterations && !algorithm_converged; i++)
	for (int i=0; i< num_iterations; i++)
	{

//	    cout << "\nregion_shrinking_placement: iteration " << i << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: region shrinking iteration " << i << endl;
	
		/* Evaluate region areas and identify the region with whose triangulation we are going to shrink, recording its vertex 
		   coordinates in region_coords.  The region selected depends on the variant of the region shrinking algorithm used, 
		   in the original, default version it is simply the region with maximal area.  In the small version, it is the
		   region of maximal area that shares an edge with the region of smallest area.
		   
		   evaluate_region_areas records the region_coords so that shrink triangulation doesn't need to work around the maximal
		   cycle again.  The coordinates in region_coords are stored in the order of type 2 and type 1 vertices encountered by
		   tracing the turning cycle bounding the maximal region.		   
		*/
		vector<double> region_area(num_cycles);
		
		matrix<double> region_coords(0,0);
		evaluate_region_areas (region_area, region_coords, cycle, num_cycles, infinite_region, num_crossings, type_3_vertex, vcoords);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement: region areas " << endl;
	for (int i=0; i< num_cycles; i++)
		debug << "region_shrinking_placement:   region " << i << ": " << region_area[i] << endl;
}

		/* identify the regions or maximal and minimal area */
		int max_region = 0;
		int min_region = 0;

		for (int i=0; i< num_cycles; i++)
		{
			if (region_area[min_region] < 0 || (region_area[i] > 0 && region_area[i] < region_area[min_region]))
				min_region = i;

			if (region_area[i] > region_area[max_region])
				max_region = i;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: max_region = " << max_region << ", min_region = " << min_region << endl;

		/* if we're doing small region placement, identify the region of maximal size adjacent to (i.e. sharing
		   an edge with) min_region and replace max_region with that region, since that is the region we wish 
		   to shrink.  We record the shared edge between min_region and max_region in shared_edge and its (2nd) 
		   index in the maximal region's row of cycle in shared_edge_index
		*/
		int shared_edge = abs(cycle[min_region][1]);
		int shared_edge_index = 0;
		if (USE_SMALL_REGION_SHRINKING_PLACEMENT)
		{
			max_region = adjacent_region(cycle, num_cycles, min_region, abs(cycle[min_region][1]), shared_edge_index);			
			double max_adjacent_region_weight = edge_weight(cycle, num_cycles, min_region, abs(cycle[min_region][1]), region_area);
			
			for (int i = 2; i <= cycle[min_region][0]; i++)
			{
				double weight = edge_weight(cycle, num_cycles, min_region, abs(cycle[min_region][i]), region_area);
				
				if (weight > max_adjacent_region_weight)
				{
					shared_edge = abs(cycle[min_region][i]);
					max_adjacent_region_weight = weight;
					max_region = adjacent_region(cycle, num_cycles, min_region, abs(cycle[min_region][i]), shared_edge_index);
				}
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "region_shrinking_placement: small region shrinking: largest region adjacent to min_region " << min_region << " is " << max_region << endl;
	debug << "region_shrinking_placement: shared edge = " << shared_edge << ", index of shared_edge in boundary cycle of maximal region = " << shared_edge_index << endl;
}	
		}

        get_region_coords (max_region, region_coords, cycle, num_cycles, infinite_region, num_crossings, type_3_vertex, vcoords);


		int num_max_region_edges = cycle[max_region][0];	
		int num_triangles = (num_max_region_edges == 2? 2: 2*num_max_region_edges);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement: number of triangles in maximal region =  " << num_triangles << endl;
	debug << "region_shrinking_placement: region_coords:"<< endl;
	for (int i=0; i < (num_triangles == 2? 4 : num_triangles+1); i++)
		debug << "region_shrinking_placement:   (" << region_coords[i][0] << "," << region_coords[i][1] << ")" << endl;
}
	
		/* Check for convergence of the algorithm due to the regions all being sufficiently close to each other in terms
		   of their area.  In the default version, we check whether the maximal region's area differs significantly from 
		   the average area of the compact regions.  In the small variant, we compare the min_region.
		   
		   If the test region differs significantly from the average, we proceed to shrink the tringulation of the maximal 
		   region. If it does not, the algorithm has converged and we break out of the main loop.
		*/
		double average_region_area = 0;
		
		for (int j=0; j< num_cycles; j++)
		{
			if (region_area[j] > 0)
				average_region_area += region_area[j];
		}
		average_region_area /= (num_cycles-1); // not the infinite cycle
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: average area of compact regions = " << average_region_area << endl;

	
		if (    (USE_SMALL_REGION_SHRINKING_PLACEMENT && average_region_area/region_area[min_region] <= region_shrinking_area_factor)
		    || (!USE_SMALL_REGION_SHRINKING_PLACEMENT && region_area[max_region]/average_region_area <= region_shrinking_area_factor)
		   )
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: region shrinking algorithm converged after " << i << " iterations, all compact regions' area within " << region_shrinking_area_factor << " of the average." << endl;

			cout << "region shrinking algorithm converged after " << i << " iterations, all compact regions have area within " << region_shrinking_area_factor << " of the average" << endl;
	
			algorithm_converged = true;
			break;
		}

		/* Proceed by shrinking each of the triangles in the maximal region towards their barycentre by 
		   region_shrinking_factor. We have to store three vertices in the shrunken triangulation for each triangle 
		   in the maximal region. The triangles are ordered in pairs according to the order of type 1 vertices in 
		   the cycle bounding the maximal region, with the triangle including the type 1 vertex followed by the 
		   corresponding triangle incident with the barycentre.
		   
		   For each triangle, the shrunken vertices are stored in the order corresponding to: type 2 vertex, type 1 
		   vertex or barycentre, type 2 vertex.  So for each triangle, the first pair of coordinates (locations 0 and 1) 
		   provide the relative position in the shrunken triangulation of the type 1 vertex at the end of the 
		   corresponding edge and type 2 vertex within that edge, or the barycentre and that type 2 vertex.
		*/
		
		matrix<double> shrunken_triangulation(num_triangles,6);
		
		shrink_triangulation(shrunken_triangulation, num_triangles, region_coords);
		
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement: vertices in maximal region after shrinking:" << endl;
	for (int i=0; i< num_triangles; i++)
	{
		debug << "region_shrinking_placement: triangle " << i << ": (";
		debug << shrunken_triangulation[i][0] << "," << shrunken_triangulation[i][1] << "), (";
		debug << shrunken_triangulation[i][2] << "," << shrunken_triangulation[i][3] << "), (";
		debug << shrunken_triangulation[i][4] << "," << shrunken_triangulation[i][5] << ")" << endl;
	}
}

		/* write the shrunken triangulation vertex data to the triangulation_shrink_file so we can draw it if required */
		ofstream shrink(triangulation_shrink_file);
		
		if (!shrink)
		{
			cout << "\nError opening triangulation_shrink_file file\n";
	if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "region_shrinking_placement: could not open " << triangulation_shrink_file << endl;
			exit(0);
		}		
		
		shrink << num_triangles << endl;
		for (int j=0; j< num_triangles; j++)
		{
			shrink << setprecision(6) << setw(15) << shrunken_triangulation[j][0] << " "		
			       << setprecision(6) << setw(15) << shrunken_triangulation[j][1] << " "		
			       << setprecision(6) << setw(15) << shrunken_triangulation[j][2] << " "		
			       << setprecision(6) << setw(15) << shrunken_triangulation[j][3] << " "		
			       << setprecision(6) << setw(15) << shrunken_triangulation[j][4] << " "		
			       << setprecision(6) << setw(15) << shrunken_triangulation[j][5] << endl;		
		}
	
		shrink.close();

		/* For the default algorithm, when deciding from which vertex to reconstruct the triangulation, 
		   and for all algorithm variants during its reconstruction, we need to understand if the maximal 
		   region's boundary contains any type 1 or type 2 vertices lying in the infinite turning cycle, 
		   so we note these vertices in the infinite sycle first.
		*/
		vector<bool> infinite_type_12_vertex_flag(3*num_crossings); // initializes to false
		for (int j=1; j <= cycle[infinite_region][0]; j++)
		{
			int edge = cycle[infinite_region][j];				
			infinite_type_12_vertex_flag[abs(edge)+num_crossings] = true;
			int vertex = (edge < 0? (abs(edge)-1)/2: edge/2);
			infinite_type_12_vertex_flag[vertex] = true;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "region_shrinking_placement: type 1 & 2 vertices lying in the infinite turning cycle flags: ";
	for (int j=0; j< num_crossings; j++)
		debug << infinite_type_12_vertex_flag[j] << " ";
	debug << endl;
}

		/* If the maximal region includes an edge in the infinite turning cycle, we shall reconstruct the triangulation rooted at the
		   type 2 vertex in that edge.  If there is no such edge but there is a type 1 vertex in the region boundary that lies in the
		   infinite turning cycle, we reassemble from that type 1 vertex. Otherwise, we root the reconstruction at the type 2 vertex 
		   in the edge having minimum opposite weight.
		   
		   First, we determine which of these cases we are working with and as we do so, we calculate the individual edge weights
		   in the maximal region.
		*/
		vector<double> opposite_weight(num_max_region_edges);
		int minimum_opposite_weight_edge = 0;
		int min_opposite_edge_index = 0;
		int base_vertex = -1;  // used as a flag, initially
		double minimum_opposite_weight = 0;


if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement: calculating opposite weights" << endl;
		
		for (int j=0; j < num_max_region_edges; j++)
		{
			int edge = cycle[max_region][j+1];				

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "region_shrinking_placement:   maximal region edge " << edge;
			
			/* check if moving along this edge brings us to a vertex in the infinite turning cycle, if it does, 
			   it is a candidate base_vertex but we will overwrite this data if there is an edge of the maximal 
			   region's boundary in the infinite turning cycle, having evaluated the weight of that edge below. 
			   
			   Note that if the current edge is in the infinite turning cycle, then the following conditional 
			   clause and the one after (considering the weight) both set minimum_opposite_weight_edge and 
			   min_opposite_edge_index to the same values but set base_vertex differently, the second (type 2 
			   vertex) value being correct in this case.
			*/
			int vertex = (edge < 0? (abs(edge)-1)/2: edge/2);
			if (infinite_type_12_vertex_flag[vertex])
			{
				/* here we use minimum_opposite_weight_edge and min_opposite_edge_index to identify the
				   edge along which we arrive at the base vertex
				*/
				minimum_opposite_weight_edge = abs(edge);
				min_opposite_edge_index = j;
				base_vertex = vertex;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", arrives at infinite turning cycle vertex " << vertex << endl;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", does not arrive at infinite turning cycle vertex" << endl;
			}

			double weight = edge_weight(cycle, num_cycles, max_region, abs(edge), region_area);
			
			if (weight == -1)
			{
				/* this edge lies in the infinite turning cycle, there may be other edges in the boundary also in the
				   infinite turning cycle but we will use the type 2 vertex in this first one as the base vertex.
				*/
				
				minimum_opposite_weight_edge = abs(edge);
				min_opposite_edge_index = j;
				minimum_opposite_weight = -1; // indicates edge in infinite turning cycle

				/* set the base vertex, at which the reconstructed triangulation will 
				   be rooted to the type 2 vertex corresponding to this edge
				*/
				base_vertex = minimum_opposite_weight_edge + num_crossings;

				break;
			}
			
			int opposite_edge = abs(cycle[max_region][(j + num_max_region_edges/2) % num_max_region_edges +1]);
			opposite_weight[j] = edge_weight(cycle, num_cycles, max_region, opposite_edge, region_area); 
			
			if (num_max_region_edges % 2)
			{
				/* we need to check whether the second edge bounds the infinite region, so we can set opposite_weight[j] correctly. */
				opposite_edge = abs(cycle[max_region][(j + num_max_region_edges/2 + 1) % num_max_region_edges +1]);
				
				weight = edge_weight(cycle, num_cycles, max_region, opposite_edge, region_area);
				if (weight == -1)
					opposite_weight[j] = -1; // i.e. infinite
				else
					opposite_weight[j] += weight;
			}
		}

		/* If minimum_opposite_weight has not been set, we have not encountered an edge of the maximal region in the 
		   infinite turning cycle.  In that case, if base_vertex has been set, we are reconstructing from 
		   that type 1 vertex.  Otherwise we need to consider the opposite weights we've just calculated.
		*/
		if (!minimum_opposite_weight)
		{
			if (base_vertex != -1)
			{
				/* We have a type 1 vertex in the maximal region boundary that lies in the infinite 
				   turning cycle, so reconstruct from this point.
				*/
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement: using the type 1 vertex " << base_vertex << " in the maximal region boundary that lies in the infinite turning cycle as the base vertex." << endl;
	debug << "region_shrinking_placement: arrived at the base vertex on edge " << minimum_opposite_weight_edge << endl;
	debug << "region_shrinking_placement: index of " << minimum_opposite_weight_edge << " in turning cycle bounding maximal region = " << min_opposite_edge_index << endl;
}
			}
			else
			{
				/* the boundary of the maximal region does not intersect the infinite turning cycle, so 
				   reconstruct the triangulation at the type 2 vertex in the edge of the maximal region 
				   having minimum opposite weight.  In the case of the small region shrinking variant, we 
				   require that the edge containing the base vertex is opposite the shared edge, meaning
				   there is either a unique choice (if the maximal region has an even number of edges), or
				   is one of only two possibilities (when it has an odd number of edges).
				   
				   In the default variant, we must consider all edges in the maximal region.
				*/
				
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement: opposite_weight: ";
	for (int k=0; k< num_max_region_edges; k++)
		debug << opposite_weight[k] << " ";
	debug << endl;
}

				if (USE_SMALL_REGION_SHRINKING_PLACEMENT)
				{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "region_shrinking_placement: selecting the base_vertex using small region shrinking procedures" << endl;
	
					if (num_max_region_edges % 2)
					{
						/* we consider the two edges opposite the shared edge in the maximal region and select the one 
						   with minimal opposite weight
						*/
						int first_opposite_edge_index = (shared_edge_index + num_max_region_edges/2) % num_max_region_edges;
						int second_opposite_edge_index = (shared_edge_index + num_max_region_edges/2 +1) % num_max_region_edges;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "region_shrinking_placement:   maximal region has an odd number of edges, first_opposite_edge_index = " << first_opposite_edge_index
	      << ", second_opposite_edge_index = " << second_opposite_edge_index << endl;
	debug << "region_shrinking_placement:   first index edge opposite weight = " << opposite_weight[first_opposite_edge_index]
	      << ", second index edge opposite weight = " << opposite_weight[second_opposite_edge_index] << endl;
}
						
						if (opposite_weight[first_opposite_edge_index] < opposite_weight[second_opposite_edge_index])
							min_opposite_edge_index = first_opposite_edge_index;
						else
							min_opposite_edge_index = second_opposite_edge_index;

						minimum_opposite_weight_edge = abs(cycle[max_region][min_opposite_edge_index+1]);
						minimum_opposite_weight = opposite_weight[min_opposite_edge_index];
					}
					else
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "region_shrinking_placement:   maximal region has an even number of edges" << endl;

						/* we select the edge opposite the shared_edge in the maximal region*/
						min_opposite_edge_index = (shared_edge_index + num_max_region_edges/2) % num_max_region_edges;
						minimum_opposite_weight_edge = abs(cycle[max_region][min_opposite_edge_index+1]);
						minimum_opposite_weight = opposite_weight[min_opposite_edge_index];
					}
				}
				else
				{
					minimum_opposite_weight_edge = abs(cycle[max_region][1]);
					min_opposite_edge_index = 0; // used to adjust the max region barycentre
					minimum_opposite_weight = opposite_weight[0];
					
					for (int j=1; j < num_max_region_edges; j++)
					{
						/* we might have initialized minimum_opposite_weight based on the infinite region weight of -1, so for j=1 always want to replace this value */
						if (minimum_opposite_weight < 0 || ( opposite_weight[j] > 0 && opposite_weight[j] < minimum_opposite_weight))
						{
							minimum_opposite_weight_edge = abs(cycle[max_region][j+1]);
							minimum_opposite_weight = opposite_weight[j];
							min_opposite_edge_index = j;
						}
					}
				}
				
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement: minimum opposite weight = " << minimum_opposite_weight << " corresponding to edge " << minimum_opposite_weight_edge << endl;
	debug << "region_shrinking_placement: minimum opposite weight edge index in turning cycle bounding maximal region = " << min_opposite_edge_index << endl;
}

				/* set the base vertex, at which the reconstructed triangulation will be rooted */
				base_vertex = minimum_opposite_weight_edge + num_crossings;

			}
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: base vertex for reconstructing triangulation = " << base_vertex << endl;

		/* To reposition the triangulation of the maximal region at the base vertex we have to calculate the new location
		   of each vertex, which we do by considering the relative position of the corresponding vertices in the shrunken
		   triangulation.  We use the base vertex as the starting reference.  We have to handle the case that the base 
		   vertex is a type 1 vertex separately.
		*/
	
		if (base_vertex < num_crossings)
		{
			/* There is no edge in the boundary of the maximal region that lies in the infinite turning cycle but there
			   is a vertex, our base vertex.
			   
			   If the maximal region has more than two edges then each of the type 2 vertices in the region boundary are joined 
			   by an edge to the barycentre.  We therefore start by adjusting the type 2 vertex on the edge along which we arrive 
			   at the base_vertex when traversing the boundary of the maximal region and then adjust the barycentre using the edge
			   connecting to this type 2 vertex.  Then, we work around the turning cycle bounding the maximal region.  For each edge, 
			   other than the one taking us to the base vertex, we adjust the location of the type 2 vertex it contains using the 
			   barycentre.  We then adjust the next type 1 vertex using the type 2 vertex in the previous edge.
			    
			   If the maximal region has just two edges, then it has just two triangles one of which has base_vertex in its 
			   boundary.  We adjust the two type 2 vertices using this triangle, then adjust the other type 1 vertex
			*/
			if (num_max_region_edges > 2)
			{
			
				/* move the type 2 vertex on the edge taking us to the base vertex towards the base_vertex, the index of this 
				   edge in the turning cycle determines which two triangles in shrunken_triangulation correspond to this edge 
				   and we need the first of these two triangles to identify the appropriate shrunken_triangulation coordinates
				*/
				int edge_u = minimum_opposite_weight_edge; // in this case the edge taking us to the base vertex
				int vertex_u = edge_u+num_crossings;
						
				double u_x = shrunken_triangulation[2*min_opposite_edge_index][0];
				double u_y = shrunken_triangulation[2*min_opposite_edge_index][1];
				double b_x = shrunken_triangulation[2*min_opposite_edge_index][2];
				double b_y = shrunken_triangulation[2*min_opposite_edge_index][3];
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 2 vertex " << vertex_u << " to type 1 vertex " << base_vertex << endl;
				
				double shift_x = u_x - b_x;
				double shift_y = u_y - b_y;
	
				vcoords[vertex_u][0] = vcoords[base_vertex][0] + shift_x;
				vcoords[vertex_u][1] = vcoords[base_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_2_vertex (" << u_x << "," << u_y << "), base_vertex (" << b_x << "," << b_y 
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_2 vertex coordinates = (" << vcoords[vertex_u][0] << "," << vcoords[vertex_u][1] << ")" << endl;
}									
				
				/* adjust the barycentre of the maximal region to vertex_u.  We use the second triangle identified by
				   min_opposite_edge_index associated with this edge to adjust the barycentre's position
				*/
				int barycentre_vertex = type_3_vertex[max_region];

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust barycentre vertex " << barycentre_vertex << " to " << vertex_u << endl;

				double barycentre_x = shrunken_triangulation[2*min_opposite_edge_index+1][2];
				double barycentre_y = shrunken_triangulation[2*min_opposite_edge_index+1][3];
				shift_x = barycentre_x - u_x;
				shift_y = barycentre_y - u_y;
				
				vcoords[barycentre_vertex][0] = vcoords[vertex_u][0] + shift_x;
				vcoords[barycentre_vertex][1] = vcoords[vertex_u][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     shrunken barycentre ( " << barycentre_x << "," << barycentre_y << "), vertex_u (" << u_x << "," << u_y
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new barycentre coordinates = (" << vcoords[barycentre_vertex][0] << "," << vcoords[barycentre_vertex][1] << ")" << endl;
}				
			
				for (int j=0; j < num_max_region_edges; j++)
				{
					/* adjust the type 2 vertex of the edge, if required */
					int edge = cycle[max_region][j+1];				
					int type_2_vertex = abs(edge)+num_crossings;
					double type_2_x;
					double type_2_y;
					
					if (abs(edge) == minimum_opposite_weight_edge)
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   edge " << abs(edge) << " takes us to the base vertex, already adjusted" << endl;
	
						continue;
					}
					else
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust vertex " << type_2_vertex << " to barycentre vertex " << barycentre_vertex << endl;
	
						type_2_x = shrunken_triangulation[2*j+1][0];
						type_2_y = shrunken_triangulation[2*j+1][1];
						barycentre_x = shrunken_triangulation[2*j+1][2];
						barycentre_y = shrunken_triangulation[2*j+1][3];
						shift_x = type_2_x - barycentre_x;
						shift_y = type_2_y - barycentre_y;
						
						vcoords[type_2_vertex][0] = vcoords[barycentre_vertex][0] + shift_x;
						vcoords[type_2_vertex][1] = vcoords[barycentre_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_2_vertex (" << type_2_x << "," << type_2_y << "), shrunken barycentre ( " 
	      << barycentre_x << "," << barycentre_y << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_2 vertex coordinates = (" << vcoords[type_2_vertex][0] << "," << vcoords[type_2_vertex][1] << ")" << endl; 
}			

					}
					
					/* adjust the type 1 vertex at the end of the edge, this vertex might still be in the infinite turning cycle, 
					   even though there is no edge of the maximal regions' boundary in the infinite turning cycle.
					*/
					int type_1_vertex = (edge < 0? (abs(edge)-1)/2: edge/2);

					if (!infinite_type_12_vertex_flag[type_1_vertex])
					{
						
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 1 vertex " << type_1_vertex << " to type 2 vertex " << type_2_vertex << endl;
	
						type_2_x = shrunken_triangulation[2*j][0];
						type_2_y = shrunken_triangulation[2*j][1];
						double type_1_x = shrunken_triangulation[2*j][2];
						double type_1_y = shrunken_triangulation[2*j][3];
						shift_x = type_1_x - type_2_x;
						shift_y = type_1_y - type_2_y;
					
						vcoords[type_1_vertex][0] = vcoords[type_2_vertex][0] + shift_x;
						vcoords[type_1_vertex][1] = vcoords[type_2_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_1_vertex (" << type_1_x << "," << type_1_y << "), type_2_vertex (" << type_2_x << "," << type_2_y
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_1 vertex coordinates = (" << vcoords[type_1_vertex][0] << "," << vcoords[type_1_vertex][1] << ")" << endl; 
}			
					}
					else
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   type 1 vertex " << type_1_vertex << " lies in the infinite turning cycle, no adjustment needed" << endl;
					}
				}
				
			}
			else
			{
				/* move the two type 2 vertices towards the base_vertex */
				int edge_u = cycle[max_region][1];				
				int vertex_u = abs(edge_u)+num_crossings;
				int edge_w = cycle[max_region][2];				
				int vertex_w = abs(edge_w)+num_crossings;
				
				/* identify which of the two triangles in the region contains the base vertex */
				int triangle = (min_opposite_edge_index == 0? 0: 1);
				double b_x = shrunken_triangulation[triangle][2];
				double b_y = shrunken_triangulation[triangle][3];
		
				/* we can get the coordinates of the shurunken type 2 vertices from either triangle */
				double u_x = shrunken_triangulation[0][0];
				double u_y = shrunken_triangulation[0][1];
				double w_x = shrunken_triangulation[0][4];
				double w_y = shrunken_triangulation[0][5];

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 2 vertex " << vertex_u << " to type 2 vertex " << base_vertex << endl;
				
				double shift_x = u_x - b_x;
				double shift_y = u_y - b_y;
	
				vcoords[vertex_u][0] = vcoords[base_vertex][0] + shift_x;
				vcoords[vertex_u][1] = vcoords[base_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_2_vertex (" << u_x << "," << u_y << "), base_vertex (" << b_x << "," << b_y 
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_2 vertex coordinates = (" << vcoords[vertex_u][0] << "," << vcoords[vertex_u][1] << ")" << endl;
}			

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 2 vertex " << vertex_w << " to type 2 vertex " << base_vertex << endl;
				
				shift_x = w_x - b_x;
				shift_y = w_y - b_y;
	
				vcoords[vertex_w][0] = vcoords[base_vertex][0] + shift_x;
				vcoords[vertex_w][1] = vcoords[base_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_2_vertex (" << w_x << "," << w_y << "), base_vertex (" << b_x << "," << b_y 
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_2 vertex coordinates = (" << vcoords[vertex_w][0] << "," << vcoords[vertex_w][1] << ")" << endl;
}			

				/* move the other type 1 vertices towards the type 2 vertex vertex_u (we could use either vertex_u or vertex _v as a 
				   reference for the other type 1 vertex, provided that the other type 1 vertex is not in the infinite turning cycle 
				*/
				int edge = cycle[max_region][(min_opposite_edge_index == 0? 1: 0)];				
				int type_1_vertex = (edge < 0? (abs(edge)-1)/2: edge/2);

				if (!infinite_type_12_vertex_flag[type_1_vertex])
				{

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 1 vertex " << type_1_vertex << " to type 2 vertex " << vertex_u << endl;	

	
					double type_1_x = shrunken_triangulation[(triangle == 0? 1:0)][2];
					double type_1_y = shrunken_triangulation[(triangle == 0? 1:0)][3];
					shift_x = type_1_x - u_x;
					shift_y = type_1_y - u_y;
					
					vcoords[type_1_vertex][0] = vcoords[base_vertex][0] + shift_x;
					vcoords[type_1_vertex][1] = vcoords[base_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_1_vertex (" << type_1_x << "," << type_1_y << "), vertex_u (" << u_x << "," << u_y
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_1 vertex coordinates = (" << vcoords[type_1_vertex][0] << "," << vcoords[type_1_vertex][1] << ")" << endl; 
}			
				}
				else
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   type 1 vertex " << type_1_vertex << " lies in the infinite turning cycle, no adjustment needed" << endl;
				}
			}
		}
		else
		{
					
			/* If the maximal region has more than two edges then each of the type 2 vertices in the region boundary are joined 
			   by an edge to the barycentre.  We therefore start by adjusting the barycentre using the edge connecting the base 
			   vertex and the barycentre, then work around the turning cycle bounding the maximal region.  For each edge, if
			   the type 2 vertex is not the base vertex, and does not lie in the infinite turning cycle, we adjust its location 
			   using the barycentre.  We then adjust the next type 1 vertex using the type 2 vertex in the previous edge.
			    
			   If the maximal region has just two edges, then it has just two triangles and so both have base_vertex in their 
			   boundary.
			*/
			if (num_max_region_edges > 2)
			{
				/* adjust the barycentre of the maximal region: min_opposite_edge_index tells us which rows of shrunken_triangulation
				   correspond to the edge of minimum opposite weight.  We use the second triangle associated with this edge to adjust
				   the barycentre's position
				*/
				int barycentre_vertex = type_3_vertex[max_region];

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust barycentre vertex " << barycentre_vertex << " to " << base_vertex << endl;

				double base_x = shrunken_triangulation[2*min_opposite_edge_index+1][0];
				double base_y = shrunken_triangulation[2*min_opposite_edge_index+1][1];
				double barycentre_x = shrunken_triangulation[2*min_opposite_edge_index+1][2];
				double barycentre_y = shrunken_triangulation[2*min_opposite_edge_index+1][3];
				double shift_x = barycentre_x - base_x;
				double shift_y = barycentre_y - base_y;
				
				vcoords[barycentre_vertex][0] = vcoords[base_vertex][0] + shift_x;
				vcoords[barycentre_vertex][1] = vcoords[base_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     shrunken barycentre ( " << barycentre_x << "," << barycentre_y << "), base_vertex (" << base_x << "," << base_y
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new barycentre coordinates = (" << vcoords[barycentre_vertex][0] << "," << vcoords[barycentre_vertex][1] << ")" << endl;
}			
	
			
				for (int j=0; j < num_max_region_edges; j++)
				{
					/* adjust the type 2 vertex of the edge, if required */
					int edge = cycle[max_region][j+1];				
					int type_2_vertex = abs(edge)+num_crossings;
					double type_2_x;
					double type_2_y;
					
					if (type_2_vertex == base_vertex)
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   edge " << abs(edge) << " contains the base vertex, no adjustment needed" << endl;
					}
					else if (!infinite_type_12_vertex_flag[type_2_vertex])
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust vertex " << type_2_vertex << " to barycentre vertex " << barycentre_vertex << endl;
	
						type_2_x = shrunken_triangulation[2*j+1][0];
						type_2_y = shrunken_triangulation[2*j+1][1];
						barycentre_x = shrunken_triangulation[2*j+1][2];
						barycentre_y = shrunken_triangulation[2*j+1][3];
						shift_x = type_2_x - barycentre_x;
						shift_y = type_2_y - barycentre_y;
						
						vcoords[type_2_vertex][0] = vcoords[barycentre_vertex][0] + shift_x;
						vcoords[type_2_vertex][1] = vcoords[barycentre_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_2_vertex (" << type_2_x << "," << type_2_y << "), shrunken barycentre ( " 
	      << barycentre_x << "," << barycentre_y << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_2 vertex coordinates = (" << vcoords[type_2_vertex][0] << "," << vcoords[type_2_vertex][1] << ")" << endl; 
}			

					}
					else
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   type 2 vertex " << type_2_vertex << " lies in the infinite turning cycle, no adjustment needed" << endl;
					}
					
					/* adjust the type 1 vertex at the end of the edge */
					int type_1_vertex = (edge < 0? (abs(edge)-1)/2: edge/2);
					
					if (!infinite_type_12_vertex_flag[type_1_vertex])
					{

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 1 vertex " << type_1_vertex << " to type 2 vertex " << type_2_vertex << endl;
	
						type_2_x = shrunken_triangulation[2*j][0];
						type_2_y = shrunken_triangulation[2*j][1];
						double type_1_x = shrunken_triangulation[2*j][2];
						double type_1_y = shrunken_triangulation[2*j][3];
						shift_x = type_1_x - type_2_x;
						shift_y = type_1_y - type_2_y;
					
						vcoords[type_1_vertex][0] = vcoords[type_2_vertex][0] + shift_x;
						vcoords[type_1_vertex][1] = vcoords[type_2_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_1_vertex (" << type_1_x << "," << type_1_y << "), type_2_vertex (" << type_2_x << "," << type_2_y
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_1 vertex coordinates = (" << vcoords[type_1_vertex][0] << "," << vcoords[type_1_vertex][1] << ")" << endl; 
}			
					}
					else
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   type 1 vertex " << type_1_vertex << " lies in the infinite turning cycle, no adjustment needed" << endl;
					}
				}
			}
			else
			{
				/* move the type 2 vertex that is not the base vertex towards the base_vertex */
				int edge_u = cycle[max_region][1];				
				int vertex_u = abs(edge_u)+num_crossings;
				int edge_w = cycle[max_region][2];				
				int vertex_w = abs(edge_w)+num_crossings;

				if ((min_opposite_edge_index == 0 && !infinite_type_12_vertex_flag[vertex_w]) || (min_opposite_edge_index == 1 && !infinite_type_12_vertex_flag[vertex_u])) 
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 2 vertex " << (min_opposite_edge_index == 0?vertex_w: vertex_u) << " to type 2 vertex " << base_vertex << endl;
			
					double u_x = shrunken_triangulation[0][0];
					double u_y = shrunken_triangulation[0][1];
					double w_x = shrunken_triangulation[0][4];
					double w_y = shrunken_triangulation[0][5];
					
					double shift_x = (min_opposite_edge_index == 0? w_x - u_x: u_x - w_x);
					double shift_y = (min_opposite_edge_index == 0? w_y - u_y: u_y - w_y);
		
					vcoords[(min_opposite_edge_index == 0?vertex_w: vertex_u)][0] = vcoords[base_vertex][0] + shift_x;
					vcoords[(min_opposite_edge_index == 0?vertex_w: vertex_u)][1] = vcoords[base_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_2_vertex (" << (min_opposite_edge_index == 0? u_x: w_x) << "," << (min_opposite_edge_index == 0? u_y: w_y) 
	      << "), base_vertex (" << (min_opposite_edge_index == 0? w_x: u_x) << "," << (min_opposite_edge_index == 0? w_y: u_y) 
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_2 vertex coordinates = (" << vcoords[(min_opposite_edge_index == 0?vertex_w: vertex_u)][0] << "," 
	      << vcoords[(min_opposite_edge_index == 0?vertex_w: vertex_u)][1] << ")" << endl;
}			
				}
				else
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   type 2 vertex " << (min_opposite_edge_index == 0?vertex_w: vertex_u)  << " lies in the infinite turning cycle, no adjustment needed" << endl;
				}


				/* move the type 1 vertices at the end of each edge towards the base vertex */
				for (int j=0; j<2; j++)
				{ 
					int edge = cycle[max_region][j+1];				
					int type_1_vertex = (edge < 0? (abs(edge)-1)/2: edge/2);
	
					if (!infinite_type_12_vertex_flag[type_1_vertex])
					{

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   adjust type 1 vertex " << type_1_vertex << " to type 2 vertex " << base_vertex << endl;
	
						double type_1_x = shrunken_triangulation[j][2];
						double type_1_y = shrunken_triangulation[j][3];
						double shift_x = type_1_x - shrunken_triangulation[0][(min_opposite_edge_index == 0? 0: 4)];
						double shift_y = type_1_y - shrunken_triangulation[0][(min_opposite_edge_index == 0? 1: 5)];
						
						vcoords[type_1_vertex][0] = vcoords[base_vertex][0] + shift_x;
						vcoords[type_1_vertex][1] = vcoords[base_vertex][1] + shift_y;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "region_shrinking_placement:     type_1_vertex (" << type_1_x << "," << type_1_y << "), base_vertex (" 
	      << shrunken_triangulation[0][(min_opposite_edge_index == 0? 0: 4)] << "," << shrunken_triangulation[0][(min_opposite_edge_index == 0? 1: 5)]
	      << "), shift = (" << shift_x << "," << shift_y << ")" << endl;
	debug << "region_shrinking_placement:     new type_1 vertex coordinates = (" << vcoords[type_1_vertex][0] << "," << vcoords[type_1_vertex][1] << ")" << endl; 
}			
					}
					else
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "region_shrinking_placement:   type 1 vertex " << type_1_vertex << " lies in the infinite turning cycle, no adjustment needed" << endl;
					}
				}
			}
		}
	
		if (TRACK_PLACEMENT_ITERATION && i % placement_iteration_tracking_step == placement_iteration_tracking_step-1)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: tracking region shrinking placement at iterative step i = " << i << endl;
	
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

	/* write the new location of the vertices to the output file */
	ofstream output(circlepack_output_file);
	
	if (!output)
	{
		cout << "\nError opening output file\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "region_shrinking_placement: could not open " << circlepack_output_file << endl;
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

	if (!algorithm_converged)
		cout << "region shrinking algorithm terminated after " << num_iterations << " iterations" << endl;
	
	output.close();
}

/* area_of_triangle uses Heron's rule for calculating the area of a triangle from the length of its sides,
   determined from the coordinates of vertices u,v and w
*/
double area_of_triangle (double u_x, double u_y, double v_x, double v_y, double w_x, double w_y)
{
	double uv = sqrt ((v_x - u_x)*(v_x - u_x) + (v_y - u_y)*(v_y - u_y));
	double vw = sqrt ((w_x - v_x)*(w_x - v_x) + (w_y - v_y)*(w_y - v_y));
	double wu = sqrt ((u_x - w_x)*(u_x - w_x) + (u_y - w_y)*(u_y - w_y));
	
	double s = (uv + vw + wu) /2;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "area_of_triangle: length of uv = " << uv << ", vw = " <<  vw << ", wu = " << wu << ", s = " << s << endl;

	double area =  sqrt(s*(s-uv)*(s-vw)*(s-wu));

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "area_of_triangle: area = " << area << endl;
	
    return area;
}

/* evaluate_region_areas records the total area of each region's sub-trianglation in region_area. */

void evaluate_region_areas (vector<double>& region_area, matrix<double>& region_coords, matrix<int>& cycle, 
                            int num_cycles, int infinite_region, int num_crossings, vector<int>& type_3_vertex, matrix<double>& vcoords)
{
	for (int i=0; i< num_cycles; i++)
	{
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "evaluate_region_areas: cycle " << i << " of length " <<  cycle[i][0] << endl;

		if (i == infinite_region)
		{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "evaluate_region_areas: cycle " << i << " bounds infinite region, setting area to -1" << endl;
			region_area[i] = -1; // infinite
		}
		else if ( cycle[i][0] == 2)
		{
			matrix<double> coords(4,2);
			
			int edge_u = cycle[i][1];	
			int vertex_u = abs(edge_u)+num_crossings;
			coords[0][0] = vcoords[vertex_u][0];
			coords[0][1] = vcoords[vertex_u][1];
	

			int vertex_v = (edge_u < 0? (abs(edge_u)-1)/2: edge_u/2);
			coords[1][0] = vcoords[vertex_v][0];
			coords[1][1] = vcoords[vertex_v][1];

			int edge_w = cycle[i][2];		
			int vertex_w = abs(edge_w)+num_crossings;
			coords[2][0] = vcoords[vertex_w][0];
			coords[2][1] = vcoords[vertex_w][1];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	
	debug << "evaluate_region_areas:   first triangle" << endl;
	debug << "evaluate_region_areas:     edge_u = " << edge_u << " vertex_u = " <<  vertex_u << " u_coords = (" << coords[0][0] << "," << coords[0][1] << ")" << endl;
	debug << "evaluate_region_areas:     vertex_v = " << vertex_v << " v_coords = (" << coords[1][0] << "," << coords[1][1] << ")" << endl;
	debug << "evaluate_region_areas:     edge_w = " << edge_w << " vertex_w = " <<  vertex_w << " w_coords = (" << coords[2][0] << "," << coords[2][1] << ")" << endl;
}
			
			region_area[i] = area_of_triangle(coords[0][0],coords[0][1],coords[1][0],coords[1][1],coords[2][0],coords[2][1]);
	
			vertex_v = (edge_w < 0? (abs(edge_w)-1)/2: edge_w/2);
			coords[3][0] = vcoords[vertex_v][0];
			coords[3][1] = vcoords[vertex_v][1];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	
	debug << "evaluate_region_areas:   second triangle" << endl;
	debug << "evaluate_region_areas:     edge_u = " << edge_u << " vertex_u = " <<  vertex_u << " u_coords = (" << coords[0][0] << "," << coords[0][1] << ")" << endl;
	debug << "evaluate_region_areas:     vertex_v = " << vertex_v << " v_coords = (" << coords[3][0] << "," << coords[3][1] << ")" << endl;
	debug << "evaluate_region_areas:     edge_w = " << edge_w << " vertex_w = " <<  vertex_w << " w_coords = (" << coords[2][0] << "," << coords[2][1] << ")" << endl;
}
			
			region_area[i] += area_of_triangle(coords[0][0],coords[0][1],coords[3][0],coords[3][1],coords[2][0],coords[2][1]);			
		}
		else
		{
			int n = 2*cycle[i][0]+1; // number of vertices in region
			matrix<double> coords(n,2);

			int barycentre_vertex = type_3_vertex[i];
			coords[n-1][0] = vcoords[barycentre_vertex][0];
			coords[n-1][1] = vcoords[barycentre_vertex][1];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "evaluate_region_areas:   barycentre_vertex = " << barycentre_vertex << " barycentre coords = (" << coords[n-1][0] << "," << coords[n-1][1] << ")" << endl;
			
			int edge_u = cycle[i][1];	
			int vertex_u = abs(edge_u)+num_crossings;
			
			double u_x = coords[0][0] = vcoords[vertex_u][0];
			double u_y = coords[0][1] = vcoords[vertex_u][1];

			region_area[i] = 0;
			double v_x,v_y,w_x,w_y;
			
			for (int j=1; j <= cycle[i][0]; j++)
			{
				/* identify the coordinates of the vertex after edge_u */
				int vertex_v = (edge_u < 0? (abs(edge_u)-1)/2: edge_u/2);
				v_x = coords[2*j-1][0] = vcoords[vertex_v][0];
				v_y = coords[2*j-1][1] = vcoords[vertex_v][1];
				

				/* identify the coordinates of the midpoint of the edge after edge_u in the cycle */
				int edge_w;
				int vertex_w;
				if (j+1 <= cycle[i][0])
				{
					edge_w = cycle[i][j+1];
					vertex_w = abs(edge_w)+num_crossings;
					w_x = coords[2*j][0] = vcoords[vertex_w][0];
					w_y = coords[2*j][1] = vcoords[vertex_w][1];
				}
				else
				{
					edge_w = cycle[i][1];
					vertex_w = abs(edge_w)+num_crossings;
					w_x = coords[0][0];
					w_y = coords[0][1];
				}
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	
	debug << "evaluate_region_areas:   j = " << j << endl;
	debug << "evaluate_region_areas:     edge_u = " << edge_u << " vertex_u = " <<  vertex_u << " u_coords = (" << u_x << "," << u_y << ")" << endl;
	debug << "evaluate_region_areas:     vertex_v = " << vertex_v << " v_coords = (" << v_x << "," << v_y << ")" << endl;
	debug << "evaluate_region_areas:     edge_w = " << edge_w << " vertex_w = " <<  vertex_w << " w_coords = (" << w_x << "," << w_y << ")" << endl;
}

				/* add to the region area the area of the two triangles either side of the edge uw */
				region_area[i] += area_of_triangle(u_x,u_y,v_x,v_y,w_x,w_y);
				region_area[i] += area_of_triangle(u_x,u_y,coords[n-1][0],coords[n-1][1],w_x,w_y);
				
				
				/* move around the cycle */
				edge_u = edge_w;
				vertex_u = vertex_w;  // only needed for debugging
				u_x = w_x;
				u_y = w_y;
			}		
		}
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "area of region" << i << " = " << region_area[i] << endl;
		
	}
}

/* get_region_coords records the coordinates of the vertices in the given region in region_vcoords.
   
   The coordinates are stored in the order determined by the corresponding turning cycle, followed by the type 3 vertex 
   coordinates, if one exists within the region.
   
*/

void get_region_coords (int region, matrix<double>& region_coords, matrix<int>& cycle, int num_cycles, int infinite_region, 
                        int num_crossings, vector<int>& type_3_vertex, matrix<double>& vcoords)
{
		
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "get_region_coords: getting coordinates of region " << region << " of length " <<  cycle[region][0] << endl;

	if ( cycle[region][0] == 2)
	{
		matrix<double> coords(4,2);
		
		int edge_u = cycle[region][1];	
		int vertex_u = abs(edge_u)+num_crossings;
		coords[0][0] = vcoords[vertex_u][0];
		coords[0][1] = vcoords[vertex_u][1];


		int vertex_v = (edge_u < 0? (abs(edge_u)-1)/2: edge_u/2);
		coords[1][0] = vcoords[vertex_v][0];
		coords[1][1] = vcoords[vertex_v][1];

		int edge_w = cycle[region][2];		
		int vertex_w = abs(edge_w)+num_crossings;
		coords[2][0] = vcoords[vertex_w][0];
		coords[2][1] = vcoords[vertex_w][1];
	
		vertex_v = (edge_w < 0? (abs(edge_w)-1)/2: edge_w/2);
		coords[3][0] = vcoords[vertex_v][0];
		coords[3][1] = vcoords[vertex_v][1];
	
		region_coords = coords;
	}
	else
	{
		int n = 2*cycle[region][0]+1; // number of vertices in region
		matrix<double> coords(n,2);

		int barycentre_vertex = type_3_vertex[region];
		coords[n-1][0] = vcoords[barycentre_vertex][0];
		coords[n-1][1] = vcoords[barycentre_vertex][1];
	
		int edge_u = cycle[region][1];	
		int vertex_u = abs(edge_u)+num_crossings;
		
		coords[0][0] = vcoords[vertex_u][0];
		coords[0][1] = vcoords[vertex_u][1];

		for (int j=1; j <= cycle[region][0]; j++)
		{
			/* identify the coordinates of the vertex after edge_u */
			int vertex_v = (edge_u < 0? (abs(edge_u)-1)/2: edge_u/2);
			coords[2*j-1][0] = vcoords[vertex_v][0];
			coords[2*j-1][1] = vcoords[vertex_v][1];
			

			/* identify the coordinates of the midpoint of the edge after edge_u in the cycle */
			int edge_w;
			int vertex_w;
			if (j+1 <= cycle[region][0])
			{
				edge_w = cycle[region][j+1];
				vertex_w = abs(edge_w)+num_crossings;
				coords[2*j][0] = vcoords[vertex_w][0];
				coords[2*j][1] = vcoords[vertex_w][1];
			}
			else
			{
				edge_w = cycle[region][1];
				vertex_w = abs(edge_w)+num_crossings;
				coords[0][0];
				coords[0][1];
			}
			
			/* move around the cycle */
			edge_u = edge_w;
			vertex_u = vertex_w;  // only needed for debugging
		}	

		region_coords = coords;
	}		
}

void shrink_triangulation (matrix<double>& shrunken_triangulation, int num_triangles, matrix<double>& region_coords)
{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "shrink_triangulation: maximal region triangulation contains " << num_triangles << " triangles" <<endl;

	if (num_triangles == 2)
	{
		double u_x = region_coords[0][0];
		double u_y = region_coords[0][1];
		double v_x = region_coords[1][0];
		double v_y = region_coords[1][1];
		double w_x = region_coords[2][0];
		double w_y = region_coords[2][1];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "shrink_triangulation: first triangle: u = (" << u_x << "," << u_y << "), v = (" << v_x << "," << v_y << "), w = (" << w_x << "," << w_y << ")" << endl;
			
		shrink_triangle(shrunken_triangulation, 0, u_x, u_y, v_x, v_y, w_x, w_y);

		v_x = region_coords[3][0];
		v_y = region_coords[3][1];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "shrink_triangulation: second triangle: u = (" << u_x << "," << u_y << "), v = (" << v_x << "," << v_y << "), w = (" << w_x << "," << w_y << ")" << endl;
			
		shrink_triangle(shrunken_triangulation, 1, u_x, u_y, v_x, v_y, w_x, w_y);
			
	}
	else
	{			
		double barycentre_x = region_coords[num_triangles][0];
		double barycentre_y = region_coords[num_triangles][1];
		double u_x = region_coords[0][0];
		double u_y = region_coords[0][1];

		double v_x,v_y,w_x,w_y;
			
		for (int j=0; j < num_triangles/2; j++)
		{
			v_x = region_coords[2*j+1][0];
			v_y = region_coords[2*j+1][1];

			if (j+1 < num_triangles/2)
			{
				w_x = region_coords[2*(j+1)][0];
				w_y = region_coords[2*(j+1)][1];
			}
			else
			{
				w_x = region_coords[0][0];
				w_y = region_coords[0][1];
			}
			
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "shrink_triangulation:   first triangle corresponding to edge " << j << ": u = (" << u_x << "," << u_y << "), v = (" << v_x << "," << v_y << "), w = (" << w_x << "," << w_y << ")" << endl;

			shrink_triangle(shrunken_triangulation, 2*j, u_x, u_y, v_x, v_y, w_x, w_y);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "shrink_triangulation:   second triangle corresponding to edge " << j << ": u = (" << u_x << "," << u_y << "), v = (" << barycentre_x << "," << barycentre_y << "), w = (" << w_x << "," << w_y << ")" << endl;

			shrink_triangle(shrunken_triangulation, 2*j+1, u_x, u_y, barycentre_x, barycentre_y, w_x, w_y);
							
			/* move around the cycle */
			u_x = w_x;
			u_y = w_y;
		}
	}
}

void shrink_triangle (matrix<double>& shrunken_triangulation, int index, double u_x, double u_y, double v_x, double v_y, double w_x, double w_y)
{
	double barycentre_x = (u_x+v_x+w_x)/3;
	double barycentre_y = (u_y+v_y+w_y)/3;
	
	shrunken_triangulation[index][0] = barycentre_x + region_shrinking_factor * (u_x - barycentre_x);
	shrunken_triangulation[index][1] = barycentre_y + region_shrinking_factor * (u_y - barycentre_y);
	shrunken_triangulation[index][2] = barycentre_x + region_shrinking_factor * (v_x - barycentre_x);
	shrunken_triangulation[index][3] = barycentre_y + region_shrinking_factor * (v_y - barycentre_y);
	shrunken_triangulation[index][4] = barycentre_x + region_shrinking_factor * (w_x - barycentre_x);
	shrunken_triangulation[index][5] = barycentre_y + region_shrinking_factor * (w_y - barycentre_y);
}

/* adjacent_region identifes the region adjacent to the given edge in the given region.  It records the index of edge
   in the boundary cycle of the adjacent region in index.
*/
int adjacent_region(matrix<int>& cycle, int num_cycles, int region, int edge, int& index)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "adjacent_region: locating edge " << edge << " in region adjacent to " << region << endl;

	/* find edge in cycle but not in region */
	bool found = false;
	int adjacent_region = -1;
	
	for (int i=0; i< num_cycles && !found; i++)
	{
		if (i == region)
			continue;
			
		for (int j=1; j<= cycle[i][0] && !found; j++)
		{
			if (abs(cycle[i][j]) == edge)
			{
				adjacent_region = i;
				index = j-1;
				found = true;
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "adjacent_region: found in region " << adjacent_region << endl;

	return adjacent_region;
}

/* edge_weight returns the weight of the given edge within the given region, that is the area of the
   other turning cycle within which edge resides.
*/
double edge_weight(matrix<int>& cycle, int num_cycles, int region, int edge, vector<double>& region_area)
{
	int index; // only needed for the call to adjacent_region
	int other_region = adjacent_region(cycle, num_cycles, region, edge, index);
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_weight:   opposite edge " << edge << " also lies in cycle " << other_region << " of weight " << region_area[other_region] << endl;

	return region_area[other_region];
}

