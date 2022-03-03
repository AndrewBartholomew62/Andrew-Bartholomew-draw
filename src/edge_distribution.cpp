/**********************************************************************************
                            Edge distribution placement
                                   January 2019

This approach seeks to expand regions of a diagram interconnected by short edges of the 
trinagulation after circle packing by contracting long edges.  To identify a threshold between 
what is considered a long or short edge, we plot a histogram of the edge length distribution, 
counting the number of edges in each of a number of divisions (default 20) between the
minimum and maximum length.  The gaps in this distribution, where the edge count is zero
indicates a potential threshold between what may be considered short or long.  We select an
"active gap" and choose the middle point of that gap as our threshold.  The active gap is set
initially to the first (lowest divsion index) gap.  Unless instructed always to use the first gap,
the algorithm seeks a gap longer than the first that remains "acceptable".  By acceptable we mean that
using this gap as the active gap results in moving a vertex (described below).

Having identified a threshold value we note whether a vertex is incident with a long or short edge,
or both.  We can then speak of a long or short vertex, understanding that a vertex may be both long 
and short.  It is these latter vertices that we wish to reposition, provided they do not lie in the 
infinite turning cycle.

We move a vertex u that is both long and short, and not in the infinite turning cycle, towards the 
centre of gravity of u, together with any other vertex v in the flower of u that is long but not 
short.  Note that since v is not short the edge uv must be long.

The algorithm terminates when the edge distribution histogram contains no gaps, since then the
edges are all roughly of the same length.

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
extern bool USE_FIRST_GAP_AS_ACTIVE;

extern bool USE_ALL_LONG_GAPS;   // temporary - remove after testing

extern int metapost_coordinate_scale_factor;
extern int placement_iteration_tracking_step;
extern int plot_steps;

extern double EPSILON;
extern double ZERO_THRESHOLD;

extern double triangulation_plot_threshold;
extern float  edge_distribution_shift_factor;


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

void set_long_and_short_vertices(vector<bool>& long_vertex, vector<bool>& short_vertex,matrix<int>& flowers, matrix<double>& vcoords, int num_type123_vertices, double threshold);


void edge_distribution_placement(metapost_control& mp_control, generic_code_data& code_data, matrix<int>& cycle, int num_cycles, 
                                int num_left_cycles, int infinite_region, int num_iterations, string title)
{
	
cout << "\n\nUsing edge distribution placement" << endl;

	int num_crossings = code_data.num_crossings;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	if (USE_SMALL_REGION_SHRINKING_PLACEMENT)
			debug << "edge_distribution_placement: using small region placement variant" << endl;
	debug << "edge_distribution_placement: code_data: ";
	write_code_data(debug,code_data);
	debug << endl;
	print_code_data(code_data, debug, "edge_distribution_placement: ");
	debug << "edge_distribution_placement: num_iterations = " << num_iterations << endl;
	if (TRACK_PLACEMENT_ITERATION)
		debug << "edge_distribution_placement: placement_iteration_tracking_step = " << placement_iteration_tracking_step << endl;
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
	debug << "edge_distribution_placement: read triangulation data from " << triangulation_output_file << endl;

	ifstream triangulation(triangulation_output_file);
	if(!triangulation)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: could not open " << triangulation_output_file << endl;
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
	debug << "edge_distribution_placement: nodecount = " << nodecount << endl;
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
	debug << "edge_distribution_placement: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
	}
	
	int num_type123_vertices = beta-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "edge_distribution_placement: number of type 1,2 & 3 vertices num_type123_vertices = " << num_type123_vertices << endl;

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
	debug << "edge_distribution_placement: node " << node << ", petal count = " << count << ": " << first_petal;
				
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
    debug << "edge_distribution_placement: type_3_vertex: ";
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
	debug << "edge_distribution_placement: vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "edge_distribution_placement:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
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
	debug << "edge_distribution_placement: boundary vertex centre of gravity = (" << boundary_cog_x << "," << boundary_cog_y << ")" << endl;
	
		for (int i=num_type123_vertices; i< nodecount; i++)
		{
			vcoords[i][0] = (vcoords[i][0]-boundary_cog_x)*boundary_vertex_retraction_factor + boundary_cog_x;
			vcoords[i][1] = (vcoords[i][1]-boundary_cog_y)*boundary_vertex_retraction_factor + boundary_cog_y;
		}
	}

    matrix<double> initial_vcoords = vcoords;
	matrix<double> two_back_vcoords = vcoords;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "edge_distribution_placement: initial_vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "edge_distribution_placement:   v" << i << " = ("<< initial_vcoords[i][0] << ", " << initial_vcoords[i][1] << ")" << endl;
}


		/* When deciding whether to move a vertex, we need to understand which type 1 or type 2 
		   vertices lie in the infinite turning cycle, so we note these vertices here.
		   we make the vector long enough to store false for type 3 vertices to simplify
		   the code deciding whether to move a vertex or not.
		*/
		vector<bool> infinite_type_12_vertex_flag(num_type123_vertices); // initializes to false
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

	/* Main iterative loop */	   
	bool algorithm_converged = false;
	
	for (int i=0; i< num_iterations; i++)
	{

//	    cout << "\nedge_distribution_placement: iteration " << i << endl;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: edge distribution iteration " << i << endl;
	
		/* Create a list of the edges in the triangulation together with the relevant information about that edge, i.e
		   it's vertices and its length and the histogram division to which it belongs.
		*/
		list<double> edge_list;
		list<edge_info> info_list;
	
		matrix<int> edge (num_type123_vertices,num_type123_vertices); // records whether we've drawn an edge between neighbours
		for (int i=0; i< num_type123_vertices; i++)
		{		    
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_distribution_placement:   vertex "<< i << endl;
			
			int num_neighbours = flowers[i][0];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_distribution_placement:   num_neighbours " << num_neighbours << endl;
			
			for (int j=1; j<= num_neighbours; j++)
			{					
				int neighbour = flowers[i][j];
				neighbour--;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_distribution_placement:     neighbour "<< neighbour;
				
				if (neighbour < num_type123_vertices && edge[i][neighbour] == 0)
				{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " edge not processed"<< endl;
	
					edge_info next_edge;
					next_edge.vertex_0 = i;
					next_edge.vertex_1 = neighbour;
					
					/* evaluate the length of the edge */
					double delta_x = vcoords[i][0] - vcoords[neighbour][0];
					double delta_y = vcoords[i][1] - vcoords[neighbour][1];
					double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);
	
					edge_list.push_back(mod_delta);

					next_edge.vertex_0 = i;
					next_edge.vertex_1 = neighbour;
					next_edge.length = mod_delta;
					
					info_list.push_back(next_edge);
											
					edge[i][neighbour] = 1;
					edge[neighbour][i] = 1;
				}
				else if (neighbour < num_type123_vertices)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " already processed"<< endl;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " is a type 4 vertex"<< endl;
				}			
			}
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: triangulation contains " << edge_list.size() << " edges" << endl;
	
		list<double>::iterator eptr = edge_list.begin();
		double max_len = *eptr;
		double min_len = *eptr;
		eptr++;
		
		while (eptr != edge_list.end())
		{
			if (*eptr > max_len)
				max_len = *eptr;			
			else if (*eptr < min_len)
				min_len = *eptr;
				
			eptr++;
		}
	
		double division_width = (max_len - min_len)/plot_steps;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "edge_distribution_placement: minimum triangulation edge length = " << min_len << ", maximum triangulation edge length = " << max_len << endl;
	debug << "edge_distribution_placement: distribution division width = " << division_width << endl;
}

		/* since (max_len - min_len)/division_width = plot_steps we will always have the division corresponding 
		   to the maximal edge equal to plot_steps, so need to allow for this in the edge_distribution vector
		*/
		vector<int> edge_distribution(plot_steps+1);
		eptr = edge_list.begin();
		int edge_count = 0;

		list<edge_info>::iterator iptr = info_list.begin();

		
		while (eptr != edge_list.end())
		{
	
			int division = static_cast<int>((*eptr-min_len)/division_width);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: edge " << edge_count << ", length = " << *eptr << ", division " << *eptr/division_width << ", adjusted to " << division << endl;
		
			edge_distribution[division]++;
			eptr++;
			
			iptr->index = edge_count;
			iptr->division_index = division;
			iptr++;

			edge_count++;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "edge_distribution_placement: triangulation edge distribution: ";
	for (int i=0; i< plot_steps+1; i++)
		debug << edge_distribution[i] << " ";
	debug << endl;
}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "edge_distribution_placement: edge information list: " << endl;
	iptr = info_list.begin();
	while (iptr != info_list.end())
	{
		debug << "edge_distribution_placement:   edge " << iptr->index << endl;
		debug << "edge_distribution_placement:     vertex_0 = " << iptr->vertex_0 << endl;	
		debug << "edge_distribution_placement:     vertex_1 = " << iptr->vertex_1 << endl;	
		debug << "edge_distribution_placement:     length = " << iptr->length << endl;	
		debug << "edge_distribution_placement:     division_index = " << iptr->division_index << endl;	
		iptr++;
	}
}

		/* identify, if it exists, the "active gap" in the distribution that we will use to determine which edges 
		   are long or short, which in turn determines which vertices will be moved.  
		   
		   There are two possible approaches:
		   
		   a) always use the first gap
		   b) consider other gaps if they are longer than the first and are "acceptable"
		   
		   Here, accceptable means that there is an edge to the right of the gap for which one of its vertices is long but 
		   not short (this is a vertex that will allow a short vertex to move).  If there are no such vertices, this gap
		   will not cause any vertices to move and is therefore unacceptable.  
		   
		   It may be that the first gap is actually unacceptable, in which case the algorithm terminates.
		   
		   We start by preparing for approach b) then override the active gap if we're asked to use option a)

		   We know that entries 0 and plot_steps of edge_distribution are both non-zero, so we can count the number of gaps 
		   by looking for a zero followed by a non_zero.
		*/
		int num_gaps = 0;

		/* long_vertex and short_vertex are used to identify acceptable gaps and to move vertices once we have identified a 
		   gap to use.
		*/
		vector<bool> long_vertex(num_type123_vertices);
		vector<bool> short_vertex(num_type123_vertices);

	   
		for (int i=0; i< plot_steps; i++)
		{
			if (edge_distribution[i] == 0 && edge_distribution[i+1] != 0)
				num_gaps++;	
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: num_gaps =  " << num_gaps << endl;

		if (num_gaps > 0)
		{
			vector<int> gap_length(num_gaps);
			vector<int> gap_start(num_gaps);
			
			int index = 0;
			for (int i=0; i< num_gaps; i++)
			{
				while(edge_distribution[index] != 0)
					index++;
				
				gap_start[i] = index;
				
				int count = 0;
		
				while(edge_distribution[index] == 0)
				{
					count++;
					index++;
				}
				
				gap_length[i] = count;
			}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "edge_distribution_placement: gap lengths =  ";
	for (int i=0; i< num_gaps; i++)
		debug << gap_length[i] << " ";
	debug << endl;

	debug << "edge_distribution_placement: gap starting indices =  ";
	for (int i=0; i< num_gaps; i++)
		debug << gap_start[i] << " ";
	debug << endl;

}

			int active_gap_index = 0;

			if (!USE_FIRST_GAP_AS_ACTIVE)
			{
				for (int i=1; i< num_gaps; i++)
				{
					if (gap_length[i] > gap_length[active_gap_index])
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement: gap " << i << " length longer than gap " << active_gap_index << endl;

						/* it looks like we have a better gap to use but check that it is acceptable */
						int middle_division_index = gap_start[i] + gap_length[i]/2;
						double division_threshold = min_len + middle_division_index*division_width;
					
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "edge_distribution_placement: gap " << i << " middle division index = " << middle_division_index << endl;
	debug << "edge_distribution_placement: gap " << i << " division_threshold = " << division_threshold << endl;
}					
						set_long_and_short_vertices(long_vertex,short_vertex,flowers,vcoords,num_type123_vertices,division_threshold);
						   
						bool acceptable_edge_found = false;
						list<edge_info>::iterator iptr = info_list.begin();
						while( iptr != info_list.end() )
						{
							if (iptr->division_index > gap_start[i])
							{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement:   edge " << iptr->index << " lies to the right of gap " << i;				

//								if (!infinite_type_12_vertex_flag[iptr->vertex_0] || !infinite_type_12_vertex_flag[iptr->vertex_1])
								if ( USE_ALL_LONG_GAPS || // temporary, will be removed, used to force all long gaps to be acceptable
								(long_vertex[iptr->vertex_0] && !short_vertex[iptr->vertex_0]) || (long_vertex[iptr->vertex_1] && !short_vertex[iptr->vertex_1]) )
								{
									acceptable_edge_found = true;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " and is acceptable" << endl;				

									break;
								}
								else
								{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " but is not acceptable" << endl;				
								}
							}
							
							iptr++;
						}
						
						if (acceptable_edge_found)
						{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement: gap " << i << " is acceptable as a threshold gap" << endl;				
							active_gap_index = i;
						}
						else
						{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement: gap " << i << " not acceptable as a threshold gap" << endl;				
						}
					}
					else
					{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement: gap " << i << " length not longer than gap " << active_gap_index << endl;
					}
				}
			}
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: active_gap_index = " << active_gap_index << ", gap length " << gap_length[active_gap_index] << endl;

			int middle_division_index = gap_start[active_gap_index] + gap_length[active_gap_index]/2;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: middle_division_index = " << middle_division_index << endl;
	
			triangulation_plot_threshold = min_len + middle_division_index*division_width;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: threshold = " << triangulation_plot_threshold << endl;

		}
		else
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: edge distribution placement algorithm converged after iteration " << i << ", no gaps found in triangulation distribution." << endl;

			cout << "edge distribution placement algorithm converged after iteration " << i << ", no gaps found in triangulation distribution" << endl;
	
			algorithm_converged = true;
			break;
		}

		/* Adjust vertex coordinates.  Start by identifying long and short vertices according to the selected triangulation_plot_threshold*/
		set_long_and_short_vertices(long_vertex,short_vertex,flowers,vcoords,num_type123_vertices,triangulation_plot_threshold);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: adjusting vertex placement" << endl;
		
		int moved_vertex_count = 0;
		
		for (int i=0; i< num_type123_vertices; i++)			
		{		    
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement:   vertex "<< i;
			
			if (i < num_type123_vertices)
			{
				if (long_vertex[i])
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " is long";
	
					if (short_vertex[i])
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << ", short";
						if (!infinite_type_12_vertex_flag[i])
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " and does not lie in the infinite turning cycle." << endl;
							
							double cog_x = vcoords[i][0];
							double cog_y = vcoords[i][1];
							int count = 1;

							int num_neighbours = flowers[i][0];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement:   num_neighbours " << num_neighbours << endl;
			
							for (int j=1; j<= num_neighbours; j++)
							{					
								int neighbour = flowers[i][j];
								neighbour--;
								
								if (long_vertex[neighbour] && !short_vertex[neighbour])
								{
									cog_x += vcoords[neighbour][0];
									cog_y += vcoords[neighbour][1];
									count++;
									
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement:     neighbour " << neighbour << " contributes to COG" << endl;
								}
								else
								{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "edge_distribution_placement:     neighbour " << neighbour << " does not contribute to COG" << endl;
								}
							}

							if (count > 1)
							{
								moved_vertex_count++;

								cog_x/=count;
								cog_y/=count;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement:   count = " << count << ", move vertex towards COG (" << cog_x << "," << cog_y << ")" << endl;
								
								vcoords[i][0] = edge_distribution_shift_factor * cog_x + (1-edge_distribution_shift_factor) * vcoords[i][0];
								vcoords[i][1] = edge_distribution_shift_factor * cog_y + (1-edge_distribution_shift_factor) * vcoords[i][1];
								
							}
							else
							{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement:   cound not find any COG contributors, vertex not moved" << endl;
							}
						}
						else
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " but lies in the infinite turning cycle" << endl;
						}
					}
					else
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " but not short" << endl;
					}
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " is not long " << endl;
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " is a type 4 vertex" << endl;
			}			
		}
		
		if (TRACK_PLACEMENT_ITERATION && i % placement_iteration_tracking_step == placement_iteration_tracking_step-1)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: tracking region shrinking placement at iterative step i = " << i << endl;
	
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

		if (moved_vertex_count == 0)
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: edge distribution placement algorithm converged after iteration " << i << ", no vertices moved in this iteration" << endl;

			cout << "edge distribution placement algorithm converged after iteration " << i << ", no vertices moved in this iteration" << endl;
	
			algorithm_converged = true;
			break;
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
				cout << "Edge distribution placement stopping after " << i << " iterations, oscillation convergence detected" << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "edge_distribution_placement: stopping after " << i << " iterations, oscillation convergence detected" << endl;
				
				algorithm_converged = true;
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
	debug << "edge_distribution_placement: could not open " << circlepack_output_file << endl;
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
		cout << "edge distribution placement algorithm terminated after specified " << num_iterations << " iterations" << endl;
	
	output.close();
}

void set_long_and_short_vertices(vector<bool>& long_vertex, vector<bool>& short_vertex,matrix<int>& flowers, matrix<double>& vcoords, int num_type123_vertices, double threshold)
{

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "set_long_and_short_vertices: threshold = " << threshold << endl;
	
	matrix<int> edge (num_type123_vertices,num_type123_vertices); // records whether we've drawn an edge between neighbours

	/* clear the vertex flags */
	for (int i=0; i< num_type123_vertices; i++)			
		long_vertex[i] = short_vertex[i] = false;

	for (int i=0; i< num_type123_vertices; i++)			
	{		    
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "set_long_and_short_vertices: vertex "<< i;
			
		int num_neighbours = flowers[i][0];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << ", num_neighbours " << num_neighbours << endl;
			
		for (int j=1; j<= num_neighbours; j++)
		{					
			int neighbour = flowers[i][j];
			neighbour--;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "set_long_and_short_vertices:   neighbour "<< neighbour;
				
			if (neighbour < num_type123_vertices) // && edge[i][neighbour] == 0)
			{

//if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
//	debug << " edge not processed"<< endl;
	
				/* evaluate the length of the edge */
				double delta_x = vcoords[i][0] - vcoords[neighbour][0];
				double delta_y = vcoords[i][1] - vcoords[neighbour][1];
				double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);

				if (mod_delta > threshold)
				{
					long_vertex[i] = true;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " identifies vertex as long" << endl;
				}
				else if (mod_delta < threshold)
				{
					short_vertex[i] = true;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " identifies vertex as short" << endl;
				}
				else  // should never get here!
				{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " doesn't identify anything!" << endl;
				}
											
			}
			else
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << " is a type 4 vertex"<< endl;
			}			
		}
	}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_long_and_short_vertices:  long_vertex: ";
	for (int i=0; i< num_type123_vertices; i++)
		debug << (long_vertex[i]? "T ": "F ");
	debug << "\nset_long_and_short_vertices:  short_vertex: ";
	for (int i=0; i< num_type123_vertices; i++)
		debug << (short_vertex[i]? "T ": "F ");
	debug << endl;
}		

}
