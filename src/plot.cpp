/**********************************************************************************
                            Plot triangulation distribution
                                   January 2019

This is expoloratory code to plot the length of the triangulation edges, not including those
incident with a type 4 vertex, as a histogram.  We count the number of edges in each of 100
intervals between the minimum and maximum triangulation edge length.

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

extern double triangulation_plot_threshold;

/********************* Function prototypes ***********************/
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, matrix<double> vcoords, vector<double> vertex_radius);
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, 
                    char const* circlepack_output_file, char const* circlepack_output_savefile);
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


void plot_triangulation_distribution(int num_steps)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "plot_triangulation_distribution: plotting triangulation with " << num_steps << " steps" << endl;

	/* Read triangulation data from triangulation_output_file */
	int nodecount;
	int alpha;
	int beta;
	int gamma;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "plot_triangulation_distribution: read triangulation data from " << triangulation_output_file << endl;

	ifstream triangulation(triangulation_output_file);
	if(!triangulation)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "plot_triangulation_distribution: could not open " << triangulation_output_file << endl;
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
	debug << "plot_triangulation_distribution: nodecount = " << nodecount << endl;
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
	debug << "plot_triangulation_distribution: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
	}
	
	int num_type123_vertices = beta-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "plot_triangulation_distribution: number of type 1,2 & 3 vertices num_type123_vertices = " << num_type123_vertices << endl;

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

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "plot_triangulation_distribution: node " << node << ", petal count = " << count << ": " << first_petal;
				
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
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "plot_triangulation_distribution: vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "plot_triangulation_distribution:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
}


	list<double> edge_list;

	matrix<int> edge (num_type123_vertices,num_type123_vertices); // records whether we've drawn an edge between neighbours
	for (int i=0; i< num_type123_vertices; i++)
	{		    
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "plot_triangulation_distribution:   vertex "<< i << endl;
			
		int num_neighbours = flowers[i][0];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "plot_triangulation_distribution:   num_neighbours " << num_neighbours << endl;
			
		for (int j=1; j<= num_neighbours; j++)
		{					
			int neighbour = flowers[i][j];
			neighbour--;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "plot_triangulation_distribution:     neighbour "<< neighbour;
				
			if (neighbour < num_type123_vertices && edge[i][neighbour] == 0)
			{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " edge not processed"<< endl;
	
				/* evaluate the length of the edge */
				double delta_x = vcoords[i][0] - vcoords[neighbour][0];
				double delta_y = vcoords[i][1] - vcoords[neighbour][1];
				double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);

				edge_list.push_back(mod_delta);
										
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
	debug << "plot_triangulation_distribution: triangulation contains " << edge_list.size() << " edges" << endl;
	
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

	double division_width = (max_len - min_len)/num_steps;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "plot_triangulation_distribution: minimum triangulation edge length = " << min_len << ", maximum triangulation edge length = " << max_len << endl;
	debug << "plot_triangulation_distribution: distribution division width = " << division_width << endl;
}

	/* since (max_len - min_len)/division_width = num_steps we will always have the division corresponding 
	   to the maximal edge equal to num_steps, so need to allow for this in the edge_distribution vector
	*/
	vector<int> edge_distribution(num_steps+1);
	eptr = edge_list.begin();
	int edge_count = 0;
	
	while (eptr != edge_list.end())
	{

		int division = static_cast<int>((*eptr-min_len)/division_width);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "plot_triangulation_distribution: edge " << edge_count << ", length = " << *eptr << ", division " << *eptr/division_width << ", adjusted to " << division << endl;
		
		edge_distribution[division]++;
		eptr++;
		edge_count++;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "plot_triangulation_distribution: triangulation edge distribution: ";
	for (int i=0; i< num_steps+1; i++)
		debug << edge_distribution[i] << " ";
	debug << endl;
}

	/* calculate the largest gap, we know that entries 0 and num_steps of edge_distribution are both non-zero, so 
      can count the number of gaps by looking for a zero followed by a non_zero.
	*/
	int num_gaps = 0;
   
	for (int i=0; i< num_steps; i++)
	{
		if (edge_distribution[i] == 0 && edge_distribution[i+1] != 0)
			num_gaps++;	
	}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "plot_triangulation_distribution: num_gaps =  " << num_gaps << endl;

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
	debug << "plot_triangulation_distribution: gap lengths =  ";
	for (int i=0; i< num_gaps; i++)
		debug << gap_length[i] << " ";
	debug << endl;

	debug << "plot_triangulation_distribution: gap starting indices =  ";
	for (int i=0; i< num_gaps; i++)
		debug << gap_start[i] << " ";
	debug << endl;

}

		int max_gap_index = 0;
		for (int i=1; i< num_gaps; i++)
		{
			if (gap_length[i] > gap_length[max_gap_index])
				max_gap_index = i;
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "plot_triangulation_distribution: max_gap_index = " << max_gap_index << ", gap length " << gap_length[max_gap_index] << endl;

		int middle_division_index = gap_start[max_gap_index] + gap_length[max_gap_index]/2;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "plot_triangulation_distribution: middle_division_index = " << middle_division_index << endl;
	
		triangulation_plot_threshold = min_len + middle_division_index*division_width;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "plot_triangulation_distribution: threshold = " << triangulation_plot_threshold << endl;

	}
}

