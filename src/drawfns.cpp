/************************************************************************
                  Support functions for draw

                  
string invstr(string& str)
double badness(string vertex_file, generic_code_data& code_data)
bool first_occurrence(matrix<int>& code_table, matrix<int>& infinite_cycle_first_visit, int crossing, int position)
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, matrix<double> vcoords, vector<double> vertex_radius, vector<int>* auxiliary_data=0)
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, 
                    char const* circlepack_output_file, char const* circlepack_output_savefile, vector<int>* auxiliary_data=0)
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles)
void print (metapost_control& mp_control, ostream& os, string prefix)
void set_edge_directions(generic_code_data& code_data, matrix<double>& vcoords, matrix<double>& edge_direction)
void calculate_theta_angles(int e1, int e2, int e3, int e4, int v,double& theta_1, double& theta_2, double& theta_3, double& theta_4,matrix<double>& vcoords)
double edge_length(double px,double py,double qx,double qy)
double cosine_rule_angle (double a, double b, double c)
double edge_argument(int v, int u, matrix<double>& vcoords)
double edge_argument(double x_delta, double y_delta)
**************************************************************************/
using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <complex>

/********************* External variables ***********************/
extern string		title;

extern ifstream     input;
extern ofstream 	output;
extern ofstream     debug;

extern char const* triangulation_output_file;
extern char const* triangulation_shrink_file;

extern double badness_threshold; 
extern int metapost_coordinate_scale_factor;
extern int metapost_hyperbolic_scale_factor;
extern bool DRAW_IN_HYPERBOLIC_DISC;
extern bool INCLUDE_BOUNDARY_VERTICES;
extern bool TRACK_PLACEMENT_ITERATION;
extern bool CHECK_INNER_HULL_CALCULATION;
extern bool IMMERSION_CROSSINGS_AT_RIGHT_ANGLES;
extern bool USE_FORCE_DIRECTED_PLACEMENT;
extern bool USE_CENTRE_OF_GRAVITY_PLACEMENT;
extern bool USE_REGION_SHRINKING_PLACEMENT;
//extern bool mp_control.state_smoothed;
extern int check_inner_hull_vertex;
extern double average_triangulation_edge_length;
extern float average_triangulation_length_factor;
//extern bool MAGNIFY_SMALL_CIRCLES;
extern float magnification_factor;
extern double pi;
extern double two_pi;
extern double ZERO_THRESHOLD;
extern double triangulation_plot_threshold;
extern float metapost_path_tension;

#include <util.h>
#include <matrix.h>
#include <draw.h>
#include <gauss-orientation.h>

/********************* Function prototypes ***********************/
void hyperbolic_representation(matrix<double>& centre, vector<double>& radius);
void mobius_transform_to_origin(complex<double>& a, complex<double>& b, complex<double> z);
complex<double> mobius_transformation(complex<double> a, complex<double> b, complex<double> z);
complex<double> inv_mobius_transformation(complex<double> a, complex<double> b, complex<double> z);
void set_edge_directions(generic_code_data& code_data, matrix<double>& vcoords, matrix<double>& edge_direction);
void calculate_theta_angles(int e1, int e2, int e3, int e4, int v,double& theta_1, double& theta_2, double& theta_3, double& theta_4,matrix<double>& vcoords);
double edge_length(double px,double py,double qx,double qy);
double cosine_rule_angle (double a, double b, double c);
double edge_argument(int v, int u, matrix<double>& vcoords);
double edge_argument(double x_delta, double y_delta);
void rotate_metapost (matrix<double>& coords, int num_coords_to_rotate, metapost_control& mp_control);
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data);
void add_virtual_crossing(int location, vector<int>& fringe, vector<int>& num_virtual_crossings_on_gauss_arc, list<gc_pc_xlabels>& virtual_crossings);
void assign_gauss_arc_direction(int arc_1, int arc_2, matrix<int>& gauss_code_table, int crossing, vector<int>& gauss_arc_direction);
string direction_to_string(int edge, vector<int>& gauss_arc_direction);
void remove_fringe_edge(int edge, vector<int>& current_fringe);


/* invstr creates the inverse of src and returns it to the call.
   The string src must be composed of si, -si, and ti components only.
   The inverse of ti is again ti.
*/
string invstr(string& str)
{
    char*	src = c_string(str);
    char*	loc_buf = new char[2*str.size() + 1]; // a slight overestimate, but always enough
    char*	mark = src;
    bool	not_finished = true;

    /* start at the end of the string and work back */
    char* cptr = strchr(src,0);
    if (strlen(src))
	    cptr--;  /* if src is the null string, cptr = src */

    char* lptr = loc_buf;
    do
    {
		if (isdigit(*cptr))
	    	cptr--;
		else
		{
		    if (*cptr == 's')
	    	{
				if (cptr !=src && *(cptr-1) == '-')
			    	mark = cptr-1;
				else
				{
		    		*lptr++ = '-';
		    		mark = cptr;
				}
				*lptr++ = 's';
	    	}
	    	else if (*cptr == 't')
	    	{
				mark = cptr;
				*lptr++ = 't';
		    }

		    cptr++;
		    while (isdigit(*cptr))
				*lptr++ = *cptr++;
	    	if (mark == src)
				not_finished = false;
	    	else
				cptr = mark-1;
		}
    } while (not_finished);

    *lptr = '\0';
    return string(loc_buf);
}

/* the badness function looks at consecutive vertices z_{n+1} and z_n in the vertex sequence
   and evaluates |z_{n+1}-z_n| and returns min{100, |z_{i+1} - z_i|}.  For links this is not
   perfect, since at component boundaries the vertex sequence does not reflect the path
   drawn by metapost but the distances between the last and first vertex at a component
   boundary is typically large, so it should not affect the value returned by the function.
*/
double badness(string vertex_file, generic_code_data& code_data)
{
	
	/* create the vertex sequence for writing the metapost path using the terminating 
	   crossing information in code_data for each edge.  The sequence starts with the 
	   type 2 vertex preceeding crossing zero, which is the first type 2 vertex in the 
	   numbering scheme.	       
	*/
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	int vertex_sequence_length = 2*num_edges;
	vector<int>& term_crossing = code_data.term_crossing;
	vector<int> vertex_sequence(vertex_sequence_length);
	for (int i=0; i< num_edges; i++)
	{
		vertex_sequence[2*i] = num_crossings+(i%num_edges); //the type_2_vertex at the midpoint of edge i
		vertex_sequence[2*i+1] = term_crossing[i];
	}

	double x1, x2, y1, y2, z, w;
	double worst=badness_threshold;
	
	ifstream input;
	input.open(vertex_file.c_str());
	
	if (input)
	{

		int num_type1234_vertices;
		input >> num_type1234_vertices;
		
		matrix<double> drawcoords(vertex_sequence_length,2);
		for (int i=0; i< vertex_sequence_length; i++)
		{
			input >> drawcoords[i][0];
			input >> drawcoords[i][1];
		}

		for (int i=1; i<vertex_sequence_length; i++)
		{
			x1 = drawcoords[i-1][0] * metapost_coordinate_scale_factor;
			x2 = drawcoords[i][0] * metapost_coordinate_scale_factor;
			y1 = drawcoords[i-1][1] * metapost_coordinate_scale_factor;
			y2 = drawcoords[i][1] * metapost_coordinate_scale_factor;
			w = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			z = sqrt(w);
			
			if (z < worst) 
				worst = z;
		}
	}
	else
	{
		cout << "\nError draw::badness failed to open " << vertex_file << endl;
		exit(0);
	}
	input.close();
	return worst;
}

/* first_occurrence detemines whether the given position at the given crossing is the first occurrence
   of the infinte region at the crossing or not.  It does this by identifying the two edges incident with 
   the given position and then checks to see if both are recorded in infinite_cycle_first_visit for that 
   crossing.  We have to check both edges incident with the region to handle Reidemeister I loops.
*/
bool first_occurrence(generic_code_data& code_data, matrix<int>& infinite_cycle_first_visit, int crossing, int position)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "first_occurrence: crossing = " << crossing << " position = " << position << endl;
    
    matrix<int,int>& code_table = code_data.code_table;
    vector<int>& first_edge_on_component = code_data.first_edge_on_component;
    vector<int>& num_component_edges = code_data.num_component_edges;

	bool first;
    
	int peer_component = code_table[COMPONENT][(code_table[OPEER][crossing]-1)/2];
	int successor = 2*crossing+1;
	int peer_successor = (code_table[OPEER][crossing]-first_edge_on_component[peer_component]+1)%num_component_edges[peer_component]+first_edge_on_component[peer_component];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "first_occurrence: peer_component = " << peer_component << " peer_successor = " << peer_successor << endl;


	
	int search_edge_1;
	int search_edge_2;
	
	switch (position)
	{
		case 0: search_edge_1 = 2*crossing; // naming edge
				search_edge_2 = code_table[OPEER][crossing];
				break;
		case 1: if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
		        {
					search_edge_1 = 2*crossing; // naming edge
					search_edge_2 = peer_successor;
				}
				else
				{
					search_edge_1 = code_table[OPEER][crossing];
					search_edge_2 = successor;
				}
				break;
		case 2: search_edge_1 = successor;
				search_edge_2 = peer_successor;
				break;
		case 3: if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
		        {
					search_edge_1 = code_table[OPEER][crossing];
					search_edge_2 = successor;
				}
				else
				{
					search_edge_1 = 2*crossing; // naming edge
					search_edge_2 = peer_successor;
				}
				break;
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "first_occurrence: search_edge_1 = " << search_edge_1 << ", search_edge_2 = " << search_edge_2 << endl;

	if ( (infinite_cycle_first_visit[crossing][0] == search_edge_1 && infinite_cycle_first_visit[crossing][1] == search_edge_2) ||
	     (infinite_cycle_first_visit[crossing][0] == search_edge_2 && infinite_cycle_first_visit[crossing][1] == search_edge_1)
	   )
		first = true;
	else
		first = false;
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "first_occurrence: first_occurrence = " << (first? "true" : "false") << endl;

	return first;
}

/* The auxiliary_data parameter in write_metapost is used to convey:
     - A specific smoothed state if a single state is to be drawn, or 
     - The gauss_crossing map when a Gauss_code (or PD data) is being drawn, for use with the gauss-crossings option, or
     - An edge list for a Hamiltonian circuit, for use with the mp_control.hamiltonians option
*/
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, matrix<double> vcoords, vector<double> vertex_radius, vector<int>* auxiliary_data=0)
{
	int num_crossings = code_data.num_crossings;
	int num_components = code_data.num_components;
//	int head = code_data.head;
	int num_edges = 2*num_crossings;
	int num_type12_vertices = num_edges+num_crossings;	
	int temp; // used for reading integers from triangulation_output_file and circlepack_output_file
	
	vector<int> state;
       
	if (mp_control.state_smoothed && auxiliary_data !=0)
	{
		state = *auxiliary_data;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: state: ";
	for (unsigned int i=0; i< state.size(); i++)
		debug << state[i] << " ";
	debug << endl;
}		
	}

	vector<int>& shortcut_crossing = code_data.shortcut_crossing;
	bool ignore_shortcut = false;
	if (code_data.head != -1 && shortcut_crossing.size())
		ignore_shortcut = true;

	if (mp_control.colourmap.length() != 0)
	{
		ifstream colourfile(mp_control.colourmap.c_str());
				
		if (!colourfile)
		{
			cout << "Error opening " << mp_control.colourmap << endl;
			exit(0);
		}
		else
		{
			mp_control.draw_colour.erase(mp_control.draw_colour.begin(),mp_control.draw_colour.end());
			
			string next_colour;
			while (getline(colourfile,next_colour))
			{
				mp_control.draw_colour.push_back(next_colour);
			}
		}
	}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: mp_control: " << endl;
    print (mp_control, debug, "write_metapost:   ");
    debug << "write_metapost: TRACK_PLACEMENT_ITERATION = " << (TRACK_PLACEMENT_ITERATION? "true": "false") << endl;
}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "write_metapost: num_crossings = " << num_crossings << endl;
    debug << "write_metapost: num_edges = " << num_edges << endl;
    debug << "write_metapost: num_type12_vertices = " << num_type12_vertices << endl;
}

	matrix<int>& code_table	 = code_data.code_table;	
	vector<int>& term_crossing = code_data.term_crossing;
	vector<int>& num_component_edges = code_data.num_component_edges;
	vector<int>& first_edge_on_component = code_data.first_edge_on_component;

	/* record the first edge that enters each crossing */
	vector<int> first_edge(num_crossings,-1);
    for (int i=0; i< num_edges; i++)
    {
		if (first_edge[term_crossing[i]] == -1)
			first_edge[term_crossing[i]] = i;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: first_edge: ";
    for (int i=0; i< num_crossings; i++)
		debug << first_edge[i] << ' ';
	debug << endl;
}

	int semi_arc = -1;
	if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1)
	{
		if (code_table[LABEL][code_data.head] == generic_code_data::POSITIVE)
			semi_arc = code_table[OPEER][code_data.head];
		else if (code_table[LABEL][code_data.head] == generic_code_data::NEGATIVE)
			semi_arc = 2*code_data.head;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: head " << code_data.head << " lies on semi-arc " << semi_arc << endl;
	
	}
	else
	{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: no head semi-arc" << endl;
	}

	
	/* Read the total node count and number of type 1, 2 & 3 vertices and the number of Reidemeister I loops
	   in the diagram from triangulation_output_file  - variables initialized to -1 to catch errors
	*/
	ifstream Tinput;
	int num_vertices = -1;
	int num_type1234_vertices = -1;
	int num_Reidemeister_I_loops = -1;
	
	Tinput.open(triangulation_output_file);
	if (!Tinput)
	{
		cout << "\nError opening output file " << triangulation_output_file << endl;
		exit(0);
    }		

	string next_line;
	getline(Tinput,next_line);
	if (next_line.find("NODECOUNT") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> num_vertices;
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "write_metapost: num_vertices = " << num_vertices << endl;
	}

	getline(Tinput,next_line);
	if (next_line.find("ALPHA") != string::npos)
	{
		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> temp; // alpha, always 1, the first vertex
		iss >> temp; // beta, the first boundary (type 5) vertex
		num_type1234_vertices = temp-1;
		iss >> num_Reidemeister_I_loops; //gamma
			
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "write_metapost: number of type 1,2,3 & 4 vertices = " << num_type1234_vertices << " number of Reidemeister I loops = " << num_Reidemeister_I_loops << endl;
	}

	vector<int> Reidemeister_I_loop_edges(num_Reidemeister_I_loops);
	if (num_Reidemeister_I_loops >0)
	{
		getline(Tinput,next_line);
		
		if (next_line.find("LOOPS") != string::npos)
		{
			char c;
			istringstream iss(next_line);
			do { iss >> c; } while (c != ':');
			
			for (int i=0; i< num_Reidemeister_I_loops; i++)
				iss >> Reidemeister_I_loop_edges[i];
		}
		
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "write_metapost: Reidemeister_I_loop_edges: ";
	for (int i=0; i< num_Reidemeister_I_loops; i++)
		debug << Reidemeister_I_loop_edges[i] << ' ';
	debug << endl;
}
	}
	
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: initial vertex_coordinates: " << endl;
    for (int i=0; i< 2*num_vertices; i++)
    {
		debug << "write_metapost:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << vcoords[i][j] << ' ';
		debug << endl;
	}	
    debug << "write_metapost: initial vertex_radius: " << endl;
    for (int i=0; i< num_vertices; i++)
		debug << "write_metapost:   vertex " << i << ": " << vertex_radius[i] << endl;
}


	/* If we're drawing in the hyperbolic plane, evaluate the Euclidean centre
	   and radius of each circle.  The centre and radius given in circlepack_output_file
	   are the hyperbolic cetre and radius in this case.
	*/
	if (DRAW_IN_HYPERBOLIC_DISC)
		hyperbolic_representation(vcoords, vertex_radius);

	/* scale vertices and radii */
	
	for (int i=0; i< 2*num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			vcoords[i][j] *= metapost_coordinate_scale_factor;
	}
	for (int i=0; i< num_vertices; i++)
		vertex_radius[i] *= metapost_coordinate_scale_factor;

/*
	for (int i=0; i< num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			orig_vcoords[i][j] *= metapost_coordinate_scale_factor;
	}
*/		
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: metapost_coordinate_scale_factor = " << metapost_coordinate_scale_factor << endl;
    debug << "write_metapost: scaled vertex coordinates: " << endl;
    for (int i=0; i< 2*num_vertices; i++)
    {
		debug << "write_metapost:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << vcoords[i][j] << ' ';
		debug << endl;
	}	
    debug << "write_metapost: scaled vertex radius: " << endl;
    for (int i=0; i< num_vertices; i++)
		debug << "write_metapost:   vertex " << i << ": " << vertex_radius[i] << endl;

/*    debug << "write_metapost: scaled original coordinates: " << endl;
    for (int i=0; i< num_vertices; i++)
    {
		debug << "write_metapost:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << orig_vcoords[i][j] << ' ';
		debug << endl;
	}	
*/
}


	average_triangulation_edge_length *= metapost_coordinate_scale_factor;	
	triangulation_plot_threshold *= metapost_coordinate_scale_factor;	

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: scaled average_triangulation_edge_length = " << average_triangulation_edge_length << endl;
    debug << "write_metapost: scaled triangulation_plot_threshold = " << triangulation_plot_threshold << endl;
}

	/* determine the frame extremities */
	double frame_minx = mp_control.frame_minx;
	double frame_maxx = mp_control.frame_maxx;
	double frame_miny = mp_control.frame_miny;
	double frame_maxy = mp_control.frame_maxy;

/*
	for (int i=1; i< num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
		{
			if (j==0 && vcoords[i][j] < frame_minx)
				frame_minx = orig_vcoords[i][j];
			
			if (j==1 && vcoords[i][j] < frame_miny)
				frame_miny = orig_vcoords[i][j];
			
			if (j==0 && vcoords[i][j] > frame_maxx)
				frame_maxx = orig_vcoords[i][j];
			
			if (j==1 && vcoords[i][j] > frame_maxy)
				frame_maxy = orig_vcoords[i][j];
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_metapost: frame_minx = " << frame_minx << ", frame_maxx = " << frame_maxx << ", frame_miny = " << frame_miny << ", frame_maxy = " << frame_maxy << endl;
*/

	/* For the drawing of axes when drawing the vertices we use the type 1 & 2 vertices to give us bounds */
	double minx = vcoords[0][0];
	double maxx = vcoords[0][0];
	double miny = vcoords[0][1];
	double maxy = vcoords[0][1];
	
	for (int i=1; i< num_type12_vertices; i++)
	{
		for (int j=0; j< 2; j++)
		{
			if (j==0 && vcoords[i][j] < minx)
				minx = vcoords[i][j];
			
			if (j==1 && vcoords[i][j] < miny)
				miny = vcoords[i][j];
			
			if (j==0 && vcoords[i][j] > maxx)
				maxx = vcoords[i][j];
			
			if (j==1 && vcoords[i][j] > maxy)
				maxy = vcoords[i][j];
		}
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_metapost: minx = " << minx << ", maxx = " << maxx << ", miny = " << miny << ", maxy = " << maxy << endl;

	/* moved outside if clause to allow drawing of the boundary circles 
	double radians;
	double M00;
	double M10;
	double M01;
	double M11; */

	/* If we are translating any of the vertices, make the required adjustments */
	if (mp_control.translations.size() != 0)
	{
		int x_percentage_unit = (maxx - minx)/100;
		int y_percentage_unit = (maxy - miny)/100;
		
		for (unsigned int i=0; i< mp_control.translations.size(); i++)
		{
			unsigned int vertex = get<0>(mp_control.translations[i]);
			int x_num = get<1>(mp_control.translations[i]);
			int y_num = get<2>(mp_control.translations[i]);
			
			if (vertex < vcoords.numrows())
			{
				vcoords[vertex][0] += x_num*x_percentage_unit;
				vcoords[vertex][1] += y_num*y_percentage_unit;
			}
		}
	}

	if (mp_control.rotate)
	{
		if (mp_control.implicit_rotation_centre)
		{
			mp_control.rotation_centre_x = vcoords[mp_control.rotation_centre_z][0];
			mp_control.rotation_centre_y = vcoords[mp_control.rotation_centre_z][1];
		}
		else if (!mp_control.explicit_rotation_centre)
		{
			mp_control.rotation_centre_x = (minx+maxx)/2;
			mp_control.rotation_centre_y = (miny+maxy)/2;
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: rotating by " << mp_control.rotation_degrees << " around ";
    if (mp_control.explicit_rotation_centre)
		debug << "explicit centre (" << mp_control.rotation_centre_x << "," << mp_control.rotation_centre_y << ")" << endl;
	else if (mp_control.implicit_rotation_centre)
	{
		debug << "implicit centre z" << mp_control.rotation_centre_z << " = (";
		debug << vcoords[mp_control.rotation_centre_z][0] << "," << vcoords[mp_control.rotation_centre_z][1] << ")" << endl;
	}
	else
		debug << "default centre (" << mp_control.rotation_centre_x << "," << mp_control.rotation_centre_y << ")" << endl;
}

	
		rotate_metapost(vcoords, (TRACK_PLACEMENT_ITERATION? 2*num_vertices: num_vertices), mp_control);
		
		/* we convert to radians by calculating rotation_degrees/360 * 2 pi and use 355/113 as an approximation to pi 

		radians = mp_control.rotation_degrees * 355 / 180 / 113;
		M00 = cos(radians);
		M10 = sin(radians);
		M01 = M10 * -1;
		M11 = M00;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_metapost: cos(rotation_degrees = " << radians << ") = " << M00 << ", sin(rotation_degrees) = " << M10 << endl;

		for (int i=0; i< (TRACK_PLACEMENT_ITERATION? 2*num_vertices: num_vertices); i++)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "write_metapost: rotating z" << i << " = (" << vcoords[i][0] << "," << vcoords[i][1] << ")" << endl;
    
			double shift_x = vcoords[i][0] - mp_control.rotation_centre_x;
			double shift_y = vcoords[i][1] - mp_control.rotation_centre_y;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "write_metapost:   shifted to (" << shift_x << "," << shift_y << ")" << endl;

			vcoords[i][0] = M00*shift_x + M01*shift_y;
			vcoords[i][1] = M10*shift_x + M11*shift_y;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "write_metapost:   rotated to (" << vcoords[i][0] << "," << vcoords[i][1] << ")" << endl;
    
			vcoords[i][0] += mp_control.rotation_centre_x;
			vcoords[i][1] += mp_control.rotation_centre_y;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "write_metapost:   shifted back to (" << vcoords[i][0] << "," << vcoords[i][1] << ")" << endl;
			
		}
*/

		/* re-evaluate the min/max values after the rotation */
		minx = vcoords[0][0];
		maxx = vcoords[0][0];
		miny = vcoords[0][1];
		maxy = vcoords[0][1];
		
		for (int i=1; i< num_type12_vertices; i++)
		{
			for (int j=0; j< 2; j++)
			{
				if (j==0 && vcoords[i][j] < minx)
					minx = vcoords[i][j];
				
				if (j==1 && vcoords[i][j] < miny)
					miny = vcoords[i][j];
				
				if (j==0 && vcoords[i][j] > maxx)
					maxx = vcoords[i][j];
				
				if (j==1 && vcoords[i][j] > maxy)
					maxy = vcoords[i][j];
			}
		}

	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_metapost: after rotation minx = " << minx << ", maxx = " << maxx << ", miny = " << miny << ", maxy = " << maxy << endl;

	os << "\n\n";
	
	if (title.length())
		os << "% " << title << endl;

	matrix<int>* Cycle = new matrix<int>(num_crossings+2, num_edges+1);
	int num_left_cycles;
	int num_cycles;
	
	os << "% ";
	write_code_data(os, code_data);
	os << endl;
	calculate_turning_cycles(code_data, *Cycle, num_left_cycles, num_cycles);

	matrix<int>& cycle = *Cycle;
	
    os << "% left turning cycles:";   
	for (int i=0; i<num_left_cycles; i++)
	{
		os << "\n% cycle " << i << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			os << cycle[i][j] << " ";
	}
	
    os << "\n% right turning cycles:";
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		os << "\n% cycle " << i << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			os << cycle[i][j] << " "	;
	}

    os << "\n% infinite turning cycle = " << mp_control.infinite_cycle << endl;
		
	delete Cycle;
	
	if (mp_control.state_smoothed && auxiliary_data !=0)
	{
		os << "% state:  ";
		for (unsigned int i=0; i< state.size(); i++)
			os << state[i] << " ";
		os << endl;
	}
	
	if (mp_control.hamiltonians && auxiliary_data->size() > 1)
	{
		vector<int>& circuit = *auxiliary_data;
		os << "% hamiltonian edge circuit: ";
		for (size_t i=0; i< circuit.size(); i++)
			os << circuit[i] << " ";
		os << endl;
	}
	
	/* scale relative to a fixed size if the scale option has been provided */
	if (mp_control.scale)
	{
		os << "% x-size = " << (maxx - minx)*mp_control.unit_size/100 << ", y-size = " << (maxy-miny)*mp_control.unit_size/100 << endl;
		float x_size = (maxx - minx)*mp_control.unit_size/100;
		float y_size = (maxy-miny)*mp_control.unit_size/100;
		float diag_size = sqrt(x_size*x_size +y_size*y_size);
		float fixed_factor = mp_control.scale*100/diag_size;
		
		for (int i=0; i< 2*num_vertices; i++)
		{
			for (int j=0; j< 2; j++)
				vcoords[i][j] *= fixed_factor;
		}
		for (int i=0; i< num_vertices; i++)
			vertex_radius[i] *= fixed_factor;
	}
	
	/* controls for coordinate output */
	int num_decimal_points = 3;
	int output_field_width = 12;
	output.setf(ios::fixed,ios::floatfield);
	output.precision(num_decimal_points);

	os << "\nbeginfig(fignum);" << endl;
	os << "fignum:=fignum+1;" << endl;
	os << "numeric u,d,theta;" << endl;
	float unit_points = mp_control.unit_size;
	os << "u=" << unit_points/100 << "pt;" << endl;
	os << "d=" << mp_control.disc_size << "u;" << endl;
	os << "path p[];" << endl;
	os << "pickup pencircle scaled " << mp_control.pen_size*0.5 << "pt;" << endl;
	
	/* vertex coordinates are numbered z<vertex> */
	for (int i=0; i< num_vertices; i++)
		os << "z" << i << "=(" << setw(output_field_width) << vcoords[i][0] << "u," << setw(output_field_width) << vcoords[i][1] << "u);" << endl;
	
	if (TRACK_PLACEMENT_ITERATION)
	{
		for (int i=0; i< num_vertices; i++)
			os << "z" << i+num_vertices << "=(" << setw(output_field_width) << vcoords[i+num_vertices][0] << "u," 
			   << setw(output_field_width) << vcoords[i+num_vertices][1] << "u);" << endl;			
	}
	
		
	if (mp_control.draw_frame_corners)
	{
		int bl = 2*num_vertices;
		int br = bl+1;
		int tr = br+1;
		int tl = tr+1;
		
		/* set up variables for the frame corners (minx,miny), (maxx,miny), (maxx,maxy), (minx,maxy) */
		os << "z" << bl << "=(" << setw(output_field_width) << frame_minx << "u," 
			   << setw(output_field_width) << frame_miny << "u);" << endl;	
		os << "z" << br << "=(" << setw(output_field_width) << frame_maxx << "u," 
			   << setw(output_field_width) << frame_miny << "u);" << endl;	
		os << "z" << tr << "=(" << setw(output_field_width) << frame_maxx << "u," 
			   << setw(output_field_width) << frame_maxy << "u);" << endl;	
		os << "z" << tl << "=(" << setw(output_field_width) << frame_minx << "u," 
			   << setw(output_field_width) << frame_maxy << "u);" << endl;	

		os << "pickup pencircle scaled 1bp;" << endl;
		os << "draw z" << bl << "--0.1[z" << bl <<",z" << br << "] dashed withdots;" << endl;
		os << "draw z" << bl << "--0.1[z" << bl <<",z" << tl << "] dashed withdots;" << endl;
		os << "draw z" << br << "--0.1[z" << br <<",z" << bl << "] dashed withdots;" << endl;
		os << "draw z" << br << "--0.1[z" << br <<",z" << tr << "] dashed withdots;" << endl;
		os << "draw z" << tr << "--0.1[z" << tr <<",z" << br << "] dashed withdots;" << endl;
		os << "draw z" << tr << "--0.1[z" << tr <<",z" << tl << "] dashed withdots;" << endl;
		os << "draw z" << tl << "--0.1[z" << tl <<",z" << tr << "] dashed withdots;" << endl;
		os << "draw z" << tl << "--0.1[z" << tl <<",z" << bl << "] dashed withdots;" << endl;
		os << "pickup pencircle scaled 0.5bp;" << endl;		
	}

	if (mp_control.draw_grid)
	{
		double x_step = mp_control.grid_size*(maxx - minx)/100;
		double y_step = mp_control.grid_size*(maxy - miny)/100;
		double right_x = minx+(100/mp_control.grid_size)*x_step;
		double top_y = miny+(100/mp_control.grid_size)*y_step;
		
		os << "for i=0 upto " << 100/mp_control.grid_size << ":" << endl;
		os << "  draw (" << setw(output_field_width) << minx << "*u,(" << setw(output_field_width) << miny << "+i*" << setw(output_field_width) << y_step << ")*u)--(" << setw(output_field_width) << right_x << "*u,(" << setw(output_field_width) << miny << "+i*" << setw(output_field_width) << y_step << ")*u) dashed evenly;" << endl;
		os << "  draw ((" << setw(output_field_width) << minx << "+i*" << setw(output_field_width) << x_step << ")*u," << setw(output_field_width) << miny << "*u)--((" << setw(output_field_width) << minx << "+i*" << setw(output_field_width) << x_step << ")*u," << setw(output_field_width) << top_y << "*u) dashed evenly;" << endl;
		os << "endfor;" << endl;
	}
	
	os << endl;
	

	/* create the vertex sequence for writing the metapost path using the terminating 
	   crossing information in code_data for each edge.  The sequence starts with the 
	   type 2 vertex preceeding crossing zero, which is the first type 2 vertex in the 
	   numbering scheme.	       
	*/
	vector<int> vertex_sequence(2*num_edges);
	for (int i=0; i< num_edges; i++)
	{
//		vertex_sequence[2*i] = num_crossings+(i%num_edges); //the type_2_vertex at the midpoint of edge i
		vertex_sequence[2*i] = num_crossings+i; //changed 22/5/16, since i<num_edges, i%num_edges == i
		vertex_sequence[2*i+1] = term_crossing[i];
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: vertex_sequence: ";
    for (int i=0; i< 2*num_edges + 2*num_Reidemeister_I_loops; i++)
		debug << vertex_sequence[i] << ' ';
	debug << endl;
}

	/* We provide the option to require that the two arcs of the metapost path at a crossing meet at right angles.
	   To achieve this we provide control over the direction of the metapost path at the end of every edge; 
	   that is, as it enters a crosing.  Metapost automatically ensures that the edge leaving the crossing does so
	   in the same direction.
	   
	   The directions are specified as vectors, initialized (by default) to zero, so that if we choose not to require 
	   perpendicular arcs at crossings, the following code produces direction vectors '{0,0}' which are interpreted by 
	   metapost as <empty>.
	*/
	matrix<double> edge_direction(num_edges,2); // (xcoord,ycoord) of direction vector
	
	if ((USE_FORCE_DIRECTED_PLACEMENT || USE_CENTRE_OF_GRAVITY_PLACEMENT) && IMMERSION_CROSSINGS_AT_RIGHT_ANGLES)
	{
		set_edge_directions(code_data,vcoords,edge_direction);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_metapost: setting edge directions: " << endl;
    for (int i=0; i< num_edges; i++)
		debug << "write_metapost: setting edge directions:   vertex " << i << "(" << edge_direction[i][0] << "," << edge_direction[i][1] << ")" << endl;
}		
		
	}
	
	/* since we have T4 vertices on Reidemeister I loop edges, we need to record the offset of each edge in the metapost 
	   path to which it belongs, so we can calculate the subpaths required to draw crossing features.  That is, we record the offset of 
	   the T2 vertex of each edge in the corresponding path.  We also note the number of R1 loops on each component.
	*/
	vector<int> edge_path_offset(num_edges);
	vector<int> num_R1_loops(num_components);
	
	if (mp_control.one_metapost_path)
	{
		/*  This approach can produce inflexion points that are difficult to remove */
		for (int i=0; i< num_components; i++)
		{			
			os << "p" << i+1 << " = ";
			bool first_edge_Reidemeister_I_loop = false;
			int first_edge_first_T4 = -1;
			int cumulative_T4_vertex_count = 0;
			
			for (int j=first_edge_on_component[i]; j< first_edge_on_component[i]+num_component_edges[i]; j++)
			{
				bool Reidemeister_I_loop = false;
				int first_T4 = -1;
				
				vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),j);
				if (vptr != Reidemeister_I_loop_edges.end())
				{
					Reidemeister_I_loop = true;
					
					int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
					first_T4 = num_type1234_vertices - 2*num_Reidemeister_I_loops + 2*loop_edge_index;

					if (j == first_edge_on_component[i])
					{
						first_edge_Reidemeister_I_loop = true;
						first_edge_first_T4 = first_T4;
					}
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "write_metapost:   edge " << j << " is a Reidemeister I loop, type 4 vertices " << first_T4 << ' ' << first_T4+1 << endl;

					os << "z" << first_T4;
					cumulative_T4_vertex_count++;

					if (mp_control.tension)
						os << "..tension..";
					else
						os << "..";
				}
								
				os << "z" << vertex_sequence[2*j];
				edge_path_offset[vertex_sequence[2*j]-num_crossings] = 2*(j-first_edge_on_component[i]) + cumulative_T4_vertex_count;
				
				if (Reidemeister_I_loop)
				{
					if (mp_control.tension)
						os << "..tension..";
					else
						os << "..";

					os << "z" << first_T4+1;
					cumulative_T4_vertex_count++;
				}
				
				if (mp_control.tension)
					os << "..tension " << metapost_path_tension << "..{";
				else
					os << "..{";
				
				os << edge_direction[j][0]<< "," << edge_direction[j][1] << "}z" << vertex_sequence[2*j+1];

				if (mp_control.tension)
					os << "..tension " << metapost_path_tension << "..";
				else
					os << "..";				
			}

			if (first_edge_Reidemeister_I_loop)
				os << "z" << first_edge_first_T4;
			else
				os << "z" << vertex_sequence[2*first_edge_on_component[i]];
				
			if (mp_control.tension)
				os << "..tension " << metapost_path_tension << "..";
			else
				os << "..";
			os << "cycle;" << endl;
			
			num_R1_loops[i] = cumulative_T4_vertex_count/2;
		}		
	}
	else
	{
		for (int i=0; i< num_components; i++)
		{
			os << "p" << 2*i+1 << " = ";
			bool first_edge_Reidemeister_I_loop = false;
			int first_edge_first_T4 = -1;
			int cumulative_T4_vertex_count = 0;

			for (int j=first_edge_on_component[i]; j< first_edge_on_component[i]+num_component_edges[i]; j++)
			{
				bool Reidemeister_I_loop = false;
				int first_T4 = -1;
				
				vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),j);
				if (vptr != Reidemeister_I_loop_edges.end())
				{
					Reidemeister_I_loop = true;
					
					int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
					first_T4 = num_type1234_vertices - 2*num_Reidemeister_I_loops + 2*loop_edge_index;

					if (j == first_edge_on_component[i])
					{
						first_edge_Reidemeister_I_loop = true;
						first_edge_first_T4 = first_T4;
					}
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "write_metapost:   edge " << j << " is a Reidemeister I loop, type 4 vertices " << first_T4 << ' ' << first_T4+1 << endl;

					os << "z" << first_T4;
					cumulative_T4_vertex_count++;

					if (mp_control.tension)
						os << "..tension..";
					else
						os << "..";
				}
				
				os << "z" << vertex_sequence[2*j];
				edge_path_offset[vertex_sequence[2*j]-num_crossings] = 2*(j-first_edge_on_component[i]) + cumulative_T4_vertex_count;
				
				if (Reidemeister_I_loop)
				{
					if (mp_control.tension)
						os << "..tension..";
					else
						os << "..";

					os << "z" << first_T4+1;
					cumulative_T4_vertex_count++;
				}
				
				if (mp_control.tension)
					os << "..tension " << metapost_path_tension << "..{";
				else
					os << "..{";
				
				os << edge_direction[j][0]<< "," << edge_direction[j][1] << "}z" << vertex_sequence[2*j+1];

				if (mp_control.tension)
					os << "..tension " << metapost_path_tension << "..";
				else
					os << "..";				
			}
								
			if (first_edge_Reidemeister_I_loop)
				os << "z" << first_edge_first_T4 << ";" << endl;
			else		
				os << "z" << vertex_sequence[2*first_edge_on_component[i]] << ";" << endl;

			num_R1_loops[i] = cumulative_T4_vertex_count/2;
		}
		
		for (int i=0; i< num_components; i++)
		{
			os << "p" << 2*i+2 << " = ";
			bool first_edge_Reidemeister_I_loop = false;
			int first_edge_first_T4 = -1;
			for (int j=first_edge_on_component[i]; j< first_edge_on_component[i]+num_component_edges[i]; j++)
			{
				bool Reidemeister_I_loop = false;
				int first_T4 = -1;
				
				vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),j);
				if (vptr != Reidemeister_I_loop_edges.end())
				{
					Reidemeister_I_loop = true;
					
					int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
					first_T4 = num_type1234_vertices - 2*num_Reidemeister_I_loops + 2*loop_edge_index;

					if (j == first_edge_on_component[i])
					{
						first_edge_Reidemeister_I_loop = true;
						first_edge_first_T4 = first_T4;
					}
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "write_metapost:   edge " << j << " is a Reidemeister I loop, type 4 vertices " << first_T4 << ' ' << first_T4+1 << endl;

					os << "z" << first_T4;

					if (mp_control.tension)
						os << "..tension..";
					else
						os << "..";
				}

				os << "z" << vertex_sequence[2*j];

				if (Reidemeister_I_loop)
				{
					if (mp_control.tension)
						os << "..tension..";
					else
						os << "..";

					os << "z" << first_T4+1;
				}
				
				if (mp_control.tension)
					os << "..tension " << metapost_path_tension << "..{";
				else
					os << "..{";
				
				os << edge_direction[j][0]<< "," << edge_direction[j][1] << "}z" << vertex_sequence[2*j+1];

				if (mp_control.tension)
					os << "..tension " << metapost_path_tension << "..";
				else
					os << "..";				
			}
			
			if (first_edge_Reidemeister_I_loop)
				os << "{direction 0 of p" << 2*i+1 << "}z" << first_edge_first_T4 << ";" << endl;
			else
				os << "{direction 0 of p" << 2*i+1 << "}z" << vertex_sequence[2*first_edge_on_component[i]] << ";" << endl;
		}
	}
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "write_metapost: edge_path_offset: ";
	for (int i=0; i< num_edges; i++)
		debug << edge_path_offset[i] << ' ';
	debug << endl;
	debug << "write_metapost: num_R1_loops: ";
	for (int i=0; i< num_components; i++)
		debug << num_R1_loops[i] << ' ';
	debug << endl;
}
	
	
	if (mp_control.draw_triangulation)
	{
		/* draw triangulation */
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost: drawing triangulation" << endl;

		os << "pickup pencircle scaled 0.25bp;" << endl;		

		if (TRACK_PLACEMENT_ITERATION && mp_control.draw_triangulation_displacement)
		{
			for (int i=0; i < num_type1234_vertices; i++)
				os << "draw z" << i << "--z" << i + num_vertices << " dashed evenly withcolor red;" << endl;

			if (INCLUDE_BOUNDARY_VERTICES)
			{			
				for (int i=num_type1234_vertices; i < num_vertices; i++)
					os << "draw z" << i << "--z" << i + num_vertices << " dashed evenly withcolor red;" << endl;
			}
		}
		
		if (CHECK_INNER_HULL_CALCULATION)
		{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: check inner hull calculation for vertex " << check_inner_hull_vertex << endl;

			ifstream IHinput ("__inner_hull.out");
			if (!IHinput)
			{
				cout << "\nError opening input file __inner_hull.out\n";
				exit(0);
		    }
			
			double boundary_cog_x=0;
			double boundary_cog_y=0;		
			double boundary_vertex_circle_radius = 0;

			int vertex;
			istringstream iss;


			getline(IHinput, next_line); //	IHinput >> next_line;
			
			{
				/* With gcc 4.8.3 I could not simply say
					iss.str(next_line);
				   then read in from the iss, since doing so caused the do clause below to fail.
				   
				   By creating a new, local istringstream within the current scope, the do clause
				   works fine.
				*/
				istringstream iss(next_line);
				iss >> boundary_cog_x >> boundary_cog_y >> boundary_vertex_circle_radius;
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: boundary_cog_x = " << boundary_cog_x << ", boundary_cog_y = " << boundary_cog_y << ", boundary_vertex_circle_radius = " << boundary_vertex_circle_radius << endl;
			
			boundary_cog_x *= metapost_coordinate_scale_factor;
			boundary_cog_y *= metapost_coordinate_scale_factor;
			boundary_vertex_circle_radius *=metapost_coordinate_scale_factor;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: scaled boundary_cog_x = " << boundary_cog_x << ", boundary_cog_y = " << boundary_cog_y << ", boundary_vertex_circle_radius = " << boundary_vertex_circle_radius << endl;
			
			do
			{		
				getline(IHinput, next_line); //	IHinput >> next_line;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: next __inner_hull.out line: " << next_line << endl;

				iss.str(next_line);
				iss >> vertex;				
			} while (vertex != check_inner_hull_vertex);		
			
			int num_IH_vertices;
			iss >> num_IH_vertices;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: number of inner hull vertices = " << num_IH_vertices << endl;

			matrix<double> IH_vertex(num_IH_vertices,2);
			
			for (int i=0; i< num_IH_vertices; i++)
			{
				iss >> IH_vertex[i][0];
				iss >> IH_vertex[i][1];
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "write_metapost: inner hull vertices: " << endl;
	for (int i=0; i< num_IH_vertices; i++)
		debug << "write_metapost: (" << IH_vertex[i][0] << "," << IH_vertex[i][1] << ")" << endl;
		
	debug <<"\nmetapost_coordinate_scale_factor = " << metapost_coordinate_scale_factor << endl;
}

			for (int i=0; i< num_IH_vertices; i++)
			{
				IH_vertex[i][0] *= metapost_coordinate_scale_factor;
				IH_vertex[i][1] *= metapost_coordinate_scale_factor;
			}


if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "write_metapost: scaled inner hull vertices: " << endl;
	for (int i=0; i< num_IH_vertices; i++)
		debug << "write_metapost: (" << IH_vertex[i][0] << "," << IH_vertex[i][1] << ")" << endl;
}

			rotate_metapost(IH_vertex, num_IH_vertices,mp_control);
			
			os << "\npath IHP[];" << endl;
			os << "IHP1=(" << IH_vertex[1][0] << "u," << IH_vertex[1][1] << "u)";
			for (int i=2; i< num_IH_vertices; i++)
			{
				os << " -- (" << IH_vertex[i][0] << "u," << IH_vertex[i][1] << "u)";
			}
			os <<";" << endl;
			os << "IHP2=(" << IH_vertex[0][0] << "u," << IH_vertex[0][1] << "u) -- IHP1 -- cycle;" << endl;
			os << "fill IHP2 withcolor (.5,.5,.5);" << endl;
			os << "draw IHP2;" << endl;			
			
			os << "pickup pencircle scaled 1bp;" << endl;
			os << "draw fullcircle scaled " << boundary_vertex_circle_radius << "/15d shifted (" << boundary_cog_x << "u," << boundary_cog_y << "u) dashed withdots;"  << endl;
			os << "pickup pencircle scaled 0.25bp;" << endl;		
			
			IHinput.close();
		}

			
		string next_line;
		do
		{		
			Tinput >> next_line;
		} while (next_line.find("FLOWERS") == string::npos);
				
		matrix<int> edge (num_vertices,num_vertices); // records whether we've drawn an edge between neighbours
		for (int i=0; i< num_vertices; i++)
		{
			if (!INCLUDE_BOUNDARY_VERTICES && i >= num_type1234_vertices)			
			    break;
			    
			int vertex;
			Tinput >> vertex;
			vertex--;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost:   vertex "<< vertex << endl;
			
			int num_neighbours;
			Tinput >> num_neighbours;	

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost:   num_neighbours "<< num_neighbours << endl;
			
			/* number of neighbours in triangulation_output_file does not include the repeated node
			   for type 1, 2 and 3 vertices and is one less than the number of neighbours 
			   for type 4 vertices
			*/
//			for (int j=0; j<= num_neighbours; j++)
			for (int j=0; j<= num_neighbours; j++)
			{					
				int neighbour;
				Tinput >> neighbour;

//				if (vertex < num_type1234_vertices && j == num_neighbours)
				if (j == num_neighbours)
					break; // don't process the repeated petal

				neighbour--;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost:     neighbour "<< neighbour;
				
				if ((INCLUDE_BOUNDARY_VERTICES || neighbour < num_type1234_vertices) && edge[vertex][neighbour] == 0)
				{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " edge not drawn"<< endl;
	
					if (TRACK_PLACEMENT_ITERATION && mp_control.draw_triangulation_displacement)
						os << "draw z" << vertex+num_vertices << "--z" << neighbour+num_vertices << " dashed evenly withcolor blue;" << endl;


					/* evaluate the original length of the edge for highlighting when magnifying */
					double delta_x = vcoords[num_vertices+vertex][0] - vcoords[num_vertices+neighbour][0];
					double delta_y = vcoords[num_vertices+vertex][1] - vcoords[num_vertices+neighbour][1];
					double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);


					if (mp_control.magnify_small_circles && mod_delta < average_triangulation_length_factor * average_triangulation_edge_length)
						os << "draw z" << vertex << "--z" << neighbour << " withcolor red;" << endl;
					else if (mp_control.highlight_small_edges && mod_delta < triangulation_plot_threshold)
						os << "draw z" << vertex << "--z" << neighbour << " withcolor red;" << endl;
					else if (mp_control.draw_immersion)
					{
/*31-3-18 test code 
if(false && neighbour <= num_type12_vertices)
	os << "draw z" << vertex << "--z" << neighbour << " dashed evenly withcolor red;" << endl;
else
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " here drawing " << vertex << " to " << neighbour << endl;
*/
						os << "draw z" << vertex << "--z" << neighbour << " dashed evenly;" << endl;
					}
					else
						os << "draw z" << vertex << "--z" << neighbour << ";" << endl;
					
						
					edge[vertex][neighbour] = 1;
					edge[neighbour][vertex] = 1;
				}
				else if (INCLUDE_BOUNDARY_VERTICES || neighbour < num_type1234_vertices)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " already drawn"<< endl;
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " is a type 4 vertex"<< endl;
				}
				
			}
		}
//setw(output_field_width) << vcoords[i][0] << "u,

		os << "pickup pencircle scaled 0.5bp;" << endl;		
	}
		
	/* if we're using region shrinking placement we may want to draw the effect of shrinking the
	   triangulation in a region of maximal area
	*/
	if (USE_REGION_SHRINKING_PLACEMENT && mp_control.draw_shrink_effect)
	{
		ifstream Sinput;
	
		Sinput.open(triangulation_shrink_file);
		if (!Sinput)
		{
			cout << "\nError opening shrink file " << triangulation_shrink_file << endl;
			exit(0);
	    }		

		int num_triangles;
		Sinput >> num_triangles;
		matrix<double> shrink_coords(3*num_triangles,2);
		
		for (int i = 0; i < 3*num_triangles; i++)
		for (int j=0; j < 2; j++)
			Sinput >> shrink_coords[i][j];

			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "write_metapost: trianglulation shrink vertices before scaling:" << endl;
	for (int i=0; i< num_triangles; i++)
	{
		debug << "write_metapost: triangle << " << i << ": (";
		debug << shrink_coords[3*i][0] << "," << shrink_coords[3*i][1] << "), (";
		debug << shrink_coords[3*i+1][0] << "," << shrink_coords[3*i+1][1] << "), (";
		debug << shrink_coords[3*i+2][0] << "," << shrink_coords[3*i+2][1] << ")" << endl;
	}
}

		for (int i = 0; i < 3*num_triangles; i++)
		for (int j=0; j < 2; j++)
			shrink_coords[i][j] *= metapost_coordinate_scale_factor;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "write_metapost: trianglulation shrink vertices after scaling:" << endl;
	for (int i=0; i< num_triangles; i++)
	{
		debug << "write_metapost: triangle << " << i << ": (";
		debug << shrink_coords[3*i][0] << "," << shrink_coords[3*i][1] << "), (";
		debug << shrink_coords[3*i+1][0] << "," << shrink_coords[3*i+1][1] << "), (";
		debug << shrink_coords[3*i+2][0] << "," << shrink_coords[3*i+2][1] << ")" << endl;
	}
}

		
		rotate_metapost (shrink_coords,3*num_triangles,mp_control);

		/* draw the shrunken triangles in red */
		for (int i = 0; i < num_triangles; i++)
		{
			os << "draw (" << setw(output_field_width) << shrink_coords[3*i][0] << "u," << setw(output_field_width) << shrink_coords[3*i][1] << "u)--(" << setw(output_field_width) << shrink_coords[3*i+1][0] << "u," << setw(output_field_width) << shrink_coords[3*i+1][1] << "u) withcolor red;" << endl;
			os << "draw (" << setw(output_field_width) << shrink_coords[3*i+1][0] << "u," << setw(output_field_width) << shrink_coords[3*i+1][1] << "u)--(" << setw(output_field_width) << shrink_coords[3*i+2][0] << "u," << setw(output_field_width) << shrink_coords[3*i+2][1] << "u) withcolor red;" << endl;
			os << "draw (" << setw(output_field_width) << shrink_coords[3*i+2][0] << "u," << setw(output_field_width) << shrink_coords[3*i+2][1] << "u)--(" << setw(output_field_width) << shrink_coords[3*i][0] << "u," << setw(output_field_width) << shrink_coords[3*i][1] << "u) withcolor red;" << endl;
		}
	}
			
	/* draw the immersion */
	if (mp_control.draw_immersion)
	{
		os << "ahlength := " << mp_control.arrowhead_bp_size << "bp;" << endl;

		if (code_data.immersion != generic_code_data::character::KNOTOID && code_data.immersion != generic_code_data::character::PURE_KNOTOID && code_data.immersion != generic_code_data::character::MULTI_LINKOID)
		{
			for (int i=0; i< num_components; i++)
			{
				int path_number = (mp_control.one_metapost_path? i+1: 2*i+2);

				if (mp_control.draw_oriented)
					os << "drawarrow subpath(0,0.5) of p" << path_number << ";" << endl;
				
				os << "draw p" << path_number;
				if (mp_control.colour && mp_control.singlecolour.length() != 0)
					os << " withcolor " << mp_control.singlecolour;
				else if (mp_control.colour && i < static_cast<int>(mp_control.draw_colour.size()))		
					os << " withcolor " << mp_control.draw_colour[i];
				os << ';' << endl;
			}

			if (mp_control.hamiltonians)
			{
				/* auxiliary data is the list of edges corresponding to a Hamiltonian circuit, which we hightlight in green */
				vector<int>& circuit = *auxiliary_data;
				
				for (size_t i=0; i< circuit.size(); i++)
				{
					int edge = circuit[i];
					int component = 0;
					for (int j=1; j< num_components; j++)
					{
						if (edge >= first_edge_on_component[j])
							component++;
						else
							break;
					}
					
					int path_number = (mp_control.one_metapost_path? component+1: 2*component+2);
					
					if (edge_path_offset[edge] !=0)
					{
						os << "draw subpath(" << edge_path_offset[edge]-1 << ',' << edge_path_offset[edge]+1 << ") of p" << path_number << " withcolor "<< mp_control.hamiltonian_colour <<" ;" << endl;						
					}
					else
					{						
						os << "draw subpath(0,1) of p" << path_number << " withcolor "<< mp_control.hamiltonian_colour << ";" << endl;						
						os << "draw subpath(" << 2*num_component_edges[component]-1 << ',' << 2*num_component_edges[component] << ") of p" << path_number << " withcolor "<< mp_control.hamiltonian_colour << ";" << endl;						
					}
				}
			}								
		}
		else if (code_data.immersion == generic_code_data::character::KNOTOID || code_data.immersion == generic_code_data::character::PURE_KNOTOID)
		{
				
//			int component_end_point = (code_data.immersion == generic_code_data::character::KNOTOID? 2*num_component_edges[0]: 2*semi_arc-1);
			int component_end_point = (code_data.immersion == generic_code_data::character::KNOTOID? 2*(num_component_edges[0]+num_R1_loops[0]): 2*(semi_arc+num_R1_loops[0])-1);
			int path_number = (mp_control.one_metapost_path? 1: 2);

			os << "drawarrow subpath(0.5,0.75) of p" << path_number;  
			if (mp_control.colour && mp_control.singlecolour.length() != 0)
				os << " withcolor " << mp_control.singlecolour;
			else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
				os << " withcolor " << mp_control.draw_colour[0];
			os << ';' << endl;
			os << "draw subpath(0.75," << component_end_point << ".5) of p" << path_number;
			if (mp_control.colour && mp_control.singlecolour.length() != 0)
				os << " withcolor " << mp_control.singlecolour;
			else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
				os << " withcolor " << mp_control.draw_colour[0];
			os << ';' << endl;
			
			if (code_data.immersion == generic_code_data::character::KNOTOID)
			{
				os << "draw subpath(0,0.25) of p" << path_number;				
				if (mp_control.colour && mp_control.singlecolour.length() != 0)
					os << " withcolor " << mp_control.singlecolour;
				else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
					os << " withcolor " << mp_control.draw_colour[0];
				os << ';' << endl;
			}
			
			if (mp_control.draw_shortcut)
			{
				if (code_data.immersion == generic_code_data::character::KNOTOID)
				{
					os << "draw subpath(0.25,0.5) of p" << path_number << " dashed " << (mp_control.dash_with_dots? "withdots" : "evenly");
					
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;				
				}
				else
				{
					os << "draw subpath(0,0.5) of p" << path_number << " dashed " << (mp_control.dash_with_dots? "withdots" : "evenly");
					
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;
					
//					os << "draw subpath(" << 2*semi_arc-1 << ".5," << 2*num_component_edges[0]+1 << ") of p" << path_number << " dashed " 
					os << "draw subpath(" << 2*(semi_arc+num_R1_loops[0])-1 << ".5," << 2*(num_component_edges[0]+num_R1_loops[0])+1 << ") of p" << path_number << " dashed " 
					   << (mp_control.dash_with_dots? "withdots" : "evenly");
					   
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;					   
				}
			}
	
			for (int i=1; i< num_components; i++)
			{
				path_number = (mp_control.one_metapost_path? i+1: 2*i+2);

				if (mp_control.draw_oriented)
				{
					os << "drawarrow subpath(0,0.5) of p" << path_number;
						
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && i < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[i];
					os << ';' << endl;						
				}
				
				os << "draw p" << path_number;
				
				if (mp_control.colour && mp_control.singlecolour.length() != 0)
					os << " withcolor " << mp_control.singlecolour;
				else if (mp_control.colour && i < static_cast<int>(mp_control.draw_colour.size()))		
					os << " withcolor " << mp_control.draw_colour[i];
				os << ';' << endl;					
			}
		}
		else if (code_data.immersion == generic_code_data::character::MULTI_LINKOID)
		{
			
			for (int i=0; i< num_components; i++)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost: immersion for component " << i << endl;
	
				int path_number = (mp_control.one_metapost_path? i+1: 2*i+2);
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost:   path_number " << path_number << endl;

				if (code_data.component_type[i].type == component_character::CLOSED)
				{
					if (mp_control.draw_oriented)
					{
						os << "drawarrow subpath(0,0.5) of p" << path_number;						
						if (mp_control.colour && mp_control.singlecolour.length() != 0)
							os << " withcolor " << mp_control.singlecolour;
						else if (mp_control.colour && i < static_cast<int>(mp_control.draw_colour.size()))		
							os << " withcolor " << mp_control.draw_colour[i];
						os << ';' << endl;						
					}
					
					os << "draw p" << path_number;
					
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && i < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[i];
					os << ';' << endl;					
				}
				else if (code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG || code_data.component_type[i].type == component_character::KNOT_TYPE_END_LEG)
				{
//					int leg_point = (code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG? 0: 2*(code_data.num_component_edges[i]-1));
					int leg_point = (code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG? 0: 2*(code_data.num_component_edges[i]+num_R1_loops[i]-1));
					int main_path_start_point = (code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG? 2: 0);
//					int component_end_point = (code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG? 2*num_component_edges[i]: 2*num_component_edges[i]-2);
					int component_end_point = (code_data.component_type[i].type == component_character::KNOT_TYPE_START_LEG? 2*(num_component_edges[i]+num_R1_loops[i]): 2*(num_component_edges[i]+num_R1_loops[i])-2);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost:   leg_point = " << leg_point << ", main_path_start_point = " << main_path_start_point << ", component_end_point = " << component_end_point << endl;

					os << "drawarrow subpath(" << leg_point << ".5," << leg_point << ".75) of p" << path_number;
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;

					/* this is the end of the path in the case KNOT_TYPE_END_LEG */
					os << "draw subpath(" << leg_point << ".75," << leg_point+2 << ") of p" << path_number;
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;						

//					os << "draw subpath(" << main_path_start_point << ',' << component_end_point << ".5) of p" << path_number;
					os << "draw subpath(" << main_path_start_point << ',' << component_end_point << ") of p" << path_number;
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;

					os << "draw subpath(" << leg_point << ',' << leg_point << ".25) of p" << path_number;					
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;

					if (mp_control.draw_shortcut)
					{
						os << "draw subpath(" << leg_point << ".25," << leg_point << ".5) of p" << path_number << " dashed " << (mp_control.dash_with_dots? "withdots" : "evenly");						
						if (mp_control.colour && mp_control.singlecolour.length() != 0)
							os << " withcolor " << mp_control.singlecolour;
						else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
							os << " withcolor " << mp_control.draw_colour[0];
						os << ';' << endl;
					}

				}
				else if (code_data.component_type[i].type == component_character::PURE_START_LEG || code_data.component_type[i].type == component_character::PURE_END_LEG)
				{
					int component_head = code_data.component_type[i].head;
					if (code_table[LABEL][component_head] == generic_code_data::POSITIVE)
						semi_arc = code_table[OPEER][component_head];
					else if (code_table[LABEL][component_head] == generic_code_data::NEGATIVE)
						semi_arc = 2*component_head;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost:   head, " << component_head << ", lies on semi-arc " << semi_arc << endl;
	
//					int leg_point = (code_data.component_type[i].type == component_character::PURE_START_LEG? 0: 2*(code_data.num_component_edges[i]-1));
					int leg_point = (code_data.component_type[i].type == component_character::PURE_START_LEG? 0: 2*(code_data.num_component_edges[i]+num_R1_loops[i]-1));
					int main_path_start_point = (code_data.component_type[i].type == component_character::PURE_START_LEG? 2: 0);
//					int head_point = 2*(semi_arc - code_data.first_edge_on_component[i])-1;
					int head_point = 2*(semi_arc +num_R1_loops[i]- code_data.first_edge_on_component[i])-1;
//					int shortcut_end_point = (code_data.component_type[i].type == component_character::PURE_START_LEG? 2*num_component_edges[i]+1: 2*num_component_edges[i]-1);
					int shortcut_end_point = (code_data.component_type[i].type == component_character::PURE_START_LEG? 2*(num_component_edges[i]+num_R1_loops[i])+1: 2*(num_component_edges[i]+num_R1_loops[i])-1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "write_metapost:   leg_point = " << leg_point << ", main_path_start_point = " << main_path_start_point << ", head_point = " << head_point << ", shortcut_end_point = " << shortcut_end_point << endl;
					
					os << "drawarrow subpath(" << leg_point << ',' << leg_point << ".75) of p" << path_number;
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;

					os << "draw subpath(" << leg_point << ".75," << leg_point+2 << ") of p" << path_number;
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;						

					os << "draw subpath(" << main_path_start_point << ',' << head_point << ".5) of p" << path_number;
					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[0];
					os << ';' << endl;
					
					if (mp_control.draw_shortcut)
					{
						os << "draw subpath(" << leg_point << ',' << leg_point << ".5) of p" << path_number << " dashed " << (mp_control.dash_with_dots? "withdots" : "evenly");						
						if (mp_control.colour && mp_control.singlecolour.length() != 0)
							os << " withcolor " << mp_control.singlecolour;
						else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
							os << " withcolor " << mp_control.draw_colour[0];
						os << ';' << endl;
						
						os << "draw subpath(" << head_point << ".5," << shortcut_end_point << ") of p" << path_number << " dashed " << (mp_control.dash_with_dots? "withdots" : "evenly");						   
						if (mp_control.colour && mp_control.singlecolour.length() != 0)
							os << " withcolor " << mp_control.singlecolour;
						else if (mp_control.colour && 0 < static_cast<int>(mp_control.draw_colour.size()))		
							os << " withcolor " << mp_control.draw_colour[0];
						os << ';' << endl;					   
					}		
				}
			}
		}
	}
	
	/* draw crossing features */

	if (mp_control.draw_crossing_features)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost: drawing crossing features" << endl;

		int state_place=0;
		for (int i=0; i< num_crossings; i++)
		{
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost:   crossing " << i;
			
			if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " is virtual" << endl;
		
				os << "draw fullcircle scaled d shifted z" << vertex_sequence[2*first_edge[i]+1] << ";" << endl;
			}
			else if (mp_control.show_odd_parity_crossings && code_table[LABEL][i] == generic_code_data::ODD)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " is odd" << endl;
		
//				os << "fill fullcircle scaled 1.2d shifted z" << vertex_sequence[2*first_edge[i]+1] << ";" << endl;
				os << "fill fullcircle scaled " << mp_control.odd_parity_disc_size << "*0.1d shifted z" << vertex_sequence[2*first_edge[i]+1] << ";" << endl;
			}
			else if (code_table[LABEL][i] != generic_code_data::FLAT)
			{
				bool A_crossing;
				bool seifert_smoothed;
				
				if ( (mp_control.seifert_circles || (mp_control.state_smoothed && auxiliary_data != 0)) && !(ignore_shortcut && shortcut_crossing[i]))
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " is smoothed" << endl;

					/* identify the components and corresponding paths of the even terminating and odd terminating edges */
					int et_edge = code_table[EVEN_TERMINATING][i];
					int ot_edge = code_table[ODD_TERMINATING][i];
					int et_component = code_table[COMPONENT][(et_edge%2? (et_edge-1)/2: et_edge/2)];
					int ot_component = code_table[COMPONENT][(ot_edge%2? (ot_edge-1)/2: ot_edge/2)];					
					int et_path = 2 * (et_component+1);
					int ot_path = 2 * (ot_component+1);
					int et_point = edge_path_offset[code_table[EVEN_TERMINATING][i]]+1;
					int ot_point = edge_path_offset[code_table[ODD_TERMINATING][i]]+1;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "write_metapost:   et_edge = " << et_edge << " ot_edge = " << ot_edge << " et_component = " << et_component << " ot_component = " << ot_component << endl;
	debug << "write_metapost:   et_path = " << et_path << " ot_path = " << ot_path << " et_point = " << et_point << " ot_point = " << ot_point << endl;
}					
					if (code_table[EVEN_TERMINATING][i] == code_table[EVEN_ORIGINATING][i])
						et_point++; // move et_point past the second type 4 vertex
					else if (code_table[ODD_TERMINATING][i] == code_table[ODD_ORIGINATING][i])
						ot_point++; // move ot_point past the second type 4 vertex

					if (mp_control.uniform_smoothed_discs)
					{
						os << "p" << 2*num_components+i+1 << " = fullcircle scaled " << mp_control.smoothed_state_disc_size << "d shifted z" << i << ";" << endl;
					}
					else
					{
						os << "r:= min(arclength (z" << i << "--(point " << et_point-1 << ".5 of p" << et_path << ")),";
						os <<         "arclength (z" << i << "--(point " << ot_point-1 << ".5 of p" << ot_path << ")),";
						os <<         "arclength (z" << i << "--(point " << et_point << ".5 of p" << et_path << ")),";
						os <<         "arclength (z" << i << "--(point " << ot_point << ".5 of p" << ot_path << ")));" << endl;

						os << "p" << 2*num_components+i+1 << " = fullcircle scaled min(2r," << mp_control.smoothed_state_disc_size << "d) shifted z" << i << ";" << endl;
						os << "theta := angle((direction " << et_point << " of p" << et_path << ")+(direction " << ot_point << " of p" << ot_path << "))-90;" << endl;
					}
					
					
					os << "fill p" << 2*num_components+i+1 << " withcolor 1white;" << endl;
					
					if (   (code_table[TYPE][i] == generic_code_data::TYPE1 && code_table[LABEL][i] == generic_code_data::NEGATIVE)
						|| (code_table[TYPE][i] == generic_code_data::TYPE2 && code_table[LABEL][i] == generic_code_data::POSITIVE)
					   )
					{
						/* positive crossing */
						if (state[state_place] == 1 || state[state_place] == 2)
						{
							A_crossing = true;
							seifert_smoothed = true;
						}
						else if (state[state_place] == -1 || state[state_place] == 3)
						{
							A_crossing = false;
							seifert_smoothed = false;
						}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost:   positive crossing " << i << " state_place = " << state_place << " A_crossing = " << A_crossing << " seifert_smoothed = " << seifert_smoothed << endl;
						
					}
					else if (   (code_table[TYPE][i] == generic_code_data::TYPE1 && code_table[LABEL][i] == generic_code_data::POSITIVE)
							  || (code_table[TYPE][i] == generic_code_data::TYPE2 && code_table[LABEL][i] == generic_code_data::NEGATIVE)
							)
					{
						/* negative crossing */
						if (state[state_place] == 1 || state[state_place] == 3)
						{
							A_crossing = false;
							seifert_smoothed = true;
						}
						else if (state[state_place] == -1 || state[state_place] == 2)
						{
							A_crossing = true;
							seifert_smoothed = false;
						}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "write_metapost:   negative crossing " << i << " state_place = " << state_place << " A_crossing = " << A_crossing << " seifert_smoothed = " << seifert_smoothed << endl;
					}
					
					if (!mp_control.seifert_circles)
					{
						os << "label (btex "; 
						if (mp_control.scriptscript_labels)
							os << "\\fiverm ";
						else if (mp_control.script_labels)
							os << "\\sevenrm ";
					    os << (A_crossing?"A":"B") << " etex, z" << i << "+if(r< " << mp_control.smoothed_disc_threshold << "u):("
					       << mp_control.label_shift << "u*cosd theta, " << mp_control.label_shift << "u*sind theta) else:(0,0) fi);" << endl;				    
					}
											
					os << "z" << 2*num_vertices+4*i << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << et_point-1<< "," << et_point << ") of p" << et_path << ";" << endl;
					os << "z" << 2*num_vertices+4*i+1 << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << ot_point-1 << "," << ot_point <<  ") of p" << ot_path << ";" << endl;
					os << "z" << 2*num_vertices+4*i+2 << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << et_point << "," << et_point+1 << ") of p" << et_path << ";" << endl;
					os << "z" << 2*num_vertices+4*i+3 << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << ot_point << "," << ot_point+1 << ") of p" << ot_path << ";" << endl;
					
//					if (state[state_place] == 1) // Seifert smoothed
					if (seifert_smoothed) // Seifert smoothed
					{
						os << "draw z" << 2*num_vertices+4*i << "--z" << 2*num_vertices+4*i+3 << ";" << endl;
						os << "draw z" << 2*num_vertices+4*i+1 << "--z" << 2*num_vertices+4*i+2 << ";" << endl;
					}
					else
					{
						os << "draw z" << 2*num_vertices+4*i << "--z" << 2*num_vertices+4*i+1 << ";" << endl;
						os << "draw z" << 2*num_vertices+4*i+2 << "--z" << 2*num_vertices+4*i+3 << ";" << endl;						
						
						os << "fill fullcircle scaled (" << mp_control.cusp_disc_size << "*0.1d) shifted 0.5[z" << 2*num_vertices+4*i << ",z" << 2*num_vertices+4*i+1 << "];" << endl;
						os << "fill fullcircle scaled (" << mp_control.cusp_disc_size << "*0.1d) shifted 0.5[z"  << 2*num_vertices+4*i+2 << ",z" << 2*num_vertices+4*i+3 << "];" << endl;
					}
					state_place++;
				}
				else
				{
					os << "fill fullcircle scaled d shifted z" << vertex_sequence[2*first_edge[i]+1] << " withcolor 1white;" << endl;
				
					/* draw the over-arc back in, if the label is positive we want the naming edge otherwise the other edge.
					   If we're at crossing zero and it's positive (i.e. we calculate edge to be zero) and if this is a Reidemeister 
					   I loop edge, then the over-arc has to be drawn back in as two subpaths that wrap from the end of the path to 
					   the beginning
					*/
					int edge;
					if (code_table[LABEL][i] == generic_code_data::POSITIVE)	
					{
		
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " has positive label, ";
				
						edge = 2*i;
					}
					else
					{
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " has negative label, ";
				
						edge = code_table[OPEER][i];
					}
						
					int component;
					if ( edge%2 )
						component = code_table[COMPONENT][(edge-1)/2];
					else
						component = code_table[COMPONENT][edge/2];
					
					int over_arc_start = edge_path_offset[edge];

					if ((edge%2 == 0 && code_table[EVEN_ORIGINATING][i] == edge) || (edge%2 == 1 && code_table[ODD_ORIGINATING][i] == edge))
						over_arc_start++; // move past the second type 4 vertex
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "over-arc starts on edge " << edge << " on component " << component << " edge_path_offset = " << over_arc_start << endl;
	
					/* we are at the midpoint vertex of edge on path p at time 2*(edge-first_edge_on_component[component])
					   and at the terminating crossing at time 2*(edge-first_edge_on_component[component])+1
					*/
					os << "draw subpath(" << over_arc_start << ".5,";
					if (mp_control.one_metapost_path)
						os << over_arc_start+1 << ".5) of p" << component+1;
					else
						os << over_arc_start+1 << ".5) of p" << 2*component+2;

					if (mp_control.colour && mp_control.singlecolour.length() != 0)
						os << " withcolor " << mp_control.singlecolour;
					else if (mp_control.colour && component < static_cast<int>(mp_control.draw_colour.size()))		
						os << " withcolor " << mp_control.draw_colour[component];
					os << ';' << endl;
						
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << " is flat " << endl;

				if (mp_control.seifert_circles)
				{
					/* identify the components and corresponding paths of the even terminating and odd terminating edges */
					int et_edge = code_table[EVEN_TERMINATING][i];
					int ot_edge = code_table[ODD_TERMINATING][i];
					int et_component = code_table[COMPONENT][(et_edge%2? (et_edge-1)/2: et_edge/2)];
					int ot_component = code_table[COMPONENT][(ot_edge%2? (ot_edge-1)/2: ot_edge/2)];					
					int et_path = 2 * (et_component+1);
					int ot_path = 2 * (ot_component+1);
					int et_point = edge_path_offset[code_table[EVEN_TERMINATING][i]]+1;
					int ot_point = edge_path_offset[code_table[ODD_TERMINATING][i]]+1;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "write_metapost:   et_edge = " << et_edge << " ot_edge = " << ot_edge << " et_component = " << et_component << " ot_component = " << ot_component << endl;
	debug << "write_metapost:   et_path = " << et_path << " ot_path = " << ot_path << " et_point = " << et_point << " ot_point = " << ot_point << endl;
}					
					if (code_table[EVEN_TERMINATING][i] == code_table[EVEN_ORIGINATING][i])
						et_point++; // move et_point past the second type 4 vertex
					else if (code_table[ODD_TERMINATING][i] == code_table[ODD_ORIGINATING][i])
						ot_point++; // move ot_point past the second type 4 vertex

					if (mp_control.uniform_smoothed_discs)
					{
						os << "p" << 2*num_components+i+1 << " = fullcircle scaled " << mp_control.smoothed_state_disc_size << "d shifted z" << i << ";" << endl;
					}
					else
					{
						os << "r:= min(arclength (z" << i << "--(point " << et_point-1 << ".5 of p" << et_path << ")),";
						os <<         "arclength (z" << i << "--(point " << ot_point-1 << ".5 of p" << ot_path << ")),";
						os <<         "arclength (z" << i << "--(point " << et_point << ".5 of p" << et_path << ")),";
						os <<         "arclength (z" << i << "--(point " << ot_point << ".5 of p" << ot_path << ")));" << endl;

						os << "p" << 2*num_components+i+1 << " = fullcircle scaled min(2r," << mp_control.smoothed_state_disc_size << "d) shifted z" << i << ";" << endl;
						os << "theta := angle((direction " << et_point << " of p" << et_path << ")+(direction " << ot_point << " of p" << ot_path << "))-90;" << endl;
					}
										
					os << "fill p" << 2*num_components+i+1 << " withcolor 1white;" << endl;
																
					os << "z" << 2*num_vertices+4*i << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << et_point-1<< "," << et_point << ") of p" << et_path << ";" << endl;
					os << "z" << 2*num_vertices+4*i+1 << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << ot_point-1 << "," << ot_point <<  ") of p" << ot_path << ";" << endl;
					os << "z" << 2*num_vertices+4*i+2 << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << et_point << "," << et_point+1 << ") of p" << et_path << ";" << endl;
					os << "z" << 2*num_vertices+4*i+3 << "=p" << 2*num_components+i+1 << " intersectionpoint subpath(" << ot_point << "," << ot_point+1 << ") of p" << ot_path << ";" << endl;
					
					os << "draw z" << 2*num_vertices+4*i << "--z" << 2*num_vertices+4*i+3 << ";" << endl;
					os << "draw z" << 2*num_vertices+4*i+1 << "--z" << 2*num_vertices+4*i+2 << ";" << endl;
				}
			}
		}
	}
	
	if (mp_control.draw_labels)
	{
		int edge_label = 0;
		
		for (int i=0; i< num_edges; i++)
		{
			if (mp_control.gauss_labels && (code_table[LABEL][term_crossing[i]] == generic_code_data::VIRTUAL || (ignore_shortcut && shortcut_crossing[term_crossing[i]]))) 
				continue;
				
			os << "label(btex $";

			if (mp_control.scriptscript_labels)
				os << "\\textfont0=\\fiverm ";
			else if (mp_control.script_labels)
				os << "\\textfont0=\\sevenrm ";

			os << (mp_control.label_edges_from_one? edge_label+1: edge_label) << "$ etex, z" << vertex_sequence[2*i] << ");" << endl;
			edge_label++;
		}
	}

	if (mp_control.label_vertices)
	{
		for (int i=0; i< num_type1234_vertices; i++)
		{
			os << "label(btex $";
			if (mp_control.scriptscript_labels)
				os << "\\textfont0=\\fiverm \\textfont1=\\fivei z";
			else if (mp_control.script_labels)
				os << "\\textfont0=\\sevenrm \\textfont1=\\seveni z";
			else
				os << "z";
			
			os << i << "$ etex, z" << i << ");" << endl;
		}

		if (mp_control.show_vertex_axes)
		{
			os << "draw (0," << miny << "u)--(0," << maxy << "u);" << endl;
			os << "draw (" << minx << "u,0)--(" << maxx << "u,0);" << endl;
		}
			
		os << "fill fullcircle scaled 3pt shifted (0," << miny << "u);" << endl;
		os << "fill fullcircle scaled 3pt shifted (0," << maxy << "u);" << endl;
		os << "fill fullcircle scaled 3pt shifted (" << minx << "u,0);" << endl;
		os << "fill fullcircle scaled 3pt shifted (" << maxx << "u,0);" << endl;

		os << "label.bot(btex (0," << int(miny) << "u) etex,(0," << miny << "u));" << endl;
		os << "label.top(btex (0," << int(maxy) << "u) etex,(0," << maxy << "u));" << endl;
		os << "label.lft(btex (" << int(minx) << "u,0) etex,(" << minx << "u,0));" << endl;
		os << "label.rt(btex (" << int(maxx) << "u,0) etex,(" << maxx << "u,0));" << endl;

	}

	if (mp_control.gauss_crossings && auxiliary_data != 0)
	{
		vector<int>& gauss_crossing_map = *auxiliary_data;
		
		for (size_t i=0; i< gauss_crossing_map.size(); i++)
		{
			int crossing = gauss_crossing_map[i];

			/* identify the components and corresponding paths of the even terminating and odd terminating edges */
			int et_edge = code_table[EVEN_TERMINATING][crossing];
			int ot_edge = code_table[ODD_TERMINATING][crossing];
			int et_component = code_table[COMPONENT][(et_edge%2? (et_edge-1)/2: et_edge/2)];
			int ot_component = code_table[COMPONENT][(ot_edge%2? (ot_edge-1)/2: ot_edge/2)];					
			int et_path = 2 * (et_component+1);
			int ot_path = 2 * (ot_component+1);
			int et_point = edge_path_offset[code_table[EVEN_TERMINATING][crossing]]+1;
			int ot_point = edge_path_offset[code_table[ODD_TERMINATING][crossing]]+1;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "write_metapost:   et_edge = " << et_edge << " ot_edge = " << ot_edge << " et_component = " << et_component << " ot_component = " << ot_component << endl;
	debug << "write_metapost:   et_path = " << et_path << " ot_path = " << ot_path << " et_point = " << et_point << " ot_point = " << ot_point << endl;
}					
			if (code_table[EVEN_TERMINATING][i] == code_table[EVEN_ORIGINATING][i])
				et_point++; // move et_point past the second type 4 vertex
			else if (code_table[ODD_TERMINATING][i] == code_table[ODD_ORIGINATING][i])
				ot_point++; // move ot_point past the second type 4 vertex

			os << "theta := angle((direction " << et_point << " of p" << et_path << ")+(direction " << ot_point << " of p" << ot_path << "))-90;" << endl;

			os << "label (btex $"; 
			if (mp_control.scriptscript_labels)
				os << "\\textfont0=\\fiverm ";
			else if (mp_control.script_labels)
				os << "\\textfont0=\\sevenrm ";

		    os << i+1 << "$ etex, z" << crossing << "+(" << mp_control.label_shift << "u*cosd theta, " << mp_control.label_shift << "u*sind theta));" << endl;				    
		}
	}
	if (mp_control.circle_packing)
	{
		for (int i=0; i< num_type1234_vertices; i++)
			os << "draw fullcircle scaled " << 2*vertex_radius[i] << "u shifted z" << i << ";" << endl;
			
		if (INCLUDE_BOUNDARY_VERTICES)
		{			
			for (int i=num_type1234_vertices; i< num_vertices; i++)
				os << "z" << i << "=(" << setw(output_field_width) << vcoords[i][0] << "u," << setw(output_field_width) << vcoords[i][1] << "u);" << endl;
			
			os << endl;
			for (int i=num_type1234_vertices; i< num_vertices; i++)
				os << "draw fullcircle scaled " << 2*vertex_radius[i] << "u shifted z" << i << ";" << endl;
		}
		
		
	}

	os << "endfig;" << endl;
	Tinput.close();
}

void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, 
                    char const* circlepack_output_file, char const* circlepack_output_savefile, vector<int>* auxiliary_data=0)
{
	ifstream Cinput;
			
	Cinput.open(circlepack_output_file);
	if (!Cinput)
	{
			cout << "\nError opening output file " << circlepack_output_file << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: Error opening output file " << circlepack_output_file << endl;
		exit(0);
    }
	    
    /* If we're tracking force directed placement or doing magnification we need the vertex coordinates from both the 
       .out and .sav versions of circlepack_output_file, the .out versions are read first, then the .sav
    */
    int num_vertices;
    Cinput >> num_vertices;
    
    matrix<double> vcoords(2*num_vertices,2);
    vector<double> vertex_radius(num_vertices);

	for (int i=0; i< num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			Cinput >> vcoords[i][j];
	}
	for (int i=0; i< num_vertices; i++)
		Cinput >> vertex_radius[i];

	Cinput.close();
	
	ifstream Sinput;

	Sinput.open(circlepack_output_savefile);
	if (!Sinput)
	{
		cout << "\nError opening output file " << circlepack_output_savefile << endl;
		exit(0);
    }
    
    Sinput >> num_vertices;

	for (int i=0; i< num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			Sinput >> vcoords[num_vertices+i][j];
	}
	
	Sinput.close();

	write_metapost(output, code_data, title, mp_control, vcoords, vertex_radius, auxiliary_data);		
}

void print (metapost_control& mp_control, ostream& os, string prefix)
{
	os << prefix << "frame_minx = " << mp_control.frame_minx << endl;
	os << prefix << "frame_maxx = " << mp_control.frame_maxx << endl;
	os << prefix << "frame_miny = " << mp_control.frame_miny << endl;
	os << prefix << "frame_maxy = " << mp_control.frame_maxy << endl;
	os << prefix << "rotate = " << (mp_control.rotate? "true": "false") << endl;
	os << prefix << "explicit_rotation_centre = " << (mp_control.explicit_rotation_centre? "true": "false") << endl;
	os << prefix << "implicit_rotation_centre = " << (mp_control.implicit_rotation_centre? "true": "false") << endl;
	os << prefix << "rotation_degrees = " << mp_control.rotation_degrees << endl;
	os << prefix << "rotation_centre_x = " << mp_control.rotation_centre_x << endl;
	os << prefix << "rotation_centre_y = " << mp_control.rotation_centre_y << endl;
	os << prefix << "rotation_centre_z = " << mp_control.rotation_centre_z << endl;
	os << prefix << "unit_size = " << mp_control.unit_size << endl;
	os << prefix << "disc_size = " << mp_control.disc_size << endl;
	os << prefix << "pen_size = " << mp_control.pen_size << endl;
	os << prefix << "horizontal_units = " << mp_control.horizontal_units << endl;
	os << prefix << "vertical_units = " << mp_control.vertical_units << endl;
	os << prefix << "infinite_cycle = " << mp_control.infinite_cycle << endl;
	os << prefix << "cudgel_space = " << mp_control.cudgel_space << endl;
	os << prefix << "dash_with_dots = " << (mp_control.dash_with_dots? "true": "false") << endl;
	os << prefix << "knotoid_leg_unbounded = " << (mp_control.knotoid_leg_unbounded? "true": "false") << endl;
	os << prefix << "one_metapost_path = " << (mp_control.one_metapost_path? "true": "false") << endl;
	os << prefix << "colour = " << (mp_control.colour? "true": "false") << endl;
	os << prefix << "draw_lace_frame = " << (mp_control.draw_lace_frame? "true": "false") << endl;
	os << prefix << "draw_immersion = " << (mp_control.draw_immersion? "true": "false") << endl;
	os << prefix << "draw_crossing_features = " << (mp_control.draw_crossing_features? "true": "false") << endl;
	os << prefix << "draw_frame_corners = " << (mp_control.draw_frame_corners? "true": "false") << endl;
	os << prefix << "draw_triangulation = " << (mp_control.draw_triangulation? "true": "false") << endl;
	os << prefix << "draw_triangulation_displacement = " << (mp_control.draw_triangulation_displacement? "true": "false") << endl;
	os << prefix << "draw_labels = " << (mp_control.draw_labels? "true": "false") << endl;
	os << prefix << "label_edges_from_one = " << (mp_control.label_edges_from_one? "true": "false") << endl;
	os << prefix << "draw_oriented = " << (mp_control.draw_oriented? "true": "false") << endl;
	os << prefix << "draw_shortcut = " << (mp_control.draw_shortcut? "true": "false") << endl;
	os << prefix << "hamiltonians = " << (mp_control.hamiltonians? "true": "false") << endl;
	os << prefix << "label_vertices = " << (mp_control.label_vertices? "true": "false") << endl;
	os << prefix << "script_labels = " << (mp_control.script_labels? "true": "false") << endl;
	os << prefix << "scriptscript_labels = " << (mp_control.scriptscript_labels? "true": "false") << endl;
	os << prefix << "seifert_circles = " << (mp_control.seifert_circles? "true": "false") << endl;
	os << prefix << "show_vertex_axes = " << (mp_control.show_vertex_axes? "true": "false") << endl;
	os << prefix << "state = " << mp_control.state << endl;
	os << prefix << "state_smoothed = " << (mp_control.state_smoothed? "true": "false") << endl;
	os << prefix << "circle_packing = " << (mp_control.circle_packing? "true": "false") << endl;
	os << prefix << "draw_shrink_effect = " << (mp_control.draw_shrink_effect? "true": "false") << endl;
	os << prefix << "highlight_small_edges = " << (mp_control.highlight_small_edges? "true": "false") << endl;
	if (mp_control.tension)
		os << prefix << "tension = true, matapost_path_tension = " << metapost_path_tension << endl;
	else
		os << prefix << "tension = false" << endl;
	if (mp_control.colourmap.length())
		os << prefix << "colourmap = " << mp_control.colourmap << endl;
	else
		os << prefix << "colourmap = default" << endl;

	if (mp_control.singlecolour.length())
		os << prefix << "singlecolour = " << mp_control.singlecolour << endl;
	else
		os << prefix << "singlecolour = not specified" << endl;

	os << prefix << "draw_colour = ";
	for (unsigned int i=0; i< mp_control.draw_colour.size(); i++)
		os << mp_control.draw_colour[i] << ' ';
	os << endl;	
	os << prefix << "lace_midpoints = ";
	for (unsigned int i=0; i< mp_control.lace_midpoints.size(); i++)
		os << mp_control.lace_midpoints[i].first << '-' << mp_control.lace_midpoints[i].second << ' ';
	os << endl;
	os << prefix << "translations = ";
	for (unsigned int i=0; i< mp_control.translations.size(); i++)
		os << get<0>(mp_control.translations[i]) << '(' << get<1>(mp_control.translations[i]) << ',' << get<1>(mp_control.translations[i]) << ") ";
	os << endl;
	os << prefix << "hamiltonian-circuit = ";
	for (unsigned int i=0; i< mp_control.hamiltonian_circuit.size(); i++)
		os << mp_control.hamiltonian_circuit[i] << ' ' ;
	os << endl;
}

	
/* The function hyperbolic_representation produces the Euclidean centres and radii
   for a set of circles currently described in hyperbolic form and required to be 
   drawn in a representation of the disc model of the hyperbolic plane.

   For each circle, we determine the Mobius transformation that takes the hyperbolic 
   centre to zero.  There are two points of intersection, X and -X, of the image of
   circle under this transformation and the x-axis, and the hyperbolic distance between
   0 and X is the hyperbolic radius.
   
   From the formula log (1+t_2)(1-t_1)/(1-t_2)(1+t_1) for hyperbolic distance we can
   evaluate the Euclidean distance of X from zero as |X| = |(e^{r_h} - 1)/(e^{r_h}+1)|.
   
   The inverse of the Mobius transformation then takes X and -X to points z_1 and z_2.
   The Euclidean centre of the circle is the point (z_1+z_2)/2 and the Euclidean radius
   is |z_1-z_2|/2
   
   Note that the hyperbolic radii produced by the circle packing are recorded in the
   form e^{-2h} for the hyperbolic radius h.  The factor of 2 is based on the choice of 
   definition for the hyperbolic metric but note that 
   
   (e^{-d} - 1)/(e^{-d}+1) = -(e^{d} - 1)/(e^{d}+1)
   
   so we may calculate the value of |X| directly from radius[i].
   
*/
void hyperbolic_representation(matrix<double>& centre, vector<double>& radius)
{
	int num_circles = radius.size();
	
	for (int i=0; i< num_circles; i++)
	{
		/* set e = e^{r_h} */
//		double e = exp(radius[i]);

		/* set t to be the modulus of the hyperbolic centre */
//		double t = sqrt(centre[i][0]*centre[i][0]+centre[i][1]*centre[i][1]);
		
		complex<double> a;
		complex<double> b;                     
		complex<double> z(centre[i][0],centre[i][1]);
		
		mobius_transform_to_origin(a,b,z);
		
//		double x = (e - 1)/(e + 1)/2;
		double x = abs((radius[i] - 1)/(radius[i] + 1));
		complex<double> z1 = inv_mobius_transformation(a,b, complex<double>(x,0));
		complex<double> z2 = inv_mobius_transformation(a,b, complex<double>(-x,0));

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "hyperbolic_representation: circle " << i << ", centre " << z << ", radius " << radius[i] << endl;
	debug << "hyperbolic_representation: Mobius transformation parameters a = " << a << ", b = " << b << endl;
	debug << "hyperbolic_representation: Euclidean value of radius (X) = " << x << endl;
	debug << "hyperbolic_representation: pre-image of X = " << z1 << ", pre-image of -X = " << z2 << endl;
}    
		
		radius[i] = abs(z1-z2)/2*metapost_hyperbolic_scale_factor;

		z1 += z2;
		z1 /= 2;
		
		centre[i][0] = z1.real()*metapost_hyperbolic_scale_factor;
		centre[i][1] = z1.imag()*metapost_hyperbolic_scale_factor;
		

		/* we always evaluate tQ, even if the hyperbolic centre is the origin
		  because we need it for the Euclidean radius
		
		double tQ = e + e*t - 1 + t;
		double d  = e + e*t + 1 - t;
		tQ /= d;

		if (t != 0)
		{
			
			double tP = 1 + t - e + e*t; 
			d =         1 + t + e - e*t;
			tP /= d;
			
			double te = (tP + tQ)/2;
			
//			scale the hyperbolic centre to get the Euclidean centre
			double f = te/t;
			
			centre[i][0] *= f;
			centre[i][1] *= f;
		}
		*/
		
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "hyperbolic_representation: hyperbolic vertex_coordinates: " << endl;
    for (int i=0; i< num_circles; i++)
    {
		debug << "hyperbolic_representation:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << centre[i][j] << ' ';
		debug << endl;
	}	
    debug << "hyperbolic_representation: hyperbolic vertex_radius: " << endl;
    for (int i=0; i< num_circles; i++)
		debug << "hyperbolic_representation:   vertex " << i << ": " << radius[i] << endl;
}

}

/* The function mobius_transform_to_origin evaluates in a and b the parameters
   for the Mobius transformation that takes z to the origin
*/
void mobius_transform_to_origin(complex<double>& a, complex<double>& b, complex<double> z)
{
	double n=norm(z);
	double r_squared = 1/(1-n); // norm is the square of abs
	double r = sqrt(r_squared);
	double s = sqrt(r_squared - 1);
	b = complex<double>(s,0);
	double theta = 355/113 - arg(z);
	a = complex<double>(r*cos(theta), r*sin(theta));
	
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "mobius_transform_to_origin: norm = " << n << ", r_squared = " << r_squared << ", r = " << r << ", theta = " << theta << endl;
}

/* the function mobius_transformation evaluates the image of z under the Mobius transformation
   that fixing the unit disc determined by a and b
*/
complex<double> mobius_transformation(complex<double> a, complex<double> b, complex<double> z)
{
	return (a*z+b)/(conj(b)*z + conj(a));
}

/* the function mobius_transformation evaluates the image of z under the inverse of the
   Mobius transformation fixing the unit disc determined by a and b
*/
complex<double> inv_mobius_transformation(complex<double> a, complex<double> b, complex<double> z)
{
	return (b-conj(a)*z)/(conj(b)*z - a);
}

void set_edge_directions(generic_code_data& code_data, matrix<double>& vcoords, matrix<double>& edge_direction)
{   
	int num_crossings = code_data.num_crossings;
	matrix<int>& code_table	 = code_data.code_table;	
	double direction_length = 100; // used to evaluate direction vectors from adjusted direction angles

	for (int i = 0; i< num_crossings; i++)
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions: setting edge directions for crossing " << i << endl;

		/* identify the type-2 vertices in the mid-point of the semi-arcs entering the crossing
		   for the first and second time, v1 and v3, together with the corresponding  type-2 vertices
		   in the mid-point of the semi-arcs leaving the crossing, v2 and v4.  We evaluate the vertices as 
		   edges initially to aid debugging
		*/
		int e1, e3, v1, v2, v3, v4;
		if (code_table[ODD_TERMINATING][i] < code_table[EVEN_TERMINATING][i])
		{
			e1 = v1 = code_table[ODD_TERMINATING][i];
			v2 = code_table[EVEN_ORIGINATING][i];
			e3 = v3 = code_table[EVEN_TERMINATING][i];
			v4 = code_table[ODD_ORIGINATING][i];
		}
		else
		{
			e1 = v1 = code_table[EVEN_TERMINATING][i];
			v2 = code_table[ODD_ORIGINATING][i];
			e3 = v3 = code_table[ODD_TERMINATING][i];
			v4 = code_table[EVEN_ORIGINATING][i];
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   e1 = " << e1 << ", e2 = " << v2 << ", e3 = " << e3 << ", e4 = " << v4 << endl;
		
		v1 += num_crossings;
		v2 += num_crossings;
		v3 += num_crossings;
		v4 += num_crossings;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   v1 = " << v1 << ", v2 = " << v2 << ", v3 = " << v3 << ", v4 = " << v4 << endl;
		
		/* calculate the angles theta_i, there are four cases depending on the type of crossing and whether 
		   code_table[ODD_TERMINATING][i] < code_table[EVEN_TERMINATING][i] or not:
		
			         OT < ET                   OT < ET                   OT > ET                   OT > ET
		       
		        v1 OT \ 3 / OO v4         v3 ET \ 4 / EO v2         v3 OT \ 4 / OO v2         v1 ET \ 3 / EO v4
		               \ /                       \ /                       \ /                       \ /
		 Type 1      4  x  2       Type 2      1  x  3       Type 1      1  x  3       Type 2      4  x  2
		               / \                       / \                       / \                       / \
		        v3 ET / 1 \ EO v2         v1 OT / 2 \ OO v4         v1 ET / 2 \ EO v4         v3 OT / 1 \ OO v2
		        
		 v-cycles:  v3 v2 v4 v1               v3 v1 v4 v2               v3 v1 v4 v2               v3 v2 v4 v1
		 
		 The v-cycle shows the anti-clockwise cyclic order of the vertices vi around the crossing starting from v3,
		 as can be seen there are just two combinations.		        
		*/
		double theta_1,theta_2,theta_3,theta_4;
		
		if (   (code_table[ODD_TERMINATING][i] < code_table[EVEN_TERMINATING][i] && code_table[TYPE][i] == generic_code_data::TYPE1)
		    || (code_table[ODD_TERMINATING][i] > code_table[EVEN_TERMINATING][i] && code_table[TYPE][i] == generic_code_data::TYPE2)
		   )
		{
			calculate_theta_angles(v3, v2, v4, v1, i, theta_1, theta_2, theta_3, theta_4, vcoords);
		}
		else
		{
			calculate_theta_angles(v3, v1, v4, v2, i, theta_1, theta_2, theta_3, theta_4, vcoords);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "set_edge_directions:   theta_1 = " << theta_1 << ", theta_2 = " << theta_2 << ", theta_3 = " << theta_3 << ", theta_4 = " << theta_4 << endl;
    debug << "set_edge_directions:     check Sum theta_i = " << theta_1 + theta_2 + theta_3 + theta_4 << endl;
}
		
		/* evaluate the new direction for edge e3 */
		double chi = edge_argument(v3,i,vcoords);

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   initial chi = " << chi << endl;

		double delta = (theta_4 - theta_1)/2; // = (theta_1 + theta_4)/2 - theta_1
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   delta = " << delta << endl;
		
		double new_chi = chi - delta;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   new chi = " << new_chi << ", cos(new_chi) = " << cos(new_chi) << ", sin(new_chi) = " << sin(new_chi) << endl;

		edge_direction[e3][0] = direction_length*cos(new_chi);
		edge_direction[e3][1] = direction_length*sin(new_chi);

		if (abs(edge_direction[e3][0]) < ZERO_THRESHOLD)
			edge_direction[e3][0] = 0;
		if (abs(edge_direction[e3][1]) < ZERO_THRESHOLD)
			edge_direction[e3][1] = 0;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   set edge direction for vertex " << v3 << " = (" << edge_direction[e3][0] << "," << edge_direction[e3][1] << ")" << endl;
			
		/* evaluate angle phi_1 */
		double phi_1;
		
		if (   (code_table[ODD_TERMINATING][i] < code_table[EVEN_TERMINATING][i] && code_table[TYPE][i] == generic_code_data::TYPE1)
		    || (code_table[ODD_TERMINATING][i] > code_table[EVEN_TERMINATING][i] && code_table[TYPE][i] == generic_code_data::TYPE2)
		   )
		{
			/* e3 lies adjacent to e1 anti-clockwise around u */
			phi_1 = (theta_1+theta_4)/2;
		}
		else
		{
			/* e3 lies adjacent to e1 clockwise around u */
			phi_1 = (theta_2+theta_3)/2;
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   phi_1 = " << phi_1 << endl;

		/* evaluate the new direction for edge e1 */
		double chi_prime = edge_argument(v1,i,vcoords);

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   initial chi' = " << chi_prime << endl;

		double delta_prime = pi/2 - phi_1;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   delta' = " << delta_prime << endl;

		double new_chi_prime = chi_prime - delta_prime;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   new chi' = " << new_chi_prime << endl;

		edge_direction[e1][0] = direction_length*cos(new_chi_prime);
		edge_direction[e1][1] = direction_length*sin(new_chi_prime);
		
		if (abs(edge_direction[e1][0]) < ZERO_THRESHOLD)
			edge_direction[e1][0] = 0;
		if (abs(edge_direction[e1][1]) < ZERO_THRESHOLD)
			edge_direction[e1][1] = 0;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_edge_directions:   set edge direction for vertex " << v1 << " = (" << edge_direction[e1][0] << "," << edge_direction[e1][1] << ")" << endl;

	}
}

/* calculate_theta_angles calculates the angles between edges e1, e2, e3, e4 in the triangulation joining a type-2 vertex to the
   vertex (crossing) v, around the vertex v; that is it calculates the angles e1 v e2, e2 v e3, e3 v e4 and e4 v e1, 
   recording them in theta_1, theta_2, theta_3, theta_4, respectively.  The call to calculate_theta_angles
   should set e1, e2, e3, e4 to the number of the type-2 vertex joined by the corresponding edge to v.
*/
void calculate_theta_angles(int e1, int e2, int e3, int e4, int v,double& theta_1, double& theta_2, double& theta_3, double& theta_4, matrix<double>& vcoords)
{	

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "calculate_theta_angles: v = " << v << ", e1 = " << e1 << ", e2 = " << e2 << ", e3 = " << e3 << ", e4 = " << e4 << endl;
	
	/* calculate the arguments of the edges joining ei to v */
	double a1 = edge_argument(v,e1,vcoords);
	double a2 = edge_argument(v,e2,vcoords);
	double a3 = edge_argument(v,e3,vcoords);
	double a4 = edge_argument(v,e4,vcoords);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "calculate_theta_angles: arguments e1: " << a1 << ", e2: " << a2 << ", e3: " << a3 << ", e4: " << a4 << endl;
	
	double b1 = a2-a1;
	double b2 = a3-a2;
	double b3 = a4-a3;
	double b4 = a1-a4;

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "calculate_theta_angles: deltas " << b1 << ", " << b2 << ", " << b3 << ", " << b4 << endl;

	theta_1 = (b1<0? b1+two_pi: b1);
	theta_2 = (b2<0? b2+two_pi: b2);
	theta_3 = (b3<0? b3+two_pi: b3);
	theta_4 = (b4<0? b4+two_pi: b4);

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "calculate_theta_angles: theta_1 = " << theta_1 << ", theta_2 = " << theta_2 << ", theta_3 = " << theta_3 << ", theta_4 = " << theta_4 << endl;

}

/* edge_length determines the length of the line joining points (px,py) and (qx,qy) */
double edge_length(double px,double py,double qx,double qy)
{
	double dx = qx-px;
	double dy = qy-py;
	return sqrt(dx*dx+dy*dy);
}

/* cosine_rule_angle calculates the angle between edges of length a and b in a triangle having edge lengths a, b, and c 
*/
double cosine_rule_angle (double a, double b, double c)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "cosine_rule_angle: edge lengths a = " << a << ", b = " << b << ", c = " << c << endl;

	double cosine = (a*a+b*b-c*c)/2/a/b;

	if (cosine > 1)
	{
		cosine = 1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "cosine_rule_angle: cosine = " << cosine << ", adjusting to 1 " << endl;	
	}
	else if (cosine < -1)
	{
		cosine = -1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "cosine_rule_angle: cosine = " << cosine << ", adjusting to -1 " << endl;
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "cosine_rule_angle: cosine = " << cosine << endl;	

	double alpha = acos(cosine);  // the principle value of acos is [0,\pi], which is correct for the angles of a triangle
						
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "cosine_rule_angle: angle = " << alpha << " = " << alpha * 360 /two_pi << " degrees" << endl;
	
	return alpha;
}

/* edge argument calculates the argument of the vector u-v with respect to the
   coordinate x-axis, in the range [0,2\pi)
*/
double edge_argument(int v, int u, matrix<double>& vcoords)
{

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: v = " << v << ", u = " << u << endl;

	double x_delta = vcoords[u][0] - vcoords[v][0];
	double y_delta = vcoords[u][1] - vcoords[v][1];
	
	return edge_argument(x_delta, y_delta);
}

/* separated for a call from convex.cpp */
double edge_argument(double x_delta, double y_delta)
{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: x_delta = " << x_delta << ", y_delta = " << y_delta << endl;

	if (abs(x_delta) < ZERO_THRESHOLD && abs(y_delta) < ZERO_THRESHOLD)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: x_delta = 0, y_delta = 0, returning 0" << endl;
	
		return 0;
	}
	else if (abs(x_delta) < ZERO_THRESHOLD)
	{
		if (y_delta < 0)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: x_delta = 0, y_delta < 0, returning 3/2 pi" << endl;
			
			return 3*pi/2;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: x_delta = 0, y_delta > 0, returning pi/2" << endl;
			
			return pi/2;
		}
	}
	else if (abs(y_delta) < ZERO_THRESHOLD)
	{
		if (x_delta < 0)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: y_delta = 0, x_delta < 0, returning pi" << endl;
			
			return pi;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: y_delta = 0, x_delta > 0, returning 0" << endl;
			
			return 0;
		}
	}
	else
	{
		double alpha = atan(abs(y_delta)/abs(x_delta));
		double beta = 0;
		if (x_delta > 0 && y_delta > 0)
			beta = alpha;
		else if (x_delta < 0 && y_delta > 0)
			beta = pi - alpha;
		else if (x_delta < 0 && y_delta < 0)
			beta = pi + alpha;
		else // (x_delta > 0 && y_delta < 0)
			beta = two_pi - alpha; 

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "edge_argument: alpha = " << alpha <<  ", beta = " << beta << endl;
			
		return beta;
	}
}

void rotate_metapost (matrix<double>& coords, int num_coords_to_rotate, metapost_control& mp_control)
{
	/* we convert to radians by calculating rotation_degrees/360 * 2 pi and use 355/113 as an approximation to pi */

	double radians = mp_control.rotation_degrees * 355 / 180 / 113;
	double M00 = cos(radians);
	double M10 = sin(radians);
	double M01 = M10 * -1;
	double M11 = M00;

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "rotate_metapost: cos(rotation_degrees = " << radians << ") = " << M00 << ", sin(rotation_degrees) = " << M10 << endl;

	for (int i=0; i< num_coords_to_rotate; i++)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "rotate_metapost: rotating z" << i << " = (" << coords[i][0] << "," << coords[i][1] << ")" << endl;
    
		double shift_x = coords[i][0] - mp_control.rotation_centre_x;
		double shift_y = coords[i][1] - mp_control.rotation_centre_y;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "rotate_metapost:   shifted to (" << shift_x << "," << shift_y << ")" << endl;

		coords[i][0] = M00*shift_x + M01*shift_y;
		coords[i][1] = M10*shift_x + M11*shift_y;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "rotate_metapost:   rotated to (" << coords[i][0] << "," << coords[i][1] << ")" << endl;
    
		coords[i][0] += mp_control.rotation_centre_x;
		coords[i][1] += mp_control.rotation_centre_y;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "_metapost:   shifted back to (" << coords[i][0] << "," << coords[i][1] << ")" << endl;
			
	}
}


vector<int> gauss_parity(generic_code_data& code_data)
{
	int num_crossings = code_data.num_crossings;
	gauss_orientation_data gauss_data(code_data);
	
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "gauss_parity: presented with generic code ";
	write_code_data(debug, code_data);
	debug << endl;
	debug << "gauss_parity: gauss data ";
//	write_gauss_data(gauss_data, debug);
	print_gauss_data(gauss_data, debug,"gauss_parity: ");
}
	int num_terms = gauss_data.num_terms;
	int num_gauss_crossings = num_terms/2;
	matrix<int>& orientation_matrix = gauss_data.orientation_matrix;
	
	vector<int> parity(num_crossings);
	
	for (int i=0; i< num_gauss_crossings; i++)
	{
		int immersion_crossing = gauss_data.immersion_crossing[i];
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "gauss_parity: gauss_crossing = " << i+1 << " immersion_crossing = " << immersion_crossing << endl;
	
		int start=0;
		for (int j=0; j< num_terms; j++)
		{
			if (orientation_matrix[1][j] == i+1)
			{
				start=j;
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "gauss_parity:   first ocurrence at index " << start << endl;
				break;
			}
		}
		
		for (int k=1; k< num_terms; k++)
		{
			if (orientation_matrix[1][start+k] == i+1)
			{
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "gauss_parity:   second ocurrence at index " << start+k << endl;
				if (k%2)
					parity[immersion_crossing] = gauss_orientation_data::parity::EVEN;
				else
					parity[immersion_crossing] = gauss_orientation_data::parity::ODD;
				
				break;
			}
		}
	}

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
{
	debug << "gauss_parity: crossing number ";
	for (int i=0; i< num_crossings; i++)
		debug << i << " ";
	debug << endl;
	debug << "gauss_parity: parity          ";
	for (int i=0; i< num_crossings; i++)
	{
		if (parity[i] == gauss_orientation_data::parity::ODD)		
			debug << "O ";
		else if (parity[i] == gauss_orientation_data::parity::EVEN)		
			debug << "E ";
		else
			debug << "N ";
	}
	debug << endl;
}
	
	return parity;
}


/* identify_gauss_crossings returns a vector, gauss_crossing_map that maps Gauss crossing i to  gauss_crossing_map[i] */
vector<int> identify_gauss_crossings(generic_code_data& code_data)
{

	matrix<int>& code_table = code_data.code_table;
	int num_crossings = code_data.num_crossings;
	
	vector<int>& term_crossing = code_data.term_crossing;
	int num_edges = 2*num_crossings;
	vector<int> edge_flag(num_edges); // initialized to zero

		
	int num_classical_crossings = num_crossings;
	bool pure_knotoid_code_data = false;
	vector<int>& shortcut_crossing = code_data.shortcut_crossing;
	
	for (int i=0; i< num_crossings; i++)
	{
		if (code_table[LABEL][i] == generic_code_data::VIRTUAL)
			num_classical_crossings--;
	}
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_crossings: num_classical_crossings = " << num_classical_crossings << endl;

	if (code_data.immersion == generic_code_data::character::PURE_KNOTOID && code_data.head != -1 && shortcut_crossing.size())
	{
		pure_knotoid_code_data = true;
		for (unsigned int i=0; i< shortcut_crossing.size(); i++)
		{
			if (shortcut_crossing[i] != 0)
				num_classical_crossings--;
		}
if (debug_control::DEBUG >= debug_control::EXHAUSTIVE)
	debug << "identify_gauss_crossings: knotoid: num_classical_crossings = " << num_classical_crossings << endl;
	}

	int num_classical_crossings_visited = 0;
		
	/* gauss_crossing_map records the immersion crossings corresponding to Gauss crossings
	   so that gauss_crossing_map[i] is the immersion crossing corresponding to the ith Gauss crossing.
	   
	   crossing_visited will be a flag to indicate whether the ith immersion crossing
	   has been recorded in gauss_crossing_map or not.
	*/
	vector<int> gauss_crossing_map(num_classical_crossings);
	vector<int> crossing_visited(num_crossings);

	int start=0;
	int edge=0;
	bool complete = false;
		
	do 
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_crossings: component_start = " << start << endl;
			
		/*	trace this component */
		do
		{	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_crossings: edge = " << edge;
			edge_flag[edge] = 1;
			int next_crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", next_crossing = " << next_crossing;
		
			if (code_table[LABEL][next_crossing] == generic_code_data::VIRTUAL)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is virtual" << endl;
			}
			else if (pure_knotoid_code_data && shortcut_crossing[next_crossing])
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", is a shortcut crossing" << endl;
			}
			else if ((edge%2 != 0 && code_table[LABEL][next_crossing] == generic_code_data::POSITIVE) ||
			    (edge%2 == 0 && code_table[LABEL][next_crossing] == generic_code_data::NEGATIVE))
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", going under" << endl;	
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << ", going over" << endl;
			}
				
			if (code_table[LABEL][next_crossing] != generic_code_data::VIRTUAL && !(pure_knotoid_code_data && shortcut_crossing[next_crossing]))
			{
				if(crossing_visited[next_crossing])
				{
					for (int i=0; i< num_classical_crossings_visited; i++)
					{
						if (gauss_crossing_map[i] == next_crossing)
						{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_crossings:   second visit to Gauss crossing " << i+1<< endl;
							break;
						}
					}
				}
				else
				{

					gauss_crossing_map[num_classical_crossings_visited] = next_crossing;
					crossing_visited[next_crossing] = 1;
					num_classical_crossings_visited++;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_crossings:   first visit, becomes Gauss crossing " << num_classical_crossings_visited << endl;
				}

				if (edge%2)
					edge = code_table[EVEN_ORIGINATING][next_crossing];
				else
					edge = code_table[ODD_ORIGINATING][next_crossing];					
			}
			else
			{
				/* just move on around the component */				
				if (edge%2)
					edge = code_table[EVEN_ORIGINATING][next_crossing];
				else
					edge = code_table[ODD_ORIGINATING][next_crossing];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "identify_gauss_crossings:   doing nothing" << endl;
			}				
		} while (edge != start);

		/* look for the start of another component */
		complete = true;
		for (int i=0; i< num_edges; i++)
		{
			if (edge_flag[i] == 0)
			{
				complete = false;
				start = i;
				edge = start;
				break;
			}
		}			
	} while (!complete);
	
	return gauss_crossing_map;
}
