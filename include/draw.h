
/**********************************************************************

This header file defines local programme types used to control objects
used by the programme 
**********************************************************************/
#include <tuple>

#include <debug-control.h>
#include <generic-code.h>

#define CHANGE 

/*
#define TYPE1   -1
#define TYPE2   1
#define FLAT 2
#define POSITIVE 1
#define NEGATIVE -1
#define VIRTUAL 0

#define EPEER 0
#define OPEER 1
#define TYPE    2
#define LABEL	3
#define EVEN_TERMINATING 4
#define ODD_TERMINATING 5
#define EVEN_ORIGINATING 6
#define ODD_ORIGINATING 7
#define COMPONENT 8
#define CODE_TABLE_SIZE 9
*/


#define IMMERSION_CODE true
#define PEER_CODE false

class metapost_control
{

public:

	enum smoothed 
	{
		NON_SEIFERT_SMOOTHED = -1,
		SEIFERT_SMOOTHED = 1,
		A_SMOOTHED = 2,
		B_SMOOTHED = 3
	};


	bool 	adjacent_cudgel_midpoints;
	int 	arrowhead_bp_size;	
	bool 	circle_packing;
	bool 	colour;
	string 	colourmap;
	float 	cudgel_space;
	int 	cusp_disc_size;
	bool 	dash_with_dots;
	int 	disc_size;
	bool 	draw_crossing_features;
	bool 	draw_frame_corners;
	bool 	draw_grid;
	bool 	draw_immersion;
	bool 	draw_labels;
	bool 	draw_lace_frame;
	bool 	draw_oriented;
	bool 	draw_shortcut;
	bool 	draw_shrink_effect;
	bool 	draw_triangulation;
	bool 	draw_triangulation_displacement;
	bool 	explicit_rotation_centre;
	double 	frame_maxx;
	double 	frame_maxy;
	double 	frame_minx;
	double 	frame_miny;
	bool 	gauss_crossings;
	bool 	gauss_labels;
	int 	grid_size;
	int 	HC_include_edge;
	bool 	hamiltonians;
	string 	hamiltonian_colour;
	bool 	highlight_small_edges;
	int 	horizontal_units;
	bool 	implicit_rotation_centre;
	int 	infinite_cycle;
	bool 	knotoid_leg_unbounded;
	bool 	label_edges_from_one;
	int 	label_shift;
	bool 	label_vertices;
	bool 	left_terminating_tail_points;
	bool 	magnify_small_circles;
	bool 	midpoints_not_tail_points;
	float 	midpoint_tension;
	int 	odd_parity_disc_size;
	bool 	one_metapost_path;
	int 	pen_size;
	bool 	right_originating_tail_points;
	bool 	rotate;
	double 	rotation_centre_x;
	double 	rotation_centre_y;
	int 	rotation_centre_z;		
	double 	rotation_degrees;
	float 	scale; 
	bool 	script_labels;
	bool 	scriptscript_labels;
	bool 	seifert_circles;
	int 	seifert_edges;
	bool 	show_odd_parity_crossings;
	bool 	show_vertex_axes;
	string 	singlecolour;
	bool 	smallarrowheads;
	int 	smoothed_disc_threshold;
	int 	smoothed_state_disc_size;
	string 	state;
	bool 	state_smoothed;
	bool 	tension;
	bool 	uniform_smoothed_discs;
	int 	unit_size;
	int 	vertical_units;

	vector<string> draw_colour;
	vector<int> hamiltonian_circuit;
	vector<pair<int,int> > lace_midpoints;
	vector<tuple<int,int,int> > translations;
	
	metapost_control(): 	
		adjacent_cudgel_midpoints(true),arrowhead_bp_size(6),
		circle_packing(false),colour(false),colourmap(""),cudgel_space(1.0),cusp_disc_size(7),	
		dash_with_dots(false),disc_size(30),draw_crossing_features(true),
		draw_frame_corners(false),draw_grid(false),draw_immersion(true),
		draw_labels(false),draw_lace_frame(false),draw_oriented(false), 
		draw_shortcut(false),draw_shrink_effect(false),draw_triangulation(false),draw_triangulation_displacement(true),
		explicit_rotation_centre(false), 
		frame_maxx(0),frame_maxy(0),frame_minx(0),frame_miny(0),
		gauss_crossings(false),gauss_labels(false),grid_size(10), 
		HC_include_edge(-1), 
		hamiltonians(false),hamiltonian_colour("green"),highlight_small_edges(false),horizontal_units(5), 
		implicit_rotation_centre(false),infinite_cycle(-1),
		knotoid_leg_unbounded(false),
		label_edges_from_one(false),label_shift(50),label_vertices(false),left_terminating_tail_points(false),
		magnify_small_circles(false),midpoints_not_tail_points(false),midpoint_tension(1.0),
		odd_parity_disc_size(12),one_metapost_path(true),
		pen_size(1),
		right_originating_tail_points(true),rotate(false),rotation_centre_x(0),rotation_centre_y(0),rotation_centre_z(0),rotation_degrees(0),
		scale(0),script_labels(false), scriptscript_labels(false),seifert_circles(false),seifert_edges(0),show_odd_parity_crossings(false),show_vertex_axes(true),
		singlecolour(""),smallarrowheads(false),smoothed_disc_threshold(30),smoothed_state_disc_size(6),state(""),state_smoothed(false), 
		tension(false),
		uniform_smoothed_discs(false),unit_size(20),
		vertical_units(5),
		draw_colour({"red","blue","ForestGreen","Brown","DarkViolet","Orange"}){}	

};

void print (metapost_control& mp_control, ostream& os, string prefix);


/* the class of edge_info objects are used by edge distribution placement */
class edge_info
{
	public:
	
	int index;
	int vertex_0;
	int vertex_1;
	double length;
	int division_index;
	
	edge_info(): vertex_0(0),vertex_1(0),length(0),division_index(0) {}
};

struct lace_control 
{	
	enum polarity {LEFT,RIGHT};
	enum boundary_type {TAIL,HEAD,BONE};
};


class lace_skeleton
{
public:
	int num_cudgels;
	matrix<int> bones;
	vector<int> head_polarity;
	
	lace_skeleton(): num_cudgels(0){}
	lace_skeleton(int n): num_cudgels(n), bones(matrix<int> (n,2)), head_polarity(vector<int>(n)) {}
};

void read_lace_skeleton(lace_skeleton& skeleton, string s);
ostream& operator << (ostream& os, const lace_skeleton& skeleton);

class lace_body
{
public:
    lace_skeleton skeleton;
	vector<int> loop_location;
	vector<int> interior_1_cell_count;
	int num_loops;
	int num_boundary_1_cells;
	int num_interior_1_cells;
	int num_2_cells;
	int num_boundary_points;
	
	lace_body(): num_loops(0), num_boundary_1_cells(0), num_interior_1_cells(0), num_2_cells(0), num_boundary_points(0){}
};

void read_lace_body (lace_body& body, string s);
ostream& operator << (ostream& os, const lace_body& body);
void print_lace_body (ostream& os, const lace_body& body, string prefix);

/* these are the column indicators for the head and tail point counts and markers */
#define EVEN_TAIL_LT 1
#define EVEN_TAIL_LO 0
#define EVEN_TAIL_UT 3
#define EVEN_TAIL_UO 2

#define ODD_TAIL_LT 0
#define ODD_TAIL_LO 1
#define ODD_TAIL_UT 2
#define ODD_TAIL_UO 3

#define EVEN_HEAD_L 0
#define EVEN_HEAD_T 2
#define EVEN_HEAD_O 1

#define ODD_HEAD_L 0
#define ODD_HEAD_T 1 
#define ODD_HEAD_O 2

