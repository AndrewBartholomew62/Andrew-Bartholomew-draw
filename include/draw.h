
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
*/

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


	double frame_minx;
	double frame_maxx;
	double frame_miny;
	double frame_maxy;
	bool rotate;
	bool explicit_rotation_centre;
	bool implicit_rotation_centre;
	double rotation_degrees;
	double rotation_centre_x;
	double rotation_centre_y;
	int arrowhead_bp_size;	
	int disc_size;
	int grid_size;
	int horizontal_units;
	int infinite_cycle;
	int pen_size;
	int rotation_centre_z;	
	int smoothed_state_disc_size;
	int unit_size;
	int vertical_units;
	float cudgel_space;
	bool dash_with_dots;
	bool knotoid_leg_unbounded;
	bool one_metapost_path;
	bool colour;
	bool draw_lace_frame;
	bool draw_immersion;
	bool draw_crossing_features;
	bool draw_frame_corners;
	bool draw_grid;
	bool draw_triangulation;
	bool draw_triangulation_displacement;
	bool draw_labels;
	bool label_edges_from_one;
	bool draw_oriented;
	bool draw_shortcut;
	bool gauss_labels;
	bool label_vertices;
	bool script_labels;
	bool scriptscript_labels;
	bool show_odd_parity_crossings;
	bool show_vertex_axes;
	bool smallarrowheads;
	bool circle_packing;
	bool draw_shrink_effect;
	bool highlight_small_edges;
	bool tension;
	bool adjacent_cudgel_midpoints;
	bool right_originating_tail_points;
	bool left_terminating_tail_points;
	bool magnify_small_circles;
	bool midpoints_not_tail_points;
	float midpoint_tension;
	string colourmap;
	string singlecolour;
	string state;
	vector<string> draw_colour;
	vector<pair<int,int> > lace_midpoints;
	vector<tuple<int,int,int> > translations;
	
	metapost_control(): frame_minx(0), frame_maxx(0), frame_miny(0),frame_maxy(0),rotate(false), explicit_rotation_centre(false), implicit_rotation_centre(false),
	                    rotation_degrees(0), rotation_centre_x(0) , rotation_centre_y(0), 
	                    arrowhead_bp_size(6),disc_size(30),grid_size(10), horizontal_units(5), infinite_cycle(-1),
	                    pen_size(1),rotation_centre_z(0),smoothed_state_disc_size(6), unit_size(20), vertical_units(5),
						cudgel_space(1.0), dash_with_dots(false), knotoid_leg_unbounded(false), one_metapost_path(false),
						colour(false), draw_lace_frame(false), draw_immersion(true), draw_crossing_features(true), draw_frame_corners(false), draw_grid(false),
						draw_triangulation(false), draw_triangulation_displacement(true), draw_labels(false), label_edges_from_one(false), draw_oriented(false), 
						draw_shortcut(false), gauss_labels(false), label_vertices(false), script_labels(false), scriptscript_labels(false), show_odd_parity_crossings(false), 
						show_vertex_axes(true), smallarrowheads(false),circle_packing(false), draw_shrink_effect(false),
	                    highlight_small_edges(false),tension(false),adjacent_cudgel_midpoints(true), right_originating_tail_points(true),
						left_terminating_tail_points(false), magnify_small_circles(false),midpoints_not_tail_points(false),midpoint_tension(1.0),colourmap(""),singlecolour(""),
	                    state(""),draw_colour({"red","blue","ForestGreen","Brown","DarkViolet","Orange"}){}	

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

