/**************************************************************************
draw generates metapost code to draw diagrams from labelled peer codes, or
labelled immersion codes.  It is currently able to draw classical and virtual 
knots and knotoids.  It was motivated by the requirement to draw knotoids.

The original algorithm was taken from the drawing routine in Knotscape, by Jim Hoste and Morwen 
Thistlethwaite, and in particular Ken Stephenson's implementation of Thurston's circle 
packing algrithm.  Subsequently Roger Fenn's circle packing algorithm and an implementation of
force directed placement of vertices was added in an attempt to improve the aesthetics of the
diagrams produced by the code.

Notes to the initial version based on Knotscape
-----------------------------------------------

My code provides the environment to read peer or immersion codes from an input stream, then 
applies the algorithms of Knotscape to produce a sequence of coordinates at key 
points as we traverse the immersion.  Armed with these coordinates it is a simple 
matter to produce metapost output that can draw the immersion with the appropriate 
over, under and virtual crossings, knotoid shortcuts, and so on.

Specifically, the following key aspects have been taken from Knotscape.  I am very 
grateful to the authors of Knotscape for making this code available in the public
domain and fully acknowledge all rights they may have in that code and it's methods.

    triangulate.c - I have written code to produce the same triangulation as described in 
    this Knotscape module starting from a labelled peer or immersion code.  The output of my
    code is a file with the same file format as that produced by the Knotscape version.
    
    ken.c - this module has been used essentially as is from Knotscape.  I have renamed the
    module circle_pack.c and have modified the writepack function to include some additional
    information useful for writing the metapost code.  I have also added some debug output 
    but the algorithms used by this module are untouched and the module is compiled from C.
    
	nodeseq.c - this module has been re-written to take into account the fact that the 
	sequence of vertices I am interested in is determined by a labelled peer or immersion code;
	also, the use of metapost means I do not need all of the coordinates used by Knotscape.
	Thus, my nodeseq is a form of the corresponding Knotscape model optimized for peer and 
	immersion codes, however, as with the other modules above, the algorithms (in particular 
	the coordinate scaling algorithm) are all taken from Knotscape.
	
	badness - I have incorporated the same assessment of a set of coordinates as used by
	Knotscape via the badness function.
	
	In addition to the above, wherever the Knotscape original modules call auxiliary functions
	I have written corresponding functions in support of my versions.  These are also 
	acknowledged effctively to have been taken directly from Knotscape.

Notes to subsequent additions
-----------------------------

Roger Fenn's circle packing algorithm was implemented in circle_pack.cpp and a boolean control
variable, USE_KEN_STEPHENSON_CIRCLE_PACKING, introduced to allow the use of the Knotscape algorithm
if desired.

Similarly, the boolean variable USE_FORCE_DIRECTED_PLACEMENT was added to apply force direction to an
initial placement of triangulation vertices.  The selected circle packing algorithm is used to calculate
this initial placement and then we apply forces to the vertices.

For both circle packing and force directed placement we may choose to include the triangulation boundary 
(type 4) vertices or not.  The default is not to include them.

For force directed placement we may choose to apply repulsive forces across all pairs of vertices or only
to vertices joined by an edge of the triangulation.

                  Coding started 21st April 2010
 
   Version 1:    support for Knotoids; identified by extending peer codes and immersion codes to allow a '^'
                 character following a crossing number to indicate the first shortcut crossing.
   Version 2:    added immersions and rotation
   Version 3:    added labelled peer codes and Metapost control qualifiers
   Version 4:    added drawing of the triangulation underlying the circle packing.
   Version 5:    added Fruchterman and Reingold force directed placement 
                 added edge magnification
   Version 6:    added inner hull and centre of gravity
   Version 7:    ported Ken Stephenson's circle packing algorithm to C++ and
                 added check for connected code data.
   Version 8:    replaced Fruchterman and Reingold force directed placement 
                 with the Plestenjak modifications
                 updated the triangulation for diagrams having two edges in the infinite turning cycle
   Version 9:    added capability to draw straight-edge disc triangulations
                 from either peer-code triangulations or dedicated input files.
                 coding started 7-12-17
   Version 10:   reinstated a form of  Fruchterman and Reingold force directed placement where adjacent 
                 type 1 and trype 2 vertices that are too close to each other repel
   Version 11:   added region shrinking placement, rearranged and improved the use of long and short
                 programme options, coding started 25-10-18
   Version 12:   added the ability to the draw smoothed states for a given knot, link knotoid or multi-knotoid diagram
   Version 13:   added lace drawing capability and the ability to draw diagrams with colour, added support for diagrams containing Reidemeister I loops (December 2021)
   Version 14:   added the ability to draw knots from Gauss codes and added the smallarrowheads option  (December 2021)
   Version 15:   added ability to label Gauss-arcs rather than immersion edges (January 2022)
   Version 16:   added the ability to draw a single smoothed state rather than all of them (February 2022)
   Version 17:   added planar diagram support and optimized the Gauss to peer code conversion (March 2022)
   Version 17.1: added automatic sizing of smoothed crossings and the ability to shift smoothed crossings labels to one side (March 2022)
   Version 18:   added labelling of Gauss crossings, as specified by Gauss code or planar diagram input, or calculated from peer code input (October 2022)
   Version 18.1  code tidy-up, moved to initialize.h
   Version 19:   added the ability to draw multi-linkoids, added shorcuts to knot-type knoid drawings.  
                 added scale option to re-size diagrams so they become  proportionate to a fixed diagonal size (August 2023)
   Version 20:   added the ability to draw Seifert circles and Hamiltonian circuits (October 2023)
   Version 21:   improved the drawing of multi-linkoids, updated and made default the use of a single metapost path per component 
                 added seifert-edges option (November 2023)
   Version 22:   added support for the Gauss code format used by J. Chen's https://www.flatknotinfo.com (November 2024)
   Version 23:   Changed the syntax for multi-linkoids to use '%' for knot-type components.  Introuduced support for singular crossings
                 using the label '@' (December 2024)

**************************************************************************/
using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <stdio.h>
#include <valarray>
#include <ctime>
#include <csignal>
#include <list>
#include <iomanip>
#include <vector>

/******************** Global Variables **********************/
string 		version  = "23";
   

extern ofstream    debug;
ofstream 	output;
ifstream    input;

#include <util.h>
#include <input-control.h>
#include <matrix.h>
#include <draw.h>
#include <gauss-orientation.h>

int printf_bool=0; // debug for circle_pack.c - deprecated but retained in case of future testing

bool DRAW_IN_HYPERBOLIC_DISC = false; // used by Ken Stphenson's circle packing
bool INCLUDE_BOUNDARY_VERTICES = false;
bool IMMERSION_CROSSINGS_AT_RIGHT_ANGLES = false;
bool USE_BADNESS_OPTIMZATION = true;
bool USE_FORCE_DIRECTED_PLACEMENT = false;
bool APPLY_FORCES_TO_TYPE_12_ONLY = true;  // only applies when not doing Plestenjak force direction
bool USE_CENTRE_OF_GRAVITY_PLACEMENT = false;
bool TRACK_PLACEMENT_ITERATION = false;
bool USE_KEN_STEPHENSON_CIRCLE_PACKING = true;
//bool mp_control.magnify_small_circles = false;
bool CHECK_INNER_HULL_CALCULATION = false;
bool RETRACT_BOUNDARY_VERTICES = false;
bool PLESTENJAK_FORCE_DIRECTION = true;
bool USE_REGION_SHRINKING_PLACEMENT = false;
bool USE_SMALL_REGION_SHRINKING_PLACEMENT = false;
bool USE_EDGE_DISTRIBUTION_PLACEMENT = false;
bool PLOT_TRIANGULATION_EDGES = false;
bool USE_FIRST_GAP_AS_ACTIVE = false;
//bool mp_control.state_smoothed = false;

bool USE_ALL_LONG_GAPS = false;  // temporary - remove after testing

bool CONVEX_DISC = false; // draw a straight line triangulation of a disc rather than a knot diagram
bool LACES = false;

bool TEST_MODE = false; 
bool VALID_BOUNDARY_COUNTS = false;  // used for testing
bool SHOW_BODY_ARC_END_LABELS = false;
bool ALLOW_HEAD_HOOKS = false;

int check_inner_hull_vertex = 0;

int figure_count = 1; // used to record fignum in write_metapost
int max_circle_packing_iterations = 1000;
int max_placement_iterations = 200;
int placement_iteration_tracking_step = 1; // 100;
float average_triangulation_length_factor = 1.0;
float magnification_factor = 1.0;
float boundary_vertex_retraction_factor = 1.0;
float metapost_path_tension = 1.0;  

/* the amount by which the maximal region's area must exceed the average area of
   the compact regions 
*/
float region_shrinking_area_factor = 1.618034;

/* shrink area triangulations by this much */
float region_shrinking_factor = 0.75;  

/* plot_steps is the number of divisions within the histogram of edge lengths when plotting a triangulation */
int plot_steps = 20;

/* the amount by which vertices are moved towards the relevant COG in edge distribution placement */
float edge_distribution_shift_factor = 0.5;

/* triangulation_plot_threshold is the threshold dividing long and short edges determined by plotting a triangulation */
double triangulation_plot_threshold = 0;

double average_triangulation_edge_length = 0;

char const* circlepack_output_file =  "__circlepack.out";
char const* circlepack_output_savefile =  "__circlepack.sav";
char const* triangulation_output_file = "__triangulate.out";
char const* triangulation_shrink_file = "__shrink.out";

	
/********************* Function prototypes ***********************/
void dual_graph_placement(metapost_control& mp_control, generic_code_data& code_data, string title);
void force_directed_placement(metapost_control& mp_control, generic_code_data& code_data, int num_iterations, string title);				
void centre_of_gravity_placement(metapost_control& mp_control, generic_code_data& code_data, int num_iterations, string title);
void region_shrinking_placement(metapost_control& mp_control, generic_code_data& code_data, matrix<int>& cycle, int num_cycles, 
                                int num_left_cycles, int infinite_region, int num_iterations, string title);
void edge_distribution_placement(metapost_control& mp_control, generic_code_data& code_data, matrix<int>& cycle, int num_cycles, 
                                int num_left_cycles, int infinite_region, int num_iterations, string title);
void plot_triangulation_distribution(int num_steps);
bool get_next_input_string(string input_file,string& input_string, string& title);
bool connected_code_data(generic_code_data& code_data);
void help_info(bool exit_after_help);
void debug_help();
void set_programme_long_option(string option, string source, metapost_control& mp_control);
void set_programme_short_option(char* cptr, metapost_control& mp_control);
bool debug_setup(char* argv_ptr);
void check_debug_option_parameters(char* pptr, string option);
void read_code_data (generic_code_data& code_data, string input_string);
void triangulate (generic_code_data& code_data, matrix<int>& cycle, int num_cycles, int num_left_cycles, int infinite_region);
double badness(string vertex_file, generic_code_data& code_data);
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, matrix<double> vcoords, vector<double> vertex_radius, vector<int>* auxiliary_data=0);
void write_metapost(ofstream& os, generic_code_data& code_data, string title, metapost_control& mp_control, 
                    char const* circlepack_output_file, char const* circlepack_output_savefile, vector<int>* auxiliary_data=0);
bool calculate_turning_cycles(generic_code_data& code_data, matrix<int>& cycle, int& num_left_cycles, int& num_cycles);
bool valid_knotoid_input(generic_code_data& code_data);
void magnify(generic_code_data& code_data);
void set_frame_corners(string coordinate_file, metapost_control& mp_control);
vector<int> gauss_parity(generic_code_data& code_data);

int circle_pack (char const* inputfile, char const* outputfile);
int KS_circle_pack (char const* inputfile, char const* outputfile);

void draw_convex_triangulation (metapost_control& mp_control, const char* filename);
void draw_lace (metapost_control& mp_control, string input_string);
bool gauss_to_peer_code(generic_code_data gauss_code_data, generic_code_data& peer_code_data, bool optimal=true, vector<int>* gauss_crossing_map=0, bool evaluate_gauss_crossing_perm=false);
vector<int> identify_gauss_crossings(generic_code_data& code_data);
list<vector<int> > hamiltonian_circuit(generic_code_data& code_data, bool list_all_circuits, bool count_circuits_only, bool edge_circuit, int include_edge);

void sigfpe_handler (int sig) 
{
	cout << "Error, draw programme received SIGFPE!" << endl;
	exit(1);
}


/* scale factor applied to circle_pack output coordinates */
int metapost_coordinate_scale_factor=0;

/* additional scaling when drawing in a representation of the hyperbolic disc */
int metapost_hyperbolic_scale_factor = 1; 

/* badness controls */

/* If the minimum distance between successive vertices in the vertex sequence
   is less than badness_threshold, we look for an alternative infinite turning
   cycle, unless a specific infinite cycle has been provided.
*/
double badness_threshold = 100; 

/* candidate turning cycles must have minimum_infinite_cycle_length to 
   be considered as an infinite cycle.  We cannot set this to the value of 2
   since then the two type 4 vertices have two edges joining each other, which 
   prevents the circle packing algorithm from identifying boundary vertices 
   correctly.
*/
int minimum_infinite_cycle_length = 3;




int dual_graph_index = 1;
namespace util
{
	string itos(long n);
}

/******************* Main Function ************************/
int main (int argc, char* argv[])
{
    string 	input_file;
    string 	output_file;
	string 	input_string;
	string 	title;
    bool    input_file_provided = false;
    bool    output_file_provided = false;

	/* default_metapost_control is initialized with the default settings */
    metapost_control default_metapost_control;

	signal(SIGFPE, sigfpe_handler); // register the SIGFPE handler

    /* Determine command line options and whether input is to
       come from the screen or a file */
    if (argc == 1)
	{
		// help_info(true);
	}
	else
    {
		for (int i=1; i < argc; i++)
		{
			if (*(argv[i]) == '-')
			{
				string test(argv[i]);
				if (*(argv[i]+1) == '-')
					set_programme_long_option(argv[i], "command line",default_metapost_control);
				else 
					set_programme_short_option(argv[i], default_metapost_control);

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw::main: after argument " << argv[i] << " default_metapost_control:" << endl;
    print (default_metapost_control, debug, "main:   ");
}
					
			}
			else if (!input_file.length()) 
			{
				input_file = argv[i];
	    		input_file_provided = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw::main: input_file " << input_file << " provided" << endl;
			}
			else if (!output_file.length()) 
			{
				output_file = argv[i];
		    	output_file_provided = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw::main: ouput_file " << output_file << " provided" << endl;
			}
			else
			{
		    	cout << "Usage draw --<task> [-<options>][<infile>[<outfile>]], type draw -h for help.\n";
		    	exit(0);
			}
		}	
    }

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw::main: command line =  ";
	for (int i = 0; i< argc; i++)
		debug << argv[i] << ' ';
	debug << endl;		
}
//	print_prog_params(debug, "detail");

	/* If the metapost_coordinate_scale_factor has not been set from 
	   the command line, set it to the default
	*/
	if (!metapost_coordinate_scale_factor)
	{
		if (CONVEX_DISC)
			metapost_coordinate_scale_factor  = 600;
		else if (USE_KEN_STEPHENSON_CIRCLE_PACKING)
			metapost_coordinate_scale_factor  = 2500;
		else if (USE_FORCE_DIRECTED_PLACEMENT)
		{
			if (PLESTENJAK_FORCE_DIRECTION)
				metapost_coordinate_scale_factor  = 1000;
			else
				metapost_coordinate_scale_factor  = 100;
			
		}
		else
			metapost_coordinate_scale_factor  = 25;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw::main: metapost_coordinate_scale_factor set to default value of " << metapost_coordinate_scale_factor << "\n";
	}
	
	
	if (input_file_provided)
	{
    	input.open(input_file.c_str());
    	if(!input)
    	{
			cout << "\nError opening input file\n";
			exit(0);
    	}	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw::main: input file: " << input_file << "\n";

		/* Read any options included in the input file */
		string next_line;
		while(getline(input,next_line))
		{
			string::size_type pos1 = next_line.find(';');
			if (pos1 != string::npos)
				next_line.erase(pos1);
				
			if (next_line.find('[')!=string::npos && next_line.find('\\')==string::npos && next_line.find('/')==string::npos && next_line.find('X')==string::npos)
			{
				/* this is an options line not a line containing a peer code
				   make sure there's an end to the option definitions on this line 
				*/
				if (next_line.find(']')==string::npos)
				{
					cout << "\nError in " << input_file << ": programme options must be enclosed within []\n";
					exit(0);
				}
							
				pos1 = next_line.find('[')+1;
				string::size_type pos2;
				
				if (next_line.find(',') != string::npos)
				{
					do
					{
						pos2 = next_line.find(',',pos1);
						set_programme_long_option(next_line.substr(pos1,pos2-pos1),"input file",default_metapost_control);
						pos1 = pos2+1;
					} while (next_line.find(',',pos1) !=string::npos);
				}
				pos2 = next_line.find(']',pos1);
				set_programme_long_option(next_line.substr(pos1,pos2-pos1),"input file",default_metapost_control);
			}			
    	}
	}

    /* prepare output file */
	if (!output_file_provided)
		output_file = "draw.out";
	
    output.open (output_file.c_str());

    if (!output)
    {
		cout << "\nError opening output file " << output_file << endl;
		exit(0);
    }
    else
    {
		output << "% Output from draw v" << version << "\n";

		if (input_file_provided)
		{
			output << "% Input file: " << input_file << "\n";
		}
		
		output << "% Command line: ";
		for (int i = 0; i< argc; i++)
			output << argv[i] << ' ';
		output << endl;		
		
		output << "\ninput boxes\n" << endl;
		output << "\ninput mpcolornames\n" << endl;
		output << "\nprologues:=3;\n" << endl;		

		output << "numeric fignum; fignum=1;" << endl;
    }


    if (!input_file_provided)
    {
	
//		if (argc > 1)
		{
			cout << "\nThis is A.Bartholomew's draw programme, v" << version << endl;
		}
		
		if (LACES)
		{
			cout << "\nThe programme will draw the lace determined by a lace code\n";
			
			cout <<"\nType help at the input prompt to view the help screens, type q to exit.\n";
	
			cout << "\nEnter a lace code.";		
		}
		else if (CONVEX_DISC)
		{
			cout << "\nThe programme will draw the convex, straight line triangulation of a disc\n";
			cout << "determined by either a labelled peer code or an abstract triangulation\n";
			cout << "description contained in a file\n";
			
			cout <<"\nType help at the input prompt to view the help screens, type q to exit.\n";
	
			cout << "\nEnter a labelled peer code or filename.";		
		}
		else
		{
			cout << "\nThe programme will produce a metapost file for drawing the knot or knotoid";
			cout << " specified by a labelled peer code, Gauss code or planar diagram.\n";
			
			cout <<"\nType help at the input prompt to view the help screens, type q to exit.\n";
	
			cout << "\nEnter diagram code.";		
		}
    }
	
	if (LACES)
	{
		input_control::ACCEPT_MAP |= input_control::lace_code;
	}
	else
	{
		input_control::ACCEPT_MAP |= input_control::peer_code;
		input_control::ACCEPT_MAP |= input_control::gauss_code;
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw::main: input_control::ACCEPT_MAP set to 0x" << hex << input_control::ACCEPT_MAP << dec << endl;
	
try {

    if (input_file_provided)
    {
		/* reset input ready for read */
		input.clear(); // state flags
		input.seekg(0);
	}

	while (get_next_input_string(input_file,input_string, title))
	{		    
		metapost_control mp_control(default_metapost_control);
	
		/* the input string might be a peer code or a the name of a file containing a description of 
		   an arbitrary disc triangulation.  If it is a peer code we read it and create the standard
		   triangulation determined by that peer code before deciding how to proceed
		*/
		bool peer_code_input = false;
		bool gauss_code_input = false;
		bool planar_diagram_input = false;
		
		generic_code_data code_data;
		vector<int> gauss_crossing_map;  // used for drawing specific smoothed states of Gauss codes and adding Gauss crossing labels.
		int num_cycles = 0;
		int num_left_cycles;
		int num_crossings;
		int num_edges;
		int infinite_region;
		matrix<int> cycle(0,0);
		
		if (input_string.find('X') != string::npos)
		{
				planar_diagram_input = true;		
				gauss_code_input = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: identified input_string as a planar diagram" << endl;
		}
		else if (input_string.find('[') != string::npos)
		{		
			peer_code_input = true;		
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: identified input_string as a peer code" << endl;
		}
		else if (input_string.find('(') != string::npos && input_string.find("K(") == string::npos) // indicates input is a lace
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: identified input_string " << input_string << " as a lace code" << endl;
		}
		else
		{
			gauss_code_input = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: identified input_string " << input_string << " as a Gauss code" << endl;				
		}
		
		if (peer_code_input || gauss_code_input)
		{
			/* first remove any unwanted qualifiers from the input string */
			string::size_type pos1 = input_string.find('{');
			string::size_type pos2;
			if (pos1 != string::npos)
			{
				pos1 += 1;
				
				if (input_string.find(',',pos1) != string::npos)
				{
					do
					{
						pos2 = input_string.find(',',pos1);
						set_programme_long_option(input_string.substr(pos1,pos2-pos1),"qualifiers",mp_control);
						pos1 = pos2+1;
					} while (input_string.find(',',pos1) !=string::npos);
				}
				pos2 = input_string.find('}',pos1);
				set_programme_long_option(input_string.substr(pos1,pos2-pos1),"qualifiers",mp_control);
				input_string.erase(input_string.find('{'));
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "draw::main: metapost_control for input_string:" << endl;
    print(mp_control, debug, "draw::main:   ");
}
			if (peer_code_input)
			{
				read_peer_code(code_data,input_string);		
				
				if (mp_control.gauss_crossings)
					gauss_crossing_map = identify_gauss_crossings(code_data);
					
			}
			else
			{	
				generic_code_data gauss_code_data;
				
				if(planar_diagram_input)
					read_planar_diagram(gauss_code_data,input_string);		
				else
					read_gauss_code(gauss_code_data,input_string);		
					
				gauss_crossing_map = vector<int>(gauss_code_data.num_crossings);

				/* If we're given a Gauss code that does not start at crossing 1, then gauss_arc 0 does not terminate at crossing 1.
				   The following call to gauss_to_peer_code evaluates the gauss_crossing_map, which tells us the order in which we first
				   encounter the classical immersion crossings as we trace the diagram from gauss_arc zero.
				   gauss_crossing_map[i] is the immersion crossing corresponding to Gauss crossing i
				*/
				gauss_to_peer_code(gauss_code_data, code_data, true, &gauss_crossing_map, true); // optimal=true, evaluate_gauss_crossing_perm = true
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "draw::main: gauss_to_peer_code returned peer code: ";
    write_peer_code(debug,code_data);
    debug<< endl;
    debug << "draw::main: gauss_to_peer_code returned gauss_crossing_map: ";
    for (int i = 0; i < gauss_code_data.num_crossings; i++)
		debug << gauss_crossing_map[i] << ' ';
	debug << endl;
}
				
			}
			
			if (!connected_code_data(code_data))
			{
				
				cout << "Input code does not describe a connected diagram" << endl;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: invalid input detected, code does not describe a connected diagram" << endl;
	
				continue;
			}
			else
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: input code describes a connected diagram" << endl;
			}
	
		
			/* There cannot be more than num_crossings + 2 cycles, and they cannot 
			   be more than num_edges in length, we add 1 to the length of a 
		   	   cycle as we store the actual length in column 0 */
			num_crossings = code_data.num_crossings;
			num_edges = 2*num_crossings;
//			matrix<int> cycle(num_crossings+2, num_edges+1);
			cycle = matrix<int>(num_crossings+2, num_edges+1);
			calculate_turning_cycles(code_data, cycle, num_left_cycles, num_cycles);
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw::main: left turning cycles:" << endl;   
	for (int i=0; i<num_left_cycles; i++)
	{
		debug << "draw::main:   cycle " << i << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " ";
		debug << endl;
	}

	debug << "draw::main: right turning cycles:" << endl;   
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
		debug << "draw::main:   cycle " << i << ": ";
		for (int j=1; j<=cycle[i][0]; j++)
			debug << cycle[i][j] << " "	;
		debug << endl;
	}
}
	
			if (code_data.head != -1)
			{
				if (!valid_knotoid_input(code_data))
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: invalid knotoid input detected" << endl;
					
					cout << "Input is not a valid knotoid, indicated shortcut does not pass under the rest of the diagram" << endl;				
					continue;			
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: valid knotoid input detected" << endl;
				}
			}
	
//			int infinite_region = -1;
			infinite_region = -1;
			if (mp_control.knotoid_leg_unbounded)
			{
				/* Select the initial choice of infinite region, this will be 
				   one bounded of maximal length that includes edge zero
				*/
				for (int i=1; i< num_cycles; i++)
				{
					int length = cycle[i][0];
					
					bool contains_zero = false;
					
					for (int j=1; j<= length; j++)
					{
						if (cycle[i][j] == 0)
						{
							contains_zero = true;
							break;
						}
					}
					
					if (contains_zero)
					{
						if (infinite_region != -1)
						{
							if (length > cycle[infinite_region][0])
								infinite_region = i;
						}
						else
							infinite_region = i;
					}
				}
			}
			else if (mp_control.infinite_cycle != -1)
			{
				infinite_region = mp_control.infinite_cycle;
			}
			else
			{
				infinite_region = 0;
				
				/* Select the initial choice of infinite region, this will be 
				   one bounded by a turning cycle of maximal length
				*/
				for (int i=1; i< num_cycles; i++)
				{
					if (cycle[i][0] > cycle[infinite_region][0])
						infinite_region = i;
				}
			}
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: initial infinite_region = " << infinite_region << endl;
					
	
			triangulate(code_data, cycle, num_cycles, num_left_cycles, infinite_region);
		}

		if (LACES)
		{
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "main: laces option provided with input string: " << input_string << endl;
    
			/* first remove any unwanted qualifiers from the input string */
			string::size_type pos1 = input_string.find('{');
			string::size_type pos2;
			if (pos1 != string::npos)
			{
				pos1 += 1;
				
				if (input_string.find(',',pos1) != string::npos) 
				{
					do
					{
						pos2 = input_string.find(',',pos1);
						set_programme_long_option(input_string.substr(pos1,pos2-pos1),"qualifiers",mp_control);
						pos1 = pos2+1;
					} while (input_string.find(',',pos1) !=string::npos);
				}
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "main: after parsing qalifiers up to last comma, we're left with string:" << input_string.substr(pos1) << endl;
}				
				pos2 = input_string.find('}',pos1);
				set_programme_long_option(input_string.substr(pos1,pos2-pos1),"qualifiers",mp_control);
				input_string.erase(input_string.find('{'));
				
			}
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "main: metapost_control for input_string:" << endl;
    print(mp_control, debug, "main:   ");
}
			
			draw_lace(mp_control,input_string);		
			
		}
		else if (CONVEX_DISC)
		{
			const char* filename;
			
			if (peer_code_input)
			{
				filename = triangulation_output_file;
			}
			else
			{
				filename = c_string(input_string);
			}
			
			draw_convex_triangulation (mp_control,filename);
			
			if (!peer_code_input)
				delete [] filename;
		}
		else
		{
			if (USE_FORCE_DIRECTED_PLACEMENT)
			{
				force_directed_placement(mp_control, code_data, max_placement_iterations, title);
			}
			else if (USE_CENTRE_OF_GRAVITY_PLACEMENT)
			{
				centre_of_gravity_placement(mp_control, code_data, max_placement_iterations, title);
			}			
			else if (USE_REGION_SHRINKING_PLACEMENT)
			{
				region_shrinking_placement(mp_control, code_data, cycle, num_cycles, num_left_cycles, infinite_region, max_placement_iterations, title);
			}
			else if (USE_EDGE_DISTRIBUTION_PLACEMENT)
			{
				edge_distribution_placement(mp_control, code_data, cycle, num_cycles, num_left_cycles, infinite_region, max_placement_iterations, title);
			}
			else
			{			
				if (USE_KEN_STEPHENSON_CIRCLE_PACKING)
					KS_circle_pack (triangulation_output_file, circlepack_output_file);
				else
					circle_pack (triangulation_output_file, circlepack_output_file);


				/* Exploratory code added January 2019 to plot the distribution of triangulation edge lengths after circle packing */
				if (PLOT_TRIANGULATION_EDGES)
				{
					plot_triangulation_distribution(plot_steps);
				}
				/* End of exploratory code added January 2019 to plot the distribution of triangulation edge lengths after circle packing */
				
				double badness_rating = badness(circlepack_output_file, code_data);
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: badness_rating when using cycle " << infinite_region << " as infinite_region = " << badness_rating << endl;
	
	
				if (USE_BADNESS_OPTIMZATION && !mp_control.knotoid_leg_unbounded && badness_rating < badness_threshold && mp_control.infinite_cycle == -1)
				{
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: badness_rating less than " << badness_threshold << ", looking for other candidates" << endl;
	
					int new_infinite_region = infinite_region;
					
					for (int i=0; i< num_cycles; i++)
					{
						if (i != infinite_region && cycle[i][0] >= minimum_infinite_cycle_length)
						{
							triangulate(code_data, cycle, num_cycles, num_left_cycles,i);
							if (USE_KEN_STEPHENSON_CIRCLE_PACKING)
								KS_circle_pack (triangulation_output_file, circlepack_output_file);
							else
								circle_pack (triangulation_output_file, circlepack_output_file);
										
							double trial_badness = badness(circlepack_output_file, code_data);
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main:   badness_rating when using cycle " << i << " as infinite_region = " << badness_rating << endl;
	
							if (trial_badness > badness_rating)
							{
								badness_rating = trial_badness;
								new_infinite_region = i;
							}
						}
						else
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main:   cycle " << i << " not a viable alternative candidate as infinite_region" << endl;
						}
					}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: best infinite_region selection is cycle " << new_infinite_region << endl;
    
					mp_control.infinite_cycle = new_infinite_region;
	
					triangulate(code_data, cycle, num_cycles, num_left_cycles,new_infinite_region);
					
					if (USE_KEN_STEPHENSON_CIRCLE_PACKING)
						KS_circle_pack (triangulation_output_file, circlepack_output_file);
					else
						circle_pack (triangulation_output_file, circlepack_output_file);
				}
				else
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: badness_rating acceptable" << endl;

					mp_control.infinite_cycle = infinite_region;

				}
			}
	
			/* save a copy of the circle packing output for drawing comparative triangulations when magnifying 
			   the write_metapost needs to see this file for all drawing variants.
			*/
			string command = "cp "+string(circlepack_output_file)+" "+string(circlepack_output_savefile);
			system(command.c_str());
		
			/* For the circle packing drawing variants we suport circle magnification */
			if (mp_control.magnify_small_circles && !USE_FORCE_DIRECTED_PLACEMENT && !USE_CENTRE_OF_GRAVITY_PLACEMENT)
				magnify(code_data);
	
			/* Write metapost of the final diagram */
		
			/* record the frame corners based on the final choice of infinite cycle */
			set_frame_corners(circlepack_output_file,mp_control); 
	
			/* save then disable control setting that we don't want in the final diagram */
			bool saved_draw_shrink_effect = mp_control.draw_shrink_effect;
			mp_control.draw_shrink_effect = false;
			
			bool saved_draw_triangulation_displacement = mp_control.draw_triangulation_displacement;
			
			if (TRACK_PLACEMENT_ITERATION && mp_control.draw_triangulation_displacement)
				mp_control.draw_triangulation_displacement = false;
			
			bool saved_highlight_small_edges = mp_control.highlight_small_edges;
			if (TRACK_PLACEMENT_ITERATION && mp_control.highlight_small_edges)
				mp_control.highlight_small_edges = false;
			
			if (mp_control.show_odd_parity_crossings)
			{
				vector<int> crossing_parity = gauss_parity(code_data);

if (debug_control::DEBUG >= debug_control::BASIC)
{
//	debug << "parity: num_parity_terms = " << num_parity_terms << " crossing_parity: ";
	debug << "crossing parity: ";
//	for (int i=0; i< num_parity_terms; i++)
	for (int i=0; i< num_crossings; i++)
	{
		if (crossing_parity[i] == gauss_orientation_data::parity::ODD)		
			debug << "O ";
		else if (crossing_parity[i] == gauss_orientation_data::parity::EVEN)		
			debug << "E ";
		else
			debug << "N ";
	}
	debug << endl;
}
				for (int i=0; i< num_crossings; i++)
				{
					if (crossing_parity[i] == gauss_orientation_data::parity::ODD)		
						code_data.code_table[generic_code_data::table::LABEL][i] = generic_code_data::ODD;
				}

			}

			if (mp_control.state_smoothed)
			{
//				bool draw_labels = mp_control.draw_labels;
//				bool draw_shortcut = mp_control.draw_shortcut;
//				mp_control.draw_labels = true;
//				mp_control.draw_shortcut = true;
				if (!mp_control.seifert_circles && mp_control.state.length() == 0)  // drawing all states not just one specific state
					write_metapost(output, code_data, title, mp_control, circlepack_output_file, circlepack_output_savefile);	
//				mp_control.draw_labels = draw_labels;
//				mp_control.draw_shortcut = draw_shortcut;

				unsigned int num_state_crossings = num_crossings;
				
				for (int i=0; i< num_crossings; i++)
				{
					if (code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::VIRTUAL)
						num_state_crossings--;
				}
				
				if (code_data.head != -1 && code_data.shortcut_crossing.size() !=0)
				{
					for (unsigned int i=0; i< code_data.shortcut_crossing.size(); i++)
					{
						if (code_data.shortcut_crossing[i])
							num_state_crossings--;
					}
				}

				if (mp_control.show_odd_parity_crossings)
				{
					for (int i=0; i< num_crossings; i++)
					{
						// enum parity {NONE = 0, ODD = 1, EVEN = 2}
						if (code_data.code_table[generic_code_data::table::LABEL][i] == generic_code_data::ODD)
							num_state_crossings--;
					}
				}
				
				
				bool finished=false;
				vector<int> state(num_state_crossings);
				for (unsigned int i=0; i< num_state_crossings; i++)
					state[i] = metapost_control::smoothed::NON_SEIFERT_SMOOTHED;
				
				/* To draw all smoothings we set the state vector to 1 to indicate a Seifert smoothed crossing and -1 to indicate a non-Seifert-smoothed crossing.
				   Note that the A and B smoothing is NOT the same as Seifert smoothing; the two are related dependent on the sign of the crossing.  If we are given
				   an explicit state to draw, we will be given an A or B indication for each Gauss crossing, regardless of the sign of the crossing and regardless of 
				   the order in which the crossings appeared in the input Gauss code.
				   
				   To indicate that we wish to force a crossing to be A smoothed or B smoothed, we use the following state vector assignments
						
						1 = Seifert smoothed
					   -1 = non-Seifert smoothed
					    2 = A-smoothed
					    3 = B=smoothed
				 
				*/
				if  (mp_control.seifert_circles || mp_control.state.length() != 0)
				{

					if (!mp_control.seifert_circles && mp_control.state.length() != num_state_crossings)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "write_metapost: Error! Incorrect number of states provided for ";
	write_code_data(debug, code_data);
	debug << "\nwrite_metapost: expected " << num_state_crossings << " states, received state data " << mp_control.state << endl;
}
						cout << "Error! Incorrect number of states provided for " << input_string << endl;
						cout << "Expected " << num_state_crossings << " states, received state data " << mp_control.state << endl;
						exit(0);
					}
					
					if (mp_control.seifert_circles)
					{
						for (unsigned int i=0; i< num_state_crossings; i++)
							state[i] = metapost_control::smoothed::SEIFERT_SMOOTHED;
					}
					else if (peer_code_input)
					{
						int place = 0;
						for (unsigned int i=0; i< mp_control.state.length(); i++)
						{
							if (mp_control.state[i] == '*' || mp_control.state[i] == 'S' || mp_control.state[i] == 'O')
								continue;
								
							state[place] = (mp_control.state[i] == 'A'? metapost_control::smoothed::A_SMOOTHED: metapost_control::smoothed::B_SMOOTHED);
							place++;
						}
					}
					else if (gauss_code_input)
					{
						/* The state data is given in the order corresponding to the Gauss crossings but we may not encounter the crossings
						   sequentially as we trace a diagram dscribed by a Guass code from gauss_arc zero.  The vector gauss_crossing_map
						   tells us the order in which we encountered the Gauss crossings in the immersion:  gauss_crossing_map[i] being the 
						   immersion crossing corresponding to Gauss crossing i
						   
						   We need to pass write_metapost a state vector in the order we encounter crossings in the immersion, so we create 
						   a sorted copy of gauss_crossing_map that specifies the order in which we need to pass the state to write_metapost.
						   
						   If were are considering parity, we need to take out of the gauss_crossing_map those crossings that are ODD, since
						   they will not be smoothed and therefore we will have no state designation for them.
						   
						*/
						if (mp_control.show_odd_parity_crossings)
						{
							vector<int> gauss_parity_map(num_state_crossings);
							int place=0;
							
							for (unsigned int i=0; i< gauss_crossing_map.size(); i++)
							{
								if (code_data.code_table[generic_code_data::table::LABEL][gauss_crossing_map[i]] != generic_code_data::ODD)
									gauss_parity_map[place++] = gauss_crossing_map[i];
							}
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "draw::main: gauss_parity_map: ";
    for (unsigned int i = 0; i < num_state_crossings; i++)
		debug << gauss_parity_map[i] << ' ';
	debug << endl;
}					
							gauss_crossing_map = gauss_parity_map;
						}					
						
						vector<int> sorted_map = gauss_crossing_map;
						sort(sorted_map.begin(),sorted_map.end());
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "draw::main: sorted gauss_to_peer_code: ";
    for (unsigned int i = 0; i < sorted_map.size(); i++)
		debug << sorted_map[i] << ' ';
	debug << endl;
}					
						for (unsigned int i=0; i< mp_control.state.length(); i++)
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: mp_control.state[" << i << "] = " << mp_control.state[i] << " gauss_crossing_map[" << i << "] = " << gauss_crossing_map[i];
    
							/* find gauss_crossing_map[i] in sorted_map */
							int place = 0;
							for (unsigned int j=0; j < sorted_map.size(); j++)
							{
								if (sorted_map[j] == gauss_crossing_map[i])
								{
									place = j;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << " found at place " << j << " of sorted gauss_crossing_map" << endl;
									break;
								}
							}
							
							state[place] = (mp_control.state[i] == 'A'? metapost_control::smoothed::A_SMOOTHED: metapost_control::smoothed::B_SMOOTHED);
						}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "draw::main: drawing single state: ";
    for (unsigned int i=0; i< num_state_crossings; i++)
		debug << state[i] << ' ';
	debug << endl;
}
					}
					else
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "draw::main: state data provided when input is neither a peer code or a Gauss code" << endl;						
					}
				}
				else
				{				
					cout << "number of smoothed states = " << pow(2,num_state_crossings) << endl;
				}
				
//				bool draw_labels = mp_control.draw_labels;
//				bool draw_shortcut = mp_control.draw_shortcut;
//				mp_control.draw_labels = false;
//				mp_control.draw_shortcut = false;
				do
				{

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "draw::main: drawing state: ";
    for (unsigned int i=0; i< num_state_crossings; i++)
		debug << state[i] << ' ';
	debug << endl;
}

					write_metapost(output, code_data, title, mp_control, circlepack_output_file, circlepack_output_savefile,&state);	
						
					/* increment state */
					int place = num_state_crossings-1;
					if (!mp_control.seifert_circles && num_state_crossings > 0 && mp_control.state.length() == 0)
					{
						do
						{
							if (state[place] == metapost_control::smoothed::NON_SEIFERT_SMOOTHED)
							{
								state[place] = metapost_control::smoothed::SEIFERT_SMOOTHED;
								for (unsigned int i = place+1; i< num_state_crossings; i++)
									state[i] = metapost_control::smoothed::NON_SEIFERT_SMOOTHED;
								break;
							}
							else if (place == 0)
								finished = true;
							else
								place--;
						} while (!finished);
					}
					else
					{
						finished = true;
					}
							
				} while (!finished);
//				mp_control.draw_labels = draw_labels;
//				mp_control.draw_shortcut = draw_shortcut;
				
			}
			else if (mp_control.hamiltonians)
			{
				list<vector<int> > circuit_list;
				
				if (mp_control.hamiltonian_circuit.size() == 1)
				{
					circuit_list = hamiltonian_circuit(code_data, false, false, true, mp_control.HC_include_edge); //list_all_circuits = false, count_circuits_only = false, edge_circuit = true;
				}
				else
				{
					circuit_list = hamiltonian_circuit(code_data, true, false, true, mp_control.HC_include_edge); //list_all_circuits = true, count_circuits_only = false, edge_circuit = true;
				}
				
				if (mp_control.hamiltonian_circuit.size() > 1)
				{
					write_metapost(output, code_data, title, mp_control, circlepack_output_file, circlepack_output_savefile,&mp_control.hamiltonian_circuit);		
				}
				else
				{
					list<vector<int> >::iterator lptr = circuit_list.begin();
					
					while (lptr != circuit_list.end())
					{
						write_metapost(output, code_data, title, mp_control, circlepack_output_file, circlepack_output_savefile,&(*lptr));		
						lptr++;
					}
					
					if (circuit_list.size() == 0)
						cout << "No Hamiltonian circuit found" << endl;
					else
						cout << "Written " << circuit_list.size() << " Hamiltonian " << (circuit_list.size()>1? "circuits": "circuit") << endl;
				}
			}
			else
			{
				if (peer_code_input && !mp_control.gauss_crossings)
					write_metapost(output, code_data, title, mp_control, circlepack_output_file, circlepack_output_savefile);							
				else
					write_metapost(output, code_data, title, mp_control, circlepack_output_file, circlepack_output_savefile,&gauss_crossing_map);		
			}
			
			mp_control.draw_shrink_effect = saved_draw_shrink_effect;
			mp_control.draw_triangulation_displacement = saved_draw_triangulation_displacement;			
			mp_control.highlight_small_edges = saved_highlight_small_edges;
		}
	}
}
catch (const bad_alloc& __e) 
{
	cout << "ERROR! Out of memory, bad_alloc thown: " << __e.what() << endl;
}
catch (...)
{
	cout << "ERROR! Exception thrown that is not caught explicitly" << endl;
}


	output << "\nend" << endl;

	if (input_file_provided)
		input.close();
	
    output.flush();
    output.close();


if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug.close();

	return 0;
} /* End of main */

void set_programme_long_option(string option, string source, metapost_control& mp_control)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: provided with option " << option << " by " << source << endl;
	
	char* cptr = c_string(option); 

	char  loc_buf[40];
	
	char* c1;
	char* c2;

	c1 = cptr;

	/* take out any leading space */
	while (*c1 == ' ')
		c1++;
	
	/* take out any leading -- */
	while (*c1 == '-')
		c1++;
	
	c2 = loc_buf;
	while (isalpha(*c1) || isdigit(*c1) || *c1 == '-' || *c1 == '_')
		*c2++ = *c1++;
    *c2 = '\0';

	if (!strcmp(loc_buf,"centre"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: rotation centre read from " << source << endl;

		bool error = false;
		
		if (*c1 == '=')
		{
			c1++;
			
			if (*c1 == 'z')
			{				
				get_number(mp_control.rotation_centre_z,++c1);
				mp_control.implicit_rotation_centre = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: rotation_centre_z read from " << source 
	      << ", rotation_centre_z = " << mp_control.rotation_centre_z << endl;
}	    
			}
			else
			{
				if (*c1 == '(')
				{
					get_number(mp_control.rotation_centre_x,++c1);
					while (isdigit(*c1))
						c1++;
						
					if (*c1 == ',')
					{
						get_number(mp_control.rotation_centre_y,++c1);
						
						while (isdigit(*c1))
							c1++;
						
						if (*c1 == ')')
						{
							mp_control.explicit_rotation_centre = true;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: explicit rotation centre (" << mp_control.rotation_centre_x 
	      << "," << mp_control.rotation_centre_y << ")" << endl;
}	    
						}
						else
							error=true;
					}
					else
						error=true;
				}
				else
					error=true;				
			}
		}
		else
			error = true;
			
		if (error)
		{
			cout << "\nYou must specify the required rotation centre if you use the centre option." << endl;
			cout << "Specify the coordinates of the centre explicitly, e.g. centre=(123.45,67.89), or implicitly, e.g. centre=z21." << endl;
			exit(0);
		}
				
	}
	else if (!strcmp(loc_buf,"colour"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting mp_control.colour true as a result of " << source << " option" << endl;

		mp_control.colour = true;
		if (*c1 == '=')
		{
			mp_control.singlecolour = ++c1;
		}
	}
	else if (!strcmp(loc_buf,"colour-map"))
	{
		if (*c1 == '=')
		{
			/* move over = sign and copy the filename to the start of loc_buf 
			++c1; 
			c2 = loc_buf;
			while (*c1 != '\0')
				*c2++ = *c1++;
		    *c2 = '\0';
			mp_control.colourmap = loc_buf;*/
			mp_control.colourmap = ++c1;
		}
		else
		{
			cout << "\nYou must specify a colour map file if you use the colour-map option, e.g. 'colour-map=mycolours'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: mp_control.colourmap read from " << source << ", mp_control.colourmap = " 
	      << mp_control.colourmap << endl;
}	
	}	
	else if (!strcmp(loc_buf,"convex-disc"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting CONVEX_DISC true as a result of " << source << " option" << endl;

		CONVEX_DISC = true;
	}
	else if (!strcmp(loc_buf,"cudgel-space"))
	{
		if (*c1 == '=')
		{
			get_number(mp_control.cudgel_space,++c1);
		}
		else
		{
			cout << "\nYou must specify a cudgel spacing factor if you use the cudgel-space option, e.g. 'cudgel-space=1.5'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: cudgel_space read from " << source << ", cudgel_space = " 
	      << mp_control.cudgel_space << endl;
}	
	}
	else if (!strcmp(loc_buf,"cusp-disc-size"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting cusp_disc_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.cusp_disc_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a disc size multiplier if you use the cusp-disc-size option, e.g. 'cusp-disc-size=4'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: cusp_disc_size read from " << source << ", cusp_disc_size = " 
	      << mp_control.cusp_disc_size << endl;
}	
	}
	else if (!strcmp(loc_buf,"cycle"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting infinite region turning cycle as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.infinite_cycle,++c1);
		}
		else
		{
			cout << "\nYou must specify a turning cycle if you use the cycle option, e.g. 'cycle=2'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: infinite_cycle read from " << source << ", infinite_cycle = " 
	      << mp_control.infinite_cycle << endl;
}	
	}
	else if (!strcmp(loc_buf,"dots"))
	{
    	mp_control.dash_with_dots = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: dash_with_dots read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"disc-size"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting disc_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.disc_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a unit size if you use the disc-size option, e.g. 'disc-size=20'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: disc_size read from " << source << ", disc_size = " 
	      << mp_control.disc_size << endl;
}	
	}
	else if (!strcmp(loc_buf,"draw-lace-frame"))
	{
		mp_control.draw_lace_frame = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set draw_lace_frame true" << endl;
	}
	else if (!strcmp(loc_buf,"edge"))
	{ 
		USE_EDGE_DISTRIBUTION_PLACEMENT = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: use edge distribution placement of triangulation vertices" << endl;

		if (*c1 == '=')
	    {
			get_number(max_placement_iterations,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set max_placement_iterations = " << max_placement_iterations << endl;

	    }	
	}
	else if (!strcmp(loc_buf,"edge-factor"))
	{
	    if (*c1 == '=')
	    {
			get_number(edge_distribution_shift_factor,++c1);
	    }
		else
		{
			cout << "\nYou must specify an edge_distribution_shift_factor if you use the edge-factor option, e.g. 'edge-factor=0.2'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: region_shrinking_factor = " << region_shrinking_factor << endl;
	}
	else if (!strcmp(loc_buf,"first-gap"))
	{ 
		USE_FIRST_GAP_AS_ACTIVE = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: USE_FIRST_GAP_AS_ACTIVE read from " << source << endl;

	}
	else if (!strcmp(loc_buf,"force"))
	{ 
		USE_FORCE_DIRECTED_PLACEMENT = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: use force directed placement of triangulation vertices" << endl;

		if (*c1 == '=')
	    {
			get_number(max_placement_iterations,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set max_placement_iterations = " << max_placement_iterations << endl;

	    }	
	}
	else if (!strcmp(loc_buf,"frame-corners"))
	{
    	mp_control.draw_frame_corners = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: frame-corners read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"gauss-crossings"))
	{
    	mp_control.gauss_crossings = true;
    	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: gauss_crossings read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"gauss-labels"))
	{
    	mp_control.gauss_labels = true;
    	mp_control.draw_labels = true;
    	mp_control.label_edges_from_one = true;
    	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: gauss_labels read from " << source << endl;
	debug << "set_programme_long_option: setting draw_labels as a result of draw_labels read from " << source << endl;
	debug << "set_programme_long_option: setting label_edges_from_one as a result of draw_labels read from " << source << endl;
}
	}
	else if (!strcmp(loc_buf,"gravity"))
	{ 
		USE_CENTRE_OF_GRAVITY_PLACEMENT = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: use centre of gravity placement of triangulation vertices" << endl;

		if (*c1 == '=')
	    {
			get_number(max_placement_iterations,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set max_placement_iterations = " << max_placement_iterations << endl;

	    }	
	}
	else if (!strcmp(loc_buf,"grid"))
	{

		mp_control.draw_grid = true;

		if (*c1 == '=')
		{
			get_number(mp_control.grid_size,++c1);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: mp_control.draw_grid read from " << source << ", grid size = " << mp_control.grid_size << endl;

	}
	else if (!strcmp(loc_buf,"hamiltonians"))
	{

		mp_control.hamiltonians = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: mp_control.hamiltonians read from " << source << endl;

		mp_control.draw_crossing_features = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set draw_crossing_features false" << endl;

	}
	else if (!strcmp(loc_buf,"hamiltonian-circuit"))
	{
		mp_control.hamiltonians = true;

		mp_control.draw_crossing_features = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set draw_crossing_features false" << endl;

		if (*c1 == '=')
		{			
			istringstream iss(++c1);
			int edge;
			while (iss >> edge)
				mp_control.hamiltonian_circuit.push_back(edge);

/*			
			char ch = ' ';
			iss >> ch; // the '{'
			do 
			{
				int edge;
				iss >> edge;
				mp_control.hamiltonian_circuit.push_back(edge);
				iss >> ch;
			} while (ch != '}');
*/
		}
		else
		{
			mp_control.hamiltonian_circuit = vector<int>(1);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: hamiltonian-circuit read from " << source;
	if (mp_control.hamiltonian_circuit.size() != 0)
	{
		debug << ", circuit = ";
		for (size_t i=0; i< mp_control.hamiltonian_circuit.size(); i++)
		 debug << mp_control.hamiltonian_circuit[i] << ' ';
	}
	debug << endl;
		
	debug << "set programme long_option: mp_control.hamiltonians set as a result of receiving hamiltonian_circuit option from " << source << endl;
}
	}
	else if (!strcmp(loc_buf,"hamiltonian-colour"))
	{
		if (*c1 == '=')
		{
			mp_control.hamiltonian_colour = ++c1;
		}
		else
		{
			cout << "\nYou must specify a colour if you use the hamiltonian-colour option, e.g. 'hamiltonian-colour=blue'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set hamiltonian-colour = " << mp_control.hamiltonian_colour << endl;

	}
	else if (!strcmp(loc_buf,"HC-include-edge"))
	{
		if (*c1 == '=')
		{
			get_number(mp_control.HC_include_edge,++c1);
		}
		else
		{
			cout << "\nYou must specify an edge to include if you use the HC-include-edge option, e.g. HC-include-edge=12" << endl;
			exit(0);
		}
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: read HC-include-edge read from " << source << ", HC-include-edge = " << mp_control.HC_include_edge << endl;
	}
	else if (!strcmp(loc_buf,"h-units"))
	{
		if (*c1 == '=')
		{
			get_number(mp_control.horizontal_units,++c1);
		}
		else
		{
			cout << "\nYou must specify a horizontal space size for laces in multiples of unit if you use the horizontal-units option, e.g. 'horizontal-units=10'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set horizontal_units = " << mp_control.horizontal_units << endl;

	}
	else if (!strcmp(loc_buf,"hyperbolic"))
	{
		DRAW_IN_HYPERBOLIC_DISC = true;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: drawing diagram in the hyperbolic plane" << endl;
	
	}	
	else if (!strcmp(loc_buf,"knotoid-leg-unbounded"))
	{
    	mp_control.knotoid_leg_unbounded = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: knotoid_leg_unbonded read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"labels"))
	{
    	mp_control.draw_labels = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: draw_labels read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"label-shift"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting label_shift as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.label_shift,++c1);
		}
		else
		{
			cout << "\nYou must specify the number of units u to shift labels if you use the label-shift option, e.g. 'label-shift=60'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: label_shift read from " << source << ", label_shift = " 
	      << mp_control.label_shift << endl;
}	
	}
	else if (!strcmp(loc_buf,"laces"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting LACES true as a result of " << source << " option" << endl;

		LACES = true;
	}
	else if (!strcmp(loc_buf,"left-term-tail-points"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting mp_control.left_terminating_tail_points true as a result of " << source << " option" << endl;

		mp_control.left_terminating_tail_points = true;
	}
	else if (!strcmp(loc_buf,"long"))   // temporary - remove after testing
	{
    	USE_ALL_LONG_GAPS = true;
	}
	else if (!strcmp(loc_buf,"magnify"))
	{
	    if (*c1 == '=')
	    {
			get_number(magnification_factor,++c1);
	    }
		else
		{
			cout << "\nYou must specify a magnification_factor as a percentage if you use the m option, e.g. '-m=10'" << endl;
			exit(0);
		}
		
		mp_control.magnify_small_circles = true;
		magnification_factor = 1 + magnification_factor/100;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: magnifying small circles" << endl;
	debug << "set_programme_long_option: magnification_factor = " << magnification_factor << endl;
}	
	}
	else if (!strcmp(loc_buf,"midpoint"))
	{
		// expect midpoint=1-2:3-4
	    if (*c1 == '=')
	    {
			c1++;
			istringstream ss(c1);
			bool done = false;
			do
			{
				pair<int,int> edge;
				char c = '\0';
				int n;
				ss >> n;
				edge.first = n;
				ss >> c;
				if (c == '-')
				{
					ss >> n;
					edge.second = n;
					mp_control.lace_midpoints.push_back(edge);
					
					ss >> c;
					if (c != ':')
						done = true;
				}
				else
					done=true;
			} while (!done);
	    }
		else
		{
			cout << "\nYou must specify a set of adjacent coordinates between which a mipoint should be added if you use the midpoint option, e.g. '--midpoint=1-2:3-4'" << endl;
			exit(0);
		}

	}
	else if (!strcmp(loc_buf,"midpoints-not-tail-points"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting mp_control.midpoints_not_tail_points true as a result of " << source << " option" << endl;

		mp_control.midpoints_not_tail_points = true;
	}
	else if (!strcmp(loc_buf,"midpoint-tension"))
	{
	    if (*c1 == '=')
	    {
			get_number(mp_control.midpoint_tension,++c1);
	    }
		else
		{
			cout << "\nYou must specify a midpoint tension ( > 1.0) if you use the midpoint-tension option, e.g. 'midpoint-tension=1.2'" << endl;
			exit(0);
		}
    	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: tension read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"no-adjacent-cudgel-midpoints"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: clearing mp_control.adjacent_cudgel_midpoints true as a result of " << source << " option" << endl;

		mp_control.adjacent_cudgel_midpoints = false;
	}
	else if (!strcmp(loc_buf,"no-right-orig-tail-points"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting mp_control.right_originating_tail_points false as a result of " << source << " option" << endl;

		mp_control.right_originating_tail_points = false;
	}
	else if (!strcmp(loc_buf,"no-vertex-axes"))
	{
    	mp_control.show_vertex_axes = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: no-axes read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"no-immersion"))
	{
    	mp_control.draw_immersion = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: no-immersion read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"no-crossings"))
	{
    	mp_control.draw_crossing_features = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: no-crossings read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"no-show-displacement"))
	{
    	mp_control.draw_triangulation_displacement = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: no-show-displacement read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"odd-parity-disc-size"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting odd_parity_disc_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.odd_parity_disc_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a disc size multiplier if you use the odd_parity-disc-size option, e.g. 'odd-parity-disc-size=4'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: odd_parity_disc_size read from " << source << ", odd_parity_disc_size = " 
	      << mp_control.odd_parity_disc_size << endl;
}	
	}
	else if (!strcmp(loc_buf,"oriented"))
	{
    	mp_control.draw_oriented = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: draw_oriented read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"packing"))
	{
    	mp_control.circle_packing = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: circle_packing read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"parity"))
	{
    	mp_control.show_odd_parity_crossings = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: SHOW_ODD_PARITY read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"pen-size"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting pen_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.pen_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a unit size if you use the pen-size option, e.g. 'pen-size=2'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: pen_size read from " << source << ", pen_size = " 
	      << mp_control.pen_size << endl;
}	
	}
	else if (!strcmp(loc_buf,"plot"))
	{
    	PLOT_TRIANGULATION_EDGES = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: PLOT_TRIANGULATION_EDGES read from " << source << endl;

		if (*c1 == '=')
		{
			get_number(plot_steps,++c1);
		}
		else
		{
			cout << "\nYou must specify the required number of divisions if you use the plot option, e.g. 'plot=100'" << endl;
			exit(0);
		}
	}
	else if (!strcmp(loc_buf,"rotate"))
	{

		mp_control.rotate = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: rotate read from " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.rotation_degrees,++c1);
		}
		else
		{
			cout << "\nYou must specify the required rotation in degrees if you use the rotate option, e.g. 'rotate=15'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: rotation_degrees read from " << source << ", rotation_degrees = " 
	      << mp_control.rotation_degrees << endl;
}	
	}
	else if (!strcmp(loc_buf,"scale"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting odd_parity_disc_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.scale,++c1);
		}
		else
		{
			cout << "\nYou must specify a scale size if you use the scale option, e.g. 'scale=1, or scale=2.5'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: scale read from " << source << ", scale = " << mp_control.scale << endl;
	}
	else if (!strcmp(loc_buf,"script-labels"))
	{
    	mp_control.draw_labels = true;
    	mp_control.script_labels = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: script_labels read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"seifert-circles"))
	{
    	mp_control.seifert_circles = true;
    	mp_control.state_smoothed = true;
    	mp_control.smoothed_state_disc_size = 1;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: seifert_circles read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"seifert-edges"))
	{
    	mp_control.seifert_edges = 2;
    	
		if (*c1 == '=')
		{
			/* it can be anything, "odd", "reversed", etc. */
			mp_control.seifert_edges = 1;
		}
    	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: seifert-edges read from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"scriptscript-labels"))
	{
    	mp_control.draw_labels = true;
    	mp_control.scriptscript_labels = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: scriptscript_labels read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"shortcut"))
	{
    	mp_control.draw_shortcut = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: draw_shortcut read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"shrink"))
	{ 
		USE_REGION_SHRINKING_PLACEMENT = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: use region shrinking placement of triangulation vertices following circle packing" << endl;

		if (*c1 == '=')
	    {
			get_number(max_placement_iterations,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set max_placement_iterations = " << max_placement_iterations << endl;

	    }	
	}
	else if (!strcmp(loc_buf,"shrink-factor"))
	{
	    if (*c1 == '=')
	    {
			get_number(region_shrinking_factor,++c1);
	    }
		else
		{
			cout << "\nYou must specify a region_shrinking_factor if you use the shrink-factor option, e.g. 'shrink-factor=0.5'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: region_shrinking_factor = " << region_shrinking_factor << endl;
	}
	else if (!strcmp(loc_buf,"shrink-area-factor"))
	{
	    if (*c1 == '=')
	    {
			get_number(region_shrinking_area_factor,++c1);
	    }
		else
		{
			cout << "\nYou must specify a region_shrinking_area_factor if you use the shrink-area-factor option, e.g. 'shrink-area-factor=1.2'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: region_shrinking_factor = " << region_shrinking_factor << endl;
	}
	else if (!strcmp(loc_buf,"show-shrink"))
	{
		mp_control.draw_shrink_effect = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: draw_shrink_effect read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"show-small"))
	{
		mp_control.highlight_small_edges = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: highlight_small_edges read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"small-shrink"))
	{ 
		USE_REGION_SHRINKING_PLACEMENT = true;
		USE_SMALL_REGION_SHRINKING_PLACEMENT = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: use small region shrinking placement of triangulation vertices following circle packing" << endl;

		if (*c1 == '=')
	    {
			get_number(max_placement_iterations,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set max_placement_iterations = " << max_placement_iterations << endl;

	    }		
	}
	else if (!strcmp(loc_buf,"smallarrowheads"))
	{
    	mp_control.smallarrowheads = true;
    	mp_control.arrowhead_bp_size = 3;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: smallarrowheads read from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"smoothed"))
	{
		mp_control.state_smoothed = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set programme long_option: option mp_control.state_smoothed read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"smoothed-disc-threshold"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting smoothed_state_disc_threshold as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.smoothed_disc_threshold,++c1);
		}
		else
		{
			cout << "\nYou must specify the threshold in units u for the smoothed disc for moving labels if you use the smoothed-disc-threshold option, e.g. 'smoothed-disc-threshold=70'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: smoothed_state_disc_threshold read from " << source << ", smoothed_state_disc_threshold = " 
	      << mp_control.smoothed_disc_threshold << endl;
}	
	}
	else if (!strcmp(loc_buf,"smoothed-disc-size"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting smoothed_state_disc_size as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.smoothed_state_disc_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a disc size multiplier if you use the smoothed-disc-size option, e.g. 'smoothed-disc-size=5'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: smoothed_state_disc_size read from " << source << ", smoothed_state_disc_size = " 
	      << mp_control.smoothed_state_disc_size << endl;
}	
	}
	else if (!strcmp(loc_buf,"state"))
	{
		if (*c1 == '=')
		{
			mp_control.state_smoothed = true;
			mp_control.state = ++c1;
		}
		else
		{
			cout << "\nYou must specify a state if you use the state option, e.g. 'state=AABABA'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: state read from " << source << ", state = " << mp_control.state << endl;
	debug << "set programme long_option: mp_control.state_smoothed set as a result of receiving state option " << source << endl;
}
	}
	else if (!strcmp(loc_buf,"tension"))
	{
    	mp_control.tension = true;

	    if (*c1 == '=')
	    {
			get_number(metapost_path_tension,++c1);
	    }
		else
		{
			cout << "\nYou must specify a matapost path tension ( > 1.0) if you use the tension option, e.g. 'tension=1.2'" << endl;
			exit(0);
		}
    	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: tension read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"triangulation"))
	{
    	mp_control.draw_triangulation = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: triangulation read from " << source << endl;
	}	
	else if (!strcmp(loc_buf,"translate"))
	{
		// expect translate=<vertex>x<x-shift>y<y-shift>:<vertex>x<x-shift>y<y-shift>:...
	    if (*c1 == '=')
	    {
			c1++;
			istringstream ss(c1);
			bool done = false;
			do
			{				
				char c = '\0';
				int v;
				int x;
				int y;
				ss >> v;
				ss >> c;
				if (c == 'x')
				{
					ss >> x;
					ss >> c;
					if (c == 'y')
					{
						ss >> y;
						tuple<int,int,int> translation(v,x,y);
						mp_control.translations.push_back(translation);
					}
					else
					{
						done = true;
					}
				}
				else
				{
					done=true;
				}
				ss >> c;
				if (c != ':')
					done = true;				
			} while (!done);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: added vertex translations" << c1 << " read from " << source << endl;
			
	    }
		else
		{
			cout << "\nYou must specify a set of vertex translations if you use the translate option, e.g. '--translate=24,5,5:34,10,10'" << endl;
			exit(0);
		}
	}
	else if (!strcmp(loc_buf,"uniform-smoothed-discs"))
	{
    	mp_control.uniform_smoothed_discs = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: uniform_smoothed_discs read from " << source << endl;
	}
	else if (!strcmp(loc_buf,"unit"))
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: setting unit as a result of " << source << " option" << endl;

		if (*c1 == '=')
		{
			get_number(mp_control.unit_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a unit size if you use the unit-size option, e.g. 'unit-size=30'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_long_option: unit_size read from " << source << ", unit_size = " 
	      << mp_control.unit_size << endl;
}	      
	}	
	else if (!strcmp(loc_buf,"v-units"))
	{ 
		if (*c1 == '=')
	    {
			get_number(mp_control.vertical_units,++c1);
	    }
		else
		{
			cout << "\nYou must specify a vertical space size for laces in multiples of unit if you use the vertical-units option, e.g. 'vertical-units=10'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: set vertical_units = " << mp_control.vertical_units << endl;

	}
	else if (!strcmp(loc_buf,"vertices"))
	{
    	mp_control.label_vertices = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_long_option: label_vertices read from " << source << endl;
	}
    else
    {
if (debug_control::DEBUG >= debug_control::SUMMARY)
        debug << "set_programme_long_option: invalid long option " << loc_buf << " read from " << source << endl;

        cout << "Invalid long option " << loc_buf << endl;
        exit(0);
    }

	delete[] cptr;
}

void set_programme_short_option(char* cptr, metapost_control& mp_control)
{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: provided with argument: " << cptr << endl;

	char* c1;

	if (strchr(cptr, 'H') && strchr(cptr, '!'))
	{
		if (strchr(cptr,'#'))
		{
			debug_help();
		}
		else
		{

			help_info(true);

//	    	exit(0);
		}
	}

	char* dptr = strchr(cptr, '#');
	if (dptr)
	{

    	/* establish a debug file */
    	debug.open ("draw.dbg"); 

    	if (!debug)
    	{
        	cout << "\nError opening debug file\n";
        	exit(0);
    	}
		else
			debug << boolalpha << "Debug information from draw version " << version << "\n\n";

if (!debug_setup(cptr))  // could probably be dptr, but the original code used cptr
{
	debug_control::DEBUG = debug_control::SUMMARY;
	debug << "set_programme_short_option: default debug options set" << endl;
}

	}

    c1 = strchr(cptr, 'a');
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(average_triangulation_length_factor,++c1);
	    }
		else
		{
			cout << "\nYou must specify an average_triangulation_length_factor if you use the a option, e.g. '-a=0.5'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: average_triangulation_length_factor = " << average_triangulation_length_factor << endl;
	}

	if (strchr(cptr, 'A'))
		IMMERSION_CROSSINGS_AT_RIGHT_ANGLES = true;

	if (strchr(cptr, 'b'))
		USE_BADNESS_OPTIMZATION = false;

	if (strchr(cptr, 'B'))
		INCLUDE_BOUNDARY_VERTICES = true;

	c1=strchr(cptr, 'c');
	if (c1)
	{
		bool error = false;
		
		if (*++c1 == '=')
		{
			c1++;
			
			if (*c1 == 'z')
			{				
				get_number(mp_control.rotation_centre_z,++c1);
				mp_control.implicit_rotation_centre = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: implicit rotation centre z" << mp_control.rotation_centre_z << endl;
	
			}
			else
			{
				if (*c1 == '(')
				{
					get_number(mp_control.rotation_centre_x,++c1);
					while (isdigit(*c1))
						c1++;
					if (*c1 == ',')
					{
						get_number(mp_control.rotation_centre_y,++c1);
						
						while (isdigit(*c1))
							c1++;
						
						if (*c1 == ')')
						{
							mp_control.explicit_rotation_centre = true;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_short_option: explicit rotation centre (" << mp_control.rotation_centre_x 
	      << "," << mp_control.rotation_centre_y << ")" << endl;
}	    
						}
						else
							error = true;

				}
					else
						error=true;
				}
				else
					error=true;				
			}
		}
		else
			error = true;
			
		if (error)
		{
			cout << "\nYou must specify the required rotation centre if you use the c option." << endl;
			cout << "Specify the coordinates of the centre explicitly, e.g. -c=(123.45,67.89), or implicitly, e.g. -c=z21." << endl;
			exit(0);
		}
	}

	c1 = strchr(cptr, 'C');
	if (c1)
	{ 
		if (*++c1 == '=')
	    {
			get_number(mp_control.infinite_cycle,++c1);
	    }
		else
		{
			cout << "\nYou must specify a turning cycle if you use the I option, e.g. '-I=2'" << endl;
			exit(0);
		}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set infinite_cycle = " << mp_control.infinite_cycle << endl;

	}

	c1=strchr(cptr, 'D');
	if (c1)
	{
		if (*++c1 == '=')
		{
			get_number(mp_control.smoothed_state_disc_size,++c1);
		}
		else
		{
			cout << "\nYou must specify a disc size multiplier if you use the D option, e.g. '-D=5'" << endl;
			exit(0);
		}
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: smoothed_state_disc_size = " << mp_control.smoothed_state_disc_size << endl;
	}

	c1 = strchr(cptr, 'd');
	if (c1)
	{ 
		if (*++c1 == '=')
	    {
			get_number(mp_control.disc_size,++c1);
	    }
		else
		{
			cout << "\nYou must specify a unit size if you use the d option, e.g. '-d=20'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set disc_size = " << mp_control.disc_size << endl;

	}

	if (strchr(cptr, 'f'))
	{
		mp_control.draw_lace_frame = true;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set draw_lace_frame true" << endl;
	}

	if (strchr(cptr, 'F'))
	{
		mp_control.draw_crossing_features = false;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set draw_crossing_features false" << endl;
	}

	c1 = strchr(cptr, 'i');
	if (c1)
	{ 
		if (*++c1 == '=')
	    {
			get_number(max_circle_packing_iterations,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set max_circle_packing_iterations = " << max_circle_packing_iterations << endl;

	    }	
	}

	if (strchr(cptr, 'h'))
		help_info(true);

	if (strchr(cptr, 'I'))
		mp_control.draw_immersion = false;


	if (strchr(cptr, 'E'))
	{
		/* use edge repulsion rather than Plastenjak force directed placement */
		PLESTENJAK_FORCE_DIRECTION = false; 

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "set_programme_short_option: using edge repulsion, PLESTENJAK_FORCE_DIRECTION = false, ";
	debug << "APPLY_FORCES_TO_TYPE_12_ONLY = " << (APPLY_FORCES_TO_TYPE_12_ONLY? "true": "false") << endl;
}
	}

	if (strchr(cptr, 'k'))
		mp_control.knotoid_leg_unbounded = true;

	if (strchr(cptr, 'K'))
	{
		USE_KEN_STEPHENSON_CIRCLE_PACKING = false;

if (debug_control::DEBUG >= debug_control::SUMMARY)
//	debug << "set_programme_short_option: using Ken Stephenson's circlepacking algorithm" << endl;
	debug << "set_programme_short_option: using the Fenn circlepacking algorithm, not Ken Stephenson's" << endl;
	}

	if (strchr(cptr, 'l'))
		mp_control.draw_labels = true;

	if (strchr(cptr, 'L'))
	{
		mp_control.draw_labels = true;
		mp_control.label_edges_from_one = true;
	}

    c1 = strchr(cptr, 'M');
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(metapost_coordinate_scale_factor,++c1);
	    }
		else
		{
			cout << "\nYou must specify a metapost_coordinate_scale_factor if you use the M option, e.g. '-M=50'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: metapost_coordinate_scale_factor = " << metapost_coordinate_scale_factor << endl;
	}


	if (strchr(cptr, 'o'))
		mp_control.draw_oriented = true;

	if (strchr(cptr, 'O'))
		mp_control.one_metapost_path = false;

	if (strchr(cptr, 'P'))
		mp_control.circle_packing = true;

    c1 = strchr(cptr, 'p');
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(mp_control.pen_size,++c1);
	    }
		else
		{
			cout << "\nYou must specify a unit size if you use the p option, e.g. '-p=4'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set pen_size = " << mp_control.pen_size << endl;
	}

    c1 = strchr(cptr, 'R');
	if (c1)
	{
		RETRACT_BOUNDARY_VERTICES = true;
		
	    if (*++c1 == '=')
	    {
			get_number(boundary_vertex_retraction_factor,++c1);
	    }
		else
		{
			cout << "\nYou must specify a boundary_vertex_retraction_factor if you use the R option, e.g. '-R=0.6'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: boundary_vertex_retraction_factor = " << boundary_vertex_retraction_factor << endl;
	}

	c1=strchr(cptr, 'r');
	if (c1)
	{
		mp_control.rotate = true;
	    if (*++c1 == '=')
	    {
			get_number(mp_control.rotation_degrees,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: rotation of " << mp_control.rotation_degrees << " degrees" << endl;

	    }
		else
		{
			cout << "\nYou must specify the required rotation in degrees if you use the r option, e.g. '-r=15'" << endl;
			exit(0);
		}
	}

	if (strchr(cptr, 'S'))
		mp_control.draw_shortcut = true;

	if (strchr(cptr, 't'))
    	mp_control.draw_triangulation = true;

	c1 = strchr(cptr, 'T');
	if (c1)
	{ 
    	TRACK_PLACEMENT_ITERATION = true;

		if (*++c1 == '=')
	    {
			get_number(placement_iteration_tracking_step,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set placement_iteration_tracking_step = " << placement_iteration_tracking_step << endl;

	    }	
	}

    c1 = strchr(cptr, 'u'); 
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(mp_control.unit_size,++c1);
	    }
		else
		{
			cout << "\nYou must specify a unit size if you use the u option, e.g. '-u=12'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set unit_size = " << mp_control.unit_size << endl;

	}
	
	if (strchr(cptr, 'v'))
		mp_control.label_vertices = true;

	c1 = strchr(cptr, 'V');
	if (c1)
	{ 
    	CHECK_INNER_HULL_CALCULATION = true;
    	mp_control.draw_triangulation = true;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: check inner hull calculation" << endl;

		if (*++c1 == '=')
	    {
			get_number(check_inner_hull_vertex,++c1);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set check_inner_hull_vertex = " << check_inner_hull_vertex << endl;

	    }	
	    	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set max_placement_iterations = 1" << endl;
	
	}

    c1 = strchr(cptr, 'x');
	if (c1)
	{
	    if (*++c1 == '=')
	    {
			get_number(printf_bool,++c1);
	    }
		else
		{
			cout << "\nYou must specify the circle packing debug level if you use the x option, e.g. '-x=3'" << endl;
			exit(0);
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_programme_short_option: set printf_bool = " << printf_bool << endl;
	}

}

bool get_next_input_string(string input_file,string& input_string, string& title)
{
	bool success = true;

	if (input_file.length())
	{
		/* get next word from input */
		get_input_word(input, input_string, title);
	}
	else
	{
		do
		{
			cout << "\n\ninput: ";
			getline(cin, input_string); // has to be getline so there can be spaces in the string
		} while (input_string.length() == 0);
	}

	if (input_string == "exit" || input_string == "q" || input_string =="Q")
		success = false;
	else if (input_string == "help")
	{
	    help_info(true);
	}
	else    	
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "draw::get_next_input_string: got string " << input_string << endl;
}
	}	
	
	if (success)
	{ 		
		cout << "\n";
		if (title.length() != 0)
			cout << title << endl;
		cout << input_string << endl;
	}
		
	return success;
}

void help_info(bool exit_after_help)
{
	cout << "\n\nThis is A. Bartholomew's programme for drawing knot diagrams, version " << version << "\n";
	
	cout << "\nUsage draw [--<long-option>][-<short-options>][<infile>[<outfile>]]\n";

	cout << "\nThe programme reads labelled peer codes, Gauss codes or planar diagram data from the <infile>, or from\n";
	cout << "a command prompt if no <infile> is provided, and then evaluates a triangulation of the disc determined\n";
	cout << "by the immersion underlying the input code.\n";
	cout << "\nBy default, the programme uses Ken Stephenson's circle packing algorithm to place the vertices of\n";
	cout << "the triangulation.  Alternatively one of the following long programme options may be used to select\n";
	cout << "a different vertex placement technique.\n";
	cout << "\nOnce the vertices have been positioned, the programme creates a metapost script describing the\n";
	cout << "diagram in <outfile>, or in the file draw.out if no <outfile> is provided.\n";

	cout << "\nProgramme control <long-option>\n";
	cout << "  convex-disc: draw the convex, straight line triangulation of a disc\n";
	cout << "  edge[=<max-iterations>] default 200: use edge distribution placement\n";
	cout << "  force[=<max-iterations>] default 200: use Plestenjak force directed placement\n";
	cout << "  gravity[=<max-iterations>] default 200: use centre of gravity placement\n";	
	cout << "  hamiltonians: draw all Hamiltonian circuits for a given diagram\n";	
	cout << "  hamiltonian_circuit[=e_1 ... e_n]: draw a single Hamiltonian edge circuit for a given diagram,\n";
	cout << "                                     if no explicit edge circuit is given, draw any Hamiltonian circuit\n";	
	cout << "  laces: draw lace diagrams from the programme input\n";
	cout << "  magnify[=<percentage>] default 0: magnify small circles in circle packing by specified percentage\n";
	cout << "  seifert-circles: draw the Seifert circles determined by the label orientation of the given diagram\n";
	cout << "  shrink[=<max-iterations>] default 200: use region shrinking placement after circle packing\n";
	cout << "  small-shrink[=<max-iterations>]: use small region shrinking placement after circle packing\n";
	cout << "  smoothed: draw all of the smoothed states for the given diagram\n";

	cout << "\nMetapost control <short-option>\n";
	cout << "  c=<centre>: specify the centre of rotation <centre> = (<x>,<y>)|z<n>\n";
	cout << "  C=<cycle>: specify the turning cycle to bound the infinite region\n";
	cout << "  D: set the smoothed state disc size multiplier to determine the size of smoothed crossings n*d (default n=6) \n";
	cout << "  d: set the unit multiplier to determine the crossing disc diameter (default 30)\n";
	cout << "  f: draw the entire frame when drawing laces\n";
	cout << "  F: do not draw the crossing features\n";
	cout << "  h: help screen\n";
	cout << "  I: do not draw the immersion (consider using F option also)\n";
	cout << "  k: draw knotoids with the leg in the unbounded component of the immersion's complement\n";
	cout << "  l: add edge labels to the diagram\n";
	cout << "  o: draw orientation for knots (orientation always shown for knotoids)\n";
	cout << "  p: set the pen multiplier n to determine the pencircle scale n*0.5pt (default n=1)\n";
	cout << "  P: draw the underlying circle packing that determines the diagram\n";
	cout << "  r=<degrees>: rotate the diagram anti-clockwise by the specified number of degrees\n";
	cout << "  S: draw the shortcut if the code describes a knotoid\n";
	cout << "  t: draw the triangulation of the unit disc\n";
	cout << "  u: set the unit dimension nn to determine u=0.nn points (default 20)\n";
	cout << "  v: label the vertices and draw coordinate axes\n";
	cout << endl;

	cout << "The Metapost control short options have corresponding long options, useful for including in\n";
	cout << "input files as global options or as qualifiers:\n";		
	cout << "  centre=<centre>: specify the centre of rotation <centre> = (<x>,<y>)|z<n>\n";
	cout << "  cycle=<cycle>: specify the turning cycle to bound the infinite region\n";
	cout << "  dots: draw knotoid shortcuts dashed with dots rather than dashed evenly\n";
	cout << "  disc-size: set the unit multiplier to determine the crossing disc diameter (default 30)\n";
	cout << "  draw-lace-frame: draw the entire frame when drawing laces\n";
	cout << "  knotoid-leg-unbounded: draw knotoids with the leg in the unbounded component of the immersion's complement\n";
	cout << "  labels: add immersion edge labels to the diagram: by default immersion edge labels are numbered from 0 (see the L option)\n";
	cout << "  no-crossings: do not draw the crossing features\n";
	cout << "  no-immersion: do not draw the immersion\n";
	cout << "  oriented: draw orientation for knots (orientation always shown for knotoids)\n";
	cout << "  packing: draw the underlying circle packing that determines the diagram\n";
	cout << "  pen-size: set the pen multiplier n to determine the pencircle scale n*0.5pt (default n=1)\n";
	cout << "  rotate=<degrees>: rotate the diagram anti-clockwise by the specified number of degrees\n";
	cout << "  shortcut: draw the shortcut if the code describes a knotoid\n";
	cout << "  triangulation: draw the triangulation of the unit disc\n";
	cout << "  unit: set the unit dimension nn to determine u=0.nn points (default 20)\n";
	cout << "  vertices: label the vertices and draw coordinate axes\n";
		
	cout << "\nadditional long options\n";
	cout << "  colour: draw different components with colours rather than just black lines\n";
	cout << "  colour-map=<filename>: use the colours in <filename> rather than the default colours\n";
	cout << "  cudgel-space=<float> default 1.0: amount by which the calculated inter-cudgels gap for laces is expanded\n";
	cout << "  cusp_disc_size=<n>: set the cusp disc size multiplier to determine the size of the non-Seifert-smoothed cusp indicator discs, 0.n*disc-size (default n=7)\n";
	cout << "  edge-factor=<float> default 0.5: amount by which vertices are moved towards the COG in edge distribution placement\n";
	cout << "  first-gap: always use the first gap as the active gap in edge distribution placement\n";
	cout << "  frame-corners: show frame corners when tracking placement iteration\n";
	cout << "  gauss-crossings: show labels for the Gauss crossings, as specified by Gauss code or planar diagram input, or calculated from peer code input\n";
	cout << "  gauss-labels: show labels for Gauss arcs, not immersion arcs: gauss-labels are numbered from 1\n";
	cout << "  grid=<grid-size>: draw a grid to assist with using the translate option, default 10 (percent of diagram width or height)\n";
	cout << "  hamiltonian_colour=<string>: the colour used for Hamiltonian circuits, default green\n";
	cout << "  h-units: size of horizontal spacing for laces, in multiples of unit\n";
	cout << "  hyperbolic: drawing diagram in the hyperbolic plane\n";
	cout << "  label-shift=<n>: set the number of units u by which the labels of small smoothed crossings or Gauss crossings should be moved (default n=50)\n";
	cout << "  left-term-tail-points: when drawing laces include tail points for body disc arcs terminating on the left side of a cudgel\n";
	cout << "  midpoints=<arc-list>: specify those adjacent coordinates in a lace diagram between which a mipoint should be added\n";
	cout << "                        <arc-list> =  p-q:r-s:t-u..., where midpoints are required between zp and zq, zr and zs, zt and zu etc.\n";
	cout << "  midpoints-not-tail-points: add midpoints rather than tail points when drawing laces\n";
	cout << "  midpoint-tension: metapost path tension for lace midpoints, default value 1.0, i.e. the metapost \"..\" default\n";
	cout << "  no-adjacent-cudgel-midpoints: when drawing laces do not include midpoints in body disc arcs joining adjacent cudgels\n";			
	cout << "  no-vertex-axes: do not show the axes when labelling the vertices of a diagram\n";			
	cout << "  no-right-orig-tail-points: when drawing laces do not include tail points for body disc originating on the right side of a cudgel\n";
	cout << "  no-show-displacement: do not show the triangulation displacement when tracking placement iteration\n";
	cout << "  odd-parity-disc-size=<n>: set the odd parity disc size multiplier to determine the size of the odd parity crossing indicator, n*0.2*disc-size (default n=6)\n";
	cout << "  plot=<divisions>: set the number of divisions for the histogram in edge distribution placement (default 20)\n";
	cout << "  parity: show crossings with odd parity in smoothed states\n";
	cout << "  scale=<float>: override other size settings and scale diagram to the specified multiple of a standard size\n";
	cout << "  script-labels: label vertices using TeX's script size font\n";
	cout << "  scriptscript-labels: label vertices using TeX's scriptscript size font\n";
	cout << "  seifert-edges: highlight the odd or even immersion edges, so edges belonging to the same Seifert circle have the same colour\n";
	cout << "  show-shrink: draw the effect of shrinking the triangulation when using region shrinking placement\n";
	cout << "  show-small: highlight the small edges of the triangulation when using edge_distribution placement\n";
	cout << "  shrink-factor=<float> default 0.75: amount by which region shrinking placement retracts trianglulation vertices towards the barycentre\n";
	cout << "  shrink-area-factor=<float> default 1.618034: region shrinking placement converges if all compact regions that have an area within this\n"; 
	cout << "                                               factor of the average area\n";
	cout << "  smoothed-disc-size=<n>: set the smoothed state disc size multiplier to determine the size of smoothed crossings n*disc-size (default n=6)\n";	
	cout << "  smoothed-disc-threshold=<n>: set the diameter in units u of the smoothed crossing disc below which the label of the crossings should be moved (default n=30)\n";
	cout << "  state: specify the smoothed state that you wish to draw as a string of A and B characters corresponding to the Gauss crossings of the diagram\n";
	cout << "  transate=<translation-list>: specify a list of vertices together with a translation in the form UxPyQ:VxRyS... where U and V are vertex numbers\n";
	cout << "                               and P,Q,R,S percentages of the diagram width and height e.g. 10x8y-33 indicates shifting vertex 10 +8% to the right -3% up\n";
	cout << "  tension: set the metapost path tension, default value 1.0, i.e. the metapost \"..\" default\n";
	cout << "  uniform-smoothed-discs: always draw state smoothed discs of the size specified by smoothed-disc-size\n";
	cout << "  v-units: size of vertical spacing for laces, in multiples of unit\n";
	
	cout << "\nadditional short options\n";
	cout << "  #: debug\n";
	cout << "  a=<float> default 1.0: average triangulation edge length factor\n";
	cout << "     used to magnify edge vertices closer than a*average_triangulation_edge_length\n";
	cout << "  A: adjust the angle of diagram arcs at crossings to be right angles\n";
	cout << "  b: do not use badness optimization\n";
	cout << "  B: include boundary vertices; may be used with f, g, P, or t.\n";
	cout << "     Boundary vertices are alway included in force directed placement,\n";
	cout << "     in which case this option only has an effect if the t option is also used.\n";
	cout << "  E: use edge repulsion rather than Plestenjak force directed placement\n";
	cout << "  i[=<max-iterations>] default 1000: set maximum circle packing iterations\n";
//	cout << "  H!: additional help screen\n";
	cout << "  #H!: display debug help screen\n";
	cout << "  K: use Fenn's circle packing algorithm rather than Ken Stephenson's\n";
	cout << "  L: start from one when labelling edges\n";
	cout << "  M=<scale_factor> : set metapost_coordinate_scale_factor: default 2500 with the circle packing\n";
	cout << "     option, 1000 with force or gravity placement, 600 with the convex-disc option, 25 otherwise\n";			
//	cout << "  O: create metapost with one cycled path for each component\n";
	cout << "  O: turn off default one_metapost_path and use two paths for each component with explicit direction, rather than a cycle\n";
	cout << "  P: draw circle packing\n";
	cout << "  R=<float>: retract boundary vertices radially towards their centroid by a factor <float> < 1\n";
	cout << "  T[=<track-step>] default 1: track placement iteration (force, gravity, shrink)\n";
	cout << "  V[=<vertex>] default 0: check inner hull calculation of <vertex>\n";
	cout << "  x=[0|1|2|3] default 0: set debug level for Stephenson circle packing\n";			

	cout << endl;
	
	if (exit_after_help)
		exit(0);
}

/***************  Functions required for setting up draw specific debug **************/

void set_main_debug_option_parameter(char* pptr, string option);
void set_main_debug_option(char* start, char* end)
{
	char  loc_buf[end-start+2];
	char* c1 = start;
	char* c2 = loc_buf;

	/* if both start and end are zero, display debug help information.
	   this has been included here in this manner so that each time a debug option
	   is added, the help will be updated (hopefully!)
	*/
	if (start == 0 && end == 0)
	{
		set_main_debug_option_parameter(0,"draw");
//		cout << "\t\tvogel, boolean" << endl;
		return;
	}

	do
	{
		*c2++ = *c1++;
	} while (c1 <= end);
	
	*c2 = '\0';

	char* pptr = strchr(loc_buf,'{');
	if (pptr)
	{
		*pptr++ = '\0';
	}

	/* now, even if there are parameters, the first part of loc_buf 
	   is a C-string that identifies the option */
	
	if (!strcmp(loc_buf,"draw"))
	{
		debug_control::DEBUG = debug_control::SUMMARY;
		debug << "main::set_main_debug_option: setting debug option debug_control::DEBUG = debug_control::SUMMARY\n";		
		
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, "draw");
		}
	}

	debug.flush();

}

void set_main_debug_option_parameter(char* pptr, string option)
{
	if (option == "draw")
	{
		if (!pptr)
		{
			cout << "\t\tdraw{summary|1:basic|2:intermediate|3:detail|4}, integer: default 0=off, no parameters sets summary" << endl;
		}
		else
		{
			if (!strcmp(pptr,"summary") || !strcmp(pptr,"1") )
			{
				debug_control::DEBUG = debug_control::SUMMARY;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::SUMMARY\n";		
			}
			if (!strcmp(pptr,"basic") || !strcmp(pptr,"2") )
			{
				debug_control::DEBUG = debug_control::BASIC;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::BASIC\n";		
			}
			else if (!strcmp(pptr,"intermediate") || !strcmp(pptr,"3"))
			{
				debug_control::DEBUG = debug_control::INTERMEDIATE;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::INTERMEDIATE\n";		
			}
			else if (!strcmp(pptr,"detail") || !strcmp(pptr,"4"))
			{
				debug_control::DEBUG = debug_control::DETAIL;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::DETAIL\n";		
			}
			else if (!strcmp(pptr,"exhaustive") || !strcmp(pptr,"5"))
			{
				debug_control::DEBUG = debug_control::EXHAUSTIVE;
				debug << "main::set_main_debug_option_parameter: setting debug option debug_control::DEBUG = debug_control::EXHAUSTIVE\n";		
			}
		}
	}
}

void main_display_default_options()
{
	cout << "\t\tdraw{summary}" << endl;
}


bool connected_code_data(generic_code_data& code_data)
{
	int num_components = code_data.num_components;

	if (num_components == 1)
	{
		return true;
	}
	else
	{
		int num_crossings = code_data.num_crossings;
		matrix<int>& code_table = code_data.code_table;
	
		bool connected = true;
	
		for (int i=0; i< num_components; i++)
		{
			/* If the odd peers of all the naming edges belonging to this component
			   lie within the range of edges on this component, then this component
			   is disconnected from the rest of the diagram.
			*/
			connected = false; // assume disconnected and seek evidence to the contrary.
			
			int first_edge = code_data.first_edge_on_component[i];
			int last_edge = (i==num_components-1? 2*num_crossings-1: code_data.first_edge_on_component[i+1]-1); 		
			
			for (int j=0; j < code_data.num_component_edges[i]/2; j++)
			{
				int odd_peer = code_table[generic_code_data::table::OPEER][first_edge/2+j];
				if (odd_peer < first_edge || odd_peer > last_edge)
				{
					connected = true;
					break;
				}
			}
			
			if (!connected)
				break;
		}
		
		return connected;
	}
}

void set_frame_corners(string coordinate_file, metapost_control& mp_control)
{
    
	ifstream Oinput;
	char* coords = c_string(coordinate_file);
	Oinput.open(coords);
	if (!Oinput)
	{
		cout << "\nError opening output file " << coordinate_file << endl;
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "set_frame_corners: could not open " << coordinate_file << endl;
		exit(0);
    }

	int num_vertices;
	Oinput >> num_vertices;

    /* read the vertex coordinates */
    matrix<double> vcoords(num_vertices,2);

	for (int i=0; i< num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			Oinput >> vcoords[i][j];
	}


if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "set_frame_corners: vertex_coordinates: " << endl;
    for (int i=0; i< num_vertices; i++)
    {
		debug << "set_frame_corners:   vertex " << i << ": ";
		for (int j=0; j<2; j++)
			debug << vcoords[i][j] << ' ';
		debug << endl;
	}	
}

	for (int i=0; i< num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
			vcoords[i][j] *= metapost_coordinate_scale_factor;
	}

	/* determine the frame extremities */
	double frame_minx = vcoords[0][0];
	double frame_maxx = vcoords[0][0];
	double frame_miny = vcoords[0][1];
	double frame_maxy = vcoords[0][1];

	for (int i=1; i< num_vertices; i++)
	{
		for (int j=0; j< 2; j++)
		{
			if (j==0 && vcoords[i][j] < frame_minx)
				frame_minx = vcoords[i][j];
			
			if (j==1 && vcoords[i][j] < frame_miny)
				frame_miny = vcoords[i][j];
			
			if (j==0 && vcoords[i][j] > frame_maxx)
				frame_maxx = vcoords[i][j];
			
			if (j==1 && vcoords[i][j] > frame_maxy)
				frame_maxy = vcoords[i][j];
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "set_frame_corners: frame_minx = " << frame_minx << ", frame_maxx = " << frame_maxx << ", frame_miny = " << frame_miny << ", frame_maxy = " << frame_maxy << endl;

	mp_control.frame_minx = frame_minx;
	mp_control.frame_maxx = frame_maxx;
	mp_control.frame_miny = frame_miny;
	mp_control.frame_maxy = frame_maxy;
	
	delete[] coords;
}

