/**************************************************************************************************

void draw_lace (metapost_control& mp_control, string input_string)
void read_lace_body(lace_body& body, string s)
void read_lace_skeleton(lace_skeleton& skeleton,string s)
ostream& operator << (ostream& os, const lace_body& body)
void print_lace_body (ostream& os, const lace_body& body, string prefix)
ostream& operator << (ostream& os, const lace_skeleton& skeleton)
bool valid_lace(lace_body& input_lace)
void assign_pairs(matrix<int>& edge_pairs, int& pair_marker, int start_label_1, int start_label_2, int start_label_3, int count_1, int count_2, int count_3, bool boundary_place_3)
void write_lace_metapost(ofstream& os, lace_body& input_lace, metapost_control& mp_control, matrix<int>& arc_end_labels,matrix<int>& path_body_arc)

**************************************************************************************************/

using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <complex>

/********************* External variables ***********************/
extern ofstream     debug;
extern ofstream     output;

extern bool TEST_MODE; 
extern bool VALID_BOUNDARY_COUNTS;  // used for testing
extern bool SHOW_BODY_ARC_END_LABELS;
extern bool ALLOW_HEAD_HOOKS;
extern float metapost_path_tension;

#include <util.h>
#include <matrix.h>
#include <draw.h>

/********************* Function prototypes ***********************/
bool valid_lace(lace_body& input_lace);
void assign_pairs(matrix<int>& edge_pairs, int& pair_marker, int start_label_1, int start_label_2, int start_label_3, int count_1, int count_2, int count_3, bool boundary_place_3);
void write_lace_metapost(ofstream& os, lace_body& input_lace, metapost_control& mp_control, matrix<int>& arc_end_labels,matrix<int>& path_body_arc);
void report_markers(matrix<int>& marker, matrix<int>& coordinate);
void set_out_head_coordinates(int side, int cudgel_x, float head_point_y_coord, int& z_index, matrix<int>& head_point_type_count, matrix<int>& head_point_coordinate, 
                              matrix<int>& tail_point_type_count, matrix<int>& tail_point_coordinate, ofstream& os);
void set_out_tail_coordinates(int side, int cudgel_x, int max_num_lower_bones, int& z_index, matrix<int>& head_point_type_count, matrix<int>& head_point_coordinate, 
                              matrix<int>& tail_point_type_count, matrix<int>& tail_point_coordinate, ofstream& os);
void add_lace_midpoints(ofstream& os, metapost_control& mp_control, int a, int b);


void draw_lace (metapost_control& mp_control, string input_string)
{
	lace_body input_lace;
	read_lace_body(input_lace, input_string);
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	print_lace_body(debug,input_lace,"draw_lace: ");
	
	if (valid_lace(input_lace))
	{						
		/* Then number of boundary points in the skeleton disc is 2(r+k) where r is the sum of the 
		   bone counts above and below the cudgel heads on the skeleton and k is the number of cudgels
		*/
		lace_skeleton& skeleton = input_lace.skeleton;
		int num_cudgels = skeleton.num_cudgels;	
		int num_boundary_points = input_lace.num_boundary_points;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: number of skeleton disc boundary points = " << num_boundary_points << endl;

		/* Identify the location of the heads and tails in the skeleton boundary.  For each cudgel, work
		   down the left side to the tail then back up the right side, thereby orienting the boundary disc
		   anti-clockwise from the basepoint.
		*/
		vector<int> boundary_point(num_boundary_points);
		int index=0;
		for (int i=0; i< num_cudgels; i++)
		{
			for (int j=0; j < skeleton.bones[i][1]; j++)		
				boundary_point[index++] = lace_control::boundary_type::BONE;
			
			if (skeleton.head_polarity[i] == lace_control::polarity::LEFT)
				boundary_point[index++] = lace_control::boundary_type::HEAD;

			for (int j=0; j < skeleton.bones[i][0]; j++)		
				boundary_point[index++] = lace_control::boundary_type::BONE;
				
			boundary_point[index++] = lace_control::boundary_type::TAIL;
				
			for (int j=0; j < skeleton.bones[i][0]; j++)		
				boundary_point[index++] = lace_control::boundary_type::BONE;
				
			if (skeleton.head_polarity[i] == lace_control::polarity::RIGHT)
				boundary_point[index++] = lace_control::boundary_type::HEAD;

			for (int j=0; j < skeleton.bones[i][1]; j++)		
				boundary_point[index++] = lace_control::boundary_type::BONE;
		}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw_lace: boundary points: ";
	for (int i=0; i< num_boundary_points; i++)
	{
		switch (boundary_point[i])
		{
			case lace_control::boundary_type::HEAD: debug << "H ";break;
			case lace_control::boundary_type::TAIL: debug << "T ";break;
			default: debug << "B ";
		}
	}
	debug << endl;
}

		/* determine the skeleton permutation */
		vector<int> skeleton_perm(num_boundary_points);
		int first_label = 1;
		for (int i=0; i< skeleton.num_cudgels; i++)
		{
			int num_bones =  skeleton.bones[i][0] + skeleton.bones[i][1];
			int last_label = first_label + 2 * num_bones + 1;
			int initial_last_label = last_label;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: cudgel " << i << " first label = " << first_label << ", num_bones = " << num_bones << ", last_label = " << last_label << endl;
							
			for (int j=0; j< num_bones; j++)
			{
				if (boundary_point[first_label-1] == lace_control::boundary_type::HEAD)
					first_label++;

				if (boundary_point[last_label-1] == lace_control::boundary_type::HEAD)
					last_label--;									
					
				skeleton_perm[first_label-1] = last_label;
				skeleton_perm[last_label-1] = first_label;
				first_label++;
				last_label--;								
			}
			
			first_label = initial_last_label+1;
		}
						
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw_lace: skeleton_perm: ";
	for (int i=0; i< num_boundary_points; i++)
		debug << skeleton_perm[i] << ' ';
	debug << endl;
}						

		int num_loops = input_lace.num_loops;
		
		/* There are num_loops+1 boundary 1-cells on the body disc.  The number of interior 1-cells and the number of 2-cells depends 
		   on whether there is more than one loop.  If we have just one loop, we add two boundary vertices either side of and immediately
		   adjacent to the base point, resulting in one interior 1-cell and two 2-cells.  If there are at least two loops the body disc 
		   divides naturally into 2-cells with vertices the basepoint and one in each loop bigon boundary.  Thus there are num_loops-1 
		   2-cells and num_loops-2 interior 1-cells.
		*/

		/* Identify how many points we have on each boundary 1-cell */
		vector<int> boundary_1_cell_count(input_lace.num_boundary_1_cells);
		
		int offset = 0; 
		index = 0;
		for (int i=0; i< num_loops; i++)
		{
			boundary_1_cell_count[index++] = input_lace.loop_location[i] - offset;
			offset = input_lace.loop_location[i];
		}
		boundary_1_cell_count[index] = num_boundary_points - offset;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw_lace: boundary 1-cell counts ";
	for (int i=0; i< input_lace.num_boundary_1_cells; i++)
		debug <<  boundary_1_cell_count[i] << ' ';
	debug << endl;
}
				
		
		/* We work towards tracing the arcs across the body disc to identify arc end labels.  We regard the 
		   body arcs as a graph with vertices the intersection with the 1-skeleton of the triangulated body 
		   disc and edges that cross a corner of a triangulated 2-cell.
							   
		   If the edges of a triangular 2-cell contain A, B and C points respectively, so that A+B+C=2m 
		   and if the opposite corner to A, B C has a,b,c edges across the corner respectively, then a+b+c=m.
		   Then, since A=b+c=m-a, we have a=m-A and similarly b=m-B and c=m-C.

		   We adopt a canonical order that is the mirror image of that used with tracks: thus, we label the 
		   boundary points anti-clockwise around the body disc and then along each interior 1-cell, again 
		   working from the left.  
		   
		   We then record the edges of our graph as pairs of these labels.  We work anti-clockwise through the 
		   triangular 2-cells and work anticlockwise around the corners of each one, starting from the base 
		   point.  For each corner, we work inwards towards the corner, assingning pairs of edge labels to 
		   the rows of edge_pairs so that column zero records the label on the 1-cell that appears earliest 
		   in the canonical order.
		   
		   if num_loops == 1 there are virtual additional vertices and num_2_cells = 2 with one interior 1-cell, 
		   if num_loops == 2 there is only one 2-cell and no interior 1-cells
		   
		   As noted above, each 2-cell has m edges, so the total number of edges across the disc is half the number
		   of boundary points plus the sume of the interior 1-cell counts.
		*/
		int num_pairs = num_boundary_points/2;
		for (int i=0; i< input_lace.num_interior_1_cells; i++)
			num_pairs += input_lace.interior_1_cell_count[i];							

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: total number of edge pairs =  " << num_pairs << endl;
							
		matrix<int> edge_pairs(num_pairs,2);

		int start_label_1;
		int start_label_2;
		int start_label_3;
		int count_1;
		int count_2;
		int count_3;
		int pair_marker = 0;

		if (num_loops == 1)
		{
			start_label_1 = 0;
			start_label_2 = 1;
			count_1 = 0;
			count_2 = boundary_1_cell_count[0];
			start_label_3 = 2*count_2+1;
			count_3 = count_2;
							
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: num_loops == 1, first triangular 2-cell" << endl;
							
			assign_pairs(edge_pairs, pair_marker,start_label_1,start_label_2,start_label_3,count_1,count_2,count_3,false);

			start_label_1 = start_label_3;
			start_label_2 += count_2;
			count_1 = count_3;
//			count_2 = count_2;
			start_label_3 = 0;
			count_3 = 0;
							
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: num_loops == 1, second triangular 2-cell" << endl;
							
			assign_pairs(edge_pairs, pair_marker,start_label_1,start_label_2,start_label_3,count_1,count_2,count_3,true);
		}
		else if (num_loops == 2)
		{
			start_label_1 = 1;
			start_label_2 = boundary_1_cell_count[0]+1;
			start_label_3 = start_label_2 + boundary_1_cell_count[1];
			count_1 = boundary_1_cell_count[0];
			count_2 = boundary_1_cell_count[1];
			count_3 = boundary_1_cell_count[2];
							
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: num_loops == 2, single triangular 2-cell" << endl;
							
			assign_pairs(edge_pairs, pair_marker,start_label_1,start_label_2,start_label_3,count_1,count_2,count_3,true);
		}
		else
		{  							
			bool boundary_place_3;
			for ( int i = 0; i < input_lace.num_2_cells; i++)
			{
				if (i==0) /* Type A 2-cell. */
				{
					start_label_1 = 1;
					start_label_2 = boundary_1_cell_count[0]+1;
					start_label_3 = num_boundary_points + 1;
					count_1 = boundary_1_cell_count[0];
					count_2 = boundary_1_cell_count[1];
					count_3 = input_lace.interior_1_cell_count[0];
					boundary_place_3 = false;
				}
				else if (i==input_lace.num_2_cells-1) /* Type C 2-cell. */
				{
					start_label_1 = start_label_3;
					start_label_2 += count_2;
					count_1 = count_3;
					count_2 = boundary_1_cell_count[i+1];
					start_label_3 = start_label_2 + count_2;
					count_3 = boundary_1_cell_count[i+2];
					boundary_place_3 = true;
				}
				else /* Type B 2-cell. */
				{
					start_label_1 = start_label_3;
					start_label_2 += count_2;
					start_label_3 += count_3;
					count_1 = count_3;
					count_2 = boundary_1_cell_count[i+1];
					count_3 = input_lace.interior_1_cell_count[i];
					boundary_place_3 = false;
				}				    	

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: 2-cell " << i << endl;
	
				assign_pairs(edge_pairs, pair_marker,start_label_1,start_label_2,start_label_3,count_1,count_2,count_3,boundary_place_3);
			}
		}	

		/* Create a disc permutation reflecting the edge_pairs, so that disc_perm[i]=j
		   If {i,j} is a row of edge_pairs.

		   1.  If a label in the first column of edge_pairs is <= num_boundary_points, it is 
		   the initial label of a component arc on the disc.  If it is > num_boundary_points, 
		   the pair in that row determines an edge whose initial point lies on an interior 1-cell.
		   In this case the label on this interior point will have occured exactly once, in the
		   second column of a row of "pair[i]", BEFORE the current pair.
		   
		   2.  If a label in the second column is <= num_boundary_points, it is the terminal 
		   point of a component arc on the disc.  If it is > num_boundary_points, the pair in 
		   that row determines an edge of t whose terminal point lies on an interior 1-cell.  
		   In this case the label at this terminal point will occur exactly once, in the first
		   column of a row of "pair[i]", AFTER the current pair.
		   
		   The disc_perm allows us to find the component arcs.  If a label in the second column of 
		   edge_pairs is <= num_boundary_points, it terminates a component arc and so shall be made	
		   to index zero in disc_perm to signify this.  If the label in the second column is
		   > num_boundary_points it will also occur in the first column in a subsequent row so no 
		   assignment is made.
		*/
		int total_num_points = num_pairs+num_boundary_points/2;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: total number of points =  " << total_num_points << endl;
						
		vector<int> disc_perm(total_num_points);
		for (int i=0; i< num_pairs; i++)
		{
			/* edge_pairs records labels numbered from 1 */
			disc_perm[edge_pairs[i][0]-1] = edge_pairs[i][1];						
			
			if (edge_pairs[i][1] <= num_boundary_points)
				disc_perm[edge_pairs[i][1]-1] = 0;									
		}
						
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw_lace: disc_perm: ";
	for (int i=0; i< total_num_points; i++)
		debug << setw(3) << i+1;
	debug << "\ndraw_lace:            ";
	for (int i=0; i< total_num_points; i++)
		debug << setw(3) << disc_perm[i];	
	debug << endl;
}
		/* count the number of disc arc components */
		int num_disc_arc_components = 0;
		for (int i=0; i< num_boundary_points; i++)
		{
			if (disc_perm[i] == 0)
				num_disc_arc_components++;
		}
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: number of disc arc components = " << num_disc_arc_components << endl;
	
		/* Identify the arc end labels, we trace through disc_perm, setting values negative 
		   to record that we've visited them.  This is why we record the labels from 1 not zero.
		*/
		matrix<int> arc_end_labels(num_disc_arc_components,2);
		for (int i=0; i< num_disc_arc_components; i++)
		{
			/* look for the start of the next component */
			int current_label;
			for (int j=0; j< num_boundary_points; j++)
			{
				if (disc_perm[j] > 0)
				{
					/* we want labels numbered from 1, not offsets numbered from zero */
					current_label = j+1;
					break;
				}
			}
			
			arc_end_labels[i][0] = current_label;
							
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: disc arc components starting at = " << current_label << ": ";
							
			int next_label;
			do
			{
				next_label = disc_perm[current_label-1];							
								
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << next_label << ' ';
	
				if (next_label == 0)
					arc_end_labels[i][1] = current_label;
					
				disc_perm[current_label-1] *= -1 ;
				current_label = next_label;
				
			} while (next_label != 0); 
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << endl;
	
		}
						

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "\ndraw_lace: arc_end_labels: " << endl;
	print (arc_end_labels, debug, 3,"draw_lace: ");
	debug << endl;
}

		/* Now we know the arc end labels in the body we can trace the lace paths by combining the 
		   bones of the skeleton with the arcs in the body disc. 
		   
		   We set up two permutations, skeleton_perm and body_perm that reflect the arcs across the 
		   two discs.  In the skeleton disc, tails and heads are not connected to anything, so we 
		   make these indices reference 0.
		*/
		vector<int> body_perm(num_boundary_points);
		for (int i=0; i< num_disc_arc_components; i++)
		{
			body_perm[arc_end_labels[i][0]-1] = arc_end_labels[i][1];
			body_perm[arc_end_labels[i][1]-1] = arc_end_labels[i][0];
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "draw_lace: body_perm: ";
	for (int i=0; i< num_boundary_points; i++)
		debug << body_perm[i] << ' ';
	debug << endl;
}						
						
		/* The body_perm and skeleton_perm allow us to determine the connected components obtained by identifying the boundaries
		   of the two discs. Those components may therefore be regarded as graphs on a 2-sphere whose vertices all lie on an equator and 
		   whose arcs join distinct vertices, lying alternatively in the northern and southern hemisphere.  The only vertices in this 
		   graph with valency one are the head and tail vertices, all other vertices have valency two.
		
		   Thus, starting at each tail, we traverse a body arc and check whether it reaches a head, a tail or a bone.  We set the body_perm
		   elements corresponding to the arc endpoints to be negative to indicate that the arc has been traversed.  If a body arc reaches a
		   head, we have traced this path and move onto the next tail.  If we reach another tail, we do not have a valid lace. If we reach a 
		   bone, we look at the skeleton perm to identify the start of the next body arc.
		   
		   If each tail lies in the same component as a head, we check that there are no additional components that are simple closed curves
		   on the sphere.  That is, are there body arcs we have not yet traversed.  If there are, then again we do not have a valid lace.
		   
		   We record the sequence of labels at the start of each body arc in path_body_arc.
		*/
		bool valid_lace = true;
		int start_label = 1;
		
		matrix<int> path_body_arc(num_cudgels, num_disc_arc_components+1);

		for (int i=0; i< num_cudgels && valid_lace; i++)
		{
			int num_bones = skeleton.bones[i][0] + skeleton.bones[i][1];
			int tail_label = start_label + num_bones + (skeleton.head_polarity[i] == lace_control::polarity::LEFT?1:0);

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace: tail " << i << endl;
							
			int body_arc_initial_label = tail_label;						
			int index = 1;
			do
			{
				path_body_arc[i][index++] = body_arc_initial_label;
				int body_arc_terminal_label = body_perm[body_arc_initial_label-1];
				body_perm[body_arc_initial_label-1] *= -1;
				body_perm[body_arc_terminal_label-1] *= -1;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_lace:   body_arc_initial_label = " << body_arc_initial_label << ", body_arc_terminal_label = " << body_arc_terminal_label;
								
				if (boundary_point[body_arc_terminal_label-1] == lace_control::boundary_type::HEAD)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " arrives at a head" << endl;
	
					break;
				}
				else if (boundary_point[body_arc_terminal_label-1] == lace_control::boundary_type::TAIL)
				{
						
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " arrives at a tail, invalid lace" << endl;
	
					valid_lace = false;
				}
				else
				{
					body_arc_initial_label = skeleton_perm[body_arc_terminal_label-1];
						
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << " arrives at a bone, connects to " << body_arc_initial_label << endl;
									
				}
			} while (valid_lace);
			
			path_body_arc[i][0] = index-1;

			start_label += 2*num_bones+2;
		}
		
		if (valid_lace)
		{
			write_lace_metapost(output, input_lace, mp_control, arc_end_labels, path_body_arc);
		}
	}
}

void read_lace_body(lace_body& body, string s)
{	
	char* line_buf = c_string(s);
	char* c1 = line_buf;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "read_lace_body: line_buf =  " << line_buf << endl;
	
	bool separator_found = false;
	int num_sides = 0;
	
	/* count the sides and prepare for reading the skeleton */
	for (unsigned int i=0; i< s.length(); i++)
	{
		/* move to start of next digit */
		while (!isdigit(*c1))
		{
			if (*c1 == '|')
			{
				separator_found = true;
				break;
			}
			c1++;
		}
		
		if (separator_found)
			break;
						
		/* find the end of the next number */
		char* c2 = c1;
		while (isdigit(*c2))
			c2++;
			
		num_sides++;
		
		c1 = c2;
	}

	if (separator_found == false)
	{
		cout << "Error!  Attempted to read lace code missing a '|' separator." << endl;
		exit(0);
	}
	
	if (num_sides %2 == 1)
	{
		cout << "Error!  Attempted to read lace code without a complete skeleton specification." << endl;
		exit(0);
	}
	
	body.skeleton = lace_skeleton(num_sides/2);
	read_lace_skeleton(body.skeleton,s);
	
	c1 = strchr(line_buf,'|');
	bool read_loop_locations = false;
	bool end_of_string = false;
	
	for (unsigned int i=0; i< s.length(); i++) // enough iterations but avoiding an infinite loop in the event of an error
	{
		int temp;
			
		/* move to start of next digit */
		while (!isdigit(*c1))
		{
			if (*c1 == 0)
			{
				end_of_string = true;
				break;
			}
			
			if (*c1 == ',')		
				read_loop_locations = true;
			c1++;
		}
		
		if (end_of_string)
			break;
			
		/* find the end of the next number */
		char* c2 = c1;
		while (isdigit(*c2))
			c2++;

		get_number(temp,c1);
		
		if (read_loop_locations)
		{
			body.interior_1_cell_count.push_back(temp);
			body.num_interior_1_cells++;
		}
		else
		{
			body.loop_location.push_back(temp);
			body.num_loops++;
		}
		
		c1 = c2;
	}
		
	body.num_boundary_1_cells = body.num_loops+1;
	
	if (body.num_loops == 1)
	{
		body.num_interior_1_cells = 1;
		body.num_2_cells = 2;
	}
	else
	{
		body.num_interior_1_cells = body.num_loops-2;
		body.num_2_cells = body.num_loops-1;
	}

	body.num_boundary_points = body.skeleton.num_cudgels;
	for (int i=0; i< body.skeleton.num_cudgels; i++)
		body.num_boundary_points = body.num_boundary_points + body.skeleton.bones[i][0] + body.skeleton.bones[i][1];
	body.num_boundary_points *= 2;
	
    delete[] line_buf;
}

void read_lace_skeleton(lace_skeleton& skeleton,string s)
{	
	char* line_buf = c_string(s);
	char* c1 = line_buf;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "read_lace_skeleton: line_buf =  " << line_buf << endl;
	debug << "read_lace_skeleton: num_cudgels = " << skeleton.num_cudgels << endl;
}
	
	for (int i=0; i< 2*skeleton.num_cudgels; i++)
	{
		/* move to start of next digit */
		while (!isdigit(*c1))
			c1++;
			
		/* find the end of the next number */
		char* c2 = c1;
		while (isdigit(*c2))
			c2++;
		
		get_number(skeleton.bones[i/2][i%2],c1);
		c1 = c2;
	}
	
	c1 = strchr(line_buf,'|');
	c1++;
	for (int i=0; i< skeleton.num_cudgels; i++)
	{
		/* move to start of next digit */
		while (!isalpha(*c1))
			c1++;

		if (*c1 == 'L')
			skeleton.head_polarity[i] = lace_control::polarity::LEFT;
		else
			skeleton.head_polarity[i] = lace_control::polarity::RIGHT;

		c1++;
	}
	
    delete[] line_buf;
}

ostream& operator << (ostream& os, const lace_body& body)
{
	os << body.skeleton << ' ';
	for (unsigned int i=0; i< body.loop_location.size(); i++)
		os << body.loop_location[i] << (i< body.loop_location.size()-1? " " : ", ");
	for (unsigned int i=0; i< body.interior_1_cell_count.size(); i++)
	{
		os << body.interior_1_cell_count[i];
		if (i< body.interior_1_cell_count.size() - 1)
			os << ' ';
	}
	
	return os;
}

void print_lace_body (ostream& os, const lace_body& body, string prefix)
{
	os << prefix << "num_cudgels = " << body.skeleton.num_cudgels << endl;
	os << prefix << "bones: ";
	for (int i=0; i< body.skeleton.num_cudgels; i++)
		os << '(' << body.skeleton.bones[i][0] << ',' << body.skeleton.bones[i][1] << ")" << ' ';
	os << endl;
	os << prefix << "head polarity: ";
	for (int i=0; i< body.skeleton.num_cudgels; i++)
		os << (body.skeleton.head_polarity[i] == lace_control::polarity::LEFT? 'L':'R');
	os << endl;
	os << prefix << "num_loops = " << body.num_loops << endl;
	os << prefix << "loop locations: ";
	for (unsigned int i=0; i< body.loop_location.size(); i++)
		os << body.loop_location[i] << " ";
	os << endl;
	os << prefix << "num_interior_1_cells = " << body.num_interior_1_cells << endl;
	os << prefix << "interior_1_cell_count: ";
	for (unsigned int i=0; i< body.interior_1_cell_count.size(); i++)
		os << body.interior_1_cell_count[i] << ' ';	
	os << endl;
	os << prefix << "num_boundary_1_cells = " << body.num_boundary_1_cells << endl;
	os << prefix << "num_2_cells = " << body.num_2_cells << endl;
	os << prefix << "num_boundary_points = " << body.num_boundary_points << endl;
}

ostream& operator << (ostream& os, const lace_skeleton& skeleton)
{
	for (int i=0; i< skeleton.num_cudgels; i++)
	{
		os << '(' << skeleton.bones[i][0] << ',' << skeleton.bones[i][1] << ')';
	}
	os << '|';
	
	for (int i=0; i< skeleton.num_cudgels; i++)
	{
		if (skeleton.head_polarity[i] == lace_control::polarity::LEFT)
			os << 'L';
		else
			os << 'R';
	}
	
	return os;
}

bool valid_lace(lace_body& input_lace)
{	
	/* check we've been given the correct number of interior 1-cell counts */
	if (input_lace.num_interior_1_cells != static_cast<int>(input_lace.interior_1_cell_count.size()))
	{
		cout << "\nInvalid lace input, expected " << input_lace.num_interior_1_cells << " interior 1-cell counts." << endl;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "valid_lace:   input_lace.num_interior_1_cells = " << input_lace.num_interior_1_cells << " not equal to interior_1_cell_count.size() = " << input_lace.interior_1_cell_count.size()<< endl;
		
		return false;
	}

	/* Check the pattern condition for each of the 2-cells.  Since we require that the body respect the loops
	   we are currently considering, the other points on the 1-cells must respect the pattern condition within
	   themselves.  Thus, we reduce the 1-cell counts to accommodate the loops.
	*/
	
	vector<int> boundary_1_cell_count(input_lace.num_loops+1);
	int offset = 0; 
	int index = 0;
	for (int i=0; i< input_lace.num_loops; i++)
	{
		boundary_1_cell_count[index++] = input_lace.loop_location[i] - offset;
		offset = input_lace.loop_location[i];
	}
	boundary_1_cell_count[index] = input_lace.num_boundary_points - offset;

	bool valid = true;
	for (int i=0; i<input_lace.num_2_cells; i++)
	{							
		int one_cell_1_count;
		int one_cell_2_count;
		int one_cell_3_count;
		
		if (input_lace.num_2_cells == 1) // num_loops == 2
		{
			one_cell_1_count = boundary_1_cell_count[0]-1;
			one_cell_2_count = boundary_1_cell_count[1]-2;
			one_cell_3_count = boundary_1_cell_count[2]-1;
		}
		else if (i==0 && input_lace.num_loops == 1)
		{
			one_cell_1_count = 0;
			one_cell_2_count = boundary_1_cell_count[0]-1;
			one_cell_3_count = input_lace.interior_1_cell_count[0]-1;
		}
		else if (i == 0)
		{
			one_cell_1_count = boundary_1_cell_count[0]-1;
			one_cell_2_count = boundary_1_cell_count[1]-2;
			one_cell_3_count = input_lace.interior_1_cell_count[0]-1;
		}
		else if (i==1 && input_lace.num_loops == 1)
		{
			one_cell_1_count = input_lace.interior_1_cell_count[0]-1;
			one_cell_2_count = boundary_1_cell_count[1]-1;
			one_cell_3_count = 0;
		}
		else if (i == input_lace.num_2_cells-1)
		{
			one_cell_1_count = input_lace.interior_1_cell_count[i-1]-1;
			one_cell_2_count = boundary_1_cell_count[i+1]-2;
			one_cell_3_count = boundary_1_cell_count[i+2]-1;
		}
		else
		{
			one_cell_1_count = input_lace.interior_1_cell_count[i-1]-1;
			one_cell_2_count = boundary_1_cell_count[i+1]-2;
			one_cell_3_count = input_lace.interior_1_cell_count[i]-1;
		}
											
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "valid_lace:   2-cell " << i << ", one_cell_1_count = " << one_cell_1_count << ", one_cell_2_count = " << one_cell_2_count << ", one_cell_3_count = " << one_cell_3_count << endl;
	
		int two_cell_boundary_points = one_cell_1_count + one_cell_2_count + one_cell_3_count; // if this is zero all 1-cell counts are zero
		
		if (two_cell_boundary_points %2 == 1)
		{
			cout << "\nInvalid lace input, 2-cell " << i << " does not have an even number of points in its boundary." << endl;
			valid = false;
							
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "valid_lace: invalid interior 1-cell options, 2-cell " << i << " does not have an even number of points in its boundary." << endl;
			break;
		}
		
		int one_cell_threshold = two_cell_boundary_points/2;
		
		if (one_cell_1_count > one_cell_threshold || one_cell_2_count > one_cell_threshold || one_cell_3_count > one_cell_threshold)
		{
			cout << "\nInvalid lace input, 2-cell " << i << " does not satisfy the pattern condition." << endl;
			valid = false;
							
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "valid_lace: invalid interior 1-cell options, 2-cell " << i << " does not satisfy the pattern condition." << endl;
			break;
		}		
	}			
	
	return valid;
}

void assign_pairs(matrix<int>& edge_pairs, int& pair_marker, int start_label_1, int start_label_2, int start_label_3, int count_1, int count_2, int count_3, bool boundary_place_3)
{

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "assign_pairs: start_label_1 = " << start_label_1 << ", count_1 = " << count_1 << endl;
	debug << "assign_pairs: start_label_2 = " << start_label_2 << ", count_2 = " << count_2 << endl;
	debug << "assign_pairs: start_label_3 = " << start_label_3 << ", count_3 = " << count_3 << endl;
	debug << "assign_pairs: boundary_place_3 = " << boundary_place_3 << endl;
}
								
	int num_2_cell_edges = (count_1 + count_2 + count_3)/2;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "assign_pairs: num_2_cell_edges = " << num_2_cell_edges << endl;
								
	/* Place 1. */
	int num_edges = num_2_cell_edges - count_2;      
								
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "assign_pairs: place 1 num_edges = " << num_edges << endl;
	
	for (int j=1; j<= num_edges; j++)
	{
		edge_pairs[pair_marker][0] = start_label_1+num_edges-j;
		if (boundary_place_3)
			edge_pairs[pair_marker][1] = start_label_3+count_3-num_edges-1+j;
		else
			edge_pairs[pair_marker][1] = start_label_3+num_edges-j;
									
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "assign_pairs:   pair " << edge_pairs[pair_marker][0] << ' ' << edge_pairs[pair_marker][1] << endl;
	
		pair_marker++;
	}

	
	/* Place 2. */
	num_edges = num_2_cell_edges - count_3;      
								
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "assign_pairs: place 2 num_edges = " << num_edges << endl;

	for (int j=1; j<= num_edges; j++)
	{
		edge_pairs[pair_marker][0] = start_label_1+count_1-num_edges-1+j;
		edge_pairs[pair_marker][1] = start_label_2+num_edges-j;
									
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "assign_pairs:   pair " << edge_pairs[pair_marker][0] << ' ' << edge_pairs[pair_marker][1] << endl;
	
		pair_marker++;
	}
					
	/* Place 3. */
	num_edges = num_2_cell_edges - count_1;      

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "assign_pairs: place 3 num_edges = " << num_edges << endl;

	for (int j=1; j<= num_edges; j++)
	{
		edge_pairs[pair_marker][0] = start_label_2+count_2-num_edges-1+j;
		
		if (boundary_place_3)
			edge_pairs[pair_marker][1] = start_label_3+num_edges-j;
		else
			edge_pairs[pair_marker][1] = start_label_3+count_3-num_edges-1+j;
									
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "assign_pairs:   pair " << edge_pairs[pair_marker][0] << ' ' << edge_pairs[pair_marker][1] << endl;
	
		pair_marker++;
	}				
}
				
void write_lace_metapost(ofstream& os, lace_body& input_lace, metapost_control& mp_control, matrix<int>& arc_end_labels,matrix<int>& path_body_arc)
{
	int num_boundary_points = input_lace.num_boundary_points;
	int num_loops = input_lace.num_loops;
	vector<int>& used_labels = input_lace.loop_location;
	lace_skeleton& skeleton = input_lace.skeleton;
	
	int num_disc_arc_components = arc_end_labels.numrows();
	
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
    debug << "write_lace_metapost: mp_control: " << endl;
    print (mp_control, debug, "write_lace_metapost:   ");
}
	
	int num_cudgels = skeleton.num_cudgels;
	int max_num_lower_bones = 0;

	for (int i=0; i< num_cudgels; i++)
	{
		if (skeleton.bones[i][0] > max_num_lower_bones)
			max_num_lower_bones = skeleton.bones[i][0];
	}	

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: max_num_lower_bones = " << max_num_lower_bones << endl;
	
	vector<int> side_count(2*num_cudgels);
	for (int i=0; i< num_cudgels; i++)
	{
		side_count[2*i+1] = side_count[2*i] = skeleton.bones[i][0]+skeleton.bones[i][1];
		if (skeleton.head_polarity[i] == lace_control::polarity::LEFT)
			side_count[2*i]++;
		else
			side_count[2*i+1]++;		
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: side_count: ";
	for (int i=0; i< 2*num_cudgels; i++)
		debug << side_count[i] << ' ';
	debug << endl;
}

	/* side_start indicates the label before the first label on a side, so that it may be
	   used to calaulate offsets
	*/
	vector<int> side_start(2*num_cudgels);
	side_start[1] = side_count[0]+1;
	
	for (int i=1; i< num_cudgels; i++)
	{
		side_start[2*i] = side_start[2*i-1]+side_count[2*i-1];
		side_start[2*i+1] = side_start[2*i]+side_count[2*i]+1;
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: side_start: ";
	for (int i=0; i< 2*num_cudgels; i++)
		debug << side_start[i] << ' ';
	debug << endl;
}

	os << "\n\n";
	
//	if (title.length())
//		os << "% " << title << endl;

	
	os << "% " << input_lace << endl;
	
	/* controls for coordinate output */
	int num_decimal_points = 3;
//	int output_field_width = 12;
	os.setf(ios::fixed,ios::floatfield);
	os.precision(num_decimal_points);

	os << "\nbeginfig(fignum);" << endl;
	os << "fignum:=fignum+1;" << endl;
	os << "numeric u,h,v;" << endl;
	os << "u=" << mp_control.unit_size/10 << "pt;" << endl;
	os << "h=" << mp_control.horizontal_units << "u;" << endl;
	os << "v=" << mp_control.vertical_units << "u;" << endl;
	os << "path p[];" << endl;
	os << "pickup pencircle scaled " << mp_control.pen_size*0.5 << "pt;" << endl;

	/* identify the maximum number of tiers (bones plus head plus tail) on each cudgel */
	int max_num_tiers = 0;
	for (int i=0; i< 2*num_cudgels; i++)
	{
		if (side_count[i] > max_num_tiers)
			max_num_tiers = side_count[i];
	}	
	
	max_num_tiers++; // allow for the tail

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: max_num_tiers =  " << max_num_tiers << endl;
	
	/* identify whether a label lies above or below the head */
	vector<int> upper_bone_label(num_boundary_points);
	
	for (int i=0; i< num_cudgels; i++)
	{
		for (int j=0; j < skeleton.bones[i][1]; j++)
		{
			upper_bone_label[side_start[2*i]+j] = 1;
			upper_bone_label[side_start[2*i+1]+side_count[2*i+1]-1-j] = 1;
		}
	}	
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: upper_bone_label: " << endl;
	for (int i=0; i< num_boundary_points; i++)
		debug << upper_bone_label[i] << ' ';
	debug << endl;
}

	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: used_labels: " << endl;
	for (int i=0; i< num_loops; i++)
		debug << used_labels[i] << ' ';
	debug << endl;
}

	 /*Identify the location of each loop and whether the side of a cudgel carries a loop.  There is at most one loop on each cudgel side 
	   and at most one in between cudgels.  We therefore record the cudgel side (as positive numbers) or the cudgel after the loop (as 
	   negative numbers).
	*/
	vector<int> loop_location(num_loops);

	vector<int> side_loop(2*num_cudgels);
	
	for (int i=0; i< num_loops; i++)
	{
		int start = 0;
		for (int j=0; j< num_cudgels; j++)
		{
			if (used_labels[i] <= start+side_count[2*j])
			{
				loop_location[i] = 2*j;
				side_loop[2*j] = used_labels[i];
				break;
			}
			else if (used_labels[i] < start+side_count[2*j]+side_count[2*j+1]+1)
			{
				loop_location[i] = 2*j+1;
				side_loop[2*j+1] = used_labels[i];
				break;
			}
			else if (used_labels[i] == start+side_count[2*j]+side_count[2*j+1]+1)			
			{
				loop_location[i] = (j+1)*-1;
				break;
			}
			else
			{
				start += side_count[2*j]+side_count[2*j+1]+1;
			}
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: loop_location: ";
	for (int i=0; i< num_loops; i++)
		debug << loop_location[i] << ' ';
	debug << endl;
    debug << "write_lace_metapost: side_loop: ";
	for (int i=0; i< 2*num_cudgels; i++)
		debug << side_loop[i] << ' ';
	debug << endl;
}


	/* Determine the number of additional points needed below each cudgel required to  draw body arcs that 
	   connect to different cudgel sides.  We do this by considering the arc_end_labels: if the two ends of
	   a body disc arc lie on different cudgel sides, then we need an additional point below each intervening cudgel.
	   
	   Since the initial label is always smaller than the terminal label, we regard initial labels that are tails 
	   as residing on the right side of a cudgel and terminal labels that are tails as residing on the left side.
	   At the same time, record the location of each arc end label.  

	   Not only do we need additional points below the cudgel but also to the side.  If the side of a cudgel 
	   does not include the head but carries a loop then we need additional points on the horizontal line at the head, 
	   since all arcs connecting to that same side are concentric.  Moreover, if a label above such a set of loops 
	   connects to a different side, then it must also pass around the concentric loops, so we need a point on the head 
	   line for that path too, regardless of whether the arc passes under the tail or not.
	      
	   If the side of a cudgel including the head also carries a loop, then it must be a loop from the tail to 
	   the head.  In this case we need additional points on a horizontal line half way between the tail and the 
	   head, which will accommodate any arcs above the head on that side that need to pass under the tail.  We
	   count these points as "head points".
	   
	   If any labels below the set of concentric loops connect to a different side, they either loop around the 
	   tail of the cudgel, or they connect to one of the additional points below an adjacent cudgel.  In either 
	   case we have an additional point on a horizontal line at the tail.
	   
	*/
	vector<int> additional_points_below(num_cudgels);
	vector<int> end_label_side(num_boundary_points);

	vector<int> additional_head_points(2*num_cudgels);
	vector<int> additional_tail_points(2*num_cudgels);

	matrix<int> head_point_type_count(2*num_cudgels,3);
	matrix<int> tail_point_type_count(2*num_cudgels,4);
	

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: evaluate additional points" << endl;
	
	for (int i=0; i< num_disc_arc_components; i++)
	{
		int initial_label = arc_end_labels[i][0];
		int terminal_label = arc_end_labels[i][1];
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: disc arc " << i << ", initial_label " << initial_label << " terminal_label " << terminal_label << endl;
		
		int initial_side = -1;
		for (int j=0; j< num_cudgels; j++)
		{
			if (initial_label <= side_start[2*j]+side_count[2*j]) // don't include the tail for initial labels
			{
				initial_side = 2*j;
				end_label_side[initial_label-1] = 2*j;
				break;
			}
			else if (initial_label <= side_start[2*j+1]+side_count[2*j+1]) // initial tails reside on the right
			{
				initial_side = 2*j+1;
				end_label_side[initial_label-1] = 2*j+1;
				break;
			}
		}
		
		int terminal_side = -1;
		for (int j=0; j< num_cudgels; j++)
		{
			if (terminal_label <= side_start[2*j]+side_count[2*j]+1) // include the tail for terminal labels
			{
				terminal_side = 2*j;
				end_label_side[terminal_label-1] = 2*j;
				break;
			}
			else if (terminal_label <= side_start[2*j+1]+side_count[2*j+1])
			{
				terminal_side = 2*j+1;
				end_label_side[terminal_label-1] = 2*j+1;
				break;
			}
		}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: initial label on side " << initial_side << endl;
    debug << "write_lace_metapost: terminal label on side " << terminal_side << endl;
}		

		if (terminal_side == initial_side)
		{	
			head_point_type_count[initial_side][(initial_side%2==1? ODD_HEAD_L:EVEN_HEAD_L)]++;
			additional_head_points[initial_side]++;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: loop head point on side " << initial_side << endl;
		}
		else if (initial_side % 2 == 0 || (initial_side % 2 == 1 && terminal_side != initial_side+1))
		{
			if (upper_bone_label[initial_label-1])
			{
				head_point_type_count[initial_side][(initial_side%2==1? ODD_HEAD_O:EVEN_HEAD_O)]++;
				additional_head_points[initial_side]++;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: upper originating head point on side " << initial_side << endl;
    
				if (initial_side % 2 == 0 || (initial_side%2 ==1 && mp_control.right_originating_tail_points))
				{
					tail_point_type_count[initial_side][(initial_side%2==1? ODD_TAIL_UO:EVEN_TAIL_UO)]++;
					additional_tail_points[initial_side]++;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: upper originating tail point on side " << initial_side << endl;
				}
    
			}
			else if (initial_side % 2 == 0 || (initial_side % 2 == 1 && mp_control.right_originating_tail_points))
			{
				/* includes arc originating at tails on odd sides */
				tail_point_type_count[initial_side][(initial_side%2==1? ODD_TAIL_LO:EVEN_TAIL_LO)]++;
				additional_tail_points[initial_side]++;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: lower originating tail point on side " << initial_side << endl;
			}

				
			if (upper_bone_label[terminal_label-1])
			{
				head_point_type_count[terminal_side][(terminal_side%2==1? ODD_HEAD_T:EVEN_HEAD_T)]++;
				additional_head_points[terminal_side]++;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: upper terminating head point on side " << terminal_side << endl;
    
				if (terminal_side % 2 == 1 || (terminal_side % 2 == 0 && mp_control.left_terminating_tail_points))
				{
					tail_point_type_count[terminal_side][(terminal_side%2==1? ODD_TAIL_UT:EVEN_TAIL_UT)]++;
					additional_tail_points[terminal_side]++;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: upper terminating tail point on side " << terminal_side << endl;
				}
    
			}			
			else if (terminal_side % 2 == 1 || (terminal_side % 2 == 0 && mp_control.left_terminating_tail_points && terminal_label < side_start[terminal_side]+side_count[terminal_side]+1 ))  // i.e not a tail on an even side
			{
				tail_point_type_count[terminal_side][(terminal_side%2==1? ODD_TAIL_LT:EVEN_TAIL_LT)]++;
				additional_tail_points[terminal_side]++;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: lower terminating tail point on side " << terminal_side << endl;
			}				
			
			int first_cudgel = (initial_side%2==1?(initial_side+1)/2: initial_side/2);
			int last_cudgel = (terminal_side%2==1?terminal_side/2: (terminal_side-1)/2);
			for (int c=first_cudgel; c <= last_cudgel; c++)
				additional_points_below[c]++;
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: head_point_type_count" << endl;
	for (int i=0; i< 2* num_cudgels; i++)
	{
		debug << "write_lace_metapost: side " << i;
		if (i%2==1)
			debug << " L = " << head_point_type_count[i][ODD_HEAD_L] << " T = " << head_point_type_count[i][ODD_HEAD_T] << " O = " << head_point_type_count[i][ODD_HEAD_L] << endl;
		else
			debug << " L = " << head_point_type_count[i][EVEN_HEAD_L] << " T = " << head_point_type_count[i][EVEN_HEAD_T] << " O = " << head_point_type_count[i][EVEN_HEAD_O] << endl;
	}
    debug << "write_lace_metapost: tail_point_type_count" << endl;
	for (int i=0; i< 2* num_cudgels; i++)
	{
		debug << "write_lace_metapost:   side " << i;
		if (i%2==1)
			debug << " LT = " << tail_point_type_count[i][ODD_TAIL_LT] << " LO = " << tail_point_type_count[i][ODD_TAIL_LO] << " UT = " << tail_point_type_count[i][ODD_TAIL_UT] << " UO = " << tail_point_type_count[i][ODD_TAIL_UO] << endl;
		else
			debug << " LT = " << tail_point_type_count[i][EVEN_TAIL_LT] << " LO = " << tail_point_type_count[i][EVEN_TAIL_LO] << " UT = " << tail_point_type_count[i][EVEN_TAIL_UT] << " UO = " << tail_point_type_count[i][EVEN_TAIL_UO] << endl;
	}
}	
	int max_num_additional_points = *max_element(additional_points_below.begin(), additional_points_below.end());
	int max_num_head_points = *max_element(additional_head_points.begin(), additional_head_points.end());
	int max_num_tail_points = *max_element(additional_tail_points.begin(), additional_tail_points.end());
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: additional points at the head ";
	for (int i=0; i< 2* num_cudgels; i++)
		debug << additional_head_points[i] << ' ';
	debug << endl;
    debug << "write_lace_metapost: additional points at the tail ";
	for (int i=0; i< 2* num_cudgels; i++)
		debug << additional_tail_points[i] << ' ';
	debug << endl;
    debug << "write_lace_metapost: additional points below cudgels ";
	for (int i=0; i< num_cudgels; i++)
		debug << additional_points_below[i] << ' ';
	debug << endl;
    debug << "write_lace_metapost: max_num_additional_points = " <<  max_num_additional_points << endl;
    debug << "write_lace_metapost: max_num_head_points = " <<  max_num_head_points << endl;
    debug << "write_lace_metapost: max_num_tail_points = " <<  max_num_tail_points << endl;
}
	

	/* work out how far apart the cudgels need to be, we have to accommodate additional head and tail points on facing sides */
	int cudgel_spacing = 0;
	for (int i=1; i< num_cudgels-1; i++)
	{
		int space_required = max(additional_head_points[2*i-1]+additional_head_points[2*i],additional_tail_points[2*i-1]+additional_tail_points[2*i])+1;
		if (space_required > cudgel_spacing)
			cudgel_spacing = space_required;
	}
	
	/* check the two outer sides to give a more balanced picture */
	if (max(additional_head_points[0],additional_tail_points[0]) > cudgel_spacing)
		cudgel_spacing = max(additional_head_points[0],additional_tail_points[0]);
	
	if (max(additional_head_points[2*num_cudgels-1],additional_tail_points[2*num_cudgels-1]) > cudgel_spacing)
		cudgel_spacing = max(additional_head_points[2*num_cudgels-1],additional_tail_points[2*num_cudgels-1]);

	/* Set out the coordinates of the z points required to draw the lace.  We regard the head as lying on the line y=0 so that we may assign points 
	   evenly arount the head more easily. The cudgels are spaced max_num_concentric_loops+max_num_additional_points units apart.
	   At the same time, record the initial z_coordinate associated with each label.
	*/
	
	cudgel_spacing *= 1.3;
	cudgel_spacing *= mp_control.cudgel_space;
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost: cudgel_spacing = " << cudgel_spacing << endl;
	
	int cudgel_x = 0;
	int z_index = 1;
	
	/* label_z_coordinate records the active z coordinate associated with each label.  If we need to extend paths between cudgels, we add loops and 
	   segments to extend paths to the additional points and record in label_z_coordinate the last z_coordinate to which a path from each label has
	   been extended.
	*/
	vector<int> label_z_coordinate(num_boundary_points); 
	vector<int> path_endpoint_coordinate(2*num_cudgels);
	
	for (int i=0; i< num_cudgels; i++)
	{
		bool left_head = (skeleton.head_polarity[i] == lace_control::polarity::LEFT);
		
		/* work down the left side of the cudgel */
		for (int j=0; j < skeleton.bones[i][1]; j++)
		{
			os << "z" << z_index << "=(" << cudgel_x << "h," << skeleton.bones[i][1]-j << "v);" << endl;
			label_z_coordinate[side_start[2*i]+j] = z_index;
			label_z_coordinate[side_start[2*i+1]+side_count[2*i+1]-1-j] = z_index;
			z_index++;
		}
		
		/* now the head */
		os << "z" << z_index << "=(" << cudgel_x << "h,0);" << endl;
		path_endpoint_coordinate[2*i] = z_index;
		
		if (left_head)
			label_z_coordinate[side_start[2*i]+skeleton.bones[i][1]] = z_index;
		else
			label_z_coordinate[side_start[2*i+1]+skeleton.bones[i][0]] = z_index;
		
		z_index++;

		for (int j=0; j< skeleton.bones[i][0]; j++)
		{
			os << "z" << z_index << "=(" << cudgel_x << "h," << -(j+1) << "v);" << endl;			
			label_z_coordinate[side_start[2*i]+skeleton.bones[i][1]+(left_head?1:0)+j] = z_index;
			label_z_coordinate[side_start[2*i+1]+skeleton.bones[i][0]-j-1] = z_index;
			z_index++;
		}
		
		/* now the tail*/
		os << "z" << z_index << "=(" << cudgel_x << "h," << -(max_num_lower_bones+1) << "v);" << endl;
		label_z_coordinate[side_start[2*i+1]-1] = z_index;		
		path_endpoint_coordinate[2*i+1] = z_index;
		
		z_index++;
		
		cudgel_x += cudgel_spacing;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: label_z_coordinate ";
	for(int i=0; i<num_boundary_points; i++)
		debug << label_z_coordinate[i] << ' ';
	debug << endl;
}
	
	/* set up the z_coordinates for the additional points below the cudgels, we record in the rows of additional_point_coordinate 
	   the coordinate index associated with each additional point below the cudgels
	*/
	matrix<int> additional_point_coordinate(num_cudgels,max_num_additional_points);
	
	cudgel_x = 0;
	for (int i=0; i< num_cudgels; i++)
	{
		for (int j=0; j< additional_points_below[i]; j++)
		{
			os << "z" << z_index << "=(" << cudgel_x << "h," << -(max_num_lower_bones+j+2) << "v);" << endl;
			additional_point_coordinate[i][j] = z_index;
			z_index++;
		}
		
		cudgel_x += cudgel_spacing;
	}	

	/* set up the z_coordinates for the head and tail points taking note of cudgel sides that have loops involving the head 
	   record in the rows of head_point_coordinate the coordinate index of each additional head point for each cudgel side
	*/
	matrix<int> head_point_coordinate(2*num_cudgels,max_num_head_points);
	matrix<int> tail_point_coordinate(2*num_cudgels,max_num_tail_points);

	cudgel_x = 0;
	for (int i=0; i< 2*num_cudgels; i++)
	{
		float head_point_y_coord = 0;
				
		/* if this side of the current cudgel contains the head and it carries a loop then it must be loop connecting the tail to the head */
		if ( 
			(side_loop[i] != 0 && ((i%2 == 0 && skeleton.head_polarity[i/2] == lace_control::polarity::LEFT) || (i%2 == 1 && skeleton.head_polarity[i/2] == lace_control::polarity::RIGHT)))
		   )
			head_point_y_coord = -1*(static_cast<float>(max_num_lower_bones+1))/2;

		if (i%2 ==1)
		{
			set_out_tail_coordinates(i, cudgel_x, max_num_lower_bones, z_index, head_point_type_count, head_point_coordinate, tail_point_type_count, tail_point_coordinate, os);
			set_out_head_coordinates(i, cudgel_x, head_point_y_coord, z_index, head_point_type_count, head_point_coordinate, tail_point_type_count, tail_point_coordinate, os);
		}
		else
		{
			set_out_head_coordinates(i, cudgel_x, head_point_y_coord, z_index, head_point_type_count, head_point_coordinate, tail_point_type_count, tail_point_coordinate, os);
			set_out_tail_coordinates(i, cudgel_x, max_num_lower_bones, z_index,  head_point_type_count, head_point_coordinate, tail_point_type_count, tail_point_coordinate, os);
		}

		if (i%2 ==1)
			cudgel_x += cudgel_spacing;
	}	

	vector<int> upper_bone_coordinate(z_index-1);
	for (int i=0; i< num_cudgels; i++)
	{
		/* work down the left side of the cudgel */
		for (int j=0; j < skeleton.bones[i][1]; j++)
			upper_bone_coordinate[label_z_coordinate[side_start[2*i]+j]] = 1;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: upper bone coordinates ";
	for(int i=0; i< z_index-1; i++)
	{
		if (upper_bone_coordinate[i] != 0)
			debug << i << ' ';
	}
	debug << endl;
}
	
	if (mp_control.draw_lace_frame)
	{
		for (int i=0; i< z_index-1; i++)
		{
			os << "fill fullcircle scaled 1mm shifted z" << i+1 << ";" << endl;
			os << "label.top (btex "; 
			if (true || mp_control.scriptscript_labels)
				os << "\\fiverm ";
			else if (mp_control.script_labels)
				os << "\\sevenrm ";
		    os << i+1 << " etex, z" << i+1 << ");" << endl;
		}
	}
	else
	{
		for (int i=0; i< 2*num_cudgels; i++)
		{
			os << "fill fullcircle scaled 1mm shifted z" << path_endpoint_coordinate[i] << ";" << endl;
		}
	}
	
	/* set up markers to hold the index into additional_point_coordinate, head_point_coordinate and tail_point_coordinate for the 
	   next additional point we should be using.  We work in towards or out away from a cudgel, depending on the nature of the path.
	   
	   Body disc arcs are oriented implicitly from arc_end_labels[i][0] to arc_end_labels[i][1].  That is, from the first end label 
	   encountered as we follow the boundary of the disc anticlockwise from the basepoint.  When attached to a skeleton the disc arcs
	   are oriented from "left to right": from an originating boundary label having a lower value than the corresponding terminating label.
	   
	   Body disc arc do not cross.  If the arcs are a_i, with originating and terminating boundary points u_i and v_i respectively, then
	   for any two arcs a_1 and a_2, if one of u_2, v_2 lies between u_1 and v_1 when moving anticlockwise from u_1 to v_1, both u_2 and 
	   v_2 must lie between u_1 and v_1.
	   
	   The above two observations, together with the fact that each side of a cudgel contains at most one (bigon) loop, show that at the 
	   bottom of the right side of a cudgel (below any arcs whose originating and terminating endpoints both lie on that side), we cannot 
	   have an originating arc below a terminating arc.  Similarly, on the left side of a cudgel, we cannot have a terminating arc below an
	   originating arc.  See the left two diagrams below.
	   
	                                        
	                         |  |                      |  |                +-------          -------+
	                         |  |-<-----+      +-----<-|  |                |-<--?              ?--<-|
	                         |  |-->--? |      | ? -->-|  |                |-->-----        ------>-|
	                         +--+       |      |       +--+                |                        |
	                               -----+      +----	                               
	   
	   At the top of a cudgel (above any arcs whose originating and terminating endpoints both lie on that side), the same restrictions 
	   apply on the right and left side of a cudgel as at the bottom. See the right two diagrams above.
	   
	   The boundary points on the side of a cudgel are therefore grouped into five (possibly empty) sets: lower terminaing (LT), 
	   lower originating (LO), loops (L), upper terminating (UT) and upper originating (UO), arranged as follows on each side of a 
	   cudgel
	   
	                          +-------------+  
	                       |  |UO         UT|  |
	                       |  |UT         UO|  | 
                           |  |L           L|  |
                           |  |LO         LT|  |
                           |  |LT         LO|  |
                           +--+             +--+

	   Note that by similar arguments to the above, on the right side of a cudgel at most one of LO and UT may be non-empty and on
	   the left side, at most one of LT and UO.
	   
	   For these purposes, the cudgel tail is regarded as part of the left (even) side of a cudgel if it originates a body disc arc
	   and on the right (odd) side if it terminates a disc arc.  When we determined the number of additional points required below 
	   each cudgel, we took the view that "Since the initial label is always smaller than the terminal label, we regard initial 
	   labels that are tails as residing on the right side of a cudgel and terminal labels that are tails as residing on the left side."
	   Thus, when a cudgel tail originates an arc, the corresponding right side of the cudgel has LT the empty set and when a cudgel tail
	   terminates an arc, that side has LO the empty set. 
	   
	   The coordinates for the additional points are recorded in the rows of head_point_coordinate and tail_point_coordinate in the order
	   corresponding to moving outwards from the cudgel.  However, as we consider the body disc arcs in order, we have to pick up these
 	   additional points in the correct order to reflect the relative position of parallel originating or terminating arcs on a given side 
	   of a cudgel.  The order for the tail points is as shown in the following diagrams:
	   
				   |  |----------------UO-+                 +-UT----------------|  |
				   |  |-----------UT-+    |                 |    +-UO-----------|  | 
				   |  |L             |    |                 |    |             L|  |
				   |  |------LO-+    |    |                 |    |    +-LT------|  |
				   |  |-LT-+    |    |    |                 |    |    |    +-LO-|  |
				   +--+----|----|----|----|----->     <-----|----|----|----|----+--+  
	                   <--   --> <--  -->                    <--  -->  <--  -->
	   
	   The order for the head points is as follows:

				   |  |-----------O--+                      +-T------------|  | 
				   |  |------T--+    |                      |    +-O-------|  |
				   |  |-L--+    |    |                      |    |    +-L--|  |
				   +--+----|----|----|----->          <-----|----|----|----+--+
	               |  |<--  <--  -->                         <--   -->  -->	   	  
	        
	   To accommodate the above sets of additional head and tail points, for each row of tail_point_coordinates we maintain four markers 
	   for each row of head_point_coordinates, three markers.  On odd sides, the markers for terminating arcs decrement towards the start 
	   of the coordinate rows and those for originating arcs increment towards the end.  For even sides the oposite is true.  Loops 
	   connecting the same side decrement towards the start.

	   Note that in the head and tail coordinate rows for the left and right sides the additional points corresponding to the sets 
	   LT,LO,UT & UO (for tails) and L, O & T (for heads) appear in a different order: e.g. for right side tails they are recorded in the 
	   order LT LO UT UO and for left side tails in the order LO LT UO UT

	   Before we can initialize the markers, we have to count the number of boundary arcs of each type on the cudgel sides.
	   
	   Below the cudgels, all the body disc arcs are oriented left to right and so we always work inwards towards the cudgel when picking
	   up additional points below.
	   
	   As an aside, note that the body disc is attached to the skeleton disc as the southern hemisphere of a 2-sphere, the skeleton disc being the
	   northern hemisphere, viewed from above the north pole.  With this perspective, when considered by itself, the body disc is viewed
	   from within the 2-sphere.  The final lace diagram lies in the plane and is obtained by an isotopy of S^2 that slides the body disc 
	   arcs around to the front of the 2-sphere.  This isotopy has the effect of reversing the notions of left and right when comparing 
	   the body disc by itself and the final lace diagram in the plane.
	*/
	vector<int> additional_point_marker(num_cudgels);
	for (int i=0; i< num_cudgels; i++)
		additional_point_marker[i] = additional_points_below[i]-1;

	matrix<int> head_point_marker(2*num_cudgels,3);
	for (int i=0; i< 2*num_cudgels; i++)
	{
		/* for both odd and even sides, the first two markers decrement towards the start of the row and the third increments towards the end */		
		head_point_marker[i][0] = head_point_type_count[i][0]-1; 
		head_point_marker[i][1] = head_point_type_count[i][0]+head_point_type_count[i][1]-1; 
		head_point_marker[i][2] = head_point_type_count[i][0]+head_point_type_count[i][1]; 
		
		if (i%2 == 1)
				{
			if (head_point_type_count[i][ODD_HEAD_L] == 0)
				head_point_marker[i][ODD_HEAD_L] = -1;

			if (head_point_type_count[i][ODD_HEAD_O] == 0)
				head_point_marker[i][ODD_HEAD_O] = -1;

			if (head_point_type_count[i][ODD_HEAD_T] == 0)
				head_point_marker[i][ODD_HEAD_T] = -1;
		}		
		else
		{
			if (head_point_type_count[i][EVEN_HEAD_L] == 0)
				head_point_marker[i][EVEN_HEAD_L] = -1;

			if (head_point_type_count[i][EVEN_HEAD_O] == 0)
				head_point_marker[i][EVEN_HEAD_O] = -1;

			if (head_point_type_count[i][EVEN_HEAD_T] == 0)
				head_point_marker[i][EVEN_HEAD_T] = -1;
		}		
	}	
	
	matrix<int> tail_point_marker(2*num_cudgels,4);
	for (int i=0; i< 2*num_cudgels; i++)
	{
		/* for both odd and even sides, the first and third markers decrement towards the start of the row and the 
		   second and fourth increment towards the end */		
		tail_point_marker[i][0] = tail_point_type_count[i][0]-1; 
		tail_point_marker[i][1] = tail_point_type_count[i][0]; 
		tail_point_marker[i][2] = tail_point_type_count[i][0]+tail_point_type_count[i][1]+tail_point_type_count[i][2]-1; 
		tail_point_marker[i][3] = tail_point_type_count[i][0]+tail_point_type_count[i][1]+tail_point_type_count[i][2]; 
		
		if (i%2 == 1)
		{
			if (tail_point_type_count[i][ODD_TAIL_LT] == 0)
				tail_point_marker[i][ODD_TAIL_LT] = -1;

			if (tail_point_type_count[i][ODD_TAIL_LO] == 0)
				tail_point_marker[i][ODD_TAIL_LO] = -1;

			if (tail_point_type_count[i][ODD_TAIL_UT] == 0)
				tail_point_marker[i][ODD_TAIL_UT] = -1;
				
			if (tail_point_type_count[i][ODD_TAIL_UO] == 0)
				tail_point_marker[i][ODD_TAIL_UO] = -1;
		}
		else
		{
			if (tail_point_type_count[i][EVEN_TAIL_LT] == 0)
				tail_point_marker[i][EVEN_TAIL_LT] = -1;

			if (tail_point_type_count[i][EVEN_TAIL_LO] == 0)
				tail_point_marker[i][EVEN_TAIL_LO] = -1;

			if (tail_point_type_count[i][EVEN_TAIL_UT] == 0)
				tail_point_marker[i][EVEN_TAIL_UT] = -1;
				
			if (tail_point_type_count[i][EVEN_TAIL_UO] == 0)
				tail_point_marker[i][EVEN_TAIL_UO] = -1;
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: initial additional_point_marker ";
	for (int i=0; i< num_cudgels; i++)
		debug << additional_point_marker[i] << ' ';
	debug << endl;
	report_markers(head_point_marker,head_point_coordinate);
	report_markers(tail_point_marker,tail_point_coordinate);
	
}	

	/* we build up a vector of z coordinate indices for each body disc arc.  Simple arcs between adjacent cudgels just have two
	   indices, indicating the indices corresponding to the arc end labels.  More complicated arcs involve the indices of the 
	   additional points through which the path passes.  
	   
	   We store the indices in a vector of vectors, since we want to push_back the indices as we work through the disc arcs.
	*/
	vector< vector<int> >disc_arc_sequence(num_disc_arc_components);

	for (int i=0; i< num_disc_arc_components; i++)
	{
		int initial_label = arc_end_labels[i][0];
		int terminal_label = arc_end_labels[i][1];
		int initial_side = end_label_side[initial_label -1];
		int terminal_side = end_label_side[terminal_label -1];
		vector<int> disc_arc = {label_z_coordinate[initial_label-1]};
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost: disc arc " << i << ", initial_label " << initial_label << " terminal_label " << terminal_label << endl;
	debug << "write_lace_metapost: initial label on side " << initial_side << endl;
    debug << "write_lace_metapost: terminal label on side " << terminal_side << endl;
}		

		if (terminal_side == initial_side)
		{			  
			int marker = (initial_side%2==1? ODD_HEAD_L:EVEN_HEAD_L);
			int coordinate = head_point_coordinate[initial_side][head_point_marker[initial_side][marker]];
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add loop head point " << coordinate << " to disc_arc" << endl;
    
			disc_arc.push_back(coordinate);
			head_point_marker[initial_side][marker]--;						
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost:   head loop point marker for side " << initial_side << " reduced to " << head_point_marker[initial_side][marker] << endl;
	report_markers(head_point_marker,head_point_coordinate);
}
		}
		else if (initial_side % 2 == 0 || (initial_side % 2 == 1 && terminal_side != initial_side+1))
		{
			int coordinate; 
			int marker; 
			
			if (upper_bone_label[initial_label-1])
			{
				marker = (initial_side%2==1? ODD_HEAD_O:EVEN_HEAD_O);
				coordinate = head_point_coordinate[initial_side][head_point_marker[initial_side][marker]];
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add upper originating head point " << coordinate << " to disc_arc" << endl;
	
				disc_arc.push_back(coordinate);
				
				if (initial_side %2 == 1)
					head_point_marker[initial_side][marker]++;
				else
					head_point_marker[initial_side][marker]--;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost:   head_point_marker for side " << initial_side << " adjusted to " << head_point_marker[initial_side][marker] << endl;				
	report_markers(head_point_marker,head_point_coordinate);
}
    
				if (initial_side % 2 == 0 || (initial_side%2 ==1 && mp_control.right_originating_tail_points))
				{
					marker = (initial_side%2==1? ODD_TAIL_UO:EVEN_TAIL_UO);
					coordinate = tail_point_coordinate[initial_side][tail_point_marker[initial_side][marker]];
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add upper originating tail point " << coordinate << " to disc_arc" << endl;
					disc_arc.push_back(coordinate);
				
					if (initial_side %2 == 1)
						tail_point_marker[initial_side][marker]++;
					else
						tail_point_marker[initial_side][marker]--;
					
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost:   tail_point_marker for side " << initial_side << " adjusted to " << tail_point_marker[initial_side][marker] << endl;	
	report_markers(tail_point_marker,tail_point_coordinate);
}	
				}
				else if (initial_side%2 == 1 && mp_control.midpoints_not_tail_points)
				{
					disc_arc.push_back(0);
				}
			}
			else if (initial_side % 2 == 0 || (initial_side % 2 == 1 && mp_control.right_originating_tail_points))
			{
				marker = (initial_side%2==1? ODD_TAIL_LO:EVEN_TAIL_LO);
				coordinate = tail_point_coordinate[initial_side][tail_point_marker[initial_side][marker]];
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add lower originating tail point " << coordinate << " to disc_arc" << endl;
				disc_arc.push_back(coordinate);

				if (initial_side %2 == 1)
					tail_point_marker[initial_side][marker]++;
				else
					tail_point_marker[initial_side][marker]--;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost:   tail_point_marker for side " << initial_side << " adjusted to " << tail_point_marker[initial_side][marker] << endl;
	report_markers(tail_point_marker,tail_point_coordinate);
}    
			}
			else if (initial_side%2 == 1 && mp_control.midpoints_not_tail_points)
			{
				disc_arc.push_back(0);
			}
			
			int first_cudgel = (initial_side%2==1?(initial_side+1)/2: initial_side/2);
			int last_cudgel = (terminal_side%2==1?terminal_side/2: (terminal_side-1)/2);
			for (int c=first_cudgel; c <= last_cudgel; c++)
			{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add additional point " << additional_point_coordinate[c][additional_point_marker[c]] << " to disc_arc" << endl;
    
				disc_arc.push_back(additional_point_coordinate[c][additional_point_marker[c]]);
				additional_point_marker[c]--;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   additional_point_marker for cudgel " << c << " reduced to " << additional_point_marker[c] << endl;
			}
			
			if (upper_bone_label[terminal_label-1])
			{
				if (terminal_side % 2 == 1 || (terminal_side % 2 == 0 && mp_control.left_terminating_tail_points))
				{
					marker = (terminal_side%2==1? ODD_TAIL_UT:EVEN_TAIL_UT);
					coordinate = tail_point_coordinate[terminal_side][tail_point_marker[terminal_side][marker]];

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add upper terminating tail point " << coordinate << " to disc_arc" << endl;
    
					disc_arc.push_back(coordinate);
				
					if (terminal_side % 2 == 1)
						tail_point_marker[terminal_side][marker]--;
					else
						tail_point_marker[terminal_side][marker]++;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost:   tail_point_marker for side " << terminal_side << " adjusted to " << tail_point_marker[terminal_side][marker] << endl;
	report_markers(tail_point_marker,tail_point_coordinate);
}
				}
				else if (terminal_side % 2 == 0 && mp_control.midpoints_not_tail_points)
				{
					disc_arc.push_back(0);
				}
	
				marker = (terminal_side%2==1? ODD_HEAD_T:EVEN_HEAD_T);
				coordinate = head_point_coordinate[terminal_side][head_point_marker[terminal_side][marker]];

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add upper terminating head point " << coordinate << " to disc_arc" << endl;
    
				disc_arc.push_back(coordinate);
			
				if (terminal_side % 2 == 1)
					head_point_marker[terminal_side][marker]--;
				else
					head_point_marker[terminal_side][marker]++;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost:   head_point_marker for side " << terminal_side << " adjusted to " << head_point_marker[terminal_side][marker] << endl;	
	report_markers(head_point_marker,head_point_coordinate);
}    
			}
			else if (terminal_side % 2 == 1 || (terminal_side % 2 == 0 && mp_control.left_terminating_tail_points && terminal_label < side_start[terminal_side]+side_count[terminal_side]+1 ))  // i.e not a tail on an even side
			{
			
				marker = (terminal_side%2==1? ODD_TAIL_LT:EVEN_TAIL_LT);
				coordinate = tail_point_coordinate[terminal_side][tail_point_marker[terminal_side][marker]];
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add lower terminating tail point " << coordinate << " to disc_arc" << endl;
    
				disc_arc.push_back(coordinate);
				
				if (terminal_side % 2 == 1)
					tail_point_marker[terminal_side][marker]--;
				else
					tail_point_marker[terminal_side][marker]++;
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "write_lace_metapost:   tail_point_marker for side " << terminal_side << " adjusted to " << tail_point_marker[terminal_side][marker] << endl;				
	report_markers(tail_point_marker,tail_point_coordinate);
}    
			}							
			else if (mp_control.midpoints_not_tail_points)// tail on even side or non-tail lower terminating arc on even side when mp_control.left_terminating_tail_points is false
			{
				disc_arc.push_back(0);
			}			
		}
	
		disc_arc.push_back(label_z_coordinate[terminal_label-1]);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   add " << terminal_label << " to disc_arc" << endl;    
    
		disc_arc_sequence[i] = disc_arc;
	}


if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	for (int i=0; i< num_disc_arc_components; i++)
	{
		vector<int>& disc_arc = disc_arc_sequence[i];
		debug << "write_lace_metapost: disc_arc_sequence[" << i << "], size " << disc_arc.size() << ": ";
		for (unsigned int j=0; j< disc_arc.size(); j++)
			debug << disc_arc[j] << ' ';
		debug << endl;
	}
}

	/* now we can write the lace paths using path_body_arc */
	for (int i=0; i< num_cudgels; i++)
	{
		os << "p[" << i << "] = ";
		for (int j=1; j<= path_body_arc[i][0]; j++)
		{
			int label = path_body_arc[i][j];
			bool initial_label;
			int body_arc;
			
			int last_coordinate = -1; // to trap errors
			
			
			/* find label in arc_end_labels */
			for (int k=0; k< num_disc_arc_components; k++)
			{
				if (arc_end_labels[k][0] == label)
				{
					initial_label = true;
					body_arc = k;
					break;
				}
				else if (arc_end_labels[k][1] == label)
				{
					initial_label = false;
					body_arc = k;
					break;
				}
			}
			
			if (initial_label)
			{
				vector<int>& disc_arc = disc_arc_sequence[body_arc];
				
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "write_lace_metapost: process disc_arc_sequence[" << body_arc << "]" << endl;
		
				if (j != 0 && mp_control.lace_midpoints.size() > 0)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   check for midpoints between " << last_coordinate << " and " << disc_arc[0] << endl;    
					add_lace_midpoints(os,mp_control,last_coordinate,disc_arc[0]);
				}
						
				for (unsigned int k=0; k < disc_arc.size()-1; k++)
				{
					if (disc_arc[k] == 0)
					{
						os << "0.5[z" << disc_arc[k-1] << ",z" << disc_arc[k+1] << "]..";
					}
					else
					{
						os << 'z' << disc_arc[k];
						
						if (mp_control.tension)
							os << "..tension " << metapost_path_tension << "..";
						else
							os << "..";
						
						if (mp_control.adjacent_cudgel_midpoints && disc_arc.size() == 2 && !(upper_bone_coordinate[disc_arc[0]] && upper_bone_coordinate[disc_arc[1]]))
							os << "0.5[z" << disc_arc[0] << ",z" << disc_arc[1] << "]..";
							
						else if (mp_control.lace_midpoints.size() > 0)
						{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   check for midpoints between " << disc_arc[k] << " and " << disc_arc[k+1]<< endl;    
							add_lace_midpoints(os,mp_control,disc_arc[k],disc_arc[k+1]);
						}
					}
				}
				
				last_coordinate = disc_arc[disc_arc.size()-1];
				
				if (j == path_body_arc[i][0])
					os << 'z' << disc_arc[disc_arc.size()-1] << ";" << endl;
			}
			else
			{
				vector<int>& disc_arc = disc_arc_sequence[body_arc];
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "write_lace_metapost: process disc_arc_sequence[" << body_arc << "]" << endl;
						
				if (j != 0 && mp_control.lace_midpoints.size() > 0)
				{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   check for midpoints between " << last_coordinate << " and " << disc_arc.size()-1 << endl;    
					add_lace_midpoints(os,mp_control,last_coordinate,disc_arc.size()-1);
				}
				
				for (unsigned int k=disc_arc.size()-1; k>0 ; k--)
				{
					if (disc_arc[k] == 0)
					{
						os << "0.5[z" << disc_arc[k-1] << ",z" << disc_arc[k+1] << "]..";
					}
					else
					{
						os << 'z' << disc_arc[k];
						
					if (mp_control.tension)
						os << "..tension " << metapost_path_tension << "..";
					else
						os << "..";
						
					if (mp_control.adjacent_cudgel_midpoints && disc_arc.size() == 2 && !(upper_bone_coordinate[disc_arc[0]] && upper_bone_coordinate[disc_arc[1]]))
						os << "0.5[z" << disc_arc[0] << ",z" << disc_arc[1] << "]..";
					else if (mp_control.lace_midpoints.size() > 0)
					{
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "write_lace_metapost:   check for midpoints between " << disc_arc[k-1] << " and " << disc_arc[k]<< endl;    
						add_lace_midpoints(os,mp_control,disc_arc[k-1],disc_arc[k]);
					}
						
					}
				}
				
				last_coordinate = disc_arc[0];
				
				if (j == path_body_arc[i][0])
					os << 'z' << disc_arc[0] << ";" << endl;
			}						
		}
		
		os << "draw p[" << i << ']';
		if (mp_control.colour && i < static_cast<int>(mp_control.draw_colour.size()))		
			os << " withcolor " << mp_control.draw_colour[i];
		os << ';' << endl;
	}
	
	os << "endfig;" << endl;
}


void report_markers(matrix<int>& marker, matrix<int>& coordinate)
{
	int rows = marker.numrows();
	int cols = marker.numcols();
	
	for (int i=0; i< rows; i++)
	{
		if (cols == 3)
		{
			debug << "report_markers: head_point_markers for side " << i;
			if (i%2==1)
			{	
			    if (marker[i][ODD_HEAD_L] != -1)
					debug << " L = " << coordinate[i][marker[i][ODD_HEAD_L]];
					
				if (marker[i][ODD_HEAD_T] != -1)
					debug << " T = " << coordinate[i][marker[i][ODD_HEAD_T]];

				if (marker[i][ODD_HEAD_O] != -1)
					debug << " O = " << coordinate[i][marker[i][ODD_HEAD_L]];				
			}
			else
			{
			    if (marker[i][EVEN_HEAD_L] != -1)
					debug << " L = " << coordinate[i][marker[i][EVEN_HEAD_L]];

			    if (marker[i][EVEN_HEAD_T] != -1)
					debug << " T = " << coordinate[i][marker[i][EVEN_HEAD_T]];

			    if (marker[i][EVEN_HEAD_O] != -1)
					debug << " O = " << coordinate[i][marker[i][EVEN_HEAD_O]];
			}
		}
		else
		{
			debug << "report_markers: tail_point_markers for side " << i;
			if (i%2==1)
			{
				if (marker[i][ODD_TAIL_LT] != -1)
					debug << " LT = " << coordinate[i][marker[i][ODD_TAIL_LT]];

				if (marker[i][ODD_TAIL_LT] != -1)
					debug << " LO = " << coordinate[i][marker[i][ODD_TAIL_LO]];

				if (marker[i][ODD_TAIL_UT] != -1)					
					debug << " UT = " << coordinate[i][marker[i][ODD_TAIL_UT]];
					
				if (marker[i][ODD_TAIL_UO] != -1)
					debug << " UO = " << coordinate[i][marker[i][ODD_TAIL_UO]];
			}
			else
			{
				if (marker[i][EVEN_TAIL_LT] != -1)
					debug << " LT = " << coordinate[i][marker[i][EVEN_TAIL_LT]];
					
				if (marker[i][EVEN_TAIL_LO] != -1)
					debug << " LO = " << coordinate[i][marker[i][EVEN_TAIL_LO]];
					
				if (marker[i][EVEN_TAIL_UT] != -1)
					debug << " UT = " << coordinate[i][marker[i][EVEN_TAIL_UT]];
					
				if (marker[i][EVEN_TAIL_UO] != -1)
					debug << " UO = " << coordinate[i][marker[i][EVEN_TAIL_UO]];
			}
		}
		debug << endl;
	}
}

void set_out_head_coordinates(int side, int cudgel_x, float head_point_y_coord, int& z_index, matrix<int>& head_point_type_count, matrix<int>& head_point_coordinate, 
                              matrix<int>& tail_point_type_count, matrix<int>& tail_point_coordinate, ofstream& os)
{
//	int last_x_coordinate = 0;
	bool set_first_x_coordinate = false;
	int x_coordinate;
	int index = 0;	
	
	for (int i=0; i<3; i++)
	{
		int count;
		switch (i)
		{
			case 0: count = head_point_type_count[side][(side%2==1?ODD_HEAD_L:EVEN_HEAD_L)]; break;		
			case 1: count = head_point_type_count[side][(side%2==1?ODD_HEAD_T:EVEN_HEAD_O)]; break;
			case 2: count = head_point_type_count[side][(side%2==1?ODD_HEAD_O:EVEN_HEAD_T)]; break;
		}
		
		for (int j=0; j < count; j++)
		{

			bool relative_x_coordinate = true;
			bool relative_to_last = false;
			
			if (j==0 && !set_first_x_coordinate)
			{
				if (side%2 == 1 && tail_point_type_count[side][ODD_TAIL_LT] != 0)
				{
					x_coordinate = tail_point_coordinate[side][tail_point_type_count[side][ODD_TAIL_LT]-1];
				}
				else
				{
					x_coordinate = cudgel_x + (side%2==1?1:-1);
					relative_x_coordinate = false;
				}
					
				set_first_x_coordinate = true;
			}
			else
			{
				relative_to_last = true;
				x_coordinate = head_point_coordinate[side][index-1];
			}


/*
			if (j==0 && !set_first_x_coordinate)
			{
				x_coordinate = cudgel_x + (side%2==1?1:-1);
				
				set_first_x_coordinate = true;
			}
			else
			{
				x_coordinate = last_x_coordinate + (side%2==1?1:-1);
			}
*/


			os << "z" << z_index << "=(";
			if (relative_x_coordinate)
			{
				if (relative_to_last)
					os << "x" << x_coordinate << (side%2==1?"+h,":"-h,");
				else
					os << "x" << x_coordinate << ',';
			}
			else
				os << x_coordinate << "h,";
				
			os << head_point_y_coord << (head_point_y_coord != 0? "v);":");")  << endl;

				
//			os << "z" << z_index << "=(" << x_coordinate << "h," << head_point_y_coord << (head_point_y_coord != 0? "v);":");")  << endl;
				
    
			head_point_coordinate[side][index++] = z_index;
//			last_x_coordinate = x_coordinate;
			z_index++;
		}
	}
}

void set_out_tail_coordinates(int side, int cudgel_x, int max_num_lower_bones, int& z_index, matrix<int>& head_point_type_count, matrix<int>& head_point_coordinate, 
                              matrix<int>& tail_point_type_count, matrix<int>& tail_point_coordinate, ofstream& os)
{
	bool set_first_x_coordinate = false;
	int x_coordinate;
	int index = 0;	

	for (int i=0; i<4; i++)
	{
		int count;
		switch (i)
		{
			case 0: count = tail_point_type_count[side][(side%2==1?ODD_TAIL_LT:EVEN_TAIL_LO)]; break;
			case 1: count = tail_point_type_count[side][(side%2==1?ODD_TAIL_LO:EVEN_TAIL_LT)]; break;
			case 2: count = tail_point_type_count[side][(side%2==1?ODD_TAIL_UT:EVEN_TAIL_UO)]; break;
			case 3: count = tail_point_type_count[side][(side%2==1?ODD_TAIL_UO:EVEN_TAIL_UT)]; break;
		}
		
		for (int j=0; j < count; j++)
		{
			bool relative_x_coordinate = true;
			bool relative_to_last = false;
			
			if (j==0 && !set_first_x_coordinate)
			{
				if (side%2 == 0 && head_point_type_count[side][EVEN_HEAD_L] != 0)
				{
					x_coordinate = head_point_coordinate[side][head_point_type_count[side][EVEN_HEAD_L]-1];
				}
				else
				{
					x_coordinate = cudgel_x + (side%2==1?1:-1);
					relative_x_coordinate = false;
				}
					
				set_first_x_coordinate = true;
			}
			else
			{
				relative_to_last = true;
				x_coordinate = tail_point_coordinate[side][index-1];
			}
			
			os << "z" << z_index << "=(";
			if (relative_x_coordinate)
			{
				if (relative_to_last)
					os << "x" << x_coordinate << (side%2==1?"+h,":"-h,");
				else
					os << "x" << x_coordinate << ',';
			}
			else
				os << x_coordinate << "h,";
				
			os << -(max_num_lower_bones+1) << "v);" << endl;
			
			tail_point_coordinate[side][index++] = z_index;
			z_index++;
		}
	}			
}

void add_lace_midpoints(ofstream& os, metapost_control& mp_control, int a, int b)
{
if (debug_control::DEBUG >= debug_control::SUMMARY)
		debug << "add_lace_midpoints: given a = " << a << ", b = " << b << endl;
		
	pair<int,int> ab = {a, b};
	pair<int,int> ba = {b, a};

	vector<pair<int,int> >& edges = mp_control.lace_midpoints; 
	
	bool ab_found = find(edges.begin(),edges.end(),ab) != edges.end();
	bool ba_found = find(edges.begin(),edges.end(),ba) != edges.end();

	if (ab_found || ba_found)
	{
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	if (ab_found) 
		debug << "add_lace_midpoints: found (" << a << ',' << b << ") in midpoints" << endl;
	else if (ba_found)
		debug << "add_lace_midpoints: found (" << b << ',' << a << ") in midpoints" << endl;
}
		if (mp_control.midpoint_tension != 1.0)
			os << "tension " <<  mp_control.midpoint_tension << "..";
		else
			os << "0.5[z" << a << ",z" << b << "]..";
	}
}
