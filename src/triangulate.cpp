/*********************************************************************************************************************
Triangulate calculates the triangulation of the disc determined by an immersion code.  This triangulation is based on 
the one used by Knotscape: a set of triangulated sub-discs centred at the immersion crossings and barycentres of those 
components of the immersion's complement having more than two edes in their boundary.  The output of the triangulation 
is written in the same form as that produced by the knotscape function triang.

The triangulation is constructed from a set of distinct left and right turning cycles, that are in 1-1 correspondance 
with the components of the immersion's complement in S^2.  These turning cycles must be calculated by the calling code
and are passed as a parameter to the function defined here.  A turning cycle is selected to bound the infinte region of
the complement when it is considered to be in R^2; that is, the choice of this cycle fixes the point at infinity.  The 
infinite region must be selected by the calling code.  In R2, the "infinite" turning cycle contains in the compact 
component of its complement the remainder of the immersion.  

Note that the order in which the turning cycles are enumerated determines a nummbering for the regions of the immersion's 
complement.

There are five types of vertex in the triangulation:

1. Each crossing of the immersion corresponds to a vertex, crossing i is numbered as vertex i.
   
2. The midpoint of each edge of the immersion is a vertex.  The midpoint of edge i is vertex number (i + num_crossings)
   
3. The barycentry of each region in the immersion's complement having more than two edges in it's boundary is a vertex.  
   These regions are identified by turning cycles of length greater than two and their corresponding vertices are 
   numbered based on the ordering of the turning cycles.

4. If there are any Reidemeister I loops present in the immersion, then a pair of type 4 vertices are added either side 
   of the type 2 vertex corresponding to the loop edge.  Thus the Reidemeister I loop boundary comprises four vertices:
   {type 1, type 5, type 2, type 5}

5. There is a ring of vertices constructed in the infinte region in 1-1 correspondance with the edges in the turning cycle
   bounding the infinte region.  These vertices are placed radially outwards from the mid-point of the edges in this 
   turning cycle.  These vertices are numbered sequentially according to the vertex that corresponds to the mid-point of 
   each edge as we proceed around the turning cycle.
   
   If the infinite turning cycle includes an edge that forms a Reidemeister I loop then for each such loop we add two 
   additional type 5 vertices, radially associated with the two type 4 vertices on the loop.
   
   If the infinite turning cycle contains only two edges then we add two additional type 5 vertices, radially associated 
   with the two type 1 vertices in the infinite turning cycle, so that each type 5 vertex is connected to one type 1 
   vertex, one type 2 vertex and two other type 5 vertices.

If the immersion contains Reidemeister I loops, then the type 1 vertex corresponding to the loop crossings are connected by 
an edge to each of the type 4 vertices added to the loop and the type 2 midpoints of the other two edges incident with the 
crossing.  If, furthermore, the Reidemeister I loop lies in the infinite turning cycle, then the type 1 vertex is connected
to a pair of type 5 vertices, one radially associated with the ingress (non-loop) edge of the infinite turning cycle at that 
crossing with respect to the orientation of the infinite turning cycle and the other radially associated with the second type
4 vertex associated with the loop (second with respect to the orientation of the infinite turning cycle).

Otherwise, type 1 vertices are connected by an edge of the triangulation to the midpoint of the four edges in the 
immersion incident with the corresponding crossing.  If such a type 1 vertex lies in the infinite turning cycle, it is also 
connected to a type 5 vertex, as described below.

The compact regions of the immersion's complement are triangluated as follows.  For a region having two sides the midpoint 
vertices of the two edges are connected by an edge.  For a region having more than two sides each midpoint is connected by 
an edge to the midpoint of the adjacent edges in the turning cycle, or to a type 4 vertex if the adjacent edge is a 
Reidemeister I loop, and to the barycentre type 3 vertex.

Type 4 vertices are connected to the Reidemsieter I loop's type 1 vertex, its type 2 vertex, the other type 4 vertex belonging 
to the same loop

If the Reidemeister I loop does not form part of the infinite turning cycle, then it lies "inside" the immersion, in which 
case, the type 4 vertices are treated much like additional type 2 vertices.  Thus they are also connected to the barycentre 
of the adjacent non-loop region and to one of the type 2 vertices belonging to the edges either side of the Reicemeister I 
loop crossing.  If the Reidemeister I loop forms part of the infinite turning cycle then the type 4 vertices are also connected 
to two type 5 vertices, as described below.

If the infinite turning cycle contains at least three edges but no Reidemeister I loops, then the type 5 vertices are connected by 
edges of the triangulation to one type 1 vertex, two type 2 vertices and two other type 5 vertices, as follows:

 a) Each type 5 vertex is connected radially to the corresponding mid-point in an edge of the infinite region turning cycle. 
 b) A type 5 vertex is connected to the type 1 vertex reached by following the infinite turning cycle from the type 2 vertex in a)
 c) A type 5 vertex is connected to the type 2 vertex in the next edge of the infinite turning cycle following the type 1 vertex in b)
 d) The type 5 vertices are connected in a ring (the boundary of our triangulated disc)

If the infinite turning cycle contains at least three edges and a Reidemeister I loop then the type 5 vertices radially associated
with the first type 4 vertex of the loop and the type 2 vertex of the loop connect to only four other vertices: a type 2 vertex
a type 4 vertex and two other type 5 vertices.

Note that the choice of turning cycle to bound the infinite region may be a left or right turning cycle.  If it is a left turning 
cycle, as we follow the turning cycle we shall be moving in a clockwise direction, and if it is a right turning cycle an 
anti-clockwise direction.  This is the case regardless of which region is selected to contain the point at infinity.

July 2017: to aid force_directed placement, the orientation of the infinite region's turning cycle was added to the triangulation 
output file, using the format BOUNDARY: clockwise if it is a left turning cycle or BOUNDARY: anticlockwise otherwise. This enables 
the correct placement of the boundary vertices by force directed placement
*********************************************************************************************************************************/
using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <iomanip>

/********************* External variables ***********************/

extern ofstream     debug;
extern char const* triangulation_output_file;

extern bool USE_FORCE_DIRECTED_PLACEMENT;


#include <util.h>
#include <matrix.h>
#include <draw.h>


/********************* Function prototypes ***********************/
//pair<int,int> adjacent_edges(matrix<int>& cycle, int num_left_cycles, int region, int edge);
pair<int,int> adjacent_edges(matrix<int>& cycle, int num_left_cycles, int region, int edge, vector<int>&type_2_vertex, vector<int>&type_4_vertex, vector<int>& Reidemeister_I_loop_edges);
bool first_occurrence(generic_code_data& code_data, matrix<int>& infinite_cycle_first_visit, int crossing, int position);

void triangulate (generic_code_data& code_data, matrix<int>& cycle, int num_cycles, int num_left_cycles, int infinite_region)
{
	int num_crossings = code_data.num_crossings;
	int num_edges = 2*num_crossings;
	
	matrix<int>& code_table = code_data.code_table;	
	vector<int>& orig_crossing = code_data.orig_crossing;
	vector<int>& term_crossing = code_data.term_crossing;

	/* It is possible that the infinite region meets a crossing in two positions (note that only the
	   infinite region can do this since the components of the immersion's complement are connected).  
	   In this case we shall want to know for each crossing in the infinite region's turning cycle 
	   which edges (ingress and egress as we traverse the infinte regions bounding turning cycle)
	   meet these crossings first and which meet them second.  To determine this we record the 
	   ingress and egress edges of the cycle that first meet each crossing.  
	*/
	matrix<int> infinite_cycle_first_visit(num_crossings,2,-1);
	
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: calculating infinite_cycle_first_visit:" << endl;
    
	for (int i=1; i<= cycle[infinite_region][0]; i++)
	{
		int this_edge = cycle[infinite_region][i];
		int next_edge = (i<cycle[infinite_region][0]? cycle[infinite_region][i+1]: cycle[infinite_region][1]);

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   this edge  = " << this_edge << ", next_edge = " << next_edge << endl;

		int crossing = (this_edge < 0 ? orig_crossing[abs(this_edge)] : term_crossing[this_edge]);

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   intermediate crossing  = " << crossing << endl;
	
		if (infinite_cycle_first_visit[crossing][0] < 0)
		{
			infinite_cycle_first_visit[crossing][0] = abs(this_edge);
			infinite_cycle_first_visit[crossing][1] = abs(next_edge);

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   setting first visit edges" << endl;
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   first visit edges already set" << endl;
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: infinite_cycle_first_visit:" << endl;
    for (int i=0; i< num_crossings; i++)
    {
		debug << "triangulate:   ";
		debug << infinite_cycle_first_visit[i][0] << ' ' << infinite_cycle_first_visit[i][1] << endl;
	}
}

	/* count the number of Reidemeister I loops present in the diagram: each such loop corresponds to a
	   turning cycle of length 1
	*/
	int num_Reidemeister_I_loops=0;
	
	for (int i=0; i< num_cycles; i++)
	{
		if (cycle[i][0] == 1)
			num_Reidemeister_I_loops++;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "triangulate: num_Reidemeister_I_loops = " << num_Reidemeister_I_loops << endl;

	int index = 0;
	vector<int> Reidemeister_I_loop_edges(num_Reidemeister_I_loops);
	for (int i=0; i< num_cycles; i++)
	{
		if (cycle[i][0] == 1)
		{
			Reidemeister_I_loop_edges[index++] = abs(cycle[i][1]);
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: Reidemeister_1_loop_edges: ";
    for (int i=0; i < num_Reidemeister_I_loops; i++)
		debug << Reidemeister_I_loop_edges[i] << ' ';
	debug << endl;
}
	

	/* Determine the regions that meet at each crossing.  We shall number the
	   regions according to their position in cycle, since each turning cycle
	   corresponds to a region.
	   
	   The regions meeting at a crossing will be stored as a row of crossing_region
	   as follows:
	             
	     \ 3 /        
          \ / 
        0  X  2
          / \ 
         / 1 \
              
         
     where the crosing is drawn in the usual manner, so that the region between 
     the two ingress edges is stored in position 0, and so on.
     
	 We complete crossing_region by considering the left and right turning cycles in turn.
	 Note that we may tell whether an edge in a turning cycle is odd or even by its sign, 
	 since we always traverse even edges following the orientation and odd edges against
	 the orientation.  Also note that a crossing has an odd and even ingress (egress) edge,
	 so if for a left(right) turning cycle we know the type of the crossing and whether
	 we're arriving on an odd or even edge, we know which region around the crossing we're
	 following.
   */
	
	matrix<int> crossing_region(num_crossings,4);
	
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: calculate crossing regions" << endl;
    
	/* First the left turning cycles */
	for (int i=0; i<num_left_cycles; i++)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   (left) turning cycle " << i << endl;

		for (int j=1; j<=cycle[i][0]; j++)
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     edge " << cycle[i][j] << endl;

			int crossing;
			
			if (cycle[i][j] > 0)
				crossing = cycle[i][j]/2;
			else
				crossing = (abs(cycle[i][j])-1)/2; // we're going backwards along odd edges

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     takes us to crossing " << crossing;
    
			if (cycle[i][j] < 0)
			{
				if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
				{
					crossing_region[crossing][2] = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 2" << endl;
				}
				else
				{
					crossing_region[crossing][1] = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 1" << endl;
				}
			}
			else
			{
				if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
				{
					crossing_region[crossing][0] = i;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 0" << endl;
				}
				else
				{
					crossing_region[crossing][3] = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 3" << endl;
				}
			}
		}
	}

	/* Now the right turning cycles */
	for (int i=num_left_cycles; i<num_cycles; i++)
	{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   (right) turning cycle " << i << endl;

		for (int j=1; j<=cycle[i][0]; j++)
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     edge " << cycle[i][j] << endl;

			int crossing;
			
			if (cycle[i][j] > 0)
				crossing = cycle[i][j]/2;
			else
				crossing = (abs(cycle[i][j])-1)/2; // we're going backwards along odd edges

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     takes us to crossing " << crossing;
    
			if (cycle[i][j] < 0)
			{
				if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
				{
					crossing_region[crossing][3] = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 3" << endl;
				}
				else
				{
					crossing_region[crossing][2] = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 2" << endl;
				}
			}
			else
			{
				if (code_table[TYPE][crossing] == generic_code_data::TYPE1)
				{
					crossing_region[crossing][1] = i;
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 1" << endl;
				}
				else
				{
					crossing_region[crossing][0] = i;
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", position 0" << endl;
				}
			}
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: crossing_region: " << endl;
    for (int i=0; i< num_crossings; i++)
    {
		debug << "triangulate:  ";
		for (int j=0; j<4; j++)
			debug << crossing_region[i][j] << ' ';
		debug << endl;
	}
}

	/* note the type 2 vertices corresponding to the mid-point of edges,
	   see notes at the top of the file
	*/
	vector<int> type_2_vertex(num_edges);
	for (int i=0; i< num_edges; i++)
		type_2_vertex[i] = i+num_crossings;

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: type_2_vertex: ";
    for (int i=0; i< num_edges; i++)
		debug << type_2_vertex[i] << ' ';
	debug << endl;
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

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: type_3_vertex: ";
    for (int i=0; i< num_cycles; i++)
		debug << type_3_vertex[i] << ' ';
	debug << endl;
}

	/* note the type 4 vertices added to Reidemeister I loops */
	int num_type_4_vertices = 2*num_Reidemeister_I_loops; 
	vector<int> type_4_vertex(num_type_4_vertices);
	for (int i=0; i< num_Reidemeister_I_loops; i++)
	{
		type_4_vertex[2*i] = temp++;
		type_4_vertex[2*i+1] = temp++;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: num_type_4_vertices = " << num_type_4_vertices << endl;
    debug << "triangulate: type_4_vertex: ";
    for (int i=0; i< num_type_4_vertices; i++)
		debug << type_4_vertex[i] << ' ';
	debug << endl;
}


	int num_type_1234_vertices = temp;
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "triangulate: num_type_1234_vertices = " << num_type_1234_vertices << endl;

	/* finally note the type 5 vertices corresponding to the mid-point of the 
	   edges in the turning cycle bounding the infinite region and the inverse
	   mapping of type 5 vertex to type 2 vertex
	*/
	int num_type_5_vertices = cycle[infinite_region][0];
	
	int num_vertices = temp + num_type_5_vertices;	

	int num_Reidemeister_loop_edges_in_infinite_turning_cycle = 0;
	for (int i=1; i<= cycle[infinite_region][0]; i++)
	{		
		if (find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),abs(cycle[infinite_region][i])) != Reidemeister_I_loop_edges.end())
			num_Reidemeister_loop_edges_in_infinite_turning_cycle++;
	}
	
	if (num_Reidemeister_loop_edges_in_infinite_turning_cycle > 0)
	{
		num_type_5_vertices += 2*num_Reidemeister_loop_edges_in_infinite_turning_cycle;
		num_vertices += 2*num_Reidemeister_loop_edges_in_infinite_turning_cycle;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: adding " << 2*num_Reidemeister_loop_edges_in_infinite_turning_cycle 
          << " type_5_vertices to allow for Reidemeister I loops in the infinite turning cycle" << endl;	
}    
	}
	else if (num_type_5_vertices == 2)  // can only have two edges AND RI loops in the infinite turning cycle if the diagram is a figure 8
	{
		num_type_5_vertices += 2;
		num_vertices += 2;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "triangulate: infinite turning cycle has length two,  adding two type_5_vertices " << endl;
		
	}
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: num_vertices = " << num_vertices << endl;
    debug << "triangulate: num_type_5_vertices = " << num_type_5_vertices << endl;
}

	vector<int> type_5_vertex(num_type_5_vertices);

	for (int i=0; i< num_type_5_vertices; i++)
		type_5_vertex[i] = temp++;
		
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: type_5_vertex: ";
    for (int i=0; i< num_type_5_vertices; i++)
		debug << type_5_vertex[i] << ' ';
	debug << endl;
}

	/* create a matrix to record the neighbours of each vertex; that is, those other vertices to which it is 
	   connected by an edge of the triangulation.  
	   
	   Type 1 vertices are connected to at most five neighbours, 
	   type 2 to at most eight (if they sit between turning cycles of length at least 3, containing a barycentre
	   type 3 vertices are connected to at most the length of the longest turning cycle plus twice the 
	   number of Reidemeister I loops,
	   type 4 vertices always connect to five neighbours
	   (standard) type 5 vertices are connected to exactly five neighbours.  
	   
	   We therefore have the maximum number of neighbours being max(8,length of longest turning cycle + 2*num_Reidemeister_I_loops).  
	   We allow an additional column (0) to record the number of neighbours.
	   
	   July 2017: we include an additional two rows in the matrix to handle non-standard type 5
	   vertices (with four neighbours) that may be added to assit in the force directed placement
	   drawing of diagrams where we want to use a turning cycle bounding the infinite region that
	   has only two edges.
	*/
	int max_neighbours = 8;
	for ( int i=0; i< num_cycles; i++)
	{
		if (max_neighbours < cycle[i][0]+2*num_Reidemeister_I_loops)
			max_neighbours = cycle[i][0]+2*num_Reidemeister_I_loops;
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "triangulate: max_neighbours = " << max_neighbours << endl;
    
    matrix<int> neighbour(num_vertices+2,max_neighbours+1);
    
    /* Start by noting the neighbours of the type 5 vertices, since this enables us to note to which type 1 vertex each of 
       these vertices is joined to, which we shall need to complete the neighbours of the type 1 vertices.
      
       Work around the infinite turning cycle so we can deal with the exceptional cases of just two edges or Reidemeister I 
       loops in the infinite turning cycle, where we need to determine neighbours differently.
    */
    int vertex = num_type_1234_vertices;
    matrix<int> type_5_vertex_corresponding_to_type_1_vertex(num_crossings,2);
    
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: type 5 vertex neighbours:" << endl;
   
    for (int i=1; i<= cycle[infinite_region][0]; i++)
    {
		int edge = cycle[infinite_region][i];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   infinite turning cycle edge " << edge;
		
		if (cycle[infinite_region][0] == 2)
		{
			/* We have just two edges in the  infinite turning cycle, so we add a type 5 vertex radially at each type 1 vertex 
			   in the infinite turning cycle to give four in all (which allows circle packing to be performed).  The type 5
			   vertices are numbered consecutively as we proceed around the turning cycle, starting with the vertex radially 
			   associated with the first edge, then the terminating or originating crossing dependent on whether the edge is 
			   even or odd respectively, then similarly the second edge and remainign crossing.
			   
			   We have two type 5 vertices associated with each type 1 vertex in this case.  We record them so that the vertex
			   radially associated with the edge appears first and that radialy associated with the crossing appears second.
			   
			   To be consistent with the other cases, if the infinite turning cycle is a right cycle, we record the neighbours
			   in a clockwise rather than an anti-clockwise manner.
			*/
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "\ntriangulate:   infinite turning cycle contains just two edges" << endl;
    
			neighbour[vertex][0] = 4; 
			neighbour[vertex+1][0] = 4; 
			
			if (vertex==num_type_1234_vertices)
				neighbour[vertex][1] = num_vertices-1;
			else
				neighbour[vertex][1] = vertex-1;
			
			neighbour[vertex+1][1] = vertex;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 1 for type_5_vertex " << vertex << " is " << neighbour[vertex][1] 
          << ", neighbour 1 for type_5_vertex " << vertex+1 << " is " << neighbour[vertex+1][1] << endl;
}
	
			neighbour[vertex][2] = type_2_vertex[abs(edge)];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     neighbour " << 2 << " for type_5_vertex " << vertex << " is " << neighbour[vertex][2] << endl;

			if (edge < 0)
			{
				neighbour[vertex][3] = (abs(edge)-1)/2;
				neighbour[vertex+1][2] = (abs(edge)-1)/2;
				
				type_5_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][0] = vertex;
				type_5_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][1] = vertex+1;
				
			}
			else
			{
				neighbour[vertex][3] = edge/2;
				neighbour[vertex+1][2] = edge/2;
				
				type_5_vertex_corresponding_to_type_1_vertex[edge/2][0] = vertex;
				type_5_vertex_corresponding_to_type_1_vertex[edge/2][1] = vertex+1;
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 3 for type_5_vertex " << vertex << " is " << neighbour[vertex][3] 
          << ", neighbour " << 2 << " for type_5_vertex " << vertex+1 << " is " << neighbour[vertex+1][2] << endl;
}
		
			int next_edge = (i==1? cycle[infinite_region][2]: cycle[infinite_region][1]);
			
			neighbour[vertex+1][3] = type_2_vertex[abs(next_edge)];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     neighbour 3 for type_5_vertex " << vertex+1 << " is " << neighbour[vertex+1][3] << endl;
			
			if (vertex < num_vertices -2)
			{
				neighbour[vertex][4] = vertex+1;
				neighbour[vertex+1][4] = vertex+2;
			}
			else
			{
				neighbour[vertex][4] = vertex+1;
				neighbour[vertex+1][4] = num_type_1234_vertices;
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour " << 4 << " for type_5_vertex " << vertex << " is " << neighbour[vertex][4] 
          << ", neighbour " << 4 << " for type_5_vertex " << vertex+1 << " is " << neighbour[vertex+1][4] << endl;
}
			
			vertex += 2;
		}
		else if (find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),abs(edge)) != Reidemeister_I_loop_edges.end())
		{		
			/* we have to deal with the type 5 vertices radially associated with a type 4, type 2 and another type 4 vertex 
			
			   If this is the first edge in the infinite turning cycle, indicated by vertex == num_type_1234_vertices, then 
			   the type 5 vertex associated with the first type 4 vertex is the one numbered num_vertices-1, the type 5 corresponding
			   to the type 2 vertex is the one numbered num_type_1234_vertices and the one corresponding to the second type 4 vertex
			   is num_type_1234_vertices+1.  Otherwise the three type vertices are vertex, vertex+1 and vertex+2 respectively.
			*/
			int first_type_5_vertex;
			int second_type_5_vertex;
			int third_type_5_vertex;
			
            if (vertex == num_type_1234_vertices)
			{
			    first_type_5_vertex = num_vertices-1;
				second_type_5_vertex = vertex;
				third_type_5_vertex = vertex+1;
			}
			else
			{
			    first_type_5_vertex = vertex;
				second_type_5_vertex = vertex+1;
				third_type_5_vertex = vertex+2;
			}

			vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),abs(edge));
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();

			int first_type_4_vertex;
			int second_type_4_vertex;
			
			if (edge <0)
			{
				first_type_4_vertex = type_4_vertex[2*loop_edge_index +1];
				second_type_4_vertex = type_4_vertex[2*loop_edge_index];
			}
			else
			{
				first_type_4_vertex = type_4_vertex[2*loop_edge_index];
				second_type_4_vertex = type_4_vertex[2*loop_edge_index +1];
			}
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << ", associated with type 5 vertices " << first_type_5_vertex << ' ' << second_type_5_vertex << " and " << third_type_5_vertex << endl;
    debug << "triangulate:   Reidemeister I loop_edge_index = " << loop_edge_index << endl;
    debug << "triangulate:   first_type_4_vertex = " << first_type_4_vertex << " second_type_4_vertex = " << second_type_4_vertex << endl;
}    

			/* if vertex == num_type_1234_vertices, so this is the first edge in the infinite turning cycle,
			   we haven't visited the Reidemeister I crossing yet 
			*/
			if (edge < 0)
				type_5_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][(vertex == num_type_1234_vertices?0:1)] = third_type_5_vertex;
			else
				type_5_vertex_corresponding_to_type_1_vertex[edge/2][(vertex == num_type_1234_vertices?0:1)] = third_type_5_vertex;

			neighbour[first_type_5_vertex][0] = 4; 
			neighbour[second_type_5_vertex][0] = 4; 
			neighbour[third_type_5_vertex][0] = 5; 
			
			if (vertex==num_type_1234_vertices)
				neighbour[first_type_5_vertex][1] = num_vertices-(vertex == num_type_1234_vertices?2:1);
			else
				neighbour[first_type_5_vertex][1] = vertex-1;
			
			neighbour[second_type_5_vertex][1] = first_type_5_vertex;
			neighbour[third_type_5_vertex][1] = second_type_5_vertex;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 1 for type_5_vertex " << first_type_5_vertex << " is " << neighbour[first_type_5_vertex][1] 
          << ", neighbour 1 for type_5_vertex " << second_type_5_vertex << " is " << neighbour[second_type_5_vertex][1]
          << ", neighbour 1 for type_5_vertex " << third_type_5_vertex << " is " << neighbour[third_type_5_vertex][1] << endl;
}
	
			neighbour[first_type_5_vertex][2] = first_type_4_vertex;
			neighbour[second_type_5_vertex][2] = type_2_vertex[abs(edge)];
			neighbour[third_type_5_vertex][2] = second_type_4_vertex;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 2 for type_5_vertex " << first_type_5_vertex << " is " << neighbour[first_type_5_vertex][2] 
          << ", neighbour 2 for type_5_vertex " << second_type_5_vertex << " is " << neighbour[second_type_5_vertex][2]
          << ", neighbour 2 for type_5_vertex " << third_type_5_vertex << " is " << neighbour[third_type_5_vertex][2] << endl;
}

			neighbour[first_type_5_vertex][3] = type_2_vertex[abs(edge)];
			neighbour[second_type_5_vertex][3] = second_type_4_vertex;
			neighbour[third_type_5_vertex][3] = (edge < 0? (abs(edge)-1)/2: edge/2);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 3 for type_5_vertex " << first_type_5_vertex << " is " << neighbour[first_type_5_vertex][3] 
          << ", neighbour 3 for type_5_vertex " << second_type_5_vertex << " is " << neighbour[second_type_5_vertex][3]
          << ", neighbour 3 for type_5_vertex " << third_type_5_vertex << " is " << neighbour[third_type_5_vertex][3] << endl;
}
		
			neighbour[first_type_5_vertex][4] = second_type_5_vertex;
			neighbour[second_type_5_vertex][4] = third_type_5_vertex;

			if (vertex < num_vertices -1)
				edge = cycle[infinite_region][i+1];
			else
				edge = cycle[infinite_region][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     next infinite turning cycle edge = " << edge << endl;
    
			neighbour[third_type_5_vertex][4] = type_2_vertex[abs(edge)];
			
if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 4 for type_5_vertex " << first_type_5_vertex << " is " << neighbour[first_type_5_vertex][4] 
          << ", neighbour 4 for type_5_vertex " << second_type_5_vertex << " is " << neighbour[second_type_5_vertex][4]
          << ", neighbour 4 for type_5_vertex " << third_type_5_vertex << " is " << neighbour[third_type_5_vertex][4] << endl;
}

			if (vertex < num_vertices - 3)
				neighbour[third_type_5_vertex][5] = vertex+(vertex == num_type_1234_vertices?2:3);
			else
				neighbour[third_type_5_vertex][5] = num_type_1234_vertices;

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     neighbour 5 for type_5_vertex " << third_type_5_vertex << " is " << neighbour[third_type_5_vertex][5] << endl;
			
			vertex += (vertex == num_type_1234_vertices?2:3);
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", associated with standard type 5 vertex " << vertex << endl;

			neighbour[vertex][0] = 5; 
	
			if (vertex==num_type_1234_vertices)
				neighbour[vertex][1] = num_vertices-1;
			else
				neighbour[vertex][1] = vertex-1;
	
			neighbour[vertex][2] = type_2_vertex[abs(edge)];

			if (edge < 0)
			{
				neighbour[vertex][3] = (abs(edge)-1)/2;
				if (type_5_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][0])
					type_5_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][1] = vertex;
				else
					type_5_vertex_corresponding_to_type_1_vertex[(abs(edge)-1)/2][0] = vertex;
			}
			else
			{
				neighbour[vertex][3] = edge/2;
				if (type_5_vertex_corresponding_to_type_1_vertex[edge/2][0])
					type_5_vertex_corresponding_to_type_1_vertex[edge/2][1] = vertex;
				else
					type_5_vertex_corresponding_to_type_1_vertex[edge/2][0] = vertex;
			}
				
if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 1 = " << neighbour[vertex][1] 
          << " neighbour 2 = " << neighbour[vertex][2] << " neighbour 3 = " << neighbour[vertex][3] << endl;
}

			if (vertex < num_vertices -1)
				edge = cycle[infinite_region][i+1];
			else
				edge = cycle[infinite_region][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     next infinite turning cycle edge = " << edge << endl;

			vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),abs(edge));
			if (vptr != Reidemeister_I_loop_edges.end())
			{
				int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge, index " << loop_edge_index << endl;
				
				if (edge <0)
					neighbour[vertex][4] = type_4_vertex[2*loop_edge_index +1];
				else
					neighbour[vertex][4] = type_4_vertex[2*loop_edge_index];
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     type 2 vertex edge" << endl;
    
				neighbour[vertex][4] = type_2_vertex[abs(edge)];
			}
				
			if (vertex < num_vertices -1)
				neighbour[vertex][5] = vertex+1;
			else
				neighbour[vertex][5] = num_type_1234_vertices;

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     neighbour 4 = " << neighbour[vertex][4] 
          << " neighbour 5 = " << neighbour[vertex][5] << endl;
}		
			vertex++;
		}
	}
      
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: type_5_vertex_corresponding_to_type_1_vertex: " << endl;
    for (int i=0; i< num_crossings; i++)
    {
		debug << "triangulate:  ";
		for (int j=0; j<2; j++)
			debug << type_5_vertex_corresponding_to_type_1_vertex[i][j] << ' ';
		debug << endl;
	}
}
    

	/* determine the neighbours of the type 1 vertices; i.e the immersion crossings.
	   We start from the region in position 0 and work anti-clockwise around the crossing,
	   only the infinite region contributes a neighbour, and each edge contributes
	   its mid-point or a type 4 vertex, depending on whether it is a Reidemeister I loop edge
	*/
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: type 1 vertex neighbours" << endl;
	for (int i=0; i< num_crossings; i++)
	{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   vertex " << i << endl;

		int nbr = 1; // index into neighbour
		
		/* region in position 0 */
		int region = crossing_region[i][0];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     region " << region << endl;

		if (region == infinite_region)
		{
			if (cycle[infinite_region][0] == 2)
			{
				if (infinite_region < num_left_cycles)
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				}
				else
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
				}
			}
			else if (type_5_vertex_corresponding_to_type_1_vertex[i][1])
			{				
				if (first_occurrence(code_data, infinite_cycle_first_visit, i, 0))
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];								
			}
			else
			{
				neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		int edge;
		if (code_table[TYPE][i] == generic_code_data::TYPE1)
			edge = code_table[EVEN_TERMINATING][i];
		else
			edge = code_table[ODD_TERMINATING][i];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     edge" << edge << endl;

		vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge, index " << loop_edge_index << endl;
    
			neighbour[i][nbr++] = type_4_vertex[2*loop_edge_index+1];
    
		}
		else
		{
			neighbour[i][nbr++] = type_2_vertex[edge];
		}
		
		/* region in position 1 */
		
		region = crossing_region[i][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     region " << region << endl;

		if (region == infinite_region)
		{
			if (cycle[infinite_region][0] == 2)
			{
				if (infinite_region < num_left_cycles)
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				}
				else
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
				}
			}
			else if (type_5_vertex_corresponding_to_type_1_vertex[i][1])
			{
				
				if (first_occurrence(code_data, infinite_cycle_first_visit, i, 1))
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
				
			}
			else
			{
				neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		if (code_table[TYPE][i] == generic_code_data::TYPE1)
			edge = code_table[EVEN_ORIGINATING][i];
		else
			edge = code_table[ODD_ORIGINATING][i];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     edge " << edge << endl;

		vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge, index " << loop_edge_index << endl;
    
			neighbour[i][nbr++] = type_4_vertex[2*loop_edge_index];
    
		}
		else
		{
			neighbour[i][nbr++] = type_2_vertex[edge];
		}
		
		/* region in position 2 */
		region = crossing_region[i][2];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     region " << region << endl;

		if (region == infinite_region)
		{
			if (cycle[infinite_region][0] == 2)
			{
				if (infinite_region < num_left_cycles)
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				}
				else
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
				}
			}
			else if (type_5_vertex_corresponding_to_type_1_vertex[i][1])
			{
				
				if (first_occurrence(code_data, infinite_cycle_first_visit, i, 2))
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
				
			}
			else
			{
				neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		if (code_table[TYPE][i] == generic_code_data::TYPE1)
			edge = code_table[ODD_ORIGINATING][i];
		else
			edge = code_table[EVEN_ORIGINATING][i];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     edge " << edge << endl;

		vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge, index " << loop_edge_index << endl;
    
			neighbour[i][nbr++] = type_4_vertex[2*loop_edge_index];
    
		}
		else
		{
			neighbour[i][nbr++] = type_2_vertex[edge];
		}

		/* region in position 3 */
		region = crossing_region[i][3];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     region " << region << endl;

		if (region == infinite_region)
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:       type 5 vertices corresponding to crossing " << type_5_vertex_corresponding_to_type_1_vertex[i][0] << ' ' << type_5_vertex_corresponding_to_type_1_vertex[i][1] << endl;
			
			if (cycle[infinite_region][0] == 2)
			{
				if (infinite_region < num_left_cycles)
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				}
				else
				{
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
				}
			}
			else if (type_5_vertex_corresponding_to_type_1_vertex[i][1])
			{								
				if (first_occurrence(code_data, infinite_cycle_first_visit, i, 3))
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];				
				else
					neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][1];				
				
			}
			else
			{
				neighbour[i][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[i][0];
			}
		}
		
		/* following edge */
		if (code_table[TYPE][i] == generic_code_data::TYPE1)
			edge = code_table[ODD_TERMINATING][i];
		else
			edge = code_table[EVEN_TERMINATING][i];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     edge " << edge << endl;

		vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge, index " << loop_edge_index << endl;
    
			neighbour[i][nbr++] = type_4_vertex[2*loop_edge_index+1];    
		}
		else
		{
			neighbour[i][nbr++] = type_2_vertex[edge];
		}
		
		/* set the number of neighbours */
		neighbour[i][0] = nbr-1;
	}
	
	/* For type 2 vertices we start with the vertex reached by following the underlying orientation of the 
	   edge and then work anti-clockwise.  If the edge is a Reidemeister I loop edge this first vertex is 
	   a type 4 vertex, otherwise it is a type 1 vertex.
	   
	   Note that crossing_region tells us the two regions incident with the edge underlying the type 2 vertex.
	   
	   Note that if we have a Reidemeister I loop in the infinite turning cycle, since we follow even links
	   according to the underlying orientation and odd edges against the underlying orientation, we have the 
	   following order for the two type 5 vertices associated with loop edge type 2 vertex when proceeding 
	   anti-clockwise from start vertex as described above.  Note that these two type 5 vertices are numbered 
	   according to the orientation of the infinite turning cycle, we refer to them as 1 and 2 below
	   
	   Reidemeister I loop oriented anticlockwise
	    - even loop edge => right infinite turning cycle, infinite region right of loop edge, right of turning cycle,type 5 order 1, 2
	    - odd loop edge => left infinite turning cycle, infinite region right of loop edge, left of turning cycle, type 5 order 2, 1

	   Reidemeister I loop oriented clockwise
	    - even loop edge => left infinite turning cycle, infinite region left of loop edge, left of turning cycle, type 5 order 2, 1
	    - odd loop edge => right infinite turning cycle, infinite region left of loop edge, right of turning cycle, type 5 order 1, 2
	    
	*/
	
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: type 2 vertex neighbours" << endl;
	for (int edge=0; edge< num_edges; edge++)
	{
		int vertex = type_2_vertex[edge];
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: edge " << edge << ", vertex " << vertex << endl;

		int nbr = 1; // index into neighbour
			
		/* start with the next vertex in the direction of the terminating crossing */
		bool Reidemeister_I_loop_edge = false;
		int loop_edge_index = -1; // initialized to catch errors
		vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			Reidemeister_I_loop_edge = true;
			loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge, index " << loop_edge_index << endl;
    
			neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index+1];
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   terminating crossing " << term_crossing[edge] << endl;
    
			neighbour[vertex][nbr++] = term_crossing[edge];
		}		
    
		/* then the left hand region */
		int region;
		if (edge % 2)
		{			
			if (code_table[TYPE][term_crossing[edge]] == generic_code_data::TYPE1)
				region = crossing_region[term_crossing[edge]][3];
			else
				region = crossing_region[term_crossing[edge]][0];
		}
		else
		{			
			if (code_table[TYPE][term_crossing[edge]] == generic_code_data::TYPE1)
				region = crossing_region[term_crossing[edge]][0];
			else
				region = crossing_region[term_crossing[edge]][3];
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   left hand region " << region << endl;
    
		if (region == infinite_region)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     infinite region" << endl;

			if (Reidemeister_I_loop_edge)
			{
				/* There are two type 5 vertices in the infinite region.  If this loop is at the start of type 5 vertex numbering, then it must
				   be at the start of the infinite turning cycle, which means that the loop edge must be zero, by the way turning cycles are 
				   constructed.  IN this case, the first type 5 neighbour is num_vertices-1 and the second is num_type_1234_vertices.
				   Otherwise, the two type 5 neighbours start one after the type 5 vertex associated with the first occurrence of the
				   Reidemeister I loop crossing in the infinite turning cycle

				   Reidemeister I loop oriented clockwise
				    - even loop edge => left infinite turning cycle, infinite region left of loop edge, left of turning cycle, type 5 order 2, 1
				    - odd loop edge => right infinite turning cycle, infinite region left of loop edge, right of turning cycle, type 5 order 1, 2
								   
				*/
				
				/* identify the Reidemeister I loop crossings */
				int crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   Reidemeister I loop crossing =  " << crossing << " edge = " << edge << endl;
	
				if (edge == 0)
				{
					neighbour[vertex][nbr++] = num_vertices-1;								
					neighbour[vertex][nbr++] = num_type_1234_vertices;								
				}
				else if (edge%2)
				{
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+1;								
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+2;								
				}
				else
				{
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+2;								
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+1;								
				}	
			}
			else
			{				
				/* Identify the type 5 vertex associated with the edge, the offset into type_5_vertex is identified by locating the edge 
				   in the infinite turning cycle, taking into account the additional type 5 vertices added for Reidemeister I loops
				*/
				int offset = 0;
				for (int i=1; i<= cycle[infinite_region][0]; i++)
				{
					if (abs(cycle[infinite_region][i]) == edge)
						break;
					else if (find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),abs(cycle[infinite_region][i])) != Reidemeister_I_loop_edges.end())
						offset += (i==1?2:3);
					else
						offset += (cycle[infinite_region][0] == 2?2:1);
				}

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     offset in type_5_vertex of type 5 vertex corresponding to edge = " << offset << endl;
				
				/* We have two neighbours in the infinite region, the type 5 vertex associated with the  edge and the type 5 
				   vertex associated with a crossing.  If the edge is odd it is the type 5 vertex corresponding to the crossing 
				   we reach along this edge .  If the edge is even we have as a neighbour the type 5 vertex corresponding to 
				   the crossing we've just left.
				   
				   If there are two type 5 vertices connected to the crossing we have to select the appropriate occurrence, 
				   determined either by having just two edges in the infinite cycle, where we recorded the type 5 vertex
				   radialy associated with the crossing second, or by checking infinite_cycle_first_visit as described below.
				   
				   We are enumerating our neighbours in an anti-clockwise manner, so if the turning cycle bounding the infinite 
				   region is a left turning cycle we want the type 5 vertex corresponding to the edge first, and otherwise the 
				   type 5 vertex corresponding to the crossing first.
				*/
				int crossing = (edge % 2? term_crossing[edge]: orig_crossing[edge]);
				int occurrence;
			
				if (cycle[infinite_region][0] == 2)
				{
					occurrence = 1;
				}
				else
				{
					/* infinite_cycle_first_visit records the ingress and egress edges in the infinite turning cycle of
					   the first encounter at a particular crossing. Since we may be traversing a left or right turning 
					   cycle as we go around the infinite region, and we may be going forwards or backwards along an even 
					   or odd edge, then  we check both the ingress and egress turning cycle edges to see if edge is part 
					   of the first or second occurrence of the vertex in the infinite turning cycle.
					   
					*/
					if (infinite_cycle_first_visit[crossing][0] == edge || infinite_cycle_first_visit[crossing][1] == edge)
						occurrence = 0;
					else
						occurrence = 1;
				}
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     crossing = " << crossing << ", occurrence = " << occurrence << endl;
			

				if (infinite_region < num_left_cycles)
				{
					neighbour[vertex][nbr++] = type_5_vertex[offset];			
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
				}
				else
				{					
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
					neighbour[vertex][nbr++] = type_5_vertex[offset];			
				}
			}
		}
		else if (cycle[region][0] == 1)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     region is a Reidemeister I loop, no neighbours of the loop type 2 vertex " << endl;
		}
		else if (cycle[region][0] == 2)
		{
			/* identify the other edge bounding this region, since the neighbour
			   we want is the type 2 vertex on that edge
			*/
			int peer_edge;
			if (abs(cycle[region][1]) == edge)
				peer_edge = abs(cycle[region][2]);
			else
				peer_edge = abs(cycle[region][1]);

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     cycle_length 2, peer edge " << peer_edge << endl;

			neighbour[vertex][nbr++] = type_2_vertex[peer_edge];			
		}
		else
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     non-infinite region of cycle_length > 2" << endl;

			if (Reidemeister_I_loop_edge)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge neighbour is the barycentre type 3 vertex " << type_3_vertex[region] << endl;
    
				neighbour[vertex][nbr++] = type_3_vertex[region];
			}
			else
			{
				/* identify neighbouring type 2 vertices in the turning cycle */
				pair<int,int> type_24_neighbours = adjacent_edges(cycle, num_left_cycles, region, edge, type_2_vertex, type_4_vertex, Reidemeister_I_loop_edges);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     type 2 or 4 neighbours in turning cycle: " << 
          type_24_neighbours.first << ' ' << type_24_neighbours.second << endl;
}
			
				neighbour[vertex][nbr++] = type_24_neighbours.second;			

				neighbour[vertex][nbr++] = type_3_vertex[region];			

				neighbour[vertex][nbr++] = type_24_neighbours.first;			
			}
		}
	
		/* next the next vertex in the direction of the originating crossing */
		if (Reidemeister_I_loop_edge)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     first type 4 vertex on Reidemeister I loop edge is " << type_4_vertex[2*loop_edge_index] << endl;
    
			neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index];
		}
		else
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   originating crossing " << orig_crossing[edge] << endl;
    
			neighbour[vertex][nbr++] = orig_crossing[edge];
		}

		/* finally the right hand region */
		if (edge % 2)
		{
			if (code_table[TYPE][term_crossing[edge]] == generic_code_data::TYPE1)
				region = crossing_region[term_crossing[edge]][0];
			else
				region = crossing_region[term_crossing[edge]][1];
		}
		else
		{			
			if (code_table[TYPE][term_crossing[edge]] == generic_code_data::TYPE1)
				region = crossing_region[term_crossing[edge]][1];
			else
				region = crossing_region[term_crossing[edge]][0];
		}

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   right hand region " << region << endl;
    
		if (region == infinite_region)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     infinite region" << endl;

			if (Reidemeister_I_loop_edge)
			{
				/* There are two type 5 vertices in the infinite region.  If this loop is at the start of type 5 vertex numbering, then it must
				   be at the start of the infinite turning cycle, which means that the loop edge must be zero, by the way turning cycles are 
				   constructed.  In this case, the first type 5 neighbour is num_vertices-1 and the second is num_type_1234_vertices.
				   Otherwise, the two type 5 neighbours start one after the type 5 vertex associated with the first occurrence of the
				   Reidemeister I loop crossing in the infinite turning cycle
				   
				   Reidemeister I loop oriented anticlockwise
				    - even loop edge => right infinite turning cycle, infinite region right of loop edge, right of turning cycle,type 5 order 1, 2
				    - odd loop edge => left infinite turning cycle, infinite region right of loop edge, left of turning cycle, type 5 order 2, 1
				*/
				
				/* identify the Reidemeister I loop crossings */
				int crossing = term_crossing[edge];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   Reidemeister I loop crossing =  " << crossing << endl;

				if (edge == 0)
				{
					neighbour[vertex][nbr++] = num_vertices-1;								
					neighbour[vertex][nbr++] = num_type_1234_vertices;								
				}
				else if (edge%2)
				{
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+2;								
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+1;								
				}
				else
				{
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+1;								
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][0]+2;								
				}
				
			}
			else
			{    
				/* Identify the type 5 vertex associated with the edge, the offset into type_5_vertex is identified by locating the edge 
				   in the infinite turning cycle, taking into account the additional type 5 vertices added for Reidemeister I loops
				*/
				int offset = 0;
				for (int i=1; i<= cycle[infinite_region][0]; i++)
				{
					if (abs(cycle[infinite_region][i]) == edge)
						break;
					else if (find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),abs(cycle[infinite_region][i])) != Reidemeister_I_loop_edges.end())
						offset += (i==1?2:3);
					else
						offset += (cycle[infinite_region][0] == 2?2:1);
				}

				/* We have two neighbours in the infinite region, the type 5 vertex associated with the  edge and the type 5 
				   vertex associated with a crossing.  If the edge is odd it is the type 5 vertex corresponding to the crossing 
				   we reach along this edge .  If the edge is even we have as a neighbour the type 5 vertex corresponding to 
				   the crossing we've just left.
				   
				   If there are two type 5 vertices connected to the crossing we have to select the appropriate occurrence, 
				   determined either by having just two edges in the infinite cycle, where we recorded the type 5 vertex
				   radialy associated with the crossing second, or by checking infinite_cycle_first_visit as described below.
				   
				   We are enumerating our neighbours in an anti-clockwise manner, so if the turning cycle bounding the infinite 
				   region is a left turning cycle we want the type 5 vertex corresponding to the edge first, and otherwise the 
				   type 5 vertex corresponding to the crossing first.
				*/
				int crossing = (edge % 2? term_crossing[edge]: orig_crossing[edge]);
				int occurrence;
				
				if (cycle[infinite_region][0] == 2)
				{
					occurrence = 1;
				}
				else
				{
					if (infinite_cycle_first_visit[crossing][0] == edge || infinite_cycle_first_visit[crossing][1] == edge)
						occurrence = 0;
					else
						occurrence = 1;
				}
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     crossing = " << crossing << ", occurrence = " << occurrence << endl;
			
				if (infinite_region < num_left_cycles)
				{
					neighbour[vertex][nbr++] = type_5_vertex[offset];			
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
				}
				else
				{					
					neighbour[vertex][nbr++] = type_5_vertex_corresponding_to_type_1_vertex[crossing][occurrence];								
					neighbour[vertex][nbr++] = type_5_vertex[offset];			
				}
			}
		}
		else if (cycle[region][0] == 1)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     region is a Reidemeister I loop, no neighbours of the loop type 2 vertex " << endl;
		}
		else if (cycle[region][0] == 2)
		{
			/* identify the other edge bounding this region, since the neighbour
			   we want is the type 2 vertex on that edge
			*/
			int peer_edge;
			if (abs(cycle[region][1]) == edge)
				peer_edge = abs(cycle[region][2]);
			else
				peer_edge = abs(cycle[region][1]);

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     cycle_length 2, peer edge " << peer_edge << endl;

			neighbour[vertex][nbr++] = type_2_vertex[peer_edge];			
		}
		else
		{

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     non-infinite region of cycle_length > 2" << endl;

			if (Reidemeister_I_loop_edge)
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     Reidemeister I loop edge neighbour is the barycentre type 3 vertex " << type_3_vertex[region] << endl;
    
				neighbour[vertex][nbr++] = type_3_vertex[region];
			}
			else
			{

				/* identify neighbouring type 2 vertices in the turning cycle */
				pair<int,int> type_24_neighbours = adjacent_edges(cycle, num_left_cycles, region, edge, type_2_vertex, type_4_vertex, Reidemeister_I_loop_edges);

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     type 2 or 4 neighbours in turning cycle: " << 
          type_24_neighbours.first << ' ' << type_24_neighbours.second << endl;
}
				neighbour[vertex][nbr++] = type_24_neighbours.second;			
	
				neighbour[vertex][nbr++] = type_3_vertex[region];			
	
				neighbour[vertex][nbr++] = type_24_neighbours.first;			
			}
		}

		/* set the number of neighbours */
		neighbour[vertex][0] = nbr-1;
	}
	
	/* Now the type 3 vertices, whose neighbours we enumerate anti-clockwise around the vertex noting the mid-point 
	   vertex for each edge in the corresponding region's turning cycle.  Enumerating anti-clockwise means that for a 
	   right turning cycles we have to work backwards along the turning cycle
	   
	   Moving anti-clockwise means that the arrangement of any type 4 vertices we might encounter are as follows:
	   
	  Reidemeister I loop oriented anticlockwise
	    - even loop edge => right turning cycle, 2nd T4 T2 1st T4
	    - odd loop edge => left turning cycle, 2nd T4 T2 1st T4 
	
	   Reidemeister I loop oriented clockwise
	    - even loop edge => left turning cycle, 1st T4 T2 2nd T4 
	    - odd loop edge => right turning cycle, 1st T4 T2 2nd T4 
	*/

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: type 3 vertex neighbours" << endl;
	for (int i=0; i< num_left_cycles; i++)
	{
		int vertex = type_3_vertex[i];
		
		if (vertex)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: vertex " << vertex << endl;
    
			int nbr = 1; // index into neighbour

			for (int j=1; j<= cycle[i][0]; j++)
			{
				int edge = abs(cycle[i][j]);
				
				vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),edge);
				if (vptr != Reidemeister_I_loop_edges.end())
				{
					int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: edge " << edge << " is a Reidemeister I loop, type 4 vertices " << type_4_vertex[2*loop_edge_index] << ' ' << type_4_vertex[2*loop_edge_index+1] << endl;

					if (edge%2)
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index+1];
					else
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index];
						
					neighbour[vertex][nbr++] = type_2_vertex[edge];			
					
					if (edge%2)
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index];
					else
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index+1];
					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: edge " << edge << " type 2 vertex = " << type_2_vertex[edge] << endl;

					neighbour[vertex][nbr++] = type_2_vertex[edge];			
				}
			}			

			/* set the number of neighbours */
			neighbour[vertex][0] = nbr-1;
		}
	}

	for (int i=num_left_cycles; i< num_cycles; i++)
	{
		int vertex = type_3_vertex[i];
		
		if (vertex)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: vertex " << vertex << endl;
   
			int nbr = 1; // index into neighbour

			for (int j=cycle[i][0]; j>= 1; j--)
			{
				int edge = abs(cycle[i][j]);

				vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),edge);
				if (vptr != Reidemeister_I_loop_edges.end())
				{
					int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
					
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: edge " << edge << " is a Reidemeister I loop, type 4 vertices " << type_4_vertex[2*loop_edge_index] << ' ' << type_4_vertex[2*loop_edge_index+1] << endl;

					if (edge%2)
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index];
					else
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index+1];
						
					neighbour[vertex][nbr++] = type_2_vertex[edge];			
					
					if (edge%2)
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index+1];
					else
						neighbour[vertex][nbr++] = type_4_vertex[2*loop_edge_index];					
				}
				else
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate: edge " << edge << " vertex = " << vertex << " nbr = " << nbr << " type 2 vertex = " << type_2_vertex[edge] << endl;

					neighbour[vertex][nbr++] = type_2_vertex[edge];						
				}
			}			
			
			/* set the number of neighbours */
			neighbour[vertex][0] = nbr-1;

		}
	}
	
	/* Finally, the type 4 vertices, which all have exactly five neighbours.  We enumerate them anti-clockwise 
	   starting from either the Reidemeister I loop's type 1 vertex or its type 2 vertex so that the first three
	   neighbours are of type 1,4 and 2, or 2, 4, and 1, regardless of the orientation of the loop.  If the loop 
	   is part of the infinite turning cycle, then the remaining two neighbours are type 5 vertices.  Otherwise,
	   they are a type 2 and a type 3 vertex, taken in the appropriate order.
	*/
	
	for (int i=0; i< num_Reidemeister_I_loops; i++)
	{
		
if (debug_control::DEBUG >= debug_control::DETAIL && i==0)
    debug << "triangulate: type 4 vertex neighbours" << endl;
    
		int edge = Reidemeister_I_loop_edges[i];
		int first_T4 = type_4_vertex[2*i];
		int second_T4 = type_4_vertex[2*i+1];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   Reidemeister I loop edge " << edge << " first type 4 vertex =  " << first_T4 << " second type 4 vertex = " << second_T4 << endl;
		
		/* find the edge in the turning cycle that is not the loop itself */
		int region = -1;
		int offset;
		
		for (int j=0; j< num_cycles && region == -1; j++)
		{			
			if( cycle[j][0] == 1)
				continue;
				
			for (int k=1; k<= cycle[j][0]; k++)
			{
				if (abs(cycle[j][k]) == edge)
				{
					region = j;
					offset = k;
					break;
				}
			}
		}
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   Reidemeister I loop edge " << edge << " lies in the non-loop turning cycle " << region << " at offset " << offset << endl;
   
		if (region == infinite_region)
		{
			/* Each Reidemeister I loop in the infinite turning cycle has three type 5 vertices radially associated with the type 2
			   and the two type 4 vertices on the loop.  These type 5 vertices are numbered based on the orientation of the infinite
			   turning cycle and are referred to here as 1 2 and 3.  The type 5 vertex preceding these three, as determined by the
			   orientation of the infinite turning cycle, is referred to here as 0.
			   
			   If this loop is at the start of type 5 vertex numbering, then it must be at the start of the infinite turning cycle, 
			   which means that the loop edge must be zero, by the way turning cycles are constructed. In this case, the preceding
			   type 5 vertex (type_5_0) is the one associated with the second occurence of the loop crossing in the infinite turning
			   cycle, which is necessarily numbered num_vertices-2.  Otherwise, the type_5_0 vertex is the one associated with the 
			   first occurence of the loop crossing in the infinite turning cycle.				   
				   			   
			   The first and second type 4 vertices on the loop edge are numbered based on the orientation of the underlying immersion
			   and the loop may be oriented anti-clockwise or clockwise with respect to the infinte region.  Thus, if we consider the
			   neighbours anti-clockwise "within the loop" first, we have the following arrangement of neighbours:

			  Reidemeister I loop oriented anticlockwise
			    - 1st T4 vertex neighbours: T2 T4 T1 T5 T5    2nd T4 vertex neighbours: T1 T4 T2 T5 T5    
			    - even loop edge => right turning cycle, 1st T4 type 5 neighbours 0,1 2nd T4 type 5 neighbours 2,3
			    - odd loop edge => left turning cycle, 1st T4 type 5 neighbours 3,2 2nd T4 type 5 neighbours 1,0
			
			   Reidemeister I loop oriented clockwise
			    - 1st T4 vertex neighbours: T1 T4 T2 T5 T5    2nd T4 vertex neighbours: T2 T4 T1 T5 T5    
			    - even loop edge => left turning cycle, 1st T4 type 5 neighbours 1,0 2nd T4 type 5 neighbours 3,2
			    - odd loop edge => right turning cycle, 1st T4 type 5 neighbours 2,3 2nd T4 type 5 neighbours 0,1			   			   
			    
			   The three type 5 vertices are those following the first type 5 vertex corresponding to the Reidemeister I 
			   loop crossing.  Since these type 5 vertices were numbered together when we considered the loop edge in the
			   infinite turning cycle, we know that they do not wrap back to the first type 5 vertex.
			*/
			
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     non-loop region " << region << " is the infinite region" << endl;
    
			int crossing_type_5_vertex;
			int type_5_0;
			int type_5_1;
			int type_5_2;
			int type_5_3;

			if (edge == 0)
			{
				crossing_type_5_vertex = type_5_vertex_corresponding_to_type_1_vertex[term_crossing[edge]][1]; // calculated as a check
				type_5_0 = num_vertices-2;
				type_5_1 = num_vertices-1;
				type_5_2 = num_type_1234_vertices;
				type_5_3 = num_type_1234_vertices +1;
			}
			else
			{
				crossing_type_5_vertex = type_5_vertex_corresponding_to_type_1_vertex[term_crossing[edge]][0];
				type_5_0 = crossing_type_5_vertex;
				type_5_1 = crossing_type_5_vertex +1;
				type_5_2 = crossing_type_5_vertex +2;
				type_5_3 = crossing_type_5_vertex +3;
			}
			

if (debug_control::DEBUG >= debug_control::DETAIL)
{
    debug << "triangulate:     crossing_type_5_vertex = " << crossing_type_5_vertex << endl;
    debug << "triangulate:     type_5_0 = " << type_5_0 << " type_5_1 = " << type_5_1 << " type_5_2 = " << type_5_2 << " type_5_3 = " << type_5_3 << endl;
}

		    if (region < num_left_cycles)
		    {
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   left turning cycle " << endl;
		
				if (edge%2)
				{
					/* left turning cycle, odd edge => anti-clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = type_2_vertex[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = term_crossing[edge];
					neighbour[first_T4][4] = type_5_3;
					neighbour[first_T4][5] = type_5_2;

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = term_crossing[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = type_2_vertex[edge];
					neighbour[second_T4][4] = type_5_1;
					neighbour[second_T4][5] = type_5_0;
				}
				else
				{
					/* left turning cycle, even edge => clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = term_crossing[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = type_2_vertex[edge];
					neighbour[first_T4][4] = type_5_1;
					neighbour[first_T4][5] = type_5_0;

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = type_2_vertex[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = term_crossing[edge];
					neighbour[second_T4][4] = type_5_3;
					neighbour[second_T4][5] = type_5_2;
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   right turning cycle " << endl;
    
				if (edge%2)
				{
					/* right turning cycle, odd edge => clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = term_crossing[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = type_2_vertex[edge];
					neighbour[first_T4][4] = type_5_2;
					neighbour[first_T4][5] = type_5_3;

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = type_2_vertex[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = term_crossing[edge];
					neighbour[second_T4][4] = type_5_0;
					neighbour[second_T4][5] = type_5_1;
				}
				else
				{
					/* right turning cycle, even edge => anti-clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = type_2_vertex[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = term_crossing[edge];
					neighbour[first_T4][4] = type_5_0;
					neighbour[first_T4][5] = type_5_1;

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = term_crossing[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = type_2_vertex[edge];
					neighbour[second_T4][4] = type_5_2;
					neighbour[second_T4][5] = type_5_3;
				}
			}        
		}
		else
		{
			/* If the Reidemeister I loop does not lie in the infinite turning cycle then the order of the final two neighbours
			   depends only on the orientation of the Reidemeister I loop and not on whether the adjacent non-loop region 
			   corresponds to a left or right turning cycle.
			
			   Thus, if we consider the neighbours anti-clockwise "within the loop" first, we have the following arrangement of 
			   neighbours:

			   Reidemeister I loop oriented anticlockwise
			    - 1st T4 vertex neighbours: T2 T4 T1 T2 T3    2nd T4 vertex neighbours: T1 T4 T2 T3 T2    
			
			   Reidemeister I loop oriented clockwise
			    - 1st T4 vertex neighbours: T1 T4 T2 T3 T2    2nd T4 vertex neighbours: T2 T4 T1 T2 T3    
			*/
			
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:     non-loop region " << region << " is not the infinite region" << endl;
    
			int preceding_edge;
			int succeeding_edge;

			int cycle_length = cycle[region][0];
						
			if (offset > 1)
				preceding_edge = abs(cycle[region][offset-1]);
			else
				preceding_edge = abs(cycle[region][cycle_length]);

			if (offset < cycle_length)
				succeeding_edge = abs(cycle[region][offset+1]);
			else
				succeeding_edge = abs(cycle[region][1]);		
				
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << ", preceding_edge = " << preceding_edge << " succeeding_edge = " << succeeding_edge << endl;

			

		    if (region < num_left_cycles)
		    {
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   left turning cycle " << endl;
						
				if (edge%2)
				{
					/* left turning cycle, odd edge => anti-clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = type_2_vertex[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = term_crossing[edge];
					neighbour[first_T4][4] = type_2_vertex[succeeding_edge];
					neighbour[first_T4][5] = type_3_vertex[region];

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = term_crossing[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = type_2_vertex[edge];
					neighbour[second_T4][4] = type_3_vertex[region];
					neighbour[second_T4][5] = type_2_vertex[preceding_edge];
				}
				else
				{
					/* left turning cycle, even edge => clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = term_crossing[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = type_2_vertex[edge];
					neighbour[first_T4][4] = type_3_vertex[region];
					neighbour[first_T4][5] = type_2_vertex[preceding_edge];

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = type_2_vertex[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = term_crossing[edge];
					neighbour[second_T4][4] = type_2_vertex[succeeding_edge];
					neighbour[second_T4][5] = type_3_vertex[region];
				}
			}
			else
			{
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "triangulate:   right turning cycle " << endl;
    
				if (edge%2)
				{
					/* right turning cycle, odd edge => clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = term_crossing[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = type_2_vertex[edge];
					neighbour[first_T4][4] = type_3_vertex[region];
					neighbour[first_T4][5] = type_2_vertex[succeeding_edge];

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = type_2_vertex[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = term_crossing[edge];
					neighbour[second_T4][4] = type_2_vertex[preceding_edge];
					neighbour[second_T4][5] = type_3_vertex[region];
				}
				else
				{
					/* right turning cycle, even edge => anti-clockwise loop */
					neighbour[first_T4][0] = 5;
					neighbour[first_T4][1] = type_2_vertex[edge];
					neighbour[first_T4][2] = type_4_vertex[2*i+1];
					neighbour[first_T4][3] = term_crossing[edge];
					neighbour[first_T4][4] = type_2_vertex[preceding_edge];
					neighbour[first_T4][5] = type_3_vertex[region];

					neighbour[second_T4][0] = 5;
					neighbour[second_T4][1] = term_crossing[edge];
					neighbour[second_T4][2] = type_4_vertex[2*i];
					neighbour[second_T4][3] = type_2_vertex[edge];
					neighbour[second_T4][4] = type_3_vertex[region];
					neighbour[second_T4][5] = type_2_vertex[succeeding_edge];
				}
			}    
		}
	}

	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "triangulate: neighbour: " << endl;
    for (int i=0; i< num_vertices; i++)
    {
		debug << "triangulate:   vertex " << i << " num neighbours = " << neighbour[i][0] << ": ";
		for (int j=1; j<= neighbour[i][0]; j++)
			debug << neighbour[i][j] << ' ';
		debug << endl;
	}
}

	/* Write out the neighbours to outputfile.  The output file is required to number
	   the vertices from 1 not from zero as we have done internally.  
	   
	   The output file should record the neighbours in an anti-clockwise manner around the 
	   vertex.   We have done this in our calculation of neighbour except in the case that the 
	   turning cycle bounding the infinite region is a right turning cycle.  In this case
	   we have recorded the type 5 vertex neighbours in a clockwise manner.
	   
	   Finally note that the output file requires that for type 1,2, and 3 vertices the 
	   first neighbour appears at the end of the list of neighbours as well as at the beginning.  
	   This is to allow the circle packing code to calculate angle sums more easily.
	   
	   The format of the output is based on that used by knotscape triang.c.  
	   
	   In March 2016 the output format was modified by adding the type 5 vertices themselves to the 
	   end of their list of neighbours.  This is to simplify the calculation of inner hulls, since 
	   type 5 vertices lie in their own link. 

	   In July 2017, to aid force_directed placement, the orientation of the infinite region's turning 
	   cycle was added to the triangulation output file, using the format BOUNDARY: clockwise if it is 
	   a left turning cycle or BOUNDARY: anticlockwise otherwise. 
	   
	   In December 2021, to facilitate drawing diagrams that include Reidemeister I loops, the use of 
	   the variable GAMMA was modified.  In knotscape triang.c, GAMMA is set to vertex[a[2]]; that is, 
	   the terminating vertex (numbered from 1) of the edge paired in the Dowker code with edge 2 
	   (edges numbered from 1), or in our terminology the terminating vertex of edge 1.  This value was 
	   never used by my code, so GAMMA was repurposed to identify the number of Reidemeister I loops
	   contained in the diagram.
	*/
	ofstream output;
	output.open(triangulation_output_file);
	if (!output)
	{
		cout << "\nError opening output file " << triangulation_output_file << endl;
		exit(0);
    }
    else
    {
	
		output << "NODECOUNT: " << num_vertices << endl;
		output << "ALPHA/BETA/GAMMA: 1 " << num_type_1234_vertices+1 << ' ' << num_Reidemeister_I_loops << endl;
//        output << "GEOMETRY: hyperbolic" << endl;
        output << "LOOPS: ";
        for (int i=0; i< num_Reidemeister_I_loops; i++)
			output << Reidemeister_I_loop_edges[i] << ' ';
		output << endl;
        output << "BOUNDARY: " << (infinite_region < num_left_cycles ? "clockwise" : "anticlockwise") << endl;
        output << "FLOWERS:\n" << endl;
        
//		for (int i=0; i< (infinite_region < num_left_cycles? num_vertices : num_type_1234_vertices); ++i)
		for (int i=0; i< num_type_1234_vertices; ++i)
		{
			output << i+1 << " " << neighbour[i][0] << ' ';
            for (int j=1; j<=neighbour[i][0]; j++) 
				output << neighbour[i][j]+1 << ' ';
			output << neighbour[i][1]+1;
			output << endl; 
		}
		
		if (infinite_region < num_left_cycles)
		{
			for (int i=num_type_1234_vertices; i < num_vertices; ++i)
			{
				output << i+1 << " " << neighbour[i][0] << ' ';
				for (int j=1; j<=neighbour[i][0]; j++) 
					output << neighbour[i][j]+1 << ' ';
				output << i+1 << endl;
			}
		}
		else
		{
			for (int i=num_type_1234_vertices; i < num_vertices; ++i)
			{
				output << i+1 << " " << neighbour[i][0] << ' ';
				for (int j=neighbour[i][0]; j>=1; j--) 
					output << neighbour[i][j]+1 << ' ';

				output << i+1 << endl;

			}
		}
		
		output << "END" << endl;
		output.close();  	
	}
}

/* adjacent_edges identifies the pair of edges adjacent to edge in the turning cycle corresponding to region and then 
   determines the type 2 or type 4 vertices adjacent to the type 2 vertex corresponding to edge.  The adjacent vertices 
   are listed so that the sequence first, edge, second appears anti-clockwise (with respect to the infinite region).
   
   Here region is not the infinite region and therefore any adjacent Reidemeister I loop edges are not incident with
   the infinite region.  Since we follow even links according to the underlying orientation and odd edges against the 
   underlying orientation, we have the arrangement of first and second type 4 vertex on the Reidemeister I loop edge
   
   Reidemeister I loop oriented anticlockwise
    - even loop edge => right turning cycle, preceding_edge == 1st type 4, succeeding_edge == 2nd type 4
    - odd loop edge => left turning cycle, preceding_edge == 1st type 4, succeeding_edge == 2nd type 4

   Reidemeister I loop oriented clockwise
    - even loop edge => left turning cycle, preceding_edge == 2nd type 4, succeeding_edge == 1st type 4
    - odd loop edge => right turning cycle, preceding_edge == 2nd type 4, succeeding_edge == 1st type 4
*/

pair<int,int> adjacent_edges(matrix<int>& cycle, int num_left_cycles, int region, int edge, vector<int>& type_2_vertex, 
                             vector<int>& type_4_vertex, vector<int>& Reidemeister_I_loop_edges)
{
	pair<int,int> neighbours;
	
	int offset=1;
	int cycle_length = cycle[region][0];

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: cycle_length = " << cycle_length << endl;

	for (int i=1; i<= cycle_length; i++)
	{
		if (abs(cycle[region][i]) == edge)
			break;
		else
			offset++;
	}

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: offset of edge in turning cycle = " << offset << endl;

    if (region < num_left_cycles)
    {
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: left turning cycle" << endl;

		int preceding_edge;
		int succeeding_edge;
		
		if (offset > 1)
			preceding_edge = abs(cycle[region][offset-1]);
		else
			preceding_edge = abs(cycle[region][cycle_length]);

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: preceding_edge = " << preceding_edge << endl;

		vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),preceding_edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			/* Reidemeister I loop oriented anticlockwise
			    - odd loop edge => left turning cycle, preceding_edge == 1st type 4, succeeding_edge == 2nd type 4
			   Reidemeister I loop oriented clockwise
			    - even loop edge => left turning cycle, preceding_edge == 2nd type 4, succeeding_edge == 1st type 4
			*/
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
			if (preceding_edge %2)
				neighbours.first = type_4_vertex[2*loop_edge_index];
			else
				neighbours.first = type_4_vertex[2*loop_edge_index+1];
		}
		else
		{
			neighbours.first = type_2_vertex[preceding_edge];
		}

		
		if (offset < cycle_length)
			succeeding_edge = abs(cycle[region][offset+1]);
		else
			succeeding_edge = abs(cycle[region][1]);		

if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: succeeding_edge = " << succeeding_edge << endl;

		vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),succeeding_edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
			if (succeeding_edge %2)
				neighbours.second = type_4_vertex[2*loop_edge_index+1];
			else
				neighbours.second = type_4_vertex[2*loop_edge_index];
		}
		else
		{
			neighbours.second = type_2_vertex[succeeding_edge];
		}
	}
    else
    {
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: right turning cycle" << endl;

		int preceding_edge;
		int succeeding_edge;

		if (offset > 1)
			succeeding_edge = abs(cycle[region][offset-1]);
		else
			succeeding_edge = abs(cycle[region][cycle_length]);
		
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: succeeding_edge = " << succeeding_edge << endl;

		vector<int>::iterator vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),succeeding_edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			/* Reidemeister I loop oriented anticlockwise
				- odd loop edge => right turning cycle, preceding_edge == 2nd type 4, succeeding_edge == 1st type 4
			   Reidemeister I loop oriented clockwise
				- even loop edge => right turning cycle, preceding_edge == 1st type 4, succeeding_edge == 2nd type 4
			*/
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
			if (succeeding_edge %2)
				neighbours.second = type_4_vertex[2*loop_edge_index];
			else
				neighbours.second = type_4_vertex[2*loop_edge_index+1];
		}
		else
		{
			neighbours.second = type_2_vertex[succeeding_edge];
		}
		
		if (offset < cycle_length)
			preceding_edge = abs(cycle[region][offset+1]);
		else
			preceding_edge = abs(cycle[region][1]);		   
    
if (debug_control::DEBUG >= debug_control::DETAIL)
    debug << "adjacent_edges: preceding_edge = " << preceding_edge << endl;

		vptr = find(Reidemeister_I_loop_edges.begin(),Reidemeister_I_loop_edges.end(),preceding_edge);
		if (vptr != Reidemeister_I_loop_edges.end())
		{
			int loop_edge_index = vptr - Reidemeister_I_loop_edges.begin();
			if (preceding_edge %2)
				neighbours.first = type_4_vertex[2*loop_edge_index+1];
			else
				neighbours.first = type_4_vertex[2*loop_edge_index];
		}
		else
		{
			neighbours.first = type_2_vertex[preceding_edge];
		}    
    
	}
	
	return neighbours;
}
