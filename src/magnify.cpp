/****************************************************************************************
                                      Magnify
                  
Magnify looks for the set of edges in the triangulation whose length is less than 
average_triangulation_length_factor * average_triangulation_edge_length.  It identifies
the vertices incident with any such edge and increases the radius of the circles
constructed at these vertices by magnification_factor, recomputing the centres of the
vertices so the islands of magnified vertices remain circle-packed, even though they will 
overlap adjacent circles.

Magnify assumes that only type 1, 2, and 3 vertices will ever need magnifying and includes 
a check that no type 4 vertices are involved, aborting its operation if an example is encountered 
where this is not the case.  We do condiser type 4 vertices for the purposes of evaluating the
average_triangulation_edge_length, if required to do so by the command line.

The magnified vertices and the edges joining them form a (possibly disconnected) sub-graph
of the triangulation.  We call each component of this subgraph an island of magnified vertices.
Note that there is no particular structure to the magnified sub-graph: some components may 
contain no uni-valent vertices, some may be trees, some may fall somewhere between these two
extremes.  We have to iterate through each component magnifying the radii so the circles
corresponding to the vertices in the magnified sub-graph remain properly packed.

We select a starting vertex for each component of the magnified sub-graph by choosing one 
incident with the maximum number of edges belonging to that component.  We prefer "internal"
vertices, that is, ones for whom all their incident edges belong to the magnified sub-graph.

By definition, a vertex v to be magnified will always have one other neighbour that is 
also magnified.  However, the set of magnified neighbours may not be contiguous around v.
We move anticlockwise around the starting vertex repositioning magnified vertices in the 
same way as the initial placement does but we only re-position vertices
relative to other magnified vertices, not to non-magnified vertices.  Therefore, if the
set of magnified neighbours is not contiguous around the starting vertex we treat each contiguous
set separately and make sure we start at the beginning of a contiguous set, with respect to our
anticlockwise motion.  We push the first neighbour in a contiguous set away from the starting vertex 
by the appropriate amount and then work around the vertex re-positioning the neighbours as we go.  

As we work through the component of the magnified subgraph we add magnified vertices to the back
of a list when their centre has been repositioned.  When we pop the next vertex from the front of 
the list we look for a magnified neighbour whose centre has already been recomputed.  
Since we may identify a neighbour in the middle of one of the vertex's sets of contiguous magnified 
neighbours we start by working clockwise from the first neighbour repositioning neighbours until
we reach the end of the current contiguous set. Then we work anti-clockwise from the first neighbour
and reposition all other magnified neighbours, across all remaining contiguous sets.

*****************************************************************************************/
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

#include <util.h>
#include <matrix.h>
#include <draw.h>


extern bool INCLUDE_BOUNDARY_VERTICES;
extern double average_triangulation_edge_length;
extern float average_triangulation_length_factor;
extern float magnification_factor;
extern double two_pi;

extern char const* circlepack_output_file;
extern char const* triangulation_output_file;

/********************* Function prototypes ***********************/
double alpha_value(double a_x, double a_y, double b_x, double b_y, double x, double y);
void push_away(int c1, int c2, matrix<double>& vcoords, vector<double>& radius);
void reposition_anticlockwise(int vertex, int first_neighbour, matrix<int>& flowers, matrix<double>& vcoords, vector<double>& radius, 
                              list<int>& vertex_list, matrix<bool>& short_edge, vector<bool>& recomputed_centre);
void reposition_clockwise(int vertex, int first_neighbour, matrix<int>& flowers, matrix<double>& vcoords, vector<double>& radius, 
                              list<int>& vertex_list, matrix<bool>& short_edge, vector<bool>& recomputed_centre);

void magnify(generic_code_data& code_data)
{

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "magnify: code_data: ";
	write_code_data(debug,code_data);
	debug << endl;
	print_code_data(debug, code_data, "magnify: ");
}

	/* Read triangulation data from triangulation_output_file */
	int nodecount;
	int alpha;
	int beta;
	int gamma;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "magnify: read triangulation data from " << triangulation_output_file  << endl;

	ifstream triangulation(triangulation_output_file);
	if(!triangulation)
	{
		cout << "\nError opening triangulation data\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: could not open " << triangulation_output_file << endl;
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
	debug << "magnify: nodecount = " << nodecount << endl;
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
	debug << "magnify: alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
	}
	
	int num_type123_vertices = beta-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "magnify: number of type 1,2 & 3 vertices = " << num_type123_vertices << endl;

	/* store the flowers for each node in a matrix with the first column set to the number of neighbours
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
	debug << "magnify: node " << node << ", petal count = " << count << ": " << first_petal;
				
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
if (debug_control::DEBUG >= debug_control::SUMMARY)
    debug << "magnify: Error opening output file " << circlepack_output_file << endl;
		exit(0);
    }
    
    int num_vertices;
    vertices >> num_vertices;
    

    matrix<double> vcoords(nodecount,2);
    vector<double> radius(nodecount);

	for (int i=0; i< nodecount; i++)
	{
		for (int j=0; j< 2; j++)
			vertices >> vcoords[i][j];
	}
	for (int i=0; i< nodecount; i++)
		vertices >> radius[i];
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "magnify: vcoords: " << endl;
	for (int i=0; i< nodecount; i++)
		debug << "magnify:   v" << i << " = ("<< vcoords[i][0] << ", " << vcoords[i][1] << ")" << endl;
    debug << "magnify: initial radii: " << endl;
    for (int i=0; i< nodecount; i++)
		debug << "magnify:   vertex " << i << ": " << radius[i] << endl;
}

	/* identify the number of vertices to be considered for the evaluation of the average
	   triangulation edge length.
	*/
	if (INCLUDE_BOUNDARY_VERTICES)
	{
		num_vertices = nodecount;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: including boundary vertices in triangulation edge length analysis, num_vertices = " << num_vertices << endl;
	}
	else
	{
		num_vertices = num_type123_vertices;

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: no boundary vertices in triangulation edge length analysis, num_vertices = " << num_vertices << endl;
	}
	   
	
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: analysing triangulation edge lengths:" << endl;

	matrix<double> displacement(nodecount,2);

	/* determine the average length of an edge in the triangulation	*/


	int num_triangulation_edges = 0;
	average_triangulation_edge_length = 0; // external variable available to write_metapost
	
	for (int u=0; u < num_vertices; u++)
	{
		/* there is always one extra vertex recorded in the rows of flower, if it's a type 1, 2, 
		   or 3 vertex, there is a repeated petal, if it is a type 4 vertex, the vertex itself
		   appears at the end of the row.
		*/
		for (int j=1; j <= flowers[u][0]; j++)
		{				
			int v = flowers[u][j] - 1; // flowers records vertices numbered from 1
			
			if (v >= num_type123_vertices && !INCLUDE_BOUNDARY_VERTICES)
				continue; // type 4 vertex

			double delta_x = vcoords[v][0] - vcoords[u][0];
			double delta_y = vcoords[v][1] - vcoords[u][1];

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "magnify: vertices " << u  << ", " << v << endl;
	debug << "magnify:   v = (" << vcoords[v][0] << "," << vcoords[v][1] << ")" << endl;
	debug << "magnify:   u = (" << vcoords[u][0] << "," << vcoords[u][1] << ")" << endl;
	debug << "magnify:   delta = (" << delta_x << "," << delta_y << ")" << endl;
}
			
			double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify:   mod_delta =  " << mod_delta << endl;
	
			average_triangulation_edge_length += mod_delta;
			num_triangulation_edges++;

		}
	}

	/* The above loop has counted every edge twice and accumulated
	   every edge length twice, so the factor of two cancels out 
	*/
	average_triangulation_edge_length /= num_triangulation_edges;
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: average_triangulation_edge_length =  " << average_triangulation_edge_length << endl;

	/* now identify those vertices joined by an edge less than average_triangulation_length_factor * average_triangulation_edge_length */
	vector<bool> magnify_vertex(num_vertices);
	matrix<bool> short_edge(num_vertices,num_vertices);
	
	for (int u=0; u < num_vertices; u++)
	{
		for (int j=1; j <= flowers[u][0]; j++)
		{				
			int v = flowers[u][j] - 1; // flowers records vertices numbered from 1
			
			if (v >= num_type123_vertices && !INCLUDE_BOUNDARY_VERTICES)
				continue; // type 4 vertex

			double delta_x = vcoords[v][0] - vcoords[u][0];
			double delta_y = vcoords[v][1] - vcoords[u][1];	
			double mod_delta = sqrt(delta_x*delta_x + delta_y*delta_y);
			
			if (mod_delta < average_triangulation_length_factor * average_triangulation_edge_length)
			{
				magnify_vertex[u] = magnify_vertex[v] = true;
				short_edge[u][v] = short_edge[v][u] = true;
			}
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "magnify: magnify_vertex:" << endl;
	for (int i=0; i< num_vertices; i++)
		debug << "magnify:   vertex " << i << ": " << (magnify_vertex[i]? "true" : "false") << endl;
}

	/* check that no type 4 vertex has been considered to need magnifying, and abort if such an example is encountered */
	for (int i= num_type123_vertices; i<num_vertices; i++)
	{
		if (magnify_vertex[i])
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: error? magnify has encountered a type 4 vertex that needs magnifying, aborting magnification" << endl;
			return;
		}
	}

	/* increase the radii of those vertices to be magnified */
	for (int i=0; i< num_type123_vertices; i++)
	{
		if (magnify_vertex[i])
		{
			radius[i] *= magnification_factor;
		}
	}

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
    debug << "magnify: magnified radii: " << endl;
    for (int i=0; i< nodecount; i++)
		debug << "magnify:   vertex " << i << ": " << radius[i] << endl;
}

	/* So, from here we assume that magnified vertices form islands amongst type 1, 2, and 3 vertices.
	
	   Repeatedly look for the magnified vertex with the greatest number of edges connected to other magnified vertices
	   and use that vertex as the start of the repositioning.  we look repeatedly as there may be more than one island of 
	   magnified vertices.
	*/
	
	vector<bool> recomputed_centre(num_vertices);
	bool found;
	do
	{

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: looking for new starting vertex" << endl;

		int max_edges = 0;
		int starting_vertex=-1;
		found = false;
		for (int i=0; i< num_vertices; i++)
		{
			if (magnify_vertex[i] && !recomputed_centre[i])
			{
				/* count the number of edges connecting i to another vertex in the same island of magnified vertices */
				int count=0;
				for (int j=1; j <= flowers[i][0]; j++)
				{
					if (magnify_vertex[flowers[i][j]-1] && short_edge[i][flowers[i][j]-1])
						count++;												
				}
				
				if (count > max_edges)
				{
					/* if all the edges of the current starting vertex conect to other magnified vertices, 
					   prefer the current starting edge unles the new candidate also has the same property.
					*/
					if (starting_vertex != -1 && max_edges == flowers[starting_vertex][0])
					{
						if (count == flowers[i][0])
						{
							starting_vertex = i;
							max_edges = count;
						}
					}
					else
					{
						starting_vertex = i;
						max_edges = count;
					}
				}
					
				found = true;
			}
		}		

		

		if (found)
		{			
			/* we have found a starting vertex for an island of magnified vertices and proceed to iterate 
			   through the component.  We start by identifying the start of a contiguous set of magnified 
			   neighbours.  We look around the neighbours of the starting vertex to see if there's a 
			   neighbour that is not magnified or in a different island and then proceed anticlockwise 
			   until we encounter one that is.
			*/
			
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: starting_vertex = " << starting_vertex << ", num edges to other magnified vertices in the same island = " << max_edges << endl;

			list<int> vertex_list;

			/* The limits in the for loop here make use of the assumption that type 4 vertices are never magnified */
			int unmagnified_neighbour = -1;
			int first_neighbour = -1;
			bool unmagnified_neighbour_found = false;

			for (int i=1; i <= flowers[starting_vertex][0]; i++)
			{
				if (magnify_vertex[flowers[starting_vertex][i]-1] == false || !short_edge[starting_vertex][flowers[starting_vertex][i]-1])
				{
					unmagnified_neighbour_found = true;
					unmagnified_neighbour = i;
					break;
				}
			}
							
			if (unmagnified_neighbour_found)
			{
				for (int i=unmagnified_neighbour+1; i <= flowers[starting_vertex][0]; i++)
				{
					if (magnify_vertex[flowers[starting_vertex][i]-1] && short_edge[starting_vertex][flowers[starting_vertex][i]-1])
					{
						first_neighbour = i;
						break;
					}				
				}

				if (first_neighbour == -1) // i.e not found yet
				{
					for (int i=1; i < unmagnified_neighbour; i++)
					{
						if (magnify_vertex[flowers[starting_vertex][i]-1] && short_edge[starting_vertex][flowers[starting_vertex][i]-1])
						{
							first_neighbour = i;
							break;
						}				
					}
				}
			}
			else
			{
				first_neighbour = 1;
			}

if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: first_neighbour is " << flowers[starting_vertex][first_neighbour]-1 << ", index = " << first_neighbour << endl;

			/* since we have chosen the first petal to be the start of a contiguous set of neighbours 
			   we need only work anticlockwise around the starting vertex, once the first neighbour
			   has been repositioned.  The starting vertex will not be moved, but we need to note that 
			   it's current location is correct under the magnification.
			*/
			push_away(starting_vertex, flowers[starting_vertex][first_neighbour]-1, vcoords, radius);
			vertex_list.push_back(flowers[starting_vertex][first_neighbour]-1);
			recomputed_centre[starting_vertex] = true;
			recomputed_centre[flowers[starting_vertex][first_neighbour]-1] = true;				

			reposition_anticlockwise(starting_vertex, first_neighbour, flowers, vcoords, radius, vertex_list, short_edge, recomputed_centre);


			while (vertex_list.size() != 0)
			{
				int vertex = *vertex_list.begin();
				vertex_list.pop_front();

if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "magnify: pop " << vertex << " from front of vertex_list" << endl;
	debug << "magnify: petals: ";
	for (int i=1; i<= flowers[vertex][0]; i++)
			debug << flowers[vertex][i] << ' ';
	debug << endl;
}	
			
				/* we have already recomputed the centre of at least one of the neighbours of
				   this vertex but we need to identify where it is in the corresponding row
				   of flowers.  We store the index of the first neighbour for which the centre
				   has been recomputed in first_neighbour.
				*/
				for (int i=1; i<=flowers[vertex][0]; i++)
				{
					if (recomputed_centre[flowers[vertex][i]-1] == true)
					{
						first_neighbour = i;
if (debug_control::DEBUG >= debug_control::SUMMARY)
{
	debug << "magnify: first neighbour for which we already have recomputed the centre is " << flowers[vertex][i]-1
	      << ", index " << first_neighbour << endl;
}	      	      
						break;
					}
				}
				
				reposition_clockwise(vertex, first_neighbour, flowers, vcoords, radius, vertex_list, short_edge, recomputed_centre);
				reposition_anticlockwise(vertex, first_neighbour, flowers, vcoords, radius, vertex_list, short_edge, recomputed_centre);
			}
		}
		else
		{
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: no new starting vertex found" << endl;
		}
	} while (found);
	
	/* write the new location of the vertices and the expanded radii to the output file */
	ofstream output(circlepack_output_file);
	
	if (!output)
	{
		cout << "\nError opening output file\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "magnify: could not open " << circlepack_output_file << endl;
		exit(0);
	}		

	output << nodecount << endl;
	output << fixed;
	
	for (int i=0; i< nodecount; i++)
	{
		output << setprecision(6) << setw(10) << vcoords[i][0] 		
		       << setprecision(6) << setw(10) << vcoords[i][1] << endl;		
	}

	for (int i=0; i< nodecount; i++)
		output << setprecision(6) << setw(10) << radius[i] << endl;		
	
	output.close();
	
}

/* push_away pushes the circle c2, away from circle c1 along the line through their centres so that
   the two circles, whose radii have already been magnified, become adjacent once again.
   The position of c1 will remain unchanged, its radius will just be increased.
   We simply push c2 away from the c1 to become c3.  We do this by translating c1 to the origin and 
   considering the image, c2', of c2 under this translation.  The coordinates of c2' alow us to 
   determine the angle \alpha between the positive x-axis and the line from the origin through c2', 
   which in turn allows us to compute the image, c3', of c3 under the above translation.
*/
void push_away(int c1, int c2, matrix<double>& vcoords, vector<double>& radius)
{			

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "push_away: pushing " << c2 << " away from " << c1 << endl;

	/* evaluate c2', the image of c2 after translating c1 to the origin */
	double c2_prime_x = vcoords[c2][0] - vcoords[c1][0];
	double c2_prime_y = vcoords[c2][1] - vcoords[c1][1];
	
	/* the angle beta rotates c2' to the x-axis, the quadrant containing
	   c2' determines whether this is added or subtracted from 0, 90, 180
	   or 270 degrees.
	*/
	double beta = atan(abs(c2_prime_y)/abs(c2_prime_x));
	double alpha = 0;
	if (c2_prime_x > 0 && c2_prime_y > 0)
		alpha = beta;
	else if (c2_prime_x < 0 && c2_prime_y > 0)
		alpha = two_pi/2 - beta;
	else if (c2_prime_x < 0 && c2_prime_y < 0)
		alpha = two_pi/2 + beta;
	else // (c2_prime_x > 0 && c2_prime_y < 0)
		alpha = -beta;
	
	/* now work out where c3' is with the new alpha */
	double a = radius[c1]+radius[c2];
	double c3_prime_x = a*cos(alpha);
	double c3_prime_y = a*sin(alpha);
	
	/* finally, set the vertex coordinate of c2 to be the translate back (to c3) from c3' */
	vcoords[c2][0] = c3_prime_x + vcoords[c1][0];
	vcoords[c2][1] = c3_prime_y + vcoords[c1][1];
}



/* reposition_anticlockwise works around the vertex from the adjacent vertex indexed by first_neighbour recomputing 
   the centres of other magnified vertices not already recomputed.  The vertex indexed by first_neighbour has
   already had its centre recomputed and lies within a contiguous set of magnified neighbours.  As we work anticlockwise
   around vertex, if we come to the end of a contiguous set of magnified neighbours we step over the unmagnified neighbours
   and start the next set by pushung the first circle away from vertex.  
    
   If there are n neighbours of the vertex, we need to cycle the columns of flowers indexed 1 to n starting from 
   first_neighbour+1, which we do by cycling the indices modulo n allowing for the fact that column 0 stores the value n.
   Note that we do not want to return to first_neighbour, so we cycle through vertices 1,...,n-1 mod n.
*/
void reposition_anticlockwise(int vertex, int first_neighbour, matrix<int>& flowers, matrix<double>& vcoords, vector<double>& radius, 
                              list<int>& vertex_list, matrix<bool>& short_edge, vector<bool>& recomputed_centre)
{	
	int n = flowers[vertex][0];
	for (int i=1; i<n; i++)
	{
		int c2_index = ((first_neighbour-1)+i-1)%n +1;
		int c3_index = ((first_neighbour-1)+i)%n +1;
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      c2_index = " << c2_index << ", c3_index = " << c3_index << endl;
		
		/* we use c1, c2, and c3 to identify the vertices with the corresponding centres */
		int c1 = vertex;
		int c2 = flowers[vertex][c2_index]-1;
		int c3 = flowers[vertex][c3_index]-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      next neighbour c3 = " << c3 << ", c2 = " << c2 << ", c1 = " << c1 << endl;

		if (!short_edge[c1][c3])
		{
			/* step over any unmagnified vertices and push the start of the next contiguous magnified set
			   away from c1, if there's another one to consider.
			*/
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      next neighbour is unmagnified, moving to next contiguous set of magnified neighbours" << endl;

			while (++i < n)
			{
				c3_index = ((first_neighbour-1)+i)%n +1;
				c3 = flowers[vertex][c3_index]-1;
				
				if (short_edge[c1][c3] && !recomputed_centre[c3])
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      another contiguous set of magnified neighbours found" << endl;
					push_away(vertex, c3, vcoords, radius);
					break;
				}
			}
		}
		else
		{			
			if (recomputed_centre[c3] == false)
			{
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      computing centre for vertex c3 = " << c3 << ", c2 = " << c2 << ", c1 = " << c1 << endl;
				
				/* evaluate c2', the image of c2 after translating c1 to the origin */
				double c2_prime_x = vcoords[c2][0] - vcoords[c1][0];
				double c2_prime_y = vcoords[c2][1] - vcoords[c1][1];

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "reposition_anticlockwise:      c1 centre = (" << vcoords[c1][0] << "," << vcoords[c1][1] << ")" << endl;
	debug << "reposition_anticlockwise:      c2 centre = (" << vcoords[c2][0] << "," << vcoords[c2][1] << ")" << endl;
	debug << "reposition_anticlockwise:      c2_prime = (" << c2_prime_x << "," << c2_prime_x << ")" << endl;
}	
				
				/* the angle beta rotates c2' to the x-axis, the quadrant containing
				   c2' determines whether this is added or subtracted from 0, 90, 180
				   or 270 degrees.
				*/
				double beta = atan( abs(c2_prime_y)/abs(c2_prime_x));

				/* alpha is the angle c2c1c3 */
				double a = radius[c1]+radius[c3];
				double b = radius[c1]+radius[c2];
				double c = radius[c2]+radius[c3];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "reposition_anticlockwise:      a = " << a << ", b = " << b << ", c = " << c << endl;

				double cosine = (a*a+b*b-c*c)/2/a/b;

				if (cosine > 1)
				{
					cosine = 1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to 1 " << endl;	
				}
				else if (cosine < -1)
				{
					cosine = -1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to -1 " << endl;
				}

				double alpha = acos(cosine);
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "reposition_anticlockwise:      cosine = " << cosine << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      beta = " << beta << ", alpha = " << alpha << endl;

				/* adjust alpha by beta before calculating c3' (thus our c3' is really the image of c3'
				   as described above after the rotation by beta)
				   
				   We have to consider each quadrant due to the ambiguity of the acos function.
				*/
				if (c2_prime_x > 0 && c2_prime_y > 0)
					alpha += beta;
				else if (c2_prime_x < 0 && c2_prime_y > 0)
					alpha = alpha + two_pi/2 - beta; //alpha + \pi - beta
				else if (c2_prime_x < 0 && c2_prime_y < 0)
					alpha = alpha + two_pi/2 + beta; //alpha + \pi + beta
				else // (c2_prime_x > 0 && c2_prime_y < 0)
					alpha -= beta;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      alpha after beta adjustment = " << alpha << endl;

				
				/* now work out where c3' is with the new alpha */
				double c3_prime_x = a*cos(alpha);
				double c3_prime_y = a*sin(alpha);
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      c3_prime = (" << c3_prime_x << "," << c3_prime_x << ")" << endl;
				
				/* finally, translate back to c3 */
				vcoords[c3][0] = c3_prime_x + vcoords[c1][0];
				vcoords[c3][1] = c3_prime_y + vcoords[c1][1];
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      c3 centre = (" << vcoords[c3][0] << "," << vcoords[c3][1] << ")" << endl;
				
				/* push vertex c3 onto the vertex_list */
				vertex_list.push_back(flowers[vertex][c3_index]-1);
				recomputed_centre[flowers[vertex][c3_index]-1] = true;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:    push " << flowers[vertex][c3_index]-1 << " onto back of vertex_list" << endl;

			}
			else
			{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_anticlockwise:      already computed centre for vertex c3 = " << flowers[vertex][c3_index]-1 << endl;				
			}
		}
	}
}

/* reposition_clockwise works around the vertex from the adjacent vertex indexed by first_neighbour recomputing 
   the centres of other magnified vertices not already recomputed.  The vertex indexed by first_neighbour has
   already had its centre recomputed and lies within a contiguous set of magnified neighbours.  As we work clockwise
   around vertex, if we come to the end of a contiguous set of magnified neighbours we stop.  
    
   If there are n neighbours of the vertex, we need to cycle backwards through the columns of flowers indexed 1 to n 
   starting from first_neighbour-1, which we do by cycling the indices modulo n allowing for the fact that column 0 stores the value n.
   Note that we do not want to return to first_neighbour, so we cycle through vertices 1,...,n-1 mod n.
*/
void reposition_clockwise(int vertex, int first_neighbour, matrix<int>& flowers, matrix<double>& vcoords, vector<double>& radius, 
                              list<int>& vertex_list, matrix<bool>& short_edge, vector<bool>& recomputed_centre)
{	
	int n = flowers[vertex][0];
	for (int i=1; i<n; i++)
	{
		int c2_index = first_neighbour-i+1;
		if (c2_index <=0)
			c2_index += n;
			
		int c3_index = first_neighbour-i;
		if (c3_index <=0)
			c3_index += n;
		
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      c2_index = " << c2_index << ", c3_index = " << c3_index << endl;
		
		/* we use c1, c2, and c3 to identify the vertices with the corresponding centres */
		int c1 = vertex;
		int c2 = flowers[vertex][c2_index]-1;
		int c3 = flowers[vertex][c3_index]-1;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      next neighbour c3 = " << c3 << ", c2 = " << c2 << ", c1 = " << c1 << endl;

		if (!short_edge[c1][c3])
		{
			/* we've reached the end of the current contiguous magnified set, so we stop */

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      next neighbour is unmagnified, clockwise repositioning complete" << endl;
	
			return;
		}
		else
		{			
			if (recomputed_centre[c3] == false)
			{
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      computing centre for vertex c3 = " << c3 << ", c2 = " << c2 << ", c1 = " << c1 << endl;
				
				/* evaluate c2', the image of c2 after translating c1 to the origin */
				double c2_prime_x = vcoords[c2][0] - vcoords[c1][0];
				double c2_prime_y = vcoords[c2][1] - vcoords[c1][1];

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "reposition_clockwise:      c1 centre = (" << vcoords[c1][0] << "," << vcoords[c1][1] << ")" << endl;
	debug << "reposition_clockwise:      c2 centre = (" << vcoords[c2][0] << "," << vcoords[c2][1] << ")" << endl;
	debug << "reposition_clockwise:      c2_prime = (" << c2_prime_x << "," << c2_prime_x << ")" << endl;
}	
				
				/* the angle beta rotates c2' to the x-axis, the quadrant containing
				   c2' determines whether this is added or subtracted from 0, 90, 180
				   or 270 degrees.
				*/
				double beta = atan( abs(c2_prime_y)/abs(c2_prime_x));

				/* alpha is the angle c2c1c3 */
				double a = radius[c1]+radius[c3];
				double b = radius[c1]+radius[c2];
				double c = radius[c2]+radius[c3];

if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "reposition_clockwise:      a = " << a << ", b = " << b << ", c = " << c << endl;

				double cosine = (a*a+b*b-c*c)/2/a/b;

				if (cosine > 1)
				{
					cosine = 1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to 1 " << endl;	
				}
				else if (cosine < -1)
				{
					cosine = -1;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "calculate_angle_sum: cosine = " << cosine << ", adjusting to -1 " << endl;
				}

				double alpha = acos(cosine);
						
if (debug_control::DEBUG >= debug_control::INTERMEDIATE)
	debug << "reposition_clockwise:      cosine = " << cosine << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      beta = " << beta << ", alpha = " << alpha << endl;

				/* adjust alpha by beta before calculating c3' (thus our c3' is really the image of c3'
				   as described above after the rotation by beta)
				   
				   We have to consider each quadrant due to the ambiguity of the acos function.
				*/
				if (c2_prime_x > 0 && c2_prime_y > 0)
					alpha = beta - alpha;
				else if (c2_prime_x < 0 && c2_prime_y > 0)
					alpha = two_pi/2 - beta - alpha; // \pi - beta - alpha
				else if (c2_prime_x < 0 && c2_prime_y < 0)
					alpha = two_pi/2 + beta - alpha; // \pi + beta - alpha
				else // (c2_prime_x > 0 && c2_prime_y < 0)
					alpha = -beta - alpha;
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      alpha after beta adjustment = " << alpha << endl;

				
				/* now work out where c3' is with the new alpha */
				double c3_prime_x = a*cos(alpha);
				double c3_prime_y = a*sin(alpha);
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      c3_prime = (" << c3_prime_x << "," << c3_prime_x << ")" << endl;
				
				/* finally, translate back to c3 */
				vcoords[c3][0] = c3_prime_x + vcoords[c1][0];
				vcoords[c3][1] = c3_prime_y + vcoords[c1][1];
	
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      c3 centre = (" << vcoords[c3][0] << "," << vcoords[c3][1] << ")" << endl;
				
				/* push vertex c3 onto the vertex_list */
				vertex_list.push_back(flowers[vertex][c3_index]-1);
				recomputed_centre[flowers[vertex][c3_index]-1] = true;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:    push " << flowers[vertex][c3_index]-1 << " onto back of vertex_list" << endl;

			}
			else
			{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "reposition_clockwise:      already computed centre for vertex c3 = " << flowers[vertex][c3_index]-1 << endl;				
			}
		}
	}
}
