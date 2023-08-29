/****************************************************************************************
                                      Convex
                  
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
extern ofstream 	output;
extern ofstream     debug;

#include <util.h>
#include <matrix.h>
#include <draw.h>

extern int max_circle_packing_iterations;
extern int metapost_coordinate_scale_factor;
extern float magnification_factor;
extern double two_pi;

extern char const* circlepack_output_file;
extern char const* triangulation_output_file;

/********************* Function prototypes ***********************/
double alpha_value(double a_x, double a_y, double b_x, double b_y, double x, double y);
void write_convex_metapost(ofstream& os, vector<double> vertex_radius, vector<double> vertex_argument, string filename,
     matrix<int>& edge_matrix, metapost_control& mp_control);
double edge_argument(double x_delta, double y_delta);
double vector_modulus (double x, double y);
pair<double,double> midpoint_between_edges(int v1, int v2, int bl, int bm, int other_vertex, matrix<double>& vcoords);
pair<double,double> intersection_point (int v1, int v2, int other_vertex, matrix<double>& vcoords);

double ZERO_ERROR = 0.00000001;

void draw_convex_triangulation (metapost_control& mp_control, const char* filename)
{
	int num_vertices=0;
	int num_edges=0;
	int num_faces=0;
	matrix<int> edge_matrix(0,0);
	matrix<int> face_matrix(0,0);
		
	cout << "\ndrawing convex triangulation described by " << filename << endl;
	
	/* The Euler characteristic, v-e+f, of a sphere is 2, that of a disc is 1 */

	ifstream input(filename);
	if(!input)
	{
		cout << "\nError opening " << filename << "\n";
if (debug_control::DEBUG >= debug_control::SUMMARY)
	debug << "draw_convex_triangulation: could not open " << filename << endl;
		exit(0);
	}		
	
	string next_line;
	getline(input,next_line);
	if (next_line.find("NODECOUNT") != string::npos)
	{
		/* filename contains triangulation data from a peer code */

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: " << filename << " contains triangulation data from a peer code" << endl;

		char c;
		istringstream iss(next_line);
		do { iss >> c; } while (c != ':');
		iss >> num_vertices;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: num_vertices = " << num_vertices << endl;

		/* store the flowers for each node in a matrix with the first column set to the number of neighbours
		   of that vertex.  No vertex is adjacent to itself, nor is any vertex adjacent to all other vertices
		   so the maximum number of neighbours is <= nodecount -2.  Recall that for interior vertices the first
		   petal appears at the end of the list as well, to simplify the calculation of angle sums.
		*/
		matrix<int> flowers(num_vertices,num_vertices);
	
		while (getline(input,next_line))
		{
			if (next_line.find("FLOWERS") != string::npos)
			{
				getline(input,next_line); // read the empty line
				
				/* now read the flowers */
				while(getline(input,next_line))
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
	
					for (int i=0; i< count; i++)
					{
						iss >> petal;
						flowers[node][i+2] = petal;						
					}	
				}		
			}
		}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: flowers" << endl;					
	print(flowers,debug,4,"draw_convex_triangulation: ");
}

		edge_matrix = matrix<int>(num_vertices, num_vertices);

		for (int i=0; i< num_vertices; i++)
		{
			for (int j=1; j<= flowers[i][0] ; j++)	
				edge_matrix[i][flowers[i][j]-1] = 1;
		}

		for (int i=0; i< num_vertices; i++)
		for (int j=0; j<i; j++)
			num_edges += edge_matrix[i][j];
		
		num_faces = num_edges-num_vertices+1; // v-e+f=1
		face_matrix = matrix<int>(num_faces,3);
		
		/* evaluate the face_matrix from the flowers. */
		int row = 0;
		for (int i=0; i< num_vertices; i++)
		{
			for (int j=1; j<= flowers[i][0] ; j++)	
			{
				if (flowers[i][j] > i+1 && flowers[i][j+1] > i+1)
				{
					face_matrix[row][0] = i+1;
					face_matrix[row][1] = flowers[i][j];
					face_matrix[row][2] = flowers[i][j+1];
					row++;
				}
			}
		}
		
	}
	else
	{
		/* filename contains native triangulation data */

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: " << filename << " contains native triangulation data" << endl;

		/* read the triangles from the file into a list of faces to count them */
		list<vector<int> > face_list;
		do 
		{
			vector<int> face(3);
			istringstream iss(next_line);
			
			for (int i=0; i< 3; i++)	
				iss >> face[i];
				
			face_list.push_back(face);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: adding face to face_list: ";
	for (int i=0; i< 3; i++)
		debug << face[i] << ' ';
	debug << endl;
}
		} while (getline(input,next_line));
		
		num_faces = face_list.size();		

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: read " << num_faces << " from " << filename << endl;

		
		/* assign the face_matrix from the list */
		face_matrix = matrix<int>(num_faces,3);
		
		list<vector<int> >::iterator lptr = face_list.begin();
		num_faces=0;
		while (lptr != face_list.end())
		{
			vector<int>& face = *lptr;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: processing face from face_list: ";
	for (int i=0; i< 3; i++)
		debug << face[i] << ' ';
	debug << endl;
}

			for (int j=0; j< 3; j++)	
				face_matrix[num_faces][j] = face[j];			
			
			num_faces++;
			lptr++;
		}
		
		/* count the number of edges in the face_matrix into a local edge matrix, 
		   sized basedon the number of faces provided		   
		*/
		matrix<int> local_edge_matrix(3*num_faces, 3*num_faces);
		int max_face_vertex = 0;
		
		for (int i=0; i< num_faces; i++)
		{
			int v1 = face_matrix[i][0];
			int v2 = face_matrix[i][1];
			int v3 = face_matrix[i][2];
		
			if (v1 > max_face_vertex)
				max_face_vertex = v1;
			if (v2 > max_face_vertex)
				max_face_vertex = v2;
			if (v3 > max_face_vertex)
				max_face_vertex = v3;
			
			v1--;
			v2--;
			v3--;
			
			local_edge_matrix[v1][v2] = local_edge_matrix[v2][v1] = 1;
			local_edge_matrix[v1][v3] = local_edge_matrix[v3][v1] = 1;
			local_edge_matrix[v2][v3] = local_edge_matrix[v3][v2] = 1;
		}
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: local_edge_matrix:" << endl;					
	print(local_edge_matrix,debug,3,"draw_convex_triangulation: ");
}
		for (int i=0; i< 3*num_faces; i++)
		for (int j=0; j< 3*num_faces; j++)
			num_edges += local_edge_matrix[i][j];
			
		num_edges/=2;
		
		/* evaluate the number of vertices from the number of faces and edges provided
		   and check that correlates with the largest vertex specified in the input file
		*/
		num_vertices = num_edges-num_faces+1; // v-e+f=1

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: num_vertices = " << num_vertices << endl;					
	debug << "draw_convex_triangulation: num_edges = " << num_edges << endl;					
	debug << "draw_convex_triangulation: num_faces = " << num_faces << endl;
}
		
		if (num_vertices != max_face_vertex)
		{
			cout << "Error! Inconsistent faces specified in " << filename <<  ".\nMaximum vertex amongst triangles = " 
			     << max_face_vertex << ", number of vertices calculated from given faces = " << num_vertices << endl;

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: Error! Inconsistent faces specified in " << filename 
	      <<  ", max_face_vertex = " << max_face_vertex << ", num_vertices calculated from faces = " << num_vertices << endl;
}
			exit(0);
		}
		
		/* Assign edge_matrix from local_edge_matrix, omitting the unnecessary rows and columns */
		edge_matrix = matrix<int>(num_vertices,num_vertices);
		for (int i=0; i< num_vertices; i++)
		for (int j=0; j< num_vertices; j++)
			edge_matrix[i][j] = local_edge_matrix[i][j];
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: Euler characteristic = " << num_vertices-num_edges+num_faces << endl;					
	debug << "draw_convex_triangulation: edge_matrix" << endl;					
	print(edge_matrix,debug,3,"draw_convex_triangulation: ");
	debug << "draw_convex_triangulation: face_matrix" << endl;					
	print(face_matrix,debug,3,"draw_convex_triangulation: ");
}

	/* identify the initial boundary edges, which are those that lie in just one face */
	matrix<int> boundary_edges(num_edges,2);
	int num_boundary_edges=0;
	
	for (int i=0; i< num_vertices; i++)
	for (int j=0; j<i; j++)
	{
		if (edge_matrix[i][j] == 1)
		{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "draw_convex_triangulation: looking for edge " << i+1 << ' ' << j+1 << endl;					
			/* count how many times edge i+1, j+1 appears in face_matrix */
			int count = 0;
			
			for (int r=0; r< num_faces; r++)
			{
				bool i_found=false;
				bool j_found=false;

				for (int c=0; c<3; c++)
				{
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "draw_convex_triangulation: face_matrix " << r << ' ' << c << "= " << face_matrix[r][c] << endl;					
					
					if (face_matrix[r][c] == i+1)
						i_found = true;

					if (face_matrix[r][c] == j+1)
						j_found = true;
				}
				
if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "draw_convex_triangulation: i_found = " << i_found << " j_found = " << j_found << endl;					

				if (i_found && j_found)
				{
					if (++count == 2)
						break;
				}
			}

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "draw_convex_triangulation: count = " << count << endl;					
			
			if (count == 1)
			{
				boundary_edges[num_boundary_edges][0] = i+1;
				boundary_edges[num_boundary_edges][1] = j+1;
				num_boundary_edges++;
			}
		}
	}
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: num_boundary_edges = " << num_boundary_edges << endl;					
	debug << "draw_convex_triangulation: boundary_edges" << endl;					
	print(boundary_edges,debug,3,"draw_convex_triangulation: ");
}
	
	
	/* the number of initial booundary vertices is equal to the number of initial boundary edges */
	int num_interior_vertices = num_vertices - num_boundary_edges;

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: num_interior_vertices = " << num_interior_vertices << endl;					
	
	/* set epsilon, the delta on the vertex radius to be 1/(num_interior vertices+1) so they are
	   equally spaced throught the unit circle
	*/
	double epsilon = double(1)/(num_interior_vertices+1);


if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: epsilon = " << epsilon << endl;					
	
	
	/* Create a copy of edge_matrix to record those edges that are removed as the shelling proceeds.
	   This is used to identify interiod edges that have to be considered when placing internal vertices
	*/
	matrix<int> local_edge_matrix(edge_matrix);
	
	/* create a list to hold the boundary as we carry out the shelling, so we may insert additional 
	   vertices as a result of type I collapses
	*/
	vector<int> boundary;
	boundary.push_back(boundary_edges[0][0]);
	boundary.push_back(boundary_edges[0][1]);
	int last_boundary_vertex = boundary_edges[0][1];
	for (int i=0; i< num_boundary_edges-2; i++)
	{
		/* find last_boundary_edge in boundary_edge and push it's peer onto the list,
		   setting its peer to be the new last_boundary_vertex.  Set vertices to zero
		   to stop us finding them again.
		*/
		   
		for (int j=1; j< num_boundary_edges; j++)
		{
			if (boundary_edges[j][0] == last_boundary_vertex)
			{
				boundary.push_back(boundary_edges[j][1]);
				last_boundary_vertex = (boundary_edges[j][1]);
				boundary_edges[j][0]=boundary_edges[j][1]=0;
				break;
			}
			else if(boundary_edges[j][1] == last_boundary_vertex)
			{
				boundary.push_back(boundary_edges[j][0]);
				last_boundary_vertex = (boundary_edges[j][0]);
				boundary_edges[j][0]=boundary_edges[j][1]=0;
				break;
			}
		}
	}
	
	
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: initial boundary: ";					
	for (unsigned int b=0; b< boundary.size(); b++)
		debug << boundary[b] << ' ';
	debug << endl;
}	

	/* we record the radius and argument of each vertex, the boundary vertices are evenly spaced
	   around the unit sphere
	*/
	vector<double> vertex_radius(num_vertices);
	vector<double> vertex_argument(num_vertices);
	vector<int> interior_vertex(num_vertices,1); // records remaining interior vertices

	/* we record the vertex coordinates for those vertices that have been placed to assist with the
	   comparison with interior edges 
	*/
	matrix<double> vcoords(num_vertices,2);

	double argument = 0;
	for (unsigned int i=0; i< boundary.size(); i++)
	{		
		int v = boundary[i]-1;
		vertex_radius[v] = 1;
		vertex_argument[v] = argument;
		vcoords[v][0] = vertex_radius[v]*cos(vertex_argument[v]);
		vcoords[v][1] = vertex_radius[v]*sin(vertex_argument[v]);
		if (abs(vcoords[v][0]) < ZERO_ERROR)
			vcoords[v][0] = 0;
		if (abs(vcoords[v][1]) < ZERO_ERROR)
			vcoords[v][1] = 0;

		argument += two_pi/num_boundary_edges;
		interior_vertex[v] = 0;
	}	

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "draw_convex_triangulation: initial vertex_coordinates: " << endl;
    print(vcoords,debug,15,"draw_convex_triangulation ");
    debug << "draw_convex_triangulation: max_circle_packing_iterations = " << max_circle_packing_iterations << endl;
}

//	for (int i=0; i< num_interior_vertices; i++)
	for (int i=0; i < (max_circle_packing_iterations < 1000? max_circle_packing_iterations:num_interior_vertices); i++)
	{
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: search for interior vertex " << i << ", interior_vertex flags = ";					
	for (int j=0; j<num_vertices; j++)
		debug << interior_vertex[j] << ' ';
	debug << endl;
}
		bool next_interior_vertex_found = false;
		/* find a boundary edge connected to an interior vertex */
		for (unsigned int j=0; j< boundary.size() && !next_interior_vertex_found; j++)
		{
			int v1 = boundary[j];
			int v2 = (j == boundary.size()-1? boundary[0]: boundary[j+1]);

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: check boundary vertices " << v1 << " and " << v2 << endl;					
			
			for (int k=0; k< num_faces && !next_interior_vertex_found; k++)
			{
				int boundary_edge_vertex_count=0;
				if (face_matrix[k][0] == v1 || face_matrix[k][0] == v2)
					boundary_edge_vertex_count++;
				if (face_matrix[k][1] == v1 || face_matrix[k][1] == v2)
					boundary_edge_vertex_count++;
				if (face_matrix[k][2] == v1 || face_matrix[k][2] == v2)
					boundary_edge_vertex_count++;
					
				if (boundary_edge_vertex_count == 2)
				{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: found boundary vertices " << v1 << " and " << v2 << " in row " << k << " of face_matrix" << endl;					

					/* is the other vertex in the face an interior vertex? */
					int other_vertex = -1;
					
					if (face_matrix[k][0] != v1 && face_matrix[k][0] != v2)
						other_vertex = face_matrix[k][0];
					else if (face_matrix[k][1] != v1 && face_matrix[k][1] != v2)
						other_vertex = face_matrix[k][1];
					else // if (face_matrix[k][2] != v1 && face_matrix[k][2] != v2)
						other_vertex = face_matrix[k][2];
						
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: other_vertex = " << other_vertex << endl;					
	
					if (interior_vertex[other_vertex-1] == 1)
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: other_vertex is the next interior vertex" << endl;					
	
						next_interior_vertex_found = true;
						
						/* We are going to perform a type 1 collapse on edge v1 v2, so record this
						   by clearing the corresponding elements of local_edge_matrix
						*/
						local_edge_matrix[v1-1][v2-1]=0;
						local_edge_matrix[v2-1][v1-1]=0;

						/* If the number of edges joining the new interior vertex to boundary vertices is greater 
						   than 2, we set the initial placement of the new vertex to be the centre of gravity of
						   these boundary vertices.  Otherwise, we set it's argument to be midway between that of 
						   v1 and v2, and it's radius to be smaller than the minimum radius of v1 and v2.
						*/
						int boundary_valency=0;
						pair <double,double> cog(0.0,0.0);
						for (unsigned int l=0; l< boundary.size(); l++)
						{
							if (edge_matrix[other_vertex-1][boundary[l]-1] == 1)
							{
								cog.first +=vcoords[boundary[l]-1][0];
								cog.second +=vcoords[boundary[l]-1][1];								
								boundary_valency++;
							}
						}						

if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: boundary valence of vertex " << other_vertex << " is " << boundary_valency << endl;
	
/* DISABLE THIS FOR NOW */
						if (false && boundary_valency > 2 && boundary_valency != boundary.size())
						{
							vcoords[other_vertex-1][0] = cog.first/boundary_valency;
							vcoords[other_vertex-1][1] = cog.second/boundary_valency;
							/* we scale the location for the call to edge_argument to be consistent with its use by force_direction */
							vertex_argument[other_vertex-1] = edge_argument(vcoords[other_vertex-1][0]*metapost_coordinate_scale_factor,
							                                                vcoords[other_vertex-1][1]*metapost_coordinate_scale_factor);
							vertex_radius[other_vertex-1] = sqrt(vcoords[other_vertex-1][0]*vcoords[other_vertex-1][0] + vcoords[other_vertex-1][1]*vcoords[other_vertex-1][1]);
						}
						else
						{
							
							/* set the argument of the other_vertex to be half way between v1 and v2. */
							vertex_argument[other_vertex-1] = (vertex_argument[v1-1]+vertex_argument[v2-1])/2;
							
							/* set the initial radius of the other vertex to be min(v1 v2) - epsilon */

//When there are a large number if interior vertice, using a multiplicative epsilon here forces vertices to the middle too soon
//							vertex_radius[other_vertex-1] = min(vertex_radius[v1-1],vertex_radius[v2-1])*epsilon; 

							vertex_radius[other_vertex-1] = min(vertex_radius[v1-1],vertex_radius[v2-1])-epsilon; 
	
							vcoords[other_vertex-1][0] = vertex_radius[other_vertex-1]*cos(vertex_argument[other_vertex-1]);
							vcoords[other_vertex-1][1] = vertex_radius[other_vertex-1]*sin(vertex_argument[other_vertex-1]);
							if (abs(vcoords[other_vertex-1][0]) < ZERO_ERROR)
								vcoords[other_vertex-1][0] = 0;
							if (abs(vcoords[other_vertex-1][1]) < ZERO_ERROR)
								vcoords[other_vertex-1][1] = 0;
								
							/* check that this placement of other vertex isn't outside the line v1 v2 and adjust it if it is */
							pair<double,double> check = intersection_point (v1, v2, other_vertex, vcoords);						
							double check_radius = vector_modulus(check.first,check.second);
							
							if (vertex_radius[other_vertex-1] > check_radius)
								vertex_radius[other_vertex-1] = 0.8*check_radius;
							
						}
						
if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: argument of vertex " << other_vertex << " = " << vertex_argument[other_vertex-1] << endl;
	debug << "draw_convex_triangulation: initial radius of vertex " << other_vertex << " set to " << vertex_radius[other_vertex-1] << endl;
	debug << "draw_convex_triangulation: initial coordinates (" << vcoords[other_vertex-1][0] << ',' << vcoords[other_vertex-1][1] << ')' << endl;
}	      
						
						/* The radius of the next interior vertex has to be chosen so that its placement does 
						   not violate the integrity of the triangulation.  Specifically, it must lie in the 
						   same closed half space as v1 and v2 for each half space determined by an interior 
						   edge or the line v1 v2.  An interior edge is an edge joining a vertex to a 
						   non-adjacent vertex in the boundary.  Note we must consider all interior edges, 
						   not just those attached to v1 and v2.  We check for interior using local_edge_matrix, 
						   since this has had edges removed that used to be in the boundary but have subsequently 
						   been collapsed.
						   						   
						*/
//bool adjustment_required = true;						
						for (unsigned int l=0; l < boundary.size()-2; l++)
						for (unsigned int m=l+2; m < (l==0? boundary.size()-1:boundary.size()); m++)
						{
							int bl = boundary[l];
							int bm = boundary[m];
							if (local_edge_matrix[bl-1][bm-1] == 1)
							{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: interior edge " << bl << ' ' << bm << " found" << endl;			
	
								/* Given two points a=(a_x,a_y) and b=(b_x,b_y) in R^2 the straight line joining a and b is given 
								   by the equation
								   
								           y = g (x-a_x) + a_y where g = (b_y - a_y)/(b_x - a_x) 
								   
								   (shift the line parallel to the y-axis by a_y, then parallel to the x-axis by a_x to translate 
								   the point a to the origin.  This moves the point b to b' = (b_x - a_x, b_y - a_y) whence the 
								   gradient of the line may be seen to be g = (b_y - a_y)/(b_x - a_x).  
								   Thus the equation of the translated line is y = gx.)
								   
								   For u=(x,y) put alpha_u = y - g (x-a_x) - a_y, then for all points u lying on the line, 
								   alpha_u = 0 and for points u, v not lying on the line alpha_u and alpha_v are non zero 
								   and have the same sign iff u and v lie in the same component of the line's 
								   complement in R^2.  Thus, u and v lie on different sides of the line iff 
								   alpha_u * alpha_v < 0.
								   
								   In our case, we know v1 and v2 lie in the same closed half space of the complement of the line
								   through bl and bm, so their alpha values are either of the same sign, or one of them is zero.
								   We may therefore add their alpha values for the purpose of comparison with the placement of 
								   other_vertex.
								*/	
								
								double alpha_v1 = alpha_value(vcoords[bl-1][0], vcoords[bl-1][1],vcoords[bm-1][0], vcoords[bm-1][1],
								                              vcoords[v1-1][0], vcoords[v1-1][1]);
								double alpha_v2 = alpha_value(vcoords[bl-1][0], vcoords[bl-1][1],vcoords[bm-1][0], vcoords[bm-1][1],
								                              vcoords[v2-1][0], vcoords[v2-1][1]);

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: alpha values of v1 and v2 with respect to interior edge " << bl << ' ' << bm 
	      << " are " << alpha_v1 << ' ' << alpha_v2 << endl;			
}
								bool adjustment_required = true;
//								adjustment_required = true;
								do
								{
									double alpha_other = alpha_value(vcoords[bl-1][0], vcoords[bl-1][1],vcoords[bm-1][0], vcoords[bm-1][1],
								                              vcoords[other_vertex-1][0], vcoords[other_vertex-1][1]);
								                              
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: alpha value for vertex " << other_vertex << " = " << alpha_other << endl;

								    if (alpha_other *(alpha_v1+alpha_v2) < 0)
								    {
										/* current placement of other_vertex is the wrong side of the interior edge 
										   change r to r + (1-r)/3
										vertex_radius[other_vertex-1] = (2*vertex_radius[other_vertex-1]+1)/3;
										vcoords[other_vertex-1][0] = vertex_radius[other_vertex-1]*cos(vertex_argument[other_vertex-1]);
										vcoords[other_vertex-1][1] = vertex_radius[other_vertex-1]*sin(vertex_argument[other_vertex-1]);
										if (abs(vcoords[other_vertex-1][0]) < ZERO_ERROR)
											vcoords[other_vertex-1][0] = 0;
										if (abs(vcoords[other_vertex-1][1]) < ZERO_ERROR)
											vcoords[other_vertex-1][1] = 0;
										*/

										/* current placement of other_vertex is the wrong side of the interior edge, change it to the 
										   mid point between the intersections of the ray through other_vertex and the two edges bl bm
										   and v1 v2
										*/
										pair<double,double> new_location = midpoint_between_edges(v1, v2, bl, bm, other_vertex, vcoords);

										vcoords[other_vertex-1][0] = new_location.first;
										vcoords[other_vertex-1][1] = new_location.second;
										/* we scale the location for the call to edge_argument to be consistent with its use by force_direction */
										vertex_argument[other_vertex-1] = edge_argument(vcoords[other_vertex-1][0]*metapost_coordinate_scale_factor,
										                                                vcoords[other_vertex-1][1]*metapost_coordinate_scale_factor);
										vertex_radius[other_vertex-1] = sqrt(vcoords[other_vertex-1][0]*vcoords[other_vertex-1][0] + vcoords[other_vertex-1][1]*vcoords[other_vertex-1][1]);
										

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: vertex " << other_vertex << " currently placed the wrong side of the interior edge " 
	      << bl << ' ' << bm << endl;
	debug << "draw_convex_triangulation: new radius = " << vertex_radius[other_vertex-1] << endl;
	debug << "draw_convex_triangulation: new coordinates (" << vcoords[other_vertex-1][0] << ',' << vcoords[other_vertex-1][1] << ')' << endl;

 alpha_other = alpha_value(vcoords[bl-1][0], vcoords[bl-1][1],vcoords[bm-1][0], vcoords[bm-1][1],
								                              vcoords[other_vertex-1][0], vcoords[other_vertex-1][1]);
	debug << "draw_convex_triangulation: alpha_other = " << alpha_other << endl;								                              
}	      

//adjustment_required = false;

									}
									else
									{
										/* current placement of other_vertex is OK */
										adjustment_required = false;
									}
									
								} while (adjustment_required);
								
			
							}
						}












if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "draw_convex_triangulation: argument for vertex " << other_vertex << " is " << vertex_argument[other_vertex-1]
          << " radius is " << vertex_radius[other_vertex-1] << endl;
    debug << "draw_convex_triangulation: updated vertex_coordinates: " << endl;
    print(vcoords,debug,15,"draw_convex_triangulation ");
}
						
						/* insert other vertex into the boundary between v1 and v2 */
						vector<int>::iterator bptr = boundary.begin();
						while (bptr != boundary.end())
						{
							if (*bptr == v1 || *bptr == v2)
							{
								bptr++;
								boundary.insert(bptr,other_vertex); // before bptr

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: boundary updated to: ";					
	for (unsigned int b=0; b< boundary.size(); b++)
		debug << boundary[b] << ' ';
	debug << endl;
}	
								interior_vertex[other_vertex-1] = 0;
								break;
							}
							else
							{
								bptr++;
							}
						}
					}
					else
					{
if (debug_control::DEBUG >= debug_control::BASIC)
	debug << "draw_convex_triangulation: other_vertex has already been placed" << endl;					
					}
				}
			}			
		}
	}
	

if (debug_control::DEBUG >= debug_control::BASIC)
{
	debug << "draw_convex_triangulation: final vertex radii ";					
	for (int i=0; i< num_vertices; i++)
		debug << vertex_radius[i] << ' ';
	debug << endl;
	debug << "draw_convex_triangulation: final vertex arguments ";					
	for (int i=0; i< num_vertices; i++)
		debug << vertex_argument[i] << ' ';
	debug << endl;
}	
	
	write_convex_metapost(output, vertex_radius, vertex_argument, filename, edge_matrix, mp_control);
}


void write_convex_metapost(ofstream& os, vector<double> vertex_radius, vector<double> vertex_argument, string filename,
     matrix<int>& edge_matrix, metapost_control& mp_control)
{

	os << "\n\n";
	
//	if (title.length())
	os << "% convex triangulation from " << filename << endl;

	
	/* controls for coordinate output */
	int num_decimal_points = 3;
	int output_field_width = 12;
	os.setf(ios::fixed,ios::floatfield);
	os.precision(num_decimal_points);

	os << "\nbeginfig(fignum);" << endl;
	os << "fignum:=fignum+1;" << endl;
	os << "numeric u,d;" << endl;
	float unit_points = mp_control.unit_size;
	os << "u=" << unit_points/100 << "pt;" << endl;
	os << "d=" << mp_control.disc_size << "u;" << endl;
	os << "path p[];" << endl;
	os << "pickup pencircle scaled " << mp_control.pen_size*0.5 << "pt;" << endl;

	
	/* vertex coordinates are numbered z<vertex> */
	int num_vertices = vertex_radius.size();
	
	matrix<int> vcoords(num_vertices,2);
		
	for (int i=0; i< num_vertices; i++)
	{
		vcoords[i][0] = vertex_radius[i]*cos(vertex_argument[i])*metapost_coordinate_scale_factor;
		vcoords[i][1] = vertex_radius[i]*sin(vertex_argument[i])*metapost_coordinate_scale_factor;
	}

if (debug_control::DEBUG >= debug_control::BASIC)
{
    debug << "write_convex_metapost: vertex_coordinates: " << endl;
    print(vcoords,debug,6,"write_convex_metapost: ");
}

	os << "draw fullcircle scaled " << 2*metapost_coordinate_scale_factor << "u dashed evenly;";
	
	for (int i=0; i< num_vertices; i++)
		os << "z" << i << "=(" << setw(output_field_width) << vcoords[i][0] << "u," << setw(output_field_width) << vcoords[i][1] << "u);" << endl;

	for (int i=0; i< num_vertices; i++)
	for (int j=0; j< i; j++)
	{
		if (edge_matrix[i][j] == 1)
			os << "draw z" << i << "--z" << j << ";" << endl;
	}

	if (mp_control.label_vertices)
	{
		
		/* evaluate the min/max vcoords values */
		double minx = vcoords[0][0];
		double maxx = vcoords[0][0];
		double miny = vcoords[0][1];
		double maxy = vcoords[0][1];
		
		for (int i=1; i< num_vertices; i++)
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
		
		if (mp_control.label_vertices)
		{
			for (int i=0; i< num_vertices; i++)
			{
				if (vcoords[i][0] !=0 || vcoords[i][1] !=0)
					os << "label(btex $z" << i+1 << "$ etex, z" << i << ");" << endl;
			}
		}
		os << "draw (0," << miny << "u)--(0," << maxy << "u) dashed evenly;" << endl;
		os << "draw (" << minx << "u,0)--(" << maxx << "u,0) dashed evenly;" << endl;
/*			
		os << "fill fullcircle scaled 0.5d shifted (0," << miny << "u);" << endl;
		os << "fill fullcircle scaled 0.5d shifted (0," << maxy << "u);" << endl;
		os << "fill fullcircle scaled 0.5d shifted (" << minx << "u,0);" << endl;
		os << "fill fullcircle scaled 0.5d shifted (" << maxx << "u,0);" << endl;

		os << "label.lft(btex (0," << int(miny) << "u) etex,(0," << miny << "u));" << endl;
		os << "label.lft(btex (0," << int(maxy) << "u) etex,(0," << maxy << "u));" << endl;
		os << "label.bot(btex (" << int(minx) << "u,0) etex,(" << minx << "u,0));" << endl;
		os << "label.bot(btex (" << int(maxx) << "u,0) etex,(" << maxx << "u,0));" << endl;
*/
	}
	
	os << "endfig;" << endl;
}

/* alpha_value determines where the point (x,y) lies in relation to the straight line through
   the (a_x,a_y) and (b_x,b_y).  It evaluates the quantity:
   
   alpha := y - g (x-a_x) - a_y where g = (b_y - a_y)/(b_x - a_x) 
   
   alpha is zero iff (x,y) lies on the line and is positive or negative according to which side
   of the line (x,y) lies.
*/
double alpha_value(double a_x, double a_y, double b_x, double b_y, double x, double y)
{
	/* scale coordinates by metapost_coordinate_scale_factor 
	a_x *= metapost_coordinate_scale_factor;
	a_y *= metapost_coordinate_scale_factor;
	b_x *= metapost_coordinate_scale_factor;
	b_y *= metapost_coordinate_scale_factor;
	x *= metapost_coordinate_scale_factor;
	y *= metapost_coordinate_scale_factor;
*/	
	double gradient = (b_y-a_y)/(b_x-a_x);
	double on_line_y_value = (x-a_x)*gradient + a_y;
	double alpha = y - on_line_y_value;
	
if (debug_control::DEBUG >= debug_control::DETAIL)
{
	debug << "alpha_value: alpha value of (" << x << "," << y << ") and line joining (" << a_x << "," << a_y << ") and (" << b_x << "," << b_y << ") is " << alpha << endl;
	debug << "alpha_value:   gradient = " << gradient << ", " << "on-line y-value for x= " << x << " is " << on_line_y_value << endl;
}	
	return alpha;
}


pair<double,double> midpoint_between_edges(int v1, int v2, int bl, int bm, int other_vertex, matrix<double>& vcoords)
{
	pair<double,double> boundary_intersection_point = intersection_point (v1, v2, other_vertex, vcoords);
	pair<double,double> interior_intersection_point = intersection_point (bl, bm, other_vertex, vcoords);
	
	pair<double,double> midpoint;
	
	midpoint.first = (boundary_intersection_point.first+interior_intersection_point.first)/2;
	midpoint.second = (boundary_intersection_point.second+interior_intersection_point.second)/2;
	
	return midpoint;
}

/* intersection_point calculates the intersection between the line through v1 v2 and the radial line from the origin 
   through other_vertex.  The line between v1 and v2 is 
   
   y = g (x - v1_x) + v1_y, where g = (v2_y - v1_y)/(v2_x - v1_x)

   and the line through other_vertex, with coordinates (u,v) say, is 
   
   y = v/u x
   
   We solve for their intersection using the general matrix equation
   
   a  b   x   =   p
   c  d   y       q
   
   
   which has solution
   
   x    =      1     d -b   p
   y         ad-bc  -c  a   q
   
   In our case, we have 
   
    g   -1   x   =  g*v1_x - v1_y
   v/u  -1   y   =  0
   
   That is, b=d=-1 and q=0, so our solution is
   
   x    =    1   -1  1   p   =   1     -p  =  p    1  =  ( p/(a-c), cp/(a-c) )
   y        c-a  -c  a   0      c-a   -cp    a-c   c
   
   where a = g, c = v/u and p = g*v1_x - v1_y   

*/
pair<double,double> intersection_point (int v1, int v2, int other_vertex, matrix<double>& vcoords)
{
	double v1_x = vcoords[v1-1][0];
	double v1_y = vcoords[v1-1][1];
	double v2_x = vcoords[v2-1][0];
	double v2_y = vcoords[v2-1][1];
	double u = vcoords[other_vertex-1][0];
	double v = vcoords[other_vertex-1][1];

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "intersection_point: v1_x = " << v1_x << ", v1_y = " << v1_y << ", v2_x = " << v2_x << ", v2_y = " << v2_y << ", u = " << u << ", v = " << v << endl;
		
	double g = (v2_y - v1_y)/(v2_x - v1_x);

	double a = g;
	double c = v/u;
	double p = g*v1_x-v1_y;

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "intersection_point: a=g = " << g << ", c = " << c << ", p = " << p << endl;
	
	pair<double,double> intersection;
	
	intersection.first = p/(a-c);
	intersection.second = c*p/(a-c);

if (debug_control::DEBUG >= debug_control::DETAIL)
	debug << "intersection_point: intersection = " << intersection.first << ',' << intersection.second << endl;
	
	return intersection;
}

double vector_modulus (double x, double y)
{
	return sqrt(x*x+y*y);
}
