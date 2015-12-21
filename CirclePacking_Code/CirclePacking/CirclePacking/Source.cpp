#include <iostream>
#include <stdlib.h>
#include <GL/freeglut.h>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;
#define PI 3.1415926535897932384626433832795


class Point2D
{
public:
	double x;
	double y;
};

class Disk
{
public:
	Point2D center;
	double radius;
	int CenterIdx; //disc index
	double *ang;
	bool interior_chk; //to check whether the vertices are on boundary or interior
	bool locatedVert; //whether a vertex found its position
	int firstnbdId;
	bool nbdPosFound; //Is a vertex have found all its neighboring veritices postion(locating disc)
	vector < std::pair<double, int>> vecSort;
};

class Complex
{
public:
	Point2D center;
	double radius;
	int CenterIdx;
};

int nvertics, nrelations; //number of vertices, numberof edges 
Disk *circles;
Complex *complex;
//contains informations of relations of each vertex.(i.e, edges) 
int *a1; 
int *a2;

//to draw a 2D circle
void circle(float x, float y, float r, int segments)
{
	glBegin(GL_LINE_STRIP);
	for (int n = 0; n <= segments; ++n) {
		float const t = 2 * PI * (float)n / (float)segments;
		glVertex2f(x + sin(t) * r, y + cos(t) * r);
	}
	glEnd();
}

void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPointSize(4);
	glColor3ub(0, 0, 0);
	glBegin(GL_POINTS);
	for (int i = 0; i < nvertics; i++)
	{
		glVertex2f(circles[i].center.x, circles[i].center.y);
		//glVertex2f(complex[i].center.x, complex[i].center.y);
	}
	glEnd();
	//glBegin(GL_LINES);
	//for (int i = 0; i < nrelations; i++)
	//{
	//	glVertex2f(complex[a1[i]].center.x, complex[a1[i]].center.y);
	//	glVertex2f(complex[a2[i]].center.x, complex[a2[i]].center.y);
	//}
	//glEnd();

	for (int i = 0; i < nvertics; i++)
	{
		glColor3ub(0, 100, 200);
		circle(circles[i].center.x, circles[i].center.y, circles[i].radius, 16);
	}
	glFinish();
}

//notation abused
double innerprod(Point2D v1, Point2D v2)
{
	return(v1.x*v2.x + v1.y*v2.y);
}

//compute angle btw v1,v2 and x-axis
//it is used to sort the edges in angle increasing order
double computeAng(Point2D p1, Point2D p2)
{
	Point2D p21; p21.x = p2.x - p1.x; p21.y = p2.y - p1.y;
	Point2D p31; p31.x = 1;			  p31.y = 0;

	double norm1 = sqrt(p21.x*p21.x + p21.y*p21.y);
	double norm2 = sqrt(p31.x*p31.x + p31.y*p31.y);
	double Ndirection = (-p21.x*p31.y) - (-p21.y*p31.x);

	if (Ndirection > 0) return(acos(innerprod(p21, p31) / norm1 / norm2));
	else return(2 * PI - acos(innerprod(p21, p31) / norm1 / norm2));//to consider the angle bigger than pi
}

void ReadFile(const char *file)
{
	FILE *fp = fopen(file, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "Model Constructor: Couldn't open %s\n", file);
		exit(-1);
	}

	fscanf(fp, "%d", &nvertics);
	fscanf(fp, "%d", &nrelations);

	circles = new Disk[nvertics]; 
	complex = new Complex[nvertics];
	double i1, i2;
	
	for (int i = 0; i < nvertics; i++)
	{
		fscanf(fp, "%lf %lf", &i1, &i2);
		//used in drawing circles with UNM method
		circles[i].center.x = i1;
		circles[i].center.y = i2;
		circles[i].CenterIdx = i;
		//I store the information twice to draw complex 
		complex[i].center.x = i1;
		complex[i].center.y = i2;
	}
	 
	a1 = new int[nrelations];
	a2 = new int[nrelations];
	//read edge connection information
	for (int j = 0; j < nrelations; j++)
	{
		fscanf(fp, "%d %d", &a1[j], &a2[j]);//nbd vertex indx
	}
	fclose(fp);


	//specify rooms for each nbd vertex
	for (int i = 0; i < nvertics; i++)
	{
		circles[i].interior_chk = TRUE;
		circles[i].radius = 6;//initial label
	}

	//----------------------------------------------------------------------------
	//--------------------------sort neighbor vertex------------------------------
	//----------------------------------------------------------------------------

	for (int j = 0; j < nrelations; j++)
	{
		double a = computeAng(circles[a1[j]].center, circles[a2[j]].center);
		//storing edges with the angles between edge and x-axis
		circles[a1[j]].vecSort.push_back(std::make_pair(a, a2[j]));
		circles[a2[j]].vecSort.push_back(std::make_pair(a < PI? a+PI : a-PI, a1[j]));
	}


	//sort the neighboring vertices' index
	for (int i = 0; i < nvertics; i++)
	{
		std::sort(circles[i].vecSort.begin(), circles[i].vecSort.end());
	}

	for (int i = 0; i < nvertics; i++)
	{
		circles[i].ang = new double[circles[i].vecSort.size()];
	}
}

//check whether the vertex is interior or on boundary
//in this part, its efficiency is bad 
void interior_CHK()
{
	for (int i = 0; i < nvertics; i++)
	{
		for (int j = 0; j < circles[i].vecSort.size(); j++)
		{
			//comapring only the number of edges and the number of angles are not enough because of the special case
			double N = (circles[circles[i].vecSort[j].second].center.x - circles[i].center.x)*(circles[circles[i].vecSort[(j + 1) % circles[i].vecSort.size()].second].center.y - circles[i].center.y) - (circles[circles[i].vecSort[(j + 1) % circles[i].vecSort.size()].second].center.x - circles[i].center.x)*(circles[circles[i].vecSort[j].second].center.y - circles[i].center.y);
			bool check=FALSE;
			if (N > 0)
			{
				for (int k = 0; k < circles[circles[i].vecSort[j].second].vecSort.size(); k++)
				{
					//if two vertices are connected and n>0 then we can determine what is interior
					if (circles[circles[i].vecSort[j].second].vecSort[k].second == circles[i].vecSort[(j + 1) % circles[i].vecSort.size()].second) 
					{
						check = TRUE; 
						circles[i].firstnbdId = 0;
						break;
					}
				}
			}
			if (!check)
			{
				circles[i].interior_chk = FALSE;
				circles[i].firstnbdId = (j + 1) % circles[i].vecSort.size();; break; 
			}
		}
	}

	printf("---------Boundary Vertex index--------- \n");
	for (int i = 0; i < nvertics; i++)
	{
		if (circles[i].interior_chk == FALSE)
		{
			printf("%d", circles[i].CenterIdx);
			printf("\n");
		}
	}

}

//given circles which are tangent to each other,
//the three radii information will give an angle of the triangle
double angSum(double r, double r1, double r2)
{
	return (2 * asin(sqrt(r1*r2 / ((r + r1)*(r + r2)))));
}


//----------------------------------------------------------------------------
//-------------------------Uniform Neighbor Model-----------------------------
//----------------------------------------------------------------------------
double UNM(Disk *d)
{
	double delta = sin(PI / d->vecSort.size()); //sin(2pi/2k);
	double angleSumOfeachVert = 0.0;
	for (int i = 0; i < d->vecSort.size(); i++)
	{
		//compute angle sum of each vertex, the goal angle sum is 2pi for each vertex.
		angleSumOfeachVert += angSum(d->radius, circles[d->vecSort[i].second].radius, circles[d->vecSort[(i + 1) % d->vecSort.size()].second].radius); 
	}
	//detailed description can be found in the picture of the report
	double beta = sin(angleSumOfeachVert / (2 * d->vecSort.size()));
	double rhat = beta *d->radius / (1 - beta);

	double tol = abs(d->radius - (1 - delta)*rhat / delta);
	d->radius = (1 - delta)*rhat / delta;
	return tol;
}

void runUNM()
{
	interior_CHK();
	//boundary condition specified
	for (int i = 0; i < nvertics; i++)
	{
		if (circles[i].interior_chk == FALSE)//boundary conditions specified to give a unique circle packing
		{
			circles[i].radius = 25;
		}
		circles[i].nbdPosFound = FALSE;//initialization
	}
	double error;
	do{
		error = 0.0;
		for (int i = 0; i < nvertics; i++)
		{
			if (circles[i].interior_chk == TRUE)
			{
				error += UNM(&circles[i]); //error is defined as the summation of the difference of radius at interior vertices
			}
		}
	} while (error > 0.01);
	printf("radius \n");
	for (int i = 0; i < nvertics; i++)
	{
		printf("%f \n", circles[i].radius);
		circles[i].locatedVert = FALSE;//initialize flag. to know what vertices are located
		circles[i].center.x = 10;
		circles[i].center.y = 10;//initialize vertex center
	}
}

//input is an index of a disc
//find appropriate postions of circles iteratively (after UNM process ended)
void find_nbd_Position(int idx)
{
	int remember;//find the vertex, located in a previous time. then set, v[remember]=standard vertex
	//to locate circles
	for (remember = 0; remember < circles[idx].vecSort.size(); remember++)
	{
		if (circles[circles[idx].vecSort[remember].second].locatedVert)
		{
			break; 
		}
	}

	double angle = computeAng(circles[idx].center, circles[circles[idx].vecSort[remember].second].center);

	//I will describe it with detail and picutre in the report
	//loop through as much as the number of edges
	int cnt = 0;
	for (int j = remember; cnt < circles[idx].vecSort.size(); j = (j + 1) % circles[idx].vecSort.size())
	{
		if (circles[circles[idx].vecSort[j].second].locatedVert == FALSE )
		{
			circles[circles[idx].vecSort[j].second].center.x = circles[idx].center.x + (circles[idx].radius + circles[circles[idx].vecSort[j].second].radius)*cos(angle);
			circles[circles[idx].vecSort[j].second].center.y = circles[idx].center.y + (circles[idx].radius + circles[circles[idx].vecSort[j].second].radius)*sin(angle);
			circles[circles[idx].vecSort[j].second].locatedVert= TRUE;
		}
		angle += circles[idx].ang[j];
		cnt++;
	}
	circles[idx].nbdPosFound = TRUE;//a vertex with index [idx] found all its neighbor circles' center postion

	for (int j = 0; j < circles[idx].vecSort.size(); j++)
	{
		//if not all neighbor center positions are located. Find the positions of circles
		if (circles[circles[idx].vecSort[j].second].nbdPosFound == FALSE)
		{ 
			if (circles[circles[idx].vecSort[j].second].interior_chk == TRUE)
			{
				find_nbd_Position(circles[idx].vecSort[j].second);
			}
			else{ break; }
		}
	}
}

void LocateDisks()
{
	//angles of each vertex will be used to locate circles(above function) 
	for (int i = 0; i < nvertics; i++)
	{
		for (int j = 0; j <circles[i].vecSort.size(); j++)
		{
			circles[i].ang[j]= angSum(circles[i].radius, circles[circles[i].vecSort[j].second].radius, circles[circles[i].vecSort[(j + 1) % circles[i].vecSort.size()].second].radius);
		}
	}
	//fix first center of a vertex
	double center_x = 100;
	double center_y = 260;
	//arbitrarily selet a vertex to locate first
	int id1 = 11;//if I change this as other index. This algorithm fails.

	circles[id1].center.x = center_x;
	circles[id1].center.y = center_y;
	circles[id1].locatedVert = TRUE;
	int id2 = circles[id1].vecSort[0].second;

	double angle = 0.0;
	//First, locate first vertex and its neighboring vertices 
	for (int j = 0; j < circles[id1].vecSort.size(); j++)
	{
		circles[circles[id1].vecSort[j].second].center.x = circles[id1].center.x + (circles[id1].radius + circles[circles[id1].vecSort[j].second].radius)*cos(angle);
		circles[circles[id1].vecSort[j].second].center.y = circles[id1].center.y + (circles[id1].radius + circles[circles[id1].vecSort[j].second].radius)*sin(angle);
		angle += circles[id1].ang[j];
		circles[circles[id1].vecSort[j].second].locatedVert = TRUE;
	}
	circles[id1].nbdPosFound = TRUE;
	circles[id2].locatedVert = TRUE;
	find_nbd_Position(id2);
}


int main(int argc, char ** argv)
{
	glutInit(&argc, argv);
	ReadFile("circles1.txt");
	runUNM();
	LocateDisks();

	int width = 500;
	int height = 500;

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(180, 100);

	glutCreateWindow("Circle Packing");
	glEnable(GL_DEPTH_TEST);

	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-100, width, -100, height, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glutDisplayFunc(myDisplay);

	glutMainLoop();
	//Sleep(5000);

	return 0;
}
