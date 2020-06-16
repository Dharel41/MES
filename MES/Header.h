#include<conio.h>
#include<iostream>
using namespace std;

struct point
{
	friend ostream &operator<<(ostream &wyjscie, const point &a);
	int id;
	double x=0, y=0;
};
struct node
{
	friend ostream &operator<<(ostream &wyjscie, const node &a);
	double x, y, t0 = 100;
	int id, status = 0;
};

struct element
{
	friend ostream &operator<<(ostream &wyjscie, const element &a);
	int id;
	node points[4];
	point new_coordinates[4];
};

struct grid
{
	int nheight, nlength;
	friend ostream &operator<<(ostream &wyjscie, const grid &a);
	node *n1;
	element *e1;
	grid(node *a, element *b, int c, int d)
	{
		n1 = a;
		e1 = b;
		nheight = c;
		nlength = d;
	}
};