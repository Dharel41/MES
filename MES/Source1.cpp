#include<conio.h>
#include<iostream>
#include"Header.h"
using namespace std;

ostream nod&operator<<(ostream &wyjscie, const node &a)
{

	wyjscie << "Wezel nr." << a.id << ".Wspolrzedne wezla (" << a.x << "," << a.y << "). Temperatura w wezle " << a.t0 << ". Status: " << a.status << endl;
	return wyjscie;
}

ostream &operator<<(ostream &wyjscie, const element &a)
{
	wyjscie << "Element " << a.id << " sklada sie z wezlow: " << a.points[0].id << " " << a.points[1].id << " " << a.points[2].id << " " << a.points[3].id << " " << endl;
	return wyjscie;
}

ostream &operator<<(ostream &wyjscie, const grid &a)
{
	for (int i = 0;i<(a.nheight*a.nlength);i++)
	{
		wyjscie << a.n1[i];
	}
	for (int i = 0;i<(a.nheight - 1)*(a.nlength - 1);i++)
	{
		wyjscie << a.e1[i];
	}
	return wyjscie;
}



