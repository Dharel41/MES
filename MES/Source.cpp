#include<conio.h>
#include<iostream>
#include"Header.h"
#include <iomanip>
using namespace std;
grid mesh(int nheight, int nlength, double height, double length);
double*** f_matrix_H(int nheight, int nlength, double height, double length,grid a);
double*** f_matrix_H_BC(int nheight, int nlength, double height, double length,grid a);
double*** f_matrix_C(int nheight, int nlength, double height, double length,grid a);
double** f_vector_P(int nheight, int nlength, double height, double length, grid a);
int main()
{
	double height, length, dTau, To, simulation_time, step_time;
	double ***Matrix_h, ***Matrix_c,***Matrix_h_bc,**vector_p;
	int nheight, nlength;
	height = 0.1;
	length = 0.1;
	nheight =4;
	nlength = 4;
	dTau = 50;
	To = 100;
	simulation_time = 500;
	step_time = 50;

	if (nheight <= 1 || nlength <= 1)
	{
		cout << "Nie mozna zbudowac siatki" << endl;
		_getch();
	}
	grid a = mesh(nheight, nlength, height, length);

	Matrix_h=f_matrix_H(nheight, nlength, height, length,a);
	Matrix_h_bc = f_matrix_H_BC(nheight, nlength, height, length, a);
	Matrix_c = f_matrix_C(nheight, nlength, height, length,a);
	vector_p = f_vector_P(nheight, nlength, height, length, a);

	/*

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
	   for (int j = 0;j<4;j++)
	    {

	    cout << vector_p[i][j] << "   ";

	    }
	cout << endl << endl;
	}



	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{

				cout << Matrix_h[i][j][k] << "   ";

			}
			cout << endl;

		}
		cout << endl << endl;
	}


	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{

				cout << Matrix_c[i][j][k] << "   ";

			}
			cout << endl;

		}
		cout << endl << endl;
	}
	*/
	double **global_H= new double *[nheight*nlength];
	double **global_C= new double *[nheight*nlength];
	double **global_H_BC = new double *[nheight*nlength];
	double **global_vector_P = new double *[nheight*nlength];
	for (int i = 0; i < nheight*nlength; i++)
	{
		global_H[i] = new double[nheight*nlength];
		global_C[i] = new double[nheight*nlength];
		global_H_BC[i] = new double[nheight*nlength];
	}
	global_vector_P[0] = new double[nheight*nlength];
	for (int i = 0; i < nheight*nlength; i++)
	{
			for (int j = 0;j < nheight*nlength;j++)
			{
				global_H[i][j] = 0;
				global_C[i][j] = 0;
				global_H_BC[i][j] = 0;
			}
			global_vector_P[0][i] = 0;
	}



	for (int i = 0;i <(nheight-1)*(nlength-1);i++)
	{

		for (int j = 0; j < 4; j++)
		{
			for (int k = 0;k < 4;k++)
			{
				global_H[a.e1[i].points[j].id][a.e1[i].points[k].id] += Matrix_h[i][j][k];
				global_C[a.e1[i].points[j].id][a.e1[i].points[k].id] += Matrix_c[i][j][k];
				global_H_BC[a.e1[i].points[j].id][a.e1[i].points[k].id] += Matrix_h_bc[i][j][k];
			}
			global_vector_P[0][a.e1[i].points[j].id] += vector_p[i][j];
		}

	}


	/*

	for (int i = 0; i < nheight*nlength; i++)
	{
		for (int j = 0;j < nheight*nlength;j++)
		{
			cout << setw(9) << global_H[i][j];
		}
		cout << endl;
	}

	cout << endl << endl;
	
	for (int i = 0; i < nheight*nlength; i++)
	{
		for (int j = 0;j < nheight*nlength;j++)
		{
			cout <<setw(9)<< global_C[i][j];
		}
		cout << endl;
	}
	*/

	//testcase

	double **Matrix_h_test= new double*[nheight*nlength];
	double **vector_p_test =new double*[nheight*nlength];
	for (int i = 0; i < nheight*nlength; i++)
	{
		Matrix_h_test[i] = new double[2*nheight*nlength];
		vector_p_test[i] = new double[nheight*nlength];
	}
	double *Ct = new double[nheight*nlength];
	double *To_tab = new double[nheight*nlength];
	double *T1_tab = new double[nheight*nlength];
	for(int i=0;i<nheight*nlength;i++)
	{
		To_tab[i] = To;
		T1_tab[i] = 0;
	}

	int count = 1;
	for (dTau;dTau < 550;dTau=dTau + step_time)
	{
		cout << "Iteration " << count << "   T=" << dTau << endl;
		count++;
		for (int i = 0; i < nheight*nlength; i++)
		{
			for (int j = 0;j < nheight*nlength;j++)
			{
				Matrix_h_test[i][j]=global_H[i][j] + global_H_BC[i][j] + (global_C[i][j] / step_time);    //[H] = [H]+[C]/dT
			}
		}

		for (int i = 0; i < nheight*nlength; i++)
		{
			Ct[i] = 0;
		}

		for (int i = 0;i < nheight*nlength;i++)
		{
			for (int j = 0;j < nheight*nlength;j++)
			{
				Ct[i] += (global_C[i][j] / step_time) * To_tab[j];                  //{[C]/dT}*{T0}
			}
		}

		for (int i = 0;i < nheight*nlength;i++)
		{
			vector_p_test[0][i] = global_vector_P[0][i]+ Ct[i];             //{P} = {P}+{[C]/dT}*{T0}
		}
		/*
		for (int i = 0; i < nheight*nlength; i++)
		{
			for (int j = 0;j < nheight*nlength;j++)
			{
				cout <<setw(9)<< Matrix_h_test[i][j];
			}
			cout << endl;
		}
		cout << endl << endl;
		*/
		cout << "{P}   ";
		for (int i = 0; i < nheight*nlength; i++)
		{
			cout << vector_p_test[0][i] << "  ";
		}
		cout << endl;
		

		////////////////////////////////////macierz odwrotna///////////////////////////////////
		double a, ratio;
		for (int i = 0; i < nheight*nlength; i++) {
			for (int j = nheight*nlength; j < 2 * nheight*nlength; j++) {
				if (i == (j - (nheight*nlength)))
					Matrix_h_test[i][j] = 1.0;
				else
					Matrix_h_test[i][j] = 0.0;
			}
		}
		for (int i = 0; i < nheight*nlength; i++) {
			for (int j = 0; j < nheight*nlength; j++) {
				if (i != j) {
					ratio = Matrix_h_test[j][i] / Matrix_h_test[i][i];
					for (int k = 0; k < 2 * nheight*nlength; k++) {
						Matrix_h_test[j][k] -= ratio * Matrix_h_test[i][k];
					}
				}
			}
		}
		for (int i = 0; i <nheight*nlength; i++) {
			a = Matrix_h_test[i][i];
			for (int j = 0; j < 2 * nheight*nlength; j++) {
				Matrix_h_test[i][j] /= a;
			}
		}


	////////////////////////////////////macierz odwrotna///////////////////////////////////

		for (int i = 0;i < nheight*nlength;i++)
		{
			for (int j = 0;j < nheight*nlength;j++)
			{
				T1_tab[i] += Matrix_h_test[i][j+ nheight*nlength] * vector_p_test[0][j];                  //{[C]/dT}*{T0}
			}
		}
		cout << "Temperatury ";
		for (int i = 0;i < nheight*nlength;i++)
		{
			To_tab[i] = T1_tab[i];
			cout << To_tab[i]<<"  ";
			T1_tab[i] = 0;
		}
		cout << endl << endl;
	}

	_getch();
}