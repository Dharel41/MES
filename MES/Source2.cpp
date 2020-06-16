#include<conio.h>
#include<iostream>
#include"Header.h"
#include<cmath>
using namespace std;

ostream &operator<<(ostream &wyjscie, const point &a)
{

	wyjscie << "Wezel nr." << a.id << ".Wspolrzedne wezla (" << a.x << "," << a.y << ")" << endl;
	return wyjscie;
}

ostream &operator<<(ostream &wyjscie, const node &a)
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

grid mesh(int nheight, int nlength, double height, double length) {

	//nheight, nlength -ilosc wezlow po wysokosci i dlugosci
	int count = 0;
	double delta_y = height / (nheight - 1);
	double delta_x = length / (nlength - 1);

	node *n1 = new node[nheight*nlength];
	element *e1 = new element[(nheight - 1)*(nlength - 1)];

	for (int i = 0;i<nlength;i++)
	{
		for (int j = 0;j<nheight;j++)
		{
			n1[count].id = count;
			n1[count].x = 0 + i*delta_x;
			n1[count].y = 0 + j*delta_y;

			if (i == 0 || j == 0 || i == nlength - 1 || j == nheight - 1)
			{
				n1[count].status = 1;
			}
			count++;
		}
	}

	count = 0;
	for (int i = 0;i<nlength - 1;i++)
	{
		for (int j = 0;j<nheight - 1;j++)
		{
			e1[count].id = count;

			e1[count].points[0] = n1[count + i];
			e1[count].points[1] = n1[count + i + nheight];
			e1[count].points[2] = n1[count + i + 1 + nheight];
			e1[count].points[3] = n1[count + i+1];
			count++;
		}
	}
	return  grid(n1, e1, nheight, nlength);
}

double*** f_matrix_H(int nheight, int nlength, double height, double length,grid a)
{
	const int conductivity = 25;
	double ksi[4] = { -1 / sqrt(3),1 / sqrt(3),1 / sqrt(3),-1 / sqrt(3) };
	double eta[4] =  { -1 / sqrt(3),-1 / sqrt(3),1 / sqrt(3),1 / sqrt(3) };
	double shape_function[4][4];   //pierwsza liczba to N1,N2,N3,N4 , druga to punkty calkowania
	for (int i = 0;i<4;i++)
	{
		shape_function[0][i] = 0.25 * (1 - ksi[i])*(1 - eta[i]);
		shape_function[1][i] = 0.25 * (1 + ksi[i])*(1 - eta[i]);
		shape_function[2][i] = 0.25 * (1 + ksi[i])*(1 + eta[i]);
		shape_function[3][i] = 0.25 * (1 - ksi[i])*(1 + eta[i]);
	}


	for(int i=0;i<((nheight-1)*(nlength-1));i++)
	{
		for (int j = 0;j < 4;j++)
		{
			a.e1[i].new_coordinates[j].id = a.e1[i].points[j].id;
			for (int k = 0;k < 4;k++)
			{
				a.e1[i].new_coordinates[j].x += (shape_function[k][j] * a.e1[i].points[k].x);
				a.e1[i].new_coordinates[j].y += (shape_function[k][j] * a.e1[i].points[k].y);
			}
		}
	}


	double *** matrix_ksi = new double **[(nheight - 1)*(nlength - 1)];
	double *** matrix_eta = new double **[(nheight - 1)*(nlength - 1)];    //tablica 3d /[nr.elem][nr fun ksztaltu][nr pkt calkowania] (pochodne funkcji ksztaltu wzgledem ksi i eta)
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		matrix_ksi[i] = new double*[4];
		matrix_eta[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			matrix_ksi[i][j] = new double[4];
			matrix_eta[i][j] = new double[4];
		}
	}
	for (int i = 0;i<((nheight - 1)*(nlength - 1));i++)
	{
	
			for(int k=0;k<4;k++)
			{
				
				matrix_ksi[i][0][k] = -0.25*(1 - eta[k]);
				matrix_ksi[i][1][k] = 0.25*(1 - eta[k]);
				matrix_ksi[i][2][k] = 0.25*(1 + eta[k]);
				matrix_ksi[i][3][k] = -0.25*(1 +eta[k]);

				matrix_eta[i][0][k] = -0.25*(1 - ksi[k]);
				matrix_eta[i][1][k] = -0.25*(1 + ksi[k]);
				matrix_eta[i][2][k] = 0.25*(1 + ksi[k]);
				matrix_eta[i][3][k] = 0.25*(1 - ksi[k]);
			}
		
		
	}

	double *** J = new double **[(nheight - 1)*(nlength - 1)];    //tablica 3d /[nr.elem][nr fun ksztaltu][nr pkt calkowania]
	double *** J2 = new double **[(nheight - 1)*(nlength - 1)];   //J - Macierz Jakobiego, J2 - Macierz Jakobiego^(-1)
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		J[i] = new double*[4];
		J2[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			J[i][j] = new double[4];
			J2[i][j] = new double[4];
		}
	}

	for (int i = 0;i < ((nheight - 1)*(nlength - 1));i++)
	{
		for (int j = 0;j < 4;j++)
		{

			for (int k = 0;k < 4;k++)
			{
				J[i][j][k] = 0;
			}
		}
	}

	for (int i = 0;i < ((nheight - 1)*(nlength - 1));i++)
	{
			for (int j = 0;j < 4;j++)
			{
				for (int k = 0;k < 4;k++)
				{

					J[i][0][j] += (matrix_ksi[i][k][j] * a.e1[i].points[k].x);
					J[i][1][j] += (matrix_ksi[i][k][j] * a.e1[i].points[k].y);
					J[i][2][j] += (matrix_eta[i][k][j] * a.e1[i].points[k].x);
					J[i][3][j] += (matrix_eta[i][k][j] * a.e1[i].points[k].y);
				}
			}
	}

	double **detJ = new double*[(nheight - 1)*(nlength - 1)];   //wyznacznik Macierzy Jakobiego - Jakobian  tablica 2d[nr.elementu][pkt.calkowania]
	for (int i = 0;i < (nheight - 1)*(nlength - 1);i++)
		detJ[i] = new double[4];


	for(int i=0;i<((nheight - 1)*(nlength - 1));i++)
	{
		for(int j=0;j<4;j++)
		{
			detJ[i][j] = (J[i][0][j] * J[i][3][j]) - (J[i][1][j] * J[i][2][j]);   
		}
	}

	
	for (int i = 0;i < (nheight - 1)*(nlength - 1);i++)
	{
		for (int j = 0;j < 4;j++)
		{
					J2[i][0][j] = (J[i][3][j] / detJ[i][j]);
					J2[i][1][j] = (-J[i][1][j] / detJ[i][j]);
					J2[i][2][j] = (-J[i][2][j] / detJ[i][j]);
					J2[i][3][j] = (J[i][0][j] / detJ[i][j]);
		}
	}

	double *** matrix_dx = new double **[(nheight - 1)*(nlength - 1)];        //tablica 3d /[nr.elem][nr pkt calkowania][nr fun ksztaltu]
	double *** matrix_dy = new double **[(nheight - 1)*(nlength - 1)];        //pochodne funkcji ksztaltu wzglêdem x i y
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		matrix_dx[i] = new double*[4];
		matrix_dy[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			matrix_dx[i][j] = new double[4];
			matrix_dy[i][j] = new double[4];
		}
	}

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for(int j=0;j<4;j++)
		{
			for(int k=0;k<4;k++)
			{
				int w = 0;
				int w1 = 2;
				matrix_dx[i][k][j] = J2[i][w][k] * matrix_ksi[i][j][k] + J2[i][w + 1][k] * matrix_eta[i][j][k];
				matrix_dy[i][k][j] = J2[i][w1][k] * matrix_ksi[i][j][k] + J2[i][w1 + 1][k] * matrix_eta[i][j][k];
			}
		}

	}


	double **** matrix_dxdx = new double ***[(nheight - 1)*(nlength - 1)];   //tablica 4d /[nr.elem][nr pkt calkowania][kolumna matrix_dx][wiersz matrix_dx] 
	double **** matrix_dydy = new double ***[(nheight - 1)*(nlength - 1)];   // iloczyny {dN/dx}{dN/dx}T    {dN/dy}{dN/dy}T 
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		matrix_dxdx[i] = new double**[4];
		matrix_dydy[i] = new double**[4];
		for (int j = 0;j<4;j++)
		{
			matrix_dxdx[i][j] = new double*[4];
			matrix_dydy[i][j] = new double*[4];
			for(int k=0;k<4;k++)
			{
				matrix_dxdx[i][j][k] = new double[4];
				matrix_dydy[i][j][k] = new double[4];
			}
		}
	}

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{
				
				for(int l=0;l<4;l++)
				{
					matrix_dxdx[i][j][k][l] = matrix_dx[i][j][l] * matrix_dx[i][j][k];
					matrix_dydy[i][j][k][l] = matrix_dy[i][j][l] * matrix_dy[i][j][k];
				}
			}

		}

	}


	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)    //mnozenie przez wyznacznik
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{

				for (int l = 0;l<4;l++)
				{
					matrix_dxdx[i][j][k][l] *= detJ[i][j];
					matrix_dydy[i][j][k][l] *= detJ[i][j];
				}
			}

		}

	}


	double **** matrix_k_sum = new double ***[(nheight - 1)*(nlength - 1)]; //sumowanie macierzy * przewodnosc
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		matrix_k_sum[i] = new double**[4];
		for (int j = 0;j<4;j++)
		{
			matrix_k_sum[i][j] = new double*[4];
			for (int k = 0;k<4;k++)
			{
				matrix_k_sum[i][j][k] = new double[4];
			}
		}
	}


	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{

				for (int l = 0;l<4;l++)
				{
					matrix_k_sum[i][j][k][l] = conductivity*(matrix_dxdx[i][j][k][l] + matrix_dydy[i][j][k][l]);
				}
			}

		}

	}

	double *** matrix_H = new double **[(nheight - 1)*(nlength - 1)]; //macierz H
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		matrix_H[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			matrix_H[i][j] = new double[4];
		}
	}

	for(int i=0;i<(nheight - 1)*(nlength - 1);i++)
	{
		for(int j=0;j<4;j++)
		{
			for(int k=0;k<4;k++)
			{
				matrix_H[i][j][k] = 0;
				for (int l = 0;l < 4;l++)
				{
					matrix_H[i][j][k] += matrix_k_sum[i][l][j][k];
				}
			}
		}
	}


		return matrix_H;
}


	
double*** f_matrix_H_BC(int nheight, int nlength, double height, double length,grid a)
{

	//////// Matrix H boundary condition
	const int convection = 300;
	double **** pc1 = new double ***[(nheight - 1)*(nlength - 1)];
	double **** pc2 = new double ***[(nheight - 1)*(nlength - 1)];
	double **** pc_sum = new double ***[(nheight - 1)*(nlength - 1)];
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		pc1[i] = new double**[4];
		pc2[i] = new double**[4];
		pc_sum[i] = new double**[4];
		for (int j = 0;j<4;j++)
		{
			pc1[i][j] = new double*[4];
			pc2[i][j] = new double*[4];
			pc_sum[i][j] = new double*[4];
			for (int k = 0;k<4;k++)
			{
				pc1[i][j][k] = new double[4];
				pc2[i][j][k] = new double[4];
				pc_sum[i][j][k] = new double[4];
			}
		}
	}
	double ksi2[8] = { -1 / sqrt(3),1,1 / sqrt(3) ,-1,1 / sqrt(3) ,1,-1 / sqrt(3),-1 };
	double eta2[8] = { -1,-1 / sqrt(3),1,1 / sqrt(3),-1,1 / sqrt(3),1,-1 / sqrt(3) };
	double shape_function2[4][8];  //[N1...N4][ksi-eta]
	for (int i = 0;i<8;i++)
	{
		shape_function2[0][i] = 0.25 * (1 - ksi2[i])*(1 - eta2[i]);
		shape_function2[1][i] = 0.25 * (1 + ksi2[i])*(1 - eta2[i]);
		shape_function2[2][i] = 0.25 * (1 + ksi2[i])*(1 + eta2[i]);
		shape_function2[3][i] = 0.25 * (1 - ksi2[i])*(1 + eta2[i]);
	}

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{
				for (int l = 0;l<4;l++)
				{
					pc1[i][j][k][l] = shape_function2[l][j] * shape_function2[k][j] * convection;
					pc2[i][j][k][l] = shape_function2[l][4 + j] * shape_function2[k][j + 4] * convection;
				}
			}
		}
	}
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{
				for (int l = 0;l<4;l++)
				{
					if (j == 0 || j == 2)
					{
						pc_sum[i][j][k][l] = (pc1[i][j][k][l] + pc2[i][j][k][l])*(height/ (nheight-1) / 2);
					}
					else
					{
						pc_sum[i][j][k][l] = (pc1[i][j][k][l] + pc2[i][j][k][l])*(length / (nlength - 1) /2);
					}

				}
			}
		}
	}



	double *** final_matrix_H = new double **[(nheight - 1)*(nlength - 1)]; //macierz H
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		final_matrix_H[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			final_matrix_H[i][j] = new double[4];
		}
	}
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{
				final_matrix_H[i][j][k] = 0;
			}
		}
	}

	double ** pow = new double *[(nheight - 1)*(nlength - 1)];
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		pow[i] = new double[4];
	}

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		if (a.e1[i].points[0].status == 1 && a.e1[i].points[1].status == 1)
			pow[i][0] = 1;
		else
			pow[i][0] = 0;

		if (a.e1[i].points[1].status == 1 && a.e1[i].points[2].status == 1)
			pow[i][1] = 1;
		else
			pow[i][1] = 0;

		if (a.e1[i].points[2].status == 1 && a.e1[i].points[3].status == 1)
			pow[i][2] = 1;
		else
			pow[i][2] = 0;

		if (a.e1[i].points[3].status == 1 && a.e1[i].points[0].status == 1)
			pow[i][3] = 1;
		else
			pow[i][3] = 0;
	}

	for (int i = 0;i < (nheight - 1)*(nlength - 1);i++)
	{
		for (int j = 0;j < 4;j++)
		{

			for (int k = 0;k < 4;k++)
			{
				final_matrix_H[i][j][k] = pow[i][0] * pc_sum[i][0][j][k] + pow[i][1] * pc_sum[i][1][j][k] + pow[i][2] * pc_sum[i][2][j][k] + pow[i][3] * pc_sum[i][3][j][k];
			}

		}
	}
	return final_matrix_H;
}








double*** f_matrix_C(int nheight, int nlength, double height, double length,grid a)
{
	const int c = 700;
	const int ro = 7800;
	double ksi[4] = { -1 / sqrt(3),1 / sqrt(3),1 / sqrt(3),-1 / sqrt(3) };
	double eta[4] = { -1 / sqrt(3),-1 / sqrt(3),1 / sqrt(3),1 / sqrt(3) };
	double shape_function[4][4];   //pierwsza liczba to N1,N2,N3,N4 , druga to punkty calkowania
	for (int i = 0;i<4;i++)
	{
		shape_function[0][i] = 0.25 * (1 - ksi[i])*(1 - eta[i]);
		shape_function[1][i] = 0.25 * (1 + ksi[i])*(1 - eta[i]);
		shape_function[2][i] = 0.25 * (1 + ksi[i])*(1 + eta[i]);
		shape_function[3][i] = 0.25 * (1 - ksi[i])*(1 + eta[i]);
	}


	for (int i = 0;i<((nheight - 1)*(nlength - 1));i++)
	{
		for (int j = 0;j < 4;j++)
		{
			a.e1[i].new_coordinates[j].id = a.e1[i].points[j].id;
			for (int k = 0;k < 4;k++)
			{
				a.e1[i].new_coordinates[j].x += (shape_function[k][j] * a.e1[i].points[k].x);
				a.e1[i].new_coordinates[j].y += (shape_function[k][j] * a.e1[i].points[k].y);
			}
		}
	}


	double *** matrix_ksi = new double **[(nheight - 1)*(nlength - 1)];
	double *** matrix_eta = new double **[(nheight - 1)*(nlength - 1)];    //tablica 3d /[nr.elem][nr fun ksztaltu][nr pkt calkowania] (pochodne funkcji ksztaltu wzgledem ksi i eta)
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		matrix_ksi[i] = new double*[4];
		matrix_eta[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			matrix_ksi[i][j] = new double[4];
			matrix_eta[i][j] = new double[4];
		}
	}
	for (int i = 0;i<((nheight - 1)*(nlength - 1));i++)
	{

		for (int k = 0;k<4;k++)
		{

			matrix_ksi[i][0][k] = -0.25*(1 - eta[k]);
			matrix_ksi[i][1][k] = 0.25*(1 - eta[k]);
			matrix_ksi[i][2][k] = 0.25*(1 + eta[k]);
			matrix_ksi[i][3][k] = -0.25*(1 + eta[k]);

			matrix_eta[i][0][k] = -0.25*(1 - ksi[k]);
			matrix_eta[i][1][k] = -0.25*(1 + ksi[k]);
			matrix_eta[i][2][k] = 0.25*(1 + ksi[k]);
			matrix_eta[i][3][k] = 0.25*(1 - ksi[k]);
		}


	}

	double *** J = new double **[(nheight - 1)*(nlength - 1)];    //tablica 3d /[nr.elem][nr fun ksztaltu][nr pkt calkowania] Macierz Jakobiego
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		J[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			J[i][j] = new double[4];

		}
	}

	for (int i = 0;i < ((nheight - 1)*(nlength - 1));i++)
	{
		for (int j = 0;j < 4;j++)
		{

			for (int k = 0;k < 4;k++)
			{
				J[i][j][k] = 0;
			}
		}
	}

	for (int i = 0;i < ((nheight - 1)*(nlength - 1));i++)
	{
		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{

				J[i][0][j] += (matrix_ksi[i][k][j] * a.e1[i].points[k].x);
				J[i][1][j] += (matrix_ksi[i][k][j] * a.e1[i].points[k].y);
				J[i][2][j] += (matrix_eta[i][k][j] * a.e1[i].points[k].x);
				J[i][3][j] += (matrix_eta[i][k][j] * a.e1[i].points[k].y);
			}
		}
	}

	double **detJ = new double*[(nheight - 1)*(nlength - 1)]; //Jakobian
	for (int i = 0;i < (nheight - 1)*(nlength - 1);i++)
		detJ[i] = new double[4];


	for (int i = 0;i<((nheight - 1)*(nlength - 1));i++)
	{
		for (int j = 0;j<4;j++)
		{
			detJ[i][j] = (J[i][0][j] * J[i][3][j]) - (J[i][1][j] * J[i][2][j]);  
		}
	}

	////////////////////////////////////Matrix C////////////////////////////////////////
	double **** ctab = new double***[(nheight - 1)*(nlength - 1)];   //{N}{N}*detJ*c*ro
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		ctab[i] = new double**[4];
		for (int j = 0;j<4;j++)
		{
			ctab[i][j] = new double*[4];
			for (int k = 0;k<4;k++)
			{
				ctab[i][j][k] = new double[4];
			}
		}
	}

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{
				for (int l = 0;l<4;l++)
				{
					ctab[i][j][k][l] = shape_function[l][j] * shape_function[k][j] * detJ[i][j] * c*ro;
				}
			}
		}
	}

	double *** matrix_C = new double **[(nheight - 1)*(nlength - 1)];     //Macierz C
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		matrix_C[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			matrix_C[i][j] = new double[4];
		}
	}



	for (int i = 0;i<(nheight - 1)*(nlength - 1);i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{
				matrix_C[i][j][k] = 0;
				for (int l = 0;l < 4;l++)
				{
					matrix_C[i][j][k] += ctab[i][l][j][k];
				}
			}
		}
	}
	
	return matrix_C;
}


















double** f_vector_P(int nheight, int nlength, double height, double length, grid a)
{

	const int convection = 300;
	const int ambient_t = 1200;
	double *** pc_sum = new double **[(nheight - 1)*(nlength - 1)];
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		pc_sum[i] = new double*[4];
		for (int j = 0;j<4;j++)
		{
			pc_sum[i][j] = new double[4];
		}
	}
	double ksi2[8] = { -1 / sqrt(3),1,1 / sqrt(3) ,-1,1 / sqrt(3) ,1,-1 / sqrt(3),-1 };
	double eta2[8] = { -1,-1 / sqrt(3),1,1 / sqrt(3),-1,1 / sqrt(3),1,-1 / sqrt(3) };
	double shape_function2[4][8];  //[N1...N4][ksi-eta]
	for (int i = 0;i<8;i++)
	{
		shape_function2[0][i] = 0.25 * (1 - ksi2[i])*(1 - eta2[i]);
		shape_function2[1][i] = 0.25 * (1 + ksi2[i])*(1 - eta2[i]);
		shape_function2[2][i] = 0.25 * (1 + ksi2[i])*(1 + eta2[i]);
		shape_function2[3][i] = 0.25 * (1 - ksi2[i])*(1 + eta2[i]);
	}

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			for (int k = 0;k<4;k++)
			{
					pc_sum[i][j][k] = shape_function2[k][j] * convection*ambient_t + shape_function2[k][4 + j] * convection*ambient_t;
					pc_sum[i][j][k] = pc_sum[i][j][k] * (height / (nheight - 1) / 2);
			}
		}
	}
	
	double ** final_vector_p = new double *[(nheight - 1)*(nlength - 1)]; 
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		final_vector_p[i] = new double[4];
	}
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		for (int j = 0;j<4;j++)
		{
			final_vector_p[i][j] = 0;
		}
	}

	double ** pow = new double *[(nheight - 1)*(nlength - 1)];
	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		pow[i] = new double[4];
	}

	for (int i = 0; i < (nheight - 1)*(nlength - 1); i++)
	{
		if (a.e1[i].points[0].status == 1 && a.e1[i].points[1].status == 1)
			pow[i][0] = 1;
		else
			pow[i][0] = 0;

		if (a.e1[i].points[1].status == 1 && a.e1[i].points[2].status == 1)
			pow[i][1] = 1;
		else
			pow[i][1] = 0;

		if (a.e1[i].points[2].status == 1 && a.e1[i].points[3].status == 1)
			pow[i][2] = 1;
		else
			pow[i][2] = 0;

		if (a.e1[i].points[3].status == 1 && a.e1[i].points[0].status == 1)
			pow[i][3] = 1;
		else
			pow[i][3] = 0;
	}

	for (int i = 0;i < (nheight - 1)*(nlength - 1);i++)
	{
		for (int j = 0;j < 4;j++)
		{

				final_vector_p[i][j] = pow[i][0] * pc_sum[i][0][j] + pow[i][1] * pc_sum[i][1][j] + pow[i][2] * pc_sum[i][2][j] + pow[i][3] * pc_sum[i][3][j];
		}
	}
	return final_vector_p;
}
