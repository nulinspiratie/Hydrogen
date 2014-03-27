#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

#define N 400
#define thermloops 4000
#define measureloops 30000
#define datainterval 100
#define acceptanceinterval 50

#define minacceptance 0.2
#define maxacceptance 0.4
#define acceptancescale 1.1

double dmax=0.01; 	//Maximum discplacement

double alpha;
double acceptance=0;

double walker[N][3];

void placewalkers();
int trialmove(int j);
void savedata();

using namespace std;

int main(int argc, char * argv[])
{
	srand48((long)time(NULL));
	if (argc!=2)
	{
		cout << "There must be an input parameter for alpha\n";
		return 1;
	}
	alpha = atof(argv[1]);
	
	//Create initial data storage file
	ofstream file("data.dat");
	file.close();

	placewalkers();
	
	for (int i=0; i<thermloops;i++)
	{
		acceptance=0;
		for (int j=0;j<N;j++)
		{
			acceptance += trialmove(j);
			//if (j<10) cout << walker[j] << " \t";
		}
		//cin.ignore();
		if (i%acceptanceinterval == 0)
		{
			if (acceptance < minacceptance * N) dmax /= acceptancescale;
			else if (acceptance > maxacceptance * N) dmax *= acceptancescale;
		}
	}
	cout << "Finished Thermalization\n";

	for (int i=0; i<measureloops;i++)
	{
		acceptance=0;
		for (int j=0;j<N;j++)
		{
			acceptance += trialmove(j);
		}
			if (i%datainterval==0) savedata();
		if (i%acceptanceinterval == 0)
		{
			if (acceptance < minacceptance * N) dmax /= acceptancescale;
			else if (acceptance > maxacceptance * N) dmax *= acceptancescale;
		}

	}
	cout << "Finished simulation\n";

	return 0;
}

void savedata()
{
	double energy=0;
	for (int i=0;i<N;i++)
		energy += alpha + ( walker[i][0] * walker[i][0] + walker[i][1] * walker[i][1] + walker[i][2] * walker[i][2] )  * (0.5 -2 * alpha * alpha);
		//energy+=exp(-alpha * walker[i] * walker[i]);
	energy/=N;
	ofstream file("data.dat",fstream::app);
	file << energy << " " << dmax << " " << acceptance << endl;
	file.close();
}

int trialmove(int j)
{
	double d[3]; 	//Discplacement
	d[0] = 2 * (drand48() - 0.5) * dmax;
	d[1] = 2 * (drand48() - 0.5) * dmax;
	d[2] = 2 * (drand48() - 0.5) * dmax;

	//double psiold=exp(-alpha * walker[j]*walker[j]);
	//double psinew = exp(- alpha * (walker[j] + displacement) * (walker[j] + displacement);
	double prob = exp(-2 *alpha * ( d[0] * (2 * walker[j][0] + d[0]) + d[1] * (2 * walker[j][1] + d[1]) + d[2] * (2 * walker[j][2] + d[2]) ) );
	if (prob > 1 || prob > drand48())
	{
		walker[j][0] += d[0];
		walker[j][1] += d[1];
		walker[j][2] += d[2];
		return 1;
	}
	else
		return 0;
}

void placewalkers()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			walker[i][j] = drand48() - 0.5;
		}
	}
}
