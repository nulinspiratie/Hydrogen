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

double walker[N][2][3];

void placewalkers();
int trialmove(int i);
void savedata();
double norm(double particles[]);

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
	double dist[3];
	double dot=0;
	for (int i=0;i<N;i++)
	{
		dot=0;
		for (int j=0;j<3;j++)
		{
			dist[j] = walker[i][0][j] - walker[i][1][j];
			dot += (walker[i][0][j]/norm(walker[i][0]) - walker[i][1][j]/norm(walker[i][1])) * dist[j];
		}
		//cout << "norm:" << norm(dist) << "\tdot=" << dot << endl;
		//cin.ignore();
		double dnorm=norm(dist);
		double mfac = 1 + alpha * dnorm;
		energy += -4 + dot * 1. / (dnorm*mfac*mfac) - 1./(dnorm*mfac*mfac*mfac) - 1./(4.*mfac*mfac*mfac*mfac) + 1./dnorm;
	}
	energy/=N;
	ofstream file("data.dat",fstream::app);
	file << energy << " " << dmax << " " << acceptance << endl;
	file.close();
}

int trialmove(int i)
{
	double newpos[2][3], distold[3], distnew[3];
	for (int j=0;j<3;j++)
	{
		newpos[0][j] = walker[i][0][j] + (drand48() - 0.5) * dmax;
		newpos[1][j] = walker[i][1][j] + (drand48() - 0.5) * dmax;
	}
	for (int j=0; j<3; j++)
	{
		distold[j] = walker[i][0][j] - walker[i][1][j];
		distnew[j] = newpos[0][j] - newpos[1][j];
	}

	double psinew = exp(-2 * norm(newpos[0])) * exp(-2 * norm(newpos[1])) * exp( norm(distnew) / (2*(1 + alpha * norm(distnew))));

	double psiold=exp(-2 * norm(walker[i][0])) * exp(-2 * norm(walker[i][1])) * exp( norm(distold) / (2*(1 + alpha * norm(distold))));
	
	double prob = (psinew * psinew) / (psiold * psiold);
	if (prob > 1 || prob > drand48())
	{
		for (int j=0;j<3;j++)
		{
			walker[i][0][j] = newpos[0][j];
			walker[i][1][j] = newpos[1][j];
		}
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
			for (int k = 0; k < 3; k++)
			{
				walker[i][k][j] = drand48() - 0.5;
			}
		}
	}
}

double norm(double particle[])
{
	return sqrt(particle[0]*particle[0] + particle[1]*particle[1] + particle[2]*particle[2]);
}	
