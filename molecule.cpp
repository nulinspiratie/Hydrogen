#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

#define N 1000
#define thermloops 30000
#define measureloops 250000
#define datainterval 25
#define acceptanceinterval 50

#define minacceptance 0.2
#define maxacceptance 0.4
#define acceptancescale 1.1

double dmax=0.01; 	//Maximum displacement

double alpha = 2;	//Parameter alpha
double a=1;	//Parameter a
double beta;	//Parameter beta
double s;	//Parameter s
double acceptance=0;

double walker[N][2][3];
double wavefunction[N];

char filename[100];

double calcwavefunction(double pos[2][3]);
void placewalkers();
int trialmove(int i);
void savedata();
void determineparam();	//Determines parameter a
double norm(double particles[], double offset=0);

using namespace std;

int main(int argc, char * argv[])
{
	srand48((long)time(NULL));
	if (argc!=3)
	{
		cout << "There must be two input parameters.\nUsage: ./Molecule beta s\n";
		return 1;
	}
	beta = atof(argv[1]);
	s = atof(argv[2]);
	cout << "beta = " << beta << "\t s = " << s << endl;
	if (s < 0.1) return 1;

	//Create initial data storage file
	sprintf(filename,"data/data-%g-%g.dat",beta,s);
	ofstream file(filename);
	file.close();
	
	determineparam();
	cout << "Parameter a: " << a << endl;

	placewalkers();
	
	cout << "Initial walkers placed and associated wavefunctions calculated\n";

	for (int i=0; i<thermloops;i++)
	{
		acceptance=0;
		for (int j=0;j<N;j++)
		{
			acceptance += trialmove(j);
		}

		//Potentially rescale maximum displacement dmax
		if (i%acceptanceinterval == 0)
		{
			if (acceptance < minacceptance * N) dmax /= acceptancescale;
			else if (acceptance > maxacceptance * N) dmax *= acceptancescale;
		}
	}
	cout << "Finished Thermalization\nStarting measurement\n";
	
	//for (int i = 0; i < 10; i++)
	//	printf("(%g, %g, %g)\n(%g, %g, %g)\n\n",walker[i][0][0],walker[i][0][1],walker[i][0][2],walker[i][1][0],walker[i][1][1],walker[i][1][2]);
	//cin.ignore();

	for (int i=0; i<measureloops;i++)
	{
		acceptance=0;
		for (int j=0;j<N;j++)
		{
			acceptance += trialmove(j);
		}
			if (i%datainterval==0) savedata();
		
		//Potentially rescale maximum displacement dmax
		if (i%acceptanceinterval == 0)
		{
			if (acceptance < minacceptance * N) dmax /= acceptancescale;
			else if (acceptance > maxacceptance * N) dmax *= acceptancescale;
		}

	}
	cout << "Finished simulation\n";

	return 0;
}

void determineparam()
{
	for (int i=0;i<50;i++)
		a = 1./(1+exp(-s/a));
}

void savedata()
{
	
	double energy=0;
	double dist[2][2];	//distance between electron and nucleus
	double eldist[3];	//electron-electron distance vector
	double eldistance;	//Electron-electron distance
	double phi[2][2];
	
	for (int i=0; i<N; i++)
	{
		dist[0][0] = norm(walker[i][0], -s/2);
		dist[0][1] = norm(walker[i][0], s/2);
		dist[1][0] = norm(walker[i][1], -s/2);
		dist[1][1] = norm(walker[i][1], s/2);
		for (int j = 0; j < 3; j++)
			eldist[j] = walker[i][0][j] - walker[i][1][j];
		eldistance = norm(eldist);
	
		phi[0][0] = exp(-dist[0][0] /a );
		phi[0][1] = exp(-dist[0][1] /a );
		phi[1][0] = exp(-dist[1][0] /a );
		phi[1][1] = exp(-dist[1][1] /a );

		double en[4]; //four components contributing to the energy
		en[0] = 1./(a * (phi[0][0]+phi[0][1])) * (phi[0][0]/dist[0][0] + phi[0][1]/dist[0][1]);
		en[1] = 1./(a * (phi[1][0]+phi[1][1])) * (phi[1][0]/dist[1][0] + phi[1][1]/dist[1][1]);
		en[2] = -(1./dist[0][0] + 1./dist[0][1] + 1./dist[1][0] + 1./dist[1][1]) + 1./eldistance;
		en[3] = -((4*beta + 1)*eldistance + 4)/(4*pow(1+beta*eldistance,4)*eldistance);
		
		//Dot product component of energy
		double dot = 0;
			dot += ((phi[0][0] * (walker[i][0][0] + s/2)/dist[0][0] + phi[0][1] * (walker[i][0][0] - s/2) / dist[0][1])/(phi[0][0] + phi[0][1]) - (phi[1][0] * (walker[i][1][0] + s/2)/dist[1][0] + phi[1][1] * (walker[i][1][0] - s/2) / dist[1][1])/(phi[1][0] + phi[1][1])) * eldist[0];
			dot += ((phi[0][0] * walker[i][0][1]/dist[0][0] + phi[0][1] * walker[i][0][1] / dist[0][1])/(phi[0][0] + phi[0][1]) - (phi[1][0] * walker[i][1][1] /dist[1][0] + phi[1][1] * walker[i][1][1] / dist[1][1])/(phi[1][0] + phi[1][1])) * eldist[1];
			dot += ((phi[0][0] * walker[i][0][2]/dist[0][0] + phi[0][1] * walker[i][0][2] / dist[0][1])/(phi[0][0] + phi[0][1]) - (phi[1][0] * walker[i][1][2] /dist[1][0] + phi[1][1] * walker[i][1][2] / dist[1][1])/(phi[1][0] + phi[1][1])) * eldist[2];
		dot /= eldistance * 2 * a * (1 + beta * eldistance) * (1 + beta * eldistance);

		energy += en[0] + en[1] + en[2] + en[3] + dot;
	}
	energy/=N;
	energy -= 1./(a*a);

	ofstream file(filename,fstream::app);
	file << energy << " " << dmax << " " << acceptance << endl;
	file.close();
}

double calcwavefunction(double pos[2][3])
{
	double dist[3];
	for (int j=0; j<3; j++)
		dist[j] = pos[0][j] - pos[1][j];
	double distance = norm(dist);
	double phi[2];
	phi[0] = exp(-norm(pos[0], -s/2) /a ) + exp(-norm(pos[0], s/2) /a );
	phi[1] = exp(-norm(pos[1], -s/2) /a ) + exp(-norm(pos[1], s/2) /a );
	double f = exp(norm(dist)/(alpha * (1 + beta * norm(dist))));
	return phi[0] * phi[1] * f;
}

int trialmove(int i)
{
	double newpos[2][3];
	for (int j=0;j<3;j++)
	{
		newpos[0][j] = walker[i][0][j] + (drand48() - 0.5) * 2 * dmax;
		newpos[1][j] = walker[i][1][j] + (drand48() - 0.5) * 2 * dmax;
	}
	double psi = calcwavefunction(newpos);
	double prob = (psi * psi) / (wavefunction[i] * wavefunction[i]);
	if (prob > 1 || prob > drand48())
	{
		for (int j=0;j<3;j++)
		{
			walker[i][0][j] = newpos[0][j];
			walker[i][1][j] = newpos[1][j];
		}
		wavefunction[i] = psi;
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
			for (int k = 0; k < 2; k++)
			{
				walker[i][k][j] = drand48() - 0.5;
			}
		}
		wavefunction[i]=calcwavefunction(walker[i]);
	}
}

double norm(double particle[], double offset)
{
	return sqrt((particle[0]-offset)*(particle[0]-offset) + particle[1]*particle[1] + particle[2]*particle[2]);
}	
