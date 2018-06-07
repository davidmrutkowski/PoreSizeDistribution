//============================================================================
// Name        : psd.cpp
// Author      : David Rutkowski
// Description : Calculates the pore size distribution (psd) based originally on code from Poreblazer (Lev Sarkisov group)
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

struct beadType {
	string id;
	double size;
};

struct bead {
	string id;
	double x;
	double y;
	double z;
};

struct cubelete {
	int x;
	int y;
	int z;
	
	double dist;
};


double periodicWrap(double value, double boxLength)
{
	double wrappedValue = value - (boxLength*(double)((int)(value/(double)(boxLength) + (0.5*copysign(1.0,value)))));
	return wrappedValue;
}

void exch(int i, int right, struct cubelete *cubeleteList)
{
	struct cubelete tempCubelete = cubeleteList[i];
	cubeleteList[i] = cubeleteList[right];
	cubeleteList[right] = tempCubelete;
}

int partition(int left, int right, struct cubelete *cubeleteList)
{
	int i = left;
	double p = cubeleteList[right].dist;
	int j = right - 1;
	int outer = 1;

	while (i <= j)
	{
		//while the values coming from the left are less than the furthest right value (the pivot point) then continue		
		while(cubeleteList[i].dist < p)
		{
			i = i + 1;
		}
		
		//while the values coming from the right are greater than the furthest right value (the pivot point) then continue
		while( (cubeleteList[j].dist > p) && (j != left))
		{
			j = j - 1;
			
			// if the j and i meet then break out of this loop (j == i?)
			if( j == left)
			{
				break;
			}
		}
		
		//break out of main loop if i is to the right of j
		if(i >= j)
		{
			break;
		}
		
		//exchange i and j
		exch(i, j, cubeleteList);
		i = i + 1;
		j = j - 1;
	}
	
	//exchange i and right once i is greater than or equal to j
	exch(i, right, cubeleteList);

	return i;
}


int quicksort(int left, int right, struct cubelete *cubeleteList)
{
	if(left >= right)
	{
		return 0;
	}

	int m = partition(left, right, cubeleteList);
	//calls quicksort on left side since the pivot is now in the correct relative position (ones that are lower than it are to the left, higher are to the right)
	int d = quicksort(left, m-1, cubeleteList);
	//calls quicksort on the right side
	d = quicksort(m+1, right, cubeleteList);

	return d;
}


int main() 
{
	cout << "begin" << endl;
	
	string atomfilename;
	string xyzfilename;
	
	// read in input.dat
	ifstream inputfile ("input.dat");
	if (inputfile.is_open())
	{
		if(inputfile >> atomfilename and inputfile >> xyzfilename)
		{
			
		}
		inputfile.close();
	}
	else
	{
		cout << "Could not find DREIDING.atoms" << endl;
		return 0;
	}
	
	// read in a file with bead diameters
	int numBeadTypes = 0;
	std::vector<beadType> beadSizes;
	
	ifstream beadtypefile (atomfilename.c_str());
	if (beadtypefile.is_open())
	{
		if(beadtypefile >> numBeadTypes)
		{
			string tempString;
			double tempSize;
			while(beadtypefile >> tempString >> tempSize)
			{
				if(tempString.compare("!") != 0)
				{
					cout << tempString << endl;
				
					beadSizes.push_back(beadType());
					beadSizes[beadSizes.size() - 1].id = tempString;
					beadSizes[beadSizes.size() - 1].size = tempSize;
				}
			}
		}
		beadtypefile.close();
	}
	else
	{
		cout << "Could not find " << atomfilename << endl;
		return 0;
	}
  
	// read in an .xyz file with system coordinates
	int numBeads;
	std::vector<bead> systemBeads;
	ifstream xyzfile (xyzfilename.c_str());
	if (xyzfile.is_open())
	{
		if(xyzfile >> numBeads)
		{
			string tempString;
			double tempX, tempY, tempZ;
			
			while(xyzfile >> tempString >> tempX >> tempY >> tempZ)
			{				
				systemBeads.push_back(bead());
				systemBeads[systemBeads.size() - 1].id = tempString;
				systemBeads[systemBeads.size() - 1].x = tempX;
				systemBeads[systemBeads.size() - 1].y = tempY;
				systemBeads[systemBeads.size() - 1].z = tempZ;
			}
		}
		xyzfile.close();
	}
	else
	{
		cout << "Could not find " << xyzfilename << endl;
		return 0;
	}
	
	
	double boxlx = 25.83200;
	double boxly = 25.83200;
	double boxlz = 25.83200;
	
	double gridSize = 0.5;
	int numGridX = (int)(boxlx / gridSize);
	int numGridY = (int)(boxly / gridSize);
	int numGridZ = (int)(boxlz / gridSize);
	
	
	double cubeleteListize = 0.01;
	
	int numCubeX = (int)(boxlx / cubeleteListize);
	int numCubeY = (int)(boxly / cubeleteListize);
	int numCubeZ = (int)(boxlz / cubeleteListize);
	
	//stores the distances to closest system bead from position associated with [i*numCubeX^2 + j*numCubeY + k]
	struct cubelete *cubeleteList = new struct cubelete[numCubeX*numCubeY*numCubeZ];
	
	// omp parallel statement here
	for(int i = 0; i < numCubeX; i++)
	{
		cout << i << " " << numCubeX << endl;
		double tempX = i * cubeleteListize;
		for(int j = 0; j < numCubeY; j++)
		{
			double tempY = j * cubeleteListize;
			for(int k = 0; k < numCubeZ; k++)
			{
				double tempZ = k * cubeleteListize;
				
				int index = i*numCubeX*numCubeX + j*numCubeY + k;
				// find distance from [i*numCubeX^2 + j*numCubeY + k] to nearest system bead
				cubeleteList[index].x = i;
				cubeleteList[index].y = j;
				cubeleteList[index].z = k;
						
				//for(int b = 0; b < systemBeads.size(); b++)
				for(int b = 0; b < 1; b++)
				{
					// need to periodically wrap these coordinates?
					// can speed this up by griding system originally?
					double xDist = tempX - systemBeads[b].x;
					double yDist = tempY - systemBeads[b].y;
					double zDist = tempZ - systemBeads[b].z;
					
					double dist = xDist*xDist + yDist*yDist + zDist*zDist;
					dist = sqrt(dist);
					
					if (dist < cubeleteList[index].dist || b == 0)
					{
						cubeleteList[index].dist = dist;
					}
				}
				//cubeleteList[i*numCubeX*numCubeX + j*numCubeY + k] = dist;
			}
		}
	}
	
	// sort cubelets using quicksort
	quicksort(0, numCubeX*numCubeY*numCubeZ, cubeleteList); 
	
	cout << numCubeX*numCubeY*numCubeZ << endl;
	/*for(int i = 0; i < numCubeX*numCubeY*numCubeZ; i++)
	{
		cout << cubeleteList[i].dist << endl;
	}*/
	// randomly pick cubelete
	
	cout << "end" << endl;
	return 0;
}