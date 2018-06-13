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
	
	double size;
};

struct cubelete {
	double x;
	double y;
	double z;
	
	double dist;
};

int determineBin(bead b, double gridSize, int numGridX, int numGridY, int numGridZ)
{
	int cellx = (int)((b.x) / gridSize);
	int celly = (int)((b.y) / gridSize);
	int cellz = (int)((b.z) / gridSize);
	
	// temporary fix for if bead is at 0.50
	if(cellx == numGridX)
	{
		cellx = 0;
	}
	if(celly == numGridY)
	{
		celly = 0;
	}
	if(cellz == numGridZ)
	{
		cellz = 0;
	}
	
	int bin = (cellx*numGridY*numGridX) + (celly)*numGridY + cellz;
	
	return bin;
}

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
		while(cubeleteList[i].dist > p)
		{
			i = i + 1;
		}
		
		//while the values coming from the right are greater than the furthest right value (the pivot point) then continue
		while( (cubeleteList[j].dist < p) && (j != left))
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
				
				int tempIndex = 0;
				while(tempString.compare(beadSizes[tempIndex].id) != 0)
				{
					tempIndex++;
					if(tempIndex >= beadSizes.size())
					{
						cout << "Couldn't find " << tempString << " type in " << atomfilename << endl;
						return 1;
					}
				}
				
				systemBeads[systemBeads.size() - 1].size = beadSizes[tempIndex].size;
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
	
	std::vector<bead*> *grid = new std::vector<bead*>[numGridX*numGridY*numGridZ];
	
	// put system beads into a grid
	for(int i = 0; i < systemBeads.size(); i++)
	{		
		int bin = determineBin(systemBeads[i], gridSize, numGridX, numGridY, numGridZ);
		
		//cout << bin << " " << numGridX*numGridY*numGridZ << endl;
		
		(grid[bin]).push_back(&(systemBeads[i]));
	}
	
	double cubeleteSize = 0.2;
	
	int numCubeX = (int)(boxlx / cubeleteSize);
	int numCubeY = (int)(boxly / cubeleteSize);
	int numCubeZ = (int)(boxlz / cubeleteSize);
	
	//stores the distances to closest system bead from position associated with [i*numCubeX^2 + j*numCubeY + k]
	int cubeleteListSize = numCubeX*numCubeY*numCubeZ;
	struct cubelete *cubeleteList = new struct cubelete[cubeleteListSize];
	
	int maxHistDist = 40.0;
	double histStep = 0.1;
	int histSize = (int)(maxHistDist / histStep);
	int *hist = new int[histSize];
	for(int h = 0; h < histSize; h++)
	{
		hist[h] = 0;
	}
	
	// omp parallel statement here
	#pragma omp parallel for
	for(int i = 0; i < numCubeX; i++)
	{
		cout << i << " " << numCubeX << endl;
		double tempX = (i + 0.5) * cubeleteSize;
		for(int j = 0; j < numCubeY; j++)
		{
			double tempY = (j + 0.5) * cubeleteSize;
			for(int k = 0; k < numCubeZ; k++)
			{
				double tempZ = (k + 0.5) * cubeleteSize;
				
				int index = i*numCubeX*numCubeX + j*numCubeY + k;
				// find distance from [i*numCubeX^2 + j*numCubeY + k] to nearest system bead
				cubeleteList[index].x = tempX;
				cubeleteList[index].y = tempY;
				cubeleteList[index].z = tempZ;
					
				//for(int b = -)
				for(int b = 0; b < systemBeads.size(); b++)
				//for(int b = 0; b < 1; b++)
				{
					// need to periodically wrap these coordinates?
					// can speed this up by griding system originally?
					double xDist = tempX - systemBeads[b].x;
					xDist = periodicWrap(xDist, boxlx);
					double yDist = tempY - systemBeads[b].y;
					yDist = periodicWrap(yDist, boxly);
					double zDist = tempZ - systemBeads[b].z;
					zDist = periodicWrap(zDist, boxlz);
					
					double dist = xDist*xDist + yDist*yDist + zDist*zDist;
					dist = sqrt(dist);
					
					dist = dist - 0.5*systemBeads[b].size;
					
					//dist = dist * 2.0;
					
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
	quicksort(0, cubeleteListSize, cubeleteList); 
	
	cout << cubeleteListSize << endl;
	/*for(int i = 0; i < numCubeX*numCubeY*numCubeZ; i++)
	{
		cout << i << " " << cubeleteList[i].dist << endl;
		cout << cubeleteList[i].x << " " << cubeleteList[i].y << " " << cubeleteList[i].z << endl;
		exit(0);
	}*/
	
	
	// randomly pick cubelete
	int numTrials = 1000;
	
	#pragma omp parallel for
	for(int i = 0; i < numTrials; i++)
	{
		double tempRand = rand() / 32768.0;
		cout << "i: " << i << " " << tempRand << endl;
		int randomCubeIndex = (int)(cubeleteListSize * tempRand);
		
		//cout << "after rand()" << endl;
		//cout << randomCubeIndex << endl;
		
		double tempX = cubeleteList[randomCubeIndex].x;
		double tempY = cubeleteList[randomCubeIndex].y;
		double tempZ = cubeleteList[randomCubeIndex].z;
		
		//cout << "before loop" << endl;
		
		for(int j = 0; j < cubeleteListSize; j++)
		{
			//cout << "j: " << j << endl;
			double distX = cubeleteList[j].x - tempX;
			distX = periodicWrap(distX, boxlx);
			double distY = cubeleteList[j].y - tempY;
			distY = periodicWrap(distY, boxly);
			double distZ = cubeleteList[j].z - tempZ;
			distZ = periodicWrap(distZ, boxlz);
			
			distX = distX*distX;
			distY = distY*distY;
			distZ = distZ*distZ;
			
			double dist = distX + distY + distZ;
			dist = sqrt(dist);
			if(dist <= cubeleteList[j].dist)
			{
				// this is the largest sphere which randomCubeIndex can fit inside
				//cout << "j: " << j << endl;
				
				int tempHistPos = (int)(cubeleteList[j].dist * 2.0 / histStep);
				
				for(int h = 0; h < tempHistPos; h++)
				{
					hist[h] += 1;
				}
				break;
			}
				
			
		}
	}
	
	for (int h = 0; h < histSize; h++)
	{
		cout << h*histStep + 0.5*histStep << " " << hist[h] / hist[0] << endl;
	}
	cout << "end" << endl;
	return 0;
}