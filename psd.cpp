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
	int x;
	int y;
	int z;
	
	double dist;
};

int determineBin(double x, double y, double z, double gridSize, int numGridX, int numGridY, int numGridZ)
{
	int cellx = (int)((x) / gridSize);
	int celly = (int)((y) / gridSize);
	int cellz = (int)((z) / gridSize);
	
	// temporary fix for if bead is at boundary
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

void determineCells(double x, double y, double z, double gridSize, int numGridX, int numGridY, int numGridZ, int cells[3])
{
	int cellx = (int)((x) / gridSize);
	int celly = (int)((y) / gridSize);
	int cellz = (int)((z) / gridSize);
	
	// temporary fix for if bead is at boundary
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
	
	cells[0] = cellx;
	cells[1] = celly;
	cells[2] = cellz;
}

/*void makelinks(int **link, int numGridX, int numGridY, int numGridZ)
{
	int rxlink[27], rylink[27], rzlink[27];
	int icell[numGridX*numGridY*numGridZ], jcell[numGridX*numGridY*numGridZ], kcell[numGridX*numGridY*numGridZ];

	unsigned int count = 0;
	
	for(int i = -1; i <= 1; i++)
	{
		for(int j = -1; j <= 1; j++)
		{
			for(int k = -1; k <= 1; k++)
			{
				rxlink[count] = i;
				rylink[count] = j;
				rzlink[count] = k;
				
				count++;
			}
		}
	}

	count = 0;
	for(int i = 0; i < numGridX; i++)
	{
		for(int j = 0; j < numGridY; j++)
		{
			for (int k = 0; k < numGridZ; k++)
			{
				icell[count] = i;
				jcell[count] = j;
				kcell[count] = k;
				count++;
			}
		}
	}
	
	//cout << "C: " << count << endl;
	
	for(int i = 0; i < count; i++)
	{
		for(int j = 0; j < 27; j++)
		{
			int cellx = icell[i] + rxlink[j];
			int celly = jcell[i] + rylink[j];
			int cellz = kcell[i] + rzlink[j];

			cellx = cellx - numGridX*(int)floor((double)(cellx)/(double)(numGridX));
			celly = celly - numGridY*(int)floor((double)(celly)/(double)(numGridY));
			cellz = cellz - numGridZ*(int)floor((double)(cellz)/(double)(numGridZ));

			int b = (cellx*numGridY*numGridX) + (celly)*numGridY + cellz;

			//cout << b << endl;
			link[i][j] = b;
		}
		//exit(0);
	}
}*/


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
	srand48((long int)time(NULL));
	
	string atomfilename;
	string xyzfilename;
	
	double boxlx = 109.966122807;
	double boxly = 109.966122807;
	double boxlz = 109.966122807;
	
	// read in input.dat
	ifstream inputfile ("psd.in");
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
	//exit(0);
	// read in a file with bead diameters
	double largestBeadDiameter = 0.0;
	
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
					//cout << tempString << endl;
				
					beadSizes.push_back(beadType());
					beadSizes[beadSizes.size() - 1].id = tempString;
					beadSizes[beadSizes.size() - 1].size = tempSize;
					
					if(tempSize > largestBeadDiameter)
					{
						largestBeadDiameter = tempSize;
					}
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
				//cout << tempString << " " << tempX << " " << tempY << " " << tempZ << endl;
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
	
	//exit(0);
	
	double gridSize = 5;
	int numGridX = (int)(boxlx / gridSize);
	int numGridY = (int)(boxly / gridSize);
	int numGridZ = (int)(boxlz / gridSize);
	

	std::vector<int> *grid = new std::vector<int>[numGridX*numGridY*numGridZ];
	
	// put system beads into a grid
	for(int i = 0; i < systemBeads.size(); i++)
	{		
		//cout << i << endl;
		int bin = determineBin(systemBeads[i].x, systemBeads[i].y, systemBeads[i].z, gridSize, numGridX, numGridY, numGridZ);
		
		//cout << bin << " " << numGridX*numGridY*numGridZ << endl;
		
		(grid[bin]).push_back(i);
	}
	
	/*int **gridLinks = (int**)malloc(numGridX*numGridY*numGridZ*sizeof(int*));
	for(int i = 0; i < numGridX*numGridY*numGridZ; i++)
	{
		gridLinks[i] = (int*)malloc(27 * sizeof(int));
	} 
	
	makelinks(gridLinks, numGridX, numGridY, numGridZ);*/
	
	//exit(0);
	double cubeleteSize = 0.2;
	
	int numCubeX = (int)(boxlx / cubeleteSize);
	int numCubeY = (int)(boxly / cubeleteSize);
	int numCubeZ = (int)(boxlz / cubeleteSize);
	

	//stores the distances to closest system bead from position associated with [i*numCubeX^2 + j*numCubeY + k]
	int cubeleteListSize = numCubeX*numCubeY*numCubeZ;
	struct cubelete *cubeleteList = new struct cubelete[cubeleteListSize];
	for(int i = 0; i < cubeleteListSize; i++)
	{
		cubeleteList[i].dist = boxlx;
	}
	
	cout << "Number of Cubeletes: " << cubeleteListSize << endl;
	cout << "Memory for Cubelete List (GB): " << cubeleteListSize / 8.0 / 1000.0 / 1000.0 / 1000.0 * (64 + 3*2) << endl;
	
	int maxHistDist = 90.0;
	double histStep = 1.0;
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
			//cout << "j: " << j << endl;
			double tempY = (j + 0.5) * cubeleteSize;
			for(int k = 0; k < numCubeZ; k++)
			{
				//cout << "k: " << k << endl;
				double tempZ = (k + 0.5) * cubeleteSize;
				
				int index = i*numCubeX*numCubeX + j*numCubeY + k;
				// find distance from [i*numCubeX^2 + j*numCubeY + k] to nearest system bead
				cubeleteList[index].x = i;
				cubeleteList[index].y = j;
				cubeleteList[index].z = k;
				
				bool notDone = true;
				int maxGrid = 0;
				
				int cells[3];
				
				determineCells(tempX, tempY, tempZ, gridSize, numGridX, numGridY, numGridZ, cells);
				
				while(notDone)
				{
					for(int m = -maxGrid; m <= maxGrid; m++)
					{
						for(int n = -maxGrid; n <= maxGrid; n++)
						{
							for(int o = -maxGrid; o <= maxGrid; o++)
							{
								if(abs(m) == maxGrid || abs(n) == maxGrid || abs(o) == maxGrid)
								{
									int cellx = cells[0] + m;
									int celly = cells[1] + n;
									int cellz = cells[2] + o;
									
									cellx = cellx - numGridX*(int)floor((double)(cellx)/(double)(numGridX));
									celly = celly - numGridY*(int)floor((double)(celly)/(double)(numGridY));
									cellz = cellz - numGridZ*(int)floor((double)(cellz)/(double)(numGridZ));
															
									int tempGrid = (cellx*numGridY*numGridX) + (celly)*numGridY + cellz;
									
									// now search through all systemBeads in tempGrid
									std::vector<int> currGrid = grid[tempGrid];
									for(int v = 0; v < currGrid.size(); v++)
									{
										int tempBead = currGrid[v];
										
										double xDist = tempX - systemBeads[tempBead].x;
										xDist = periodicWrap(xDist, boxlx);
										double yDist = tempY - systemBeads[tempBead].y;
										yDist = periodicWrap(yDist, boxly);
										double zDist = tempZ - systemBeads[tempBead].z;
										zDist = periodicWrap(zDist, boxlz);
										
										double dist = xDist*xDist + yDist*yDist + zDist*zDist;
										dist = sqrt(dist);
										
										dist = dist - 0.5*systemBeads[tempBead].size;
										
										if(dist < cubeleteList[index].dist)
										{
											cubeleteList[index].dist = dist;
										}
									}
								}
							}
						}
					}
					
					/* if the minimum distance stored is less than the possible distance from cubelette center to surface of a system bead
					   lying outside the already considered grids then don't need to search further grids*/
					if(cubeleteList[index].dist <= gridSize*maxGrid - 0.5*(largestBeadDiameter))
					{
						notDone = false;
					}
					else
					{
						//add another layer
						maxGrid += 1;
					}
				}
			}
		}
	}
	
	// sort cubelets using quicksort
	quicksort(0, cubeleteListSize, cubeleteList); 
	
	cout << cubeleteListSize << endl;
	
	
	// randomly pick cubeletes
	int numTrials = 5000;
	
	int maxCubelete = 0;
	
	#pragma omp parallel for
	for(int i = 0; i < numTrials; i++)
	{
		//double tempRand = rand() / 32768.0;
		double tempRand = drand48();
		cout << "i: " << i << " " << tempRand << endl;
		int randomCubeIndex = (int)(cubeleteListSize * tempRand);
		
		//cout << "after rand()" << endl;
		//cout << randomCubeIndex << endl;
		
		double tempX = (cubeleteList[randomCubeIndex].x + 0.5) * cubeleteSize;
		double tempY = (cubeleteList[randomCubeIndex].y + 0.5) * cubeleteSize;
		double tempZ = (cubeleteList[randomCubeIndex].z + 0.5) * cubeleteSize;
		
		//cout << "before loop" << endl;
		
		for(int j = 0; j < cubeleteListSize; j++)
		{
			//cout << "j: " << j << endl;
			double distX = (cubeleteList[j].x + 0.5) * cubeleteSize - tempX;
			distX = periodicWrap(distX, boxlx);
			double distY = (cubeleteList[j].y + 0.5) * cubeleteSize - tempY;
			distY = periodicWrap(distY, boxly);
			double distZ = (cubeleteList[j].z + 0.5) * cubeleteSize - tempZ;
			distZ = periodicWrap(distZ, boxlz);
			
			distX = distX*distX;
			distY = distY*distY;
			distZ = distZ*distZ;
			
			double dist = distX + distY + distZ;
			dist = sqrt(dist);
			if(dist <= cubeleteList[j].dist)
			{
				// this is the largest sphere which randomCubeIndex can fit inside
				if(j > maxCubelete)
				{
					maxCubelete = j;
				}
				
				int tempHistPos = (int)(cubeleteList[j].dist * 2.0 / histStep);
				
				for(int h = 0; h < tempHistPos; h++)
				{
					hist[h] += 1;
				}
				break;
			}
		}
	}
	
	cout << (maxCubelete + 1.0) / (double)cubeleteListSize * 100.0 << "% cubelete list utilization" << endl;
	
	ofstream psdCummulative;
	psdCummulative.open("psd_cummulative.out");
	if(psdCummulative.is_open())
	{
		for (int h = 0; h < histSize; h++)
		{
			psdCummulative << h*histStep + 0.5*histStep << " " << (double)hist[h] / (double)hist[0] << endl;
		}
		
		psdCummulative.close();
	}
	
	// derivative
	ofstream psd;
	psd.open("psd.out");
	if(psd.is_open())
	{
		for (int h = 1; h < histSize - 1; h++)
		{
			double x1 = (h-1)*histStep + 0.5*histStep;
			double x2 = (h+1)*histStep + 0.5*histStep;
			psd << (h)*histStep + 0.5*histStep << " " << -((double)hist[h+1] - (double)hist[h-1] ) / (x2 - x1) / (double)hist[0] << endl;
		}
		
		psd.close();
	}

	return 0;
}