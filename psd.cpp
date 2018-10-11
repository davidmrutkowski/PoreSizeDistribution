//============================================================================
// Name        : psd.cpp
// Author      : David Rutkowski
// Description : Calculates the pore size distribution (psd) based originally on code from Poreblazer (Lev Sarkisov group)
// Compile with "icpc -fopenmp -O2 -o psd psd.cpp" or equivalent g++ expression
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <set>
#include <ctime>
#include <chrono>

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

struct point {
	int x;
	int y;
	int z;
};

int determineBin(double x, double y, double z, double gridSizeX, double gridSizeY, double gridSizeZ, int numGridX, int numGridY, int numGridZ)
{
	int cellx = (int)(x / gridSizeX);
	int celly = (int)(y / gridSizeY);
	int cellz = (int)(z / gridSizeZ);
	
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
	
	int bin = (cellx*numGridY*numGridZ) + (celly)*numGridZ + cellz;
	
	return bin;
}

void determineCells(double x, double y, double z, double gridSizeX, double gridSizeY, double gridSizeZ, int numGridX, int numGridY, int numGridZ, int cells[3])
{
	int cellx = (int)(x / gridSizeX);
	int celly = (int)(y / gridSizeY);
	int cellz = (int)(z / gridSizeZ);
	
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
	string atomfilename;
	string xyzfilename;
	
	double boxlx = 10.0;
	double boxly = 10.0;
	double boxlz = 10.0;
	
	// doesnt do anything with these values currently, always assumes box is rectangular
	double boxAngleX, boxAngleY, boxAngleZ;
	
	int numTrials = 1295029;
	
	// minHistDist doesn't do anything currently and is effectively 0.0 regardless of input file
	double minHistDist = 0.0;
	double maxHistDist = 100.0;
	double histStep = 0.1;
		
	double cubeleteSize = 1.0;
	
	long int randomSeed = 0;
	
	bool outputXYZ = false;
	
	// read in input.dat
	int fileLineCount = 0;
	std::string line;
	char chars[] = "'";
	ifstream inputfile ("psd.in");
	if(inputfile.is_open())
	{
		while(std::getline(inputfile,line))
		{
			if(!(line.empty() || line.substr(0,1).compare("!") == 0))
			{
				std::istringstream iss(line);
			
				if(fileLineCount == 0)
				{
					iss >> atomfilename;
				}
				else if(fileLineCount == 1)
				{
					iss >> xyzfilename;
					//https://stackoverflow.com/questions/5891610/how-to-remove-certain-characters-from-a-string-in-c
					for (unsigned int i = 0; i < strlen(chars); i++)
					{
						xyzfilename.erase (std::remove(xyzfilename.begin(), xyzfilename.end(), chars[i]), xyzfilename.end());
					}
				}
				else if(fileLineCount == 2)
					iss >> numTrials;
				else if(fileLineCount == 3)
					iss >> minHistDist;
				else if(fileLineCount == 4)
					iss >> histStep;
				else if(fileLineCount == 5)
					iss >> maxHistDist;
				else if(fileLineCount == 6)
					iss >> boxlx >> boxly >> boxlz;
				else if(fileLineCount == 7)
					iss >> boxAngleX >> boxAngleY >> boxAngleZ;
				else if(fileLineCount == 8)
				{
					iss >> randomSeed;
				}
				else if(fileLineCount == 9)
				{
					string tempString;
					iss >> tempString;
					
					if(tempString.compare("T") == 0 || tempString.compare("1") == 0)
					{
						outputXYZ = true;
					}
					
					break;
				}
					
				fileLineCount++;
			}
		}
		inputfile.close();
	}
	else
	{
		cout << "Could not find psd.in" << endl;
		return 0;
	}
		
	//srand48((long int)time(NULL));
	srand48(randomSeed);
	
	// read in a file with bead diameters
	double largestBeadDiameter = 0.0;
	
	int numBeadTypes = 0;
	std::vector<beadType> beadSizes;
	
	fileLineCount = 0;
	
	ifstream beadtypefile (atomfilename.c_str());
	if (beadtypefile.is_open())
	{
		while(std::getline(beadtypefile,line))
		{	
			if(!(line.empty() || line.substr(0,1).compare("!") == 0))
			{
				std::istringstream iss(line);
				
				if(fileLineCount == 0)
				{
					iss >> numBeadTypes;
					fileLineCount++;
				}
				else
				{
					string tempString;
					double tempSize;
					
					if(iss >> tempString >> tempSize)
					{
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
		}
		
		beadtypefile.close();
	}
	else
	{
		cout << "Could not find " << atomfilename << endl;
		return 0;
	}
	
	
	double minX = 0.0;
	double minY = 0.0;
	double minZ = 0.0;
	
	// read in an .xyz file with system coordinates
	fileLineCount = 0;
	int numBeads;
	std::vector<bead> systemBeads;
	ifstream xyzfile (xyzfilename.c_str());
	if (xyzfile.is_open())
	{
		while(std::getline(xyzfile,line))
		{
			if(!(line.empty() || line.substr(0,1).compare("!") == 0))
			{
				std::istringstream iss(line);
				
				if (fileLineCount == 0)
				{
					iss >> numBeads;
					
					fileLineCount++;
				}
				else
				{
					string tempString;
					double tempX, tempY, tempZ;
					
					if(iss >> tempString >> tempX >> tempY >> tempZ)
					{					
						systemBeads.push_back(bead());
						systemBeads[systemBeads.size() - 1].id = tempString;
						systemBeads[systemBeads.size() - 1].x = tempX;
						systemBeads[systemBeads.size() - 1].y = tempY;
						systemBeads[systemBeads.size() - 1].z = tempZ;
						
						if(tempX < minX)
							minX = tempX;
						if(tempY < minY)
							minY = tempY;
						if(tempZ < minZ)
							minZ = tempZ;
						
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
					else
					{
						cout << "Ignoring line: " << line << endl;
					}	
				}
			}
		}
		xyzfile.close();
	}
	else
	{
		cout << "Could not find " << xyzfilename << endl;
		return 0;
	}

	// done reading from files
	
	// shift simulation box if necessary (any negative coordinates)
	for(int i = 0; i < systemBeads.size(); i++)
	{
		systemBeads[i].x -= minX;
		systemBeads[i].y -= minY;
		systemBeads[i].z -= minZ;
	}
	
	// should determine if grid is needed based on # of system beads (<1000 dont bother with grid?)
	// should automatically determine this based on numbered of particles evenly dispersed in system
	int desiredNumParticlesInGrid = 20;
	
	double systemDensity = systemBeads.size() / boxlx / boxly / boxlz;
	double gridSize = (float)desiredNumParticlesInGrid / systemDensity;
	gridSize = pow(gridSize, 1.0/3.0);
	cout << "System gridsize is: " << gridSize << endl;
	
	int numGridX = (int)(boxlx / gridSize);
	int numGridY = (int)(boxly / gridSize);
	int numGridZ = (int)(boxlz / gridSize);
	
	if(numGridX < 3)
		numGridX = 3;
	if(numGridY < 3)
		numGridY = 3;
	if(numGridZ < 3)
		numGridZ = 3;
	
	double gridSizeX = boxlx / (double)numGridX;
	double gridSizeY = boxly / (double)numGridY;
	double gridSizeZ = boxlz / (double)numGridZ;
	
	double minGridSize = gridSizeX;
	if(gridSizeY < minGridSize)
		minGridSize = gridSizeY;
	if(gridSizeZ < minGridSize)
		minGridSize = gridSizeZ;
	
	std::vector<int> *grid = new std::vector<int>[numGridX*numGridY*numGridZ];
	
	// put system beads into a grid
	for(int i = 0; i < systemBeads.size(); i++)
	{
		int bin = determineBin(systemBeads[i].x, systemBeads[i].y, systemBeads[i].z, gridSizeX, gridSizeY, gridSizeZ, numGridX, numGridY, numGridZ);
		
		if(bin >= numGridX*numGridY*numGridZ)
		{
			cout << "Error bin too big: " << bin << endl;
			cout << systemBeads[i].x << " " << systemBeads[i].y << " " << systemBeads[i].z << endl;
			cout << boxlx << " " << boxly << " " << boxlz << endl;
		}
		
		(grid[bin]).push_back(i);
	}
	
	int numCubeX = (int)(boxlx / cubeleteSize);
	int numCubeY = (int)(boxly / cubeleteSize);
	int numCubeZ = (int)(boxlz / cubeleteSize);
	
	if(numCubeX < 1)
		numCubeX = 1;
	if(numCubeY < 1)
		numCubeY = 1;
	if(numCubeZ < 1)
		numCubeZ = 1;
	
	//stores the distances to closest system bead from position associated with [i*numCubeX^2 + j*numCubeY + k]
	int cubeleteListSize = numCubeX*numCubeY*numCubeZ;
	struct cubelete *cubeleteList = new struct cubelete[cubeleteListSize];
	for(int i = 0; i < cubeleteListSize; i++)
	{
		cubeleteList[i].dist = boxlx*boxly*boxlz;
	}
	
	cout << "Number of Cubeletes: " << cubeleteListSize << endl;
	cout << "Memory for Cubelete List (GB): " << cubeleteListSize / 8.0 / 1000.0 / 1000.0 / 1000.0 * (64 + 3*8) << endl;
	
	// generate grid neighbor mapping	
	int xMin = (int)(numGridX*0.5);
	int yMin = (int)(numGridY*0.5);
	int zMin = (int)(numGridZ*0.5);
	
	int xMax = (int)round(numGridX*0.5);
	int yMax = (int)round(numGridY*0.5);
	int zMax = (int)round(numGridZ*0.5);
	
	int maxGridNeighbor = xMin;
	if(yMin > maxGridNeighbor)
		maxGridNeighbor = yMin;
	if(zMin > maxGridNeighbor)
		maxGridNeighbor = zMin;
	maxGridNeighbor = maxGridNeighbor + 1;
	
	
	vector<vector<struct point> > gridNeighbors(maxGridNeighbor);
	
	for(int m = -xMin; m < xMax; m++)
	{
		for(int n = -yMin; n < yMax; n++)
		{
			for(int o = -zMin; o < zMax; o++)
			{
				int largest = abs(m);
				if(largest < abs(n))
					largest = abs(n);
				if(largest < abs(o))
					largest = abs(o);
				
				struct point tempPoint;
				tempPoint.x = m;
				tempPoint.y = n;
				tempPoint.z = o;
				
				gridNeighbors[largest].push_back(tempPoint);
			}
		}
	}
	
	//auto start = std::chrono::system_clock::now();
	
	#pragma omp parallel for
	for(int i = 0; i < numCubeX; i++)
	{		
		double tempX = (i + 0.5) * cubeleteSize;
		for(int j = 0; j < numCubeY; j++)
		{
			double tempY = (j + 0.5) * cubeleteSize;
			for(int k = 0; k < numCubeZ; k++)
			{
				double tempZ = (k + 0.5) * cubeleteSize;
				
				int index = i*numCubeY*numCubeZ + j*numCubeZ + k;
				
				// find distance from [i*numCubeX^2 + j*numCubeY + k] to nearest system bead
				cubeleteList[index].x = i;
				cubeleteList[index].y = j;
				cubeleteList[index].z = k;
				
				bool done = false;
				int maxGridX = 0;
				int maxGridY = 0;
				int maxGridZ = 0;
				
				int cells[3];
				
				determineCells(tempX, tempY, tempZ, gridSizeX, gridSizeY, gridSizeZ, numGridX, numGridY, numGridZ, cells);
				
				int layer = 0;
				// need to ensure that maxGrid does not go larger than it can possibly be!
				while(done == false && layer < maxGridNeighbor)
				{
					for(int gn = 0; gn < gridNeighbors[layer].size(); gn++)
					{
						int m = gridNeighbors[layer][gn].x;
						int n = gridNeighbors[layer][gn].y;
						int o = gridNeighbors[layer][gn].z;
						
						int cellx = cells[0] + m;
						int celly = cells[1] + n;
						int cellz = cells[2] + o;
						
						cellx = cellx - numGridX*(int)floor((double)(cellx)/(double)(numGridX));
						celly = celly - numGridY*(int)floor((double)(celly)/(double)(numGridY));
						cellz = cellz - numGridZ*(int)floor((double)(cellz)/(double)(numGridZ));
												
						int tempGrid = (cellx*numGridY*numGridZ) + (celly)*numGridZ + cellz;
						
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
					
					// if the minimum distance stored is less than the possible distance from cubelette center to surface of a system bead
					// lying outside the already considered grids then don't need to search further grids
					if(cubeleteList[index].dist <= minGridSize*layer - 0.5*(largestBeadDiameter))
					{
						//cout << "layer: " << layer << endl;
						done = true;
					}
					else
					{
						//add another grid layer
						layer += 1;
					}
				}
			}
		}
	}
	
	// sort cubelets using quicksort
	quicksort(0, cubeleteListSize-1, cubeleteList); 	
	
	// randomly pick cubeletes	
	int histSize = (int)((maxHistDist)/ histStep);
	int *hist = new int[histSize];
	for(int h = 0; h < histSize; h++)
	{
		hist[h] = 0;
	}
	
	int maxCubelete = 0;
	
	int zeroCubeletePosition = cubeleteListSize - 1;
	// search through cubeleteList backwards to find when the position where values start being positive or zero
	for(int i = cubeleteListSize - 1; i >= 0; i--)
	{
		if(cubeleteList[i].dist >= 0.0)
		{
			zeroCubeletePosition = i;
			break;
		}
	}
	
	//std::set<int> outputCubeletes;
	int *tempHistArray = new int[numTrials];
	int *outputCubeletes = new int[numTrials];
	
	//cout << zeroCubeletePosition << endl;
	
	#pragma omp parallel for
	for(int i = 0; i < numTrials; i++)
	{
		//double tempRand = rand() / 32768.0;
		double tempRand = drand48();
		//cout << "i: " << i << " " << tempRand << endl;
		int randomCubeIndex = (int)(zeroCubeletePosition * tempRand);
		
		double tempX = (cubeleteList[randomCubeIndex].x + 0.5) * cubeleteSize;
		double tempY = (cubeleteList[randomCubeIndex].y + 0.5) * cubeleteSize;
		double tempZ = (cubeleteList[randomCubeIndex].z + 0.5) * cubeleteSize;
		
		for(int j = 0; j <= zeroCubeletePosition; j++)
		{
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
				
				// this is not thread safe!
				//outputCubeletes.insert(j);
				
				// this is thread safe!
				outputCubeletes[i] = randomCubeIndex;
				
				if(j > maxCubelete)
				{
					maxCubelete = j;
				}
				
				int tempHistPos = (int)((cubeleteList[j].dist * 2.0) / histStep);
				
				
				tempHistArray[i] = tempHistPos;
				
				break;
			}
		}
	}
	
	//auto end = std::chrono::system_clock::now();

	//auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	//std::cout << elapsed.count() << endl;
	
	for(int i = 0; i < numTrials; i++)
	{
		int tempHistPos = tempHistArray[i];
		
		for(int h = 0; h < tempHistPos; h++)
		{
			hist[h] += 1;
		}
	}
	
	cout << (maxCubelete + 1.0) / (double)zeroCubeletePosition * 100.0 << "% cubelete list utilization" << endl;
	
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
	
	// output .xyz
	std::set<int> poreCubeletes;
	for(int i = 0; i < numTrials; i++)
	{
		poreCubeletes.insert(outputCubeletes[i]);
	}
	
	if(outputXYZ == true)
	{
		ofstream xyzFile;
		xyzFile.open("Pores.xyz");
		if(xyzFile.is_open())
		{
			for(std::set<int>::iterator it=poreCubeletes.begin(); it != poreCubeletes.end(); it++)
			{
				int currIndex = *it;
				
				double tempX = (cubeleteList[currIndex].x + 0.5) * cubeleteSize;
				double tempY = (cubeleteList[currIndex].y + 0.5) * cubeleteSize;
				double tempZ = (cubeleteList[currIndex].z + 0.5) * cubeleteSize;
				
				double tempRadius = cubeleteList[currIndex].dist;
				
				xyzFile << "Ar" << " " << tempX << " " << tempY << " " << tempZ << " " << tempRadius << endl;
			}
			
			xyzFile.close();
		}
	}
	
	// output shifted coordinated, new.xyz
	ofstream newXYZFile;
	newXYZFile.open("new.xyz");
	if(newXYZFile.is_open())
	{
		newXYZFile << systemBeads.size() << endl << endl;
		
		for(int i = 0; i < systemBeads.size(); i++)
		{			
			newXYZFile << systemBeads[i].id << " " << systemBeads[i].x << " " << systemBeads[i].y << " " << systemBeads[i].z << endl;
		}
		
		newXYZFile.close();
	}

	return 0;
}