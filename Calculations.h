#pragma once

#include <algorithm>
#include <fstream>
#include "PoscarInfo.h"

//Header file for function definitions of various manipulations of POSCAR and CONTCAR files

//Globals:  stuff to change paramaters as needed
double GENERAL_BOND_DISTANCE = 3; ///(angstroms) distance between any two atoms for them to be considered bonded (works well for POSCARS with only a few atom types)

struct bondInfo
{
	std::string type;
	double physBondDist;
};

//Struct to make life easier:  specifically for writing bond info
struct atom
{
	int id;
	std::vector <atomPair> involvedIn;
};

//rules for sorting vector of bondinfo:  sorts by type, then id, then 
bool compareBondInfoRules(bondInfo a, bondInfo b)
{
	return a.type < b.type;
}

//rules for sorting atom pairs
bool compareAtomPair(atomPair a, atomPair b)
{
	return (a.pairedAtoms[0].atomType + a.pairedAtoms[1].atomType) < (b.pairedAtoms[0].atomType + b.pairedAtoms[1].atomType);
}

//Comparison rule for sorting strings
bool stringSort(std::string a, std::string b)
{
	return a < b;
}

//Distance between two points in R3
double dist(const Coords coordA, const Coords coordB)
{
	return sqrt((((coordB.a - coordA.a)*(coordB.a - coordA.a)))+(((coordB.b - coordA.b)*(coordB.b - coordA.b)))+(((coordB.c - coordA.c)*(coordB.c - coordA.c))));
}

//TODO: make this work for any shape - not just orthogonal cartesian coords
//Extends supercell in potentilly three directions (along the supercell unit vectors).  Marks original atoms and new atoms respectivly
//Sytax: negative, then positive for x, y, z each in that order
void Poscar::extendSupercell(unsigned int negX, unsigned int posX, unsigned int negY, unsigned int posY, unsigned int negZ, unsigned int posZ)
{
	//Vectors for holding new values
	std::vector <Coords> negXCrds, posXCrds, negYCrds, posYCrds, negZCrds, posZCrds;

	//Convert to cartesian for easier algorythm writing.  Also takes care of universal scale constant.  Find original atoms
	convertToCartesian();
	for (int i = 0; i < atomCoords.size(); i++) ///mark original atoms as such
		atomCoords[i].extraInfo = "original";
	

	//Filling all vectors with their respective coordinates...

	//X coords:
	if (negX != 0)
		for (int i = 0; i < negX; i++)
			for (int j = 0; j < atomCoords.size(); j++)
			{
				Coords tmp = atomCoords[j]; ///copy each atom in original coords (including the ID - this is important)
				tmp.a = tmp.a - (superCellVectorA[0] * (i + 1)); ///changing corresponding value
				tmp.extraInfo[0] = 'n'; ///tagging atoms
				negXCrds.push_back(tmp);
			}
	//filling atomCoords
	for (int i = 0; i < negXCrds.size(); i++)
		atomCoords.push_back(negXCrds[i]);

	if (posX != 0)
		for (int i = 0; i < posX; i++)
			for (int j = 0; j < atomCoords.size(); j++)
			{
				Coords tmp = atomCoords[j]; ///copy each atom in original coords (including the ID - this is important)
				tmp.a = tmp.a + (superCellVectorA[0] * (i + 1)); ///changing corresponding value
				tmp.extraInfo[0] = 'p'; ///tagging atoms
				posXCrds.push_back(tmp);
			}
	//filling atomCoords
	for (int i = 0; i < posXCrds.size(); i++)
		atomCoords.push_back(posXCrds[i]);

	//Y coords:
	if (negY != 0)
		for (int i = 0; i < negY; i++)
			for (int j = 0; j < atomCoords.size(); j++)
			{
				Coords tmp = atomCoords[j]; ///copy each atom in original coords (including the ID - this is important)
				tmp.b = tmp.b - (superCellVectorB[1] * (i + 1)); ///changing corresponding value
				tmp.extraInfo[1] = 'n'; ///tagging atoms
				negYCrds.push_back(tmp);
			}
	//filling atomCoords
	for (int i = 0; i < negYCrds.size(); i++)
		atomCoords.push_back(negYCrds[i]);

	if (posY != 0)
		for (int i = 0; i < posY; i++)
			for (int j = 0; j < atomCoords.size(); j++)
			{
				Coords tmp = atomCoords[j]; ///copy each atom in original coords (including the ID - this is important)
				tmp.b = tmp.b + (superCellVectorB[1] * (i + 1)); ///changing corresponding value
				tmp.extraInfo [1] = 'p'; ///tagging atoms
				posYCrds.push_back(tmp);
			}
	//filling atomCoords
	for (int i = 0; i < posYCrds.size(); i++)
		atomCoords.push_back(posYCrds[i]);

	//Z coords:
	if (negZ != 0)
		for (int i = 0; i < negZ; i++)
			for (int j = 0; j < atomCoords.size(); j++)
			{
				Coords tmp = atomCoords[j]; ///copy each atom in original coords (including the ID - this is important)
				tmp.c = tmp.c - (superCellVectorC[2] * (i + 1)); ///changing corresponding value
				tmp.extraInfo[2] = 'n'; ///tagging atoms
				negZCrds.push_back(tmp);
			}
	//filling atomCoords
	for (int i = 0; i < negZCrds.size(); i++)
		atomCoords.push_back(negZCrds[i]);

	if (posZ != 0)
		for (int i = 0; i < posZ; i++)
			for (int j = 0; j < atomCoords.size(); j++)
			{
				Coords tmp = atomCoords[j]; ///copy each atom in original coords (including the ID - this is important)
				tmp.c = tmp.c + (superCellVectorC[2] * (i + 1)); ///changing corresponding value
				tmp.extraInfo[2] = 'p'; ///tagging atoms
				posZCrds.push_back(tmp);
			}
	//filling atomCoords
	for (int i = 0; i < posZCrds.size(); i++)
		atomCoords.push_back(posZCrds[i]);

	//extending the supercell vectors:
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = superCellVectorA[i] * ((negX + posX) + 1);
		superCellVectorB[i] = superCellVectorB[i] * ((negY + posY) + 1);
		superCellVectorC[i] = superCellVectorC[i] * ((negZ + posZ) + 1);
	}

	//remove duplicate atoms that inevitably exist because I wrote this in like 20 minutes
	//removeDuplicates();

	//Updating atom type counts
	updateAtomCounts();
}

//Important:  this function assumes:  the shape is two-dimensional parallelogram in the xy plane, supercell vector A sits along the x axis
double Poscar::grapheneAreaApprox()
{
	convertToCartesian();
	
	if (modelType == "Bulk")
	{
		std::vector <double> vectA, vectB, vectC;

		for (int i = 0; i < 3; i++)
		{
			vectA.push_back(superCellVectorA[i]);
			vectB.push_back(superCellVectorB[i]);
			vectC.push_back(superCellVectorC[i]);
		}

		double a = sqrt((superCellVectorA[0] * superCellVectorA[0]) + (superCellVectorA[1] * superCellVectorA[1]) + (superCellVectorA[2] * superCellVectorA[2]));
		double b = sqrt((superCellVectorB[0] * superCellVectorB[0]) + (superCellVectorB[1] * superCellVectorB[1]) + (superCellVectorB[2] * superCellVectorB[2]));
		double c = sqrt((superCellVectorC[0] * superCellVectorC[0]) + (superCellVectorC[1] * superCellVectorC[1]) + (superCellVectorC[2] * superCellVectorC[2]));
		double alpha = angleBetween(vectB, vectC); //the angle between b and c
		double beta = angleBetween(vectA, vectC); //the angle between a and c
		double gamma = angleBetween(vectA, vectB); //the angle between a and b
		double volume = a*b*c*sqrt(1 - (cos(alpha) * cos(alpha)) - (cos(beta) * cos(beta)) - (cos(gamma) * cos(gamma)) + 2 * cos(alpha)*cos(beta)*cos(gamma));
		return volume / (superCellVectorC[2]);
	}

	else
		if (modelType == "Molecular")
		{
			double minX, minY, maxX, maxY;
			minX = minY = 1000;
			maxX = maxY = -1000;

			for (int i = 0; i < atomCoords.size(); i++)
				if (atomCoords[i].atomType == "C")
				{
					if (atomCoords[i].a < minX)
						minX = atomCoords[i].a;
					if (atomCoords[i].b < minY)
						minY = atomCoords[i].b;
					if (atomCoords[i].a > maxX)
						maxX = atomCoords[i].a;
					if (atomCoords[i].b > maxY)
						maxY = atomCoords[i].b;
				}

			return (maxX - minX)*(maxY - minY);
		}

	return -5;
}

//really, just writing this to use when i figure out a way to determine verticies for POSCAR files
//better approximation for area.  but needs verticies (see below syntax).  Still assumes it is a planar 4-vertexed shape.  Assumes coordinate pairs are in cartesian
double areaGivenVerticies(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
	double area = 0;         
	double X [4] = {x1, x2, x3, x4}; ///arrays to make writing the algorythm less annoying
	double Y [4] = {y1, y2, y3, y4}; ///Y is parallel to X

	for (int i = 0, j = 3; i < 4; j = i, i++)
		area = area + (X[i] + X[j]) * (Y[j] - Y[i]); ///I'll have the figuring for this algorythm in my research notebook in case there is debugging to do, 
													///but for the purposes of the Graphene stuff, it should work perfectly fine
	return area / 2;
}

//Shifts the atoms in the supercell
void Poscar::translateAtoms(double xTrans, double yTrans, double zTrans)
{
	convertToCartesian();

	for (int i = 0; i < atomCoords.size(); i++)
	{
		atomCoords[i].a += xTrans;
		atomCoords[i].b += yTrans;
		atomCoords[i].c += zTrans;
	}
}

//Function to apply strain to a model.  Note: even numbers of steps are preferred, the middle number will be the exact same as the input file
void Poscar::applyStrain(bool strainX, bool strainY, bool strainZ, int nSteps, double maxStrainMagnitude)
{
	Poscar origFile("readAll", infilePath);
	
	//Preliminary step:  if molecular...
	double minX, minY, minZ, maxX, maxY, maxZ;
	if (modelType == "Molecular")
	{
		convertToCartesian();

		//Shift the model to the bottom-left of the (hopefully rectangular) supercell - this is equivalent to 'fixing' the edges constant
		//(This is outside the main loop to increase preformance)
		///First, find how far to shift TODO: (should really have a seperate function for this)
		
		minX = minY = minZ = 1000;
		maxX = maxY = maxZ = -1000;
		///find max and mins
		for (int i = 0; i < origFile.atomCoords.size(); i++)
		{
			if (origFile.atomCoords[i].a < minX)
				minX = origFile.atomCoords[i].a;
			if (origFile.atomCoords[i].b < minY)
				minY = origFile.atomCoords[i].b;
			if (origFile.atomCoords[i].c < minZ)
				minZ = origFile.atomCoords[i].c;
			if (origFile.atomCoords[i].a > maxX)
				maxX = origFile.atomCoords[i].a;
			if (origFile.atomCoords[i].b > maxY)
				maxY = origFile.atomCoords[i].b;
			if (origFile.atomCoords[i].c > maxZ)
				maxZ = origFile.atomCoords[i].c;
		}
		///Move the atoms
		origFile.translateAtoms(-minX, -minY, -minZ);
	}

	//MAIN LOOP
	for (int i = 0; i <= nSteps; i++)
	{
		//Copy contents of original poscar into new one	
		Poscar newFile = origFile;

		//Bulk Case:
		if (newFile.modelType == "Bulk")
		{
			newFile.convertToDirect();
			if (strainX)
			{
				newFile.superCellVectorA[0] = origFile.superCellVectorA[0] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				newFile.superCellVectorB[0] = origFile.superCellVectorB[0] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				newFile.superCellVectorC[0] = origFile.superCellVectorC[0] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
			}
			if (strainY)
			{
				newFile.superCellVectorA[1] = origFile.superCellVectorA[1] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				newFile.superCellVectorB[1] = origFile.superCellVectorB[1] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				newFile.superCellVectorC[1] = origFile.superCellVectorC[1] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
			}
			if (strainZ)
			{
				newFile.superCellVectorA[2] = origFile.superCellVectorA[2] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				newFile.superCellVectorB[2] = origFile.superCellVectorB[2] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				newFile.superCellVectorC[2] = origFile.superCellVectorC[2] * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
			}
		}

		else

		//Molecular Case:
		if (newFile.modelType == "Molecular")
		{
			newFile.convertToCartesian();

			//Change new file coords based on original file coords.  The use of the same index SHOULD be fine, since the vectors are parallel
			for (int j = 0; j < newFile.atomCoords.size(); j++)
			{
				if (strainX)
					newFile.atomCoords[j].a = origFile.atomCoords[j].a * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				if (strainY)
					newFile.atomCoords[j].b = origFile.atomCoords[j].b * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
				if (strainZ)
					newFile.atomCoords[j].c = origFile.atomCoords[j].c * (1 + (-maxStrainMagnitude + (2 * i * maxStrainMagnitude / nSteps)));
			}
		}

		else
		{
			std::cout << "model type not defined." << std::endl;
			break;
		}

		//Write the edit number to a new file.  Need to do the following workaround because the g++ compiler complains about doing it the easy way
		std::stringstream name;
		name << i;
		std::string istr = name.str();
		newFile.write(("Edit" + istr).c_str());
	}

	///Translate the atoms back to where they should go, if the model type was molecular:
	if (origFile.modelType == "Molecular")
		origFile.translateAtoms(minX, minY, minZ);
}

void Poscar::fetchAtomPairs()
{
	if (!atomPairs.empty())
		return;

	convertToCartesian();

	//make a copy of the input file so you dont change the original input instance
	Poscar extendedPoscar = *this;

	//Extend inpt file in all directions
	extendedPoscar.extendSupercell(1, 1, 1, 1, 1, 1); ///this could take a long time to do

	//'shaving' coordinates so that the rest of the calculations go faster:  if any coords are more than (a bit more than) a bond distance away from the supercell, disregard them
	//X direction
	std::vector <Coords> tmpx;
	for (int i = 0; i < extendedPoscar.atomCoords.size(); i++) ///putting atoms worth keeping in a vector
	{
		///positive x check
		if (extendedPoscar.atomCoords[i].extraInfo[0] == 'p')
		{
			if (extendedPoscar.atomCoords[i].a - GENERAL_BOND_DISTANCE < superCellVectorA[0] + 0.1)
				tmpx.push_back(extendedPoscar.atomCoords[i]);
		}
		else
			///negative x check
			if (extendedPoscar.atomCoords[i].extraInfo[0] == 'n')
			{
				if (abs(extendedPoscar.atomCoords[i].a) < GENERAL_BOND_DISTANCE +  0.1)
					tmpx.push_back(extendedPoscar.atomCoords[i]);
			}
			else
				///original check 
				if (extendedPoscar.atomCoords[i].extraInfo == "original")
					tmpx.push_back(extendedPoscar.atomCoords[i]);
	}
	//Y direction
	std::vector <Coords> tmpy;
	for (int i = 0; i < tmpx.size(); i++) ///putting atoms worth keeping in a vector
	{
		///positive y check
		if (tmpx[i].extraInfo[1] == 'p')
		{
			if (tmpx[i].b - GENERAL_BOND_DISTANCE < superCellVectorB[1] + 0.1)
				tmpy.push_back(tmpx[i]);
		}
		else
			///negative y check
			if (tmpx[i].extraInfo[1] == 'n')
			{
				if (abs(tmpx[i].b) < GENERAL_BOND_DISTANCE + 0.1)
					tmpy.push_back(tmpx[i]);
			}
			else
				///original check 
				if (tmpx[i].extraInfo == "original")
					tmpy.push_back(tmpx[i]);
	}
	//Z direction
	std::vector <Coords> tmpz;
	for (int i = 0; i < tmpy.size(); i++) ///putting atoms worth keeping in a vector
	{
		///positive z check
		if (tmpy[i].extraInfo[2] == 'p')
		{
			if (tmpy[i].c - GENERAL_BOND_DISTANCE < superCellVectorC[2] + 0.1)
				tmpz.push_back(tmpy[i]);
		}
		else
			///negative z check
			if (tmpy[i].extraInfo[2] == 'n')
			{
				if (abs(tmpy[i].c) < GENERAL_BOND_DISTANCE + 0.1)
					tmpz.push_back(tmpy[i]);
			}
			else
				///original check 
				if (tmpy[i].extraInfo == "original")
					tmpz.push_back(tmpy[i]);
	}

	//Updating atom type counts
	extendedPoscar.atomCoords = tmpz;
	extendedPoscar.updateAtomCounts(tmpz);

	int count = 0;
	std::vector <atomPair> tmp2;
	for (int i = 0; i < extendedPoscar.atomCoords.size(); i++)
		if (extendedPoscar.atomCoords[i].extraInfo == "original")
		{
			count++;
			for (int j = 0; j < extendedPoscar.atomCoords.size(); j++)
				if (extendedPoscar.atomCoords[i] != extendedPoscar.atomCoords[j])
					if ((new atomPair(extendedPoscar.atomCoords[i], extendedPoscar.atomCoords[j]))->lenBetween < GENERAL_BOND_DISTANCE)
						if ((new atomPair(extendedPoscar.atomCoords[i], extendedPoscar.atomCoords[j]))->lenBetween > 0.1)
							tmp2.push_back(*new atomPair(extendedPoscar.atomCoords[i], extendedPoscar.atomCoords[j]));
		}
	/*
	//Get rid of double counting
	{
	std::vector <atomPair> tmp1;
	for (int i = 0; i < tmp2.size(); i++)
		for (int j = 0; j < tmp2.size(); j++)
			if (i != j)
				if (tmp2[i] == tmp2[j] && tmp2[i].extraInfo != "cp")
					tmp2[j].extraInfo = "cp";

	for (int i = 0; i < tmp2.size(); i++)
		if (tmp2[i].extraInfo != "cp")
			tmp1.push_back(tmp2[i]);

	tmp2 = tmp1;
	}
	
	//TODO Initialize fraction of each bond that is in the unit cell

	//Remove any bonds that have negative components
	std::vector <atomPair> tmp1;
	for (int i = 0; i < tmp2.size(); i++)
		if (tmp2[i].pairedAtoms[0].extraInfo[0] == 'n' || tmp2[i].pairedAtoms[1].extraInfo[0] == 'n' || tmp2[i].pairedAtoms[0].extraInfo[1] == 'n' || tmp2[i].pairedAtoms[1].extraInfo[1] == 'n' || tmp2[i].pairedAtoms[0].extraInfo[2] == 'n' || tmp2[i].pairedAtoms[1].extraInfo[2] == 'n')
					tmp2[i].extraInfo[4] = 'r';

	for (int i = 0; i < tmp2.size(); i++)
		if (tmp2[i].extraInfo[4] != 'r')
			tmp1.push_back(tmp2[i]);

	tmp2 = tmp1;
	*/
	atomPairs = tmp2;
	updateAtomCounts();
}

void sortAtomPairs(std::vector <atomPair> &a)
{
	for (int i = 0; i > a.size(); i++)
		if (a[i].pairedAtoms[0].atomType < a[i].pairedAtoms[1].atomType)
		{
			Coords tmp = a[i].pairedAtoms[1];
			a[i].pairedAtoms[1] = a[i].pairedAtoms[0];
			a[i].pairedAtoms[0] = tmp;
		}
}

//REQUIRES info for bond lengths
void Poscar::fetchAtomBonds()
{
	if (!atomBonds.empty())
		return;

	fetchAtomPairs();

	std::ifstream infile("C:\\Users\\baron\\Desktop\\bondDataInfo");
	if (infile.fail())
	{
		std::cout << "could not open 'bondDataInfo' for reading.  Make sure a file named 'bondDataInfo' is present in the current directory.\n";
		std::cout << "the file should be in the following format:\nATOM1 ATOM2 (maximum bonding dist between ATOM1 and ATOM2) `newline`\netc.";
	}

	std::vector <bondInfo> bondInfoVect;

	//Fill array of bond info
	while (!infile.eof())
	{
		bondInfo bnd;
		std::string tmp1, tmp2;

		infile >> tmp1 >> tmp2;
		bnd.type = (tmp1 + tmp2).c_str();
		infile >> bnd.physBondDist;

		bondInfoVect.push_back(bnd);
	}

	sortAtomPairs(atomPairs);
	std::sort(atomPairs.begin(), atomPairs.end(), compareAtomPair);
	std::sort(bondInfoVect.begin(), bondInfoVect.end(), compareBondInfoRules);

	//Filling atomBonds vector
	for (int i = 0; i < atomPairs.size(); i++)
		for (int j = 0; j < bondInfoVect.size(); j++)
			if ((atomPairs[i].pairedAtoms[0].atomType + atomPairs[i].pairedAtoms[1].atomType).c_str() == bondInfoVect[j].type
				|| (atomPairs[i].pairedAtoms[1].atomType + atomPairs[i].pairedAtoms[0].atomType).c_str() == bondInfoVect[j].type)
				if (atomPairs[i].lenBetween <= bondInfoVect[j].physBondDist)
					atomBonds.push_back(atomPairs[i]);


	infile.close();

	//properly fil the types for atomBonds
	for (int i = 0; i < atomBonds.size(); i++)
		atomBonds[i].type = (atomBonds[i].pairedAtoms[0].atomType + atomBonds[i].pairedAtoms[1].atomType).c_str();
	atomPair tmp;
	tmp.type = "END";
	atomBonds.push_back(tmp);
}

//REQUIRES info for bond lengths
void Poscar::fetchAtomBonds(std::string path)
{
	if (!atomBonds.empty())
		return;

	fetchAtomPairs();

	std::ifstream infile((path).c_str());
	if (infile.fail())
	{
		std::cout << "could not open 'bondDataInfo' for reading.  Make sure a file named 'bondDataInfo' is present in the current directory.\n";
		std::cout << "the file should be in the following format:\nATOM1 ATOM2 (maximum bonding dist between ATOM1 and ATOM2) `newline`\netc.";
		return;
	}

	std::vector <bondInfo> bondInfoVect;

	//Fill array of bond info
	while (!infile.eof())
	{
		bondInfo bnd;
		std::string tmp1, tmp2;

		infile >> tmp1 >> tmp2;
		bnd.type = (tmp1 + tmp2).c_str();
		infile >> bnd.physBondDist;

		bondInfoVect.push_back(bnd);
	}

	sortAtomPairs(atomPairs);
	std::sort(atomPairs.begin(), atomPairs.end(), compareAtomPair);
	std::sort(bondInfoVect.begin(), bondInfoVect.end(), compareBondInfoRules);

	//Filling atomBonds vector
	for (int i = 0; i < atomPairs.size(); i++)
		for (int j = 0; j < bondInfoVect.size(); j++)
			if ((atomPairs[i].pairedAtoms[0].atomType + atomPairs[i].pairedAtoms[1].atomType).c_str() == bondInfoVect[j].type
				|| (atomPairs[i].pairedAtoms[1].atomType + atomPairs[i].pairedAtoms[0].atomType).c_str() == bondInfoVect[j].type)
				if (atomPairs[i].lenBetween <= bondInfoVect[j].physBondDist && atomPairs[i].lenBetween > 0.1)
					atomBonds.push_back(atomPairs[i]);

	//Remove duplicate bonds
	for (int i = 0; i < atomBonds.size(); i++){
		for (int j = 0; j < atomBonds.size(); j++)
			if (i != j)
				if ((atomBonds[i] &= atomBonds[j]) && (atomBonds[i].extraInfo != "remo"))
					atomBonds[j].extraInfo = "remo";
	}

	std::vector <atomPair> tmp1;
	for (int i = 0; i < atomBonds.size(); i++)
		if (atomBonds[i].extraInfo != "remo")
			tmp1.push_back(atomBonds[i]);

	atomBonds = tmp1;
	

	infile.close();

	//properly fill the types for atomBonds
	for (int i = 0; i < atomBonds.size(); i++)
		atomBonds[i].type = (atomBonds[i].pairedAtoms[0].atomType + atomBonds[i].pairedAtoms[1].atomType).c_str();
	atomPair tmp;
	tmp.type = "END";
	atomBonds.push_back(tmp);
}

//Writes bond densities to a csv file
void Poscar::writeBondInfo(std::string path)
{
	std::ofstream outfile((path).c_str());
	std::vector <atom> typesAndCounts; ///new vector to hold the types and their counts.

	//find number of atoms total
	int totalAtoms = 0;
	for (int i = 0; i < atomTypeNums.size(); i++)
		totalAtoms += atomTypeNums[i];

	//Initialize vector of atoms before filling with information
	for (int i = 0; i < totalAtoms; i++)
	{
		atom tmp;
		tmp.id = i;
		typesAndCounts.push_back(tmp);
	}
	
	//fill the typesAndCounts vector: by filling each ID's vector with the atoms it is involved in
	for (int i = 0; i < atomBonds.size(); i++)
		for (int n = 0; n < 2; n++)
			for (int j = 0; j < typesAndCounts.size(); j++)
				if (atomBonds[i].pairedAtoms[n].id == typesAndCounts[j].id)
					typesAndCounts[j].involvedIn.push_back(atomBonds[i]);

	
	//initialize the counts for the involvedIn bonds, for each atom seperatly
	for (int i = 0; i < typesAndCounts.size(); i++)
		for  (int k = 0; k < typesAndCounts[i].involvedIn.size(); k++)
			for (int m = 0; m < typesAndCounts[i].involvedIn.size(); m++)
				if (typesAndCounts[i].involvedIn[k].type == typesAndCounts[i].involvedIn[m].type)
					typesAndCounts[i].involvedIn[k].genNum++;
				

	//get rid of extras
	for (int i = 0; i < typesAndCounts.size(); i++)
		for (int j = 0; j < typesAndCounts[i].involvedIn.size(); j++)
			for (int k = 0; k < typesAndCounts[i].involvedIn.size(); k++)
				if (typesAndCounts[i].involvedIn[j].type == typesAndCounts[i].involvedIn[k].type && typesAndCounts[i].involvedIn[j].type != "na" && j != k)
					typesAndCounts[i].involvedIn[j].type = "na";
			


	std::vector <std::string> toWrite;
	//change the names of the types of bonds (add a number to show how many of each type an atom has), an prepare to write
	for (int i = 0; i < typesAndCounts.size(); i++)
		for (int j = 0; j < typesAndCounts[i].involvedIn.size(); j++)
			if (typesAndCounts[i].involvedIn[j].type != "na")
			{
				std::stringstream newNameNum;
				newNameNum << typesAndCounts[i].involvedIn[j].genNum;
				typesAndCounts[i].involvedIn[j].type =(newNameNum.str() + "!" + typesAndCounts[i].involvedIn[j].type).c_str();
			}


	//fill toWrite vector
	for (int i = 0; i < typesAndCounts.size(); i++)
		for (int j = 0; j < typesAndCounts[i].involvedIn.size(); j++)
			if (typesAndCounts[i].involvedIn[j].type != "na")
				toWrite.push_back(typesAndCounts[i].involvedIn[j].type);
	toWrite.push_back(*new std::string = "END");
	toWrite.push_back(*new std::string = "endf");

	//get volume, since we want to print the density 
	//TODO:  write a way to determine if the material is 2D, in which case you should get the area (since density isnt meaningful for 2D materials
	std::vector <double> vectA, vectB, vectC;

	for (int i = 0; i < 3; i++)
	{
		vectA.push_back(superCellVectorA[i]);
		vectB.push_back(superCellVectorB[i]);
		vectC.push_back(superCellVectorC[i]);
	}

	double a = sqrt((superCellVectorA[0] * superCellVectorA[0]) + (superCellVectorA[1] * superCellVectorA[1]) + (superCellVectorA[2] * superCellVectorA[2]));
	double b = sqrt((superCellVectorB[0] * superCellVectorB[0]) + (superCellVectorB[1] * superCellVectorB[1]) + (superCellVectorB[2] * superCellVectorB[2]));
	double c = sqrt((superCellVectorC[0] * superCellVectorC[0]) + (superCellVectorC[1] * superCellVectorC[1]) + (superCellVectorC[2] * superCellVectorC[2]));
	double alpha = angleBetween(vectB, vectC); //the angle between b and c
	double beta = angleBetween(vectA, vectC); //the angle between a and c
	double gamma = angleBetween(vectA, vectB); //the angle between a and b
	double volume = a*b*c*sqrt(1 - (cos(alpha) * cos(alpha)) - (cos(beta) * cos(beta)) - (cos(gamma) * cos(gamma)) + 2 * cos(alpha)*cos(beta)*cos(gamma));

	//sort the vector for easier writing
	std::sort(toWrite.begin(), toWrite.end(), stringSort);

	//write to outfile csv
	outfile << "name: num!bond type,count,density (1 / angst^3),";
	double count = 0;///keep this as double since c++ dosent like dividing ints by floats
	for (int i = 0; toWrite[i+1] != "endf"; i++)
	{
		count++;
		if (toWrite[i+1] != toWrite[i])
		{
			outfile << toWrite[i] << "," << count << "," << count / volume << ",end,";
			count = 0;
		}
	}

	outfile.close();
}
