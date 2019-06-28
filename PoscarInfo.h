#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>

//Header file for reading in POSCAR and CONTCAR files.  
//As of right now, file formats that include atomic velocities at the bottom are not supported, and will likley cause errors

//Globals:  stuff to change paramaters as needed
double GENERAL_BOND_DISTANCE_ = 2.5; ///(angstroms) distance between any two atoms for them to be considered bonded (works well for POSCARS with only a few atom types)
double CtoC_BOND_DISTANCE_ = 1.6; ///(angstroms) distance between two carbon bonds for them to be considered bonded (normally 1.42 but should be bigger to account for uncertainty)

//Main variable type for this class.  Holds information about individual atoms
struct Coords
{
	std::string atomType; ///holds the type of atom being described
	double a, b, c; ///atomic coordinates 
	std::string flags; ///string for the optional molecular dynamics flags after the coordinates.  Can also hold a blank string to avoid annoying errors
	int id; ///useful way for IDing atoms in the case that an iterative loop is not enough

	std::string extraInfo; ///string for tags that may be useful for other calculations.  Not used in this header file

	//Default constructor
	Coords()
	{
		atomType = "undf"; ///can probably be left alone, but may be useful
		a = -99;
		b = -99;
		c = -99;
		flags = " Uninitialized instance...you should not see this. : Coords - Default Constructor";
		id = -1;
		extraInfo = "undf : Coords: - Default Constructor";
	}

	//Operator declerations
	Coords& operator= (const Coords &coords);
	bool operator== (const Coords &) const;
	bool operator!= (const Coords &) const;
};

double dist_(const Coords, const Coords);

//Operators
Coords& Coords::operator= (const Coords &coords)
{
	atomType = coords.atomType;
	a = coords.a;
	b = coords.b;
	c = coords.c;
	flags = coords.flags;
	id = coords.id;
	extraInfo = coords.extraInfo;
	return *this;
}

bool Coords::operator ==(const Coords &coords) const
{
	return (dist_(*this, coords) < 0.01); ///TODO:  try to find a better way to determine if two atoms are close to eachother
}

bool Coords::operator !=(const Coords & coords) const
{
	return !(*this == coords);
}

//Extra variable type.  Holds information about pairs (not necessairily bonds) of atoms.  Max supported is two atoms per pair
struct atomPair
{
	Coords pairedAtoms[2]; ///info about the two related atoms
	double lenBetween; ///length between the two paired atoms
	double genNum; ///general number for various calculations.  should be reset at the end of every function, since it is used for multiple things
	std::string type; ///type of pair, not type of atom. for example, two atoms A and B should create a type of AB
	std::string extraInfo;

	//Default class constructor.  Not overly useful, but may be needed in some implementations
	atomPair()
	{
		Coords a, b; ///leave these uninitialized -- Poscar class constructor, or a call to the initialization, should set values automatically
		pairedAtoms[0] = a;
		pairedAtoms[1] = b;

		genNum = 0;
		lenBetween = 0;
		type = "uninitialized...you should not see this. : Atompair - Default Constructor";
		extraInfo = "undf : AtomPair - Default Constructor";
	}

	//Paramaterized class constructor.  Uses the info from the two input atoms, ***under the assumption that they are initialized already***
	atomPair(const Coords atomA, const Coords atomB)
	{
		pairedAtoms[0] = atomA;
		pairedAtoms[1] = atomB;

		genNum = 0;
		lenBetween = dist_(atomA, atomB);
		type = (atomA.atomType + atomB.atomType).c_str(); ///should add a sort by alphabetical order so C-H bonds and H-C bonds are considered the same thing
		extraInfo = "undf";
	}

	//Operator declerations
	atomPair& operator= (const atomPair &pair);
	bool operator== (const atomPair &) const;
	bool operator&= (const atomPair &) const;
	bool operator!= (const atomPair &) const;
};

//Operators
atomPair& atomPair::operator= (const atomPair &pair)
{
	pairedAtoms[0] = pair.pairedAtoms[0];
	pairedAtoms[1] = pair.pairedAtoms[1];
	lenBetween = pair.lenBetween;
	genNum = pair.genNum;
	type = pair.type;
	extraInfo = pair.extraInfo;
	return *this;
}

///Compare by Position
bool atomPair::operator ==(const atomPair &atompair) const
{
	return (((pairedAtoms[0] == atompair.pairedAtoms[0]) && (pairedAtoms[1] == atompair.pairedAtoms[1])) || ((pairedAtoms[0] == atompair.pairedAtoms[1]) && (pairedAtoms[1] == atompair.pairedAtoms[0])));
}

///Compare by ID
bool atomPair::operator &=(const atomPair &atompair) const
{
	return (((pairedAtoms[0].id == atompair.pairedAtoms[0].id) && (pairedAtoms[1].id == atompair.pairedAtoms[1].id)) || ((pairedAtoms[0].id == atompair.pairedAtoms[1].id) && (pairedAtoms[1].id == atompair.pairedAtoms[0].id)));
}

bool atomPair::operator !=(const atomPair & atompair) const
{
	return !(*this == atompair);
}

//Comparison rule for sorting strings
bool stringSort_(std::string a, std::string b)
{
	return a < b;
}

//-------------------------------------------------------------------------------------------------------------------

//Class Decleration
class Poscar
{
	public:
		//Variables
		std::string defaultPath; //= DEFAULTOPENPATH.c_str();
		std::string defaultWritePath; //= (DEFAULTOPENPATH + "write").c_str();
		std::string infilePath;
		std::string fileTitle; ///first line of file that is ignored by VASP
		std::string modelType; ///defines if model is bulk or molecular
		double universalScaleFactor; ///all atomic coordinates and lattice vectors are scaled by this.  A negative value is read as a predefined cell volume
		double superCellVectorA [3];  ///-----      
		double superCellVectorB [3];  ///      > defines the supercell
		double superCellVectorC [3];  ///-----
		std::vector <std::string> atomTypes; ///vector to hold arbitrary number of atom types ('C', 'Ag', etc.)
		std::vector <int> atomTypeNums; ///vector (parllel to above) to hold the number of atoms of each type
		bool selectiveDynamicsTag; ///false if the selective dynamics tag is not there
		bool directTag, cartesianTag; ///true if the system is in direct or cartesian coordinates, respectivly
		std::vector <Coords> atomCoords; ///vector of all Coords (defined before the class)
		std::vector <atomPair> atomPairs; ///vector of every atom's relation to every other atom.  get rid of this if the code takes too long to run
		std::vector <atomPair> atomBonds; ///similar to atomPairs, but with a distance restriction on what atoms are related to one another

		//Constructor Declerations
		Poscar();
		Poscar(std::string);
		Poscar(std::string, std::string);

		//Functions Declerations included in this header file
		void fetchFileTitle();
		void fetchModelType();
		void fetchUniversalScaleFactor();
		void fetchSuperCellVectors();
		void fetchAtomTypes();
		void fetchAtomTypeNums();
		void fetchSelectiveDynamicsTag();
		void fetchDirectTag();
		void fetchCartesianTag();
		void fetchAtomCoords();
		void convertToDirect();
		void convertToCartesian();
		void removeTaggedAtoms(std::string);
		void removeDuplicates(); ///note:  NOT automatically called upon construction of class instance.  need to manually call this if you want it to be used
		void removeDuplicates(std::string); ///paramatarized version.  can removes by coords ("coords") or by atom ID ("id") or non-original atoms [for use with supercell extend] ("unoriginal")
		void updateAtomTypes();
		void updateAtomCounts();
		void updateAtomCounts(std::vector <Coords>);
		void print();
		void print(std::string, std::string);
		void write();
		void write(std::string);
		void clearEmptyValues();

		Poscar& operator= (const Poscar &poscar);

		//Function Declerations included in the 'Calculations' header file
		void extendSupercell(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int);
		double grapheneAreaApprox(); ///assuming the POSCAR describes a 4-vertex plane, calculates the surface area of the carbons
		void translateAtoms(double, double, double);
		void applyStrain(bool, bool, bool, int, double); ///strains models
		void fetchAtomPairs(); ///NOT automatically called in class constructor
		void fetchAtomBonds(); ///NOT automatically called in class constructor
		void fetchAtomBonds(std::string);
		void writeBondInfo(std::string);
};

//Class Constructors-----------------------------------------------------------------------
//Default class constructor; user will have to supply most of the info somehow
Poscar::Poscar()
{
	defaultPath = "C:\\Users\\baron\\Desktop\\POSCAR";
	defaultWritePath = (defaultPath + "write").c_str();

	infilePath = defaultPath;
	std::string str = "nothing read in: Poscar - Default Constructor";
	fileTitle = str;
	modelType = "undf";
	universalScaleFactor = 0;
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = -1;
		superCellVectorB[i] = -1;
		superCellVectorC[i] = -1;
	}
	atomTypes.push_back(str);
	atomTypeNums.push_back(-1);
	selectiveDynamicsTag = 0;
	directTag = 0;
	cartesianTag = 0;
	atomCoords;
	atomPairs;
	atomBonds;
}

//Default class constructor; user will have to supply most of the info somehow; but with the filepath defined
Poscar::Poscar(std::string path)
{
	infilePath = path;
	defaultWritePath = (infilePath + "write").c_str();

	std::string str = "nothing read in: Poscar - string Constructor";
	fileTitle = str;
	modelType = "undf";
	universalScaleFactor = 0;
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = -1;
		superCellVectorB[i] = -1;
		superCellVectorC[i] = -1;
	}
	atomTypes.push_back(str);
	atomTypeNums.push_back(-1);
	selectiveDynamicsTag = 0;
	directTag = 0;
	cartesianTag = 0;
	atomCoords;
	atomPairs;
	atomBonds;
}

//Special (useful) class constructor.  Input: (string) instructions["readHead"] or ["readAll"], (string) file path
Poscar::Poscar(std::string instructions, std::string path)
{
	//Reads only the head of the file.  Does NOT return an iterator to the beginning of atomic coordinates; user will have to do this before reading those in
	if (instructions == "ReadHead" || instructions == "readHead")
	{
		infilePath = path;
		defaultWritePath = (infilePath + "write").c_str();

		fetchFileTitle();
		fetchUniversalScaleFactor();
		fetchSuperCellVectors();
		fetchAtomTypes();
		fetchModelType();
		fetchAtomTypeNums();
		fetchSelectiveDynamicsTag();
		fetchDirectTag();
		fetchCartesianTag();
	}

	//Reads the entire file, storing atomic coordinates in a vector
	if (instructions == "ReadAll" || instructions == "readAll")
	{
		infilePath = path;
		defaultWritePath = (infilePath + "write").c_str();

		fetchFileTitle();
		fetchUniversalScaleFactor();
		fetchSuperCellVectors();
		fetchAtomTypes();
		fetchModelType();
		fetchAtomTypeNums();
		fetchSelectiveDynamicsTag();
		fetchDirectTag();
		fetchCartesianTag();
		fetchAtomCoords();
		atomPairs;
		atomBonds;
	}
}


//Destructor (not needed unless dynamic memory allocation is used.  As of right now, it isnt.)
//~Poscar() {}

//-----------------------------------------------------------------------------------------
//Overloaded = sign.  This "should" be very useful
Poscar& Poscar::operator= (const Poscar &poscar)
{
	infilePath = poscar.infilePath;
	defaultWritePath = poscar.defaultWritePath;
	fileTitle = poscar.fileTitle;
	universalScaleFactor = poscar.universalScaleFactor;
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = poscar.superCellVectorA[i];
		superCellVectorB[i] = poscar.superCellVectorB[i];
		superCellVectorC[i] = poscar.superCellVectorC[i];
	}
	atomTypes = poscar.atomTypes;
	atomTypeNums = poscar.atomTypeNums;
	modelType = poscar.modelType;
	selectiveDynamicsTag = poscar.selectiveDynamicsTag;
	directTag = poscar.directTag;
	cartesianTag = poscar.cartesianTag;
	atomCoords = poscar.atomCoords;
	atomPairs = poscar.atomPairs;
	atomBonds = poscar.atomBonds;
	return *this;
}

//----------------------------------------------------------------------------------------------------------------------

//Non-Member Function Definitions----------------------------------------------------------------------------------------

std::string readNthLine(const std::string filename, int n)
{
		std::ifstream infile(filename.c_str());
		std::string s;

		//Make this faster
		s.reserve(1024);

		//skip n lines
		for (int i = 0; i < n; ++i)
			std::getline(infile, s);

		std::getline(infile, s);
		return s;
}

double angleBetween (std::vector <double> a, std::vector <double> b)
{
	return acos(((a[0]*b[0]) + (a[1]*b[1]) + (a[2] * b[2]))/(sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]))*sqrt((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]))));
}

//Distance between two points in R3.  Why can I not use the one defined in Calculations.h?  Even though i'm trying to include it?  I dont know
double dist_(const Coords coordA, const Coords coordB)
{
	return sqrt((((coordB.a - coordA.a)*(coordB.a - coordA.a))) + (((coordB.b - coordA.b)*(coordB.b - coordA.b))) + (((coordB.c - coordA.c)*(coordB.c - coordA.c))));
}

//Formula shamelessly taken from wikipedia: https://en.wikipedia.org/wiki/Fractional_coordinates
void setTransformMatrix (std::vector <double> vectA, std::vector <double> vectB, std::vector <double> vectC, std::vector <double> &toReturn)
{
	toReturn.clear();
	//syntax for vector positions in reference to the matrix:
	//        | [0]  [1]  [2] |
	//        | [3]  [4]  [5] |
	//        | [6]  [7]  [8] |

	//Definitions (if you are debugging this: 1) i'm sorry.  2) you should look at the wiki link for the matrix formula used there)
	double a = sqrt((vectA[0] * vectA[0]) + (vectA[1] * vectA[1]) + (vectA[2] * vectA[2]));
	double b = sqrt((vectB[0] * vectB[0]) + (vectB[1] * vectB[1]) + (vectB[2] * vectB[2]));
	double c = sqrt((vectC[0] * vectC[0]) + (vectC[1] * vectC[1]) + (vectC[2] * vectC[2]));
	double alpha = angleBetween(vectB, vectC); //the angle between b and c
	double beta = angleBetween(vectA, vectC); //the angle between a and c
	double gamma = angleBetween(vectA, vectB); //the angle between a and b
	double volume = a*b*c*sqrt(1 - (cos(alpha) * cos(alpha)) - (cos(beta) * cos(beta)) - (cos(gamma) * cos(gamma)) + 2*cos(alpha)*cos(beta)*cos(gamma));
	
	//Filling vector / matrix
	toReturn.push_back(1 / a); //[0]
	toReturn.push_back(-cos(gamma) / (a * sin(gamma))); //[1]
	toReturn.push_back(b * c * ((cos(alpha) * cos(gamma) - cos(beta))/(volume * sin(gamma)))); //[2]
	toReturn.push_back(0); //[3]
	toReturn.push_back(1 / (b * sin(gamma))); //[4]
	toReturn.push_back(a * c * ((cos(beta) * cos(gamma) - cos(alpha)) / (volume * sin(gamma)))); //[5]
	toReturn.push_back(0); //[6]
	toReturn.push_back(0); //[7]
	toReturn.push_back(a * b * sin(gamma)/(volume)); //[8]
}

//(Almost) shamelessly taken from: https://stackoverflow.com/questions/8425214/splitting-string-into-a-vectorstring-of-words
std::vector<std::string> split(const std::string& s)
{
	std::vector<std::string> ret;
	typedef std::string::size_type string_size;
	string_size i = 0;

	// invariant: we have processed characters [original value of i, i) 
	while (i != s.size()) {
		// ignore leading blanks
		// invariant: characters in range [original i, current i) are all spaces
		while (i != s.size() && isspace(s[i]))
			++i;

		// find end of next word
		string_size j = i;
		// invariant: none of the characters in range [original j, current j)is a space
		while (j != s.size() && !isspace(s[j]))
			j++;

		// if we found some nonwhitespace characters 
		if (i != j) {
			// copy from s starting at i and taking j - i chars
			ret.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return ret;
}

bool empty(const std::string s)
{
	std::string s_  = s.c_str(); ///c_str() because unix dosent like "normal" strings
	if (s_ == "" || s_ == " " || s_ == "\t" || s_.length() == 0)
		return 1;
	else return 0;
}


//Member Function Definitions:-------------------------------------------------------------------------------------------------------------------

void Poscar::fetchFileTitle()
{
	//Open file:
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";

	//Get first line, update the class title
	getline(infile, fileTitle);
	fileTitle = fileTitle.c_str(); ///c_str() because unix dosent like "normal" strings
	
	infile.close();
}

//Warning: works based off of whether the input file has carbons (molecular model) or no carbons (bulk model).  Don't rely on this to work in general.
void Poscar::fetchModelType()
{
	if (atomTypes.size() <= 1)
		modelType = "Bulk";
	else
		modelType = "Molecular";
}

void Poscar::fetchUniversalScaleFactor()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";
	
	//move to line before line to read in
	int n = 1; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for(;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update class param
	infile >> universalScaleFactor;

	infile.close();
}

void Poscar::fetchSuperCellVectors()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 2; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update class params
	for (int i = 0; i < 3; i++)
		infile >> superCellVectorA[i];
	for (int i = 0; i < 3; i++)
		infile >> superCellVectorB[i];
	for (int i = 0; i < 3; i++)
		infile >> superCellVectorC[i];

	infile.close();
}

void Poscar::fetchAtomTypes()
{
	if (!atomTypes.empty())
		return;

	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchAtomTypes(): File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 5; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update params:
		//Get the non-blank line that holds all of the atom types
	std::string atomTypeLine;
	for (;;)
	{
		getline(infile, atomTypeLine);
		if (!empty(atomTypeLine.c_str()))
			break;
	}
	atomTypes = split(atomTypeLine);

	infile.close();
}

void Poscar::fetchAtomTypeNums()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 6; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update params
	for (int i = 0; i < atomTypes.size(); i++)
	{	
		int n;
		infile >> n; ///should be in line with the atom types
		atomTypeNums.push_back(n);
	}

	infile.close();
}

void Poscar::fetchSelectiveDynamicsTag()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 7; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update param
	std::string s;
	infile >> s;
	if ((s[0] == 'S' || s[0] == 's') && (s[1] == 'E' || s[1] == 'e'))
		selectiveDynamicsTag = 1;
	else selectiveDynamicsTag = 0;

	infile.close();
}

void Poscar::fetchDirectTag()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";

	//need to skip different number of lines depending on whether the 'Sel' tag is there or not, since it takes up an entire line
	int n;
	if (selectiveDynamicsTag)
		 n = 8;
	else
		 n = 7;
	
	//move to line before line to read in
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update param
	std::string s;
	infile >> s;
	if ((s[0] == 'D' || s[0] == 'd') && (s[1] == 'i' || s[1] == 'i'))
		directTag = 1;
	else directTag = 0;

	infile.close();
}

void Poscar::fetchCartesianTag()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";

	//need to skip different number of lines depending on whether the 'Sel' tag is there or not, since it takes up an entire line
	int n;
	if (selectiveDynamicsTag)
		n = 8;
	else
		n = 7;

	//move to line before line to read in
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update param
	std::string s;
	infile >> s;
	if ((s[0] == 'C' || s[0] == 'c') && (s[1] == 'A' || s[1] == 'a'))
		cartesianTag = 1;
	else cartesianTag = 0;

	infile.close();
}

void Poscar::fetchAtomCoords()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "File could not be found...expect errors\n";

	//need to skip different number of lines depending on whether the 'Sel' tag is there or not, since it takes up an entire line
	int n;
	if (selectiveDynamicsTag)
		n = 9;
	else
		n = 8;

	//move to line before line to read in
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Over-complicated reading in of atomic coordinates so that each atom can be assigned its correct atom type
	for (int i = 0; i < atomTypes.size(); i++)
	{
		for (int j = atomTypeNums[i]; j > 0; j--) ///if there will be any error in this code, it will be here, with this for loop's termination rule
		{
			Coords crds;
			infile >> crds.a >> crds.b >> crds.c;
			getline(infile, crds.flags);
			crds.atomType = atomTypes[i];
			atomCoords.push_back(crds);
		}
	}

	//And setting ID numbers for the atoms
	for (int i = 0; i < atomCoords.size(); i++)
		atomCoords[i].id = i;

	infile.close();
}

void Poscar::convertToDirect()
{
	if (directTag)  ///leaving this here just in case you'd want to do something special if file is already in cartesian
	{
		return;
	} 
	else
		if (cartesianTag)
		{
			std::vector <double> transformMatrix;
			//Somehow managed to set the supercell vectors as arrays, which is annoying to work with.  Put them in vectors:
			//May be worth fixing this later to improve performance...for now, i just have the vectors in a limited scope

			//Fill transform matrix; syntax is in the above function comment
			{	
				std::vector <double> vA, vB, vC;
				for (int i = 0; i < 3; i++)
				{
					vA.push_back(universalScaleFactor * superCellVectorA[i]);
					vB.push_back(universalScaleFactor * superCellVectorB[i]);
					vC.push_back(universalScaleFactor * superCellVectorC[i]);
				}
				setTransformMatrix(vA, vB, vC, transformMatrix);
			}
			//Conversion math
			for (int i = 0; i < atomCoords.size(); i++)
			{
				atomCoords[i].a = universalScaleFactor * ((transformMatrix[0] * atomCoords[i].a) + (transformMatrix[1] * atomCoords[i].b) + (transformMatrix[2] * atomCoords[i].c));
				atomCoords[i].b = universalScaleFactor * ((transformMatrix[3] * atomCoords[i].a) + (transformMatrix[4] * atomCoords[i].b) + (transformMatrix[5] * atomCoords[i].c));
				atomCoords[i].c = universalScaleFactor * ((transformMatrix[6] * atomCoords[i].a) + (transformMatrix[7] * atomCoords[i].b) + (transformMatrix[8] * atomCoords[i].c));
			}

			//Update tags
			for (int i = 0; i < 3; i++)
			{
				superCellVectorA[i] = universalScaleFactor * superCellVectorA[i];
				superCellVectorB[i] = universalScaleFactor * superCellVectorB[i];
				superCellVectorC[i] = universalScaleFactor * superCellVectorC[i];
			}
			universalScaleFactor = 1.0;
			directTag = 1;
			cartesianTag = 0;
		}
}

//Converts all atoms to cartesian.  Could do just one atom, but I cant see a reason that you'd ever want to mix the two coordinate systems in a single file
void Poscar::convertToCartesian()
{
	if (cartesianTag) ///leaving this here just in case you'd want to do something special if file is already in cartesian
	{
		return;
	}
	else 
		if (directTag)
		{
		///conversion math
			for (int i = 0; i < atomCoords.size(); i++)
			{
				atomCoords[i].a = universalScaleFactor * ((superCellVectorA[0] * atomCoords[i].a) + (superCellVectorB[0] * atomCoords[i].b) + (superCellVectorC[0] * atomCoords[i].c));	
				atomCoords[i].b = universalScaleFactor * ((superCellVectorA[1] * atomCoords[i].a) + (superCellVectorB[1] * atomCoords[i].b) + (superCellVectorC[1] * atomCoords[i].c));
				atomCoords[i].c = universalScaleFactor * ((superCellVectorA[2] * atomCoords[i].a) + (superCellVectorB[2] * atomCoords[i].b) + (superCellVectorC[2] * atomCoords[i].c));
			}

			//Updating tags
			for (int i = 0; i < 3; i++)
			{
				superCellVectorA[i] = universalScaleFactor * superCellVectorA[i];
				superCellVectorB[i] = universalScaleFactor * superCellVectorB[i];
				superCellVectorC[i] = universalScaleFactor * superCellVectorC[i];
			}
			universalScaleFactor = 1.0;
			directTag = 0;
			cartesianTag = 1;
		}
}

//Removes atoms that have the same R3 coordinates.  does NOT take any other parts of the Coords class into account (id, extraInfo, etc.)
void Poscar::removeTaggedAtoms(std::string tag)
{
	std::vector <Coords> tmp;

	for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo == tag)
					tmp.push_back(atomCoords[i]);

	atomCoords = tmp;
	updateAtomCounts(tmp);
}

//Removes atoms that have the same R3 coordinates.  does NOT take any other parts of the Coords class into account (id, extraInfo, etc.)
void Poscar::removeDuplicates()
{
	std::vector <Coords> tmp;

	for (int i = 0; i < atomCoords.size(); i++)
		for (int j = 0; j < atomCoords.size(); j++)
			if (i != j)
				if (atomCoords[i] == atomCoords[j] && atomCoords[i].extraInfo != "cp")
					atomCoords[j].extraInfo = "cp";
				
	for (int i = 0; i < atomCoords.size(); i++)
		if (atomCoords[i].extraInfo != "cp")
			tmp.push_back(atomCoords[i]);

	atomCoords = tmp;
	updateAtomCounts(tmp);
}

void Poscar::removeDuplicates(std::string instruct)
{
	std::vector <Coords> tmp;

	if (instruct == "coords")
	{
		for (int i = 0; i < atomCoords.size(); i++)
			for (int j = 0; j < atomCoords.size(); j++)
				if (i != j)
					if (atomCoords[i] == atomCoords[j] && atomCoords[i].extraInfo != "cp")
						atomCoords[j].extraInfo = "cp";

		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "cp")
				tmp.push_back(atomCoords[i]);
	}

	if (instruct == "id")
	{
		for (int i = 0; i < atomCoords.size(); i++)
			for (int j = 0; j < atomCoords.size(); j++)
				if (i != j)
					if (atomCoords[i].id == atomCoords[j].id && atomCoords[i].extraInfo != "cp")
						atomCoords[j].extraInfo = "cp";

		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "cp")
				tmp.push_back(atomCoords[i]);
	}

	if (instruct == "unoriginal")
	{
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "original")
				atomCoords[i].extraInfo = "cp";

		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "cp")
				tmp.push_back(atomCoords[i]);
	}

	atomCoords = tmp;
	updateAtomCounts(tmp);
}

void Poscar::updateAtomTypes()
{
	bool ok = true;
	std::vector <std::string> tmp, tmp2;

	for (int i = 0; i < atomCoords.size(); i++)
		if (atomCoords[i].atomType != "undf" && atomCoords[i].atomType != "END" && (atomCoords[i].atomType[0] != 'u' && atomCoords[i].atomType[1] != 'n'))
			tmp.push_back(atomCoords[i].atomType);

	std::sort(tmp.begin(), tmp.end(), stringSort_);
	tmp.push_back("END");

	for (int i = 0; i < tmp.size(); i++)
		if (tmp[i] != "END")
			if (tmp[i] != tmp[i+1])
				tmp2.push_back(tmp[i]);

	atomTypes = tmp2;
}

void Poscar::updateAtomCounts()
{
	std::vector <int> tmp;
	for (int i = 0; i < atomTypes.size(); i++)
		tmp.push_back(0);


	for (int i = 0; i < tmp.size(); i++)
		for (int j = 0; j < atomCoords.size(); j++)
			if (atomCoords[j].atomType == atomTypes[i])
				tmp[i] ++;

	atomTypeNums = tmp;
}

//TEMPORARY in case this small change breaks everything:  the original: vv The new: ^^
/*
void Poscar::updateAtomCounts()
{
	//Update params
	for (int i = 0; i < atomTypeNums.size(); i++)
		atomTypeNums[i] = 0;

	for (int i = 0; i < atomTypeNums.size(); i++)
		for (int j = 0; j < atomCoords.size(); j++)
			if (atomCoords[j].atomType == atomTypes[i])
				atomTypeNums[i] ++;
}
*/

void Poscar::updateAtomCounts(std::vector <Coords> crds)
{
	//Update params
	for (int i = 0; i < atomTypeNums.size(); i++)
		atomTypeNums[i] = 0;

	for (int i = 0; i < atomTypeNums.size(); i++)
		for (int j = 0; j < crds.size(); j++)
			if (crds[j].atomType == atomTypes[i])
				atomTypeNums[i] ++;
}

void Poscar::print()
{
	//Print original input file
	using namespace std;
	cout << "Printing input file..." << endl << endl << endl;

	//Open file
	ifstream infile(infilePath.c_str());
	if (infile.fail())
		cout << "File could not be found...expect errors\n";

	//Print line by line
	string s;
	while (!infile.eof())
	{
		getline(infile, s);
		cout << s << endl;
	}

	infile.close();

	//Print class contents (mostly for debugging)
	cout << endl << endl << "Printing class interpretation of input file..." << endl << endl << endl;
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);

	cout << fileTitle << endl;
	cout << universalScaleFactor << endl;
	cout << left << setfill('0') << setprecision(10) << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
	cout << left << setfill('0') << setprecision(10) << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
	cout << left << setfill('0') << setprecision(10) << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
	for (int i = 0; i < atomTypes.size(); i++)
		cout << atomTypes[i] << "   ";
	cout << endl;
	for (int i = 0; i < atomTypeNums.size(); i++)
		cout << atomTypeNums[i] << "   ";
	cout << endl;
	if (selectiveDynamicsTag)
		cout << "Sel" << endl;
	if (directTag)
		cout << "Direct   ";
	if (cartesianTag)
		cout << "Cartesian   ";
	cout << endl;
	for (int j = 0; j < atomTypes.size(); j++)
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].atomType == atomTypes[j])
				cout << left << setfill('0') << setprecision(10) << atomCoords[i].a << "\t" << atomCoords[i].b << "\t" << atomCoords[i].c << "\t" << atomCoords[i].flags << endl;
	cout << "FILE END.  THERE SHOULD BE NO SPACE BETWEEN THIS MESSAGE AND THE FINAL LINE OF THE FILE." << endl;
}

void Poscar::print(std::string headOrFull, std::string normalOrExtra)
{
	using namespace std;
	//Always print head of file
	cout << endl << endl << "Printing class interpretation of input file..." << endl << endl << endl;
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	
	if (normalOrExtra != "extra" && normalOrExtra != "Extra")
	{
		cout << fileTitle << endl;
		cout << universalScaleFactor << endl;
		cout << left << setfill('0') << setprecision(10) << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
		cout << left << setfill('0') << setprecision(10) << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
		cout << left << setfill('0') << setprecision(10) << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
		for (int i = 0; i < atomTypes.size(); i++)
			cout << atomTypes[i] << "   ";
		cout << endl;
		for (int i = 0; i < atomTypeNums.size(); i++)
			cout << atomTypeNums[i] << "   ";
		cout << endl;
		if (selectiveDynamicsTag)
			cout << "Sel" << endl;
		if (directTag)
			cout << "Direct   ";
		if (cartesianTag)
			cout << "Cartesian   ";
		cout << endl;
	}

	int linenum = 1;
	if (normalOrExtra == "extra" || normalOrExtra == "Extra")
	{
		cout << "Linenum " << linenum << "\t" << fileTitle << endl;
		linenum ++;
		cout << "Linenum " << linenum << "\t" << universalScaleFactor << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t" << left << setfill('0') << setprecision(10) << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t" << left << setfill('0') << setprecision(10) << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t" << left << setfill('0') << setprecision(10) << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t";
		linenum++;
		for (int i = 0; i < atomTypes.size(); i++)
			cout << atomTypes[i] << "   ";
		cout << endl;
		cout << "Linenum " << linenum << "\t";
		for (int i = 0; i < atomTypeNums.size(); i++)
			cout << atomTypeNums[i] << "   ";
		cout << endl;
		if (selectiveDynamicsTag)
		{
			linenum++;
			cout << "Linenum " << linenum << "\tSel" << endl;
		}
		linenum++;
		cout << "Linenum " << linenum << "\t";
		if (directTag)
			cout << "Direct   ";
		if (cartesianTag)
			cout << "Cartesian   ";
		cout << endl;
	}

	//if user wants the rest of file, print that.  check for extra options as well
	if (headOrFull == "Full" || headOrFull == "full")
	{
		for (int j = 0; j < atomTypes.size(); j++)
		{
			for (int i = 0; i < atomCoords.size(); i++)
				if (atomCoords[i].atomType == atomTypes[j])
				{
					if (normalOrExtra == "extra" || normalOrExtra == "Extra")
					{
						linenum++;
						cout << "Linenum " << linenum << "\tAtomnum" << i + 1 << "\t";
					}
					cout << left << setfill('0') << setprecision(10) << atomCoords[i].a << "\t" << atomCoords[i].b << "\t" << atomCoords[i].c << "\t" << atomCoords[i].flags;
					if (normalOrExtra == "extra" || normalOrExtra == "Extra")
						cout << "\tAtom Type: " << atomCoords[i].atomType << "   extra:  " << atomCoords[i].extraInfo << endl;
					cout << endl;
				}
		}
		cout << "FILE END.  THERE SHOULD BE NO SPACE BETWEEN THIS MESSAGE AND THE FINAL LINE OF THE FILE." << endl;
	}

}

void Poscar::write()
{
	//Open file
	using namespace std;
	ofstream outfile(defaultWritePath.c_str());


	//Write file.  literally the exact same as print, except 'cout' has been replaced with 'outfile'
	outfile << fileTitle << endl;

	outfile << fixed;
	outfile << setprecision(3);

	outfile << left << universalScaleFactor << endl;
	outfile << setprecision(10);
	outfile << left << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
	outfile << left << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
	outfile << left << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
	for (int i = 0; i < atomTypes.size(); i++)
		outfile << atomTypes[i] << "   ";
	outfile << endl;
	for (int i = 0; i < atomTypeNums.size(); i++)
		outfile << atomTypeNums[i] << "   ";
	outfile << endl;
	if (selectiveDynamicsTag)
		outfile << "Sel" << endl;
	if (directTag)
		outfile << "Direct   ";
	if (cartesianTag)
		outfile << "Cartesian   ";
	outfile << endl;
	//Print atoms in order, just in case they got mixed up somehow
	for (int j = 0; j < atomTypes.size(); j++)
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].atomType == atomTypes[j])
				outfile << left << setfill('0') << setprecision(10) << atomCoords[i].a << "\t" << atomCoords[i].b << "\t" << atomCoords[i].c << "\t" << atomCoords[i].flags << endl;
	
	outfile.close();
}

void Poscar::write(std::string userPath)
{
	//Open file
	using namespace std;
	ofstream outfile(userPath.c_str());

	//Write file.  literally the exact same as print, except 'cout' has been replaced with 'outfile'
	outfile << fileTitle << endl;

	outfile << fixed;
	outfile << setprecision(3);

	outfile << left << universalScaleFactor << endl;
	outfile << setprecision(10);
	outfile << left << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
	outfile << left << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
	outfile << left << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
	for (int i = 0; i < atomTypes.size(); i++)
		outfile << atomTypes[i] << "   ";
	outfile << endl;
	for (int i = 0; i < atomTypeNums.size(); i++)
		outfile << atomTypeNums[i] << "   ";
	outfile << endl;
	if (selectiveDynamicsTag)
		outfile << "Sel" << endl;
	if (directTag)
		outfile << "Direct   ";
	if (cartesianTag)
		outfile << "Cartesian   ";
	outfile << endl;
	//Print atoms in order, just in case they got mixed up somehow
	for (int j = 0; j < atomTypes.size(); j++)
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].atomType == atomTypes[j])
				outfile << left << setfill('0') << setprecision(10) << atomCoords[i].a << "\t" << atomCoords[i].b << "\t" << atomCoords[i].c << "\t" << atomCoords[i].flags << endl;

	outfile.close();
}

void Poscar::clearEmptyValues()
{
	for (int i = 0; i < atomTypes.size(); i++)
		if (atomTypes[i] == "" || atomTypes[i] == " ")
			atomTypes[i] = "del";
	std::vector <std::string> tmp;

	for (int i = 0; i < atomTypes.size(); i++)
		if (atomTypes[i] != "del")
			tmp.push_back(atomTypes[i].c_str());

	atomTypes = tmp;
}
