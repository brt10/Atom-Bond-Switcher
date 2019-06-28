#include <iostream>
#include <ctime>    
#include <cstdlib>  
#include "PoscarInfo.h"
#include "Calculations.h"



std::string infileLoc = "/storage/work/vkb5066/scriptTests/tst13/POSCAR";
std::string bondInfoLoc = "/storage/work/vkb5066/scriptTests/tst13/bondDataInfo";
std::string writeLoc = "/storage/work/vkb5066/scriptTests/tst13/POSCARnew";

bool verifyBonds(Poscar, std::string, std::string);

int main()
{
	//Read in poscar, initialize vector of bonds
	Poscar POSCAR("ReadAll", infileLoc.c_str());
	POSCAR.fetchAtomBonds(bondInfoLoc.c_str());

	srand(time(NULL));

	//Loop through asking user what they want to replace until they want to end program
	for(;;)
	{
		//User input
		std::string A1, A2, B1, B2;
		double percentToChange = 100;
		{
		std::cout << "Replace bond type A with type B:\nType 'EXIT' at any time to leave.\nEnter Bond type A as two strings seperated by a space or a return...\n";
		std::cin >> A1;
		if (A1 == "EXIT")
			break;
		std::cin >> A2;
		if (A2 == "EXIT")
			break;

		std::cout << "Replace with what?\n";
		std::cin >> B1;
		if (B1 == "EXIT")
			break;
		std::cin >> B2;
		if (B2 == "EXIT")
			break;

		std::cout << "What percentage of the " << A1 << "-" << A2 << " bonds should be changed?\n";
		std::cin >> percentToChange;
		if (percentToChange > 100)
			percentToChange = 100;
		if (percentToChange < 0)
			percentToChange = 0;
		std::cout << "Replacing " << percentToChange << "% of bonds\n";
		}

		//Make sure the bond type to replace actually exists in instance of Poscar
		if (!verifyBonds(POSCAR, A1, A2))
		{
			std::cout << "No " << A1 << "-" << A2 << " bonds found in the input file\n";
			break;
		}

		//Send warning if user is trying to replace both existing atoms with two different atoms
		if (((A1 != A2) && (A1 != B1) && (A1 != B2) && (A2 != B1) && (A2 != B2) && (B1 != B2)) || ((A1 == A2) && (B1 != B2)))
		{
			std::cout << "You are attempting to do an illegal operation:  At least one argument in the input leaves me unsure of what atoms to change\n";
			break;
		}

		std::cout << "Replacing " << A1 << "-" << A2 << " bonds with " << B1 << "-" << B2 << " bonds...\n";

		//Create new vector of atomPairs that will eventually replace old atomCoords vector.  Initialize
		std::vector <atomPair> newAtomPairVect;
		for (int i = 0; i < POSCAR.atomBonds.size(); i++)
		{
			int randomNum = rand() % 100;
			if ((A1 == POSCAR.atomBonds[i].pairedAtoms[0].atomType && A2 == POSCAR.atomBonds[i].pairedAtoms[1].atomType) || (A1 == POSCAR.atomBonds[i].pairedAtoms[1].atomType && A2 == POSCAR.atomBonds[i].pairedAtoms[0].atomType))
			{///Consider this bond for editing --- TODO:  there is probably a way to write the below in a more general way; figure it out
				Coords tmp0 = POSCAR.atomBonds[i].pairedAtoms[0];
				Coords tmp1 = POSCAR.atomBonds[i].pairedAtoms[1];
				//case where B1 == B2 (easy)
				if (B1 == B2)
				{
					tmp0.atomType = B1;
					tmp1.atomType = B2;
				}
				else ///B1 must not equal B2, and one of them must match with A1 or A2
				{    ///Replace the A that does not have a matching B with the B that does not have a matching A
					if (B1 != tmp0.atomType && B1 != tmp1.atomType)
					{
						if (tmp0.atomType != B2)
							tmp0.atomType = B1;
						else 
							if (tmp1.atomType != B2)
								tmp1.atomType = B1;
					}
					else
						if (B2 != tmp0.atomType && B2 != tmp1.atomType)
						{
							if (tmp0.atomType != B1)
								tmp0.atomType = B2;
							else
								if (tmp1.atomType != B1)
									tmp1.atomType = B2;
						}
						else
							std::cout << "uh oh\n";
				}
				//Factoring in random amount of switched bonds
				if (randomNum <= percentToChange)
					newAtomPairVect.push_back(*new atomPair(tmp0, tmp1));
				else
					newAtomPairVect.push_back(*new atomPair(POSCAR.atomBonds[i].pairedAtoms[0], POSCAR.atomBonds[i].pairedAtoms[1]));
			}
			else ///else the bond is perfectly fine; keep it
				newAtomPairVect.push_back(POSCAR.atomBonds[i]);
		}

		//Set up new coords:  put all instances of atoms that are not already ID'd into a coords vector
		std::vector <Coords> newCoordsVect;
		for (int i = 0; i < newAtomPairVect.size(); i++)
			for (int n = 0; n < 2; n++)
				newCoordsVect.push_back(newAtomPairVect[i].pairedAtoms[n]);
				
		//Create new instance of Poscar, exactly the same as the input instance, but change its atomCoords vector
		Poscar POSCARnew = POSCAR;
		POSCARnew.atomCoords = newCoordsVect;
		POSCARnew.removeDuplicates();
		POSCARnew.removeDuplicates("unoriginal");
		POSCARnew.updateAtomTypes();
		POSCARnew.updateAtomCounts(); ///update atom type counts, just in case

		//Finially, write the new Poscar to a file
		POSCARnew.write((writeLoc + A1 + A2 + "_to_" + B1 + B2).c_str());

		std::cout << "\n\n";
	}

	return 0;
}





//Checks if bond type a + b (or b + a) exists anywhere in the 'instances' stored bonds
bool verifyBonds (Poscar instance, std::string a, std::string b)
{
	for (int i = 0; i < instance.atomBonds.size(); i++)
		if ((a + b).c_str() == instance.atomBonds[i].type.c_str() || (b + a).c_str() == instance.atomBonds[i].type.c_str())
			return false;
	return true;
}
