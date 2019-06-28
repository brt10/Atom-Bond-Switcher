Instructions: 
	1. Edit the file paths at the top of the driver .cpp file to appropriate strings
	2. Edit your bond data file to include info about the atoms that you will be using (included in this version are Si, C, Fe, and H)
	3. Be sure to have, at least, g++ (GCC) 5.3.1 20160406 (Red Hat 5.3.1-6) installed before compiling with 
	   " g++ AtomSwitcherDrive.cpp "
	4. Run with " ./a.out "
	5. Follow instructions as written to console by the code
		There are a couple of atom bond switches that you can not do:
			Switching bond A-A with bond A-B
			Switching bond A-B with bond C-D
		These are not allowed because the computer would not know what atom to switch (the program should stop the user from trying to make these switches, but it might not)


Input:
	1. Correctly formatted VASP input file
	2. Bond data info file (formatted in the same way as the included 'bondDataInfo' file is)
Output:
	1. VASP output file with appropriate atom switches
		Output will probably be written in cartesian coordinates, despite its input formatting
		(to change this, add " POSCARnew.convertToDirect(); " before the line that writes the output in the .cpp driver file)
