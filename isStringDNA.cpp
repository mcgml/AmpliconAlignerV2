/*
* Filename : isStringDNA.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Checks supplied string is DNA.
* Status: Release
*/

#include <string>

using namespace std;

bool isStringDNA(const string& str) { //check DNA sequence input

	for (unsigned short n = 0; n < str.length(); n++) {

		if (str[n] != 'A' && str[n] != 'T' && str[n] != 'G' && str[n] != 'C') {
			return 1; //ERROR: contains non-standard base
		}

	}

	return 0;
}