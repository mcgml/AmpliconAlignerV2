/*
* Filename : isReadNMasked.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Returns true if all bases a N.
* Status: Release
*/

#include <string>
#include "AmpliconAlignerV2.h"

using namespace std;

bool isReadNMasked(const string& read) {

	for (unsigned n = 0; n < read.size(); ++n) {
		if (read[n] != 'N') {
			return false;
		}
	}

	return true;
}
