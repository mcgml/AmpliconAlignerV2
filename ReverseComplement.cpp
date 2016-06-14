/*
* Filename : ReverseComplement.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Returns the reverse complement of a DNA sequence
* Status: Release
*/

#include <string>

using namespace std;

string ReverseComplement(const string& DNA) {

	string revcomp;

	for (short n = (DNA.length() - 1); n > -1; --n) {

		if (DNA[n] == 'A') {
			revcomp += 'T';
		} else if (DNA[n] == 'T') {
			revcomp += 'A';
		} else if (DNA[n] == 'G') {
			revcomp += 'C';
		} else if (DNA[n] == 'C') {
			revcomp += 'G';
		} else if (DNA[n] == 'a') {
			revcomp += 't';
		} else if (DNA[n] == 't') {
			revcomp += 'a';
		} else if (DNA[n] == 'g') {
			revcomp += 'c';
		} else if (DNA[n] == 'c') {
			revcomp += 'g';
		} else {
			revcomp += DNA[n];
		}

	}

	return revcomp;
}