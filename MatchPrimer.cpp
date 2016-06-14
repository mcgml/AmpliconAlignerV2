/*
* Filename : MatchPrimer.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Uses gapless alignement of primer to read
* Status: Release
*/

#include <string>

using namespace std;

bool MatchPrimer(const string& Seq, const string& Primer) //iterate over bases of primer and match to seq
{
	float BasesMatched = 0;
	unsigned MaxMismatchLen = 3, PrimerLen = Primer.length(); //no mismatches in the last 3bp -- prevents indels through phase shift and reduced off-target reads

	for (unsigned base = 0; base < PrimerLen; ++base) {

		if (Seq[base] == Primer[base]) {
			BasesMatched++;
		} else if (base > (PrimerLen - MaxMismatchLen)) {
			return 0;
		}

	}

	if (BasesMatched / (PrimerLen - MaxMismatchLen) > 0.8) { //check if match is acceptable
		return 1;
	}

	return 0;
}