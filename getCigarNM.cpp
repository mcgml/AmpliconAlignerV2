/*
* Filename : getCigarNM.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Calculates read alignment cigar and edit distance.
* Status: Release
*/

#include <string>
#include <vector>
#include <seqan/Align.h>
#include <boost/lexical_cast.hpp>
#include "AmpliconAlignerV2.h"

using namespace std;

bool getCigarNM(TRow& row1, TRow& row2, const unsigned LeftPrimerLengthStrandConverted, const unsigned RightPrimerLengthStrandConverted, pair<string, unsigned>& CigarNM, unsigned& SingleBaseMisMatchFrequency) {
	
	TRowIterator Rit = seqan::begin(row1), RitEnd = seqan::end(row1), Qit = seqan::begin(row2), QitEnd = seqan::end(row2);
	unsigned j = 0;
	vector<pair<char, unsigned>> CigarPairs;
	pair<char, unsigned> Cigartemp;
	CigarNM.first = "";
	CigarNM.second = 0;
	SingleBaseMisMatchFrequency = 0;

	//iterate over Query Alignment
	for (; Qit != QitEnd; ++Qit, ++Rit) {
		j++; //1-based

		if (isGap(Qit)) { //Query deleteion

			if (Qit == seqan::begin(row2)) {
				Cigartemp.first = 'D';
				Cigartemp.second = 1;
			} else if (Cigartemp.first != 'D') {
				CigarPairs.push_back(Cigartemp);
				Cigartemp.first = 'D';
				Cigartemp.second = 1;
			} else { //contig base
				Cigartemp.second++;
			}

		} else if (isGap(Rit)) { //Query insertion

			if (Qit == seqan::begin(row2)) {
				Cigartemp.first = 'I';
				Cigartemp.second = 1;
			} else if (Cigartemp.first != 'I') {
				CigarPairs.push_back(Cigartemp);
				Cigartemp.first = 'I';
				Cigartemp.second = 1;
			} else { //contig base
				Cigartemp.second++;
			}

		} else { //Query match/mismatch

			if (Qit == seqan::begin(row2)) {
				Cigartemp.first = 'M';
				Cigartemp.second = 1;
			} else if (Cigartemp.first != 'M') {
				CigarPairs.push_back(Cigartemp);
				Cigartemp.first = 'M';
				Cigartemp.second = 1;
			} else { //contig base
				Cigartemp.second++;
			}

			if (value(Qit) != value(Rit) && j > LeftPrimerLengthStrandConverted && j <= (QitEnd - seqan::begin(row2)) - RightPrimerLengthStrandConverted) { //base mismatch not in primer
				CigarNM.second++; //edit distance
				SingleBaseMisMatchFrequency++;
			}

		}

	} //end iterating over Alignemnt

	//load last Cigar
	CigarPairs.push_back(Cigartemp);

	//softclip primers
	if (CigarPairs.size() == 1) { //all bases match/mismatch

		if (CigarPairs[0].first == 'M' && CigarPairs[0].second > (LeftPrimerLengthStrandConverted + RightPrimerLengthStrandConverted)) { //matching bases exceed primer length
			CigarNM.first = boost::lexical_cast<string>(LeftPrimerLengthStrandConverted);
			CigarNM.first += 'S';
			CigarNM.first += boost::lexical_cast<string>(CigarPairs[0].second - (LeftPrimerLengthStrandConverted + RightPrimerLengthStrandConverted));
			CigarNM.first += 'M';
			CigarNM.first += boost::lexical_cast<string>(RightPrimerLengthStrandConverted);
			CigarNM.first += 'S';
		} else {
			return 1; //just primer; skipped
		}

	} else if (CigarPairs.size() > 2) {

		if (CigarPairs[0].first == 'M' && CigarPairs[0].second >= LeftPrimerLengthStrandConverted && CigarPairs[CigarPairs.size() - 1].first == 'M' && CigarPairs[CigarPairs.size() - 1].second >= RightPrimerLengthStrandConverted) { //matching bases exceed primer length
			
			CigarNM.first = boost::lexical_cast<string>(LeftPrimerLengthStrandConverted);
			CigarNM.first += 'S';

			if (CigarPairs[0].second > LeftPrimerLengthStrandConverted) {
				CigarNM.first += boost::lexical_cast<string>(CigarPairs[0].second - LeftPrimerLengthStrandConverted);
				CigarNM.first += 'M';
			}

			for (j = 1; j < CigarPairs.size() - 1; j++) { //iterate over Cigars ignore first and last
				CigarNM.first += boost::lexical_cast<string>(CigarPairs[j].second);
				CigarNM.first += CigarPairs[j].first;

				//region is del or ins add to edit distance
				if (CigarPairs[j].first != 'M') {
					CigarNM.second += CigarPairs[j].second;
				}

			}

			if (CigarPairs[CigarPairs.size() - 1].second > RightPrimerLengthStrandConverted) {
				CigarNM.first += boost::lexical_cast<string>(CigarPairs[CigarPairs.size() - 1].second - RightPrimerLengthStrandConverted);
				CigarNM.first += 'M';
			}

			CigarNM.first += boost::lexical_cast<string>(RightPrimerLengthStrandConverted);
			CigarNM.first += 'S';

		} else {
			return 1;
		}

	} else {
		return 1; //0 or 2 skipped
	}

	return 0;
}