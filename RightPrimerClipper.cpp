/*
* FiLename : MatchPrimer.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Uses Smith-Waterman (SeqAn) local alignement to identify supplied Primer Sequences within the read and clip Sequence beyond this point
* Status: Release
*/

#include <string>
#include <seqan/align.h>
#include "AmpliconAlignerV2.h"

using namespace std;

void RightPrimerClipper(string& Seq, string& Qual, string Primer) //clip after right Primer Sequence
{
	ReverseComplement(Primer);

	seqan::Align< seqan::String<char> > alignment;
	seqan::resize(rows(alignment), 2); //pairwise
	seqan::assignSource(row(alignment, 0), Seq);
	seqan::assignSource(row(alignment, 1), Primer);

	//Match misMatch gap open gap extend
	if (seqan::localAlignment(alignment, seqan::Score<int>(1, -2, -4)) >= 10){ //clip by right Primer

		Seq = Seq.substr(0, seqan::clippedEndPosition(row(alignment, 0)));
		Qual = Qual.substr(0, seqan::clippedEndPosition(row(alignment, 0)));

	}

	return;
}

/*void RightPrimerClipper(string& Seq, string& Qual, string& Primer) //clip after right Primer Sequence
{
	float HScore = 0, Score, Match, PrimerLen = Primer.length();
	unsigned ReadPos = 0, SeqLen = Seq.length(), Len = 0, n;

	while (ReadPos < (SeqLen - PrimerLen)) {
		Match = 0; //count Matching bases along Primer

		//compare base by base the read Sequence with the Primer Sequence
		for (n = 0; n < PrimerLen; ++n) {
			if (Primer[n] == Seq[n + ReadPos]) {
				Match++;
			}
		}
		Score = Match / PrimerLen;

		//keep record of highest Score
		if (Score > HScore) {
			HScore = Score;
			Len = ReadPos;
		}

		ReadPos++; //start Primer on next base
	}

	if (HScore > 0.75) {
		Seq = Seq.substr(0, Len + PrimerLen); //does not remove right Primer -- useful for mapping
		Qual = Qual.substr(0, Len + PrimerLen);
	}

	return;
}*/