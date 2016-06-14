/*
* Filename : ReadMerger.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Merges overlapping paired-end reads using gapless alignment taking the highest quality bases and recalibrating the Scores across the overlap. 
* Status: Release
*/

#include <string>
#include <algorithm>
#include "AmpliconAlignerV2.h"

using namespace std;

bool ReadMerger(const string& SeqR1, const string& QualR1, string SeqR2, string QualR2,
	const unsigned MaxQScore, const unsigned QScorePhredOffset, pair<string, string>& MergedRead) {

	/*									Method
	R1 ---->	R1 ---->		B1 R1 ----> B2 R1 ---->  B3 R1 ---->   B4 R1 ---->
	R2   <----	R2   ----> (RC) B1 R2 ----> B2 R2  ----> B3 R2   ----> B4 R2    ----> etc
	*/

	//Parameters
	unsigned MinScore = 15, MismatchPenalty = 4, MatchAward = 1; //use positive values
	unsigned MisMatchDenominator = 20; //overlap length / MisMatchDenominatorless; than 5% MisMatches

	unsigned ReadPos = 0, SeqR1Len = SeqR1.length(), SeqR2Len = SeqR2.length(), n, BestPos, MisMatches;
	int Score, Q1, Q2, BestScore = 0, SecondBestScore = 0;

	//convert R2 orientation and complement
	SeqR2 = ReverseComplement(SeqR2);
	reverse(QualR2.begin(), QualR2.end());

	//match base by base reads and Score
	while (ReadPos < SeqR1Len) { //iterate over SeqR1
		Score = 0;
		MisMatches = 0;

		//Fix R1 in place, start R1 base 1 at R2 base 1, move R2 left to right one base at a time and check for matches/MisMatches against R1

		for (n = 0; n + ReadPos < SeqR1Len && n < SeqR2Len; ++n) { //stop loop when SeqR2 (ReadPos) extends beyond the length of the SeqR1 OR get to the end of R2

			if (SeqR1[n + ReadPos] == SeqR2[n]) {
				Score += MatchAward;
			} else {
				Score -= MismatchPenalty;
				MisMatches++;
			}

			if (MisMatches > (float)((SeqR1Len - ReadPos) / MisMatchDenominator)) {  //length of potential overlap over maxmismatchdenominator
				break; //stop checking if read exceeds maximum MisMatches for the whole overlap; improves preformance and accuracy
			}

		}

		if (Score > BestScore) {
			SecondBestScore = BestScore;
			BestScore = Score;
			BestPos = ReadPos;
		}

		ReadPos++; //start read on next base
	}

	//check best alignment & merge
	if (BestScore > MinScore && (float) SecondBestScore / BestScore < 0.9 && SeqR2Len + BestPos >= SeqR1Len) {
		//Score is adequate, Score is sufficently higher than the second best & R1 does not have adapter

		//attach start of SeqR1
		MergedRead.first = SeqR1.substr(0, BestPos);
		MergedRead.second = QualR1.substr(0, BestPos);

		//take consensus across overlap
		for (n = 0; n + BestPos < SeqR1Len; ++n) {

			Q1 = QualR1[n + BestPos] - QScorePhredOffset;
			Q2 = QualR2[n] - QScorePhredOffset;

			if (SeqR1[n + BestPos] == SeqR2[n]) { //base is the same; match
				MergedRead.first += SeqR1[n + BestPos];

				if (Q1 + Q2 > MaxQScore) { //add quality Scores together; cap at 40
					MergedRead.second += toascii(MaxQScore + QScorePhredOffset);
				} else {
					MergedRead.second += toascii(Q1 + Q2 + QScorePhredOffset);
				}

			} else { //bases not the same; mismatch

				/*CLC BIO:If the two Scores of the input reads are approximately equal, the resulting Score will be very low which will reflect the fact that it is a very unreliable base.
				On the other hand, if one Score is very low and the other is high, it is likely that the base with the high quality Score is indeed correct,
				and this will be reflected in a relatively high quality Score.*/

				if (Q1 >= Q2){ //R1 is more likely to be correct even if QScores are the same
					MergedRead.first += SeqR1[n + BestPos]; //use highest scoring base
					MergedRead.second += toascii((Q1 - Q2) + QScorePhredOffset); //
				} else {
					MergedRead.first += SeqR2[n]; //use highest scoring base
					MergedRead.second += toascii((Q2 - Q1) + QScorePhredOffset);
				}

			}

		}

		//attach end of SeqR2
		MergedRead.first += SeqR2.substr(SeqR1Len - BestPos, string::npos);
		MergedRead.second += QualR2.substr(SeqR1Len - BestPos, string::npos);

		return true;

	} else {
		return false; //no merge occured
	}

}