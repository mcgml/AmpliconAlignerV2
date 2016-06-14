/*
* Filename : GetAmplicons.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Extract fields from supplied amplicon file.
* Status: Release
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "AmpliconAlignerV2.h"

using namespace std;

bool GetAmplicons(ifstream& Amplicons_in, vector<AmpliconRecord>& AmpliconRecords, vector<string>& SamHeaders){

	string AmpliconLine;
	vector<string> AmpliconFields;
	AmpliconRecord TempRecord;

	//Parse primer file
	if (Amplicons_in.is_open()) {
		while (Amplicons_in.good()) {
			getline(Amplicons_in, AmpliconLine);

			boost::trim(AmpliconLine); //remove whitespace at either end of line; MS excel likes putting this in.

			//skip empty lines and headers
			if (AmpliconLine == "" || AmpliconLine[0] == '#') {
				continue;
			} else if (AmpliconLine[0] == '@') {
				SamHeaders.push_back(AmpliconLine);
			} else {

				//tokenize string
				boost::split(AmpliconFields, AmpliconLine, boost::is_any_of("\t"), boost::token_compress_on); //substrings in elements

				if (AmpliconFields.size() != 7) {

					std::cerr << "ERROR: Amplicon list improperly formatted." << endl;
					std::cerr << "AmpliconName Chr Start RefSequence LeftPrimerLength RightPrimerLength Strand(+/-)" << endl;
					return 1;

				} else {

					boost::to_upper(AmpliconFields[3]); //convert sequence to upper-case

					if (isStringDNA(AmpliconFields[3]) == 1) { //check sequence contains only ACTG in upper-case
						std::cerr << "ERROR: " << AmpliconFields[0] << " sequence contains non-standard bases." << endl;
						return 1;
					}

					TempRecord.ID = AmpliconFields[0]; //ampliconname
					TempRecord.Chrom = AmpliconFields[1]; //chromosome
					TempRecord.RefSeq = AmpliconFields[3]; //read 1 sequence

					//left primer
					TempRecord.LeftPrimerLen = boost::lexical_cast<unsigned>(AmpliconFields[4]);
					TempRecord.LeftPrimer = AmpliconFields[3].substr(0, TempRecord.LeftPrimerLen);

					//right primer
					AmpliconFields[3] = ReverseComplement(AmpliconFields[3]);
					TempRecord.RightPrimerLen = boost::lexical_cast<unsigned>(AmpliconFields[5]);
					TempRecord.RightPrimer = AmpliconFields[3].substr(0, TempRecord.RightPrimerLen);

					if (AmpliconFields[6] == "+") { //is+strand?
						TempRecord.Strand = true; //needed to output SAM correctly
						TempRecord.Pos = boost::lexical_cast<unsigned>(AmpliconFields[2]) + TempRecord.LeftPrimerLen; //1-based left coordinate //add primer length to coordinate; after soft-clipping read must be shifted
					} else if (AmpliconFields[6] == "-") { //revcomp needed
						TempRecord.Strand = false; //needed to output SAM correctly
						TempRecord.Pos = boost::lexical_cast<unsigned>(AmpliconFields[2]) + TempRecord.RightPrimerLen;
						TempRecord.RefSeq = ReverseComplement(TempRecord.RefSeq);
					} else {
						std::cerr << "ERROR: " << AmpliconFields[0] << " strand field must be + or -" << endl;
						return 1;
					}

					AmpliconRecords.push_back(TempRecord);
				}

				//delete elements (tokens) for next line
				AmpliconFields.clear();

			}
		}

		Amplicons_in.close();

	} else {
		std::cerr << "ERROR: Unable to open amplicon file" << endl;
		return 1;
	}

	return 0;
}