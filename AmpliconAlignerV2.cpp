/*
* Filename : AmpliconAlignerV2.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Globally aligns paired-end reads to amplicon sequences and outputs in SAM format.
* Status: Release
*/

/*
TODO: need to iterate over all possilble Reference sequences including off-target
TODO: calculate mapping quality score: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577856/
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <seqan/Align.h>
#include "AmpliconAlignerV2.h"

using namespace std;

int main(int argc, char* argv[]) {

	float Version = 2.1;

	//check argument number is correct; print usage
	if (argc != 5) { //program ampliconlist r1 r2 prefix
		std::cerr << "\nProgram: AmpliconAligner v" << Version << endl;
		std::cerr << "Contact: Matthew Lyon, Wessex Regional Genetics Lab (matthew.lyon@salisbury.nhs.uk)\n" << endl;
		std::cerr << "Usage: AmpliconAligner <AmpliconList> <Read1.fastq.gz> <Read2.fastq.gz> <OutputFilenamePrefix>\n" << endl;
		std::cerr << "AmpliconID Chr Start RefSequence LeftPrimerLength RightPrimerLength Strand(+/-)\n" << endl;
		return -1;
	}

	//parameters
	const unsigned minIsize = 5;
	const unsigned MaxQScore = 40;
	const unsigned QScorePhredOffset = 33;
	const float MaxSingleBaseMisMatch = 0.05; //maximum fraction of mismatching bases relative to the wildtype length
	
	unsigned LineNo = 0, Start, TotalReads = 0, nMaskedReads = 0, PrimerMatchedReads = 0, TotalUsableReads = 0, n, 
		LeftPrimerLengthStrandConverted, RightPrimerLengthStrandConverted, TotalMappedReads = 0, TotalNotMergedReads = 0, SingleBaseMisMatchFrequency;
	string Read1Line, Read2Line, Header, Header1, Seq1, Qual1, Header2, Seq2, Qual2, Index, FlowCellID, R1FASTQ = argv[2], R2FASTQ = argv[3];
	pair<string, string> MergedRead;
	pair<string, unsigned> CigarNM;
	vector<string> SamHeaders;
	vector<AmpliconRecord> AmpliconRecords;
	unordered_map <string, Stat> Stats;
	TSequence Ref, Query;
	TAlign Align;
	TRowIterator Rit, RitEnd, Qit, QitEnd;
	int NWScore;
	seqan::resize(rows(Align), 2); //always pairwise
	boost::iostreams::filtering_stream<boost::iostreams::input> R1FilterStream, R2FilterStream;

	//is FASTQ input gziped?
	if (R1FASTQ.substr(R1FASTQ.find_last_of('.'), string::npos) != ".gz" || R2FASTQ.substr(R2FASTQ.find_last_of('.'), string::npos) != ".gz" ) { //?FASTQ is gzipped
		cerr << "ERROR: FASTQ files must be unmodified and gzipped" << endl;
		return -1;
	}

	//filstreams
	ifstream Amplicons_in(argv[1]);
	ifstream R1_in(R1FASTQ, ios_base::in | ios_base::binary);
	ifstream R2_in(R2FASTQ, ios_base::in | ios_base::binary);
	ofstream SAM_out((string) argv[4] + ".sam");
	ofstream STATS_out((string) argv[4] + "_MappingStats.txt");

	//populate amplicon records
	if (GetAmplicons(Amplicons_in, AmpliconRecords, SamHeaders) == 1) {
		return -1;
	}

	try 
	{
		//prepare gzip decompression streams for file reading
		R1FilterStream.push(boost::iostreams::gzip_decompressor());
		R1FilterStream.push(R1_in);
		R2FilterStream.push(boost::iostreams::gzip_decompressor());
		R2FilterStream.push(R2_in);

		//parse FASTQs
		if (R1_in.is_open() && R2_in.is_open()) {
			while (R1FilterStream.good() && R2FilterStream.good()) {
				getline(R1FilterStream, Read1Line);
				getline(R2FilterStream, Read2Line);

				//skip empty lines
				if (Read1Line == "" || Read2Line == "") {
					continue;
				}

				LineNo++;

				if (LineNo == 1) {

					Header = Read1Line.substr(1, Read1Line.find_first_of(' ') - 1);
					Header1 = Read1Line;
					Header2 = Read2Line;

					TotalReads++;

					if (TotalReads < 15) { //check the first few reads

						//check read Headers are the same in both files
						if (Header != Read2Line.substr(1, Read2Line.find_first_of(' ') - 1)) {
							std::cerr << "ERROR: Read header " << Header << " does not match " << Read2Line.substr(0, Read2Line.find_first_of(' ')) << endl;
							std::cerr << "ERROR: FASTQ read headers are not synchronised" << endl;
							return -1;
						}

						//check read no is correct
						if (Read1Line.substr(Read1Line.find_first_of(' ') + 1, 1) != "1") {
							std::cerr << "ERROR: " << argv[2] << " contains R" << Read1Line.substr(Read1Line.find_first_of(' ') + 1, 1) << " reads" << endl;
							return -1;
						}
						if (Read2Line.substr(Read2Line.find_first_of(' ') + 1, 1) != "2") {
							std::cerr << "ERROR: " << argv[3] << " contains R" << Read2Line.substr(Read2Line.find_first_of(' ') + 1, 1) << " reads" << endl;
							return -1;
						}

						if (TotalReads == 1) {

							//check the Index number is the same across FASTQs and reads
							Index = Read1Line.substr(Read1Line.find_last_of(':') + 1, std::string::npos); //set Index no

							if (Index != Read2Line.substr(Read2Line.find_last_of(':') + 1, std::string::npos)) {
								std::cerr << "ERROR: Index in read Headers " << Read1Line << " and " << Read2Line << " do not match" << endl;
								return -1;
							}

							FlowCellID = GetFlowCellID(Header); //set flowcell ID

							//write SAM Headers to file
							if (SAM_out.is_open()) {

								if (SamHeaders.size() == 0) {
									std::cerr << "ERROR: No SAM Headers were provided in the reference file. You must apply these manually to pass Picard validation." << endl;
								} else {

									//print sam Headers
									for (n = 0; n < SamHeaders.size(); ++n) {
										SAM_out << SamHeaders[n] << "\012";
									}

								}

								SAM_out << "@RG\tID:" << argv[4] << '_' << FlowCellID << "\tSM:" << argv[4] << "\tPL:ILLUMINA\tLB:" << argv[4] << "\012";
								SAM_out << "@PG\tID:IndelAmpliconAligner\tPN:IndelAmpliconAligner\tCL:";
								STATS_out << "#ID:" << argv[4] << '_' << FlowCellID << "\n";
								STATS_out << "#CL:";

								//print command line arguments
								for (n = 0; n < argc; ++n) {
									if (n == 0) {
										SAM_out << argv[n];
										STATS_out << argv[n];
									} else {
										SAM_out << ' ' << argv[n];
										STATS_out << ' ' << argv[n];
									}
								}

								SAM_out << "\tVN:" << Version << "\012";
								STATS_out << "\n#PG:IndelAmpliconAligner v" << Version << "\n";
								SAM_out << "@CO\tReads were globally Aligned using amplicon specific reference sequences\012";

							} else {
								std::cerr << "ERROR: Could not write headers to SAM file. Check file is not in use." << endl;
								return -1;
							}


						} else if (Index != Read1Line.substr(Read1Line.find_last_of(':') + 1, std::string::npos)) {
							std::cerr << "ERROR: Read headers contain mixed indexes" << endl;
							return -1;
						} else if (Index != Read2Line.substr(Read2Line.find_last_of(':') + 1, std::string::npos)) {
							std::cerr << "ERROR: Read headers contain mixed indexes" << endl;
							return -1;
						}

					}

				} else if (LineNo == 2) {
					Seq1 = Read1Line;
					Seq2 = Read2Line;
				} else if (LineNo == 4) {
					Qual1 = Read1Line;
					Qual2 = Read2Line;

					LineNo = 0;

					//skip N masked reads
					if (isReadNMasked(Seq1) == true || isReadNMasked(Seq2) == true) {
						nMaskedReads++;
						continue;
					}

					//iterate over amplicons
					for (n = 0; n < AmpliconRecords.size(); ++n) {

						if (MatchPrimer(Seq1, AmpliconRecords[n].LeftPrimer) == 1) { //smith-waterman alignment to position 0
							if (MatchPrimer(Seq2, AmpliconRecords[n].RightPrimer) == 1) {
								//read matches to this amplicon

								PrimerMatchedReads++; //total number of ontarget reads

								//trim adapter
								RightPrimerClipper(Seq1, Qual1, AmpliconRecords[n].RightPrimer);
								RightPrimerClipper(Seq2, Qual2, AmpliconRecords[n].LeftPrimer);

								//reduce primer dimer; insert size less than minIsize ignored
								if (Seq1.length() >= AmpliconRecords[n].LeftPrimerLen + AmpliconRecords[n].RightPrimerLen + minIsize &&
									Seq2.length() >= AmpliconRecords[n].LeftPrimerLen + AmpliconRecords[n].RightPrimerLen + minIsize) {

									Stats[AmpliconRecords[n].ID].Usable++; //usable reads by amplicon
									TotalUsableReads++;

									//merge reads into 1 contig
									if (ReadMerger(Seq1, Qual1, Seq2, Qual2, MaxQScore, QScorePhredOffset, MergedRead) == 1) { //MergedRead contains seq and qual merged
										Stats[AmpliconRecords[n].ID].Merged++; //merged reads

										/*TODO: need to iterate over all possilble reference sequences including off-target here*/

										//convert read and reference sequence to + strand for Alignment
										if (AmpliconRecords[n].Strand == false) { //is+Strand

											//Ref and Query must be reverse ConvertDNAComplemented to Reflect + strand
											MergedRead.first = ReverseComplement(MergedRead.first);
											reverse(MergedRead.second.begin(), MergedRead.second.end());

											LeftPrimerLengthStrandConverted = AmpliconRecords[n].RightPrimerLen;
											RightPrimerLengthStrandConverted = AmpliconRecords[n].LeftPrimerLen;

										} else {
											LeftPrimerLengthStrandConverted = AmpliconRecords[n].LeftPrimerLen;
											RightPrimerLengthStrandConverted = AmpliconRecords[n].RightPrimerLen;
										}

										//load sequences into Alignment matrix
										Query = MergedRead.first;
										Ref = AmpliconRecords[n].RefSeq;

										seqan::assignSource(seqan::row(Align, 1), Query);
										seqan::assignSource(seqan::row(Align, 0), Ref);

										//global pairwise Alignment
										NWScore = seqan::globalAlignment(Align, seqan::Score<int, seqan::Simple>(1, -3, -1, -8)); //match mismatch gapextend gapopen

										if (NWScore < 0) {
											break; //poor Alignment; discard this read and proceed to next
										}

										//calculate cigar string and edit distance for SAM output 
										TRow &row1 = seqan::row(Align, 0);
										TRow &row2 = seqan::row(Align, 1);

										if (getCigarNM(row1, row2, LeftPrimerLengthStrandConverted, RightPrimerLengthStrandConverted, CigarNM, SingleBaseMisMatchFrequency) == 1) {
											break;
										} else if ((float)SingleBaseMisMatchFrequency / AmpliconRecords[n].RefSeq.length() > MaxSingleBaseMisMatch) { //too many single base mismatches (false alignment)
											break;
										}

										//write Alignment to SAM
										if (SAM_out.is_open()) {

											if (AmpliconRecords[n].Strand == true) { //is+Strand
												SAM_out << Header << "\t0\t";
											} else if (AmpliconRecords[n].Strand == false) {
												SAM_out << Header << "\t16\t"; //read was reverse ConvertDNAComplemented
											}

											//downscale mapping score in acceptable range
											if (NWScore > 60) {
												SAM_out << AmpliconRecords[n].Chrom << "\t" << AmpliconRecords[n].Pos << "\t" << 60 << "\t" << CigarNM.first << "\t*\t0\t0\t" << MergedRead.first << "\t" << MergedRead.second;
											} else {
												SAM_out << AmpliconRecords[n].Chrom << "\t" << AmpliconRecords[n].Pos << "\t" << NWScore << "\t" << CigarNM.first << "\t*\t0\t0\t" << MergedRead.first << "\t" << MergedRead.second;
											}

											//optional fields
											SAM_out << "\tRG:Z:" << argv[4] << '_' << FlowCellID;
											SAM_out << "\tNM:i:" << CigarNM.second; //edit distance- including every base of an indel
											SAM_out << "\tAS:i:" << NWScore; //true alignment score
											SAM_out << "\tCO:Z:" << AmpliconRecords[n].ID << "\012"; //amplicon name

											Stats[AmpliconRecords[n].ID].Mapped++; //mapped reads by amplicon
											TotalMappedReads++;

										} else {
											std::cerr << "ERROR: Could not ouput Alignments to SAM file. Check file is not in use." << endl;
											return -1;
										}

									} else {
										TotalNotMergedReads++;
									}

								}

							} //check if right primer matches R2

							break; //if R1 primer matches do not continue looking for matches

						} //check if left primer matches R1
					} //end iterating over primer pairs

				}

			}

		} else {
			std::cerr << "ERROR: Unable to open FASTQ file(s)" << endl;
			return -1;
		}

	} catch (boost::iostreams::gzip_error& e) {
		std::cerr << "ERROR: Problem with gzip extraction: " << e.what() << endl;
		return -1;
	} catch (exception& e) {
		std::cerr << "ERROR: " << e.what() << endl;
		return -1;
	}

	//mapping stats
	STATS_out << "#TotalReads:" << TotalReads << "\n";
	STATS_out << "#PrimerMatchedPairs:" << PrimerMatchedReads << "\n";
	STATS_out << "#UsablePairs:" << TotalUsableReads << "\n";
	STATS_out << "#UnmergedPairs:" << TotalNotMergedReads << ' ' << (float)TotalNotMergedReads / TotalUsableReads * 100 << "%\n";
	STATS_out << "#TotalAlignedPairs:" << TotalMappedReads << ' ' << (float)TotalMappedReads / TotalUsableReads * 100 << "%\n";
	
	STATS_out << "#Amplicon\tUsableReads\tMergedReads\tMappedReads\n";
	for (n = 0; n < AmpliconRecords.size(); ++n) {
		STATS_out << AmpliconRecords[n].ID << "\t" << Stats[AmpliconRecords[n].ID].Usable << "\t" << Stats[AmpliconRecords[n].ID].Merged << "\t" << Stats[AmpliconRecords[n].ID].Mapped << "\n";
	}

	Amplicons_in.close();
	R1_in.close();
	R2_in.close();
	SAM_out.close();
	STATS_out.close();

	return 0;
}