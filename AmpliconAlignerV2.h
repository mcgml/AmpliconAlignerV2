/*
* Filename : AmpliconAlignerV2.h
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Globally aligns paired-end reads to amplicon sequences and outputs in SAM format.
* Status: Release
*/

#include <string>
#include <seqan/Align.h>

using namespace std;

typedef struct {
	string ID;
	string Chrom;
	string RefSeq;
	unsigned Pos;
	string LeftPrimer;
	string RightPrimer;
	unsigned LeftPrimerLen;
	unsigned RightPrimerLen;
	bool Strand; //is+Strand
} AmpliconRecord;

typedef struct {
	unsigned Usable;
	unsigned Merged;
	unsigned Mapped;
} Stat;

typedef seqan::String<char> TSequence;                 // sequence type
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;      // align type
typedef seqan::Row<TAlign>::Type TRow;
typedef seqan::Iterator<TRow>::Type TRowIterator;

string GetFlowCellID(const string& header);
bool MatchPrimer(const string& Seq, const string& Primer);
string GetFlowCellID(const string& header);
void RightPrimerClipper(string& Seq, string& Qual, string Primer);
string ReverseComplement(const string& DNA);
bool ReadMerger(const string& SeqR1, const string& QualR1, string SeqR2, string QualR2,
	const unsigned MaxQScore, const unsigned QScorePhredOffset, pair<string, string>& MergedRead);
bool GetAmplicons(ifstream& Amplicons_in, vector<AmpliconRecord>& AmpliconRecords, vector<string>& SamHeaders);
bool isStringDNA(const string& str);
bool getCigarNM(TRow& row1, TRow& row2, const unsigned LeftPrimerLengthStrandConverted, const unsigned RightPrimerLengthStrandConverted, pair<string, unsigned>& CigarNM, unsigned& SingleBaseMisMatchFrequency);
bool isReadNMasked(const string& read);
