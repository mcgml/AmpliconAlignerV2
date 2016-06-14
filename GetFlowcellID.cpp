/*
* Filename : GetFlowcellID.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Extract flowcellID from read header
* Status: Release
*/

#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;

string GetFlowCellID(const string& header) {

	vector<string> SplitVec; // Search for tokens
	boost::split(SplitVec, header, boost::is_any_of(":"), boost::token_compress_on); //substrings in elements

	return SplitVec[2].substr(SplitVec[2].find_first_of('-') + 1);

}