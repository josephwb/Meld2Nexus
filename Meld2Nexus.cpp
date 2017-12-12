/*

To compile, type the following in a unix prompt:

make

*/

// TODO: add fasta option ***

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

#include "General.h"

void usage();
void processCommandLineArguments (int argc, char *argv[], vector <string> & listFileNames, string & outName);
void getNumTaxaChar (string & seqFileName, int & numTaxa, int & numChar, bool & interleavedData);
vector < vector <string> > collectTaxaAlignment (string & seqFileName, int const& numTaxa,
	int const& numChar, bool const& interleavedData);
void addFile(string & fileName, vector <int> & geneLengths, vector <string> & geneNames,
	vector < vector <string> > & taxaAlignment);
string emptySequence (int const& length);
int sumLengths (vector <int> const& geneLengths);
void printNewNexus (vector <int> const& geneLengths, vector <string> const& geneNames,
	vector < vector <string> > const& taxaAlignment, string & outName);

bool DEBUG = false;

int main (int argc, char *argv[]) {
	vector <int> geneLengths;
	vector <string> geneNames;
	vector <string> taxonNames;
	vector <string> fileNames;
	vector < vector <string> > taxaAlignment;
	string outName = "Merged.NEX";
	
	processCommandLineArguments(argc, argv, fileNames, outName);
	
	for (int i = 0; i < (int)fileNames.size(); i++) {
		addFile(fileNames[i], geneLengths, geneNames, taxaAlignment);
	}
	
// write it already!
	printNewNexus (geneLengths, geneNames, taxaAlignment, outName);
	
	cout << endl << "Successfully created merged Nexus file '" << outName 
		<< "' from " << fileNames.size() << " alignments." << endl;
	
	cout << endl << "Fin." << endl;
	return 0;
}

void usage() {
	cout
	<< "Usage:" << endl
	<< endl
	<< "Call as either:" << endl
	<< endl
	<< "   ./Meld2Nexus -c file1.NEX file2.NEX file3.NEX ... [-o outFileName]" << endl
	<< endl
	<< "or" << endl
	<< endl
	<< "   ./Meld2Nexus -f fileContainingNexusFileNames [-o outFileName]" << endl
	<< endl;
}

void processCommandLineArguments (int argc, char *argv[], vector <string> & listFileNames, string & outName)
{
	if (argc == 1) {
		cout << "No arguments given." << endl << endl;
		usage();
		exit(0);
	} else {
		for (int i = 1; i < argc; i++) {
			string temp = argv[i];
			
			if (temp == "-h" || temp == "-help") {
				cout << endl
				<< "Program description: merge Nexus alignment files, preserving partition information." << endl
				<< endl
				<< "To compile, type the following in a unix prompt:" << endl << endl
				<< "   make" << endl << endl;
				usage();
				exit(0);
			} else if (temp == "-f") { // read from file
				i++;
				listFileNames = readFileList(argv[i]);
				continue;
			} else if (temp == "-c") { // read from commandline
				i++;
				while (i < argc) {
					char c = argv[i][0];
					if (c == '-') {
						i--;
						break;
					}
					checkValidInputFile(argv[i]);
					listFileNames.push_back(argv[i]);
					i++;
				}
			} else if (temp == "-o") {
				i++;
				outName = argv[i];
			} else {
				cout
				<< "*** Unknown command-line argument '" << argv[i] << "' encountered. ***" << endl << endl;
				usage();
				cout << endl << "Exiting because of improper command-line argument '" << argv[i] << "' encountered." << endl << endl;
				exit(0);
			}
		}
	}
	checkValidOutputFile(outName, 0);
}

void getNumTaxaChar (string & seqFileName, int & numTaxa, int & numChar, bool & interleavedData) {
	ifstream inputUserFile;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	bool numTaxaEncountered = false;
	bool numCharEncountered = false;
	bool equalSignEncountered = false;
	bool semicolonEncountered = false;
	bool matrixEncountered = false; // know to stop looking
	
	inputUserFile.open(seqFileName.c_str());
	string line;
	
// Looking for pattern like 'dimensions ntax=53 nchar=16620;'
// 	- can be in either order, but must be stated on same line (for now)
// 	- no spaces allowed next to equal sign (for now)

// Looking for pattern like 'Format datatype=dna [gap=-] [missing=?] {[interleave=yes] or [interleave]};'
// 	- no spaces allowed next to equal sign (for now)

	while (getline(inputUserFile,line) && !matrixEncountered) {
		int stringPosition = 0;
		commentLine = checkCommentLineNexus(line);
		whiteSpaceOnly = checkWhiteSpaceOnly(line);
		if (line.empty() || commentLine || whiteSpaceOnly) {
			continue;
		} else {
			if (checkStringValue(line, "matrix", stringPosition)) {	// Done - won't find the information further down; really only used for 'interleave'
				matrixEncountered = true;
				continue;
			} else if (checkStringValue(line, "dimensions", stringPosition)) {
				stringPosition = 0;
				while (!numTaxaEncountered || !numCharEncountered) {
					stringPosition++;
					string tempString = removeStringSuffix(parseString(line,stringPosition), '=', equalSignEncountered);
					if (checkStringValue(tempString, "ntax", 0)) {
						tempString = removeStringPrefix(parseString(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', numTaxaEncountered);
						numTaxa = convertStringtoInt(tempString);
						if (DEBUG) {cout << "NTax = " << numTaxa << endl;}
						numTaxaEncountered = true;
					}
					if (checkStringValue(tempString, "nchar", 0)) {
						tempString = removeStringPrefix(parseString(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', numCharEncountered);
						numChar = convertStringtoInt(tempString);
						if (DEBUG) {cout << "NChar = " << numChar << endl;}
						numCharEncountered = true;
					}
				}
			} else if (checkStringValue(line, "format", stringPosition)) {
				stringPosition = 0;
				while (!semicolonEncountered) {
					stringPosition++;
					
					string tempString = removeStringSuffix(parseString(line,stringPosition), '=', equalSignEncountered); // where '=' is used i.e. 'interleave=yes;'
					if (checkStringValue(tempString, "interleave", 0)) {
						tempString = removeStringPrefix(parseString(line,stringPosition), '=');
						tempString = removeStringSuffix(tempString, ';', semicolonEncountered);
						if (checkStringValue(tempString, "yes", 0)) {
							interleavedData = true;
							semicolonEncountered = true;
							cout << "Data are in interleaved format." << endl;
							continue;
						} else if (checkStringValue(tempString, "no", 0)) {
							interleavedData = false;
							semicolonEncountered = true;
							cout << "Data are not in interleaved format." << endl;
							continue;
						}
					}
					tempString = removeStringSuffix(parseString(line,stringPosition), ';', semicolonEncountered); // where '=' is NOT used i.e. 'interleave;'
					if (checkStringValue(tempString, "interleave", 0)) {										  // OR a space follows interleave declaration i.e. 'interleave [=];'
						interleavedData = true;
						semicolonEncountered = true;
						cout << "Data are in interleaved format." << endl;
						continue;
					}
				}
			}
		}
	}
	inputUserFile.close();
}

vector < vector <string> > collectTaxaAlignment (string & seqFileName, int const& numTaxa,
	int const& numChar, bool const& interleavedData)
{
// PLEASE NOTE: search strategy below uses very strict format assumptions - that which is exported by PAUP*
// 	- do not be surprised if this fucks up - it is probably a simple rearrangement of terms

	vector < vector <string> > taxaAlignment;
	ifstream inputAlignment;
	vector <string> tempStringVector;
	bool commentLine = false;
	bool whiteSpaceOnly = false;
	bool matrixEncountered = false;
	bool allCharacterRead = false;
	int numCharRead = 0;
	
	inputAlignment.open(seqFileName.c_str());
	string line;
	
	if (!interleavedData) {
		while (!matrixEncountered) { // Ignore lines until 'matrix' is encountered
			getline(inputAlignment,line);
			commentLine = checkCommentLineNexus(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly) {
				continue;
			} else if (checkStringValue(line, "matrix", 0)) {
				matrixEncountered = true;
			} else {
				continue;
			}
		}
		numCharRead = 0;
// Read in every non-empty (or non-whitespace), non-commented-out line
		for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++) {
			getline(inputAlignment,line);
			commentLine = checkCommentLineNexus(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly) {
				taxonIter--;
				continue;
			} else {
// First string is taxon name, second is sequence
				if (DEBUG) {cout << "Reading in taxon '" << parseString(line, 0) << "'..." << endl;}
				tempStringVector.push_back(parseString(line, 0));
				tempStringVector.push_back(parseString(line, 1));
				
				taxaAlignment.push_back(tempStringVector);
				
				tempStringVector.clear();
				
// Count sites encountered - for error-checking; only checking first sequence (for now)
				if (taxonIter == 0) {
					int charCounter = 0;
					for (string::size_type iterCharacters = 0; iterCharacters < parseString(line, 1).size(); iterCharacters++) {
						charCounter++;
					}
					numCharRead += charCounter;
				}
			}
		}
		if (numCharRead == numChar) {
			allCharacterRead = true;
			if (DEBUG) {cout << "numCharRead (" << numCharRead << ") == numChar (" << numChar << ") declared in Nexus file. Woo-hoo!" << endl;}
		}
	} else if (interleavedData) {
// Ignore lines until 'matrix' is encountered
		while (!matrixEncountered) {
			getline(inputAlignment,line);
			commentLine = checkCommentLineNexus(line);
			whiteSpaceOnly = checkWhiteSpaceOnly(line);
			
			if (line.empty() || commentLine || whiteSpaceOnly) {
				continue;
			} else if (checkStringValue(line, "matrix", 0)) {
				matrixEncountered = true;
			} else {
				continue;
			}
		}
		bool firstPass = true;
		numCharRead = 0;
		while (!allCharacterRead) {
// Read in every non-empty (or non-whitespace), non-commented-out line
			for (int taxonIter = 0; taxonIter < numTaxa; taxonIter++) {
				getline(inputAlignment,line);
				commentLine = checkCommentLineNexus(line);
				whiteSpaceOnly = checkWhiteSpaceOnly(line);
				
				if (line.empty() || commentLine || whiteSpaceOnly) {
					taxonIter--;
					continue;
				} else {
// First string is taxon name, second is sequence
					string taxonName = (parseString(line, 0));
					string taxonSequence = (parseString(line, 1));
					
					if (firstPass) {
						if (DEBUG) {cout << "Reading in taxon '" << parseString(line, 0) << "'..." << endl;}
						tempStringVector.push_back(taxonName);		// Taxon name
						tempStringVector.push_back(taxonSequence);	// sequence
						taxaAlignment.push_back(tempStringVector);
					}
					if (!firstPass) {
						taxaAlignment[taxonIter][1] += taxonSequence;
					}
					
// Count sites encountered - for error-checking; only checking first sequence (for now)
					if (taxonIter == 0) {
						int charCounter = 0;
						for (string::size_type iterCharacters = 0; iterCharacters < taxonSequence.size(); iterCharacters++) {
							charCounter++;
						}
						numCharRead += charCounter;
					}
					tempStringVector.clear();
				}
			}
			firstPass = false;
			if (numCharRead == numChar) {
				allCharacterRead = true;
				if (DEBUG) {cout << "numCharRead (" << numCharRead << ") == numChar (" << numChar << ") declared in Nexus file. Woo-hoo!" << endl;}
			}
		}
	}
	return taxaAlignment;
}

void addFile (string & fileName, vector <int> & geneLengths, vector <string> & geneNames, vector < vector <string> > & taxaAlignment) {
	vector < vector <string> > curAlignment;
	vector <string> tempStringVector;
	string missingData;
	
	int sumLength = sumLengths(geneLengths);
	int sumTaxa = 0;
	if (sumLength != 0) {
		sumTaxa = taxaAlignment.size();
	}
	
// for current file
	int numTaxa = 0;
	int numChar = 0;
	bool interleavedData = false;
	
	getNumTaxaChar(fileName, numTaxa, numChar, interleavedData);
	curAlignment = collectTaxaAlignment(fileName, numTaxa, numChar, interleavedData);
	
	geneLengths.push_back(numChar);
	geneNames.push_back(getRootName(fileName));
	
	if (taxaAlignment.empty()) { // first file
		taxaAlignment = curAlignment;
	} else {
		bool matchTaxon = false;
		missingData = emptySequence(sumLength);
		for (int i = 0; i < numTaxa; i++) { // loop over taxa in current file
			matchTaxon = false;
			string taxon = curAlignment[i][0];
			string seq = curAlignment[i][1];
			for (int j = 0; j < sumTaxa; j++) { // loop over taxa in master
				if (taxon == taxaAlignment[j][0]) {
					matchTaxon = true;
					taxaAlignment[j][1] += seq;
					break;
				}
			}
			if (!matchTaxon) { // new taxon
				if (DEBUG) {cout << "New taxon '" << taxon << "' encountered" << endl;}
				tempStringVector.push_back(taxon); // name
				seq = missingData + seq;
				tempStringVector.push_back(seq);
				taxaAlignment.push_back(tempStringVector);
			}
			tempStringVector.clear();
		}
// fill in for taxa in master file that were not present in current file. need only consider original sumTaxa.
// not most efficient way
		sumLength += numChar;
		missingData = emptySequence(numChar);
		for (int i = 0; i < sumTaxa; i++) {
			if ((int)taxaAlignment[i][1].size() != sumLength) {
				taxaAlignment[i][1] += missingData;
			}
		}
	}
}

void printNewNexus (vector <int> const& geneLengths, vector <string> const& geneNames,
	vector < vector <string> > const& taxaAlignment, string & outName)
{
	ofstream NEXUS;
	
	int numTaxa = taxaAlignment.size();
	int start = 1;
	int stop = 0;
	
	vector <string> tempVect;
	for (int i = 0; i < numTaxa; i++) {
		tempVect.push_back(taxaAlignment[i][0]);
	}
	
	cout << endl << "Combined data involves " << numTaxa << " unique taxon names." << endl;
	string longest = getLongestName(tempVect);
	
	//cout << "Longest name = " << longest << endl;
	
	NEXUS.open(outName.c_str());
	
	NEXUS
	<< "#NEXUS" << endl
	<< endl
	<< "begin data;" << endl
	<< "	dimensions ntax=" << numTaxa << " nchar=" << taxaAlignment[0][1].size() << ";" << endl
	<< "	format datatype=rna missing=? gap=-;" << endl
	<< "matrix" << endl;
	
	tempVect.clear();
	for (int i = 0; i < numTaxa; i++) {
		string formatted = addFormattingSpaces(longest, taxaAlignment[i][0]);
		formatted = formatted + "	" + taxaAlignment[i][1];
		tempVect.push_back(formatted);
	}
	
	// sort
	sort(tempVect.begin(), tempVect.end());
	
	// print alignment
	for (int i = 0; i < numTaxa; i++) {
		NEXUS << tempVect[i] << endl;
	}
	NEXUS << ";" << endl << endl
	<< "End;" << endl << endl;
	
	// print out CHARSET information;
	NEXUS << "Begin paup;" << endl << endl;
	
	for (int i = 0; i < (int)geneNames.size(); i++) {
		stop = start + geneLengths[i] - 1;
		NEXUS << "	CHARSET " << geneNames[i] << " = " << start << "-" << stop
			<< "; [" << geneLengths[i] << " characters]" << endl;
		start = stop + 1;
	}
	
	NEXUS << endl << "End;" << endl;
	
	NEXUS.close();
}

int sumLengths (vector <int> const& geneLengths) {
	int sum = 0;
	if (!geneLengths.empty()) {
		for (int i = 0; i < (int)geneLengths.size(); i++) {
			sum += geneLengths[i];
		}
	}
	return sum;
}

string emptySequence (int const& length) { // create a string of ?s
	string res = "-";
	for (int i = 1; i < length; i++) {
		res += "-";
	}
	return res;
}
