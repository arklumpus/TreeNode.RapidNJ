#ifndef DATALOADER_POINTER_H
#define DATALOADER_POINTER_H

#include "stdinclude.h"
#include "dataloader.hpp"
#include "bitStringUtils.hpp"

class dataloaderPointer : public dataloader {

public:
	dataloaderPointer(InputType sequenceType, int inputSequenceCount, int inputSequenceLength, int* inputSequenceNamesLengths, char** inputSequenceNames, char** inputSequenceData) {
		sequences = NULL;
		bitStrings = NULL;
		gapFilters = NULL;
		sequenceLength = 0;
		sequenceCount = 0;
		sequenceNames = new vector<string>;

		type = sequenceType;
		sequenceLength = inputSequenceLength;

		if (fastdist) {
			bitStrings = new vector<unsigned int*>;

			if (type == DNA) {
				bitStringsCount = sequenceLength / 64 + 6;
				paddingLength = bitStringsCount * 64 - sequenceLength;
				gapFilters = new vector<unsigned int*>;
			}
			else {
				bitStringsCount = sequenceLength / 16 + 8;
				paddingLength = bitStringsCount * 16 - sequenceLength;
			}
		}
		else {
			sequences = new vector<char*>;
		}

		for (int i = 0; i < inputSequenceCount; i++)
		{
			int nameLength = inputSequenceNamesLengths[i];

			string name(inputSequenceNames[i], inputSequenceNamesLengths[i]);

			storeSequence(&name, inputSequenceData[i]);
		}
	}

	void storeSequence(string* name, char* characters)
	{
		if (fastdist) {
			unsigned int* bitString = (unsigned int*)_mm_malloc(bitStringsCount * 4 * sizeof(unsigned int), 16);
			if (type == DNA) {
				unsigned int* gapFilter;
				gapFilter = (unsigned int*)_mm_malloc(bitStringsCount * 4 * sizeof(unsigned int), 16);
				encodeDNASequence(bitString, gapFilter, characters);
				gapFilters->push_back(gapFilter);
			}
			else {
				encodeProteinSequence(bitString, characters);
			}
			bitStrings->push_back(bitString);
			sequenceNames->push_back(*name);
			sequenceCount++;
		}
		else
		{
			char* data = new char[sequenceLength];
			for (unsigned int i = 0; i < sequenceLength; i++) {
				data[i] = resolveChar(characters[i]);
			}
			sequences->push_back(data);
			sequenceNames->push_back(*name);
			sequenceCount++;
		}
	}

	void load(string filename) {

	}

	inline char resolveChar(char c) {
		//ignore positions with ambigious/unknown nucleotides
		if (type == DNA) {
			if (!(c == 'a' || c == 'A' || c == 'c' || c == 'C' || c == 'g' || c == 'G' || c == 't' || c == 'T' || c == 'u' || c == 'U')) {
				return '-';
			}
			else {
				return c;
			}
		}
		else {
			if (c == '-' || c == '.' || c == 'X' || c == 'x' || c == 'z' || c == 'Z' || c == 'b' || c == 'B' || c == 'J' || c == 'j' || c == '?') {
				return '-';
			}
			else {
				return c;
			}
		}
	}

	inline void encodeProteinSequence(unsigned int* bitString, char* data) {
		for (unsigned int i = 0; i < sequenceLength; i++) {
			int offset = i % 4;
			unsigned int bitStringIdx = i / 4;
			if (offset == 0) {
				bitString[bitStringIdx] = 0;
			}
			char c = resolveChar(data[i]);
			bitString[bitStringIdx] += (c << (offset * 8));
		}

		for (unsigned int i = sequenceLength; i < paddingLength + sequenceLength; i++) {
			int offset = i % 4;
			unsigned int bitStringIdx = i / 4;
			if (offset == 0) {
				bitString[bitStringIdx] = 0;
			}
			char c = '-';
			bitString[bitStringIdx] += (c << (offset * 8));
		}
	}

	inline void encodeDNASequence(unsigned int* bitString, unsigned int* gapFilter, char* data) {
		for (unsigned int i = 0; i < sequenceLength; i++) {
			int offset = i % (BLOCK_SIZE / 2);
			unsigned int bitStringIdx = i / (BLOCK_SIZE / 2);
			if (offset == 0) {
				bitString[bitStringIdx] = 0;
				gapFilter[bitStringIdx] = 0;
			}
			char c = data[i];
			switch (c) {
				//Nothing is done for gaps and ambigious nucleotides
			case 'A':
			case 'a':
				//A has value 00 so nothing is added to the bitStrings
				gapFilter[bitStringIdx] += (Gbin << (offset * 2));
				break;
			case 'C':
			case 'c':
				bitString[bitStringIdx] += (Cbin << (offset * 2));
				gapFilter[bitStringIdx] += (Gbin << (offset * 2));
				break;
			case 'G':
			case 'g':
				bitString[bitStringIdx] += (Gbin << (offset * 2));
				gapFilter[bitStringIdx] += (Gbin << (offset * 2));
				break;
			case 'T':
			case 't':
				bitString[bitStringIdx] += (Tbin << (offset * 2));
				gapFilter[bitStringIdx] += (Gbin << (offset * 2));
				break;
			default:
				break;
			}
		}

		for (unsigned int i = sequenceLength; i < sequenceLength + paddingLength; i++) {
			int offset = i % (BLOCK_SIZE / 2);
			unsigned int bitStringIdx = i / (BLOCK_SIZE / 2);
			if (offset == 0) {
				bitString[bitStringIdx] = 0;
				gapFilter[bitStringIdx] = 0;
			}
		}
	}

	unsigned int** getBitStrings() {
		return &(*bitStrings)[0];
	}

	unsigned int** getGapFilters() {
		return &(*gapFilters)[0];
	}

	unsigned int getSequenceCount() {
		return sequenceCount;
	}

	unsigned int getSequenceLength() {
		return sequenceLength;
	}

	unsigned int getBitStringsCount() {
		return bitStringsCount;
	}

	vector<string>* getSequenceNames() {
		return sequenceNames;
	}

	vector<char*>* getSequences() {
		return sequences;
	}

	void setSequences(vector<char*>* val) {
		sequences = val;
	}

	~dataloaderPointer() {
		if (sequenceNames != NULL) {
			delete sequenceNames;
		}

		if (sequences != NULL) {
			for (unsigned int i = 0; i < sequences->size(); i++) {
				delete[] sequences->at(i);
			}
			delete sequences;
		}
		if (bitStrings != NULL) {
			for (unsigned int i = 0; i < bitStrings->size(); i++) {
				_mm_free(bitStrings->at(i));
			}
			delete bitStrings;
		}
		if (gapFilters != NULL) {
			for (unsigned int i = 0; i < gapFilters->size(); i++) {
				_mm_free(gapFilters->at(i));
			}
			delete gapFilters;
		}
	}

private:
	vector<string>* sequenceNames;
};

#endif