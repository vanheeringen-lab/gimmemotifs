
// Copyright 2009, 2010, 2011, 2012 Christopher Benner <cbenner@gmail.com>
// 
// This file is part of HOMER
//
// HOMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HOMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.


#include "Hashtable.h"

LinkedListItem::LinkedListItem() {
	next = NULL;
	obj = NULL;
}
LinkedList::LinkedList() {
	firstItem = NULL;
	lastItem = NULL;
	total = 0;
}
LinkedList::~LinkedList() {
	while (firstItem != NULL) {
		LinkedListItem* next = firstItem->next;
		delete firstItem;
		firstItem = next;
	}
	lastItem = NULL;
	total = 0;
}
void LinkedList::add(void* obj) {
	LinkedListItem* item = new LinkedListItem();
	item->obj = obj;
	if (lastItem == NULL) {
		firstItem = item;
	} else {
		lastItem->next = item;
	}
	lastItem = item;
	total++;
}
void* LinkedList::get(int index) {
	LinkedListItem* item = firstItem;
	for (int i=1;i<index;i++) {
		if (item == NULL) return NULL;
		item = item->next;
	}
	if (item == NULL) return NULL;
	return item->obj;
}
void* LinkedList::remove(int index) {
	if (firstItem == NULL) return NULL;

	LinkedListItem* item = firstItem;
	LinkedListItem* last = NULL;
	for (int i=1;i<index;i++) {
		if (item == NULL) return NULL;
		last = item;
		item = item->next;
	}
	if (item == NULL) return NULL;
	if (last == NULL) {//first item;
		firstItem = item->next;
	}
	if (item->next == NULL) {
		lastItem = NULL;
		if (last != NULL) {
			lastItem = last;
			last->next = NULL;
		}
	} else {
		if (last != NULL) {
			last->next = item->next;
		}
	}
	void *obj = item->obj;
	delete item;
	total--;
	return obj;
}
void** LinkedList::toArray(int &numberOfItems) {
	if (total == 0) return NULL;
	void** array = new void*[total];
	for (int i=0;i<total;i++) array[i] = NULL;
	numberOfItems = 0;
	LinkedListItem* item = firstItem;
	while (item != NULL) {
		array[numberOfItems++] = item->obj;
		item = item->next;
	}
	return array;
}


Hashlist::Hashlist(void* obj, char* str) {
	elements = new void*[1];
	hashstrs = new char*[1];
	lowMemFlag = 0;
	elements[0] = obj;
	hashstrs[0] = new char[strlen(str)+1];
	strcpy(hashstrs[0], str);
	length = 1;
}
Hashlist::Hashlist() {
}
Intlist::Intlist() {
}
LongIntlist::LongIntlist() {
}
Floatlist::Floatlist() {
}
Doublelist::Doublelist() {
}
Intlist::Intlist(int I, char* str) {
	elements = new int[1];
	hashstrs = new char*[1];
	elements[0] = I;
	hashstrs[0] = new char[strlen(str)+1];
	strcpy(hashstrs[0], str);
	length = 1;
}
LongIntlist::LongIntlist(long long int I, char* str) {
	elements = new long long int[1];
	hashstrs = new char*[1];
	elements[0] = I;
	hashstrs[0] = new char[strlen(str)+1];
	strcpy(hashstrs[0], str);
	length = 1;
}
Floatlist::Floatlist(float I, char* str) {
	elements = new float[1];
	hashstrs = new char*[1];
	elements[0] = I;
	hashstrs[0] = new char[strlen(str)+1];
	strcpy(hashstrs[0], str);
	length = 1;
}
Doublelist::Doublelist(double I, char* str) {
	elements = new double[1];
	hashstrs = new char*[1];
	elements[0] = I;
	hashstrs[0] = new char[strlen(str)+1];
	strcpy(hashstrs[0], str);
	length = 1;
}


Hashlist::~Hashlist() {
	if (hashstrs != NULL) {
		for (int i=0;i<length;i++) {
			if (hashstrs[i] != NULL) {
				if (lowMemFlag==0) {
					delete[] (hashstrs[i]);
					hashstrs[i] = NULL;
				} else {
					hashstrs[i] = NULL;
				}
			}
		}
		delete[] hashstrs;
		hashstrs = NULL;
	}
	if (elements != NULL) {
		delete[] elements;
		elements = NULL;
	}
}
Intlist::~Intlist() {
	if (hashstrs != NULL) {
		for (int i=0;i<length;i++) {
			if (hashstrs[i] != NULL) {
				delete[] (hashstrs[i]);
				hashstrs[i] = NULL;
			}
		}
		delete[] hashstrs;
		hashstrs= NULL;
	}
	if (elements != NULL) {
		delete[] elements;
		elements = NULL;
	}
}
LongIntlist::~LongIntlist() {
	if (hashstrs != NULL) {
		for (int i=0;i<length;i++) {
			if (hashstrs[i] != NULL) {
				delete[] (hashstrs[i]);
				hashstrs[i] = NULL;
			}
		}
		delete[] hashstrs;
		hashstrs= NULL;
	}
	if (elements != NULL) {
		delete[] elements;
		elements = NULL;
	}
}
Floatlist::~Floatlist() {
	if (hashstrs != NULL) {
		for (int i=0;i<length;i++) {
			if (hashstrs[i] != NULL)
				delete[] (hashstrs[i]);
				hashstrs[i] = NULL;
		}
		delete[] hashstrs;
		hashstrs = NULL;
	}
	if (elements != NULL) {
		delete[] elements;
		elements = NULL;
	}
}
Doublelist::~Doublelist() {
	if (hashstrs != NULL) {
		for (int i=0;i<length;i++) {
			if (hashstrs[i] != NULL)
				delete[] (hashstrs[i]);
				hashstrs[i] = NULL;
		}
		delete[] hashstrs;
		hashstrs = NULL;
	}
	if (elements != NULL) {
		delete[] elements;
		elements = NULL;
	}
}


void* Hashlist::search(char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i]) == 0) {
			return elements[i];
		}
	}
	return NULL;
}
int Intlist::search(char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i]) == 0) {
			return elements[i];
		}
	}
	return EMPTY_INT;
}
long long int LongIntlist::search(char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i]) == 0) {
			return elements[i];
		}
	}
	return EMPTY_INT;
}
float Floatlist::search(char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i]) == 0) {
			return elements[i];
		}
	}
	return EMPTY_FLOAT;
}
double Doublelist::search(char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i]) == 0) {
			return elements[i];
		}
	}
	return EMPTY_DOUBLE;
}




int Hashlist::insert(void* obj, char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			elements[i] = obj;
			return 0;
		}
	}
	void** newelements = new void*[length+1];
	char** newnames = new char*[length+1];
	for (int i=0;i<length;i++) {
		newelements[i] = elements[i];
		newnames[i] = hashstrs[i];
	}
	delete []elements;
	delete []hashstrs;
	elements = newelements;
	hashstrs = newnames;
	elements[length] = obj;
	hashstrs[length] = new char[strlen(str)+1];
	strcpy(hashstrs[length], str);
	length++;
	return 1;
}
int Intlist::insert(int I, char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			elements[i] = I;
			return 0;
		}
	}
	int* newelements = new int[length+1];
	char** newnames = new char*[length+1];
	for (int i=0;i<length;i++) {
		newelements[i] = elements[i];
		newnames[i] = hashstrs[i];
	}
	delete []elements;
	delete []hashstrs;
	elements = newelements;
	hashstrs = newnames;
	elements[length] = I;
	hashstrs[length] = new char[strlen(str)+1];
	strcpy(hashstrs[length], str);
	length++;
	return 1;
}
int LongIntlist::insert(long long int I, char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			elements[i] = I;
			return 0;
		}
	}
	long long int* newelements = new long long int[length+1];
	char** newnames = new char*[length+1];
	for (int i=0;i<length;i++) {
		newelements[i] = elements[i];
		newnames[i] = hashstrs[i];
	}
	delete []elements;
	delete []hashstrs;
	elements = newelements;
	hashstrs = newnames;
	elements[length] = I;
	hashstrs[length] = new char[strlen(str)+1];
	strcpy(hashstrs[length], str);
	length++;
	return 1;
}
int Floatlist::insert(float I, char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			elements[i] = I;
			return 0;
		}
	}
	float* newelements = new float[length+1];
	char** newnames = new char*[length+1];
	for (int i=0;i<length;i++) {
		newelements[i] = elements[i];
		newnames[i] = hashstrs[i];
	}
	delete []elements;
	delete []hashstrs;
	elements = newelements;
	hashstrs = newnames;
	elements[length] = I;
	hashstrs[length] = new char[strlen(str)+1];
	strcpy(hashstrs[length], str);
	length++;
	return 1;
}
int Doublelist::insert(double I, char* str) {
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			elements[i] = I;
			return 0;
		}
	}
	double* newelements = new double[length+1];
	char** newnames = new char*[length+1];
	for (int i=0;i<length;i++) {
		newelements[i] = elements[i];
		newnames[i] = hashstrs[i];
	}
	delete []elements;
	delete []hashstrs;
	elements = newelements;
	hashstrs = newnames;
	elements[length] = I;
	hashstrs[length] = new char[strlen(str)+1];
	strcpy(hashstrs[length], str);
	length++;
	return 1;
}

void* Hashlist::remove(char* str) {
	int found = 0;
	void* returnObj=NULL;
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			found = 1;
			returnObj = elements[i];
			if (lowMemFlag==0) {
				delete[] (hashstrs[i]);
			}
		}
		if (found) {
			if (i==length-1) {
				elements[i]=NULL;
				hashstrs[i]=NULL;
			} else {
				elements[i] = elements[i+1];
				hashstrs[i] = hashstrs[i+1];
			}
		}
	}
	if (found) length--;
	return returnObj;
}
	
int Intlist::remove(char* str) {
	int found = 0;
	int returnObj=EMPTY_INT;
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			found = 1;
			returnObj = elements[i];
			delete[] (hashstrs[i]);
		}
		if (found) {
			if (i==length-1) {
				elements[i]=EMPTY_INT;
				hashstrs[i]=NULL;
			} else {
				elements[i] = elements[i+1];
				hashstrs[i] = hashstrs[i+1];
			}
		}
	}
	if (found) length--;
	return returnObj;
}
int LongIntlist::remove(char* str) {
	int found = 0;
	long long int returnObj=EMPTY_INT;
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			found = 1;
			returnObj = elements[i];
			delete[] (hashstrs[i]);
		}
		if (found) {
			if (i==length-1) {
				elements[i]=EMPTY_INT;
				hashstrs[i]=NULL;
			} else {
				elements[i] = elements[i+1];
				hashstrs[i] = hashstrs[i+1];
			}
		}
	}
	if (found) length--;
	return returnObj;
}
	
float Floatlist::remove(char* str) {
	int found = 0;
	float returnObj=EMPTY_FLOAT;
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			found = 1;
			returnObj = elements[i];
			delete[] (hashstrs[i]);
		}
		if (found) {
			if (i==length-1) {
				elements[i]=EMPTY_FLOAT;
				hashstrs[i]=NULL;
			} else {
				elements[i] = elements[i+1];
				hashstrs[i] = hashstrs[i+1];
			}
		}
	}
	if (found) length--;
	return returnObj;
}
double Doublelist::remove(char* str) {
	int found = 0;
	double returnObj=EMPTY_DOUBLE;
	for (int i=0;i<length;i++) {
		if (strcmp(str,hashstrs[i])==0) {
			found = 1;
			returnObj = elements[i];
			delete[] (hashstrs[i]);
		}
		if (found) {
			if (i==length-1) {
				elements[i]=EMPTY_DOUBLE;
				hashstrs[i]=NULL;
			} else {
				elements[i] = elements[i+1];
				hashstrs[i] = hashstrs[i+1];
			}
		}
	}
	if (found) length--;
	return returnObj;
}



Hashtable::Hashtable() {
	maxsize = DEFAULT_SIZE;
	init();
}
Inttable::Inttable() {
	maxsize = DEFAULT_SIZE;
	init();
}
LongInttable::LongInttable() {
	maxsize = DEFAULT_SIZE;
	init();
}
Floattable::Floattable() {
	maxsize = DEFAULT_SIZE;
	init();
}
Doubletable::Doubletable() {
	maxsize = DEFAULT_SIZE;
	init();
}
Hashtable::Hashtable(int size) {
	maxsize = size;
	init();
}
Inttable::Inttable(int size) {
	maxsize = size;
	init();
}
LongInttable::LongInttable(int size) {
	maxsize = size;
	init();
}
Floattable::Floattable(int size) {
	maxsize = size;
	init();
}
Doubletable::Doubletable(int size) {
	maxsize = size;
	init();
}
void Hashtable::init() {
	hash = NULL;
	lowMemFlag = 0;
	hash = new Hashlist*[maxsize];
	for (int i=0;i<maxsize;i++) {
		hash[i] = NULL;
	}
	total = 0;
}
void Inttable::init() {
	hash = NULL;
	hash = new Intlist*[maxsize];
	for (int i=0;i<maxsize;i++) {
		hash[i] = NULL;
	}
	total = 0;
}
void LongInttable::init() {
	hash = NULL;
	hash = new LongIntlist*[maxsize];
	for (int i=0;i<maxsize;i++) {
		hash[i] = NULL;
	}
	total = 0;
}
void Floattable::init() {
	hash = NULL;
	hash = new Floatlist*[maxsize];
	for (int i=0;i<maxsize;i++) {
		hash[i] = NULL;
	}
	total = 0;
}
void Doubletable::init() {
	hash = NULL;
	hash = new Doublelist*[maxsize];
	for (int i=0;i<maxsize;i++) {
		hash[i] = NULL;
	}
	total = 0;
}

Hashtable::~Hashtable() {
	if (hash != NULL) {
		for (int i=0;i<maxsize;i++) {
			if (hash[i] != NULL) {
				delete hash[i];
				hash[i] = NULL;
			}
		}
		delete []hash;
		hash = NULL;
	}
}
Inttable::~Inttable() {
	if (hash != NULL) {
		for (int i=0;i<maxsize;i++) {
			if (hash[i] != NULL) {
				delete hash[i];
				hash[i] = NULL;
			}
		}
		delete []hash;
		hash = NULL;
	}
}
LongInttable::~LongInttable() {
	if (hash != NULL) {
		for (int i=0;i<maxsize;i++) {
			if (hash[i] != NULL) {
				delete hash[i];
				hash[i] = NULL;
			}
		}
		delete []hash;
		hash = NULL;
	}
}
Floattable::~Floattable() {
	if (hash != NULL) {
		for (int i=0;i<maxsize;i++) {
			if (hash[i] != NULL) {
				delete hash[i];
				hash[i] = NULL;
			}
		}
		delete []hash;
		hash = NULL;
	}
}
Doubletable::~Doubletable() {
	if (hash != NULL) {
		for (int i=0;i<maxsize;i++) {
			if (hash[i] != NULL) {
				delete hash[i];
				hash[i] = NULL;
			}
		}
		delete []hash;
		hash = NULL;
	}
}
void Hashtable::setLowMem(int a) {
	lowMemFlag = a;
	for (int i=0;i<maxsize;i++) {
		if (hash[i] != NULL) {
			hash[i]->lowMemFlag = a;
		}
	}
}
void* Hashtable::search(char* hashstr) {
	if (hashstr == NULL) return NULL;
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return NULL;
	} else {
		return hash[bin]->search(hashstr);
	}
}
int Inttable::search(char* hashstr) {
	if (hashstr == NULL) return EMPTY_INT;
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_INT;
	} else {
		return hash[bin]->search(hashstr);
	}
}
long long int LongInttable::search(char* hashstr) {
	if (hashstr == NULL) return EMPTY_INT;
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_INT;
	} else {
		return hash[bin]->search(hashstr);
	}
}
float Floattable::search(char* hashstr) {
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_FLOAT;
	} else {
		return hash[bin]->search(hashstr);
	}
}
double Doubletable::search(char* hashstr) {
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_DOUBLE;
	} else {
		return hash[bin]->search(hashstr);
	}
}

void Hashtable::insert(void* obj, char* hashstr) {
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	int add = 1;
	if (hash[bin] == NULL) {
		hash[bin] = new Hashlist(obj,hashstr);
	} else {
		add = hash[bin]->insert(obj,hashstr);
	}
	total+=add;
}
void Inttable::insert(int I, char* hashstr) {
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	int add = 1;
	if (hash[bin] == NULL) {
		hash[bin] = new Intlist(I,hashstr);
	} else {
		add = hash[bin]->insert(I,hashstr);
	}
	total+=add;
}
void LongInttable::insert(long long int I, char* hashstr) {
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	int add = 1;
	if (hash[bin] == NULL) {
		hash[bin] = new LongIntlist(I,hashstr);
	} else {
		add = hash[bin]->insert(I,hashstr);
	}
	total+=add;
}
void Floattable::insert(float I, char* hashstr) {
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	int add = 1;
	if (hash[bin] == NULL) {
		hash[bin] = new Floatlist(I,hashstr);
	} else {
		add = hash[bin]->insert(I,hashstr);
	}
	total+=add;
}
void Doubletable::insert(double I, char* hashstr) {
	unsigned int hashval = hashfunc(hashstr);
	unsigned int bin = hashval % maxsize;
	int add = 1;
	if (hash[bin] == NULL) {
		hash[bin] = new Doublelist(I,hashstr);
	} else {
		add = hash[bin]->insert(I,hashstr);
	}
	total+=add;
}

void* Hashtable::remove(char* str) {
	unsigned int hashval = hashfunc(str);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return NULL;
	} else {
		void* val = NULL;
		val = hash[bin]->remove(str);
		if (val != NULL) total--;
		return val;
	}
}

int Inttable::remove(char* str) {
	unsigned int hashval = hashfunc(str);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_INT;
	} else {
		int val = 0;
		val = hash[bin]->remove(str);
		if (val != EMPTY_INT) total--;
		return val;
	}
}
int LongInttable::remove(char* str) {
	unsigned int hashval = hashfunc(str);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_INT;
	} else {
		int val = 0;
		val = hash[bin]->remove(str);
		if (val != EMPTY_INT) total--;
		return val;
	}
}
float Floattable::remove(char* str) {
	unsigned int hashval = hashfunc(str);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_FLOAT;
	} else {
		float val = 0.0;
		val = hash[bin]->remove(str);
		if (val < EMPTY_FLOAT-1e3 || val > EMPTY_FLOAT+1e3) total--;
		return val;
	}
}
double Doubletable::remove(char* str) {
	unsigned int hashval = hashfunc(str);
	unsigned int bin = hashval % maxsize;
	if (hash[bin] == NULL) {
		return EMPTY_DOUBLE;
	} else {
		double val = 0.0;
		val = hash[bin]->remove(str);
		if (val != EMPTY_DOUBLE) total--;
		return val;
	}
}

char** Hashtable::keys() {
	char** allkeys = new char*[total];
	int index = 0;
	for (int i=0;i<maxsize;i++) {
		if (hash[i] != NULL) {
			for (int j=0;j<hash[i]->length;j++) {
				if (lowMemFlag) {
					allkeys[index++] = hash[i]->hashstrs[j];
				} else {
					copystrhash(allkeys[index++], hash[i]->hashstrs[j]);
				}
			}
		}
	}
	return allkeys;
}
char** Inttable::keys() {
	char** allkeys = new char*[total];
	int index = 0;
	for (int i=0;i<maxsize;i++) {
		if (hash[i] != NULL) {
			for (int j=0;j<hash[i]->length;j++) {
				copystrhash(allkeys[index++], hash[i]->hashstrs[j]);
			}
		}
	}
	return allkeys;
}
char** LongInttable::keys() {
	char** allkeys = new char*[total];
	int index = 0;
	for (int i=0;i<maxsize;i++) {
		if (hash[i] != NULL) {
			for (int j=0;j<hash[i]->length;j++) {
				copystrhash(allkeys[index++], hash[i]->hashstrs[j]);
			}
		}
	}
	return allkeys;
}
char** Floattable::keys() {
	char** allkeys = new char*[total];
	int index = 0;
	for (int i=0;i<maxsize;i++) {
		if (hash[i] != NULL) {
			for (int j=0;j<hash[i]->length;j++) {
				copystrhash(allkeys[index++], hash[i]->hashstrs[j]);
			}
		}
	}
	return allkeys;
}
char** Doubletable::keys() {
	char** allkeys = new char*[total];
	int index = 0;
	for (int i=0;i<maxsize;i++) {
		if (hash[i] != NULL) {
			for (int j=0;j<hash[i]->length;j++) {
				copystrhash(allkeys[index++], hash[i]->hashstrs[j]);
			}
		}
	}
	return allkeys;
}

void copystrhash(char* &dest, char* src) {
    if (src == NULL) {
        dest = NULL;
        return;
    }
    int length = strlen(src)+1;
    dest = new char[length];
    strcpy(dest, src);
}

unsigned int Hashtable::hashfunc(char* str) {
	unsigned int hashNum = 5381;
	int c = *str++;
	while (c) {
		hashNum = ((hashNum << 5) + hashNum*33) + c;
		c=*str++;
	}
	return hashNum;
}

unsigned int Inttable::hashfunc(char* str) {
	unsigned int hashNum = 5381;
	int c = *str++;
	while (c) {
		hashNum = ((hashNum << 5) + hashNum*33) + c;
		c=*str++;
	}
	return hashNum;
}
unsigned int LongInttable::hashfunc(char* str) {
	unsigned int hashNum = 5381;
	int c = *str++;
	while (c) {
		hashNum = ((hashNum << 5) + hashNum*33) + c;
		c=*str++;
	}
	return hashNum;
}

unsigned int Floattable::hashfunc(char* str) {
	unsigned int hashNum = 5381;
	int c = *str++;
	while (c) {
		hashNum = ((hashNum << 5) + hashNum*33) + c;
		c=*str++;
	}
	return hashNum;
}
unsigned int Doubletable::hashfunc(char* str) {
	unsigned int hashNum = 5381;
	int c = *str++;
	while (c) {
		hashNum = ((hashNum << 5) + hashNum*33) + c;
		c=*str++;
	}
	return hashNum;
}


void Floattable::stats() {
	int a=0,b=0,c=0,d=0;
	for (int i=0;i<maxsize;i++ ){
		if (hash[i] == NULL) {
			a++;
		} else {
			if (hash[i]->length < 2) {
				b++;
			} else if (hash[i]->length < 5) {
				c++;
			} else {
				d++;
			}
		}
	}
	fprintf(stderr, "Stats:\nNULL: %d\n1: %d\n2-4: %d\n>4: %d\n",a,b,c,d);
}
