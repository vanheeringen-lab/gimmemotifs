
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


#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <limits.h>


#ifndef HASHTABLE_H
#define HASHTABLE_H


#define EMPTY_INT -1080706050
#define EMPTY_FLOAT -1.23e10
#define EMPTY_FLOAT_CHECK -1.20e10
#define EMPTY_DOUBLE -1.23e100
#define EMPTY_DOUBLE_CHECK -1.20e100

#define DEFAULT_SIZE 100000

void copystrhash(char* &dest, char* src);

class Hashlist {
public:
	int length;
	void** elements;
	char** hashstrs;
	int lowMemFlag;
	Hashlist();
	Hashlist(void* obj, char* str);
	~Hashlist();
	void* search(char*);
	int insert(void* obj, char* str);
	void* remove(char* str);
	char** keys();
};

class Intlist {
public:
	int length;
	int* elements;
	char** hashstrs;
	Intlist();
	Intlist(int I, char* str);
	~Intlist();
	int search(char*);
	int insert(int I, char* str);
	int remove(char* str);
	char** keys();
};

class LongIntlist {
public:
	int length;
	long long int* elements;
	char** hashstrs;
	LongIntlist();
	LongIntlist(long long int I, char* str);
	~LongIntlist();
	long long int search(char*);
	int insert(long long int I, char* str);
	int remove(char* str);
	char** keys();
};

class Floatlist {
public:
	int length;
	float* elements;
	char** hashstrs;
	Floatlist();
	Floatlist(float F, char* str);
	~Floatlist();
	float search(char*);
	int insert(float F, char* str);
	float remove(char* str);
	char** keys();
};

class Doublelist {
public:
	int length;
	double* elements;
	char** hashstrs;
	Doublelist();
	Doublelist(double F, char* str);
	~Doublelist();
	double search(char*);
	int insert(double F, char* str);
	double remove(char* str);
	char** keys();
};


class Hashtable {
public:
	Hashtable();
	Hashtable(int size);
	~Hashtable();
	void init();
	void insert(void*, char*);
	void* search(char*);
	void* remove(char*);
	unsigned int hashfunc(char*);
	int maxsize;
	int total;
	int lowMemFlag;
	void setLowMem(int a);
	Hashlist** hash;
	char** keys();
};

class Inttable {
public:
	Inttable();
	Inttable(int);
	~Inttable();
	void init();
	void insert(int, char*);
	int search(char*);
	int remove(char*);
	int maxsize;
	int total;
	Intlist** hash;
	char** keys();
	unsigned int hashfunc(char*);
};

class LongInttable {
public:
	LongInttable();
	LongInttable(int);
	~LongInttable();
	void init();
	void insert(long long int, char*);
	long long int search(char*);
	int remove(char*);
	int maxsize;
	int total;
	LongIntlist** hash;
	char** keys();
	unsigned int hashfunc(char*);
};

class Floattable {
public:
	Floattable();
	Floattable(int);
	~Floattable();
	void init();
	void insert(float, char*);
	float search(char*);
	float remove(char*);
	int maxsize;
	int total;
	Floatlist** hash;
	char** keys();
	unsigned int hashfunc(char*);
	void stats();
};

class Doubletable {
public:
	Doubletable();
	Doubletable(int);
	~Doubletable();
	void init();
	void insert(double, char*);
	double search(char*);
	double remove(char*);
	int maxsize;
	int total;
	Doublelist** hash;
	char** keys();
	unsigned int hashfunc(char*);
	void stats();
};

class LinkedListItem;
class LinkedList {
public:
	LinkedListItem* firstItem;
	LinkedListItem* lastItem;
	int total;

	LinkedList();
	~LinkedList();
	void add(void*);
	void* get(int);
	void* remove(int);
	void** toArray(int& numberOfItems);
};
class LinkedListItem {
public:
	LinkedListItem* next;
	void* obj;
	LinkedListItem();
};

#endif
