/*
 * Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
 * 
 * This module is free software. You can redistribute it and/or modify it under 
 * the terms of the MIT License, see the file COPYING included with this 
 * distribution.
 *
 * This module contains all the code to compare motifs
 *
 */

#include <Python.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#if PY_MAJOR_VERSION >= 3
    #define PyInt_FromLong PyLong_FromLong
    #define PyString_Check PyUnicode_Check
    #define PyString_Size PyUnicode_GET_LENGTH
#endif


void fill_matrix(double matrix[][4], PyObject *matrix_o) {
	// Parse a PyObject matrix into a C array	
	// Structure of the matrix is as follows: [[freqA, freqC, freqG, freqT], ...]
	int matrix_len = PyList_Size(matrix_o);
	int i,j;
	PyObject *item;
	for (i = 0; i < matrix_len; i++) {
		
		PyObject *row = PyList_GetItem(matrix_o, i);
		if (!PyList_Check(row)) {
			PyErr_SetString( PyExc_TypeError, "wrong matrix structure");
			return;
		}	
		if (PyList_Size(row) != 4) {
			PyErr_SetString( PyExc_TypeError, "4 nucleotides!");
			return;
		}		
		for (j = 0; j < 4; j++) {
			item = PyList_GetItem(row, j);
			matrix[i][j] = PyFloat_AsDouble(item);
		}
	}
}


double mean(double array[], double length) {
	double sum = 0;
	int i;
	for (i = 0; i < length; i++) {
		sum += array[i];
	}
	return sum / length;
}

double sum(double array[], double length) {
	double sum = 0;
	int i;
	for (i = 0; i < length; i++) {
		sum += array[i];
	}
	return sum;
}

double wic(double col1[], double col2[]) {
		int n;
		double x,y;
		double log2 = log(2.0);
		double factor = 2.5;
		double score = 0;;
		double score_a = 0.0;
		double score_b = 0.0;
		double pseudo = 0.0000001;
		double sum1 = 0;	
		double sum2 = 0;	
		
		for (n = 0; n < 4; n++) {
			sum1 += col1[n] + pseudo;
			sum2 += col2[n] + pseudo;
		}

		for (n = 0; n < 4; n++) {
			x = (col1[n] + pseudo) / sum1;
			y = (col2[n] + pseudo) / sum2;
			x = x * log(x / 0.25) / log2;
			y = y * log(y / 0.25) / log2;
			score += fabs(x - y);
			score_a += x;
			score_b += y;
		}
		return sqrt(score_a * score_b) - factor * score;
}

double matrix_wic_mean(double matrix1[][4], double matrix2[][4], int length) {
	// Return the mean pcc of two matrices 
	int i;
	double result[length];

	for (i = 0; i < length; i++ ) {
		result[i] = wic(matrix1[i], matrix2[i]);
	}
	return mean(result, length);
}

double matrix_wic_sum(double matrix1[][4], double matrix2[][4], int length) {
	// Return the mean pcc of two matrices 
	int i;
	double result[length];

	for (i = 0; i < length; i++ ) {
		result[i] = wic(matrix1[i], matrix2[i]);
	}
	return sum(result, length);
}

double pcc(double col1[], double col2[]) {
	// Return the Pearson correlation coefficient of two motif columns
	int n;
	double sum1 = 0;
	double sum2 = 0;
	double mean1, mean2;

	for (n = 0; n < 4; n++) {
		sum1 += col1[n];
		sum2 += col2[n];
	}
	mean1 = sum1 / 4;
	mean2 = sum2 / 4;

	float a = 0;
	float x = 0;
	float y = 0;
	for (n = 0; n < 4; n++) {
		if (((col1[n] - mean1) != 0) && ((col2[n] - mean2) != 0)) {
			a += (col1[n] - mean1) * (col2[n] - mean2);
			x += pow((col1[n] - mean1), 2);
			y += pow((col2[n] - mean2), 2);
		}
		else {
			return 0;
		}
	}

	return  a / sqrt(x * y);
}

double matrix_pcc_mean(double matrix1[][4], double matrix2[][4], int length) {
	// Return the mean pcc of two matrices 
	int i;
	double result[length];

	for (i = 0; i < length; i++ ) {
		result[i] = pcc(matrix1[i], matrix2[i]);
	}
	return mean(result, length);
}

double ed(double col1[], double col2[]) {
	// Return the euclidian distance of two motif columns
	int n;
	double score = 0;

	for (n = 0; n < 4; n++) {
		score = score + pow((col1[n] - col2[n]), 2);
	}

	return  -sqrt(score);
}

double matrix_ed_mean(double matrix1[][4], double matrix2[][4], int length) {
	// Return the mean euclidian distance of two matrices 
	int i;
	double result[length];

	for (i = 0; i < length; i++ ) {
		result[i] = ed(matrix1[i], matrix2[i]);
	}
	return mean(result, length);
}

double distance(double col1[], double col2[]) {
	// Return the distance between two motifs (Harbison et al.) 
	int n;
	double d = 0;

	for (n = 0; n < 4; n++) {
		d = d + pow(col1[n] - col2[n], 2);
	}
	d = d * 1 / sqrt(2);

	return  -d;
}

double matrix_distance_mean(double matrix1[][4], double matrix2[][4], int length) {
	// Return the mean distance (Harbison et al.) of two matrices 
	int i;
	double result[length];

	for (i = 0; i < length; i++ ) {
		result[i] = distance(matrix1[i], matrix2[i]);
	}
	return mean(result, length);
}

static PyObject * c_metrics_score(PyObject *self, PyObject * args)
{
	PyObject *matrix1_o;
	PyObject *matrix2_o;
	const char *metric;
	const char *combine;
	int matrix1_len;
	int matrix2_len;

	if (!PyArg_ParseTuple(args, "OOss", &matrix1_o, &matrix2_o, &metric, &combine)) {
		PyErr_SetString( PyExc_TypeError, "Error parsing arguments");
		return NULL;
	}

	// Retrieve frequency matrix
	if (!PyList_Check(matrix1_o)) {
		PyErr_SetString( PyExc_TypeError, "Error: missing required argument matrix1");
		return NULL;
	}

	if (!PyList_Check(matrix2_o)) {
		PyErr_SetString( PyExc_TypeError, "Error: missing required argument matrix2");
		return NULL;
	}
	// Fill matrix1
	matrix1_len = PyList_Size(matrix1_o);
	double matrix1[matrix1_len][4];
	fill_matrix(matrix1, matrix1_o);
	
	// Fill matrix2
	matrix2_len = PyList_Size(matrix2_o);
	double matrix2[matrix2_len][4];
	fill_matrix(matrix2, matrix2_o);

	if (matrix1_len != matrix2_len) {
		PyErr_SetString( PyExc_TypeError, "Error: matrices have different sizes");
		return NULL;
	}
	int i;
	double result[matrix1_len];
	if (!strcmp(metric, "wic")) {
		for (i = 0; i < matrix1_len; i++ ) {
			result[i] = wic(matrix1[i], matrix2[i]);
		}
	}
	else if (!strcmp(metric, "ed")) {
		for (i = 0; i < matrix1_len; i++ ) {
			result[i] = ed(matrix1[i], matrix2[i]);
		}
	}
	else if (!strcmp(metric, "pcc")) {
		for (i = 0; i < matrix1_len; i++ ) {
			result[i] = pcc(matrix1[i], matrix2[i]);
		}
	}
	else if (!strcmp(metric, "distance")) {
		for (i = 0; i < matrix1_len; i++ ) {
			result[i] = distance(matrix1[i], matrix2[i]);
		}
	}
	else {
			PyErr_SetString( PyExc_TypeError, "Unknown metric");
			return NULL;
	}
	
	if (!strcmp(combine, "mean")) {
		return Py_BuildValue("f", mean(result, matrix1_len));
	}
	else if (!strcmp(combine, "sum")) {
		return Py_BuildValue("f", sum(result, matrix1_len));
	}
	else {
			PyErr_SetString( PyExc_TypeError, "Unknown combine");
			return NULL;
	}

}

int get_truncate_len(len1, len2, pos) {
	// 
	
	if (pos < 0) {
		len2 += pos;
	}
	else if (pos > 0) {
		len1 -= pos;
	}

	if (len1 > len2) {
		return len2;
	}
	else {
		return len1;
	}
	
}

void fill_tmp_matrices(double matrix1[][4], double matrix2[][4], int pos, int l, double tmp_matrix1[][4], double tmp_matrix2[][4]) {
	int start1 = 0;
	int start2 = 0;
	int i, n;

	if (pos < 0) {
		start2 = -pos;
	}
	else {
		start1 = pos;
	}
	
	for (i = 0; i < l; i++) {
		for (n = 0; n < 4; n++) {
			tmp_matrix1[i][n] = matrix1[start1 + i][n];
			tmp_matrix2[i][n] = matrix2[start2 + i][n];
		}
	}

}

int index_at_max(double scores[], int len) {
	double max = scores[0];
	int i_at_max = 0;
	int i;
	
	if (len == 1) {
		return 0;
	}
	
	for (i = 1; i < len; i++) {
		if (scores[i] > max) {
			max = scores[i];
			i_at_max = i;
		}
	}
	return i_at_max;
}

void fill_rc_matrix(double matrix[][4], int len, double rc_matrix[][4]) {
	int i = 0;

	for (i = 0; i < len; i++) {
		rc_matrix[i][0] = matrix[len - i - 1][3];
		rc_matrix[i][1] = matrix[len - i - 1][2];	
		rc_matrix[i][2] = matrix[len - i - 1][1];
		rc_matrix[i][3] = matrix[len - i - 1][0];
	}

}


static PyObject * c_metrics_max_subtotal(PyObject *self, PyObject * args)
{
	PyObject *matrix1_o;
	PyObject *matrix2_o;
	const char *metric;
	const char *combine;
	int matrix1_len;
	int matrix2_len;
	if (!PyArg_ParseTuple(args, "OOss", &matrix1_o, &matrix2_o, &metric, &combine)) {
		PyErr_SetString( PyExc_TypeError, "Error parsing arguments");
		return NULL;
	}

	// Retrieve frequency matrix
	if (!PyList_Check(matrix1_o)) {
		PyErr_SetString( PyExc_TypeError, "Error: missing required argument matrix1");
		return NULL;
	}

	if (!PyList_Check(matrix2_o)) {
		PyErr_SetString( PyExc_TypeError, "Error: missing required argument matrix2");
		return NULL;
	}
	// Fill matrix1
	matrix1_len = PyList_Size(matrix1_o);
	double matrix1[matrix1_len][4];
	fill_matrix(matrix1, matrix1_o);
	
	// Fill matrix2
	matrix2_len = PyList_Size(matrix2_o);
	double matrix2[matrix2_len][4];
	fill_matrix(matrix2, matrix2_o);

	// Fill rc matrix2
	double rc_matrix2[matrix2_len][4];
	fill_rc_matrix(matrix2, matrix2_len, rc_matrix2);

	// Assign correct function
	
	int pos = -2;
	int l;
	float max_score;
	int min_overlap = 6;
	int nr_matches;
	int start_pos = -(matrix2_len - min_overlap);
	nr_matches = matrix1_len + matrix2_len - 2 * min_overlap + 1;
	int i;	
	int positions[nr_matches * 2];
	double scores[nr_matches * 2];
	int orients[nr_matches * 2];

	double (*ptr_metric_function)(double[][4], double[][4], int) = NULL; 
		
	if ((!strcmp(metric, "wic")) && (!strcmp(combine, "mean"))) {
		ptr_metric_function = &matrix_wic_mean;
	}
	else if ((!strcmp(metric, "wic")) && (!strcmp(combine, "sum"))) {
		ptr_metric_function = &matrix_wic_sum;
	}
	else if ((!strcmp(metric, "pcc")) && (!strcmp(combine, "mean"))) {
		ptr_metric_function = &matrix_pcc_mean;
	}
	else if ((!strcmp(metric, "ed")) && (!strcmp(combine, "mean"))) {
		ptr_metric_function = &matrix_ed_mean;
	}
	else if ((!strcmp(metric, "distance")) && (!strcmp(combine, "mean"))) {
		ptr_metric_function = &matrix_distance_mean;
	}
	else {
		PyErr_SetString( PyExc_TypeError, "Unknown metric or combination");
		return NULL;
	}
	
	for (i = 0; i < nr_matches; i++) {
		pos = i + start_pos;
		l = get_truncate_len(matrix1_len, matrix2_len, pos);
		
		// Normal matrix
		double tmp_matrix1[l][4];
		double tmp_matrix2[l][4];

		fill_tmp_matrices(matrix1, matrix2, pos, l, tmp_matrix1, tmp_matrix2);
		max_score = (*ptr_metric_function) (tmp_matrix1, tmp_matrix2, l);
		positions[i] = pos;
		scores[i] = max_score;
		orients[i] = 1;
		//printf("CC 1 Len %i score %f\n", l, max_score);

		// Reverse complement matrix
		fill_tmp_matrices(matrix1, rc_matrix2, pos, l, tmp_matrix1, tmp_matrix2);
		max_score = (*ptr_metric_function) (tmp_matrix1, tmp_matrix2, l);
		positions[i + nr_matches] = pos;
		scores[i + nr_matches] = max_score;
		orients[i+ nr_matches] = -1;

		//printf("CC -1 Len %i score %f\n", l, max_score);
	}
	
	i = index_at_max(scores, nr_matches * 2);
	return Py_BuildValue("fii", scores[i], positions[i], orients[i]);
	
}

static PyObject * c_metrics_pwmscan(PyObject *self, PyObject * args)
{
	
        PyObject *pwm_o;
        PyObject *cutoff_o;
        char *seq;
        int seq_len;
        int n_report;
        int pwm_len;
        int i, j;
        int scan_rc;
        int return_all = 0;

        if (!PyArg_ParseTuple(args, "sOOii|i", &seq, &pwm_o, &cutoff_o, &n_report, &scan_rc, &return_all))
                return NULL;

        seq_len = strlen(seq);

        // Retrieve frequency matrix
        if (!PyList_Check(pwm_o))
                return NULL;

        // Weight matrices
        pwm_len = PyList_Size(pwm_o);
        double pwm[pwm_len][4];
        fill_matrix(pwm, pwm_o);
        //Py_DECREF(pwm_o);

        // Cutoff for every spacer length
        double cutoff;
        cutoff = PyFloat_AsDouble(cutoff_o);

        // Scan sequence
        int j_max = seq_len - pwm_len + 1;
        double score_matrix[j_max];
        double rc_score_matrix[j_max];
        double score, rc_score;

        if (j_max < 0) { j_max = 0;}
        int m;
	int c;
	double pwm_min = -50;	
	for (j = 0; j < j_max; j++) {
		score = 0;
		rc_score = 0;
		for (m = 0; m < pwm_len; m++) {
			switch(seq[j + m]) {
                                case 'A':
                                        score += pwm[m][0];
                                        rc_score += pwm[pwm_len - m - 1][3];
                                        break;
                                case 'C':
                                        score += pwm[m][1];
                                        rc_score += pwm[pwm_len - m - 1][2];
                                        break;
                                case 'G':
                                        score += pwm[m][2];
                                        rc_score += pwm[pwm_len - m - 1][1];
                                        break;
                                case 'T':
                                        score += pwm[m][3];
                                        rc_score += pwm[pwm_len - m - 1][0];
                                        break;
				case 'N':
					pwm_min = pwm[m][0];
					for (c = 1; c < 4; c++) {
						if (pwm[m][c] < pwm_min) {
							pwm_min = pwm[m][c];
						}
					}
					score += pwm_min;
					rc_score += pwm_min; 
					break;
			}
		
		}
		score_matrix[j] = score;
		rc_score_matrix[j] = rc_score;
	}
	
	if (return_all) {
		PyObject *return_list = PyList_New(j_max);

		for (j = 0; j < j_max; j++) {
    			PyList_SetItem(return_list, j, PyFloat_FromDouble(score_matrix[j]));
		}
	    	return return_list;
	}

	// Initialize matrices of n_report highest scores and corresponding positions + strands
	double maxScores[n_report];
	double maxPos[n_report];
	int maxStrand[n_report];
	for (j = 0; j < n_report; j++) {
		maxScores[j] = -100;
		maxPos[j] = -1;
		maxStrand[j] = 1;
	}
	PyObject*  return_list = PyList_New(0);



	int p,q;
	PyObject *x;
	for (j = 0; j < j_max; j++) {
		score = score_matrix[j];
		if (n_report > 0) {
			if (score >= cutoff) {
				p = n_report - 1;
				while ((p >= 0) && (score > maxScores[p])) {
					p--;
				}
				if (p < (n_report-1)) {
					for (q = n_report - 1; q > (p + 1); q--) {
						maxScores[q] = maxScores[q - 1];
						maxPos[q] = maxPos[q - 1];
						maxStrand[q] = maxStrand[q - 1];
					}
					maxScores[p + 1] = score;
					maxPos[p + 1] = j;
					maxStrand[p + 1] = 1;
				}
			}
		}
		else {
			if (score >= cutoff) {
				PyObject* row = PyList_New(0);
				x = PyFloat_FromDouble(score);
				PyList_Append(row, x); Py_DECREF(x);
				x = PyInt_FromLong((long) j);
				PyList_Append(row, x); Py_DECREF(x);
				x = PyInt_FromLong((long) 1);
				PyList_Append(row, x); Py_DECREF(x);
				PyList_Append(return_list, row);
				Py_DECREF(row);
			}
		}
	}

	if (scan_rc) {
		for (j = 0; j < j_max; j++) {
			score = rc_score_matrix[j];
			if (n_report > 0) {
				if (score >= cutoff) {
					p = n_report - 1;
					while ((p >= 0) && (score > maxScores[p])) {
						p--;
					}
					if (p < (n_report-1)) {
						for (q = n_report - 1; q > (p + 1); q--) {
							maxScores[q] = maxScores[q - 1];
							maxPos[q] = maxPos[q - 1];
							maxStrand[q] = maxStrand[q - 1];
						}
						maxScores[p + 1] = score;
						maxPos[p + 1] = j;
						maxStrand[p + 1] = -1;
					}
				}
			}
			else {
				if (score >= cutoff) {
					PyObject* row = PyList_New(0);
					x = PyFloat_FromDouble(score);
					PyList_Append(row, x); Py_DECREF(x);
					x = PyInt_FromLong((long) j);
					PyList_Append(row, x); Py_DECREF(x);
					x =  PyInt_FromLong((long) 1);
					PyList_Append(row, x); Py_DECREF(x);
					PyList_Append(return_list, row);
					Py_DECREF(row);
				}
			}
		}
	}

	for (i = 0; i < n_report; i++) {
		if (maxPos[i] > - 1) {
			PyObject* row = PyList_New(0);
			x = PyFloat_FromDouble(maxScores[i]);
			PyList_Append(row, x);  Py_DECREF(x);
			x = PyInt_FromLong((long)maxPos[i]);
			PyList_Append(row, x); Py_DECREF(x);
			x = PyInt_FromLong((long)maxStrand[i]);
			PyList_Append(row, x); Py_DECREF(x);
			PyList_Append(return_list, row);
			Py_DECREF(row);
		}
	}

	return return_list;
}




static PyObject * c_metrics_pfmscan(PyObject *self, PyObject * args)
{
	
        PyObject *pfm_o;
        PyObject *cutoff_o;
        char *seq;
        int seq_len;
        int n_report;
        int pwm_len;
        int i, j;
        int scan_rc;
        int return_all = 0;

        if (!PyArg_ParseTuple(args, "sOOii|i", &seq, &pfm_o, &cutoff_o, &n_report, &scan_rc, &return_all))
                return NULL;

        seq_len = strlen(seq);

        // Retrieve frequency matrix
        if (!PyList_Check(pfm_o))
                return NULL;

        // Weight matrices
        pwm_len = PyList_Size(pfm_o);
        double pfm[pwm_len][4];
        double pwm[pwm_len][4];
        fill_matrix(pfm, pfm_o);
        //Py_DECREF(pwm_o);

        double g = 0.25;
        double z = 0.01;
        for (i = 0; i < pwm_len; i++) {
                for (j = 0; j < 4; j++) {
                        pwm[i][j] = log(pfm[i][j] / g + z);
                }
        }

        // Cutoff for every spacer length
        double cutoff;
        cutoff = PyFloat_AsDouble(cutoff_o);

        // Scan sequence
        int j_max = seq_len - pwm_len + 1;
        double score_matrix[j_max];
        double rc_score_matrix[j_max];
        double score, rc_score;

        if (j_max < 0) { j_max = 0;}
        int m;
	int c;
	double pwm_min = -50;	
	for (j = 0; j < j_max; j++) {
		score = 0;
		rc_score = 0;
		for (m = 0; m < pwm_len; m++) {
			switch(seq[j + m]) {
                                case 'A':
                                        score += pwm[m][0];
                                        rc_score += pwm[pwm_len - m - 1][3];
                                        break;
                                case 'C':
                                        score += pwm[m][1];
                                        rc_score += pwm[pwm_len - m - 1][2];
                                        break;
                                case 'G':
                                        score += pwm[m][2];
                                        rc_score += pwm[pwm_len - m - 1][1];
                                        break;
                                case 'T':
                                        score += pwm[m][3];
                                        rc_score += pwm[pwm_len - m - 1][0];
                                        break;
				case 'N':
					pwm_min = pwm[m][0];
					for (c = 1; c < 4; c++) {
						if (pwm[m][c] < pwm_min) {
							pwm_min = pwm[m][c];
						}
					}
					score += pwm_min;
					rc_score += pwm_min; 
					break;
			}
		
		}
		score_matrix[j] = score;
		rc_score_matrix[j] = rc_score;
	}
	
	if (return_all) {
		PyObject *return_list = PyList_New(j_max);

		for (j = 0; j < j_max; j++) {
    			PyList_SetItem(return_list, j, PyFloat_FromDouble(score_matrix[j]));
		}
	    	return return_list;
	}

	// Initialize matrices of n_report highest scores and corresponding positions + strands
	double maxScores[n_report];
	double maxPos[n_report];
	int maxStrand[n_report];
	for (j = 0; j < n_report; j++) {
		maxScores[j] = -100;
		maxPos[j] = -1;
		maxStrand[j] = 1;
	}
	PyObject*  return_list = PyList_New(0);



	int p,q;
	PyObject *x;
	for (j = 0; j < j_max; j++) {
		score = score_matrix[j];
		if (n_report > 0) {
			if (score >= cutoff) {
				p = n_report - 1;
				while ((p >= 0) && (score > maxScores[p])) {
					p--;
				}
				if (p < (n_report-1)) {
					for (q = n_report - 1; q > (p + 1); q--) {
						maxScores[q] = maxScores[q - 1];
						maxPos[q] = maxPos[q - 1];
						maxStrand[q] = maxStrand[q - 1];
					}
					maxScores[p + 1] = score;
					maxPos[p + 1] = j;
					maxStrand[p + 1] = 1;
				}
			}
		}
		else {
			if (score >= cutoff) {
				PyObject* row = PyList_New(0);
				x = PyFloat_FromDouble(score);
				PyList_Append(row, x); Py_DECREF(x);
				x = PyInt_FromLong((long) j);
				PyList_Append(row, x); Py_DECREF(x);
				x = PyInt_FromLong((long) 1);
				PyList_Append(row, x); Py_DECREF(x);
				PyList_Append(return_list, row);
				Py_DECREF(row);
			}
		}
	}

	if (scan_rc) {
		for (j = 0; j < j_max; j++) {
			score = rc_score_matrix[j];
			if (n_report > 0) {
				if (score >= cutoff) {
					p = n_report - 1;
					while ((p >= 0) && (score > maxScores[p])) {
						p--;
					}
					if (p < (n_report-1)) {
						for (q = n_report - 1; q > (p + 1); q--) {
							maxScores[q] = maxScores[q - 1];
							maxPos[q] = maxPos[q - 1];
							maxStrand[q] = maxStrand[q - 1];
						}
						maxScores[p + 1] = score;
						maxPos[p + 1] = j;
						maxStrand[p + 1] = -1;
					}
				}
			}
			else {
				if (score >= cutoff) {
					PyObject* row = PyList_New(0);
					x = PyFloat_FromDouble(score);
					PyList_Append(row, x); Py_DECREF(x);
					x = PyInt_FromLong((long) j);
					PyList_Append(row, x); Py_DECREF(x);
					x =  PyInt_FromLong((long) 1);
					PyList_Append(row, x); Py_DECREF(x);
					PyList_Append(return_list, row);
					Py_DECREF(row);
				}
			}
		}
	}

	for (i = 0; i < n_report; i++) {
		if (maxPos[i] > - 1) {
			PyObject* row = PyList_New(0);
			x = PyFloat_FromDouble(maxScores[i]);
			PyList_Append(row, x);  Py_DECREF(x);
			x = PyInt_FromLong((long)maxPos[i]);
			PyList_Append(row, x); Py_DECREF(x);
			x = PyInt_FromLong((long)maxStrand[i]);
			PyList_Append(row, x); Py_DECREF(x);
			PyList_Append(return_list, row);
			Py_DECREF(row);
		}
	}

	return return_list;
}



static PyMethodDef CoreMethods[] = {
	{"score", c_metrics_score, METH_VARARGS,"Test"},
	{"c_max_subtotal", c_metrics_max_subtotal, METH_VARARGS,"Test"},
	{"pfmscan", c_metrics_pfmscan, METH_VARARGS,"Test"},
	{"pwmscan", c_metrics_pwmscan, METH_VARARGS,"Test"},
	{NULL, NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "c_metrics",     /* m_name */
        "This is a module",  /* m_doc */
        -1,                  /* m_size */
        CoreMethods,         /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };
#endif

static PyObject * moduleinit(void) {
    PyObject *m;
#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule3("c_metrics", CoreMethods, "c_metrics_module");
#endif
    return m;
}

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
    PyInit_c_metrics(void)
    {
	return moduleinit();
    };
#else
    initc_metrics(void)
    {
	moduleinit();
    };
#endif

