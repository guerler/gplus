/* Copyright (C) 2010 Aysam Guerler and Ernst-Walter Knapp 

This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details. 

The GNU General Public License details are available <www.gnu.org/licenses/>. 
*/ 

// show detailed information on screen 
//#define LG_ON 

// make PDB-file and PyMol-script 
#define PDB_ON 

// product details 
#define PRODUCT_NAME "GANGSTA+ - Non-Sequential Protein Structure Alignment" 
#define PRODUCT_COPY "Freie Universitaet Berlin (Biochemistry, AG Knapp)" 
#define PRODUCT_AUTHORS "Aysam Guerler (aysam.guerler@gmail.com)" 
#define PRODUCT_VERSION "3.08 - Release 22/09/2010" 

// const 
#define TOLERANCE 0.01 
#define EPSILON 0.00000000001 
#define MINUSEPSILON -0.00000000001 
#define LARGE 999999999 
#define BSIZE 1048576 

// parameters 
class Config 
{ 
	public: 
	// constructor 
	Config() 
	{ 
		// default alignment detection parameters 
		alCoreDistance = 4; // maximal core distance 
		alResidueDistance = 6; // maximal residue assignment cut off 
		alInversion = false; // allow inversion of SSEs
	} 
	
	// alignment 
	int alCoreDistance; 
	double alResidueDistance; 
	bool alInversion;
}; 

// parameters 
class GPlusConfig 
{ 
	
	public: 
	// constructor 
	GPlusConfig() 
	{ 
		evDepth = 5000; // evaluation depth 
		coResults = 200; // number of results to be refined 
		
		// prepare 
		ppCoreDelta = 7; // maximal sse length difference 
		
		// info object parameters 
		inDistMax = 11; 
		inRescale = 5; 
	} 
	
	// general configuration 
	int coResults; 
	
	// info 
	double inDistMax, inRescale; 
	
	// prepare 
	int ppCoreDelta; 
	
	// evaluation 
	int evDepth; 
}; 

// generic string buffer 
char strbuf[2048]; 

// std includes 
using namespace std; 
#include <assert.h> 
#include <unistd.h> 
#include <sys/types.h> 
#include <sys/stat.h> 
#include <time.h> 
#include <map> 
#include <string.h> 
#include <errno.h> 
#include <string> 
#include <list> 
#include <float.h> 
#include <time.h> 
#include <unistd.h> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <iomanip> 
#include <vector> 
#include <stdlib.h> 
#include <cstdarg> 
#include <stdio.h> 
#include <algorithm> 
#include <math.h> 
#include <stddef.h> 
#include <dirent.h> 

// 
// external 
// 

struct Kabsch 
{ 
	
	typedef struct 
	{ 
		double m[4][4]; 
	} MATRIX; 
	
	#define vdiff2(a,b) ( ((a)[0]-(b)[0]) * ((a)[0]-(b)[0]) + ((a)[1]-(b)[1]) * ((a)[1]-(b)[1]) + ((a)[2]-(b)[2]) * ((a)[2]-(b)[2]) ) 
	
	/* 
	calculate rmsd between two structures 
	Params: v1 - first set of points 
	v2 - second set of points 
	N - number of points 
	mtx - return for transfrom matrix used to align structures 
	Returns: rmsd score 
	Notes: mtx can be null. Transform will be rigid. Inputs must 
	be previously aligned for sequence alignment*/ 
	
	double rmsd(double *v1, double *v2, int N, double *mtx) 
	{ 
		double cent1[3]; 
		double cent2[3]; 
		MATRIX tmtx; 
		MATRIX tempmtx; 
		MATRIX move1; 
		MATRIX move2; 
		int i; 
		double answer; 
		double *temp1 = 0; 
		double *temp2 = 0; 
		int err; 
		
		assert(N > 3); 
		
		temp1 = (double*) malloc(N * 3 * sizeof(double)); 
		temp2 = (double*) malloc(N * 3 * sizeof(double)); 
		if(!temp1 || !temp2) 
		goto error_exit; 
		
		centroid(cent1, v1, N); 
		centroid(cent2, v2, N); 
		for(i=0;i<N;i++) 
		{ 
			temp1[i*3+0] = v1[i*3+0] - cent1[0]; 
			temp1[i*3+1] = v1[i*3+1] - cent1[1]; 
			temp1[i*3+2] = v1[i*3+2] - cent1[2]; 
			
			temp2[i*3+0] = v2[i*3+0] - cent2[0]; 
			temp2[i*3+1] = v2[i*3+1] - cent2[1]; 
			temp2[i*3+2] = v2[i*3+2] - cent2[2]; 
		} 
		
		err = getalignmtx(temp1, temp2, N, &tmtx); 
		if(err == -1) 
		goto error_exit; 
		
		mtx_trans(&move1, -cent2[0], -cent2[1], -cent2[2]); 
		mtx_mul(&tempmtx, &move1, &tmtx); 
		mtx_trans(&move2, cent1[0], cent1[1], cent1[2]); 
		mtx_mul(&tmtx, &tempmtx, &move2); 
		memcpy(temp2, v2, N * sizeof(double) * 3); 
		for(i=0;i<N;i++) 
		mulpt(&tmtx, temp2 + i * 3); 
		answer = alignedrmsd(v1, temp2, N); 
		free(temp1); 
		free(temp2); 
		if(mtx) 
		memcpy(mtx, &tmtx.m, 16 * sizeof(double)); 
		
		return answer; 
		error_exit: 
		free(temp1); 
		free(temp2); 
		if(mtx) 
		{ 
			for(i=0;i<16;i++) 
			mtx[i] = 0; 
		} 
		return sqrt(-1.0); 
	} 
	
	/* 
	calculate rmsd between two aligned structures (trivial) 
	Params: v1 - first structure 
	v2 - second structure 
	N - number of points 
	Returns: rmsd 
	*/ 
	static double alignedrmsd(double *v1, double *v2, int N) 
	{ 
		double answer =0; 
		int i; 
		
		for(i=0;i<N;i++) 
		answer += vdiff2(v1 + i *3, v2 + i * 3); 
		return sqrt(answer/N); 
	} 
	
	/* 
	compute the centroid 
	*/ 
	static void centroid(double *ret, double *v, int N) 
	{ 
		int i; 
		
		ret[0] = 0; 
		ret[1] = 0; 
		ret[2] = 0; 
		for(i=0;i<N;i++) 
		{ 
			ret[0] += v[i*3+0]; 
			ret[1] += v[i*3+1]; 
			ret[2] += v[i*3+2]; 
		} 
		ret[0] /= N; 
		ret[1] /= N; 
		ret[2] /= N; 
	} 
	
	/* 
	get the matrix needed to align two structures 
	Params: v1 - reference structure 
	v2 - structure to align 
	N - number of points 
	mtx - return for rigid body alignment matrix 
	Notes: only calculates rotation part of matrix. 
	assumes input has been aligned to centroids 
	*/ 
	static int getalignmtx(double *v1, double *v2, int N, MATRIX *mtx) 
	{ 
		MATRIX A = { {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}} }; 
		MATRIX At; 
		MATRIX Ainv; 
		MATRIX temp; 
		double tv[3]; 
		double tw[3]; 
		double tv2[3]; 
		double tw2[3]; 
		int k, i, j; 
		int flag = 0; 
		double correction; 
		
		correction = absmaxv(v1, N * 3) * absmaxv(v2, N * 3); 
		
		for(k=0;k<N;k++) 
		for(i=0;i<3;i++) 
		for(j=0;j<3;j++) 
		A.m[i][j] += (v1[k*3+i] * v2[k*3+j])/correction; 
		
		while(flag < 3) 
		{ 
			for(i=0;i<4;i++) 
			for(j=0;j<4;j++) 
			At.m[i][j] = A.m[j][i]; 
			
			memcpy(&Ainv, &A, sizeof(MATRIX)); 
			/* this will happen if all points are in a plane */ 
			if( mtx_invert((double *) &Ainv, 4) == -1) 
			{ 
				if(flag == 0) 
				{ 
					crossproduct(tv, v1, v1+3); 
					crossproduct(tw, v2, v2+3); 
				} 
				else 
				{ 
					crossproduct(tv2, tv, v1); 
					crossproduct(tw2, tw, v2); 
					memcpy(tv, tv2, 3 * sizeof(double)); 
					memcpy(tw, tw2, 3 * sizeof(double)); 
				} 
				for(i=0;i<3;i++) 
				for(j=0;j<3;j++) 
				A.m[i][j] += tv[i] * tw[j]; 
				
				flag++; 
			} 
			else 
			flag = 5; 
		} 
		if(flag != 5) 
		return -1; 
		
		mtx_mul(&temp, &At, &A); 
		mtx_root(&temp); 
		mtx_mul(mtx, &temp, &Ainv); 
		return 0; 
	} 
	
	/* 
	get the crossproduct of two vectors. 
	Params: ans - return pinter for answer. 
	pt1 - first vector 
	pt2 - second vector. 
	Notes: crossproduct is at right angles to the two vectors. 
	*/ 
	static void crossproduct(double *ans, double *pt1, double *pt2) 
	{ 
		ans[0] = pt1[1] * pt2[2] - pt1[2] * pt2[1]; 
		ans[1] = pt1[0] * pt2[2] - pt1[2] * pt2[0]; 
		ans[2] = pt1[0] * pt2[1] - pt1[1] * pt2[0]; 
	} 
	
	/* 
	Denman-Beavers square root iteration 
	*/ 
	static void mtx_root(MATRIX *mtx, int limit = 1000) 
	{ 
		MATRIX Y = *mtx; 
		MATRIX Z; 
		MATRIX Y1; 
		MATRIX Z1; 
		MATRIX invY; 
		MATRIX invZ; 
		MATRIX Y2; 
		int iter = 0; 
		int i, ii; 
		
		mtx_identity(&Z); 
		
		do 
		{ 
			invY = Y; 
			invZ = Z; 
			if( mtx_invert((double *) &invY, 4) == -1) 
			return; 
			if( mtx_invert((double *) &invZ, 4) == -1) 
			return; 
			for(i=0;i<4;i++) 
			for(ii=0;ii<4;ii++) 
			{ 
				Y1.m[i][ii] = 0.5 * (Y.m[i][ii] + invZ.m[i][ii]); 
				Z1.m[i][ii] = 0.5 * (Z.m[i][ii] + invY.m[i][ii]); 
			} 
			Y = Y1; 
			Z = Z1; 
			
			mtx_mul(&Y2, &Y, &Y); 
		} 
		while(!almostequal(&Y2, mtx) && iter++ < limit); 
		
		*mtx = Y; 
	} 
	
	/* 
	Check two matrices for near-enough equality 
	Params: a - first matrix 
	b - second matrix 
	Returns: 1 if almost equal, else 0, epsilon 0.0001f. 
	*/ 
	static int almostequal(MATRIX *a, MATRIX *b) 
	{ 
		int i, ii; 
		double d = 0.0; 
		
		for(i=0;i<4;i++) 
		for(ii=0;ii<4;ii++) 
		d += fabs(a->m[i][ii] - b->m[i][ii]); 
		
		if (d < 0.00000000001) 
		return 1; 
		else 
		return 0; 
	} 
	
	/*static int almostequal(MATRIX *a, MATRIX *b) 
	{ 
		int i, ii; 
		double epsilon = 0.001f; 
		
		for(i=0;i<4;i++) 
		for(ii=0;ii<4;ii++) 
		if(fabs(a->m[i][ii] - b->m[i][ii]) > epsilon) 
		return 0; 
		return 1; 
	}*/ 
	
	
	/* 
	multiply a point by a matrix. 
	Params: mtx - matrix 
	pt - the point (transformed) 
	*/ 
	static void mulpt(MATRIX *mtx, double *pt) 
	{ 
		double ans[4] = {0}; 
		int i; 
		int ii; 
		
		for(i=0;i<4;i++) 
		{ 
			for(ii=0;ii<3;ii++) 
			{ 
				ans[i] += pt[ii] * mtx->m[ii][i]; 
			} 
			ans[i] += mtx->m[3][i]; 
		} 
		pt[0] = ans[0]; 
		pt[1] = ans[1]; 
		pt[2] = ans[2]; 
	} 
	
	/* 
	multiply two matrices. 
	Params: ans - return pointer for answer. 
	x - first matrix 
	y - second matrix. 
	Notes: ans may not be equal to x or y. 
	*/ 
	static void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y) 
	{ 
		int i; 
		int ii; 
		int iii; 
		
		for(i=0;i<4;i++) 
		for(ii=0;ii<4;ii++) 
		{ 
			ans->m[i][ii] = 0; 
			for(iii=0;iii<4;iii++) 
			ans->m[i][ii] += x->m[i][iii] * y->m[iii][ii]; 
		} 
	} 
	
	
	/* 
	create an identity matrix. 
	Params: mtx - return pointer. 
	*/ 
	static void mtx_identity(MATRIX *mtx) 
	{ 
		int i; 
		int ii; 
		
		for(i=0;i<4;i++) 
		for(ii=0;ii<4;ii++) 
		{ 
			if(i==ii) 
			mtx->m[i][ii] = 1.0f; 
			else 
			mtx->m[i][ii] = 0; 
		} 
	} 
	
	/* 
	create a translation matrix. 
	Params: mtx - return pointer for matrix. 
	x - x translation. 
	y - y translation. 
	z - z translation 
	*/ 
	static void mtx_trans(MATRIX *mtx, double x, double y, double z) 
	{ 
		mtx->m[0][0] = 1; 
		mtx->m[0][1] = 0; 
		mtx->m[0][2] = 0; 
		mtx->m[0][3] = 0; 
		
		mtx->m[1][0] = 0; 
		mtx->m[1][1] = 1; 
		mtx->m[1][2] = 0; 
		mtx->m[1][3] = 0; 
		
		mtx->m[2][0] = 0; 
		mtx->m[2][1] = 0; 
		mtx->m[2][2] = 1; 
		mtx->m[2][3] = 0; 
		
		mtx->m[3][0] = x; 
		mtx->m[3][1] = y; 
		mtx->m[3][2] = z; 
		mtx->m[3][3] = 1; 
	} 
	
	/* 
	matrix invert routine 
	Params: mtx - the matrix in raw format, in/out 
	N - width and height 
	Returns: 0 on success, -1 on fail 
	*/ 
	static int mtx_invert(double *mtx, int N) 
	{ 
		int indxc[100]; /* these 100s are the only restriction on matrix size */ 
		int indxr[100]; 
		int ipiv[100]; 
		int i, j, k; 
		int irow = 0, icol = 0; 
		double big; 
		double pinv; 
		int l, ll; 
		double dum; 
		double temp; 
		
		assert(N <= 100); 
		
		for(i=0;i<N;i++) 
		ipiv[i] = 0; 
		
		for(i=0;i<N;i++) 
		{ 
			big = 0.0; 
			
			/* find biggest element */ 
			for(j=0;j<N;j++) 
			if(ipiv[j] != 1) 
			for(k=0;k<N;k++) 
			if(ipiv[k] == 0) 
			if(fabs(mtx[j*N+k]) >= big) 
			{ 
				big = fabs(mtx[j*N+k]); 
				irow = j; 
				icol = k; 
			} 
			
			ipiv[icol]=1; 
			
			if(irow != icol) 
			for(l=0;l<N;l++) 
			{ 
				temp = mtx[irow * N + l]; 
				mtx[irow * N + l] = mtx[icol * N + l]; 
				mtx[icol * N + l] = temp; 
			} 
			
			indxr[i] = irow; 
			indxc[i] = icol; 
			
			
			/* if biggest element is zero matrix is singular, bail */ 
			if(mtx[icol* N + icol] == 0) 
			goto error_exit; 
			
			pinv = 1.0/mtx[icol * N + icol]; 
			
			mtx[icol * N + icol] = 1.0; 
			
			for(l=0;l<N;l++) 
			mtx[icol * N + l] *= pinv; 
			
			for(ll=0;ll<N;ll++) 
			if(ll != icol) 
			{ 
				dum = mtx[ll * N + icol]; 
				mtx[ll * N + icol] = 0.0; 
				for(l=0;l<N;l++) 
				mtx[ll * N + l] -= mtx[icol * N + l]*dum; 
			} 
		} 
		
		
		/* unscramble matrix */ 
		for (l=N-1;l>=0;l--) 
		{ 
			if (indxr[l] != indxc[l]) 
			for (k=0;k<N;k++) 
			{ 
				temp = mtx[k * N + indxr[l]]; 
				mtx[k * N + indxr[l]] = mtx[k * N + indxc[l]]; 
				mtx[k * N + indxc[l]] = temp; 
			} 
		} 
		
		return 0; 
		
		error_exit: 
		return -1; 
	} 
	
	/* 
	get the asolute maximum of an array 
	*/ 
	static double absmaxv(double *v, int N) 
	{ 
		double answer = 0.0; 
		int i; 
		
		for(i=0;i<N;i++) 
		if(answer < fabs(v[i])) 
		answer = fabs(v[i]); 
		return answer; 
	} 
}; 

// 
// internal 
// 

// general library 
// vector 
template < typename T > 
struct Vec : public vector<T> 
{ 
	// constructor 
	Vec() {} 
	
	// fill 
	Vec(int n) 
	{ 
		vector<T>::resize(n); 
	} 
	
	// fill 
	inline void fill(T elem) 
	{ 
		for (int i = 0; i < vector<T>::size(); i++) 
		vector<T>::operator[](i) = elem; 
	} 
	
	inline T* ptr() 
	{ 
		return &(*this)[0]; 
	} 
	
	T operator[] (int i) const 
	{ 
		assert (i < vector<T>::size()); 
		return vector<T>::operator[](i); 
	} 
	
	T& operator[] (int i) 
	{ 
		//assert (i < vector<T>::size()); 
		return vector<T>::operator[](i); 
	} 
	
	inline void append(T element) 
	{ 
		vector<T>::push_back(element); 
	} 
}; 

// vector for bool 
template < > 
struct Vec < bool > : public vector<bool> 
{ 
	// fill 
	inline void fill(bool elem) 
	{ 
		for (unsigned int i = 0; i < vector<bool>::size(); i++) 
		vector<bool>::operator[](i) = elem; 
	} 
}; 

// default vector class 
struct Vec3 
{ 
	// 
	// data 
	// 
	double x, y, z; 
	
	// 
	// construction 
	// 
	Vec3(double x, double y, double z) 
	{ 
		initialize (x, y, z); 
	} 
	
	Vec3(double x = 0.0) 
	{ 
		initialize (x, x, x); 
	} 
	
	int size() 
	{ 
		return 3; 
	} 
	
	// 
	// operators 
	// 
	Vec3 operator+(const Vec3 &v) const 
	{ 
		Vec3 v0; 
		v0.x = x + v.x; 
		v0.y = y + v.y; 
		v0.z = z + v.z; 
		return v0; 
	} 
	
	Vec3 operator-(const Vec3& v) const 
	{ 
		Vec3 v0; 
		v0.x = x - v.x; 
		v0.y = y - v.y; 
		v0.z = z - v.z; 
		return v0; 
	} 
	
	void operator+=(const Vec3 &v) 
	{ 
		x += v.x; 
		y += v.y; 
		z += v.z; 
	} 
	
	void operator-=(const Vec3 &v) 
	{ 
		x -= v.x; 
		y -= v.y; 
		z -= v.z; 
	} 
	
	void operator/=(const double &d) 
	{ 
		x /= d; 
		y /= d; 
		z /= d; 
	} 
	
	void operator*=(const double &d) 
	{ 
		x *= d; 
		y *= d; 
		z *= d; 
	} 
	
	void operator-=(const double &d) 
	{ 
		x -= d; 
		y -= d; 
		z -= d; 
	} 
	
	void operator+=(const double &d) 
	{ 
		x += d; 
		y += d; 
		z += d; 
	} 
	
	void operator*=(const Vec3 &v) 
	{ 
		x = y * v.z - z * v.y; 
		y = z * v.x - x * v.z; 
		z = x * v.y - y * v.x; 
	} 
	
	
	double operator*(const Vec3& v) const 
	{ 
		return dot(v); 
	} 
	
	Vec3 operator*(const double& d) const 
	{ 
		return Vec3 (x*d, y*d, z*d); 
	} 
	
	bool operator==(const Vec3 &v) const 
	{ 
		return (length() == v.length()); 
	} 
	
	bool operator==(const double &d) const 
	{ 
		return (x == d && y == d && z == d); 
	} 
	
	bool operator!=(const double &d) const 
	{ 
		return !operator==(d); 
	} 
	
	// unary minus 
	friend Vec3 operator-(const Vec3& v) 
	{ 
		return Vec3(-v.x, -v.y, -v.z); 
	} 
	
	// inefficient for compatibility only 
	inline double& operator[](int i) 
	{ 
		assert(i < 3); 
		assert(i >= 0); 
		return (&x)[i]; 
	} 
	
	// 
	// methods 
	// 
	
	// distance 
	double dist(const Vec3& v) const 
	{ 
		return (Vec3 ((*this) - v).length()); 
	} 
	
	// dot product 
	double dot(const Vec3& v) const 
	{ 
		return (x*v.x + y*v.y + z*v.z); 
	} 
	
	// cross product 
	Vec3 cross(const Vec3& v) const 
	{ 
		return Vec3 (y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); 
	} 
	
	// one cross product self = cross(v1, v2) 
	template < typename TA, typename TB > 
	void cross(TA& v1, TB& v2) 
	{ 
		x= v1.y*v2.z-v2.y*v1.z; 
		y= v2.x*v1.z-v1.x*v2.z; 
		z= v1.x*v2.y-v2.x*v1.y; 
	} 
	
	// multiplying a cross product by a scalar is very common 
	// one cross product v3 = k*cross(v1, v2) 
	template < typename TA, typename TB > 
	static Vec3 cross (double k, TA& v1, TB& v2) 
	{ 
		return Vec3( k*(v1.y*v2.z-v2.y*v1.z), 
		k*(v2.x*v1.z-v1.x*v2.z), 
		k*(v1.x*v2.y-v2.x*v1.y) ); 
	} 
	
	// cross 
	void neg() 
	{ 
		x = -x; 
		y = -y; 
		z = -z; 
	} 
	
	int equals(Vec3& v) 
	{ 
		return ( (x == v.x) && (y == v.y) && (z == v.z) ); 
	} 
	
	// set all 
	inline void set(double x, double y, double z) 
	{ 
		this->x = x; 
		this->y = y; 
		this->z = z; 
	} 
	
	// get length 
	double length() const 
	{ 
		double factor = x*x + y*y + z*z; 
		if (factor > 0) 
		return sqrt(factor); 
		else 
		return 0.0; 
	} 
	
	// is normalized 
	bool isNormalized() const 
	{ 
		return (length() == 1); 
	} 
	
	// normalize 
	void normalize() 
	{ 
		(*this) /= length(); 
	} 
	
	// scalar 
	double scalar(const Vec3 &v) const 
	{ 
		return (x * v.x) + (y * v.y) + (z * v.z); 
	} 
	
	
	private: 
	// initialize 
	void initialize (double x, double y, double z) 
	{ 
		this->x = x; 
		this->y = y; 
		this->z = z; 
	} 
}; 

// output 
ostream& operator<<(ostream& stream, const Vec3& v) 
{ 
	stream << setw(15) << v.x << setw(15) << v.y << setw(15) << v.z; 
	return stream; 
} 


// 4d vector 
struct Vec4 
{ 
	// values 
	double w, x, y, z; 
	
	// construct 
	Vec4() 
	{ 
		initialize(0, 0, 0, 0); 
	} 
	
	Vec4(double w, double x, double y, double z) 
	{ 
		initialize(w, x, y, z); 
	} 
	
	// operator 
	void operator*=(const double &d) 
	{ 
		w *= d; 
		x *= d; 
		y *= d; 
		z *= d; 
	} 
	
	// operator 
	void operator/=(const double &d) 
	{ 
		w /= d; 
		x /= d; 
		y /= d; 
		z /= d; 
	} 
	
	// get length 
	double length() const 
	{ 
		double factor = w*w + x*x + y*y + z*z; 
		if (factor > 0) 
		return sqrt(factor); 
		else 
		return 0; 
	} 
	
	// inefficient for compatibility only 
	inline double& operator[](int i) 
	{ 
		assert(i < 4); 
		assert(i >= 0); 
		return (&w)[i]; 
	} 
	private: 
	void initialize(double w, double x, double y, double z) 
	{ 
		this->w = w; 
		this->x = x; 
		this->y = y; 
		this->z = z; 
	} 
}; 

// log file 
struct Msg 
{ 
	// write title 
	static void title() 
	{ 
		#ifdef LG_ON 
		cout << PRODUCT_NAME << endl; 
		cout << PRODUCT_COPY << endl; 
		cout << PRODUCT_AUTHORS << endl; 
		next(); 
		#endif 
	} 
	
	static void line() 
	{ 
		#ifdef LG_ON 
		cout << endl; 
		#endif 
	} 
	
	static void write(const char* pszFormat, ...) 
	{ 
		#ifdef LG_ON 
		// buffer and timer 
		static char buf[2048]; 
		
		// write the formated log string to szLine 
		va_list argList; 
		va_start( argList, pszFormat ); 
		vsprintf( buf, pszFormat, argList ); 
		va_end( argList ); 
		cout << buf << endl; 
		#endif 
	} 
	
	static void rewrite(const char* pszFormat, ...) 
	{ 
		#ifdef LG_ON 
		// buffer and timer 
		static char buf[2048]; 
		
		// write the formated log string to szLine 
		va_list argList; 
		va_start( argList, pszFormat ); 
		vsprintf( buf, pszFormat, argList ); 
		va_end( argList ); 
		cout << buf << "          \r" << flush; 
		#endif 
	} 
	
	// error handle 
	template < typename T > 
	static void error (const char* title, T msg) 
	{ 
		cout << endl << "================================================================================" << endl << endl; 
		cout << "Error occured : " << title << endl; 
		cout << endl << "================================================================================" << endl << endl; 
		cout << " " << msg; 
		cout << endl << "================================================================================" << endl << endl; 
		
		// exit 
		exit(0); 
	} 
	
	static void next() 
	{ 
		#ifdef LG_ON 
		cout << endl << "================================================================================" << endl << endl; 
		#endif 
	} 
}; 

// random value generator 
struct Random 
{ 
	inline static void init(unsigned int seed) 
	{ 
		srand (seed); 
	} 
	
	inline static void init() 
	{ 
		srand (time(NULL)); 
	} 
	
	// random 
	static double get () 
	{ 
		return rand() / double (RAND_MAX); 
	} 
	
	inline static int get(int max) 
	{ 
		return (int) (max * get()); 
	} 
	
	inline static double get(double max) 
	{ 
		return max * get(); 
	} 
	
	inline static double sign() 
	{ 
		if (get() > 0.5) 
		return 1.0; 
		else 
		return -1.0; 
	} 
	
	// random pi 
	template < typename TVec3 > 
	static void get (TVec3& vec) 
	{ 
		// get 
		for (int i = 0; i < 3; i++) 
		vec[i] = Random::get(1.0) * Random::sign(); 
		
		// unit vector 
		vec /= vec.length(); 
	} 
}; 


// local library 
struct Convert 
{ 
	// int to string 
	inline static string toString (int i) 
	{ 
		static char buffer[16]; 
		snprintf(buffer, sizeof(buffer), "%d", i); 
		return string(buffer); 
	} 
	
	// anything to string 
	template < typename TAny > 
	inline static string toString(TAny t) 
	{ 
		std::ostringstream sstrm; 
		sstrm << t; 
		return sstrm.str(); 
	} 
	
	// format 
	inline static double toDbl(string str) 
	{ 
		return atof(str.c_str()); 
	} 
	
	inline static char toChar (string str) 
	{ 
		return str.c_str()[0]; 
	} 
	
	inline static int toInt (string str) 
	{ 
		return atoi(str.c_str()); 
	} 
}; 

// local library 
struct FName 
{ 
	// returns path/path/name.typ.typ -> name.typ.typ 
	static string getfullname(string path) 
	{ 
		return path.substr(path.rfind("/") + 1); 
	} 
	
	// returns path/path/name.typ.typ -> name 
	static string getname(string path) 
	{ 
		// find 
		int start = path.rfind("/") + 1; 
		if (start != -1) 
		path = path.substr(start); 
		
		// find 
		int end = path.find("."); 
		if (end != -1) 
		path = path.substr(0, end); 
		
		// return 
		return path; 
	} 
	
	static string reduceprefix(string path) 
	{ 
		return path.substr(path.find("_") + 1); 
	} 
	
	// returns name.typ.typ -> name.typ 
	static string reduce(string fname) 
	{ 
		int pos = fname.rfind("."); 
		if (pos != -1) 
		return fname.substr(0, fname.rfind(".")); 
		else 
		return fname; 
	} 
	
	// returns name.typ.typ -> .typ 
	static string gettype(string fname) 
	{ 
		int pos = fname.rfind("."); 
		if (pos != -1) 
		return fname.substr(pos); 
		else 
		return ""; 
	} 
}; 

// local library 
struct Lib 
{ 
	// str clean
	static string strclean(string str)
	{
        string output = "";
		for (int i = 0; i < str.length(); i++)
		{
			if (str[i] != ' ' && str[i] != '\n' && str[i] != '\r' && str[i] != '\t' && str[i] != 0)
				output += str[i];
		} 
		return output;
	}
	
	/** 
	tokenizer 
	*/ 
	template < typename TChar > 
	static void split(Vec < string >& data, TChar* ptr) 
	{ 
		// reset 
		data.clear(); 
		string str = ""; 
		
		// loop 
		bool valid = false; 
		while (*ptr != '\0') 
		{ 
			// add data 
			if (*ptr == '\n' || *ptr == ' ' || *ptr == '\r' || *ptr == '\t') 
			{ 
				if (valid) 
				{ 
					data.push_back(str); 
					str = ""; 
				} 
				
				valid = false; 
			} else { 
				str += *ptr; 
				valid = true; 
			} 
			
			// next character 
			ptr++; 
		} 
		
		// add final item 
		if (str != "") 
		data.push_back(str); 
	} 
	
	// bubble sort 
	template < typename TA, typename TB > 
	static void sort (Vec < TA >& data, Vec < TB >& index, int size = 0) 
	{ 
		// get size 
		int n = size; 
		if (n == 0) 
		n = data.size(); 
		if (index.size() < n) 
		Msg::error("Lib::sort()", "Size of index array is too small!"); 
		
		// sort data 
		TA temp; 
		TB tint; 
		for (int i = 0; i < n-1; i++) 
		for (int j = i+1; j < n; j++) 
		{ 
			if (data[i] > data[j]) 
			{ 
				temp = data[i]; 
				data[i] = data[j]; 
				data[j] = temp; 
				
				tint = index[i]; 
				index[i] = index[j]; 
				index[j] = tint; 
			} 
		} 
	} 
	
	// shuffle 
	template < typename TArray > 
	static void shuffle (TArray& v) 
	{ 
		// parameters 
		int value, index; 
		
		// loop 
		for (int i = 0; i < v.size(); i++) 
		{ 
			// get random index 
			index = Random::get((int) v.size() - 1); 
			
			// swap 
			value = v[i]; 
			v[i] = v[index]; 
			v[index] = value; 
		} 
	} 
	
	// sign 
	template < typename TItem > 
	inline static int sign (TItem f) 
	{ 
		// setup sign 
		if (f > 0) 
		return 1; 
		else 
		if (f < 0) 
		return -1; 
		else 
		return 0; 
	} 
	
	// square 
	inline static double sqr(double val) 
	{ 
		return pow(val, 2); 
	} 
}; 

// parameters 
char buffer[BSIZE]; 

// file handle 
class File 
{ 
	// file stream 
	ifstream* ifs; 
	
	// maxsize 
	int ndata; 
	public: 
	// data 
	Vec < string > data; 
	
	// real construct 
	File() 
	{ 
		ifs = 0; 
		ndata = 0; 
	} 
	
	// construct 
	bool construct(string filename) 
	{ 
		return open (filename); 
	} 
	
	// open file 
	bool open (string filename) 
	{ 
		// file stream 
		if (ifs != 0) 
		delete ifs; 
		
		// open stream 
		ifs = new ifstream; 
		
		// open filestream 
		ifs->open (filename.c_str(), ios::in); 
		if (!ifs->good()) 
		return false; 
		else 
		return true; 
	} 
	
	// destructor 
	~File() 
	{ 
		close(); 
	} 
	
	// close 
	void close() 
	{ 
		if (ifs != 0) 
		{ 
			delete ifs; 
			ifs = 0; 
		} 
	} 
	
	// open file 
	static bool exists (string filename) 
	{ 
		// open file 
		ifstream ifs; 
		ifs.open (filename.c_str(), ios::in); 
		if (ifs.good()) 
		return true; 
		else 
		return false; 
	} 
	
	// good 
	inline bool good() 
	{ 
		return ifs->good(); 
	} 
	
	// default get 
	inline string get(int i = 0) 
	{ 
		if (i < ndata) 
		return data[i]; 
		else 
		return ""; 
	} 
	
	inline string get(int i, int j) 
	{ 
		if (i >= 0 && ndata > 0) 
		return string(data[0], i, min (j, ndata - i)); 
		else 
		return ""; 
	} 
	
	// get 
	inline char getChar (int i = 0) 
	{ 
		return Convert::toChar(get(i)); 
	} 
	
	inline int getInt (int i = 0) 
	{ 
		return Convert::toInt(get(i)); 
	} 
	
	inline double getDbl(int i = 0) 
	{ 
		return Convert::toDbl(get(i)); 
	} 
	
	inline double getDouble(int i = 0) 
	{ 
		return getDbl(i); 
	} 
	
	// reader 
	inline bool readLine() 
	{ 
		// read 
		ifs->getline(buffer, BSIZE, '\n'); 
		data.resize(1); 
		data[0] = string (buffer); 
		ndata = data[0].length(); 
		
		// return 
		if (ndata > 0) 
		return true; 
		else 
		return false; 
	} 
	
	inline bool read() 
	{ 
		// load 
		ifs->getline(buffer, BSIZE, '\n'); 
		
		// split 
		Lib::split (data, buffer); 
		ndata = data.size(); 
		
		// return 
		if (ndata > 0) 
		return true; 
		else 
		return false; 
	} 
	
	// return size of entries 
	inline int size() 
	{ 
		return ndata; 
	} 
}; 

// secondary structure element 
struct SSE 
{ 
	int start, end; 
	int n; 
	char type; 
	
	// copy constructor 
	SSE(SSE *sse) 
	{ 
		this->start = sse->start; 
		this->end = sse->end; 
		this->type = sse->type; 
		this->n = sse->n; 
	} 
	
	SSE(char type, int start = 0, int end = 0) 
	{ 
		this->start = start; 
		this->end = end; 
		this->type = type; 
		n = end - start + 1; 
	} 
	
	// recieve 
	inline int length() const { return n; } 
	inline int getStart() const { return start; } 
	inline int getEnd() const { return end; } 
	inline char getType() const { return type; } 
	
	// alter 
	void setStart(int start) 
	{ 
		this->start = start; 
		this->n = end - start + 1; 
	} 
	
	void setEnd(int end) 
	{ 
		this->end = end; 
		this->n = end - start + 1; 
	} 
}; 



struct SpecDetails 
{ 
	// types 
	static const int AMINO = 20; 
	static const int ATOMS = 37; 
	
	// setup atom 
	inline static int getResCode (string res) 
	{ 
		// set amino acid code sorted by single code 
		if (res == "ALA") return 0; // 'A' 
		else if (res == "CYS") return 1; // 'C' 
		else if (res == "ASP") return 2; // 'D' 
		else if (res == "GLU") return 3; // 'E' 
		else if (res == "PHE") return 4; // 'F' 
		else if (res == "GLY") return 5; // 'G' 
		else if (res == "HIS") return 6; // 'H' 
		else if (res == "ILE") return 7; // 'I' 
		else if (res == "LYS") return 8; // 'K' 
		else if (res == "LEU") return 9; // 'L' 
		else if (res == "MET") return 10; // 'M' 
		else if (res == "ASN") return 11; // 'N' 
		else if (res == "PRO") return 12; // 'P' 
		else if (res == "GLN") return 13; // 'Q' 
		else if (res == "ARG") return 14; // 'R' 
		else if (res == "SER") return 15; // 'S' 
		else if (res == "THR") return 16; // 'T' 
		else if (res == "VAL") return 17; // 'V' 
		else if (res == "TRP") return 18; // 'W' 
		else if (res == "TYR") return 19; // 'Y' 
		return 20; 
	} 
}; 

// atom
struct SpecAtom
{   
    // const
    static const char SSECOIL   = '?';
	static const int  SSEEMPTY  = -1;
	
	// parameters
	SpecDetails	mdetail;
	
	// data
	Vec3       pos;
	string     name;
	string     res;
	int        rescode;	
	int        resno;
	int        respdb;
	int        sseno;	
	char       sse;

    // additional pdb
	double     occupancy;
    double     bfactor;
    
	// construct
	SpecAtom(Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, char sse = SSECOIL, int sseno = SSEEMPTY)
	{
	        set (pos, name, res, resno, respdb, occupancy, bfactor, sse, sseno);      
	}
    
	// setup atom
	inline void set (Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, char sse = SSECOIL, int sseno = SSEEMPTY)
	{
		// data
		this->pos        = pos;
		this->name       = Lib::strclean(name);
		this->res        = res;
		this->resno      = resno;
		this->occupancy  = occupancy;
		this->bfactor    = bfactor;

        // original residue identifier from PDB
        if (respdb == LARGE)
            this->respdb = resno;
		else
            this->respdb = respdb;
	
        // secondary information
		this->sse        = sse;
		this->sseno      = sseno;

        // get residue code
		this->rescode    = mdetail.getResCode(res);
	}          
};

// molecule object
class SpecMolecule
{
	// parameters
	int 		n;		
	public:
	
	// molecule name or pdbcode
	string 		name, info;

	// lists
	Vec < SpecAtom >   latom;
	Vec < SpecAtom >   lhetatom;
	Vec < SSE  >       lsse;
	Vec < int >        lcalpha;
	Vec < char >       lssebyres;
	Vec < int >        lssecenter;

	// tables
	Vec < Vec < SpecAtom* > > tresidue;

	// type tags
	static const char SSECOIL   = '?';
	static const char SSEALPHA  = 'H';
	static const char SSEBETA   = 'S';
	static const char SSETURN   = 'T';	
	static const int  SSEEMPTY  = -1;
	string ATOMCA;

	// construct
	SpecMolecule()
	{
		construct();
	}

	// construction
	void construct (string name = "", string info = "")
	{
		// initialize
		initialize (info, name);
	}

	// open initialize
	void construct(SpecMolecule* mol, string name = "")
	{
		// change name
		if (name == "")
			name = mol->name;
			
		// initialize
		initialize (mol->info, name);

		// setup size and center
		n = mol->size();
	
		// copy lists
		latom.insert(latom.begin(), mol->latom.begin(), mol->latom.end());
		lhetatom.insert(lhetatom.begin(), mol->lhetatom.begin(), mol->lhetatom.end());
		lsse.insert(lsse.begin(), mol->lsse.begin(), mol->lsse.end() );
		lcalpha.insert(lcalpha.begin(), mol->lcalpha.begin(), mol->lcalpha.end() );
		lssebyres.insert(lssebyres.begin(), mol->lssebyres.begin(), mol->lssebyres.end() );
		lssecenter.insert(lssecenter.begin(), mol->lssecenter.begin(), mol->lssecenter.end() );
	}     
	
	// copy coordinates
	void readCoordinates(SpecMolecule* mol)
	{
		// assert
		assert (n == mol->n);
		
		// copy
		for (int i = 0; i < n; i++)
		    latom[i].pos = mol->latom[i].pos;
	}

	// copy coordinates
	void readCAlpha(SpecMolecule* mol)
	{
		// assert
		assert (lcalpha.size() == mol->lcalpha.size());
		
		// copy
		for (int i = 0; i < lcalpha.size(); i++)
		    latom[lcalpha[i]].pos = mol->latom[mol->lcalpha[i]].pos;
	}

	// load molecule
	void readStructure(SpecMolecule* mol)
	{	
		// initialize
		int cres = 0;
		if (latom.size() > 0)
			cres = latom[latom.size() - 1].resno;

		// add atoms		
		int nres = 0;        
		for (int i = 0; i < mol->latom.size(); i++)
		{
		    // atom
		    SpecAtom* a = &mol->latom[i];

		    // check
		    if(nres != a->resno)
		    {
                nres = a->resno;
                cres++;
		    }

		    // append
		    append(a->pos, a->name, a->res, cres);
		}

		// prepare
		lssebyres.insert(lssebyres.end(), mol->lssebyres.begin(), mol->lssebyres.end() );
	}

	// copy calphas
	void copyCAlpha(SpecMolecule* selection)
	{
		// setup molecule
		selection->construct(name);
	
		// setup atoms
		SpecAtom* a;
		for (int i = 0; i < lcalpha.size(); i++)
		{
		    a = &latom[lcalpha[i]];
		    selection->append (a->pos, a->name, a->res, a->resno, a->respdb, a->occupancy, a->bfactor, a->sse);
		}

		// update
		selection->finalize();
	}

	// copy calphas
	void copyBone(SpecMolecule* selection)
	{
		// setup molecule
		selection->construct(name);
	
		// setup atoms
		SpecAtom* a;
		for (int i = 0; i < latom.size(); i++)
		{
		    a = &latom[i];
		    if (a->name == "CA" || a->name == "CB" || a->name == "N" || a->name == "O" || a->name == "C")
    		    selection->append (a->pos, a->name, a->res, a->resno, a->respdb, a->occupancy, a->bfactor, a->sse);
		}

		// update
		selection->finalize();
		selection->lg();
	}

	// copy calphas
	void copyHeavy(SpecMolecule* selection)
	{
		// setup molecule
		selection->construct(name);
	
		// setup atoms
		SpecAtom* a;
		for (int i = 0; i < latom.size(); i++)
		{
		    a = &latom[i];
            if (a->name.size() > 0)
		    if (a->name[0] != 'H')
    		    selection->append (a->pos, a->name, a->res, a->resno, a->respdb, a->occupancy, a->bfactor, a->sse);
		}

		// update
		selection->finalize();
		selection->lg();
	}

	// number of elements
	inline int size()
	{
		return n;
	}

	// residues
	inline int sizeRes()
	{
		return (int) lcalpha.size();
	}

	// append atom
	inline void append (SpecAtom& a)
	{
		append (a.pos, a.name, a.res, a.resno, a.respdb, a.occupancy, a.bfactor, a.sse);
	}

	// append atom
	inline void append (Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, int sse = SSECOIL)
	{
		// increase counter
		n++;
		    	
		// resize
		latom.append (SpecAtom(pos, name, res, resno, respdb, occupancy, bfactor, sse));
	}

	// append hetatom
	inline void appendhet (Vec3 pos, string name, string res, int resno = 0, int respdb = LARGE)
	{
		lhetatom.append (SpecAtom(pos, name, res, resno, respdb));
	}

	// set
	inline void set (int i, Vec3 pos, string name, string res, int resno, int respdb = LARGE, double occupancy = 0.0, double bfactor = 0.0, int sse = SSECOIL)
	{
		latom[i].set(pos, name, res, resno, respdb, occupancy, bfactor, sse);
	}

	// center of gravity
	inline Vec3 cog()
	{
		// mean vector
		Vec3 tmp (0, 0, 0);
		int cx = 0;
		for (int i = 0; i < n; i++)
		{
   			tmp += latom[i].pos;
   	   		cx++;
        }
		tmp /= cx;
			
		// return 
		return tmp;
	}

	// center of gravity
	inline void setcenter()
	{
		// mean vector
		Vec3 tmp (0, 0, 0);
		for (int i = 0; i < n; i++)
			tmp += latom[i].pos;
		tmp /= n;
		
        for (int i = 0; i < n; i++)
            latom[i].pos += tmp;
	}
	
	// log
	void lg(bool details = false)
	{
		// sse seq
		string seq;        
		for (int i = 0; i < lsse.size(); i++)
		    seq += lsse[i].getType();

		// title
		Msg::write ("Molecule %s %s [%i, %i]", name.c_str(), seq.c_str(), lsse.size(), lcalpha.size());

		if (details)
		{
		    for (int i = 0; i < n; i++)
			Msg::write ("ATOM %i \t %f,%f,%f,%s,%s,%i,%c", i, latom[i].pos.x, latom[i].pos.y, latom[i].pos.z, latom[i].name.c_str(), latom[i].res.c_str(), latom[i].resno, latom[i].sse);	
		    
		    for (int i = 0; i < lsse.size(); i++)
			Msg::write ("SSE  %i \t %c \t %i \t %i", i, lsse[i].getType(), lsse[i].getStart(), lsse[i].getEnd());
		}
	}

	// finalize
	inline bool finalize()
	{   
		// verify
		if (n == 0)
		{
			Msg::write("No atoms or secondary structure was defined in molecule.");
			return false;
		}

		// verify consistency
		if (latom[0].resno != 0)
		{
			Msg::write("Residue numbers inconsistent.");
			return false;
		}
		for (int i = 0; i < latom.size() - 1; i++)
		{
			int diff = latom[i+1].resno - latom[i].resno;
			if (!(diff == 0 || diff == 1))
			{
				Msg::write("Residue numbers inconsistent.");
				return false;
			}
		}

		// update secondary structure information if missing
		if (lssebyres.size() == 0)
		{
			lssebyres.resize(latom[n - 1].resno + 1);
			lssebyres.fill(SSECOIL);
		}

		// verify consistency
		if (lssebyres.size() != latom[n - 1].resno + 1)
		{
			Msg::write("Residue numbers inconsistent with secondary information per residue (%i, %i).", lssebyres.size(), latom[n - 1].resno + 1);
			return false;
		}

		// setup calpha	
		finalizeAlpha();

		// setup secondary structure details
		finalizeSecondary();
			
		// return
		return true;
	}

	// check for c-alpha atom
	bool isCA(SpecAtom& a)
	{
		// check
		if(a.name == ATOMCA)
			return true;

		// return
		return false;
	}

	// rmsd upon forcetopology
	double rmsd(SpecMolecule* mol)
	{
		// verify
		if (lcalpha.size() != mol->lcalpha.size())
			Msg::error ("SpecMolecule::rmsd", "Number of atoms not equivalent");
		
		// calculate
		double dist = 0;
		int naligned = 0;
		for(int i = 0; i < lcalpha.size(); i++)
		{
			dist += pow (latom[lcalpha[i]].pos.dist(mol->latom[mol->lcalpha[i]].pos), 2);
            naligned++;			
        }

		// return
		return sqrt ( dist / naligned );
	}

private:	
	// initialize class.x
	void initialize(string& info, string& name)
	{
		// definitions
		ATOMCA = "CA";
		 
		// initialize
		this->info  	= info;
		this->name  	= name;
		this->n     	= 0;
	
		// clear all lists
		latom.clear();
		lhetatom.clear();		
		lsse.clear();
		lssebyres.clear();
		lcalpha.clear();
		lssecenter.clear();

		// clear all tables
		tresidue.clear();
	}

	// setup calpha list
	void finalizeAlpha()
	{
		// tresidue
		tresidue.resize(lssebyres.size());

		// lcalphas
		lcalpha.resize(lssebyres.size());
		lcalpha.fill(0);

		int cres = 0; 
		if (n > 0 && lssebyres.size() > 0)
		{
			// last residue
			int last = latom[0].resno;
	    
			// loop
			for (int i = 0; i < n; i++)
			{								
				// residue counter
				if (last != latom[i].resno)
				{
					last = latom[i].resno;
					cres++;
				}

				// is c-alpha atom
				if (isCA(latom[i]))
					lcalpha[cres] = i;
			
				// append residue table
				tresidue[cres].push_back(&latom[i]);
			
				// sse
				latom[i].sse = lssebyres[cres];
			}
		}
	}

	// setup secondary structure details
	void finalizeSecondary()
	{
		char last = 0;
		int nsse = 0;
		bool valid = false;

		// over all atoms
		for (int i = 0; i < n; i++)
		{            
			// append sses
			if (last != latom[i].sse)
			{                    
				if (valid)
					lsse[nsse-1].setEnd(latom[i-1].resno);

				if (latom[i].sse == SSEALPHA || latom[i].sse == SSEBETA)
				{
					lsse.push_back(SSE (latom[i].sse, latom[i].resno));
					nsse++;
				}

				// update
				last    = latom[i].sse;
				valid   = (last == SSEALPHA || last == SSEBETA);
			}

			// backup sse index
			if (latom[i].sse == SSEALPHA || latom[i].sse == SSEBETA)
				latom[i].sseno = nsse - 1;
		}

		// finalize
		if (valid)
			lsse[nsse-1].setEnd(latom[n-1].resno);

		//
		// CENTER CALPHA IN SSE
		//
		for (int i = 0; i < lsse.size(); i++)
			lssecenter.push_back (lcalpha[ lsse[i].getStart() + lsse[i].length() / 2 ]);
	}
};

// force set 
struct ForceSet 
{ 
	// parameter set 
	Vec3 force, amom, cog; 
	double angle, energy; 
	
	// load 
	ForceSet() 
	{ 
		// reset 
		reset(); 
	} 
	
	// reset parameters 
	void reset() 
	{ 
		// results 
		cog = Vec3(0.0, 0.0, 0.0); 
		force = Vec3(0.0, 0.0, 0.0); 
		amom = Vec3(0.0, 0.0, 0.0); 
		energy = 0.0; 
		angle = 1.0; 
	} 
}; 

// force topology 
struct ForceTopo 
{ 
	// index arrays 
	Vec < Vec < int > > a; 
	Vec < Vec < int > > b; 
	
	// resize 
	template < typename TTopo > 
	void construct (TTopo* ft) 
	{ 
		// clear 
		a.clear(); 
		b.clear(); 
		
		// resize 
		a.insert(a.begin(), ft->a.begin(), ft->a.end()); 
		b.insert(b.begin(), ft->b.begin(), ft->b.end()); 
	} 
	
	// clear 
	void clear() 
	{ 
		a.clear(); 
		b.clear(); 
	} 
	
	// resize 
	void resize(int n) 
	{ 
		// clear 
		a.clear(); 
		b.clear(); 
		
		// resize 
		a.resize(n); 
		b.resize(n); 
	} 
	
	// reset 
	void reset(int n) 
	{ 
		resize(n); 
	} 
	
	// number of groups 
	int nGroups() 
	{ 
		assert (a.size() == b.size()); 
		return a.size(); 
	} 
	
	// count interaction pairs 
	int npairs() 
	{ 
		int n = 0; 
		for (int k = 0; k < a.size(); k++) 
		n += a[k].size() * b[k].size(); 
		return n; 
	} 
}; 


// rotation 
struct Rotation 
{ 
	public: 
	double q[4]; 
	
	// axis and angle to quaternion 
	Rotation (Vec3& axis, double angle) 
	{ 
		// prepare 
		axis.normalize(); 
		axis *= sin(angle / 2.0); 
		
		// set 
		q[0] = axis[0]; 
		q[1] = axis[1]; 
		q[2] = axis[2]; 
		q[3] = cos(angle / 2.0); 
	} 
	
	// axis and angle to quaternion 
	Rotation (Vec4& v) 
	{ 
		q[0] = v.w; 
		q[1] = v.x; 
		q[2] = v.y; 
		q[3] = v.z; 
	} 
}; 

// matrix 
class Mat4 
{ 
	
	// matrix 
	double m[4][4]; 
	public: 
	void identity() 
	{ 
		m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0; m[0][3] = 0.0; 
		m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0; m[1][3] = 0.0; 
		m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0; m[2][3] = 0.0; 
		m[3][0] = 0.0; m[3][1] = 0.0; m[3][2] = 0.0; m[3][3] = 1.0; 
	} 
	
	void setRotate (double* q) 
	{ 
		m[0][0] = (double)(1 - 2.0 * (q[1] * q[1] + q[2] * q[2])); 
		m[0][1] = (double)( 2.0 * (q[0] * q[1] + q[2] * q[3])); 
		m[0][2] = (double)( 2.0 * (q[2] * q[0] - q[1] * q[3])); 
		m[0][3] = 0.0; 
		m[1][0] = (double)( 2.0 * (q[0] * q[1] - q[2] * q[3])); 
		m[1][1] = (double)(1 - 2.0 * (q[2] * q[2] + q[0] * q[0])); 
		m[1][2] = (double)( 2.0 * (q[1] * q[2] + q[0] * q[3])); 
		m[1][3] = 0.0; 
		m[2][0] = (double)( 2.0 * (q[2] * q[0] + q[1] * q[3])); 
		m[2][1] = (double)( 2.0 * (q[1] * q[2] - q[0] * q[3])); 
		m[2][2] = (double)(1 - 2.0 * (q[1] * q[1] + q[0] * q[0])); 
		m[2][3] = 0.0; 
		m[3][0] = 0.0; 
		m[3][1] = 0.0; 
		m[3][2] = 0.0; 
		m[3][3] = 1.0; 
	} 
	
	void setTranslate (Vec3& t) 
	{ 
		m[0][0] = 1.0; 
		m[0][1] = 0.0; 
		m[0][2] = 0.0; 
		m[0][3] = 0.0; 
		m[1][0] = 0.0; 
		m[1][1] = 1.0; 
		m[1][2] = 0.0; 
		m[1][3] = 0.0; 
		m[2][0] = 0.0; 
		m[2][1] = 0.0; 
		m[2][2] = 1.0; 
		m[2][3] = 0.0; 
		m[3][0] = t[0]; 
		m[3][1] = t[1]; 
		m[3][2] = t[2]; 
		m[3][3] = 1.0; 
	} 
	
	void setTransform(Vec3 &t, double* q, Vec3 &c) 
	{ 
		Mat4 matrix; 
		
		// reset 
		identity(); 
		
		// encode translation 
		if (t != 0) 
		{ 
			matrix.setTranslate(t); 
			set(matrix * (*this)); 
		} 
		
		// shift to center 
		if (c != 0) 
		{ 
			matrix.setTranslate(c); 
			set(matrix * (*this)); 
		} 
		
		// rotate 
		matrix.setRotate(q); 
		set(matrix * (*this)); 
		
		// shift back 
		if (c != 0) 
		{ 
			c.neg(); 
			matrix.setTranslate(c); 
			set(matrix * (*this)); 
		} 
	} 
	
	friend Mat4 operator*(const Mat4& m1, const Mat4& m2) 
	{ 
		Mat4 out; 
		for (int i = 0; i < 4; i++) 
		{ 
			for (int j = 0; j < 4; j++) 
			{ 
				out.m[i][j] = 0.0; 
				for (int k = 0; k < 4; k++) 
				out.m[i][j] += m1.m[i][k] * m2.m[k][j]; 
			} 
		} 
		return out; 
	} 
	
	void reset() 
	{ 
		for (int i = 0; i < 4; i++) 
		for (int j = 0; j < 4; j++) 
		m[i][j] = 0.0; 
	} 
	
	void set(const Mat4& matrix ) 
	{ 
		for (int i = 0; i < 4; i++) 
		for (int j = 0; j < 4; j++) 
		m[i][j] = matrix.m[i][j]; 
	} 
	
	void multVecMatrix(Vec3 &src, Vec3 &dst) 
	{ 
		double f = 1.0 / (src[0]*m[0][3]+src[1]*m[1][3]+src[2]*m[2][3]+m[3][3]); 
		dst[0] = f*(src[0]*m[0][0]+src[1]*m[1][0]+src[2]*m[2][0]+m[3][0]); 
		dst[1] = f*(src[0]*m[0][1]+src[1]*m[1][1]+src[2]*m[2][1]+m[3][1]); 
		dst[2] = f*(src[0]*m[0][2]+src[1]*m[1][2]+src[2]*m[2][2]+m[3][2]); 
	} 
}; 



// transformation 
struct Trans 
{ 
	/** 
	apply transformation 
	**/ 
	
	// apply on whole molecule 
	static void apply (Vec<Vec3>& lpos, Vec3& force, Vec3& amom, double& angle, Vec3& cog) 
	{ 
		// parameter 
		Vec3 tmp; 
		Mat4 trans; 
		
		// rotation 
		Rotation rot (amom, -angle); 
		
		// generate transform matrix (rigid transformation) 
		trans.setTransform(force, rot.q, cog); 
		
		// compute new positions 
		for(int i = 0; i < lpos.size(); i++) 
		{ 
			// transformation 
			trans.multVecMatrix(lpos[i], tmp); 
			lpos[i] = tmp; 
		} 
	} 
	
	// apply on whole molecule 
	static void apply (SpecMolecule* mol, Vec<Vec3>& v, ForceSet& fs) 
	{ 
		// parameter 
		Vec3 tmp; 
		Mat4 trans; 
		Rotation rot (fs.amom, -fs.angle); 
		
		// generate transform matrix (rigid transformation) 
		trans.setTransform(fs.force, rot.q, fs.cog); 
		
		// compute new positions 
		for(int i = 0; i < mol->size(); i++) 
		{ 
			trans.multVecMatrix(mol->latom[i].pos, tmp); 
			mol->latom[i].pos = tmp; 
		} 
		
		// compute new positions 
		for(int i = 0; i < v.size(); i++) 
		{ 
			trans.multVecMatrix(v[i], tmp); 
			v[i] = tmp; 
		} 
	} 
	
	// apply on whole molecule 
	static void apply (SpecMolecule* mol, ForceSet& fs) 
	{ 
		// parameter 
		Vec3 tmp; 
		Mat4 trans; 
		Rotation rot (fs.amom, -fs.angle); 
		
		// generate transform matrix (rigid transformation) 
		trans.setTransform(fs.force, rot.q, fs.cog); 
		
		// compute new positions 
		for(int i = 0; i < mol->size(); i++) 
		{ 
			trans.multVecMatrix(mol->latom[i].pos, tmp); 
			mol->latom[i].pos = tmp; 
		} 
	} 
	
	// apply on whole molecule and a set of coordinates 
	static void apply (SpecMolecule* mol, Vec <Vec3>& lcoord, Vec <Vec3>& lnormals, ForceSet& fs) 
	{ 
		// parameter 
		Vec3 tmp; 
		Mat4 trans; 
		Rotation rot (fs.amom, -fs.angle); 
		
		// generate transform matrix (rigid transformation) 
		trans.setTransform(fs.force, rot.q, fs.cog); 
		
		// compute new positions 
		for(int i = 0; i < mol->size(); i++) 
		{ 
			trans.multVecMatrix(mol->latom[i].pos, tmp); 
			mol->latom[i].pos = tmp; 
		} 
		
		// compute new positions 
		for(int i = 0; i < lcoord.size(); i++) 
		{ 
			trans.multVecMatrix(lcoord[i], tmp); 
			lcoord[i] = tmp; 
		} 
		
		// compute new positions 
		for(int i = 0; i < lnormals.size(); i++) 
		{ 
			trans.multVecMatrix(lnormals[i], tmp); 
			lnormals[i] = tmp; 
		} 
	} 
	
	// apply on whole molecule and a set of coordinates 
	static void apply (SpecMolecule* mol, Vec <Vec3>& lcoord, Vec4& quat, Vec3 force = Vec3(0), Vec3 cog = Vec3(0)) 
	{ 
		// parameter 
		Vec3 tmp; 
		Mat4 trans; 
		
		// rotation 
		Rotation rot (quat); 
		
		// generate transform matrix (rigid transformation) 
		trans.setTransform(force, rot.q, cog); 
		
		// compute new positions 
		for(int i = 0; i < mol->size(); i++) 
		{ 
			trans.multVecMatrix(mol->latom[i].pos, tmp); 
			mol->latom[i].pos = tmp; 
		} 
		
		// compute new positions 
		for(int i = 0; i < lcoord.size(); i++) 
		{ 
			trans.multVecMatrix(lcoord[i], tmp); 
			lcoord[i] = tmp; 
		} 
	} 
	
	// apply on whole molecule 
	static void apply (SpecMolecule* target, SpecMolecule* mol, ForceSet& fs) 
	{ 
		// parameter 
		Vec3 tmp; 
		Mat4 trans; 
		Rotation rot (fs.amom, -fs.angle); 
		
		// generate transform matrix (rigid transformation) 
		trans.setTransform(fs.force, rot.q, fs.cog); 
		
		// compute new positions 
		for(int i = 0; i < target->size(); i++) 
		{ 
			trans.multVecMatrix(target->latom[i].pos, tmp); 
			target->latom[i].pos = tmp; 
		} 
		
		// compute new positions 
		for(int i = 0; i < mol->size(); i++) 
		{ 
			trans.multVecMatrix(mol->latom[i].pos, tmp); 
			mol->latom[i].pos = tmp; 
		} 
	} 
	
	// apply translation on whole molecule 
	template < typename TVector > 
	static void apply (SpecMolecule* mol, TVector& shift) 
	{ 
		// compute new positions 
		for(int i = 0; i < mol->size(); i++) 
		{ 
			mol->latom[i].pos.x += shift.x; 
			mol->latom[i].pos.y += shift.y; 
			mol->latom[i].pos.z += shift.z; 
		} 
	} 
	
	// transformation matrix 
	inline static void apply (Vec < Vec3>& pos, Vec < Vec3 >& m, Vec3& t) 
	{ 
		// apply translation on every atom 
		double x, y, z; 
		for (int i = 0; i < pos.size(); i++) 
		{ 
			// calculate 
			x = pos[i].x * m[0][0] + pos[i].y * m[0][1] + pos[i].z * m[0][2] + t[0]; 
			y = pos[i].x * m[1][0] + pos[i].y * m[1][1] + pos[i].z * m[1][2] + t[1]; 
			z = pos[i].x * m[2][0] + pos[i].y * m[2][1] + pos[i].z * m[2][2] + t[2]; 
			
			// update 
			pos[i].x = x; 
			pos[i].y = y; 
			pos[i].z = z; 
		} 
	} 
	
	// transformation matrix 
	inline static void apply (SpecMolecule* mol, Vec < Vec3 >& m, Vec3& t) 
	{ 
		// 
		// apply translation on every atom 
		double x, y, z; 
		for (int i = 0; i < mol->size(); i++) 
		{ 
			// calculate 
			x = mol->latom[i].pos.x * m[0][0] + mol->latom[i].pos.y * m[0][1] + mol->latom[i].pos.z * m[0][2] + t[0]; 
			y = mol->latom[i].pos.x * m[1][0] + mol->latom[i].pos.y * m[1][1] + mol->latom[i].pos.z * m[1][2] + t[1]; 
			z = mol->latom[i].pos.x * m[2][0] + mol->latom[i].pos.y * m[2][1] + mol->latom[i].pos.z * m[2][2] + t[2]; 
			
			// update 
			mol->latom[i].pos.x = x; 
			mol->latom[i].pos.y = y; 
			mol->latom[i].pos.z = z; 
		} 
		
		// all other atoms 
		for (int i = 0; i < mol->lhetatom.size(); i++) 
		{ 
			// calculate 
			x = mol->lhetatom[i].pos.x * m[0][0] + mol->lhetatom[i].pos.y * m[0][1] + mol->lhetatom[i].pos.z * m[0][2] + t[0]; 
			y = mol->lhetatom[i].pos.x * m[1][0] + mol->lhetatom[i].pos.y * m[1][1] + mol->lhetatom[i].pos.z * m[1][2] + t[1]; 
			z = mol->lhetatom[i].pos.x * m[2][0] + mol->lhetatom[i].pos.y * m[2][1] + mol->lhetatom[i].pos.z * m[2][2] + t[2]; 
			
			// update 
			mol->lhetatom[i].pos.x = x; 
			mol->lhetatom[i].pos.y = y; 
			mol->lhetatom[i].pos.z = z; 
		} 
	} 
	
	/** 
	align 
	**/ 
	
	// kabsch alignment 
	inline static void align (Vec < Vec3 >& a, Vec < Vec3 >& b, SpecMolecule* mol) 
	{ 
		// parameters 
		Vec < Vec3 > r; 
		Vec3 t; 
		
		// apply 
		align (a, b, r, t); 
		apply(mol, r, t); 
	} 
	
	// kabsch alignment 
	inline static void align (Vec < Vec3 >& a, Vec < Vec3 >& b) 
	{ 
		// parameters 
		Vec < Vec3 > r; 
		Vec3 t; 
		
		// apply 
		align (a, b, r, t); 
		apply(a, r, t); 
	} 
	
	// kabsch alignment 
	inline static void align (Vec < Vec3 >& a, Vec < Vec3 >& b, Vec < Vec3 >& rot, Vec3& trans) 
	{ 
		// amount 
		int n = a.size(); 
		
		// results 
		Kabsch::MATRIX M = { {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,0}} }; 
		
		// reset 
		if (n == b.size() && n > 3) 
		{ 
			// setup matrices a and b 
			double ref[n*3]; 
			double mov[n*3]; 
			for (int i = 0; i < n; i++) 
			for (int j = 0; j < 3; j++) 
			{ 
				mov[i*3 + j] = a[i][j]; 
				ref[i*3 + j] = b[i][j]; 
			} 
			
			// get matrix 
			Kabsch kbs; 
			kbs.rmsd(ref, mov, n, (double*) &M.m); 
		} 
		
		// copy transformation 
		rot.resize(3); 
		for (int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++) 
		rot[i][j] = M.m[j][i]; 
		for (int i = 0; i < 3; i++) 
		trans[i] = M.m[3][i]; 
		
		// fix chirality 
		double det = rot[0][0]*rot[1][1]*rot[2][2] + rot[0][1]*rot[1][2]*rot[2][0] + rot[0][2]*rot[1][0]*rot[2][1] 
		- rot[0][2]*rot[1][1]*rot[2][0] - rot[0][1]*rot[1][0]*rot[2][2] - rot[0][0]*rot[1][2]*rot[2][1]; 
		if (det < 0) 
		rot[0] *= -1.0; 
	} 
	
	/** 
	complex rmsd calculator 
	**/ 
	
	// rmsd upon forcetopology 
	static double rmsd (SpecMolecule* target, SpecMolecule* ref, SpecMolecule* mol) 
	{ 
		// verify 
		if (ref->lcalpha.size() != mol->lcalpha.size()) 
		Msg::error ("Trans<T>::rmsd", "Amount of atoms not equivalent"); 
		
		// sizes 
		int nt = target->lcalpha.size(); 
		int nr = ref->lcalpha.size(); 
		int n = nt + nr; 
		
		// vector list 
		Vec < Vec3 > va, vb; 
		va.resize(n); 
		vb.resize(n); 
		
		// generate decoy complex 
		for(int i = 0; i < nt; i++) 
		va[i] = vb[i] = target->latom[target->lcalpha[i]].pos; 
		
		for(int i = 0; i < nr; i++) 
		{ 
			va[nt + i] = mol->latom[mol->lcalpha[i]].pos; 
			vb[nt + i] = ref->latom[ref->lcalpha[i]].pos; 
		} 
		
		// parameters 
		Vec < Vec3 > r; 
		Vec3 t; 
		
		// apply 
		align (va, vb, r, t); 
		apply(va, r, t); 
		
		// calculate 
		double dist = 0; 
		for(int i = 0; i < n; i++) 
		dist += pow (va[i].dist(vb[i]), 2); 
		
		// return 
		return sqrt ( dist / n ); 
	} 
}; 

// job 
struct TransMatrix 
{ 
	// values 
	Vec < Vec3 > r, from, to; 
	Vec3 t, center; 
	
	// static size 
	static const int nalpha = 5; 
	
	// construct 
	TransMatrix() 
	{ 
		t = 0.0; 
		r.resize(3); 
		r.fill(0.0); 
	} 
	
	// set reference 
	void construct (SpecMolecule* mol) 
	{ 
		// read from 
		from.resize(nalpha); 
		for (int i = 0; i < nalpha; i++) 
		from[i] = mol->latom[mol->lcalpha[i]].pos; 
	} 
	
	// copy constructor 
	void construct (TransMatrix& matrix) 
	{ 
		from.assign(matrix.from.begin(), matrix.from.end()); 
		to.assign(matrix.to.begin(), matrix.to.end()); 
		r.assign(matrix.r.begin(), matrix.r.end()); 
		t = matrix.t; 
	} 
	
	// set orientations 
	void prepare (SpecMolecule* mol) 
	{ 
		// read to 
		to.resize(nalpha); 
		for (int i = 0; i < nalpha; i++) 
		to[i] = mol->latom[mol->lcalpha[i]].pos; 
		
		// get matrix 
		Trans::align(from, to, r, t); 
		
		// center 
		center = mol->cog(); 
	} 
	
	// set orientations 
	void prepare (Vec <Vec3 > pos) 
	{ 
		// read to 
		to.resize(nalpha); 
		for (int i = 0; i < nalpha; i++) 
		to[i] = pos[i]; 
		
		// get matrix 
		Trans::align(from, to, r, t); 
	} 
	
	// set orientations 
	void apply (SpecMolecule* mol, Vec <Vec3 >& pos) 
	{ 
		// read from 
		from.resize(nalpha); 
		for (int i = 0; i < nalpha; i++) 
		from[i] = mol->latom[mol->lcalpha[i]].pos; 
		
		// from 
		to.resize(nalpha); 
		for (int i = 0; i < nalpha; i++) 
		to[i] = pos[i]; 
		
		// get matrix 
		Trans::align(from, to, r, t); 
		Trans::apply(mol, r, t); 
	} 
	
	// apply 
	void apply (SpecMolecule* mol) 
	{ 
		Trans::apply(mol, r, t); 
	} 
	
	// apply 
	void apply (Vec <Vec3>& pos) 
	{ 
		Trans::apply(pos, r, t); 
	} 
	
	// print 
	void print(double epsilon = 0.000001) 
	{ 
		// clean up zeros 
		for (int i = 0; i < 3; i++) 
		{ 
			for (int j = 0; j < 3; j++) 
			if (fabs(r[i][j]) < epsilon) 
			r[i][j] = 0.0; 
			if (fabs(t[i]) < epsilon) 
			t[i] = 0.0; 
		} 
		
		// write 
		printf (" i          t(i)         u(i,1)         u(i,2)         u(i,3)\n"); 
		printf (" 1    %14.10f %14.10f %14.10f %14.10f\n", t[0], r[0][0], r[0][1], r[0][2]); 
		printf (" 2    %14.10f %14.10f %14.10f %14.10f\n", t[1], r[1][0], r[1][1], r[1][2]); 
		printf (" 3    %14.10f %14.10f %14.10f %14.10f\n", t[2], r[2][0], r[2][1], r[2][2]); 
	} 
	
	// reverse 
	void reverse() 
	{ 
		// swap 0 
		double s = 0; 
		s = r[1][0]; 
		r[1][0] = r[0][1]; 
		r[0][1] = s; 
		
		// swap 0 
		s = r[2][0]; 
		r[2][0] = r[0][2]; 
		r[0][2] = s; 
		
		// swap 0 
		s = r[1][2]; 
		r[1][2] = r[2][1]; 
		r[2][1] = s; 
		
		// Invert the translation part by negating it; and then rotate it by the new rotation part 
		Vec3 v = Vec3(-t.x, -t.y, -t.z); 
		t.x = v.x * r[0][0] + v.y * r[0][1] + v.z * r[0][2]; 
		t.y = v.x * r[1][0] + v.y * r[1][1] + v.z * r[1][2]; 
		t.z = v.x * r[2][0] + v.y * r[2][1] + v.z * r[2][2]; 
	} 
	
	// get string 
	string toString() 
	{ 
		return Convert::toString(r[0].x) + " " + 
		Convert::toString(r[0].y) + " " + 
		Convert::toString(r[0].z) + " " + 
		Convert::toString(r[1].x) + " " + 
		Convert::toString(r[1].y) + " " + 
		Convert::toString(r[1].z) + " " + 
		Convert::toString(r[2].x) + " " + 
		Convert::toString(r[2].y) + " " + 
		Convert::toString(r[2].z) + " " + 
		Convert::toString(t.x) + " " + 
		Convert::toString(t.y) + " " + 
		Convert::toString(t.z); 
	} 
	
	// prepare matrix 
	static double load (SpecMolecule* mol, string input, double& rmsd) 
	{ 
		// split data 
		Vec < string > data; 
		Lib::split (data, input.c_str()); 
		
		// verify 
		if (data.size() <= 12) 
		return -1.0; 
		
		// load data 
		Vec < Vec3 > r; 
		r.resize(3); 
		r.fill(0.0); 
		Vec3 t = 0.0; 
		
		// read 
		r[0].x = Convert::toDbl(data[0]); 
		r[0].y = Convert::toDbl(data[1]); 
		r[0].z = Convert::toDbl(data[2]); 
		r[1].x = Convert::toDbl(data[3]); 
		r[1].y = Convert::toDbl(data[4]); 
		r[1].z = Convert::toDbl(data[5]); 
		r[2].x = Convert::toDbl(data[6]); 
		r[2].y = Convert::toDbl(data[7]); 
		r[2].z = Convert::toDbl(data[8]); 
		t.x = Convert::toDbl(data[9]); 
		t.y = Convert::toDbl(data[10]); 
		t.z = Convert::toDbl(data[11]); 
		
		// apply 
		Trans::apply(mol, r, t); 
		
		// read last column rmsd 
		if (data.size() <= 13) 
		rmsd = (double) LARGE; 
		else 
		rmsd = Convert::toDbl(data[13]); 
		
		// return 
		return Convert::toDbl(data[12]); 
	} 
}; 

// sse 
class Secondary 
{ 
	
	public: 
	// variables 
	int j1, j2, j3, j4, j5; 
	double dis13, dis14, dis15, dis24, dis25, dis35; 
	
	// print 
	void load (SpecMolecule* mol, Vec <Vec3>& lcalpha) 
	{ 
		// prepare sse by residue list 
		int n = lcalpha.size(); 
		if (n == 0) 
		{ 
			Msg::write("No Residues found in PDB."); 
			return; 
		} 
		
		// prepare 
		mol->lssebyres.resize(n); 
		mol->lssebyres.fill(SpecMolecule::SSECOIL); 
		
		// 1->coil, 2->helix, 3->turn, 4->strand 
		for (int i = 1; i < n; i++) 
		{ 
			j1=i-2; 
			j2=i-1; 
			j3=i; 
			j4=i+1; 
			j5=i+2; 
			if(j1 > 0 && j5 < n) 
			{ 
				dis13 = lcalpha[j1].dist(lcalpha[j3]); 
				dis14 = lcalpha[j1].dist(lcalpha[j4]); 
				dis15 = lcalpha[j1].dist(lcalpha[j5]); 
				dis24 = lcalpha[j2].dist(lcalpha[j4]); 
				dis25 = lcalpha[j2].dist(lcalpha[j5]); 
				dis35 = lcalpha[j3].dist(lcalpha[j5]); 
				mol->lssebyres[i] = make_sec(dis13,dis14,dis15,dis24,dis25,dis35); 
			} 
		} 
		
		// smooth 
		smooth(mol); 
	} 
	
	// print 
	void loadgeneric (SpecMolecule* mol, Vec <Vec3>& lcalpha, int splitfactor = 30) 
	{ 
		// prepare sse by residue list 
		int n = lcalpha.size(); 
		if (n == 0) 
		{ 
			Msg::write("No Residues found in PDB."); 
			return; 
		} 
		
		int split = max (int (n / double (splitfactor)), 5); 
		
		// prepare 
		mol->lssebyres.resize(n); 
		mol->lssebyres.fill(SpecMolecule::SSEBETA); 
		
		// 1->coil, 2->helix, 3->turn, 4->strand 
		for (int i = 0; i < n; i++) 
		{ 
			if (i % split == 0) 
			mol->lssebyres[i] = SpecMolecule::SSECOIL; 
		} 
	} 
	
	// print 
	int make_sec (double dis13, double dis14, double dis15, double dis24, double dis25, double dis35) 
	{ 
		// helix 
		double delta = 2.5; //2.1 ->2.5 
		if(fabs(dis15-6.37) < delta) 
		if(fabs(dis14-5.18) < delta) 
		if(fabs(dis25-5.18) < delta) 
		if(fabs(dis13-5.45) < delta) 
		if(fabs(dis24-5.45) < delta) 
		if(fabs(dis35-5.45) < delta) 
		return SpecMolecule::SSEALPHA; 
		
		// strand 
		delta = 2.0; //1.42 ->2.0 
		if(fabs(dis15-13.0) < delta) 
		if(fabs(dis14-10.4) < delta) 
		if(fabs(dis25-10.4) < delta) 
		if(fabs(dis13-6.1) < delta) 
		if(fabs(dis24-6.1) < delta) 
		if(fabs(dis35-6.1) < delta) 
		return SpecMolecule::SSEBETA; 
		
		// turn 
		if (dis15 < 8.0) 
		return SpecMolecule::SSETURN; 
		
		// coil 
		return SpecMolecule::SSECOIL; 
	} 
	
	// smooth 
	void smooth (SpecMolecule* mol) 
	{ 
		int j; 
		
		// x%x >> xxx 
		for (int i = 0; i < mol->lssebyres.size() - 2; i++) 
		if(mol->lssebyres[i] != SpecMolecule::SSECOIL) 
		{ 
			j = mol->lssebyres[i]; 
			if(mol->lssebyres[i+2] == j) 
			mol->lssebyres[i+1] = j; 
		} 
		
		// -x- >> ----- 
		for (int i = 0; i < mol->lssebyres.size() - 2; i++) 
		if(mol->lssebyres[i] == SpecMolecule::SSECOIL) 
		if(mol->lssebyres[i+1] != SpecMolecule::SSECOIL) 
		if(mol->lssebyres[i+2] == SpecMolecule::SSECOIL) 
		mol->lssebyres[i+1] = SpecMolecule::SSECOIL; 
	} 
	
	// according to tmalign 
	static int getnumeric (char ssecode) 
	{ 
		if (ssecode == SpecMolecule::SSEALPHA) 
		return 2; 
		else if (ssecode == SpecMolecule::SSEBETA) 
		return 4; 
		else if (ssecode == SpecMolecule::SSETURN) 
		return 3; 
		else 
		return 1; 
	} 
}; 

class Storage
{

	// file
	File f;	

public:
	// molecule
	string path, pathalternative;

	// construct
	void construct (string path = "", string pathalternative = "")
	{
		// initialize
		this->path              = path;
		this->pathalternative   = pathalternative;
	}   
	
	// read
	bool read (SpecMolecule* mol, string pdbname, string name = "Molecule")
	{        
		// indices
		Vec < Vec3 > lcalpha;
	    
		// setup molecule
		mol->construct(name);

		// load pdb with/without hydrogens
		pdb(mol, pdbname, lcalpha);

		// load structure
		if (mol->size() > 0)
		{
			Secondary sec;
			sec.load(mol, lcalpha);
		}

		// finalize
		if (mol->finalize())
		{
			if (mol->lcalpha.size() < 5)
			{
				Msg::write("Molecule has less than five residues.");
				return false;
			} else
				return true;   
		} else
			return false; 
	}

	// pdb reader
	void pdb(SpecMolecule* mol, string name, Vec <Vec3>& lrespos, bool verbose = true)
	{   
		// read file     
		load (name);

        // key
		string key = "";
		
		// loop
        int nres = 0;
		int resno = 0;
		string chain = "";
		string resname = "";
		
		// residue information
        Vec3 respos;
		
		// validate residue
		SpecMolecule residue;
		residue.construct();

		// read file
		while (f.good())
		{
			// read
			f.readLine();

			// get pdb key
			key = Lib::strclean(f.get(0, 6));

			// end model
			if (key == "END" || key == "ENDMDL" || key == "TER")
				break;

			// atom
			if (key == "ATOM" || (key == "HETATM" && f.get(21, 1) == chain))
			{
				// generate index for sse information				
				if (resname == "")
				{
					// backup residue details
                    resname = f.get(17, 3);
					chain   = f.get(21, 1);
	       			resno   = Convert::toInt(f.get(22, 4));
					
                    // reset
					respos  = 0.0;
				}

				// verify new residue
				if (resname != f.get(17, 3) || chain != f.get(21, 1) || resno != Convert::toInt(f.get(22, 4)))
				{
					// verify
					if (respos == 0.0)
					{
                        if (verbose)
						  Msg::write("Residue %i incomplete in %s.", nres, name.c_str());
				    } else {
    					// copy residue into molecule
                        for (int i = 0; i < residue.size(); i++)
                            mol->append(residue.latom[i]);

                        // add calpha
                        lrespos.push_back(respos);
                        
                        // count
                        nres++;
                    }
                
					// reset residue
					residue.construct();
	
					// backup residue details
                    resname = f.get(17, 3);
					chain   = f.get(21, 1);
	       			resno   = Convert::toInt(f.get(22, 4));
					
                    // reset
					respos     = 0.0;					
				}

				// get atom name and coordinate
				string atomname = f.get(12, 4);
				Vec3 atompos = Vec3( Convert::toDbl(f.get(30, 8)) , Convert::toDbl(f.get(38, 8)), Convert::toDbl(f.get(46, 8)) );

				// check calpha coordinate
				if (Lib::strclean(atomname) == "CA")
					respos = atompos;

				// append			
				residue.append(atompos, atomname, resname, nres, resno, Convert::toDbl(f.get(54, 6)), Convert::toDbl(f.get(60, 6)));
			} else if (key == "HETATM")
				mol->appendhet(Vec3( Convert::toDbl(f.get(30, 8)) , Convert::toDbl(f.get(38, 8)), Convert::toDbl(f.get(46, 8)) ), f.get(12, 4), f.get(17, 3));
		}
        		
		// append last residue
        if (respos != 0.0)
        {
            // add atoms
            for (int i = 0; i < residue.size(); i++)
                mol->append(residue.latom[i]);

            // add calpha
            lrespos.push_back(respos);
        }

		// return
		f.close();
	}

	// save
	void save (Vec < SpecMolecule >& molb, string name)
	{
		FILE *fp;
		open (fp, name);
	
		for (int i = 0; i < molb.size(); i++)
			print (fp, &molb[i], Convert::toString(i));
	
		// finalize
		close (fp);
	}
	
	// save
	void save (Vec < SpecMolecule* >& molb, string name)
	{
		FILE *fp;
		open (fp, name);
	
		for (int i = 0; i < molb.size(); i++)
			print (fp, molb[i], Convert::toString(i));
	
		// finalize
		close (fp);
	}
	
	// save
	void save (SpecMolecule* mola, Vec < SpecMolecule* >& molb, string name)
	{
		FILE *fp;
		open(fp, name);    
	
		// print
		print (fp, mola, "A", "receptor");
		for (int i = 0; i < molb.size(); i++)
		{
			if (i == 0)
				print (fp, molb[i], "B", "ligand reference");
			else
				printMinimal (fp, molb[i], "B", "ligand " + Convert::toString(i-1));                
		}
	
		// finalize
		close (fp);
	}
		   
	void save (SpecMolecule* mola, SpecMolecule* molb, string name)
	{
		FILE *fp;
		open(fp, name);                
		print (fp, mola);
		print (fp, molb, "B");
		close (fp);
	}

	void save (SpecMolecule* mol, string name, bool minimal = false)
	{
		FILE *fp;
		open (fp, name);
		
		// print
		if (minimal)
    		printMinimal (fp, mol);		
		else
    		print (fp, mol);

		// finalize
		close (fp);
	}

private:
	// open file
	void open(FILE*& fp, string name)
	{
		string nm (name);
	    if ((fp = fopen(nm.c_str(), "wb")) == NULL)
			Msg::error("open()", "Can not open file " + string(nm.c_str()));
	}
	
	void load(string name)
	{
		// setup path		
		if (!f.open(path + name))
			if (!f.open(pathalternative + name))
				Msg::write ("Unable to open file %s.", string (path + name).c_str());
	}
            
	// close file
	void close(FILE* fp)
	{
		fflush(fp);
		fclose(fp);
	}
	
	// write molecule to file
	void print (FILE *fp, SpecMolecule* mol, string chain = "A", bool details = false)
	{	    
        // info
        if (mol->info != "")
            fprintf (fp, string (mol->info + "\n").c_str());            
            
		// print more details
		if (details)
		{
			// are there sses
			if (mol->lsse.size() > 0)
			{
				fprintf (fp, "REM Secondary structure information\n");
				fprintf (fp, "REM No. Type \t From \t To\n");
				for (int i = 0; i < mol->lsse.size(); i++)
					fprintf (fp, "REM %i \t %c \t %i \t %i\n", i, mol->lsse[i].getType(), mol->lsse[i].getStart(), mol->lsse[i].getEnd());
			}
			
			// calphas
			fprintf (fp, "REM C-Alpha table\n");
			for (int i = 0; i < mol->lcalpha.size(); i++)
				fprintf (fp, "REM %i \t %i \n", i, mol->lcalpha[i]);
	
			// sse per residues
			fprintf (fp, "REM SSE-Residue table\n");
			for (int i = 0; i < mol->lssebyres.size(); i++)
				fprintf (fp, "REM %i \t %c \n", i, mol->lssebyres[i]);
	
			// sse centers
			fprintf (fp, "REM SSE center table\n");
			for (int i = 0; i < mol->lssecenter.size(); i++)
				fprintf (fp, "REM %i \t %i \n", i, mol->lssecenter[i]);
		}
		
		// print atoms
		for (int i = 0; i < mol->size(); i++)
		{
            fprintf (fp, "ATOM  %5d  %-4s%-3s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12c    %12c\n",
                        i,
				        mol->latom[i].name.c_str(),
						mol->latom[i].res.c_str(),
						chain.c_str(),
						min (9999, mol->latom[i].resno),
						mol->latom[i].pos.x, mol->latom[i].pos.y, mol->latom[i].pos.z,
						min(99.0, mol->latom[i].occupancy), max (-99.0, min(99.0, mol->latom[i].bfactor)), mol->latom[i].name[0], mol->latom[i].sse);
        }

		// print hetatoms
		for (int i = 0; i < mol->lhetatom.size(); i++)
		{
			fprintf (fp, "HETATM%5d  %-4s%-4s%1s%4d    %8.3f%8.3f%8.3f  1.00  1.00%12c    %12c\n",
						i,
						mol->lhetatom[i].name.c_str(),
						mol->lhetatom[i].res.c_str(),
						"H",
						min (9999, mol->lhetatom[i].resno),
						mol->lhetatom[i].pos.x, mol->lhetatom[i].pos.y, mol->lhetatom[i].pos.z,
						mol->lhetatom[i].name[0], mol->lhetatom[i].sse);
        }
	}
	
	// write molecule in minimal mode
	void printMinimal (FILE *fp, SpecMolecule* mol, string chain = "A", string tag = "")
	{
		// log
		if (tag != "")
		{
			tag = "REM " + tag + "\n";			
			fprintf (fp, tag.c_str());
		}
			
		// loop
		for (int index = 0; index < mol->lcalpha.size(); index++)
		{
			int i = mol->lcalpha[index];
			fprintf(fp,"ATOM  %5d  %-4s%-4s%1s%4d    %8.3f%8.3f%8.3f\n",
						i, 	mol->latom[i].name.c_str(),
						mol->latom[i].res.c_str(),
						chain.c_str(),
						min (9999, mol->latom[i].resno),
						mol->latom[i].pos.x, mol->latom[i].pos.y, mol->latom[i].pos.z);
		}		
	}
};

/** 
Minimization Interface 
*/ 
template < typename TInput > 
class IStep 
{ 
	public : 
	// weight or step size 
	double weight; 
	
	// interval 
	double wMin, wMax; 
	
	// weights 
	Vec <double> w; 
	
	// pre forces 
	Vec <TInput> pre; 
	
	// rho 
	double rhoplus, rhominus; 
	
	// step size in percentage of weight 
	void construct (double weight, int size = 1, double max = 2.0, double min = EPSILON, double rhoplus = 1.5, double rhominus = 0.5) 
	{ 
		// init params 
		this->weight = weight; 
		this->wMin = min; 
		this->wMax = max; 
		
		// adaption parameters 
		this->rhoplus = rhoplus; 
		this->rhominus = rhominus; 
		
		// resize arrays 
		w.resize(size); 
		pre.resize(size); 
	} 
	
	// reset values 
	void initialize () 
	{ 
		w.fill(weight); 
		pre.fill(TInput (0.0)); 
	} 
}; 

// broyden-fletcher-goldfarb-shanno (bfgs) method 
class StepBFGS 
{ 
	// parameters 
	double grad, fLength; 
	
	// parent 
	IStep < Vec3 > p; 
	public: 
	
	void construct(double weight) 
	{ 
		p.construct(weight); 
	} 
	
	void initialize() 
	{ 
		p.initialize(); 
	} 
	
	inline double adjust(Vec3& f) 
	{ 
		// calc. gradient 
		grad = p.pre[0].scalar(f); 
		
		// switch 
		if (grad > 0.0) 
		p.w[0] = min(p.w[0] * p.rhoplus, p.wMax); 
		else 
		if (grad < 0.0) 
		p.w[0] = max(p.w[0] * p.rhominus, p.wMin); 
		
		// backup f 
		p.pre[0] = f; 
		
		// return 
		return p.w[0]; 
	} 
}; 

// rprop 
class StepRPROP 
{ 
	// parameters 
	double grad; 
	
	// direction 
	int f1; 
	
	// interface 
	IStep < double > p; 
	
	public : 
	void construct(double weight, int size = 1) 
	{ 
		p.construct(weight, size); 
	} 
	
	void initialize() 
	{ 
		p.initialize(); 
	} 
	
	inline double adjust(double& f, int i) 
	{ 
		// setup sign 
		f1 = Lib::sign(f); 
		
		// calc. gradient 
		grad = f1 * p.pre[i]; 
		
		// rprop 
		if (grad > 0) 
		{ 
			// increase learning rate 
			p.w[i] = min(p.w[i] * p.rhoplus, p.wMax); 
			p.pre[i] = f1; 
		} else if (grad < 0) 
		{ 
			// actual shift in direction 
			p.w[i] = max(p.w[i] * p.rhominus, p.wMin); 
			p.pre[i] = 0; 
		} else 
		{ 
			// no direction or lastly shift in direction 
			p.pre[i] = f1; 
		} 
		
		// result 
		return f1 * p.w[i]; 
	} 
}; 

// force 
template < typename TPotential > 
class ForceSingle 
{ 
	// molecule 
	SpecMolecule* mol; 
	SpecMolecule* target; 
	
	// force set 
	ForceSet fs; 
	
	// potential function 
	TPotential potential; 
	
	// minimizer 
	StepBFGS adAngle; 
	StepRPROP adForce; 
	
	// parameters 
	double epsilon, weight; 
	int maxSteps; 
	public: 
	// topologies 
	ForceTopo ft; 
	
	// constructor 
	template < typename TPar > 
	void construct (SpecMolecule* target, SpecMolecule* mol, TPar* par) 
	{ 
		construct(target, mol); 
	} 
	
	// constructor 
	void construct (SpecMolecule* target, SpecMolecule* mol, double epsilon = 0.001, double weight = 0.5, int maxSteps = 50) 
	{ 
		// parameters 
		this->target = target; 
		this->mol = mol; 
		
		// potential 
		potential.construct(&fs, target, mol); 
		
		// load parameters 
		this->epsilon = epsilon; 
		this->weight = weight; 
		this->maxSteps = maxSteps; 
		
		// construct adjuster 
		adAngle.construct(weight); 
		adForce.construct(weight, 3); 
	} 
	
	// improve positioning 
	void minimize() 
	{ 
		// initialize force 
		initialize(); 
		
		// apply force 
		int steps = maxSteps; 
		while (steps > 0 && !next()) 
		{ 
			Trans::apply (mol, fs); 
			steps--; 
		} 
	} 
	
	// improve positioning 
	void minimizeall() 
	{ 
		// initialize force 
		initialize(); 
		
		// apply force 
		int steps = maxSteps; 
		while (steps > 0 && !nextall()) 
		{ 
			Trans::apply (mol, fs); 
			steps--; 
		} 
	} 
	
	// return energy 
	double energy() 
	{ 
		return fs.energy; 
	} 
	private: 
	// get force 
	int nextall() 
	{ 
		// reset force params 
		fs.reset(); 
		
		// calc center of gravity 
		fs.cog = mol->cog(); 
		
		// evaluate 
		for(int i = 0; i < target->lcalpha.size(); i++) 
		for (int j = 0; j < mol->lcalpha.size(); j++) 
		potential.evaluate(target->lcalpha[i], mol->lcalpha[j]); 
		
		// setup 
		fs.angle = adAngle.adjust(fs.amom); 
		
		// setup force 
		for (int i = 0; i < 3; i++) 
		fs.force[i] = adForce.adjust(fs.force[i], i); 
		
		// criteria 
		if (fs.force.length() < epsilon && fs.amom.length() < epsilon) 
		return true; 
		else 
		return false; 
	} 
	
	// get force 
	int next() 
	{ 
		// reset force params 
		fs.reset(); 
		
		// calc center of gravity 
		fs.cog = mol->cog(); 
		
		// evaluate 
		for (int k = 0; k < ft.a.size(); k++) 
		for (int j = 0; j < ft.a[k].size(); j++) 
		for(int i = 0; i < ft.b[k].size(); i++) 
		potential.evaluate(ft.a[k][j], ft.b[k][i]); 
		
		// setup 
		fs.angle = adAngle.adjust(fs.amom); 
		
		// setup force 
		for (int i = 0; i < 3; i++) 
		fs.force[i] = adForce.adjust(fs.force[i], i); 
		
		// criteria 
		if (fs.force.length() < epsilon && fs.amom.length() < epsilon) 
		return true; 
		else 
		return false; 
	} 
	
	void initialize() 
	{ 
		// initialize step manager 
		adForce.initialize(); 
		adAngle.initialize(); 
	} 
}; 

// alignement potential - laurenz function 
class AlignmentPotential 
{ 
	// molecules 
	SpecMolecule *target, *mol; 
	
	// force result set 
	ForceSet* fs; 
	
	// parameters 
	Vec3 v; 
	double r; 
	public: 
	
	// constructor 
	void construct (ForceSet* fs, SpecMolecule* target, SpecMolecule* mol) 
	{ 
		// parameters 
		this->target = target; 
		this->mol = mol; 
		this->fs = fs; 
	} 
	
	// get force 
	inline void evaluate(int i, int j) 
	{ 
		// get distance 
		v = target->latom[i].pos - mol->latom[j].pos; 
		
		// force 
		r = v.length(); 
		v /= (0.01 + r*r); 
		
		// force = sum partForce 
		fs->force += v; 
		fs->amom += v.cross(mol->latom[j].pos - fs->cog); 
	} 
}; 


// alignment topology 
class AlignmentTopology 
{ 
	
	// molecule 
	SpecMolecule* target; 
	SpecMolecule* mol; 
	double cutdist; 
	bool allowinversion; 
	
	// force topology 
	ForceTopo* ft; 
	public: 
	// constructor 
	void construct (SpecMolecule* target, SpecMolecule* mol, ForceTopo* ft, double cutdist, bool allowinversion) 
	{ 
		// initialize 
		this->target = target; 
		this->mol = mol; 
		this->ft = ft; 
		this->cutdist = cutdist; 
		this->allowinversion = allowinversion; 
	} 
	
	// generate topology upon alignment 
	template < typename TVec > 
	double alignmentRestrictive(TVec& map, double threshold = 3.0) 
	{ 
		// mapped on target index 
		int index; 
		
		// counter 
		int i, j, k; 
		
		// assignment parameters 
		int lai, smi; 
		bool go; 
		
		// distance parameters 
		double mind, r; 
		int mini = 0; 
		int minj = 0; 
		
		// small and large 
		SSE *smSSE, *laSSE; 
		
		// ensures unique residues 
		Vec < bool > molresidues, targetresidues; 
		molresidues.resize(mol->lcalpha.size()); 
		molresidues.fill(false); 
		targetresidues.resize(target->lcalpha.size()); 
		targetresidues.fill(false); 
		
		// resize 
		Vec < int > pairMol; 
		Vec < int > pairTarget; 
		
		// direction 
		double direction; 
		int inverse = 0; 
		int sgn; 
		
		// over all genes 
		for (i = 0 ; i < map.size(); i++) 
		{ 
			// check mapping 
			index = map[i]; 
			
			// skip gaps 
			if (index < 0 || index >= target->lsse.size()) 
			continue; 
			
			// select sse 
			smSSE = &mol->lsse[i]; 
			laSSE = &target->lsse[index]; 
			
			// 
			// determine direction 
			// 
			
			// verify direction 
			sgn = 1; 
			
			// sequential alignment direction 
			direction = (target->latom[target->lcalpha[laSSE->getStart()]].pos.dist(mol->latom[mol->lcalpha[smSSE->getStart()]].pos) + 
			target->latom[target->lcalpha[laSSE->getEnd()]].pos.dist(mol->latom[mol->lcalpha[smSSE->getEnd()]].pos)) - 
			(target->latom[target->lcalpha[laSSE->getEnd()]].pos.dist(mol->latom[mol->lcalpha[smSSE->getStart()]].pos) + 
			target->latom[target->lcalpha[laSSE->getStart()]].pos.dist(mol->latom[mol->lcalpha[smSSE->getEnd()]].pos)); 
			
			// direction 
			if (direction > threshold) 
			{ 
				if (allowinversion) 
				{ 
					sgn = -1; 
					inverse++; 
				} else 
				continue; 
			} 
			
			// 
			// get crossing point 
			// 
			
			// get first coordinate 
			mind = LARGE; 
			for (j = laSSE->getStart(); j <= laSSE->getEnd(); j++) 
			for (k = smSSE->getStart(); k <= smSSE->getEnd(); k++) 
			{ 
				r = target->latom[target->lcalpha[j]].pos.dist(mol->latom[mol->lcalpha[k]].pos); 
				if (mind > r) 
				{ 
					mind = r; 
					mini = j; 
					minj = k; 
				} 
			} 
			
			// 
			// assignment 
			// 
			
			// get indices 
			lai = mini; 
			smi = minj; 
			
			// find backward fitting elements 
			while (((lai - sgn) >= 0) && ((lai - sgn) < target->sizeRes()) && ((smi - 1) >= 0)) 
			{ 
				// no double 
				if (molresidues[smi - 1] || targetresidues[lai - sgn]) 
				break; 
				
				// check distance 
				if (target->latom[target->lcalpha[lai - sgn]].pos.dist(mol->latom[mol->lcalpha[smi - 1]].pos) > cutdist) 
				break; 
				
				// update counter 
				lai-=sgn; 
				smi--; 
			} 
			
			// add interaction for sse elements 
			go = true; 
			while(go) 
			{ 
				if (!molresidues[smi] && !targetresidues[lai]) 
				{ 
					molresidues[smi] = true; 
					targetresidues[lai] = true; 
					
					pairTarget.push_back (target->lcalpha[lai]); 
					pairMol.push_back (mol->lcalpha[smi]); 
				} else 
				break; 
				
				// update counter 
				lai+=sgn; 
				smi++; 
				
				// criteria 
				if (lai >= 0 && lai < target->sizeRes() && smi < mol->sizeRes()) 
				{ 
					if (target->latom[target->lcalpha[lai]].pos.dist(mol->latom[mol->lcalpha[smi]].pos) < cutdist) 
					go = true; 
					else 
					go = false; 
				} else 
				go = false; 
			} 
		} 
		
		// setup force topology 
		ft->resize(pairMol.size()); 
		for (i = 0; i < ft->a.size(); i++) 
		{ 
			ft->a[i].push_back( pairTarget[i] ); 
			ft->b[i].push_back( pairMol[i] ); 
		} 
		
		// return inverse 
		return double (inverse) / double (map.size()); 
	} 
	
	// generic topo map 
	template < typename TVec > 
	void map(TVec& map) 
	{ 
		// parameters 
		int index, s, e; 
		
		// resize 
		ft->resize(mol->lsse.size()); 
		
		// over all genes 
		for (int i = 0 ; i < min (map.size(), mol->lsse.size()); i++) 
		{ 
			// check mapping 
			index = map[i]; 
			
			// skip gaps 
			if (index < 0 || index >= target->lsse.size()) 
			continue; 
			
			// add interaction partners in target 
			s = target->lsse[index].getStart(); 
			e = target->lsse[index].getEnd(); 
			for (int j = s; j <= e; j++) 
			ft->a[i].push_back(target->lcalpha[j]); 
			
			// add interaction partners in molecule 
			s = mol->lsse[i].getStart(); 
			e = mol->lsse[i].getEnd(); 
			for (int j = s; j <= e; j++) 
			ft->b[i].push_back(mol->lcalpha[j]); 
		} 
	} 
	
	// get number of pairs 
	int npairs() 
	{ 
		return ft->npairs(); 
	} 
	
	// rmsd upon forcetopology 
	double rmsd() 
	{ 
		double dist = 0; 
		for (int k = 0; k < ft->a.size(); k++) 
		for (int j = 0; j < ft->a[k].size(); j++) 
		for(int i = 0; i < ft->b[k].size(); i++) 
		dist += pow (target->latom[ft->a[k][j]].pos.dist(mol->latom[ft->b[k][i]].pos), 2); 
		
		// return 
		return sqrt (dist / ft->npairs()); 
	} 
}; 

struct CubeItem 
{ 
	int map, index; 
	CubeItem(int map, int index) 
	{ 
		this->map = map; 
		this->index = index; 
	} 
}; 

class Cube 
{ 
	// const 
	static const int GAP = -2; 
	
	// 3d world projection 
	Vec < Vec < Vec < Vec < int > > > > matrix; 
	
	// contact matrix 
	Vec < Vec < int > > contacts; 
	
	// detection 
	int detection; 
	
	// temporary indices 
	int i, j, k, x, y, z; 
	
	// temporary map variables 
	int tSSE, mSSE, v, c; 
	
	// coordinate 
	int coord[3]; 
	
	// boundary 
	double box[6]; 
	
	// dimensions 
	int dim[3]; 
	
	// molecules 
	SpecMolecule* mol; 
	SpecMolecule* target; 
	
	public: 
	// pair list 
	Vec < int > tList, mList; 
	
	// construct 3d cube 
	void construct (SpecMolecule *target, SpecMolecule *mol, int detection) 
	{ 
		// initialize 
		this->target = target; 
		this->mol = mol; 
		this->detection = detection; 
		
		// reset boundary 
		for (i = 0; i < 6; i++) 
		box[i] = 0; 
		
		// get boundery coordinates 
		SpecAtom* a; 
		for (i = 0; i < target->lcalpha.size(); i++) 
		if (i < target->lcalpha.size()) 
		{ 
			// get atom 
			a = &target->latom[target->lcalpha[i]]; 
			
			// get values 
			box[0] = min (box[0], a->pos.x - 0.5); 
			box[1] = max (box[1], a->pos.x + 0.5); 
			box[2] = min (box[2], a->pos.y - 0.5); 
			box[3] = max (box[3], a->pos.y + 0.5); 
			box[4] = min (box[4], a->pos.z - 0.5); 
			box[5] = max (box[5], a->pos.z + 0.5); 
		} 
		
		// get dimensions 
		for (i = 0; i < 3; i++) 
		dim[i] = (int) (box[2 * i + 1] - box[2 * i]) + 1; 
		
		// setup cube 
		matrix.resize(dim[0]); 
		for (x = 0; x < dim[0]; x++) 
		{ 
			matrix[x].resize(dim[1]); 
			for (y = 0; y < dim[1]; y++) 
			{ 
				matrix[x][y].resize(dim[2]); 
				for (z = 0; z < dim[2]; z++) 
				matrix[x][y][z].clear(); 
			} 
		} 
		
		// initialize 
		for (i = 0; i < target->lcalpha.size(); i++) 
		{ 
			// get atom 
			a = &target->latom[target->lcalpha[i]]; 
			
			// verify 
			if (a->sseno == SpecMolecule::SSEEMPTY) 
			continue; 
			
			// select coordinate 
			for (j = 0; j < 3; j++) 
			coord[j] = (int) (a->pos[j] - box[2 * j]); 
			
			// write values in cut off 
			for (x = coord[0] - detection; x < coord[0] + detection; x++) 
			for (y = coord[1] - detection; y < coord[1] + detection; y++) 
			for (z = coord[2] - detection; z < coord[2] + detection; z++) 
			{ 
				if (validate (x, y, z)) 
				matrix[x][y][z].push_back(a->sseno); 
			} 
		} 
		
		// setup contact matrix 
		contacts.resize(mol->lsse.size()); 
		for (i = 0; i < contacts.size(); i++) 
		contacts[i].resize(target->lsse.size()); 
	} 
	
	// check if index is valid 
	inline bool validate (int x, int y, int z) 
	{ 
		return (x >= 0 && y >= 0 && z >= 0) && (x < dim[0] && y < dim[1] && z < dim[2]); 
	} 
	
	// evaluate alignment 
	template < typename TVec > 
	bool evaluate(TVec& map, bool update) 
	{ 
		// reset contacts 
		for (i = 0; i < contacts.size(); i++) 
		contacts[i].fill(0); 
		
		// projection of target backbone into cube 
		for (i = 0; i < mol->lcalpha.size(); i++) 
		{ 
			// get atom 
			SpecAtom* a = &mol->latom[mol->lcalpha[i]]; 
			
			// select coordinate 
			c = -1; 
			for (j = 0; j < 3; j++) 
			{ 
				c = (int) (a->pos[j] - box[2 * j]); 
				if ( (c < 0) || (c >= dim[j]) ) 
				{ 
					c = -1; 
					break; 
				} 
				
				// backup 
				coord[j] = c; 
			} 
			
			// skip 
			if (c == -1) 
			continue; 
			
			// copy coordinates 
			x = coord[0]; 
			y = coord[1]; 
			z = coord[2]; 
			
			// is inside cube 
			if (validate (x, y, z)) 
			{ 
				for (int k = 0; k < matrix[x][y][z].size(); k++) 
				{ 
					// indices 
					tSSE = matrix[x][y][z][k]; 
					mSSE = mol->latom[mol->lcalpha[i]].sseno; 
					
					// is assigned 
					if (mSSE == SpecMolecule::SSEEMPTY) 
					continue; 
					
					// core type 
					if (target->lsse[tSSE].getType() != mol->lsse[mSSE].getType()) 
					continue; 
					
					// count 
					contacts[mSSE][tSSE]++; 
				} 
			} 
		} 
		
		// 
		// create map 
		// 
		
		// setup map 
		map.resize(mol->lsse.size()); 
		map.fill(GAP); 
		
		// setup map 
		for (i = 0; i < map.size(); i++) 
		{ 
			// find maximum 
			tSSE = mSSE = c = 0; 
			for (j = 0; j < contacts.size(); j++) 
			{ 
				// search in target 
				Vec < int >::const_iterator iter = std::max_element (contacts[j].begin(), contacts[j].end()); 
				
				// backup maximum 
				if (*iter > c) 
				{ 
					mSSE = j; 
					tSSE = iter - contacts[j].begin(); 
					c = *iter; 
				} 
			} 
			
			// finalized 
			if (c < 3) 
			break; 
			
			// remove conflicting entries 
			for (j = 0; j < mol->lsse.size(); j++) 
			contacts[j][tSSE] = 0; 
			for (j = 0; j < target->lsse.size(); j++) 
			contacts[mSSE][j] = 0; 
			
			// append 
			map[mSSE] = tSSE; 
		} 
		
		// update map 
		if (update) 
		mapupdate(map); 
		
		// verify map 
		for (i = 0; i < map.size(); i++) 
		if (map[i] >= 0) 
		return true; 
		
		// return status 
		return false; 
	} 
	
	// update map 
	void mapupdate (Vec<int>& map) 
	{ 
		// copy gapless 
		list < CubeItem > gapless; 
		for (i = 0; i < map.size(); i++) 
		if (map[i] >= 0) 
		gapless.push_back(CubeItem(map[i], i)); 
		
		// verify 
		bool jp, jn, jg; 
		Vec < int > lgaps; 
		
		// generate iterators 
		list<CubeItem>::iterator it, ip, in, idelete; 
		
		// list 
		it = gapless.begin(); 
		while (true) 
		{ 
			// next 
			it++; 
			
			// verify 
			if (it == gapless.end()) 
			break; 
			
			// get previous item 
			ip = it; 
			ip--; 
			
			// get next item 
			in = it; 
			in++; 
			
			// verify 
			if (in == gapless.end()) 
			break; 
			
			// get flags 
			jp = (ip->map > it->map) ? true : false; 
			jn = (it->map > in->map) ? true : false; 
			jg = (ip->map > in->map) ? true : false; 
			
			// is break? 
			if ((jp || jn) && !jg) 
			{ 
				lgaps.push_back(it->index); 
				idelete = it; 
				it--; 
				gapless.erase(idelete); 
			} 
		} 
		
		// copy back 
		for (i = 0; i < lgaps.size(); i++) 
		map[lgaps[i]] = GAP; 
	} 
}; 

// pymol script writer 
class PyMol 
{ 
	// parameters 
	SpecMolecule *target, *mol; 
	
	// pymol 
	fstream pms; 
	string path; 
	string fname; 
	string coldefa, coldefb; 
	string colmarka, colmarkb; 
	
	// code 
	char chain1; 
	char chain2; 
	public: 
	// colors: lightblue, density, lightorange, orange 
	// constructor 
	void construct (SpecMolecule* target, SpecMolecule* mol, string path = "", bool numeric = false) 
	{ 
		this->target= target; 
		this->mol = mol; 
		this->path = path; 
		this->fname = target->name + "_" + mol->name; 
		
		// change color 
		coldefa = "density"; 
		colmarka = "orange"; 
		coldefb = "lightblue"; 
		colmarkb = "lightorange"; 
		
		// set chain names 
		if (!numeric) 
		{ 
			chain1 = 'A'; 
			chain2 = 'B'; 
		} else { 
			chain1 = '0'; 
			chain2 = '1'; 
		} 
	} 
	
	// returns the full filename of the pdb file 
	string getfilename() 
	{ 
		return path + fname + ".pdb"; 
	} 
	
	// write all 
	void writeinterface (Vec<double>& lrefa, Vec<double>& lrefb, Vec<Vec3>& lpreda, double tresh = 4.5) 
	{ 
		// open 
		open(); 
		reset(); 
		
		// verify size 
		if (lrefa.size() != target->lcalpha.size() || lrefb.size() != mol->lcalpha.size() || lpreda.size() > target->lcalpha.size()) 
		Msg::error("PyMol::writeinterface()", "Size mismatch!"); 
		
		// mark predicted ones 
		Vec < bool > lflags; 
		lflags.resize(lrefa.size()); 
		lflags.fill(0.0); 
		for (int i = 0; i < target->lcalpha.size(); i++) 
		for (int j = 0; j < lpreda.size(); j++) 
		if (target->latom[target->lcalpha[i]].pos == lpreda[j]) 
		lflags[i] = true; 
		
		// select 
		for (int i = 0; i < lrefa.size(); i++) 
		{ 
			string strcol = ""; 
			
			// it is an interface residue 
			if (lrefa[i] < tresh) 
			{ 
				// true positive or false negative 
				if (lflags[i]) 
				strcol = "green"; 
				else 
				strcol = colmarka; 
			} else { 
				// false positive 
				if (lflags[i]) 
				strcol = "red"; 
			} 
			residue(i, i, chain1, "color " + strcol); 
		} 
		
		// select 
		for (int i = 0; i < lrefb.size(); i++) 
		if (lrefb[i] < tresh) 
		residue(i, i, chain2, "color " + colmarkb); 
		
		// close 
		close(); 
	} 
	
	// write all 
	void writecolors () 
	{ 
		open(); 
		reset(); 
		close(); 
	} 
	
	// write all 
	template < typename TVec > 
	void writemap (TVec& m) 
	{ 
		// initialize 
		open(); 
		
		// over map 
		int s, e; 
		for (int k = 0; k < m.size(); k++) 
		{ 
			// skip gaps 
			if (! (m[k] >= 0) ) 
			continue; 
			
			s = target->lsse[m[k]].getStart(); 
			e = target->lsse[m[k]].getEnd(); 
			residue(s, e, chain1, "color " + colmarka); 
			sse(s, e-1, target->lsse[m[k]].getType(), chain1); 
			
			s = mol->lsse[k].getStart(); 
			e = mol->lsse[k].getEnd(); 
			residue(s, e, chain2, "color " + colmarkb); 
			sse(s, e-1, mol->lsse[k].getType(), chain2); 
		} 
		
		// close file 
		close(); 
	} 
	
	// pymol output only interaction partners 
	template < typename TTopo > 
	void write (TTopo* ft) 
	{ 
		// initialize 
		open(); 
		
		// define all sses 
		reset(); 
		
		// over all a 
		int index, first, j, i; 
		for (i = 0 ; i < ft->a.size(); i++) 
		{ 
			j = 0; 
			while (j < ft->a[i].size()) 
			{ 
				// first index 
				first = ft->a[i][j]; 
				index = first; 
				
				// write sse 
				while(target->latom[first].sse == target->latom[ ft->a[i][j] ].sse) 
				{ 
					// get index 
					index = ft->a[i][j]; 
					residue (target->latom[index].resno, chain1, "color " + colmarka); 
					
					// counter 
					j++; 
					if (j >= ft->a[i].size()) 
					break; 
				} 
				
				// write sse 
				sse (target->latom[first].resno, target->latom[index].resno, target->latom[first].sse, chain1); 
			} 
		} 
		
		for (i = 0 ; i < ft->b.size(); i++) 
		{ 
			j = 0; 
			while (j < ft->b[i].size()) 
			{ 
				// first index 
				first = ft->b[i][j]; 
				index = first; 
				
				// write sse 
				while(mol->latom[first].sse == mol->latom[ ft->b[i][j] ].sse) 
				{ 
					// get index 
					index = ft->b[i][j]; 
					residue (mol->latom[index].resno, chain2, "color " + colmarkb); 
					
					// counter 
					j++; 
					if (j >= ft->b[i].size()) 
					break; 
				} 
				
				// write sse 
				sse (mol->latom[first].resno, mol->latom[index].resno, mol->latom[first].sse, chain2); 
			} 
		} 
		
		// close file 
		close(); 
	} 
	
	private: 
	
	// write all 
	void reset () 
	{ 
		// over all sses 
		int s, e; 
		for (int i = 0 ; i < target->lsse.size(); i++) 
		{ 
			s = target->lsse[i].getStart(); 
			e = s + target->lsse[i].length(); 
			if (s < e) 
			{ 
				for (int j = s; j < e; j++) 
				residue(j, chain1, "color " + coldefa); 
				sse(s, e-1, target->lsse[i].getType(), chain1); 
			} 
		} 
		
		for (int i = 0 ; i < mol->lsse.size(); i++) 
		{ 
			s = mol->lsse[i].getStart(); 
			e = s + mol->lsse[i].length(); 
			if (s < e) 
			{ 
				for (int j = s; j < e; j++) 
				residue (j, chain2, "color " + coldefb); 
				sse(s, e-1, mol->lsse[i].getType(), chain2); 
			} 
		} 
	} 
	
	// 
	// pymol basics 
	// 
	void open() 
	{ 
		// write file 
		string fn = string(path + fname + ".pml"); 
		pms.open(fn.c_str(),ios_base::out); 
		pms << "load " << fname << ".pdb" << endl; 
		pms << "bg_color white" << endl; 
		pms << "hide everything" << endl; 
		pms << "show spheres, chain H" << endl; 
		pms << "show cartoon" << endl; 
		pms << "color " << coldefa << ", chain " << chain1 << endl; 
		pms << "color " << coldefb << ", chain " << chain2 << endl; 
	} 
	
	
	// close file 
	void close(double sphere = 0.5, double stick = 0.35) 
	{ 
		//pms << "show surface" << endl; 
		//pms << "set transparency = 0.75" << endl; 
		pms << "set sphere_scale=" << sphere << endl; 
		pms << "set stick_radius=" << stick << endl; 
		
		// finalize 
		pms.close(); 
	} 
	
	// colors 
	inline static string getcolor(int i) 
	{ 
		// switch 
		static string colors[] = {"yellow","marine","green","lightmagenta","cyan","density","red", 
		"forest","deepolive","lime","sand","deeppurple","palecyan","orange"}; 
		i = i % 13; 
		return colors[i]; 
	} 
	
	// write residue 
	void residue(int start, int end, char chain, string item) 
	{ 
		pms << "show cartoon, resi " << start << "-" << end << " AND chain " << chain << endl; 
		pms << item << ", resi " << start << "-" << end << " AND chain " << chain << endl; 
	} 
	
	// write residue 
	void residue(int residue, char chain, string item) 
	{ 
		pms << "show cartoon, resi " << residue << " AND chain " << chain << endl; 
		pms << item << ", resi " << residue << " AND chain " << chain << endl; 
	} 
	
	// write sse 
	void sse(int start, int end, char type, char chain) 
	{ 
		if (type == 'E') type = 'S'; 
		pms << "alter " << chain << "/" << min (start, end) << "-" << max (start, end) << "/, ss='" << type << "'" << endl; 
	} 
}; 


// analytical optimization 
class AlignmentScheme 
{ 
	
	// parameters 
	Config par; 
	
	public: 
	// molecule 
	SpecMolecule* target; 
	SpecMolecule* mol; 
	
	// topology and cube 
	AlignmentTopology topo; 
	Cube cube; 
	
	// single body translation and minimization with potential 
	ForceSingle <AlignmentPotential> force; 
	
	// count inverse elements 
	double inverse; 
	
	// constructor 
	void construct (SpecMolecule* target, SpecMolecule* mol) 
	{ 
		// initialize 
		this->target = target; 
		this->mol = mol; 
		this->inverse = 0.0; 
		
		// load force 
		force.construct(target, mol, &par); 
		
		// topology generator 
		topo.construct(target, mol, &force.ft, par.alResidueDistance, par.alInversion); 
		
		// cube constructor 
		cube.construct(target, mol, par.alCoreDistance); 
	} 
	
	// optimize directly from coordinates 
	bool optimize() 
	{ 
		// create new map 
		Vec <int> map; 
		
		// get residue assignment 
		return load(map); 
	} 
	
	// optimize by correct map 
	template < typename TVec > 
	bool optimize(TVec& map) 
	{ 
		// crude align 
		topo.map(map); 
		force.minimize(); 
		
		// return 
		if (load(map)) 
		return true; 
		else 
		return false; 
	} 
	
	// remaster template 
	template < typename TVec > 
	bool load(TVec& map, int nrepeat = 10) 
	{ 
		// repeat 
		for (int i = 0; i < nrepeat; i++) 
		{ 
			// remaster map 
			if (!cube.evaluate(map, !par.alInversion)) 
			return false; 
			
			// get residue based topology 
			topo.alignmentRestrictive(map); 
			
			// finalize alignment by kabsch 
			finalize(); 
		} 
		
		// get residue based topology 
		inverse = topo.alignmentRestrictive(map); 
		
		// return 
		return true; 
	} 
	
	private: 
	// prepare optimization 
	template < typename TVec > 
	void prepare(TVec& map) 
	{ 
		// number of pairs 
		Vec < Vec3 > tcalpha, mcalpha; 
		for (int c = 0; c < min ((int) mol->lsse.size(), (int) map.size()); c++) 
		{ 
			// skip gaps 
			if (map[c] < 0 || map[c] >= target->lsse.size()) 
			continue; 
			
			// add pairs 
			tcalpha.push_back( target->latom[ target->lssecenter[map[c]] ].pos ); 
			mcalpha.push_back( mol->latom [ mol->lssecenter[c] ].pos ); 
		} 
		
		// apply alignment 
		Trans::align(mcalpha, tcalpha, mol); 
	} 
	
	// apply final kabsch 
	void finalize() 
	{ 
		// list 
		Vec < Vec3 > lista, listb; 
		
		// copy target atoms to lista 
		lista.resize(force.ft.a.size()); 
		for (int i = 0; i < force.ft.a.size(); i++) 
		lista[i] = target->latom[force.ft.a[i][0]].pos; 
		
		// copy target atoms to listb 
		listb.resize(force.ft.b.size()); 
		for (int i = 0; i < force.ft.b.size(); i++) 
		listb[i] = mol->latom[force.ft.b[i][0]].pos; 
		
		// apply alignment 
		Trans::align(listb, lista, mol); 
	} 
}; 

// bags 
struct Bags 
{ 
	// bags 
	Vec < int > alpha, beta; 
	
	// construct 
	void construct (SpecMolecule* mol) 
	{ 
		// clear bags 
		alpha.clear(); 
		beta.clear(); 
		
		// create bags 
		for (int i = 0; i < mol->lsse.size(); i++) 
		{ 
			if (mol->lsse[i].getType() == SpecMolecule::SSEALPHA) 
			alpha.push_back(i); 
			else 
			beta.push_back(i); 
		} 
	} 
	
	// shuffle 
	inline void shuffle() 
	{ 
		Lib::shuffle(alpha); 
		Lib::shuffle(beta); 
	} 
}; 



// map features 
struct MapFeatures 
{ 
	
	// verify non gaps 
	template < typename TVec > 
	static inline bool verifyNonGaps(TVec& m, double minsse = 0.5) 
	{ 
		// count gaps 
		int nongaps = MapFeatures::countNonGaps(m); 
		if (nongaps >= int (m.size() * minsse) && nongaps > 2) 
		return true; 
		else 
		return false; 
	} 
	
	// check sequential 
	template < typename TVec > 
	static inline int countNonGaps(TVec& m) 
	{ 
		int nongaps = 0; 
		for (int i = 0; i < m.size(); i++) 
		{ 
			if (m[i] != -2) 
			nongaps++; 
		} 
		return nongaps; 
	} 
	
	// check sequential 
	template < typename TVec > 
	static inline bool makeSequential(TVec& m) 
	{ 
		// verify 
		if (m.size() == 0) 
		return false; 
		
		// copy 
		Vec < int > data, index; 
		for (int i = 0; i < m.size(); i++) 
		{ 
			// gap 
			if (m[i] < 0) 
			continue; 
			
			// copy 
			data.push_back(m[i]); 
			index.push_back(i); 
		} 
		
		// count 
		Vec < int > length, start; 
		
		// check for single break 
		for (int j = 0; j < data.size() - 1; j++) 
		{ 
			// count 
			int cx = 0; 
			for (int i = j; i < data.size() - 1; i++) 
			{ 
				// is a break 
				if (data[i] < data[i + 1]) 
				cx--; 
				else 
				break; 
			} 
			length.push_back(cx); 
			start.push_back(j); 
		} 
		
		// sort 
		Lib::sort(length, start); 
		
		if (start.size() > 0) 
		{ 
			for (int i = 0; i < start[0]; i++) 
			m[index[i]] = -2; 
			
			for (int i = start[0] + (-1) * length[0] + 1; i < index.size(); i++) 
			m[index[i]] = -2; 
			
			int cx = 0; 
			for (int i = 0; i < m.size(); i++) 
			if (m[i] >= 0) 
			cx++; 
			
			// return 
			if (cx > 0) 
			return true; 
		} 
		
		return false; 
	} 
	
	// check sequential 
	template < typename TVec > 
	static inline bool isSequential(TVec& m) 
	{ 
		// verify 
		if (m.size() == 0) 
		return false; 
		
		// copy 
		Vec < int > data; 
		for (int i = 0; i < m.size(); i++) 
		{ 
			// gap 
			if (m[i] < 0) 
			continue; 
			
			// copy 
			data.push_back(m[i]); 
		} 
		
		// check for single break 
		for (int i = 0; i < data.size() - 1; i++) 
		{ 
			// is a break 
			if (data[i] > data[i + 1]) 
			return false; 
		} 
		
		// return 
		return true; 
	} 
	
	// check circular permutation 
	template < typename TVec > 
	static inline bool isCircular(TVec& m) 
	{ 
		// verify 
		if (m.size() == 0) 
		return false; 
		
		// copy 
		Vec < int > data; 
		for (int i = 0; i < m.size(); i++) 
		{ 
			// gap 
			if (m[i] < 0) 
			continue; 
			
			// copy 
			data.push_back(m[i]); 
		} 
		
		// get size 
		int nsize = data.size(); 
		
		// check for single break 
		bool flag = false; 
		for (int i = 0; i < nsize - 1; i++) 
		{ 
			// is a break 
			if (data[i] > data[i + 1]) 
			{ 
				// verify no previous break 
				if (flag) 
				return false; 
				
				// set flag 
				flag = true; 
			} 
		} 
		
		// return 
		if (flag && data[0] > data[nsize - 1]) 
		return true; 
		else 
		return false; 
	} 
	
	// log 
	template < typename TVec > 
	static void lg(TVec& m, int count = 0, double score = 0.0) 
	{ 
		string data = ""; 
		for (int j = 0; j < m.size(); j++) 
		{ 
			// set data 
			if (m[j] >= 0) 
			data += Convert::toString (m[j]); 
			else 
			data += "G"; 
			
			// set space 
			if (j+1 < m.size()) 
			data += " "; 
		} 
		
		// write 
		if (count != 0) 
		Msg::write ("%s\t[%f, %i]", data.c_str(), score, count); 
		else 
		Msg::write ("%s", data.c_str()); 
	} 
}; 



// mapping 
struct Map 
{ 
	// const 
	static const int GAP = -2; 
	static const int ZERO = -LARGE; 
	
	// data 
	Vec < int > data; 
	
	// scores 
	double score; 
	
	// count non gaps 
	int count; 
	
	// check used values 
	Vec < bool > check; 
	
	// help to fill sequential 
	int index; 
}; 

// map link 
struct MapLink 
{ 
	// data 
	Map** ptr; 
	int start, n; 
	
	// construct 
	template < typename MapBase > 
	void construct(MapBase* mb, int start, int n) 
	{ 
		this->ptr = (Map**) &mb->maps[start]; 
		this->start = start; 
		this->n = n; 
	} 
}; 

// map ram database 
class MapBase 
{ 
	public: 
	// data 
	SpecMolecule *target, *mol; 
	
	// count valid scores 
	int counter; 
	
	// maps 
	Vec < Map* > maps; 
	
	// amount of items 
	int n, tn; 
	
	// construct 
	bool construct (SpecMolecule* target, SpecMolecule* mol) 
	{ 
		return construct(target, mol, target->lsse.size(), mol->lsse.size()); 
	} 
	
	// construct 
	bool construct (SpecMolecule* target, SpecMolecule* mol, int tn, int n) 
	{ 
		// initialize 
		this->target = target; 
		this->mol = mol; 
		
		// clear maps 
		clear(); 
		
		// amount 
		this->n = n; 
		this->tn = tn; 
		
		// reset counter 
		counter = 0; 
		
		// return 
		return true; 
	} 
	
	// clear 
	void clear() 
	{ 
		// remove mappings 
		for (int i = 0; i < maps.size(); i++) 
		delete maps[i]; 
		
		// delete index 
		maps.clear(); 
	} 
	
	// destruct 
	~MapBase() 
	{ 
		clear(); 
	} 
	
	/** 
	MAP OPERATIONS 
	**/ 
	// construct 
	inline Map* mapping() 
	{ 
		// construct 
		Map* m = new Map; 
		
		// resize data 
		m->data.resize(n); 
		
		// reset flags 
		m->check.clear(); 
		m->check.resize(tn); 
		
		// clear 
		clear (m); 
		
		// return 
		return m; 
	} 
	
	// copy 
	inline Map* mapping (Map* copy) 
	{ 
		// new map 
		Map* m = mapping(); 
		
		// setup 
		set(m, copy); 
		
		// return 
		return m; 
	} 
	
	// append mapping in new copy 
	inline Map* mapping (Map* copy, int i) 
	{ 
		// create 
		Map* m = mapping(copy); 
		
		// add new element 
		append (m, i); 
		
		// return 
		return m; 
	} 
	
	// copy 
	template < typename TVec > 
	inline void set (Map* m, TVec& copy) 
	{ 
		// clear 
		clear(m); 
		
		// setup 
		for (int i = 0; i < copy.size(); i++) 
		append (m, copy[i]); 
	} 
	
	// set 
	inline void set (Map* m, Map* copy) 
	{ 
		// initialize 
		m->count = copy->count; 
		m->score = copy->score; 
		m->index = copy->index; 
		
		// copy data 
		m->data.assign(copy->data.begin(), copy->data.end()); 
		m->check.assign(copy->check.begin(), copy->check.end()); 
	} 
	
	// clear 
	inline void clear (Map* m) 
	{ 
		// initialize 
		m->index = 0; 
		m->count = 0; 
		m->score = 0; 
		
		// resize data 
		m->data.fill(Map::GAP); 
		
		// reset flags 
		m->check.fill(false); 
	} 
	
	// add 
	inline void append (Map* m, int i) 
	{ 
		// skip gaps 
		if (i < 0) 
		m->data[m->index++] = Map::GAP; 
		else if(!m->check[i]) 
		{ 
			// alter values 
			m->data[m->index++] = i; 
			
			// update 
			m->count++; 
			m->check[i] = true; 
		} 
	} 
	
	// update 
	inline void update (Map* m, int index, int value) 
	{ 
		// clear old 
		if (m->data[index] >= 0) 
		{ 
			m->check[m->data[index]] = false; 
			m->count--; 
		} 
		
		// set new 
		if (value >= 0) 
		{ 
			m->data[index] = value; 
			m->check[value] = true; 
			m->count++; 
		} else 
		m->data[index] = Map::GAP; 
	} 
	
	// finished 
	inline bool finalize (Map* m, int limit = 2) 
	{ 
		// count values 
		if (m->count >= limit) 
		{ 
			// fill reset with gaps 
			for (int i = m->index; i < n; i++) 
			m->data[i] = Map::GAP; 
			m->index = n; 
			
			// return 
			return true; 
		} else 
		return false; 
	} 
	
	// return true if maps are equal 
	inline bool equal (Map* m, Map* m0) 
	{ 
		for (int i = 0; i < n; i++) 
		{ 
			if (m->data[i] != m0->data[i]) 
			return false; 
		} 
		
		// return 
		return true; 
	} 
	
	// number of inequal assignments 
	inline int distance (Map* m, Map* m0) 
	{ 
		int cx = 0; 
		for (int i = 0; i < n; i++) 
		{ 
			if (m->data[i] != m0->data[i]) 
			cx++; 
		} 
		
		// return 
		return cx; 
	} 
	
	/** 
	Map base 
	*/ 
	
	// add mapping to base 
	inline void add (Map* m) 
	{ 
		// counter 
		if (m->score > (double) Map::ZERO) 
		counter++; 
		
		// add to list 
		maps.push_back(m); 
	} 
	
	// size 
	inline int size() 
	{ 
		return maps.size(); 
	} 
	
	/** 
	Sort 
	*/ 
	// sort 
	inline void sort () 
	{ 
		sort (maps); 
	} 
	
	inline void sort (Vec <Map*>& maps) 
	{ 
		sort (maps, maps.size()); 
	} 
	
	inline void sort(Vec <Map*>& maps, int size) 
	{ 
		sort ( (Vec <Map*>*) &maps[0], size); 
	} 
	
	// qsort map 
	inline void sort(Vec <Map*>* maps, int size) 
	{ 
		qsort(maps, size, sizeof(Map*), &compare); 
	} 
	
	// quicksort compare by score 
	static int compare(const void* p, const void* p0) 
	{ 
		// evaluate 
		double value = (* (Map**) p0)->score - (* (Map**) p)->score; 
		if (value < MINUSEPSILON) 
		return -1; 
		else if (value > EPSILON) 
		return 1; 
		else 
		return 0; 
	} 
	
	// mapling sorter 
	inline void sort(MapLink& mapLink) 
	{ 
		// sort 
		std::sort ( maps.begin() + mapLink.start, maps.begin() + mapLink.start + mapLink.n, compareGreater); 
		
		// relink data 
		mapLink.ptr = (Map**) &maps[mapLink.start]; 
	} 
	
	// compare by score 
	static bool compareGreater(const Map* m, const Map* m0) 
	{ 
		return m->score > m0->score; 
	} 
	
	// log 
	void lg(int amount = 0) 
	{ 
		// loop 
		amount = (amount == 0) ? maps.size() : (min (amount, (int) maps.size())); 
		for (int i = 0; i < amount; i++) 
		MapFeatures::lg(maps[i]->data, maps[i]->count, maps[i]->score); 
	} 
}; 




// cache 
struct Cache 
{ 
	Vec < Vec < int > > icc, dist; 
	Vec < Vec < char > > type; 
	Vec < Vec < bool > > cmap; 
}; 

// information common among all individuals 
class Info 
{ 
	
	// molecules 
	SpecMolecule* mol; 
	
	public: 
	// map cache 
	Cache cache; 
	
	// type descriptors 
	static const char TYPENULL = '0'; 
	static const char TYPECROSS = 'X'; 
	static const char TYPEPAR = 'P'; 
	static const char TYPEANTI = 'A'; 
	
	// number of cores 
	int nCores; 
	
	// parameters 
	double indistmax, inrescale; 
	
	// construct 
	bool construct(SpecMolecule* mol) 
	{ 
		GPlusConfig par; 
		return construct (mol, par.inDistMax, par.inRescale); 
	} 
	
	// construct 
	bool construct(SpecMolecule* mol, double indistmax, double inrescale) 
	{ 
		// initialize 
		this->mol = mol; 
		this->indistmax = indistmax; 
		this->inrescale = inrescale; 
		
		// get number of cores 
		nCores = mol->lsse.size(); 
		
		// get residue contact map 
		if (!tableContacts()) 
		return false; 
		
		// create maps 
		if (!tableICC()) 
		return false; 
		
		// create maps 
		tableType(); 
		tableDist(); 
		
		// return 
		return true; 
	} 
	
	private: 
	// get contact map 
	bool tableContacts() 
	{ 
		// parameters 
		double dist; 
		int n = mol->lcalpha.size(); 
		
		// resize 
		cache.cmap.resize(n); 
		for (int i = 0; i < n; i++) 
		cache.cmap[i].resize(n); 
		
		// fill 
		int cx = 0; 
		for (int i = 0; i < n; i++) 
		for (int j = i + 1; j < n; j++) 
		{ 
			// get distance 
			dist = mol->latom[mol->lcalpha[i]].pos.dist(mol->latom[mol->lcalpha[j]].pos); 
			
			// check distance 
			if (dist <= indistmax) 
			{ 
				cache.cmap[i][j] = true; 
				cache.cmap[j][i] = true; 
				cx++; 
			} else { 
				cache.cmap[i][j] = false; 
				cache.cmap[j][i] = false; 
			} 
		} 
		
		
		// total number of contacts 
		if (cx == 0) 
		{ 
			Msg::write("Protein has no residue contacts."); 
			return false; 
		} else 
		Msg::write("Protein has %i residue contacts.", cx); 
		
		// return 
		return true; 
	} 
	
	// initialize the inter core contact map 
	int tableICC() 
	{ 
		// allocate memory and initialize 
		cache.icc.resize(nCores); 
		
		// initialize space 
		for (int i = 0; i < nCores; i++) 
		{ 
			cache.icc[i].resize(nCores); 
			for (int j = 0; j < nCores; j++) 
			cache.icc[i][j] = 0; 
		} 
		
		// fill data 
		int sumCont = 0; 
		for(int i = 0; i < nCores; i++) 
		for (int j = i+1; j < nCores; j++) 
		for(int r1= mol->lsse[i].getStart(); r1 <= mol->lsse[i].getEnd(); r1++) 
		for(int r2= mol->lsse[j].getStart(); r2 <= mol->lsse[j].getEnd(); r2++) 
		if (cache.cmap[r1][r2]) 
		{ 
			cache.icc[i][j]++; 
			cache.icc[j][i]++; 
			sumCont++; 
		} 
		
		// total number of contacts 
		if (sumCont == 0) 
		{ 
			Msg::write("Protein has no secondary structure contacts."); 
			return false; 
		} 
		
		// rescale strands 
		if (inrescale > 1.0) 
		{ 
			// for all cores 
			for (int i = 0; i < nCores; i++) 
			{ 
				// we use only beta-strands 
				if (mol->lsse[i].getType() != SpecMolecule::SSEBETA) 
				continue; 
				for (int j = i+1; j < nCores; j++) 
				cache.icc[i][j] = cache.icc[j][i] = (int) ceil(cache.icc[i][j] * inrescale); 
			} 
		} 
		
		// return 
		return true; 
	} 
	
	// create distance map 
	bool tableDist() 
	{ 
		// allocate space 
		cache.dist.resize(nCores); 
		for (int i = 0; i < nCores; i++) 
		cache.dist[i].resize(nCores); 
		
		// calculate the distances 
		for (int i = 0; i < nCores; i++) 
		{ 
			cache.dist[i][i] = 0; 
			for (int j = i+1; j < nCores; j++) 
			cache.dist[i][j] = cache.dist[j][i] = 
			(int) round(mol->latom[mol->lssecenter[i]].pos.dist(mol->latom[mol->lssecenter[j]].pos) * 10); 
		} 
		
		// return 
		return true; 
	} 
	
	// create type map 
	bool tableType() 
	{ 
		// allocate space 
		cache.type.resize(nCores); 
		for (int i = 0; i < nCores; i++) 
		cache.type[i].resize(nCores); 
		
		// for all cores 
		for (int i = 0; i < nCores; i++) 
		{ 
			// initialize 
			cache.type[i][i] = TYPENULL; 
			
			// loop 
			for (int j = i+1; j < nCores; j++) 
			{ 
				// initialize 
				int sumMax = -LARGE; 
				int sumMin = LARGE; 
				int diffMax = -LARGE; 
				int diffMin = LARGE; 
				int distType = 0; 
				
				// determine direction 
				for (int c1 = mol->lsse[i].getStart(); c1 <= mol->lsse[i].getEnd(); c1++) 
				{ 
					for (int c2 = mol->lsse[j].getStart(); c2 <= mol->lsse[j].getEnd(); c2++) 
					{ 
						if (cache.cmap[c1][c2]) 
						{ 
							int sum = c1 + c2; 
							int diff = c1 - c2; 
							
							if (sum>sumMax) sumMax = sum; 
							if (sum<sumMin) sumMin = sum; 
							if (diff>diffMax) diffMax = diff; 
							if (diff<diffMin) diffMin = diff; 
						} 
					} 
					
					// setup 
					distType = (sumMax-sumMin)-(diffMax-diffMin); 
					if (distType > 1) 
					{ 
						cache.type[i][j] = TYPEPAR; 
						cache.type[j][i] = TYPEPAR; 
					} else if (distType < -1) 
					{ 
						cache.type[i][j] = TYPEANTI; 
						cache.type[j][i] = TYPEANTI; 
					} else if (sumMin < 999999) 
					{ 
						cache.type[i][j] = TYPECROSS; 
						cache.type[j][i] = TYPECROSS; 
					} else { 
						cache.type[i][j] = TYPENULL; 
						cache.type[j][i] = TYPENULL; 
					} 
				} 
			} 
		} 
		
		// return 
		return true; 
	} 
}; 




// evaluation 
class ObjGplus 
{ 
	// data 
	SpecMolecule *target, *mol; 
	
	// parameters 
	int ldiff, distSum, numTM, numGaps, lastGene, seqGaps, numDQ, badDist; 
	double singlediff1, singlediff2, diff, lmalus, typeMatch; 
	
	// size 
	int n; 
	
	// gap penalty 
	double gapPenalty; 
	public: 
	// info 
	Info infoTarget, infoMol; 
	
	// construct 
	bool construct (SpecMolecule* target, SpecMolecule* mol) 
	{ 
		// initialize 
		this->target = target; 
		this->mol = mol; 
		
		// load info object 
		if(!infoMol.construct (mol)) 
		return false; 
		
		if (!infoTarget.construct (target)) 
		return false; 
		
		// set size 
		n = mol->lsse.size(); 
		
		// gap penalty 
		gapPenalty = 0.3 / max ( (mol->lsse.size() * 0.5), 2.0 ); 
		
		// return 
		return true; 
	} 
	
	// get score 
	inline void get (Map* m) 
	{ 
		// parameters 
		lastGene = -1; 
		distSum = ldiff = numTM = numGaps = seqGaps = numDQ = badDist = 0; 
		singlediff1 = singlediff2 = diff = lmalus = typeMatch = 0.0; 
		
		// loop over cores 
		for (int i = 0; i < n; i++) 
		{ 
			// count gaps 
			if (m->data[i] == Map::GAP) 
			{ 
				numGaps++; 
				continue; 
			} 
			
			// check sequential 
			if (lastGene > m->data[i]) 
			seqGaps++; 
			lastGene = m->data[i]; 
			
			// judge length difference 
			ldiff = mol->lsse[i].length() - target->lsse[m->data[i]].length(); 
			if (mol->lsse[i].getType() == SpecMolecule::SSEALPHA) 
			{ 
				// for alpha helices larger length differences are acceptable 
				if (ldiff > 5) 
				lmalus += 0.01; 
				else if ((ldiff < 0) && (ldiff > -7)) 
				lmalus += 0.005; 
				else if (ldiff < -7) 
				lmalus += 0.03; 
			} else { 
				if (ldiff > 3) 
				lmalus += 0.01; 
				else if ((ldiff < 0) && (ldiff > -3)) 
				lmalus += 0.005; 
				else if ((ldiff < -2) && (ldiff > -4)) 
				lmalus += 0.01; 
				else if (ldiff < -3) 
				lmalus += 0.03; 
			} 
			
			// check coremaps 
			distSum = 0; 
			numDQ = 0; 
			for (int j = i; j < n; j++) 
			{ 
				// skip gap 
				if (m->data[j] == Map::GAP) 
				continue; 
				
				// core contact differences 
				diff += abs(infoMol.cache.icc[i][j] - infoTarget.cache.icc[m->data[i]][m->data[j]]); 
				singlediff1 += infoMol.cache.icc[i][j]; 
				singlediff2 += infoTarget.cache.icc[m->data[i]][m->data[j]]; 
				distSum += abs(infoMol.cache.dist[i][j] - infoTarget.cache.dist[m->data[i]][m->data[j]]); 
				numDQ++; 
				
				// evaluate sse orientation 
				if ((infoMol.cache.icc[i][j] > 0) && (infoTarget.cache.icc[m->data[i]][m->data[j]] > 0)) 
				{ 
					numTM++; 
					if (infoMol.cache.type[i][j] == infoTarget.cache.type[m->data[i]][m->data[j]]) 
					// exact match 
					typeMatch += 1.0; 
					else if (( infoMol.cache.type[i][j] == Info::TYPEPAR || 
					infoMol.cache.type[i][j] == Info::TYPEANTI) && 
					( infoTarget.cache.type[m->data[i]][m->data[j]] == Info::TYPEPAR || 
					infoTarget.cache.type[m->data[i]][m->data[j]] == Info::TYPEANTI )) 
					// anti-parallel match 
					typeMatch += 0.7; 
					else if (( infoMol.cache.type[i][j] == Info::TYPECROSS) || 
					( infoTarget.cache.type[m->data[i]][m->data[j]] == Info::TYPECROSS )) 
					// cross-cross match 
					typeMatch += 0.5; 
					else 
					// other 
					typeMatch += 0.3; 
				} 
			} 
			
			// judge distance 
			if (distSum / numDQ > 100) 
			badDist++; 
		} 
		
		// reset score 
		m->score = 0.0; 
		
		// sequence bonus 
		m->score += 0.11 / ((double) seqGaps + 1); 
		
		// sum contributions 
		m->score -= (badDist + numGaps) * gapPenalty; 
		m->score -= lmalus; 
		
		if (singlediff1 + singlediff2 > 0) 
		m->score += (1.0 - diff / (singlediff1 + singlediff2)) * 0.6; 
		
		if (numTM > 0) 
		m->score += (typeMatch / ((double) numTM)) * 0.4; 
		
		// ensure range 
		if (m->score < 0.0) 
		m->score = (double) Map::ZERO; 
		else if (m->score > 1.0) 
		m->score = 1.0; 
	} 
}; 



// leaf 
struct Leaf 
{ 
	Vec <Leaf*> next; 
	Leaf(int n) 
	{ 
		next.resize(n); 
		for (int i = 0; i < n; i++) 
		next[i] = 0; 
	} 
}; 

// tree 
class Tree 
{ 
	// maximal distance 
	static const int MAXDIST = 5; 
	
	// leafs 
	Leaf *head, *index; 
	
	// counter 
	int n, id, cx; 
	public: 
	
	Tree() 
	{ 
		// initialize 
		head = NULL; 
		index = NULL; 
		n = 0; 
	} 
	
	~Tree() 
	{ 
		clean(head); 
	} 
	
	// clean 
	void clean(Leaf* index) 
	{ 
		// index 
		if (index != NULL) 
		{ 
			for (int i = 0; i < n; i++) 
			clean(index->next[i]); 
			delete index; 
		} 
	} 
	
	void construct (int n) 
	{ 
		// reset 
		clean(head); 
		
		// for a gap (+ 1) 
		this->n = n + 1; 
		
		// setup main leaf 
		head = new Leaf(this->n); 
	} 
	
	bool search(vector <int>& map) 
	{ 
		index = head; 
		cx = 0; 
		for (int i = 0; i < map.size(); i++) 
		{ 
			if (map[i] != Map::GAP) 
			id = map[i]; 
			else 
			// for a gap (- 1) 
			id = n - 1; 
			
			if (index->next[id] == 0) 
			{ 
				// not found 
				index->next[id] = new Leaf(n); 
				return false; 
			} else { 
				// counter 
				if (++cx >= MAXDIST) 
				return true; 
			} 
			
			// go to next element 
			index = index->next[id]; 
		} 
		
		return true; 
	} 
}; 




// prepare 
template < typename TObjective > 
class Prepare 
{ 
	// maximal number of tupel maps 
	static const int NMAX = 999999; 
	
	// parameters 
	SpecMolecule* target; 
	SpecMolecule* mol; 
	GPlusConfig par; 
	
	// map base 
	MapBase* mb; 
	
	// objective 
	TObjective* obj; 
	
	// sse bags 
	Bags bag; 
	
	// counter 
	int ntotal; 
	int nsse; 
	public: 
	
	// constructor 
	bool construct (MapBase* mb, TObjective* obj) 
	{ 
		// initialize 
		this->mb = mb; 
		this->target= mb->target; 
		this->mol = mb->mol; 
		this->obj = obj; 
		this->nsse = mb->n; 
		
		// reset size 
		ntotal = 0; 
		
		// verify 
		if (target->lsse.size() == 0 || mol->lsse.size() == 0) 
		return false; 
		
		// generate bags 
		bag.construct(target); 
		
		// create mappings 
		generate (mb->mapping()); 
		if (ntotal == 0) 
		{ 
			par.ppCoreDelta = LARGE; 
			generate (mb->mapping()); 
		} 
		
		// sort mapping 
		mb->sort(mb->maps); 
		
		// log 
		if (ntotal > 0) 
		{ 
			Msg::write ("Number of mappings %i. Maps above tresh %i. Highscore %f.", 
			ntotal, mb->counter, mb->maps[0]->score); 
			return true; 
		} else 
		return false; 
	} 
	
	// generator 
	void generate(Map* m) 
	{ 
		// size limit 
		if (ntotal > NMAX) 
		{ 
			delete m; 
			return; 
		} 
		
		// completed 
		if (mb->finalize(m)) 
		{ 
			// get score 
			obj->get(m); 
			mb->add(m); 
			
			// increase counter 
			ntotal++; 
		} else { 
			// is it possible to elongate? 
			if(m->index == nsse) 
			{ 
				delete m; 
				return; 
			} 
			
			// get length of molecule sse 
			int lng = mol->lsse[m->index].length(); 
			
			// check type 
			if (mol->lsse[m->index].getType() == SpecMolecule::SSEALPHA) 
			{ 
				// append sequence 
				for (int i = 0; i < bag.alpha.size(); i++) 
				{ 
					if (!m->check[bag.alpha[i]]) 
					{ 
						// check length 
						if (abs (lng - target->lsse[bag.alpha[i]].length()) < par.ppCoreDelta) 
						generate (mb->mapping(m, bag.alpha[i])); 
					} 
				} 
			} else { 
				// append sequence 
				for (int i = 0; i < bag.beta.size(); i++) 
				{ 
					if (!m->check[bag.beta[i]]) 
					{ 
						// check length 
						if (abs (lng - target->lsse[bag.beta[i]].length()) < par.ppCoreDelta) 
						generate (mb->mapping(m, bag.beta[i])); 
					} 
				} 
			} 
			
			// try gap 
			if (m->index != nsse) 
			generate (mb->mapping(m, Map::GAP)); 
			
			// delete 
			delete m; 
		} 
	} 
}; 




// search engine 
template < typename TObjective> 
class Evaluate 
{ 
	// parameter set 
	MapBase* mb; 
	SpecMolecule* target; 
	SpecMolecule* mol; 
	GPlusConfig par; 
	
	// new map properties 
	Map* m; 
	double s; 
	
	// population 
	Vec <Map*> *pop; 
	
	// counter 
	int nScores, nPop, nDepth; 
	
	// counter 
	int i, ia, ib, j, n; 
	
	// search tree 
	Tree tree; 
	
	// type 
	TObjective* obj; 
	public: 
	// constructor 
	bool construct(MapBase* mb, TObjective* obj) 
	{ 
		// initialize 
		this->mb = mb; 
		this->target= mb->target; 
		this->mol = mb->mol; 
		this->obj = obj; 
		
		// parameters 
		nScores = 0; 
		nPop = 0; 
		
		// initialize population 
		population(); 
		
		// initialize search tree 
		tree.construct(mb->tn); 
		
		// log 
		Msg::write("Evaluation of %i total maps, %i above thresh, %i evaluation depth.", nPop, mb->counter, nDepth); 
		
		// make optimization 
		combinatorics(); 
		
		// verify 
		if (nScores > 0) 
		return true; 
		else 
		return false; 
	} 
	
	// are they compatible 
	inline bool verify (Map* m, Map* input) 
	{ 
		// loop 
		for (ib = 0; ib < input->data.size(); ib++) 
		{ 
			// input is not GAP 
			if (input->data[ib] != Map::GAP) 
			{ 
				// value exists in both mappings at same place 
				if (m->check[input->data[ib]]) 
				if (m->data[ib] != input->data[ib]) 
				return false; 
			} 
		} 
		
		// return 
		return true; 
	} 
	private: 
	// initialize population 
	void population() 
	{ 
		// setup calculation space 
		nPop = min ((int) mb->size(), par.evDepth); 
		
		// search depth 
		nDepth = min (mb->counter, (int) (0.5 * nPop)); 
		
		// link map base 
		pop = &mb->maps; 
	} 
	
	// merge 
	inline void merge (Map* m, Map* input) 
	{ 
		// loop 
		for (ia = 0; ia < input->data.size(); ia++) 
		{ 
			// input is not GAP 
			if (input->data[ia] != Map::GAP) 
			{ 
				// map is GAP 
				if (m->data[ia] == Map::GAP) 
				{ 
					m->data[ia] = input->data[ia]; 
					m->check[m->data[ia]] = true; 
					m->count++; 
				} 
			} 
		} 
		
		// update score 
		if (tree.search(m->data)) 
		m->score = (double) Map::ZERO; 
		else 
		obj->get(m); 
	} 
	
	// next 
	inline bool combinatorics() 
	{ 
		// generate children 
		n = nPop; 
		for (i = 0; i < nDepth; i++) 
		{ 
			for (j = i + 1; j < nDepth; j++) 
			{ 
				// is compatible 
				if (verify((*pop)[i], (*pop)[j])) 
				{ 
					// link map 
					m = (*pop)[n - 1]; 
					
					// merge 
					mb->set(m, (*pop)[i]); 
					merge(m, (*pop)[j]); 
					
					// verify score 
					if (m->score > max ((*pop)[j]->score, (*pop)[i]->score)) 
					if (--n < nDepth) 
					goto iterate; 
				} 
			} 
		} 
		
		// goto tag 
		iterate: 
		
		// sort 
		mb->sort((*pop), nPop); 
		
		// break criteria 
		if (n == nPop) 
		return false; 
		else { 
			// update 
			nScores += nPop - n; 
			
			// log 
			Msg::write("%i evalutation depth, %i successful maps generated.", i, nPop - n); 
			
			// return true 
			return true; 
		} 
	} 
}; 



// result item 
struct OptimizeItem 
{ 
	// data 
	Vec <Vec3> pos; 
	Vec < int > map; 
	double rmsd, inverse, score, globalscore; 
	int pairs, nongaps; 
	bool circular, insequence; 
	ForceTopo topo; 
	
	// reset 
	void construct () 
	{ 
		// clear 
		pos.clear(); 
		map.clear(); 
		topo.clear(); 
		
		// reset 
		rmsd = 0.0; 
		pairs = 0; 
		score = 0; 
		circular = false; 
		insequence = false; 
		nongaps = 0; 
		inverse = 0.0; 
		globalscore = 0.0; 
	} 
	
	// constructor 
	template < typename TVec, typename TTopo > 
	void construct (SpecMolecule* target, SpecMolecule* mol, TVec& map, double rmsd, int pairs, double inverse, TTopo* ft, int normalization = 0, int globalpairs = 0) 
	{ 
		// assert 
		int nalpha = TransMatrix::nalpha; 
		
		// copy coordinates 
		pos.resize(nalpha); 
		for (int i = 0; i < nalpha; i++) 
		pos[i] = mol->latom[mol->lcalpha[i]].pos; 
		
		// topology 
		if (ft != NULL) 
		topo.construct(ft); 
		
		// normalize score 
		int defaultnormalization = min (mol->lcalpha.size(), target->lcalpha.size()); 
		if (normalization == 0) 
		normalization = defaultnormalization; 
		
		// score 
		this->score = double (pairs) / double (normalization); 
		this->globalscore = double (globalpairs) / double (defaultnormalization); 
		
		// link additional data 
		this->map = map; 
		this->rmsd = rmsd; 
		this->pairs = pairs; 
		this->circular = MapFeatures::isCircular(map); 
		this->insequence = MapFeatures::isSequential(map); 
		this->nongaps = MapFeatures::countNonGaps(map); 
		this->inverse = inverse / this->nongaps; 
		
		// log 
		MapFeatures::lg(map); 
		Msg::write("%f rmsd, %i pairs, %f score.", rmsd, pairs, score); 
	} 
}; 


// result collector 
struct Optimize 
{ 
	// types 
	typedef GPlusConfig TPar; 
	typedef OptimizeItem TItem; 
	typedef Map TMap; 
	
	// result pointer 
	TItem* item; 
	
	// results objects 
	TItem itemfirst; 
	TItem itemsecond; 
	
	// swap flag 
	bool swap; 
	
	// raw molecule 
	SpecMolecule raw; 
	SpecMolecule *mol; 
	
	// transformation 
	TransMatrix tmat; 
	
	// apply shift 
	void applytransformation (bool info = false) 
	{ 
		TransMatrix tmat; 
		mol->construct(&raw); 
		tmat.apply(mol, item->pos); 
	} 
	
	// apply shift 
	void gettransformation (TransMatrix& tmat) 
	{ 
		// set reference 
		tmat.construct(&raw); 
		tmat.prepare(item->pos); 
		
		// swap 
		if(!swap) 
		tmat.reverse(); 
	} 
	
	// apply shift 
	void printtransformation () 
	{ 
		TransMatrix tmat; 
		
		// set reference 
		gettransformation(tmat); 
		
		// show 
		tmat.print(); 
	} 
	
	// construct 
	TItem* construct (SpecMolecule*& target, SpecMolecule*& mol) 
	{ 
		// construct result items 
		itemfirst.construct(); 
		itemsecond.construct(); 
		
		// link first item by default 
		item = &itemfirst; 
		
		// swap for smaller maps 
		if (target->lsse.size() < mol->lsse.size()) 
		{ 
			SpecMolecule* tmp = mol; 
			mol = target; 
			target = tmp; 
			swap = true; 
		} else 
		swap = false; 
		
		// link molecule 
		this->mol = mol; 
		
		// copy reference 
		raw.construct(mol); 
		
		// parameters 
		TPar par; 
		
		// objective function 
		ObjGplus obj; 
		if(!obj.construct (target, mol)) 
		return item; 
		
		// additional modules 
		MapBase mb; // maps database 
		Prepare<ObjGplus> pp; // prepare data set tupels 
		Evaluate<ObjGplus> ev; // evaluation strategy 
		
		// maximal amount 
		int maxAmount = par.coResults; 
		
		// create map database 
		if(!mb.construct(target, mol)) 
		{ 
			Msg::write ("Mapbase construction failed."); 
			return item; 
		} 
		
		// log 
		Msg::line(); 
		
		// evaluation 
		if (!pp.construct(&mb, &obj)) 
		{ 
			Msg::write ("Preparation failed."); 
			return item; 
		} 
		
		// optimize 
		ev.construct(&mb, &obj); 
		
		// setup objects 
		AlignmentScheme oa; 
		oa.construct (target, mol); 
		
		// loop parameters 
		int k = 0; 
		int amount = 0; 
		bool valid = true; 
		
		// number of aligned residues 
		int pairs = 0; 
		
		// loop 
		while (mb.maps.size() > k && amount < maxAmount) 
		{ 
			// get map 
			TMap* m = mb.maps[k]; 
			
			// verify with filter for bindom 
			if (valid) 
			{ 
				// optimize 
				if (oa.optimize(m->data)) 
				{ 
					// get pairs 
					pairs = oa.topo.npairs(); 
					if (pairs > itemfirst.pairs) 
					itemfirst.construct (target, mol, m->data, oa.topo.rmsd(), pairs, oa.inverse, &oa.force.ft); 
					if (pairs > itemsecond.pairs && MapFeatures::isSequential(m->data)) 
					itemsecond.construct (target, mol, m->data, oa.topo.rmsd(), pairs, oa.inverse, &oa.force.ft); 
					
					// get sequential structure alignment 					
					MapFeatures::makeSequential(m->data); 
					oa.topo.alignmentRestrictive(m->data); 
					
					// get pairs 
					pairs = oa.topo.npairs(); 
					if (pairs > itemfirst.pairs) 
					itemfirst.construct (target, mol, m->data, oa.topo.rmsd(), pairs, oa.inverse, &oa.force.ft); 
					if (pairs > itemsecond.pairs && MapFeatures::isSequential(m->data)) 
					itemsecond.construct (target, mol, m->data, oa.topo.rmsd(), pairs, oa.inverse, &oa.force.ft);
				} 
				
				// counter 
				amount++; 
			} 
			
			// update counter 
			k++; 
			
			// verify 
			if (mb.maps.size() > k) 
			{ 
				if (!mb.equal(mb.maps[k], mb.maps[k-1])) 
				valid = true; 
				else 
				valid = false; 
			} 
		} 
		
		// next 
		Msg::line(); 
		
		// return 
		return item; 
	} 
}; 


// main module 
struct ModAlign 
{ 
	// storage 
	Storage store; 
	
	// setup molecules 
	SpecMolecule targetObj, molObj; 
	
	// construct 
	double construct (string ida, string idb) 
	{ 
		// load maps of ida on idb 
		store.read(&targetObj, ida, FName::getname(ida)); 
		
		// construct 
		return construct(&targetObj, idb); 
	} 
	
	// construct 
	double construct (SpecMolecule* target, string idb) 
	{ 
		// load molecule 
		store.read(&molObj, idb, FName::getname(idb)); 
		
		// link molecule 
		SpecMolecule *mol = &molObj; 
		
		// reset score 
		double score = 0.0; 
		
		// verify loaded molecules 
		if(target->lsse.size() > 0 && mol->lsse.size() > 0) 
		{ 
			// output 
			printf ("Chain 1: %15s      Size= %i\n", target->name.c_str(), (int) target->lcalpha.size()); 
			printf ("Chain 2: %15s      Size= %i\n", mol->name.c_str(), (int) mol->lcalpha.size()); 
			
			// log 
			target->lg(); 
			mol->lg(); 
			
			// optimization scheme 
			Optimize optimize; 
			OptimizeItem* o = optimize.construct(target, mol); 
			
			// verify 
			if (o->score > 0.0) 
			{ 
				// backup score 
				score = o->score; 
				
				// check sequential 
				cout << "Features: "; 
				if (o->insequence) 
				cout << "sequential, "; 
				else 
				cout << "non-sequential, "; 
				
				// is circular 
				if (o->circular) 
				cout << "circ.perm., "; 
				else 
				cout << "no-perm., "; 
				
				// inverse elements 
				if (o->inverse > 0.0) 
				{ 
					string tag = "inv." + Convert::toString(int (o->inverse * 100)) + "%"; 
					cout << tag.c_str(); 
				} else 
				cout << "no-inver."; 
				cout << endl; 
				
				// result 
				printf ("Aligned length=%4i, RMSD=%6.2f, GP-Score=%6.5f", o->pairs, o->rmsd, o->score); 
				printf ("\n"); 
				
				cout << endl << " -------- rotation matrix to rotate Chain-1 to Chain-2 ------" << endl; 
				
				// write transformation 
				optimize.printtransformation(); 
				
				#ifdef PDB_ON 
				// apply on molecule coordinates 
				optimize.applytransformation(); 
				
				// pymol output 
				PyMol pmol; 
				pmol.construct(target, mol); 
				pmol.write(&o->topo); 
				
				// save results 
				store.save(target, mol, pmol.getfilename()); 
				#endif 
			} 
			
			// completed 
			cout << endl << "GANGSTA+ Done."; 
		} else { 
			Msg::error ("Files not found", "Could not read pdb file(s)."); 
		}
		
		// log 
		Msg::line(); 
		
		// return score 
		return score; 
	} 
}; 

// main 
int main(int n, char* arg[]) 
{ 
	// show details 
	Msg::title(); 
	
	// result 
	cout << " **************************************************************************" << endl; 
	cout << " *                              GANGSTA+                                  *" << endl; 
	cout << " *                                                                        *" << endl; 
	cout << " * Reference: A. Guerler and E. W. Knapp, Protein Sci. 2008 17, 1374-82   *" << endl; 
	cout << " * Comments on the program, please email to: aysam.guerler@gmail.com      *" << endl; 
	cout << " **************************************************************************" << endl; 
	cout << endl; 
	
	// check args 
	if (n < 3) 
	Msg::error ("Missing information", "[pdb file] [pdb file]"); 
	
	// run 
	ModAlign spec; 
	spec.construct((string) arg[1], (string) arg[2]); 
	
	// return 
	return 0; 
} 
