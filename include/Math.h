#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cassert>
//#include <nrtypes.h>  

// ���������� �����:

typedef long double ld;
typedef std::vector<ld> Vector;
typedef std::vector<double> Vector1;  
typedef std::vector<Vector> Matrix;
typedef ld(*Function)(const Vector & x);
/*
template<class T>
inline const T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
/*
template<class T>
inline void SWAP(T &a, T &b)
{
	T dum = a; a = b; b = dum;
}
/*
template<class T>
inline const T MAX(const T &a, const T &b)
{
	return b > a ? (b) : (a);
}

template<class T>
inline const T SQR(const T a) { return a * a; }
*/
typedef std::pair<Vector, int>(*Method)(Function, Vector, int);
//typedef std::pair<Vector, int> powell(int *a, DP func(Vec_I_DP &x), int NDIM, int n);

// �������� ��������� ������������ �����:
const ld COMPARE_EPS = 0.0000000000000001L;

// �������������� ��������� (M_E, M_PI, M_PI_2) � ����� ����������

// ��������� ������ � �����:
std::ostream& operator<<(std::ostream& os, const Vector& v);
std::ostream& operator<<(std::ostream& os, const Matrix& m);

// �������� ��������� ������� �� ����� � ��������� �� ��������:
Vector& operator*=(Vector& v, const ld value);
Vector operator*(const ld value, Vector v);
Vector operator*(Vector v, const ld value);

// ������� ����� ��� �������:
Vector operator-(Vector v);

// �������� � ��������� ��������:
Vector& operator+=(Vector & v1, const Vector& v2);
Vector& operator-=(Vector & v1, const Vector& v2);
Vector operator+(Vector v1, const Vector& v2);
Vector operator-(Vector v1, const Vector& v2);

// ��������� ������������ ��������:
ld dot(const Vector& v1, const Vector& v2);

// ����� �������:
ld norm(const Vector& v);
/*
template<class T>
inline const T MAX(const T &a, const T &b)
{
	return b > a ? (b) : (a);
}
*/

ld MAX(const ld &a, const ld &b);

void SWAP(ld &a, ld &b);

ld SIGN(const ld &a, const ld &b);

ld SQR(const ld a);

// ��������� ������� �� ������:
Vector operator*(const Matrix& m, const Vector& v);

// �������� ��������� ������� ����
bool is_zero(const Vector& v);

// �������� ������ � ������������ R^n:
// �� i-�� ����� 1, �� ���� ��������� ����
Vector id_vect(int size, int i);

// ����� ������ ��������� �� �������:
Vector grad(Function f, const Vector& point, const ld h = 1e-4);

// ����� ���������� ������� ����� (������� ��������, ������� ������ ��� �������� �����������)
// � ����� x � ����� h � ������������ O(h)
Matrix hess(Function f, const Vector& x, const ld h = 1e-4);