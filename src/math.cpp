#include "math.hpp"

/*
// �������������� ���������:
const ld M_E = std::exp(1.0L);
const ld M_PI = std::acos(-1.0L);
const ld M_PI_2 = M_PI / 2.0L;
*/

// ��������� ������ � �����:
std::ostream& operator<<(std::ostream& os, const Vector& v) {
	for (auto& it : v) {
		os << std::setprecision(8) << std::fixed << std::setw(16) << it;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
	for (auto& it : m) {
		os << it << std::endl;
	}
	return os;
}

// �������� ��������� ������� �� ����� � ��������� �� ��������:
Vector& operator*=(Vector& v, const ld value) {
	for (auto & it : v) {
		it *= value;
	}
	return v;
}

Vector operator*(const ld value, Vector v) {
	return v *= value;
}

Vector operator*(Vector v, const ld value) {
	return v *= value;
}

// ������� ����� ��� �������:
Vector operator-(Vector v) {
	return v *= -1;
}

// �������� � ��������� ��������:
Vector& operator+=(Vector & v1, const Vector& v2) {
	assert(v1.size() == v2.size());
	for (int i = 0; i < (int)v1.size(); ++i) {
		v1[i] += v2[i];
	}
	return v1;
}

Vector& operator-=(Vector & v1, const Vector& v2) {
	return v1 += -v2;
}

Vector operator+(Vector v1, const Vector& v2) {
	return v1 += v2;
}

Vector operator-(Vector v1, const Vector& v2) {
	return v1 -= v2;
}

// ��������� ������������ ��������:
ld dot(const Vector& v1, const Vector& v2) {
	assert(v1.size() == v2.size());
	ld sum = 0;
	for (int i = 0; i < (int)v1.size(); ++i) {
		sum += v1[i] * v2[i];
	}
	return sum;
}

// ����� �������:
ld norm(const Vector& v) {
	return std::sqrt(dot(v, v));
}
//________________________________________________



ld MAX(const ld &a, const ld &b) {
	return b > a ? (b) : (a);
}

 ld SIGN(const ld &a, const ld &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

ld SQR(const ld a) { return a * a; }
//_________________________________________________

// ��������� ������� �� ������:
Vector operator*(const Matrix& m, const Vector& v) {
	int nRows = (int)m.size();
	Vector ans(nRows);
	for (int i = 0; i < nRows; ++i) {
		ans[i] = dot(m[i], v);
	}
	return ans;
}

// �������� ��������� ������� ����
bool is_zero(const Vector& v) {
	for (const auto & it : v) {
		if (std::abs(it) > COMPARE_EPS) {
			return false;
		}
	}
	return true;
}

// �������� ������ � ������������ R^n:
// �� i-�� ����� 1, �� ���� ��������� ����
Vector id_vect(int size, int i) {
	Vector answer(size);
	answer[i] = 1;
	return answer;
}

// ����� ������ ��������� �� �������:
Vector grad(Function f, const Vector& point, const ld h) {
	int n = (int)point.size();
	Vector right = point;
	Vector left = point;
	Vector answer(n);
	for (int i = 0; i < n; ++i) {
		// ������� ���������� i-�� ���������:
		right[i] += h;
		left[i] -= h;
		// �� ������� ����������� ��������� � ������������ O(h^2):
		answer[i] = (f(right) - f(left)) / (2 * h);
		// ������� ���������� �������:
		right[i] -= h;
		left[i] += h;
	}
	return answer;
}

// ����� ���������� ������� ����� (������� ��������, ������� ������ ��� �������� �����������)
// � ����� x � ����� h � ������������ O(h)
Matrix hess(Function f, const Vector& x, const ld h) {
	int n = (int)x.size();
	Matrix ans(n, Vector(n));
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			ans[i][j] = (
				f(x + (id_vect(n, i) + id_vect(n, j))*h)
				- f(x + id_vect(n, i)*h)
				- f(x + id_vect(n, j)*h)
				+ f(x)) / h / h;
		}
	return ans;
}
