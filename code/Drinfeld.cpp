#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_pE.h>

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
//#include <math.h>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace NTL;



// Generic function definitions

int ilog(double base, double targ) {
	return (int) (ilogb(targ)/ ilogb(base));
}

void operator_eval(ZZ_pE& res, vector<ZZ_pE>& op, ZZ_pE& elem, long q ) {
	clear(res);
	vector<ZZ_pE> powers;
	powers.resize(op.size());
	powers[0] = elem;
	for (int i = 0; i < op.size(); i++) {
		if (i < op.size() - 1) {
			power(powers[i+1], powers[i], q);
		}
		res += powers[i] * op[i];
	}
}


void build_poly(ZZ_pX& poly, int* arr, int len) {
	poly.SetLength(len);
	for (int i = 0; i < len; i++) {
		poly[i] = conv<ZZ_p>(ZZ(arr[i]));
	}
	return;
}

void elem_exp(mat_ZZ_pE& ret, mat_ZZ_pE& matr, double expo ) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			power(ret[i][j], matr[i][j], expo);
		}
	}
	return;
}



ZZ_pX char_poly1(ZZ_pE& g, ZZ_pE& del, long p, long q_exp, int n, ZZ_pX& P) {
	int m = n, cnum = n/2 + 1, rnum = n, q = pow(p, q_exp);
	int tot = rnum + cnum;
	ZZ_pX cons = ZZ_pX(INIT_MONO, 1, 1);
	ZZ_pE mono = conv<ZZ_pE>(cons);

	vec_ZZ_p x;
	ZZ_p det;
	vector<ZZ_pE> phit;
	phit.resize(3);

	ZZ_pE alpha = conv<ZZ_pE>(ZZ(1));
	set(alpha);


	phit[0] = mono;
	phit[1] = g;
	phit[2] = del;


	vec_ZZ_pE rvals;
	vec_ZZ_p trvals;
	rvals.SetLength(rnum);
	trvals.SetLength(rnum);

	ZZ_pX b = pow(-1,m) * norm(del) * P;

	vec_ZZ_p b_vec;
	b_vec.SetLength(rnum);


	ZZ_pE phi_b_alpha;
	clear(phi_b_alpha);
	vector<ZZ_pE> phis;
	phis.resize(tot);
	phis[0] = alpha;



	for (int i = 1; i< tot; i++) {
		operator_eval(phis[i], phit, phis[i-1], q);
	}

	for (int i = 0; i < rnum + 1; i++) {
		phi_b_alpha += phis[i]*b[i];
	}
	ZZ_pE r = alpha + phi_b_alpha;

	rvals[0] = r;
	trvals[0] = trace(r);

	for (int i = 1; i < rnum; i++) {
		operator_eval(rvals[i], phit, rvals[i-1], q);
		trvals[i] = trace(rvals[i]);
	}
	mat_ZZ_p A;
	A.SetDims(rnum, cnum);
	set(A[0][0]);

	for (int i = 0; i < cnum; i++) {
		A[0][i] = trace(phis[i]);
	}

	for (int i = 1; i < rnum; i++) {
		for (int j = 0; j < cnum - 1; j++) {
			A[i][j] = A[i-1][j + 1];
		}
		A[i][cnum - 1] = trace(phis[cnum - 1 + i]);
	}
	x.SetLength(cnum);

	for (int i = 0; i < rnum; i++) {
		b_vec[i] = trvals[i];
	}

	b_vec.SetLength(cnum);
	A.SetDims(cnum,cnum);

	double begin = GetTime();
	solve(det, A, x, b_vec);
	double end = GetTime();
	cout << "Hankel System solving time: " << end-begin << endl;

	ZZ_pX poly;
	poly.SetLength(x.length());
	vector<ZZ_pE> skewpoly;
	for (int i = 0; i < x.length(); i++) {
		poly[i] = x[i];

	}

	return poly;
}

ZZ_pEX char_poly2(ZZ_pE g, ZZ_pE del, long p, long q_exp, int n, ZZ_pX P) {

	int m = n, cnum = n/2 + 1, rnum = n, q = pow(p, q_exp);
	int tot = rnum + cnum;
	ZZ_pX cons = ZZ_pX(INIT_MONO, 1, 1);
	ZZ_pE mono = conv<ZZ_pE>(cons);
	ZZ_pE alpha, beta;

	vec_ZZ_pE ei;
	vec_ZZ_pE rvals;
	ZZ_pEX interpol;

	ZZ_pX elem;
	ei.SetLength(cnum);
	rvals.SetLength(cnum);

	int maxdeg = ilog(p, cnum);
	elem.SetLength(maxdeg);

	for (int i = 0; i < cnum; i++) {
		ei[i] = conv<ZZ_p>(ZZ(i));
		alpha = (-1)*(mono - ei[i]) / del;
		beta = (-1)*g / del;

		mat_ZZ_pE M, B, N;
		N.SetDims(2,2);
		B.SetDims(2,2);
		M.SetDims(2,2);
		M[0][0] = 0;
		M[1][0] = 1;
		M[0][1] = alpha;
		M[1][1] = beta;

		B = M;

		int ind = 0;

		while(2*ind + 1 < n) {
			elem_exp(B, M, pow(q, ind + 1));
			mul(M,M,B);
			(ind*=2)++;

		}
		elem_exp(B, B, pow(q, ind));
		while (ind + 1 < n) {
			elem_exp(B,B,q);
			mul(M,M,B);
			ind++;
		}

		if (M[1][0] != 0) {
			ZZ_pE temp;
			power(temp, M[0][0], q);
			rvals[i] = M[0][0] + temp;
			power(temp, M[1][0], q);
			rvals[i] += beta*temp;
		}
		else rvals[i] = M[0][0];


	}
	interpolate(interpol, ei, rvals);

	return interpol;

}



int main() {
	int p = 17, q_exp = 1, n = 6;


	ZZ_p::init(ZZ(p));
	ZZ_pX P;
	BuildIrred(P, n);
	ZZ_pE::init(P);

	ZZ_pE g, del;
	set(g);
	set(del);

	ZZ_pX out = char_poly1(g,del,p,q_exp,n,P);
	ZZ_pEX out2 = char_poly2(g,del,p,q_exp,n,P);

	cout << "Char poly (random): " << out << endl << "(deterministic): " << out2 << endl;
}
