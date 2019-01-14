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

void operator_eval(ZZ_pE& res, ZZ_pE& op, ZZ_pE& elem, long q ) {
	clear(res);
	vector<ZZ_pE> powers;
	ZZ_pX repp = rep(op);
	long d = deg(repp);
	if (d == -1) {
		res = 0;
		return;
	}
	powers.resize(d + 1);
	powers[0] = 1;
	for (int i = 0; i <= d; i++) {
		if (i < d) {
			mul(powers[i+1], powers[i], elem);
		}
		res += powers[i] * repp[i];
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

void elem_exp(mat_ZZ_pE& ret, mat_ZZ_pE& matr, double q, double expo ) {
	ZZ_pX cons = ZZ_pX(INIT_MONO, 1, 1);
	ZZ_pE mono = conv<ZZ_pE>(cons);
	ZZ_pE xq;
	ZZ_pE temp;
	power(xq, mono, q);
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			if (matr[i][j] == 0) {
				ret[i][j] = 0;
			}
			else {
				ret[i][j] = matr[i][j];
				for (int k = 0; k < expo; k++) {
					operator_eval(temp, ret[i][j], xq, 1 );
					ret[i][j] = temp;
				}
			}
		}
	}
	return;
}


void frob(ZZ_pE& ret, ZZ_pE& matr, double q, double expo ) {
	cout << "entering" << endl;
	ZZ_pX cons = ZZ_pX(INIT_MONO, 1, 1);
	ZZ_pE mono = conv<ZZ_pE>(cons);
	ZZ_pE xq;
	ZZ_pE temp;
	power(xq, mono, q);
			if (matr == 0) {
				ret = 0;
			}
			else {
				ret = matr;
				for (int k = 0; k < expo; k++) {
					operator_eval(temp, ret, xq, 1 );
					ret = temp;
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
	int tot = rnum + cnum, ind = 0;
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
		ei[i] = conv<ZZ_pE>(ZZ(i));
		alpha = (-1)*(mono - ei[i]) / del;
		beta = (-1)*g / del;
		ind = 0;

		mat_ZZ_pE M, B, N, temp;
		N.SetDims(2,2);
		B.SetDims(2,2);
		M.SetDims(2,2);
		M[0][0] = 0;
		M[1][0] = 1;
		M[0][1] = alpha;
		M[1][1] = beta;
		B = M;
		N = M;

		while(2*ind + 1 < n) {
			elem_exp(B, M, q, ind + 1);
			mul(temp,M,B);
			M = temp;
			(ind*=2)++;
			//cout << M << endl;
		}

		elem_exp(B, N, q, ind);

		while (ind + 1 < n) {
			elem_exp(B,B,q);
			mul(M,M,B);
			//cout << M << endl;
			ind++;
		}

		if (M[1][0] != 0) {
			ZZ_pE temp;
			power(temp, M[0][0], q);
			rvals[i] = M[0][0] + temp;
			power(temp, M[1][0], q);
			rvals[i] += beta*temp;
		}
		else {
			rvals[i] = M[0][0];
		}


	}
	interpolate(interpol, ei, rvals);

	return interpol;

}

void gekeler(ZZ_pE g, ZZ_pE del, long p, long q_exp, int n, ZZ_pX P) {

	int cnum = n/2 + 1;
	mat_ZZ_pE f;
	mat_ZZ_pE A;
	vec_ZZ_pE alphas;
	ZZ_pE det;
	alphas.SetLength(cnum);
	int q = pow(p, q_exp);
	vec_ZZ_pE x;
	x.SetLength(cnum);
	q = p;
	ZZ_pX Pn = P;
	ZZ_pX cons = ZZ_pX(INIT_MONO, 1, 1);
	ZZ_pE mono = conv<ZZ_pE>(cons);

	ZZ_pX b = pow(-1,n) * norm(del) * P;
	ZZ_pE temp;
	ZZ_pE temp2;
	ZZ_pE norm_delta_inv = inv(conv<ZZ_pE>(norm(del)));
	power(temp, del, pow(q,2*n));

	f.SetDims(2*n+1, 2*n+1);
	A.SetDims(cnum, cnum );
	f[0][0] = 1;
	//f[0][1] =
	f[1][0] = mono;
	f[1][1] = g;
	f[1][2] = del;


	for (int i = 2; i <= n; i++) {
		for (int j = 0; j < 2*i - 1; j++) {
			if (j == 0 ) {
				power(f[i][j], mono, i);
			}
			else if (j == 1) {
				power(temp, f[i-1][j-1], q);
				f[i][j] = mono*f[i-1][j] + g*temp;
			}
			else {
				power(temp, f[i-1][j-1], q);
				f[i][j] = mono*f[i-1][j] + g*temp;
				power(temp, f[i-1][j-2], q*q);
				f[i][j] += del*temp;
			}

		}
		power(temp, f[i-1][2*i-2], q);
		f[i][2*i-1] += g*temp;
		power(temp, f[i-1][2*i-3], q*q);
		f[i][2*i-1] += del*temp;


		power(temp, f[i-1][2*i-2], q*q);
		f[i][2*i] += del*temp;

	}

	cout << "f: " << f << endl;

	if (n % 2 == 0) {
		alphas[cnum-1] = 1 + pow(-1,n)*norm_delta_inv*Pn[n]*f[n][2*n];
	}

	int end;
	if (n % 2 == 1) {
		end = cnum;
	}
	else end = cnum - 1;

	for (int j = 0; j < end; j++) {
		for(int i = ((2*j+n)/2); i <= n; i++  ) {
			alphas[j] += Pn[i]*f[i][2*j+n];
			cout << "i: " << i << " p: " << Pn[i] << " f: " << f[i][2*j+n] << endl;
			cout << "j: " << j << " :: " << alphas[j] << endl;
		}
		alphas[j] *= pow(-1,n)*norm_delta_inv;

	}

	for (int i = 0; i < cnum; i++) {
		for (int j = 0; j < cnum; j++) {
			cout << f[i][2*j];
			A[j][i] = f[i][2*j];
		}
	}

	cout << "matrix: " << A << endl;
	cout << "alphas: " << alphas << endl;

	vec_ZZ_pE xx;
	xx.SetLength(cnum);


	solve(det, A, xx, alphas);

	cout << "gekeler " << xx << endl;

}


int main() {
	int p = 1163, q_exp = 1, n = 50;


	ZZ_p::init(ZZ(p));
	ZZ_pX P;
	BuildIrred(P, n);
	ZZ_pE::init(P);

	ZZ_pE g, del;
	set(g);
	set(del);
	double start1, start2, start3, end1, end2, end3;

	ZZ_pE base, res, op, elem;
	ZZ_pX cons = ZZ_pX(INIT_MONO, 1, 1);
	ZZ_pE mono = conv<ZZ_pE>(cons);
	elem = mono;
	elem = elem*mono + elem + 1;
	op = mono*mono + mono + 1;
	power(base, elem, p*p);
	frob(res, elem, p, 2);
	//cout << "ref: " << base << "res: " << res << endl;

	ZZ_pX out1 = char_poly1(g,del,p,q_exp,n,P);


	ZZ_pEX out2 = char_poly2(g,del,p,q_exp,n,P);
	//cout << "P: " << P << endl;

	//cout << "Char poly (det): " << out2 << " char poly (rand): " << out1 << endl;
}
