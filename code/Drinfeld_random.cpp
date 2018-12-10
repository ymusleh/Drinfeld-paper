



// Generic function definitions

int ilog(double base, double targ) {
	return (int) (ilogb(targ)/ ilogb(base));
}


void skew_mul(vector<ZZ_pE> & ret, vector<ZZ_pE>& a, vector<ZZ_pE>& b, int q ) {
	for (int i = 0; i < ret.size(); i++) clear(ret[i]);
	ret.resize(b.size() + a.size());
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < b.size(); j++) {
			add(ret[i + j], ret[i + j], a[i] * power(b[j], pow(q,i)) );
		}
	}
}

void skew_add(vector<ZZ_pE> & ret, vector<ZZ_pE>& a, vector<ZZ_pE>& b ) {
	int n_deg = max(b.size(), a.size());
	ret.resize(n_deg);
	for (int i = 0; i < n_deg; i++) {
		if (i < a.size()) {
			if (i < b.size()) {
				add(ret[i], b[i], a[i]);
			}
			else ret[i] = a[i];
		}
		else ret[i] = b[i];
	}
	return;

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
	int m = n;
	int cnum = n/2 + 1;
	int rnum = n;
	int q = pow(p, q_exp);
	q = p;
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
	cout << poly << std::endl;

	return poly;
}
