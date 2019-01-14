

#include <iostream>
#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZX.h>

#include "BaseChange.h"
#include "MultipointEval.h"
#include "Util.h"
#include "HasseLift.h"
#include "HasseLiftExt.h"
#include "MultiComposeMod.h"
#include "FrobComp.h"

using namespace std;

void testHasseLift() {
    ZZ_pX f, g, delta, lift1, lift2;

    long degree = 10000;

    // build a square-free modulus
    Util util;
    util.randomMonic(f, degree);
    diff(g, f);
    GCD(g, g, f);
    div(f, f, g);
    ZZ_pXModulus F;
    build(F, f);
    degree = deg(F);

    random(g, degree);
    random(delta, degree);

    HasseLiftExt hasseLiftExt(g, delta, F);
    long start = util.getTimeMillis();
    hasseLiftExt.computeNaive(lift1, degree / 2);
    cout << util.getTimeMillis() - start << endl;

    start = util.getTimeMillis();
    hasseLiftExt.compute(lift2, degree / 2, 1);
    cout << util.getTimeMillis() - start << endl;

    if (lift1 == lift2)
        cout << "OK" << endl;
    else
        cout << "Failed" << endl;
}

int main(int argc, char** argv) {

    ZZ p;
    NextPrime(p, to_ZZ(7));
    ZZ_p::init(p);
    
    testHasseLift();
    
    return 0;
}

