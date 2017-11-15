/*
 * cexp.cpp
 *
 *  Created on: Jul 18, 2013
 *      Author: terhorst
 */


#include <iostream>
#include <cstdio>
#include <unordered_map>
#include <vector>
#include <google/sparse_hash_map>

#include "cexp.h"
#include "hash.h"
#include "rhs.h"
#include "lhs.h"
#include "transition.h"
#include "mgf.h"

#include "print.h"

using google::sparse_hash_map; 

static int cond_time;
static ExpVec cond_exp;
static long cacheHits = 0;
static long cacheMisses = 0;

std::unordered_map<std::pair<int, Params>, Freq> freqMemo;
Freq deterministic(int t, const Params &p) {
    std::pair<int, Params> key{t, p};
    Params p2 = p;
    if (!freqMemo.count(key)) {
        Freq ret;
        if (t == 0) {
            ret = p.init;
        } else {
            p2.init = deterministic(t - 1, p);
            ret = transition(p2);
        }
        freqMemo[key] = ret;
    }
    return freqMemo[key];
}

double Edelta_lhs(int t, const Params &p, const ExpVec &e) {
    return Edelta(t, p, (t == cond_time) ? e + cond_exp : e);
}

// std::unordered_map<EdeltaOpts, double> EdeltaMemo;
sparse_hash_map<EdeltaOpts, double, std::hash<EdeltaOpts>> EdeltaMemo;
void initHashMap() {
    return;
    EdeltaOpts empty;
    empty.t = -9999;
    // EdeltaMemo.set_empty_key(empty);
}

double Edelta(int t, const Params &p, const ExpVec &e) {
    int ct = 0;
    ExpVec ce = Exponents::ZERO;
    if (cond_time < t) {
        ct = cond_time;
        ce = cond_exp;
    }
    EdeltaOpts key{t, p, e, ct, ce};
    if (!EdeltaMemo.count(key)) {
        cacheMisses++;
        // std::cout << t << " " << e << " " << ct << " " << ce << std::endl;
        double r, l;
        r = rhs(t, p, e);
        if (t == 1) {
            l = mgf(p, e);
        } else {
            l = lhs(t, p, e);
        }
        EdeltaMemo[key] = l - r;
        // std::cout << "t: " << t << " e:" << e << " l:" << l << " r:" << r << std::endl;
    } else {
        cacheHits++;
    }
    return EdeltaMemo[key];
}

double EZ(const Params &p, int t, const ExpVec &e) {
    return (t == 1) ? mgf(p, e) : lhs(t, p, e);
}

double Expectation(const Params &p, int t, Position pos) {
    cond_time = 0;
    cond_exp = Exponents::ZERO;
    if (t == 1) {
        return mgf(p, Exponents::Z(pos));
    } else {
        return lhs(t, p, Exponents::Z(pos));
    }
}

double Covariance(const Params &p, int t1, const Position pos1, int t2, const Position pos2) {
    if (t1 < t2) {
        return Covariance(p, t2, pos2, t1, pos1);
    }
    ExpVec e1 = Exponents::Z(pos1);
    ExpVec e2 = Exponents::Z(pos2);
    cond_time = 0;
    cond_exp = Exponents::ZERO;
    double ed1 = Edelta(t1, p, e1);
    double ed2 = Edelta(t2, p, e2);
    double ret;
    if (t1 == t2) {
        ret = Edelta(t1, p, e1 + e2);
    } else {
        cond_time = t2;
        cond_exp = e2;
        ret = Edelta(t1, p, e1);
    }
    cond_time = 0;
    cond_exp = Exponents::ZERO;
    return ret - ed1 * ed2;
}

double cov(const Params &p, int t1, const Marginal m1, int t2, const Marginal m2) {
    double ret = 0.0;
    for (Position pos1 : Exponents::marginals(m1)) {
        for (Position pos2 : Exponents::marginals(m2)) {
            ret += Covariance(p, t1, pos1, t2, pos2);
        }
    }
    return ret;
}

double E(const Params &p, int t, const Marginal m) {
    double ret = 0.0;
    for (Position pos : Exponents::marginals(m)) {
        ret += Expectation(p, t, pos);
    }
    return ret;
}

void resetMemo() {
    cacheHits = 0;
    cacheMisses = 0;
    EdeltaMemo.clear();
    freqMemo.clear();
}

void printCacheStats() {
    std::cout << "Cache hits: " << cacheHits << "\tmisses :" << cacheMisses 
        << "\t hit rate: " << (double)cacheHits / (cacheHits + cacheMisses) << std::endl;
}

// Helper functions to expose to Cython
namespace helpers {
    double varL(const Params &p, int t) {
        return cov(p, t, Marginal::L, t, Marginal::L);
    }

    double varS(const Params &p, int t) {
        return cov(p, t, Marginal::S, t, Marginal::S);
    }

    double varR(const Params &p, int t) {
        return cov(p, t, Marginal::R, t, Marginal::R);
    }

    double covLL(const Params &p, int t1, int t2) {
        return cov(p, t1, Marginal::L, t2, Marginal::L);
    }

    double covSS(const Params &p, int t1, int t2) {
        return cov(p, t1, Marginal::S, t2, Marginal::S);
    }

    double covRR(const Params &p, int t1, int t2) {
        return cov(p, t1, Marginal::R, t2, Marginal::R);
    }

    double covLR(const Params &p, int tL, int tR) {
        return cov(p, tL, Marginal::L, tR, Marginal::R);
    }

    double covLS(const Params &p, int tL, int tS) {
        return cov(p, tL, Marginal::L, tS, Marginal::S);
    }

    double covRS(const Params &p, int tR, int tS) {
        return cov(p, tR, Marginal::R, tS, Marginal::S);
    }

    double EL(const Params &p, int t) {
        return E(p, t, Marginal::L);
    }

    double ES(const Params &p, int t) {
        return E(p, t, Marginal::S);
    }
    double ER(const Params &p, int t) {
        return E(p, t, Marginal::R);
    }
}
