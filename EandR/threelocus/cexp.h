/*
 * cexp.h
 *
 *  Created on: Jul 18, 2013
 *      Author: terhorst
 */

#ifndef CEXP_H
#define CEXP_H

#include <vector>
#include <array>
#include <map>

typedef std::array<int, 7> ExpVec;
inline ExpVec operator+(const ExpVec &e1, const ExpVec &e2) {
    ExpVec r;
    for (int i = 0; i < 7; i++) {
        r[i] = e1[i] + e2[i]; 
    }
    return r;
}

enum class Position {L, S, R, LS, LR, RS, LRS};
enum class Marginal {L, S, R};

namespace Exponents {
    const ExpVec ZERO = {0, 0, 0, 0, 0, 0, 0};
    inline ExpVec Z(const Position &p) {
	switch (p) {
	    case Position::L:
 	        return {1, 0, 0, 0, 0, 0, 0};
            case Position::S:
                return {0, 1, 0, 0, 0, 0, 0};
	    case Position::R:
                return {0, 0, 1, 0, 0, 0, 0};
	    case Position::LS:
                return {0, 0, 0, 1, 0, 0, 0};
	    case Position::LR:
                return {0, 0, 0, 0, 1, 0, 0};
	    case Position::RS:
	        return {0, 0, 0, 0, 0, 1, 0};
	    case Position::LRS:
                return {0, 0, 0, 0, 0, 0, 1};
        }
    }
    inline std::array<Position, 4> marginals(const Marginal &m) {
    	switch (m) {
	    case Marginal::L:
	        return {Position::L, Position::LS, Position::LR, Position::LRS};
	    case Marginal::S:
                return {Position::S, Position::LS, Position::RS, Position::LRS};
	    case Marginal::R:
		return {Position::R, Position::LR, Position::RS, Position::LRS};
	}
    }
}

struct Freq {
    double zL, zS, zR, zLS, zLR, zRS, zLRS;
    bool operator==(const Freq& f) const {
        return zL==f.zL && zS==f.zS && zR==f.zR && 
            zLS==f.zLS && zLR==f.zLR && zRS==f.zRS && 
            zLRS==f.zLRS;
    }
};

struct Params {
    Freq init;
    int n;
    double h, s;
    double rLR, rLS, rRS, rLS_R, rRS_L, rLR_S;
    bool operator==(const Params& p) const {
        return init==p.init && n==p.n && h==p.h && s==p.s && rLR==p.rLR &&
               rLS==p.rLS && rRS==p.rRS && rLS_R==p.rLS_R && rRS_L==p.rRS_L && 
               rLR_S==p.rLR_S;
    }
};

struct EdeltaOpts {
    int t;
    Params p;
    ExpVec e;
    int ct;
    ExpVec ce;
    bool operator==(const EdeltaOpts& ed) const {
        return t==ed.t && p==ed.p && e==ed.e && ct==ed.ct && ce==ed.ce;
    }
};

Freq deterministic(int, const Params&);
double Edelta(int, const Params&, const ExpVec&);
double Edelta_lhs(int, const Params&, const ExpVec&);
double cov(const Params&, int, Marginal, int, Marginal);
double E(const Params&, int, Marginal);
void resetMemo();
void initHashMap();
void printCacheStats();

inline double Edelta(int t, const Params &p, 
    int eL, int eS, int eR, int eLS, int eLR, int eRS, int eLRS) {
    ExpVec e = {eL, eS, eR, eLS, eLR, eRS, eLRS};
    return Edelta(t, p, e);
}

inline double Edelta_lhs(int t, const Params &p, 
    int eL, int eS, int eR, int eLS, int eLR, int eRS, int eLRS) {
    ExpVec e = {eL, eS, eR, eLS, eLR, eRS, eLRS};
    return Edelta_lhs(t, p, e);
}

// Helper functions to expose to Cython
namespace helpers {
    double varL(const Params &p, int t);
    double varS(const Params &p, int t);
    double varR(const Params &p, int t);
    double covLL(const Params &p, int t1, int t2);
    double covSS(const Params &p, int t1, int t2);
    double covRR(const Params &p, int t1, int t2);
    double covLR(const Params &p, int tL, int tR);
    double covLS(const Params &p, int tL, int tS);
    double covRS(const Params &p, int tR, int tS);
    double EL(const Params &p, int t);
    double ES(const Params &p, int t);
    double ER(const Params &p, int t);
}
#endif
