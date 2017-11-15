#ifndef RHS_H
#define RHS_H

double rhs_0000000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double ret;
#include "include/_rhs_0000000.cii"
    return ret;
}

double rhs_0000001(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1;
    double ret;
#include "include/_rhs_0000001.cii"
    return ret;
}

double rhs_0000002(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t2, t4;
    double ret;
#include "include/_rhs_0000002.cii"
    return ret;
}

double rhs_0000010(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1;
    double ret;
#include "include/_rhs_0000010.cii"
    return ret;
}

double rhs_0000011(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0000011.cii"
    return ret;
}

double rhs_0000020(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t2, t4;
    double ret;
#include "include/_rhs_0000020.cii"
    return ret;
}

double rhs_0000100(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1;
    double ret;
#include "include/_rhs_0000100.cii"
    return ret;
}

double rhs_0000101(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0000101.cii"
    return ret;
}

double rhs_0000110(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0000110.cii"
    return ret;
}

double rhs_0000200(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t2, t4;
    double ret;
#include "include/_rhs_0000200.cii"
    return ret;
}

double rhs_0001000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1;
    double ret;
#include "include/_rhs_0001000.cii"
    return ret;
}

double rhs_0001001(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0001001.cii"
    return ret;
}

double rhs_0001010(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0001010.cii"
    return ret;
}

double rhs_0001100(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0001100.cii"
    return ret;
}

double rhs_0002000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t2, t4;
    double ret;
#include "include/_rhs_0002000.cii"
    return ret;
}

double rhs_0010000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1;
    double ret;
#include "include/_rhs_0010000.cii"
    return ret;
}

double rhs_0010001(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0010001.cii"
    return ret;
}

double rhs_0010010(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0010010.cii"
    return ret;
}

double rhs_0010100(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0010100.cii"
    return ret;
}

double rhs_0011000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0011000.cii"
    return ret;
}

double rhs_0020000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t2, t4;
    double ret;
#include "include/_rhs_0020000.cii"
    return ret;
}

double rhs_0100000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1;
    double ret;
#include "include/_rhs_0100000.cii"
    return ret;
}

double rhs_0100001(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0100001.cii"
    return ret;
}

double rhs_0100010(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0100010.cii"
    return ret;
}

double rhs_0100100(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0100100.cii"
    return ret;
}

double rhs_0101000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0101000.cii"
    return ret;
}

double rhs_0110000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_0110000.cii"
    return ret;
}

double rhs_0200000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t2, t4;
    double ret;
#include "include/_rhs_0200000.cii"
    return ret;
}

double rhs_1000000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1;
    double ret;
#include "include/_rhs_1000000.cii"
    return ret;
}

double rhs_1000001(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_1000001.cii"
    return ret;
}

double rhs_1000010(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_1000010.cii"
    return ret;
}

double rhs_1000100(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_1000100.cii"
    return ret;
}

double rhs_1001000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_1001000.cii"
    return ret;
}

double rhs_1010000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_1010000.cii"
    return ret;
}

double rhs_1100000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t4, t6;
    double ret;
#include "include/_rhs_1100000.cii"
    return ret;
}

double rhs_2000000(int t, const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    Freq f = deterministic(t, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t2, t4;
    double ret;
#include "include/_rhs_2000000.cii"
    return ret;
}

double rhs(int t, const Params &p, const ExpVec &e) {
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0000000(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return rhs_0000001(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 2)
        return rhs_0000002(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return rhs_0000010(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 1)
        return rhs_0000011(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 2 && e[6] == 0)
        return rhs_0000020(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return rhs_0000100(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 1)
        return rhs_0000101(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 1 && e[6] == 0)
        return rhs_0000110(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 2 && e[5] == 0 && e[6] == 0)
        return rhs_0000200(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0001000(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return rhs_0001001(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return rhs_0001010(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return rhs_0001100(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 2 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0002000(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0010000(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return rhs_0010001(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return rhs_0010010(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return rhs_0010100(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0011000(t, p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 2 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0020000(t, p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0100000(t, p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return rhs_0100001(t, p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return rhs_0100010(t, p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return rhs_0100100(t, p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0101000(t, p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0110000(t, p);
    if (e[0] == 0 && e[1] == 2 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_0200000(t, p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_1000000(t, p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return rhs_1000001(t, p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return rhs_1000010(t, p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return rhs_1000100(t, p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_1001000(t, p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_1010000(t, p);
    if (e[0] == 1 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_1100000(t, p);
    if (e[0] == 2 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return rhs_2000000(t, p);
    return 0.0;
}


#endif
