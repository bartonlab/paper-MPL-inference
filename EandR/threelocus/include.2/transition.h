#ifndef TRANSITION_H
#define TRANSITION_H

Freq transition(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    double t1, t10, t100, t101, t103, t105, t106, t107, t109, t11, t111, t12, t120, t125, t13, t135, t136, t137, t138, t139, t14, t140, t141, t16, t17, t18, t19, t2, t20, t21, t22, t25, t26, t27, t28, t29, t3, t30, t31, t33, t34, t35, t38, t4, t47, t5, t53, t6, t63, t68, t69, t70, t71, t72, t73, t76, t81, t86, t87, t90, t91, t92, t94, t95, t96, t97, t99;
    Freq f;
    double zL = p.init.zL;
    double zLp;
    double zS = p.init.zS;
    double zSp;
    double zR = p.init.zR;
    double zRp;
    double zLS = p.init.zLS;
    double zLSp;
    double zLR = p.init.zLR;
    double zLRp;
    double zRS = p.init.zRS;
    double zRSp;
    double zLRS = p.init.zLRS;
    double zLRSp;
#include "include/_transition.cii"
    f.zL = zLp;
    f.zS = zSp;
    f.zR = zRp;
    f.zLS = zLSp;
    f.zLR = zLRp;
    f.zRS = zRSp;
    f.zLRS = zLRSp;
    return f;
}


#endif
