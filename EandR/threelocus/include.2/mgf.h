#ifndef MGF_H
#define MGF_H

double mgf_0000000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double ret;
#include "include/_mgf_0000000.cii"
    return ret;
}

double mgf_0000001(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t101, t11, t12, t13, t14, t17, t19, t2, t22, t24, t25, t26, t29, t3, t37, t4, t40, t41, t45, t47, t48, t49, t5, t50, t51, t52, t53, t58, t59, t6, t60, t63, t64, t65, t66, t67, t69, t72, t85;
    double ret;
#include "include/_mgf_0000001.cii"
    return ret;
}

double mgf_0000002(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t103, t12, t13, t136, t137, t14, t15, t16, t19, t21, t24, t26, t27, t28, t3, t31, t39, t4, t42, t43, t47, t49, t5, t50, t51, t52, t53, t54, t55, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t8, t87;
    double ret;
#include "include/_mgf_0000002.cii"
    return ret;
}

double mgf_0000010(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t104, t105, t106, t107, t109, t11, t112, t12, t121, t125, t13, t14, t15, t16, t17, t18, t2, t22, t23, t3, t30, t33, t34, t36, t4, t40, t43, t5, t57, t58, t6, t62, t65, t7, t75, t8, t85, t9, t91, t93, t96;
    double ret;
#include "include/_mgf_0000010.cii"
    return ret;
}

double mgf_0000011(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t139, t14, t15, t16, t19, t21, t24, t26, t27, t28, t3, t31, t32, t39, t4, t42, t43, t47, t48, t49, t5, t50, t51, t52, t53, t54, t55, t57, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t79, t8, t87, t92;
    double ret;
#include "include/_mgf_0000011.cii"
    return ret;
}

double mgf_0000020(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t106, t107, t108, t109, t11, t111, t114, t12, t123, t127, t13, t137, t138, t14, t15, t16, t17, t18, t19, t20, t24, t25, t3, t32, t35, t36, t38, t4, t42, t45, t5, t59, t6, t60, t64, t67, t7, t77, t8, t87, t9, t93, t95, t98;
    double ret;
#include "include/_mgf_0000020.cii"
    return ret;
}

double mgf_0000100(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t102, t11, t12, t13, t14, t17, t19, t2, t22, t24, t25, t26, t29, t3, t38, t4, t41, t42, t46, t48, t49, t5, t50, t51, t52, t53, t54, t59, t6, t60, t61, t64, t65, t66, t67, t68, t70, t73, t76, t86;
    double ret;
#include "include/_mgf_0000100.cii"
    return ret;
}

double mgf_0000101(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t139, t14, t15, t16, t19, t21, t24, t26, t27, t28, t3, t31, t32, t39, t4, t42, t43, t47, t49, t5, t50, t51, t52, t53, t54, t55, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t77, t79, t8, t87, t92;
    double ret;
#include "include/_mgf_0000101.cii"
    return ret;
}

double mgf_0000110(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t105, t106, t107, t108, t109, t11, t111, t114, t119, t12, t123, t127, t129, t13, t132, t134, t136, t137, t14, t140, t15, t16, t17, t18, t19, t20, t24, t25, t3, t32, t35, t36, t38, t4, t42, t43, t45, t49, t5, t58, t59, t6, t60, t64, t67, t7, t77, t8, t87, t9, t92, t93, t95, t98;
    double ret;
#include "include/_mgf_0000110.cii"
    return ret;
}

double mgf_0000200(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t104, t12, t13, t137, t138, t14, t15, t16, t19, t21, t24, t26, t27, t28, t3, t31, t4, t40, t43, t44, t48, t5, t50, t51, t52, t53, t54, t55, t56, t6, t61, t62, t63, t66, t67, t68, t69, t7, t70, t72, t75, t78, t8, t88;
    double ret;
#include "include/_mgf_0000200.cii"
    return ret;
}

double mgf_0001000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t101, t11, t12, t13, t14, t17, t19, t2, t22, t23, t24, t25, t26, t29, t3, t37, t4, t40, t41, t45, t47, t48, t49, t5, t50, t51, t52, t53, t58, t59, t6, t60, t63, t64, t65, t66, t67, t69, t72, t85, t86, t88;
    double ret;
#include "include/_mgf_0001000.cii"
    return ret;
}

double mgf_0001001(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t138, t14, t15, t16, t19, t21, t24, t25, t26, t27, t28, t3, t31, t32, t39, t4, t42, t43, t47, t49, t5, t50, t51, t52, t53, t54, t55, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t79, t8, t87, t88, t90, t92;
    double ret;
#include "include/_mgf_0001001.cii"
    return ret;
}

double mgf_0001010(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t105, t106, t107, t108, t109, t11, t111, t114, t115, t117, t119, t12, t123, t127, t129, t13, t134, t135, t136, t137, t139, t14, t15, t16, t17, t18, t19, t20, t24, t25, t3, t32, t35, t36, t38, t4, t42, t43, t45, t49, t5, t58, t59, t6, t60, t64, t67, t7, t77, t8, t87, t9, t92, t93, t95, t98;
    double ret;
#include "include/_mgf_0001010.cii"
    return ret;
}

double mgf_0001100(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t102, t104, t108, t117, t12, t126, t13, t136, t137, t139, t14, t15, t16, t19, t21, t24, t25, t26, t27, t28, t3, t31, t32, t4, t40, t43, t44, t48, t5, t50, t51, t52, t53, t54, t55, t56, t6, t60, t61, t62, t63, t66, t67, t68, t69, t7, t70, t72, t75, t78, t8, t80, t88, t89, t91, t93;
    double ret;
#include "include/_mgf_0001100.cii"
    return ret;
}

double mgf_0002000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t103, t12, t13, t136, t137, t14, t15, t16, t19, t21, t24, t25, t26, t27, t28, t3, t31, t39, t4, t42, t43, t47, t49, t5, t50, t51, t52, t53, t54, t55, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t8, t87, t88, t90;
    double ret;
#include "include/_mgf_0002000.cii"
    return ret;
}

double mgf_0010000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t107, t11, t137, t14, t15, t16, t17, t18, t19, t2, t20, t21, t22, t23, t24, t25, t26, t28, t29, t3, t30, t31, t32, t33, t36, t4, t40, t45, t47, t5, t55, t57, t6, t60, t63, t64, t67, t68, t71, t82, t83, t84, t85, t87, t99;
    double ret;
#include "include/_mgf_0010000.cii"
    return ret;
}

double mgf_0010001(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t138, t14, t141, t142, t15, t16, t19, t21, t24, t26, t27, t28, t3, t30, t31, t32, t39, t4, t42, t43, t47, t48, t49, t5, t50, t51, t52, t53, t54, t55, t57, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t77, t79, t8, t87, t92;
    double ret;
#include "include/_mgf_0010001.cii"
    return ret;
}

double mgf_0010010(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t105, t106, t107, t108, t109, t11, t111, t113, t114, t119, t12, t123, t124, t127, t129, t13, t132, t134, t136, t137, t139, t14, t141, t142, t15, t16, t17, t18, t19, t20, t24, t25, t3, t32, t35, t36, t38, t4, t42, t43, t45, t49, t5, t58, t59, t6, t60, t64, t67, t7, t77, t8, t87, t9, t92, t93, t95, t98;
    double ret;
#include "include/_mgf_0010010.cii"
    return ret;
}

double mgf_0010100(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t102, t104, t108, t117, t12, t126, t13, t136, t137, t139, t14, t141, t142, t15, t16, t19, t21, t24, t26, t27, t28, t3, t30, t31, t32, t33, t4, t40, t43, t44, t48, t49, t5, t50, t51, t52, t53, t54, t55, t56, t58, t6, t60, t61, t62, t63, t66, t67, t68, t69, t7, t70, t72, t75, t78, t8, t80, t88, t93;
    double ret;
#include "include/_mgf_0010100.cii"
    return ret;
}

double mgf_0011000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t138, t14, t141, t142, t15, t16, t19, t21, t24, t25, t26, t27, t28, t3, t30, t31, t32, t39, t4, t42, t43, t47, t48, t49, t5, t50, t51, t52, t53, t54, t55, t57, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t77, t79, t8, t87, t88, t90, t92;
    double ret;
#include "include/_mgf_0011000.cii"
    return ret;
}

double mgf_0020000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t111, t113, t12, t13, t139, t14, t140, t141, t15, t16, t18, t19, t20, t21, t22, t23, t24, t27, t28, t29, t3, t30, t31, t32, t33, t35, t36, t37, t4, t40, t49, t5, t55, t6, t65, t7, t74, t75, t79, t8, t81, t82, t83, t84, t85, t88, t91, t92, t94;
    double ret;
#include "include/_mgf_0020000.cii"
    return ret;
}

double mgf_0100000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t100, t103, t11, t123, t124, t127, t13, t136, t15, t16, t17, t18, t19, t2, t23, t25, t28, t3, t32, t33, t34, t35, t36, t4, t43, t44, t47, t48, t49, t5, t50, t51, t53, t54, t55, t58, t6, t67, t73, t83, t89, t90, t91, t93, t95, t98, t99;
    double ret;
#include "include/_mgf_0100000.cii"
    return ret;
}

double mgf_0100001(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t138, t14, t140, t141, t15, t16, t19, t21, t23, t24, t25, t26, t27, t28, t3, t31, t32, t39, t4, t42, t43, t47, t48, t49, t5, t50, t51, t52, t53, t54, t55, t57, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t79, t8, t87, t88, t90, t92;
    double ret;
#include "include/_mgf_0100001.cii"
    return ret;
}

double mgf_0100010(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t105, t106, t107, t108, t109, t11, t111, t114, t115, t117, t119, t12, t122, t123, t124, t127, t129, t13, t134, t135, t136, t137, t139, t14, t140, t141, t15, t16, t17, t18, t19, t20, t24, t25, t3, t32, t35, t36, t38, t4, t42, t43, t45, t49, t5, t58, t59, t6, t60, t64, t67, t7, t77, t8, t87, t9, t92, t93, t95, t98;
    double ret;
#include "include/_mgf_0100010.cii"
    return ret;
}

double mgf_0100100(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t102, t104, t108, t117, t12, t126, t13, t136, t137, t139, t14, t141, t142, t15, t16, t19, t21, t23, t24, t25, t26, t27, t28, t3, t31, t32, t4, t40, t43, t44, t48, t49, t5, t50, t51, t52, t53, t54, t55, t56, t58, t6, t60, t61, t62, t63, t66, t67, t68, t69, t7, t70, t72, t75, t78, t8, t80, t88, t89, t91, t93;
    double ret;
#include "include/_mgf_0100100.cii"
    return ret;
}

double mgf_0101000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t138, t14, t140, t141, t15, t16, t19, t21, t23, t24, t25, t26, t27, t28, t3, t31, t32, t39, t4, t42, t43, t47, t48, t49, t5, t50, t51, t52, t53, t54, t55, t57, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t79, t8, t87, t88, t90, t92;
    double ret;
#include "include/_mgf_0101000.cii"
    return ret;
}

double mgf_0110000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t104, t109, t118, t119, t12, t121, t123, t128, t13, t138, t139, t14, t140, t142, t143, t144, t15, t16, t18, t21, t22, t23, t24, t25, t26, t27, t28, t29, t3, t30, t32, t33, t34, t35, t36, t37, t4, t40, t41, t42, t43, t44, t45, t48, t5, t57, t6, t63, t7, t73, t78, t79, t8, t81, t84, t85, t87, t90, t91, t93, t94, t96, t97, t98, t99;
    double ret;
#include "include/_mgf_0110000.cii"
    return ret;
}

double mgf_0200000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t100, t101, t102, t105, t12, t125, t126, t129, t13, t138, t139, t140, t15, t17, t18, t19, t20, t21, t25, t27, t3, t30, t34, t35, t36, t37, t38, t4, t45, t46, t49, t5, t50, t51, t52, t53, t55, t56, t57, t6, t60, t69, t7, t75, t8, t85, t91, t92, t93, t95, t97;
    double ret;
#include "include/_mgf_0200000.cii"
    return ret;
}

double mgf_1000000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t106, t108, t11, t111, t12, t13, t136, t14, t16, t17, t18, t19, t2, t20, t21, t22, t25, t26, t27, t28, t29, t3, t30, t31, t33, t34, t35, t38, t4, t47, t5, t53, t6, t63, t69, t70, t71, t72, t73, t76, t81, t87, t91, t92, t95, t99;
    double ret;
#include "include/_mgf_1000000.cii"
    return ret;
}

double mgf_1000001(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t138, t14, t140, t141, t15, t16, t18, t19, t21, t24, t25, t26, t27, t28, t3, t31, t32, t39, t4, t42, t43, t47, t49, t5, t50, t51, t52, t53, t54, t55, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t77, t79, t8, t87, t88, t90, t92;
    double ret;
#include "include/_mgf_1000001.cii"
    return ret;
}

double mgf_1000010(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t10, t105, t106, t107, t108, t109, t11, t111, t114, t115, t117, t119, t12, t123, t126, t127, t129, t13, t132, t134, t135, t136, t137, t14, t140, t141, t142, t15, t16, t17, t18, t19, t20, t24, t25, t3, t32, t35, t36, t38, t4, t42, t43, t45, t49, t5, t58, t59, t6, t60, t64, t67, t7, t77, t8, t87, t9, t92, t93, t95, t98;
    double ret;
#include "include/_mgf_1000010.cii"
    return ret;
}

double mgf_1000100(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t102, t104, t108, t117, t12, t126, t13, t136, t137, t139, t14, t140, t141, t15, t16, t18, t19, t21, t24, t25, t26, t27, t28, t3, t31, t32, t33, t4, t40, t43, t44, t48, t5, t50, t51, t52, t53, t54, t55, t56, t6, t60, t61, t62, t63, t66, t67, t68, t69, t7, t70, t72, t75, t78, t8, t80, t88, t89, t91, t93;
    double ret;
#include "include/_mgf_1000100.cii"
    return ret;
}

double mgf_1001000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t101, t103, t107, t116, t12, t125, t13, t135, t136, t139, t14, t140, t141, t15, t16, t18, t19, t21, t24, t25, t26, t27, t28, t3, t31, t32, t39, t4, t42, t43, t47, t49, t5, t50, t51, t52, t53, t54, t55, t59, t6, t60, t61, t62, t65, t66, t67, t68, t69, t7, t71, t74, t77, t79, t8, t87, t88, t90, t92;
    double ret;
#include "include/_mgf_1001000.cii"
    return ret;
}

double mgf_1010000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t100, t104, t105, t106, t107, t116, t117, t12, t126, t127, t129, t13, t131, t132, t134, t136, t137, t139, t140, t142, t143, t144, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t3, t30, t31, t32, t33, t34, t35, t38, t39, t4, t40, t41, t42, t43, t46, t5, t55, t6, t61, t7, t71, t76, t78, t79, t8, t81, t83, t87, t97, t99;
    double ret;
#include "include/_mgf_1010000.cii"
    return ret;
}

double mgf_1100000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t100, t101, t102, t105, t106, t115, t12, t124, t125, t126, t128, t129, t13, t131, t134, t136, t138, t139, t141, t143, t144, t15, t17, t18, t19, t20, t21, t24, t25, t27, t3, t30, t34, t35, t36, t37, t38, t4, t42, t45, t46, t49, t5, t50, t51, t52, t53, t55, t56, t57, t6, t60, t69, t7, t75, t8, t85, t90, t91, t92, t93, t95, t97, t99;
    double ret;
#include "include/_mgf_1100000.cii"
    return ret;
}

double mgf_2000000(const Params &p) {
    int n = p.n;
    double s = p.s;
    double h = p.h;
    double rRS = p.rRS;
    double rLR = p.rLR;
    double rLS = p.rLS;
    double rLS_R = p.rLS_R;
    double rRS_L = p.rRS_L;
    double rLR_S = p.rLR_S;
    int t = 1;
    Freq f = deterministic(t - 1, p);
    double dtL = f.zL;
    double dtS = f.zS;
    double dtR = f.zR;
    double dtLS = f.zLS;
    double dtLR = f.zLR;
    double dtRS = f.zRS;
    double dtLRS = f.zLRS;
    double t1, t111, t113, t12, t13, t138, t139, t14, t140, t15, t16, t18, t19, t20, t21, t22, t23, t24, t27, t28, t29, t3, t30, t31, t32, t33, t35, t36, t37, t4, t40, t49, t5, t55, t6, t65, t7, t71, t74, t79, t8, t80, t81, t82, t83, t84, t87, t91, t96, t97, t99;
    double ret;
#include "include/_mgf_2000000.cii"
    return ret;
}

double mgf(const Params &p, const ExpVec &e) {
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0000000(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return mgf_0000001(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 2)
        return mgf_0000002(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return mgf_0000010(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 1)
        return mgf_0000011(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 2 && e[6] == 0)
        return mgf_0000020(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return mgf_0000100(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 1)
        return mgf_0000101(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 1 && e[6] == 0)
        return mgf_0000110(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 2 && e[5] == 0 && e[6] == 0)
        return mgf_0000200(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0001000(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return mgf_0001001(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return mgf_0001010(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return mgf_0001100(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 0 && e[3] == 2 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0002000(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0010000(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return mgf_0010001(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return mgf_0010010(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return mgf_0010100(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 1 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0011000(p);
    if (e[0] == 0 && e[1] == 0 && e[2] == 2 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0020000(p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0100000(p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return mgf_0100001(p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return mgf_0100010(p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return mgf_0100100(p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0101000(p);
    if (e[0] == 0 && e[1] == 1 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0110000(p);
    if (e[0] == 0 && e[1] == 2 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_0200000(p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_1000000(p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 1)
        return mgf_1000001(p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 1 && e[6] == 0)
        return mgf_1000010(p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 1 && e[5] == 0 && e[6] == 0)
        return mgf_1000100(p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 0 && e[3] == 1 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_1001000(p);
    if (e[0] == 1 && e[1] == 0 && e[2] == 1 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_1010000(p);
    if (e[0] == 1 && e[1] == 1 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_1100000(p);
    if (e[0] == 2 && e[1] == 0 && e[2] == 0 && e[3] == 0 && e[4] == 0 && e[5] == 0 && e[6] == 0)
        return mgf_2000000(p);
    return 0.0;
}


#endif
