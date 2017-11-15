from __future__ import division
cimport numpy as np
cimport cython
import math

initHashMap()

cdef int* eta = [-1, 1, 1, -1]

cdef double EZi(double[:] pi, double r, double N, int t, int i):
    cdef double zi = pi[i - 1]
    cdef double ld = pi[0] * pi[3] - pi[1] * pi[2]
    return zi + eta[i - 1] * t * r * ld * (1.0 - (t - 1) / (2 * N))


cdef double ErZiLD(double[:] pi, double r, double N, int t, int i):
    cdef double zi = pi[i - 1]
    cdef double z1 = pi[0]
    cdef double z2 = pi[1]
    cdef double z3 = pi[2]
    cdef double z4 = pi[3]
    cdef double ld = pi[0] * pi[3] - pi[1] * pi[2]
    return (r / N) * (zi * ld * (N - 3 * t) + (t / 2) * ((1 - eta[i - 1]) * z1 * z4 - (1 + eta[i - 1]) * z2 * z3))

cdef double EZiZj(double[:] pi, double r, double N, int t, int i, int j):
    cdef double zi = pi[i - 1]
    cdef double zj = pi[j - 1]
    cdef int ei = eta[i - 1]
    cdef int ej = eta[j - 1]
    cdef double z1 = pi[0]
    cdef double z2 = pi[1]
    cdef double z3 = pi[2]
    cdef double z4 = pi[3]
    cdef double yij = -zi * zj if i != j else zi * (1.0 - zi)
    cdef double ld = pi[0] * pi[3] - pi[1] * pi[2]
    return (zi * zj + ei * ej * t * r * (ei * zi + ej * zj) * ld +
            (t / N) * (yij + r * ((t + 1 - abs(ei - ej)) / 2 * (z1 * z4 + z2 * z3) -
                ei * ej * (2 * t - 1) * ld * (ei * zi + ej * zj) - 
                abs(ei + ej) / 8 * ((ei + ej) * (t + 1) * ld * int(i != j) +
                    4 * t * ((ei + 1) * z2 * z3 + (1 - ei) * z1 * z4)))))

cdef double EZitiZjtj(double[:] pi, double r, double N, int ti, int i, int tj, int j):
    if ti < tj:
        return EZitiZjtj(pi, r, N, tj, j, ti, i)
    u = ti - tj
    return EZiZj(pi, r, N, tj, i, j) + eta[i - 1] * u * (1 - (u - 1) / (2 * N)) * ErZiLD(pi, r, N, tj, j)


cdef double covZiZj(double[:] pi, double r, double N, int ti, int i, int tj, int j):
    return EZitiZjtj(pi, r, N, ti, i, tj, j) - EZi(pi, r, N, ti, i) * EZi(pi, r, N, tj, j)

def covXY_neutral(double[:] pi, double r, double N, int tx, int ty):
    "Compute cov(X_(tx), Y_(ty)) under neutrality."
    # X = z1_tx + z2_tx
    # Y = z1_ty + z3_ty
    # Neutral equations were derived for a haploid model, so multiply
    # pop. size by 2
    return covZiZj(pi, r, 2 * N, tx, 1, ty, 1) + covZiZj(pi, r, 2 * N, tx, 2, ty, 1) + \
            covZiZj(pi, r, 2 * N, tx, 1, ty, 3) + covZiZj(pi, r, 2 * N, tx, 2, ty, 3)

def covXX_neutral(double[:] pi, double r, double N, int tx1, int tx2):
    "Compute cov(X_(tx1), X_(tx2)) under neutrality."
    # Neutral equations were derived for a haploid model, so multiply
    # pop. size by 2
    return covZiZj(pi, r, 2 * N, tx1, 1, tx2, 1) + covZiZj(pi, r, 2 * N, tx1, 2, tx2, 1) + \
            covZiZj(pi, r, 2 * N, tx1, 2, tx2, 1) + covZiZj(pi, r, 2 * N, tx1, 2, tx2, 2)

# @cython.boundscheck(False)
# @cython.cdivision(True)
def mean_cov_neutral(
        double[:, :] M,
        double[:, :, :, :] C, 
        int[:] positions,
        double[:, :, :] hap_freqs,
        int[:] times, 
        int N,
        double r):
    '''Calculate mean and covariance under a neutral model.

    Parameters:
    M -- Matrix of expected values, updated in-place (T x L)
    C -- Matrix of covariances, updated in-place (T x L x T x L)
    positions -- Vector of site positions, used for recombination 
                 distance calculation (L)
    hap_freqs -- Matrix of haplotype frequencies for each pair of
                 sites (L x L x 3)
    times -- Vector of times for each data point (T)
    N -- Effective population size.
    r -- Recombination rate (/ base-pair / generation).
    '''
    cdef int T = C.shape[0]
    cdef int L = C.shape[1]
    cdef int i_t1, i_t2, i_ell1, i_ell2, t1, t2
    cdef double[:] pi 

    for i_ell1 in range(L):
        M[:, i_ell1] = hap_freqs[i_ell1, i_ell1, 0] + hap_freqs[i_ell1, i_ell1, 1]
        for i_ell2 in range(i_ell1, L):
            pi = hap_freqs[i_ell1, i_ell2]
            rr = min(abs(positions[i_ell1] - positions[i_ell2]) * r, 0.5)
            for i_t1 in range(T):
                t1 = times[i_t1]
                for i_t2 in range(T):
                    t2 = times[i_t2]
                    if i_ell1 == i_ell2:
                        c = covXX_neutral(pi, rr, N, t1, t2)
                    else:
                        c = covXY_neutral(pi, rr, N, t1, t2)
                    C[i_t1, i_ell1, i_t2, i_ell2] = \
                            C[i_t2, i_ell2, i_t1, i_ell1] = c

def mean_cov_approx(
        double[:, :] M, 
        double[:, :, :, :] C, 
        int[:] positions,
        int i_sel,
        double[:, :, :] hap_freqs,
        int[:] times, 
        int n,
        double r,
        double h,
        double s,
        bint skip_selected_site,
        double alpha,
        double beta,
        bint debug):
    '''Calculate mean and covariance with selection.

    Parameters:
    M -- Matrix of expected values (T x L)
    C -- Matrix of covariances (T x L x T x L)
    positions -- Vector of site positions, used for recombination 
                 distance calculation (L)
    i_sel -- Index of selected site.
    hap_freqs -- Matrix of haplotype frequencies for each pair of
                 sites + selected site (L x L x 7)
    times -- Vector of times for each data point (T)
    N -- Effective population size.
    r -- Recombination rate (/ base-pair / generation).
    h -- Dominance coefficient.
    s -- Selection coefficient.
    skip_selected_site -- Do not calculate selected site.
    '''
    cdef int T = M.shape[0]
    cdef int L = M.shape[1]
    cdef Params p
    cdef double c, ex
    cdef int i_t1, i_t2, i_ell1, i_ell2, t1, t2
    cdef double pL, pR, pLR, corr

    p.n = n
    p.h = h
    p.s = s
    for i_ell1 in range(L):
        for i_ell2 in range(i_ell1, L):
            p.init.zL = hap_freqs[i_ell1, i_ell2, 0]
            p.init.zS = hap_freqs[i_ell1, i_ell2, 1]
            p.init.zR = hap_freqs[i_ell1, i_ell2, 2]
            p.init.zLS = hap_freqs[i_ell1, i_ell2, 3]
            p.init.zLR = hap_freqs[i_ell1, i_ell2, 4]
            p.init.zRS = hap_freqs[i_ell1, i_ell2, 5]
            p.init.zLRS = hap_freqs[i_ell1, i_ell2, 6]
            pL = p.init.zL + p.init.zLS + p.init.zLR + p.init.zLRS
            pR = p.init.zR + p.init.zRS + p.init.zLR + p.init.zLRS
            pLR = p.init.zLR + p.init.zLRS
            corr = (pLR - pL * pR) / math.sqrt(pL * (1. - pL) * pR * (1. - pR))
            p.rLS = min(abs(positions[i_ell1] - positions[i_sel]) * r, 0.5)
            p.rLR = min(abs(positions[i_ell1] - positions[i_ell2]) * r, 0.5)
            p.rRS = min(abs(positions[i_ell2] - positions[i_sel]) * r, 0.5)
            if i_ell1 <= i_ell2 <= i_sel:
                # LRS
                p.rLS_R = 0.0
                p.rLR_S = p.rRS
                p.rRS_L = p.rLR
            elif i_ell1 <= i_sel <= i_ell2:
                # LSR
                p.rLR_S = 0.0
                p.rLS_R = p.rRS
                p.rRS_L = p.rLS
            else:
                # SLR
                p.rRS_L = 0.0
                p.rLS_R = p.rLR
                p.rLR_S = p.rLS

            for i_t1 in range(T):
                t1 = times[i_t1]
                if i_ell1 == i_sel:
                    M[i_t1, i_ell1] = ES(p, t1)
                else:
                    M[i_t1, i_ell1] = EL(p, t1)
                for i_t2 in range(T):
                    t2 = times[i_t2]
                    c = corr * math.exp(-alpha * p.rLR) * math.exp(-beta * abs(t1 - t2))
                    if debug:
                        print(i_t1, i_ell1, i_t2, i_ell2)
                    C[i_t1, i_ell1, i_t2, i_ell2] = \
                            C[i_t2, i_ell2, i_t1, i_ell1] = c

# @cython.boundscheck(False)
# @cython.cdivision(True)
def mean_cov(
        double[:, :] M, 
        double[:, :, :, :] C, 
        int[:] positions,
        int i_sel,
        double[:, :, :] hap_freqs,
        int[:] times, 
        int n,
        double r,
        double h,
        double s,
        bint skip_selected_site):
    '''Calculate mean and covariance with selection.

    Parameters:
    M -- Matrix of expected values (T x L)
    C -- Matrix of covariances (T x L x T x L)
    positions -- Vector of site positions, used for recombination 
                 distance calculation (L)
    i_sel -- Index of selected site.
    hap_freqs -- Matrix of haplotype frequencies for each pair of
                 sites + selected site (L x L x 7)
    times -- Vector of times for each data point (T)
    N -- Effective population size.
    r -- Recombination rate (/ base-pair / generation).
    h -- Dominance coefficient.
    s -- Selection coefficient.
    skip_selected_site -- Do not calculate selected site.
    '''
    cdef int T = M.shape[0]
    cdef int L = M.shape[1]
    cdef Params p
    cdef double c, ex
    cdef int i_t1, i_t2, i_ell1, i_ell2, t1, t2

    p.n = n
    p.h = h
    p.s = s
    for i_ell1 in range(L):
        if i_ell1 == i_sel and skip_selected_site:
            continue
        for i_ell2 in range(i_ell1, L):
            if i_ell2 == i_sel and skip_selected_site:
                continue
            p.init.zL = hap_freqs[i_ell1, i_ell2, 0]
            p.init.zS = hap_freqs[i_ell1, i_ell2, 1]
            p.init.zR = hap_freqs[i_ell1, i_ell2, 2]
            p.init.zLS = hap_freqs[i_ell1, i_ell2, 3]
            p.init.zLR = hap_freqs[i_ell1, i_ell2, 4]
            p.init.zRS = hap_freqs[i_ell1, i_ell2, 5]
            p.init.zLRS = hap_freqs[i_ell1, i_ell2, 6]

            p.rLS = min(abs(positions[i_ell1] - positions[i_sel]) * r, 0.5)
            p.rLR = min(abs(positions[i_ell1] - positions[i_ell2]) * r, 0.5)
            p.rRS = min(abs(positions[i_ell2] - positions[i_sel]) * r, 0.5)
            if i_ell1 <= i_ell2 <= i_sel:
                # LRS
                p.rLS_R = 0.0
                p.rLR_S = p.rRS
                p.rRS_L = p.rLR
            elif i_ell1 <= i_sel <= i_ell2:
                # LSR
                p.rLR_S = 0.0
                p.rLS_R = p.rRS
                p.rRS_L = p.rLS
            else:
                # SLR
                p.rRS_L = 0.0
                p.rLS_R = p.rLR
                p.rLR_S = p.rLS

            for i_t1 in range(T):
                t1 = times[i_t1]
                if i_ell1 == i_sel:
                    M[i_t1, i_ell1] = ES(p, t1)
                else:
                    M[i_t1, i_ell1] = EL(p, t1)
                for i_t2 in range(T):
                    t2 = times[i_t2]
                    if i_ell1 == i_sel and i_ell2 == i_sel:
                        c = covSS(p, t1, t2)
                    elif i_ell1 == i_sel and i_ell2 != i_sel:
                        c = covRS(p, t1, t2)
                    elif i_ell1 != i_sel and i_ell2 == i_sel:
                        c = covLS(p, t1, t2)
                    elif i_ell1 == i_ell2:
                        c = covLL(p, t1, t2)
                    else:
                        c = covLR(p, t1, t2)
                    C[i_t1, i_ell1, i_t2, i_ell2] = \
                            C[i_t2, i_ell2, i_t1, i_ell1] = c

cdef inline Params _mkparams(obj):
    cdef Freq f
    cdef Params p
    n = 1.0 * obj.n
    f.zL = obj.zL
    f.zS = obj.zS
    f.zR = obj.zR
    f.zLR = obj.zLR
    f.zRS = obj.zRS
    f.zLS = obj.zLS
    f.zLRS = obj.zLRS
    p.init = f
    p.n = obj.n
    p.h = obj.h
    p.s = obj.s
    p.rLR = obj.rLR
    p.rRS = obj.rRS
    p.rLS_R = obj.rLS_R
    p.rRS_L = obj.rRS_L
    p.rLR_S = obj.rLR_S
    return p

def reset():
    resetMemo()

def cacheStats():
    printCacheStats()

class Expectation:
    def __init__(self, int n, double h, double s, 
            double zL, double zS, double zR, double zLS, double zLR, double zRS, double zLRS, 
            double rLR, double rLS, double rRS, 
            double rLS_R, double rRS_L, double rLR_S):
        '''int n, double h, double s, int zL, int zS, int zR, 
            int zLS, int zLR, int zRS, int zLRS, double rLR, double rLS, double rRS, 
            double rLS_R, double rRS_L, double rLR_S'''
        self.n = n
        self.h = h
        self.s = s
        self.zL = zL
        self.zS = zS
        self.zR = zR
        self.zLS = zLS
        self.zLR = zLR
        self.zRS = zRS
        self.zLRS = zLRS
        self.rLR = rLR
        self.rLS = rLS
        self.rRS = rRS
        self.rLS_R = rLS_R
        self.rRS_L = rRS_L
        self.rLR_S = rLR_S
        self._p = _mkparams(self)

    def EL(self, t):
        return EL(self._p, t)

    def ES(self, t):
        return ES(self._p, t)

    def ER(self, t):
        return ER(self._p, t)

    def varL(self, t):
        return varL(self._p, t)

    def varS(self, t):
        return varS(self._p, t)

    def varR(self, t):
        return varR(self._p, t)

    def covLL(self, t1, t2):
        return covLL(self._p, t1, t2)

    def covSS(self, t1, t2):
        return covSS(self._p, t1, t2)

    def covRR(self, t1, t2):
        return covRR(self._p, t1, t2)

    def covLR(self, tL, tR):
        return covLR(self._p, tL, tR)

    def covLS(self, tL, tS):
        return covLS(self._p, tL, tS)

    def covRS(self, tR, tS):
        return covRS(self._p, tR, tS)

    def deterministic(self, t):
        return deterministic(t, self._p)
