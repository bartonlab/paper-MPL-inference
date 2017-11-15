cdef extern from "cexp.h":
    cdef struct Freq:
        double zL, zS, zR, zLS, zLR, zRS, zLRS

    cdef struct Params:
        Freq init
        int n
        double h, s
        double rLR, rLS, rRS, rLS_R, rRS_L, rLR_S

    cdef Freq deterministic(int, const Params&)
    cdef void resetMemo()
    cdef void printCacheStats()
    cdef void initHashMap()

cdef extern from "cexp.h" namespace "helpers":
    cdef double varL(const Params&, int)
    cdef double varS(const Params&, int)
    cdef double varR(const Params&, int)
    cdef double covLL(const Params&, int, int)
    cdef double covSS(const Params&, int, int)
    cdef double covRR(const Params&, int, int)
    cdef double covLR(const Params&, int, int)
    cdef double covLS(const Params&, int, int)
    cdef double covRS(const Params&, int, int)
    cdef double EL(const Params&, int)
    cdef double ES(const Params&, int)
    cdef double ER(const Params&, int)
