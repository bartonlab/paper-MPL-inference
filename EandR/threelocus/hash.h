#ifndef HASH_H
#define HASH_H

#include <cmath>
#include <functional>

namespace std {

    template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

    template<typename T, size_t N> 
        struct hash<array<T, N>> : public unary_function<array<T, N>, size_t>
        {
            size_t operator()(const array<T, N> &a) const
            {
                size_t ret = 0;
                for (T i = 0; i < N; ++i)
                    hash_combine(ret, a[i]);
                return ret;
            }
        };

    template<class T1, class T2>
        struct hash<pair<T1, T2>> : public unary_function<pair<T1, T2>, size_t>
        {
            size_t operator()(pair<T1, T2> const &p) const
            {
                size_t ret = 0;
                hash_combine(ret, p.first);
                hash_combine(ret, p.second);
                return ret;
            }
        };

    template<>
        struct hash<Freq> : public unary_function<Freq, size_t>
        {
            size_t operator()(Freq const &f) const
            {
                size_t ret = 0;
                hash_combine(ret, f.zL);
                hash_combine(ret, f.zS);
                hash_combine(ret, f.zR);
                hash_combine(ret, f.zLS);
                hash_combine(ret, f.zLR);
                hash_combine(ret, f.zRS);
                hash_combine(ret, f.zLRS);
                return ret;
            }
        };

    template<>
        struct hash<Params> : public std::unary_function<Params, size_t>
        {
            size_t operator()(const Params &p) const
            {
                size_t ret = 0;
                hash_combine(ret, p.init);
                hash_combine(ret, p.n);
                hash_combine(ret, p.h);
                hash_combine(ret, p.s);
                hash_combine(ret, p.rLR);
                hash_combine(ret, p.rLS);
                hash_combine(ret, p.rLS_R);
                hash_combine(ret, p.rRS_L);
                hash_combine(ret, p.rLR_S);
                return ret;
            }
        };

    template<>
        struct hash<EdeltaOpts> : public unary_function<EdeltaOpts, size_t>
        {
            size_t operator()(const EdeltaOpts &ed) const
            {
                size_t ret = 0;
                hash_combine(ret, ed.t);
                hash_combine(ret, ed.p);
                hash_combine(ret, ed.e);
                hash_combine(ret, ed.ct);
                hash_combine(ret, ed.ce);
                return ret;
            }
        };


}

#endif
