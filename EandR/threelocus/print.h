#ifndef PRINT_H
#define PRINT_H

#include "prettyprint.hpp"

template<typename S>
S& operator<< ( S& stm, const Freq &f ) {
    return stm << "{ zL:" << f.zL << " zS:" << f.zS << " zR:" << f.zR <<
        " zLS:" << f.zLS << " zLR:" << f.zLR << " zRS:" << f.zRS <<
        " zLRS:" << f.zLRS << " }";
}


#endif
