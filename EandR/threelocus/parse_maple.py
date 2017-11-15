#!/usr/bin/env python2.7

from contextlib import contextmanager
import itertools
import sys

ORDER = int(sys.argv[1])
EXPS = "L S R LS LR RS LRS".split()
ZVARS = ["z" + s for s in EXPS]

@contextmanager
def header_file(fname):
    fname += ".h"
    f = open(fname, "wt")
    hid = fname.upper().replace("_", "__").replace(".", "_")
    f.write("#ifndef %s\n#define %s\n\n" % (hid, hid))
    yield f
    f.write("\n#endif\n")


@contextmanager
def double_func(ostream, fname, args, **kwargs):
    rtype = kwargs.get("rtype", "double")
    def stmt_write(s, indent_level=1, semi=True):
        pre = "    " * indent_level
        post = ";" if semi else ""
        ostream.write(pre + s + post + "\n")
    ostream.write("%s %s(" % (rtype, fname))
    ostream.write(", ".join(args))
    ostream.write(") {\n")
    yield stmt_write
    ostream.write("}\n\n")


def parse_inc(finc):
    with open(finc, 'rt') as f:
        vv = [line.strip().split()[0] for line in f]
        rs = [v for v in vv if v[0] == 'r']
        ts = set([v for v in vv if v[0] == 't'])
    return sorted(ts), sorted(rs)

ords = [v for v in itertools.product(range(ORDER + 1), repeat=7) if sum(v) <= ORDER]

unpacked_params = ["int n = p.n"]
unpacked_params += ["double %s = p.%s" % (name, name) for name in 
        "s h rRS rLR rLS rLS_R rRS_L rLR_S".split()]

for fname, fargs in (
    ("lhs", ["int t", "const Params &p"]),
    ("rhs", ["int t", "const Params &p"]),
    ("mgf", ["const Params &p"])):
    with header_file(fname) as h:
        for v in ords:
            vstr = "".join(map(str, v))
            incfname = "_%s_%s.cii" % (fname, vstr)
            ts, rs = parse_inc(incfname)
            with double_func(h, fname + "_" + vstr, fargs) as stmt_write:
                for p in unpacked_params:
                    stmt_write(p)
                if fname == "mgf":
                    stmt_write("int t = 1")
                stmt_write("Freq f = deterministic(%s, p)" % ("t" if fname == "rhs" else "t - 1"))
                for ev in EXPS:
                    stmt_write("double dt%s = f.z%s" % (ev, ev))
                if ts:
                    stmt_write("double %s" % ", ".join(ts))
                lrs = len(rs)
                a = "[%i]" % len(rs) if lrs > 1 else ""
                stmt_write("double ret%s" % a)
                stmt_write('#include "include/%s"' % incfname, 0, False)
                stmt_write('return %s' % " + ".join(rs))
        # Create switching function
        with double_func(h, fname, fargs + ["const ExpVec &e"]) as stmt_write:
            for v in ords:
                vstr = "".join(map(str, v))
                stmt_write("if (" + " && ".join("e[%i] == %i" % (i, v[i]) for i in range(7)) + ")", 1, False)
                stmt_write("return %s_%s(%s)" % (fname, vstr, ", ".join(a[-1] for a in fargs)), 2)
            stmt_write("return 0.0")

with header_file("transition") as h:
    with double_func(h, "transition", ["const Params &p"], rtype="Freq") as stmt_write:
        for p in unpacked_params:
            stmt_write(p)
        ts, _ = parse_inc("_transition.cii")
        stmt_write("double %s" % ", ".join(ts))
        stmt_write("Freq f")
        for zv in ZVARS:
            stmt_write("double %s = p.init.%s" % (zv, zv))
            stmt_write("double %sp" % zv);
        stmt_write('#include "include/_transition.cii"', 0, False)
        for zv in ZVARS:
            stmt_write("f.%s = %sp" % (zv, zv))
        stmt_write("return f")
