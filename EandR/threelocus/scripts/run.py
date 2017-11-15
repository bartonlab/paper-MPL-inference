#!/usr/bin/env python

import sys
import itertools as it

order = int(sys.argv[1])
iparts = (ijk for ijk in it.product(range(order + 1), repeat=7) if sum(ijk) <= order)
with open("run.sh", "w") as f:
    for ip in iparts:
        args = {('ARG%i' % i): j for i, j in enumerate(ip)}
        args['ORDER'] = order
        args = " ".join("-D '%s=%i'" % aa for aa in args.items())
        f.write("~/maple17/bin/maple %s ../3l.mpl\n" % args)
