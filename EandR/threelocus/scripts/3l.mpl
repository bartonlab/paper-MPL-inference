with(combinat, cartprod): ss := [`$`(0 .. ORDER)]: 
cp := cartprod([ss, ss, ss, ss, ss, ss, ss]):
ret := []:
while not cp[finished] do 
    ret := [op(ret), cp[nextvalue]()] 
end do:
vals := select(proc (x) options operator, arrow; convert(x, `+`) <= ORDER; end proc, ret):
valsm1 := vals:
ret := 'ret';
# read "../vals.mpl":

DLIST := [dL, dS, dR, dLS, dLR, dRS, dLRS]:

xition := solve(
    [pLp=zLp + zLSp + zLRp + zLRSp,
     pSp=zSp + zLSp + zRSp + zLRSp,
     pRp=zRp + zLRp + zRSp + zLRSp,
     cLSp=(zLSp + zLRSp) - pLp*pSp,
     cLRp=(zLRp + zLRSp) - pLp*pRp,
     cRSp=(zRSp + zLRSp) - pSp*pRp,
     cLRSp=zLRSp - pLp*cRSp - pSp*cLRp - pRp*cLSp - pLp*pSp*pRp],
     [zLp, zSp, zRp, zLSp, zLRp, zRSp, zLRSp]
)[1]:


pL := zL + zLS + zLR + zLRS;
pS := zS + zLS + zRS + zLRS;
pR := zR + zLR + zRS + zLRS;
cLS := (zLS + zLRS) - pL*pS:
cLR := (zLR + zLRS) - pL*pR:
cRS := (zRS + zLRS) - pS*pR:
cLRS := zLRS - pL*cRS - pS*cLR - pR*cLS - pL*pS*pR:

a_S_0 := s * (h + (1 - 2 * h) * pS) / (1 + s * pS * (2 * h + (1 - 2 * h) * pS)):
a_S_S := (1 - 2 * h) * s / (1 + s * pS * (2 * h + (1 - 2 * h) * pS)):
qS := 1 - pS:
qL := 1 - pL:
qR := 1 - pR:

rLRS := rLS_R + rRS_L + rLR_S:

# Deltas
delta_pS := a_S_0 * pS * qS:
delta_pL := a_S_0 * cLS:
delta_pR := a_S_0 * cRS:

g := (r) -> -r + a_S_0 * (1 - r) * (1 - 2 * pS) + a_S_S * r * pS * qS;

deltat_cLS := g(rLS) * cLS:
deltat_cRS := g(rRS) * cRS:
deltat_cLR := -rLR * cLR + a_S_0 * (1 - rLR) * cLRS + a_S_S * rLR * cLS * cRS:
deltat_cLRS := (-rLRS + a_S_0 * (1 - rLRS) * (1 - 2 * pS) + a_S_S * rLR_S * pS * qS) * cLRS -
    a_S_0 * (rLS_R + rRS_L) * pS * qS * cLR -
    (a_S_0 * (2 - rLS_R - rRS_L) - a_S_S * (1 - 2 * pS) * (rLS_R + rRS_L)) * cLS * cRS:

delta_cLS := deltat_cLS - a_S_0^2 * pS * qS * cLS:
delta_cRS := deltat_cRS - a_S_0^2 * pS * qS * cRS:
delta_cLR := deltat_cLR - a_S_0^2 * cLS * cRS:
delta_cLRS := deltat_cLRS - 
    a_S_0 * (pS * qS * deltat_cLR + cLS * deltat_cRS + cRS * deltat_cLS) +
    2 * a_S_0^3 * pS * qS * cLS * cRS:

B := eval(xition, [pLp=pL + delta_pL, pSp=pS + delta_pS, pRp=pR + delta_pR, 
     cLSp=cLS + delta_cLS, cLRp=cLR + delta_cLR, cRSp=cRS + delta_cRS, 
     cLRSp=cLRS + delta_cLRS]
):

if ARG0 + ARG1 + ARG2 + ARG3 + ARG4 + ARG5 + ARG6 = 0 then
CodeGeneration:-C(B, optimize=true, coercetypes=false, deducetypes=false, defaulttype=float, output="_transition.cii"):
end if:

mgff := (p1*exp(t1)+p2*exp(t2)+p3*exp(t3)+p4*exp(t4)+p5*exp(t5)+p6*exp(t6)+p7*exp(t7)+(1-p1-p2-p3-p4-p5-p6-p7))^n:
multimgf := proc (i1, i2, i3, i4, i5, i6, i7) options operator, arrow; 
    s := i1 + i2 + i3 + i4 + i5 + i6 + i7:
    if s = 0 then 1 
    else 
        n^(-s)*simplify(subs(t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0, 
                             diff(mgff, `$`(t1, i1), `$`(t2, i2), `$`(t3, i3),
                                        `$`(t4, i4), `$`(t5, i5), `$`(t6, i6), 
                                        `$`(t7, i7))));
    end if;
end proc:

mgfsubs := proc (iL, iS, iR, iLS, iLR, iRS, iLRS) options operator, arrow;
    eval(multimgf(iL, iS, iR, iLS, iLR, iRS, iLRS),
         [p1=zLp, p2=zSp, p3=zRp, p4=zLSp, p5=zLRp, p6=zRSp, p7=zLRSp]);
    eval(%, B):
end proc:

deltasubs := proc (iL, iS, iR, iLS, iLR, iRS, iLRS)
    ms := mgfsubs(iL, iS, iR, iLS, iLR, iRS, iLRS):
    eval(%, [zL=zL + dL, zS=zS + dS, zR=zR + dR,
             zLS=zLS + dLS, zLR=zLR + dLR, zRS=zRS + dRS,
             zLRS=zLRS + dLRS]):
end proc:


llhs := proc (aL, aS, aR, aLS, aLR, aRS, aLRS)
rary := []:
ms := deltasubs(aL, aS, aR, aLS, aLR, aRS, aLRS):
ind := 1:
for vv in valsm1 do 
   print(vv):
   coeftayl(ms, DLIST = [0, 0, 0, 0, 0, 0, 0], vv):
   cg := eval(%, [zL=dtL, zS=dtS, zR=dtR, zLS=dtLS, zLR=dtLR, zRS=dtRS, zLRS=dtLRS]):
   rary := [op(rary), ret[ind] = cg * Edelta_lhs(t - 1, p, vv[1], vv[2], vv[3], vv[4], vv[5], vv[6], vv[7])]:
   ind := ind + 1:
end do:
return rary:
end proc:

rrhs := proc(aL, aS, aR, aLS, aLR, aRS, aLRS)
    local e, c, xx; 
    xx := expand(
        (dtL + dL)^aL *
        (dtS + dS)^aS *
        (dtR + dR)^aR *
        (dtLS + dLS)^aLS *
        (dtLR + dLR)^aLR *
        (dtRS + dRS)^aRS *
        (dtLRS + dLRS)^aLRS
        - dL^aL * dS^aS * dR^aR * dLS^aLS * dLR^aLR * dRS^aRS * dLRS^aLRS);
    c := coeffs(xx, DLIST, 'e'); 
    convert(zip(proc (ee, cc) 
        options operator, arrow; 
        degs := map((term) -> degree(ee, term), DLIST):
        Edelta(t, p, degs[1], degs[2], degs[3], degs[4], degs[5], degs[6], degs[7])*cc;
    end proc, [e], [c]), `+`);
end proc:

aL := ARG0;
aS := ARG1;
aR := ARG2;
aLS := ARG3;
aLR := ARG4;
aRS := ARG5;
aLRS := ARG6;
print(aL, aS, aR, aLS, aLR, aRS, aLRS);
ll := llhs(aL, aS, aR, aLS, aLR, aRS, aLRS):
outl := cat("_lhs_", aL, aS, aR, aLS, aLR, aRS, aLRS, ".cii"):
CodeGeneration:-C(ll, optimize=true, coercetypes=false, deducetypes=false, defaulttype=float, output=outl):

rr := rrhs(aL, aS, aR, aLS, aLR, aRS, aLRS):
outl := cat("_rhs_", aL, aS, aR, aLS, aLR, aRS, aLRS, ".cii"):
CodeGeneration:-C([ret = rr], optimize=true, coercetypes=false, deducetypes=false, defaulttype=float, output=outl):

mgfsubs(aL, aS, aR, aLS, aLR, aRS, aLRS):
mm := eval(%, [zL=dtL, zS=dtS, zR=dtR, zLS=dtLS, zLR=dtLR, zRS=dtRS, zLRS=dtLRS]):
outl := cat("_mgf_", aL, aS, aR, aLS, aLR, aRS, aLRS, ".cii"):
CodeGeneration:-C([ret = mm], optimize=true, coercetypes=false, deducetypes=false, defaulttype=float, output=outl):
