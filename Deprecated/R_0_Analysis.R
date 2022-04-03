# R0 For Human Specific Population ----------------------------------------

f1 = quote(betaHH*Ish*Sh)
f2 = quote(betaHH*Irh*Sh*(1-a))

vm1 = quote((r*Ish+uh*Ish))
vm2 = quote((r*Irh+uh*Irh))

vp1 = 0
vp2 = 0

v1 = substitute(a-b, list(a = vm1, b = vp1))
v2 = substitute(a-b, list(a = vm2, b = vp2))

f11 = D(f1, "Ish"); f12 = D(f1, "Irh")
f21 = D(f2, "Ish"); f22 = D(f2, "Irh")

v11 = D(v1, "Ish"); v12 = D(v1, "Irh")
v21 = D(v2, "Ish"); v22 = D(v2, "Irh")

PARAS = list(Sh = 0.98, Ish = 0, Irh = 0,
             a = 0.3493491, betaHH = 0.00001, r =5.5^-1, uh = 28835^-1)

f = with(PARAS,
         matrix(c(eval(f11), eval(f12), eval(f21), eval(f22)),
                nrow = 2, byrow = 2))

v = with(PARAS,
         matrix(c(eval(v11), eval(v12), eval(v21), eval(v22)),
                nrow = 2, byrow = 2))

eigen(f %*% solve(v))$values

# R0 For Whole System Specific Population ----------------------------------------

fsa = quote(betaAA*Isa*Sh)
fra = quote(betaHH*Irh*Sh*(1-a))

fsh = quote(betaHH*Ish*Sh)
frh = quote(betaHH*Irh*Sh*(1-a))

vmsa = quote((r*Ish+uh*Ish))
vmra = quote((r*Irh+uh*Irh))

vpsh = quote((r*Ish+uh*Ish))
vprh = quote((r*Ish+uh*Ish))



vp1 = 0
vp2 = 0

v1 = substitute(a-b, list(a = vm1, b = vp1))
v2 = substitute(a-b, list(a = vm2, b = vp2))


PARAS = list(Sh = 0.98, Ish = 0, Irh = 0,
             a = 0.3493491, betaHH = 0.0281, r =60^-1, uh = 240^-1)