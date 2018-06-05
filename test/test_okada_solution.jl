##
include("src/RateState.jl")
importall RateState

## parameters setting
α = 0.3
x = 15.
y = 30.
z = -20.
dep = 60.
dip = 45.
al1 = -50.
al2 = 50.
aw1 = 0.
aw2 = -30.
disl1 = 1.
disl2 = 1.
disl3 = 1.

res1 = dc3d_fortran(x, y, z, α, dep, dip, al1, al2, aw1, aw2, disl1, disl2, disl3)
res2 = dc3d_wrapper([x, y, z], α, dep, dip, [al1, al2], [aw1, aw2], [disl1, disl2, disl3])
res3 = dc3d_wrapper([x y z], α, dep, dip, [al1 al2], [aw1 aw2], [disl1 disl2 disl3])
