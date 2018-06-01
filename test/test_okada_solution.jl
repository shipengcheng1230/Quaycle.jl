##
include("../src/dc3d.jl")

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

res = dc3d(x, y, z, α, dep, dip, al1, al2, aw1, aw2, disl1, disl2, disl3)
