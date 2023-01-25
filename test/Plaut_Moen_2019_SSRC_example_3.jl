using ThinWalledColumn

L = 2438.0
z = range(0.0, L, 12)
A = fill(6650.0, length(z))
Ix = fill(2.12E8, length(z))
Iy = fill(6.37E6, length(z))
Io = fill(2.19E8, length(z))
J = fill(2.11E5, length(z))
Cw = fill(3.06E11, length(z))
E = fill(200.0, length(z))
G = E ./ 2.6
xo = fill(0.0, length(z))
yo = fill(0.0, length(z))
hx = fill(0.0, length(z))
hy = fill(225.0, length(z))
kx = fill(0.0, length(z))
ky = fill(0.0, length(z))
kϕ = fill(0.0, length(z))

a1 = L/1000
a2 = L/1000
a3 = 0.00766

uo_zz = a1 .* -π^2/L^2 .* sin.(π*z/L)
vo_zz = a2 .* -π^2/L^2 .* sin.(π*z/L)
ϕo_zz = a3 .* -π^2/L^2 .* sin.(π*z/L)

#u''=v''=ϕ''=0 (simply-supported), u'=v'=ϕ'=0  (fixed), u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free, e.g., a cantilever)
end_boundary_conditions = ["simply-supported", "simply-supported"]

#supports
#z, u, v, ϕ
supports = [(0.0, "fixed", "fixed", "fixed"),
            (L, "fixed", "fixed", "fixed")]


P = fill(750.0, length(z))

model = ThinWalledColumn.solve(z, A, Ix, Iy, Io, J, Cw, E, G, hx, hy, xo, yo, kx, ky, kϕ, uo_zz, vo_zz, ϕo_zz, P, supports, end_boundary_conditions);

isapprox(maximum(model.solution.ϕ .+ a3), 0.0097, rtol=0.01)   #Figure 17

isapprox(maximum(model.solution.u .+ a1), 3.78, rtol=0.01)   #Figure 18