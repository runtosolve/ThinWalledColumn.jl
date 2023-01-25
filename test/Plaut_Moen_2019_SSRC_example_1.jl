using ThinWalledColumn

L = 2438.0
z = range(0.0, L, 12)
A = fill(272.0, length(z))
Ix = fill(363370.0, length(z))
Iy = fill(64100.0, length(z))
Io = fill(716363.0, length(z))
J = fill(188.0, length(z))
Cw = fill(122720891.0, length(z))
E = fill(200.0, length(z))
G = E ./ 2.6
xo = fill(-32.59, length(z))
yo = fill(0.0, length(z))
hx = fill(0.0, length(z))
hy = fill(0.0, length(z))
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


P = fill(5.0, length(z))

model = ThinWalledColumn.solve(z, A, Ix, Iy, Io, J, Cw, E, G, hx, hy, xo, yo, kx, ky, kϕ, uo_zz, vo_zz, ϕo_zz, P, supports, end_boundary_conditions);

isapprox(maximum(model.solution.ϕ .+ a3), 0.02, rtol=0.01)   #Figure 4


P = fill(15.0, length(z))

model = ThinWalledColumn.solve(z, A, Ix, Iy, Io, J, Cw, E, G, hx, hy, xo, yo, kx, ky, kϕ, uo_zz, vo_zz, ϕo_zz, P, supports, end_boundary_conditions);

isapprox(maximum(model.solution.v .+ a2), 3.4, rtol=0.01)   #Figure 5