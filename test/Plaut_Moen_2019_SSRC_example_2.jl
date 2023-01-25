using ThinWalledColumn

L = 6100.0
z = range(0.0, L, 12)
A = fill(6450.0, length(z))
Ix = fill(1.357E7, length(z))
Iy = fill(2.527E7, length(z))
Io = fill(4.17E7, length(z))
J = fill(624350.0, length(z))
Cw = fill(8.62E8, length(z))
E = fill(200.0, length(z))
G = E ./ 2.6
xo = fill(0.0, length(z))
yo = fill(-23.6, length(z))
hx = fill(0.0, length(z))
hy = fill(0.0, length(z))
kx = fill(0.0, length(z))
ky = fill(0.0, length(z))
kϕ = fill(0.0, length(z))

a1 = -L/1000   #tricky, not fully clear SSRC paper
a2 = L/1000
a3 = 0.0192

uo_zz = a1 .* -π^2/L^2 .* sin.(π*z/L)
vo_zz = a2 .* -π^2/L^2 .* sin.(π*z/L)
ϕo_zz = a3 .* -π^2/L^2 .* sin.(π*z/L)

#u''=v''=ϕ''=0 (simply-supported), u'=v'=ϕ'=0  (fixed), u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free, e.g., a cantilever)
end_boundary_conditions = ["simply-supported", "simply-supported"]

#supports
#z, u, v, ϕ
supports = [(0.0, "fixed", "fixed", "fixed"),
            (L, "fixed", "fixed", "fixed")]


P = fill(600.0, length(z))

model = ThinWalledColumn.solve(z, A, Ix, Iy, Io, J, Cw, E, G, hx, hy, xo, yo, kx, ky, kϕ, uo_zz, vo_zz, ϕo_zz, P, supports, end_boundary_conditions);

isapprox(maximum(model.solution.ϕ .+ a3), 0.0246, rtol=0.01)   #Figure 11


P = fill(1000.0, length(z))

model = ThinWalledColumn.solve(z, A, Ix, Iy, Io, J, Cw, E, G, hx, hy, xo, yo, kx, ky, kϕ, uo_zz, vo_zz, ϕo_zz, P, supports, end_boundary_conditions);

isapprox(minimum(model.solution.u .+ a1), -26.9, rtol=0.01)   #Figure 13