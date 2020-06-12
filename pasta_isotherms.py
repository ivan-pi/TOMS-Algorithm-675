import sympy as sp

aw, X, T = sp.symbols('a_{w} X_{db} T')

# Oswin
# a_\mathrm{w} =  \frac{(X_\mathrm{db}/a)^{1/b}}{1 + (X_\mathrm{db}/a)^{1/b}}
# 
# Henderson
# a_\mathrm{w} =  1 - \exp{\left(-a(T-b)^c {X_\mathrm{db}}^d\right)}
#
# Henderson 2
# a_\mathrm{w} =  1 - \exp{\left(-a T^{-b} {X_\mathrm{db}}^{cT^{-d}}\right)}
#
# GAB
# X_\mathrm{db} = \frac{X_\mathrm{m} C K a_\mathrm{w}}{\left(1 - K a_\mathrm{w}\right)\left(1 - K a_\mathrm{w} + C K a_\mathrm{w}\right)}

# Oswin relationship
print("\n -- Oswin --")
k0, k1, n0, n1 = sp.symbols('k_0 k_1 n_0 n_1')

oswin = (k0 + k1*T)*(aw/(1-aw))**(n0 + n1*T)
oswin_aw = sp.solve(oswin - X,aw)[0]
print(oswin_aw)

oswin_diff = sp.diff(oswin_aw,X)
print(oswin_diff)

# Modified Henderson Equation 1
print("\n -- Henderson 1 --")
A, B, C, D = sp.symbols('A B C D')
henderson1_aw = 1 - sp.exp(-A*(T-B)**C*(X*100)**D)
henderson1 = sp.solve(henderson1_aw, X)[0]
print(henderson1_aw)

henderson1_diff = sp.diff(henderson1_aw,X)
print(henderson1_diff)

# Modified Henderson Equation 2

henderson2_aw = 1 - sp.exp(-A*T**(-B)*X**(C*T**(-D)))
henderson2 = sp.solve(henderson2_aw - aw, X)[0]

# GAB equation
print("\n -- GAB --")
Xm, k = sp.symbols('X_m k')
gab = Xm*C*k*aw/(1-k*aw)/(1-k*aw+C*k*aw)
gab_aw = sp.solve(gab - X,aw)[1]
gab_aw = sp.simplify(gab_aw)
print("gab_aw = ",gab_aw)
print(sp.solve(gab - X,aw,force=True))

gab_diff = sp.simplify(sp.diff(gab_aw,X))
print(gab_diff)