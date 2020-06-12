import sympy as sp
from sympy.printing.fcode import print_fcode

# Define symbols
Xdb, T = sp.symbols('Xdb T')
a1, a2, b1, b2 = sp.symbols('a1 a2 b1 b2')

# Oswin sorption isotherm
a = a1 + a2*T
b = b1 + b2*T
aw = (Xdb/a)**(1/b)/(1 + (Xdb/a)**(1/b))

# Derivative with respect to Xdb
daw = sp.simplify(sp.diff(aw,Xdb))

# Simplify common subexpressions
rep,red_expr = sp.cse(sp.Matrix([aw,daw]))

# Print reduced variables
for var, expr in rep:
    print_fcode(expr,assign_to=var,source_format='free')

# Print reduced expressions
result = sp.MatrixSymbol('res',2,1)
print_fcode(red_expr[0],assign_to=result,source_format='free')


# !> Oswin sorption isotherm
# subroutine oswin(Xdb,T,a1,a2,b1,b2,aw,daw)
#   real(dp), intent(in) :: Xdb,T,a1,a2,b1,b2
#     !! Input parameters (moisture,temperature,coefficients)
#   real(dp), intent(out) :: aw
#   real(dp), intent(out), optional :: daw
#     !! Output (water activity, slope of water activity)
  
#   real(dp) :: x0,x1,x2    ! Temporary variables

#   x0 = 1.0_dp/(T*b2 + b1)
#   x1 = (Xdb/(T*a2 + a1))**x0
#   x2 = x1 + 1
#   aw = x1/x2
#   if (present(daw)) daw = x0*x1/(Xdb*x2**2)
# end subroutine



