from sympy import *

R, h, f, dh, df, dhh, dff, dhf = symbols('R h f dh df dhh dff dhf')

x = (R+h)*cos(f)
y = (R+h)*sin(f)

hx = cos(f)
hy = sin(f)

fx = -sin(f)/(R+h)
fy =  cos(f)/(R+h)

#(hx dh + fx df) (hy dh + fy df)

dxx = hx*hx * dhh + hx * diff(hx,h) * dh + \
      hx*fx * dhf + hx * diff(fx,h) * df + \
      fx*hx * dhf + fx * diff(hx,f) * dh + \
      fx*fx * dff + fx * diff(fx,f) * df

dyy = hy*hy * dhh + hy * diff(hy,h) * dh + \
      hy*fy * dhf + hy * diff(fy,h) * df + \
      fy*hy * dhf + fy * diff(hy,f) * dh + \
      fy*fy * dff + fy * diff(fy,f) * df

dxy = hx*hy * dhh + hx * diff(hy,h) * dh + \
      hx*fy * dhf + hx * diff(fy,h) * df + \
      fx*hy * dhf + fx * diff(hy,f) * dh + \
      fx*fy * dff + fx * diff(fy,f) * df

tr = simplify(dxx+dyy) # Trace
dethess = simplify(dxx*dyy-dxy*dxy) # Hessian determinant


print(tr)
print(dethess)





