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

print(dethess)


'''
Ahh = dhh + dh/(R+h)
Aff = dff/(R+h)**2

numer = (R+h)**4 * dethess

print(numer)

ans = collect(numer,df).coeff(df,2)
Fdf2 = simplify(factor(ans))
print('[df^2] ',Fdf2)
numer = simplify(numer - Fdf2 * df**2)
print(numer)

ans = collect(numer,dff).coeff(dff,1)
Fdff = simplify(factor(ans))
print('[dff] ',Fdff)
numer = simplify(numer - Fdff * dff)
print(numer)

ans = collect(numer,dhh).coeff(dhh,1)
Fdhh = simplify(factor(ans))
print('[dhh] ', Fdhh)
numer = simplify(numer - Fdhh * dhh)
print(numer)

ans = collect(dethess*(R+h)**4,dhf).coeff(dhf,2)
Fdhf2 = simplify(factor(ans))
print('[dhf^2] ',Fdhf2)

numer = simplify(numer - Fdhf2 * dhf**2)
print(numer)

rec = dhh*dh * (R+h)**3 + (dff*dhh-dhf**2)*(R+h)**2 + 2*df*dhf * (R+h) - df*df
'''

rec = dhh*dh/(R+h) + (dff*dhh-dhf**2)/(R+h)**2 + 2*df*dhf/(R+h)**3 - df*df/(R+h)**4

ans = simplify(dethess - rec)
print(ans)

pprint(rec)




