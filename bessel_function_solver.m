format longE;
z = 0:0.1:10;
J=-sqrt(z).*(besselj(-2/3,(2/3).*z.^(3/2))+besselj(2/3,(2/3).*z.^(3/2)))./(-besselj(-1/3,(2/3).*z.^(3/2))+besselj(1/3,(2/3).*z.^(3/2)));
plot(z,J)
x0=2.33863;
fun= @f;
y=fzero(fun,x0)

function y=f(z)
y=besselj(-1/3,(2/3).*z.^(3/2))+besselj(1/3,(2/3).*z.^(3/2));
end