format longE;
z = 0:0.1:10;
A=airy(1,-z);
plot(z,A)
x0=-1;
fun= @f;
y=fzero(fun,x0)

function y=f(z)
y=airy(1,-z);
end