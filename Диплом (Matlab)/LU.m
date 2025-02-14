function y = LU(x)
% Описание ЛУПЧ
k0 = 15;
x0 = 0.025;
al = 1.2;
y = zeros(1, length(x));
for i = 1 : length(x)
    if (x(i) >= x0)
        y(i) = k0*x0*(al*log(x(i)/x0) + 1);
    else
       y(i) = k0*x(i); 
    end
end
