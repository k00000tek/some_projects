function y = AD(x, t)
% Описание АД
tau = 0.00001;
y = zeros(1, length(x));
y(1) = x(1);
iu = 1;
xu = x(1);
for i = 2 : length(x)
    if (y(i - 1) <= x(i))
        y(i) = x(i);
        iu = i;
        xu = x(i);
    else
       y(i) = exp(-t(i - iu)/tau)*xu; 
    end
end
