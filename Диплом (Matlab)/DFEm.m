function [V, D] = DFEm(U, Un, ang)
% Вычисление смещения и ошибки
Err = Un - U(:, ang);
max_err_V = 0;
maxu=0;
tim_err = 1;
% Вычисление максимального смещенния сигнала
for i = 1:5000
    if (abs(Err(i)) > abs(max_err_V))
        max_err_V = Err(i);
        tim_err = i;
    end 
end
ym_V = Un(tim_err);
% Иллюстрация вычислений
figure(14)
hold on
plot(1:361, U(tim_err, :))
plot( ang, Un(tim_err), '*')
hold off
xlim([0 361])
ylim([-3.5 3.5])
% Вычисление ошибки пеленгации
if ang >= 181
    for j = 321 : -1 : 181
        if (U(tim_err, j) >= maxu)
            maxu = U(tim_err, j);
            maxuj = j;
        end
    end
    max_err_D = maxuj;
    for j = 361 : -1 : 181
        if (U(tim_err, j) >= ym_V)
            max_err_D = j;
        end
    end
else
    for j = 41 : 181      
        if (U(tim_err, j) <= maxu)
            maxu = U(tim_err, j);
            maxuj = j;
        end
    end
    max_err_D = maxuj;
    for j = 1 : 181
        if (U(tim_err, j) <= ym_V)
            max_err_D = j;
        end
    end
end
V = abs(max_err_V);
D = abs((max_err_D - ang) / 10);
end