% Исследование ошибок пеленгации при наличии ГБШ в приемном канале
%clear;
close all;

% 1. Формирование ДНА
teta = -pi/10 : pi/180/10 : pi/10; % диапазон -18..+18 град.
ind = 191; % индекс для визуализации отклонения 1 град.
Da = length(teta);
syms x % Определение угла разнесения ДН
equ = sin(12*x)/12/x == 1/sqrt(2);
teta0 = double(vpasolve(equ, x)); % Угол разнесения ДН
f1 = F(teta0 - teta); % Функция F - 
f2 = F(teta0 + teta); % амплитудная ДН

figure(1) % Визуализация ДН в полярной СК
teta_fc = -pi : pi/180 : pi; % для построений
f1_fc = F(teta0 - teta_fc);  % диапазон углов
f2_fc = F(teta0 + teta_fc);  % от -180 до 180
polarplot(teta_fc, f1_fc, 'm', teta_fc, f2_fc, 'b');
legend('f_{1}(\theta)', 'f_{2}(\theta)', 'FontSize', 14)
title('Диаграммы направленности антенн', 'FontSize', 14)

% 2. Формирование сигнала цели
Ds = 100000; % Кол-во временных отсчётов
t = linspace(0, 3e-6, Ds);
Ts = t(2) - t(1); % Шаг дискретизации
E = sin(2*pi*1e9*t); % Сигнал цели

% 3. Формирование аддитивной смеси с ГБШ
S = RandStream('mt19937ar', 'Seed', (5499 - count_stat));
% Значение ОСШ
SNR = -30 : 5 : 10; % диапазон -30..10 дБ
%SNR = zeros(1,9);
dB0 = 7; % индекс для визуализации шума 0 дБ
Dn = length(SNR);
for k = 1 : Dn
    En(k, :) = awgn(E, SNR(k), 'measured', S);
    reset(S);
end

figure(2) % Визуализация сигнала цели (с ГБШ и без)
plot(t, E, 'k', t, En(dB0, :), 'k:');
legend('E(t)', 'E(t)+N(t)', 'FontSize', 14)
title('Сигнал от цели (уровень шума 0 дБ)', 'FontSize', 14)
xlabel('t, с', 'FontSize', 14)
xlim([0 3e-9])
ylabel('U, В', 'FontSize', 14)
ylim([-2.5 2.5])

% 4. Формирование принятых сигналов
u1 = zeros(Ds, Da); % Выделение памяти
u2 = u1;
U1 = zeros(Ds, Da, Dn);
U2 = U1;

for k = 1 : Dn
    for j = 1 : Da
        for i = 1 : Ds
            % В отсутствие шума
            u1(i, j) = E(1, i) * f1(1, j); % в 1 канале
            u2(i, j) = E(1, i) * f2(1, j); % во 2 канале
            % При наличии шума
            U1(i, j, k) = En(k, i) * f1(1, j); % в 1 канале
            U2(i, j, k) = En(k, i) * f2(1, j); % во 2 канале
        end
    end
end

clear E;
clear En;

% 5. Формирование сигнала гетеродина
Uget = sin(2*pi*(1e6+1e9)*t);

% 6. Формирование сигналов на смесителях
usm1 = zeros(Ds, Da); % Выделение памяти
usm2 = usm1;
Usm1 = zeros(Ds, Da, Dn);
Usm2 = Usm1;

for k = 1 : Dn
    for j = 1 : Da
        for i = 1 : Ds
            % В отсутствие шума
            usm1(i, j) = Uget(1, i) * u1(i, j); % в 1 канале
            usm2(i, j) = Uget(1, i) * u2(i, j); % во 2 канале
            % При наличии шума
            Usm1(i, j, k) = Uget(1, i) * U1(i, j, k); % в 1 канале
            Usm2(i, j, k) = Uget(1, i) * U2(i, j, k); % во 2 канале
        end
    end
end

figure(3) % Визуализация сигналов на смесителях
plot(t, Usm1(:, ind, dB0), ':m', t, Usm2(:, ind, dB0), ':b', t, usm1(:, ind), '-m', t, usm2(:, ind), '-b');
legend({'U_{sm1(N)}(t, 1\circ)', 'U_{sm2(N)}(t, 1\circ)', 'U_{sm1}(t, 1\circ)', 'U_{sm2}(t, 1\circ)'}, 'NumColumns', 2, 'FontSize', 14)
title('Сигналы на смесителях', 'FontSize', 14)
xlabel('t, с', 'FontSize', 14)
xlim([0 3e-6])
ylabel('U, В', 'FontSize', 14)
ylim([-5 5])

clear u1;
clear u2;
clear U1;
clear U2;
clear Uget;

% 7. Переход в область частот
fs = 1/Ts; % частота дискретизации                      
f = (- Ds/2 : Ds/2 - 1)*fs/Ds; % сдвинутые частоты
p = [1+2i 3i];
y1 = zeros(Ds, Da, 'like', p);
y2 = y1;
Y1 = zeros(Ds, Da, Dn, 'like', p);
Y2 = Y1;
y1 = fftshift(fft(usm1), 1); % БПФ со сдвигом нулевой
y2 = fftshift(fft(usm2), 1); % частотной составляющей
Y1 = fftshift(fft(Usm1), 1); 
Y2 = fftshift(fft(Usm2), 1); 

clear usm1;
clear usm2;
clear Usm1;
clear Usm2;

figure(4) % Визуализация спектров сигналов
plot(f, abs(Y1(:, ind, dB0)), '-m*', f, abs(Y2(:, ind, dB0)), '-b*', f, abs(y1(:, ind)), '-mo', f, abs(y2(:, ind)), '-bo');
legend({'|S_{sm1(N)}(f, 1\circ)|', '|S_{sm2(N)}(f, 1\circ)|', '|S_{sm1}(f, 1\circ)|', '|S_{sm2}(f, 1\circ)|'}, 'Location', 'southeast', 'FontSize', 14)
title('Спектры сигналов', 'FontSize', 14)
xlabel('f, Гц', 'FontSize', 14)
xlim([-3e9 3e9])
ylabel('|S(f)|', 'FontSize', 14)
ylim([0 2e4])

% 8. Описание ФНЧ Бесселя
omega = 2*pi*f; 
% Фильтр Бесселя 1-го порядка, f0 = 8e6
[H, W] = besself(1, 8e6*2*pi);
[h, w] = freqs(H, W, omega);

figure(5) % Визуализация частотных характеристик
freqs(H, W);
title('Частотные характеристики ФНЧ', 'FontSize', 14)
subplot(2, 1, 1)
xlabel('\omega, рад/с', 'FontSize', 14)
ylabel('|H(\omega)|', 'FontSize', 14)
subplot(2, 1, 2)
xlabel('\omega, рад/с', 'FontSize', 14)
ylabel('\phi(\omega), \circ', 'FontSize', 14)

% 9. Выделение промежуточной частоты
uf1 = zeros(Ds, Da, 'like', p); % Выделение памяти
uf2 = uf1;
for j = 1 : Da
    for i = 1 : Ds
        % В отсутствие шума
        uf1(i, j) = y1(i, j)*h(i); % Фильтрация 1 канала
        uf2(i, j) = y2(i, j)*h(i); % Фильтрация 2 канала
    end
end
clear y1;
clear y2;

Uf1 = zeros(Ds, Da, Dn, 'like', p);
Uf2 = Uf1;
for k =1 : Dn
    for j = 1 : Da
        for i = 1 : Ds
            % При наличии шума
            Uf1(i, j, k) = Y1(i, j, k)*h(i); % Фильтрация 1 канала
            Uf2(i, j, k) = Y2(i, j, k)*h(i); % Фильтрация 2 канала
        end
    end
end
clear Y1;
clear Y2;

figure(6) % Визуализация отфильтрованных спектров
subplot(2, 1, 1)
plot(f, abs(Uf1(:, ind, dB0)), '-m*', f, abs(Uf2(:, ind, dB0)), '-b*', f, abs(uf1(:, ind)), '-mo', f, abs(uf2(:, ind)), '-bo');
legend('|S_{if1(N)}(f, 1\circ)|', '|S_{if2(N)}(f, 1\circ)|', '|S_{if1}(f, 1\circ)|', '|S_{if2}(f, 1\circ)|', 'FontSize', 14)
title('Спектры сигналов после фильтрации', 'FontSize', 14)
xlabel('f, Гц', 'FontSize', 14)
xlim([-3e9 3e9])
ylabel('|S(f)|', 'FontSize', 14)
ylim([0 2e4])
subplot(2, 1, 2)
plot(f, abs(Uf1(:, ind, dB0)), '-m*', f, abs(Uf2(:, ind, dB0)), '-b*', f, abs(uf1(:, ind)), '-mo', f, abs(uf2(:, ind)), '-bo');
xlabel('f, Гц', 'FontSize', 14)
xlim([-3e6 3e6])
ylabel('|S(f)|', 'FontSize', 14)
ylim([0 2e4])

% 10. Переход во временную область
x1 = zeros(Ds, Da, 'like', p);
x2 = x1;
x1 = ifft(ifftshift(uf1, 1)); % сдвиг, ОПФ
clear uf1;
x2 = ifft(ifftshift(uf2, 1));
clear uf2;

X1 = zeros(Ds, Da, Dn, 'like', p);
X2 = X1;
X1 = ifft(fftshift(Uf1, 1));
clear Uf1;
X2 = ifft(ifftshift(Uf2, 1));
clear Uf2;

uph1 = real(x1); % Выделение действительных частей
uph2 = real(x2);
Uph1 = real(X1); % Выделение действительных частей
Uph2 = real(X2);

clear x1;
clear x2;
clear X1;
clear X2;

figure(7) % Визуализация сигналов на промежуточных частотах
plot(t, Uph1(:, ind, dB0), '--m*', t, Uph2(:, ind, dB0), '--b*', t, uph1(:, ind), '-m', t, uph2(:, ind), '-b', 'MarkerIndices', 1:1000:Ds);
legend('U_{if1(N)}(t, 1\circ)', 'U_{if2(N)}(t, 1\circ)', 'U_{if1}(t, 1\circ)', 'U_{if2}(t, 1\circ)', 'FontSize', 14)
title('Сигналы на промежуточных частотах', 'FontSize', 14)
xlabel('t, с', 'FontSize', 14)
xlim([0 3e-6])
ylabel('U, В', 'FontSize', 14)
ylim([-0.5 0.5])

% 11. Применение ЛУПЧ
ulu1 = zeros(Ds, Da); % Выделение памяти
ulu2 = ulu1;
Ulu1 = zeros(Ds, Da, Dn);
Ulu2 = Ulu1;

for k = 1 : Dn
    for j = 1 : Da
        % В отсутствие шума
        ulu1(:, j) = LU(uph1(:, j)); % Усиление 1 канала
        ulu2(:, j) = LU(uph2(:, j)); % Усиление 2 канала
        % При наличии шума
        Ulu1(:, j, k) = LU(Uph1(:, j, k)); % Усиление 1 канала
        Ulu2(:, j, k) = LU(Uph2(:, j, k)); % Усиление 2 канала
    end
end

clear uph1;
clear uph2;
clear Uph1;
clear Uph2;

% 12. Применение амплитудных детекторов
uad1 = zeros(Ds, Da); % Выделение памяти
uad2 = uad1;
Uad1 = zeros(Ds, Da, Dn);
Uad2 = Uad1;

for k = 1 : Dn
    for j = 1 : Da
        % В отсутствие шума
        uad1(:, j) = AD(ulu1(:, j), t); % Детектирование 1 канала
        uad2(:, j) = AD(ulu2(:, j), t); % Детектирование 2 канала
        % При наличии шума
        Uad1(:, j, k) = AD(Ulu1(:, j, k), t); % Детектирование 1 канала
        Uad2(:, j, k) = AD(Ulu2(:, j, k), t); % Детектирование 2 канала
    end
end

figure(8) % Визуализация усиленных и детектированных сигналов
hold on
plot(t, Ulu1(:, ind, dB0), '-m*', t, Ulu2(:, ind, dB0), '-b*', t, ulu1(:, ind), '-m', t, ulu2(:, ind), '-b', 'MarkerIndices', 1:1000:Ds);
plot(t, Uad1(:, ind, dB0), '--mo', t, Uad2(:, ind, dB0), '--bo', t, uad1(:, ind), '--m', t, uad2(:, ind), '--b', 'MarkerIndices', 1:1000:Ds);
hold off
box on
title('Сигналы после усиления и детектирования', 'FontSize', 14)
legend('U_{lu1(N)}(t, 1\circ)', 'U_{lu2(N)}(t, 1\circ)', 'U_{lu1}(t, 1\circ)', 'U_{lu2}(t, 1\circ)', 'U_{ad1(N)}(t, 1\circ)', 'U_{ad2(N)}(t, 1\circ)', 'U_{ad1}(t, 1\circ)', 'U_{ad2}(t, 1\circ)', 'Location', 'eastoutside', 'FontSize', 14)
xlabel('t, с', 'FontSize', 14)
xlim([0 3e-6])
ylabel('U, В', 'FontSize', 14)
ylim([0 2])

clear ulu1;
clear ulu2;
clear Ulu1;
clear Ulu2;

% 13. Описание вычитающего устройства
uvu = uad1 - uad2; % В отсутствие шума
Uvu = Uad1 - Uad2; % При наличии шума

clear uad1;
clear uad2;
clear Uad1;
clear Uad2;

figure(9) % Визуализация сигнала на ВУ
hold on
plot(t, Uvu(:, ((ind - 10) : 10 : (ind + 70)), dB0), 'Color', '#D95319', 'LineWidth', 2, 'LineStyle', ':');
plot(t, uvu(:, ((ind - 10) : 10 : (ind + 70))), 'Color', '#7E2F8E', 'LineWidth', 2);
hold off
box on
grid on
title('Сигналы на выходе вычитающего устройства', 'FontSize', 14)
legend({'U_{ву(N)}(t, 0\circ)', 'U_{ву(N)}(t, +1\circ)', 'U_{ву(N)}(t, +2\circ)', 'U_{ву(N)}(t, +3\circ)', 'U_{ву(N)}(t, +4\circ)', 'U_{ву(N)}(t, +5\circ)', 'U_{ву(N)}(t, +6\circ)', 'U_{ву(N)}(t, +7\circ)', 'U_{ву(N)}(t, +8\circ)', 'U_{ву}(t, 0\circ)', 'U_{ву}(t, +1\circ)', 'U_{ву}(t, +2\circ)', 'U_{ву}(t, +3\circ)', 'U_{ву}(t, +4\circ)', 'U_{ву}(t, +5\circ)', 'U_{ву}(t, +6\circ)', 'U_{ву}(t, +7\circ)', 'U_{ву}(t, +8\circ)'}, 'Location', 'eastoutside', 'NumColumns', 2, 'FontSize', 14)
xlabel('t, с', 'FontSize', 14)
xlim([0 3e-6])
ylabel('U, В', 'FontSize', 14)
ylim([0 2])

% 14. Формирование пеленгационных характеристик
deg = teta*180/pi;
di = 61 : 10: 301; % Выборка от -12 до +12 град.

figure(10) % Визуализация ПХ
hold on
plot(deg, Uvu(1, :, dB0), 'Color', '#D95319', 'LineStyle', ':', 'LineWidth', 2);
plot(deg, uvu(1, :), 'Color', '#7E2F8E', 'LineWidth', 2);
hold off
box on
grid on
title('Пеленгационные характеристики', 'FontSize', 14)
legend({'U_{ву(N)}(\theta)', 'U_{ву}(\theta)'}, 'Location', 'southeast', 'FontSize', 14)
xlabel('\theta, \circ', 'FontSize', 14)
xlim([-15 15])
ylabel('U, В', 'FontSize', 14)
ylim([-3.5 3.5])

% 15. Вычисление ошибок пеленгации 
if (count_stat == 0)
    Dfe_V = zeros(length(di), Dn + 1); % Выделение памяти
    Dfe_D = Dfe_V;
end
De = length(di);

% Поиск максимальных значения отклонения в [В] и [град.]
for k = 1 : Dn
    for j = 1 : De
        Dfe_V(j, 1) = j - 13;
        Dfe_D(j, 1) = j - 13;
        [V, D] = DFEm(uvu(:, :), Uvu(:, di(j), k), di(j));
        Dfe_V(j, k + 1) = (Dfe_V(j, k + 1) * count_stat + V) / (count_stat + 1); 
        Dfe_D(j, k + 1) = (Dfe_D(j, k + 1) * count_stat + D) / (count_stat + 1);
    end
end
count_stat = count_stat + 1;

clear Uvu;
clear uvu;

figure(11) % Визуализация зависимости ошибки от угла рассогласования
subplot(1,2,1)
plot(Dfe_V(:, 1), Dfe_V(:, 2:end), 'LineWidth', 1);
grid on
title('Смещение выходного сигнала', 'FontSize', 14)
legend('ОСШ = -30 дБ', 'ОСШ = -25 дБ', 'ОСШ = -20 дБ', 'ОСШ = -15 дБ', 'ОСШ = -10 дБ', 'ОСШ = -5 дБ', 'ОСШ = 0 дБ', 'ОСШ = 5 дБ', 'ОСШ = 10 дБ', 'Location', 'southoutside', 'Numcolumns', 3, 'FontSize', 14)
xlabel('\theta, \circ', 'FontSize', 14)
xlim([-12 12])
ylabel('U_{err}(\theta), В', 'FontSize', 14)
ylim([0 2.5])
subplot(1,2,2)
plot(Dfe_D(:, 1), Dfe_D(:, 2:end), 'LineWidth', 1);
grid on
title('Ошибка пеленгации', 'FontSize', 14)
xlabel('\theta, \circ', 'FontSize', 14)
xlim([-12 12])
ylabel('\theta_{err}(\theta), \circ', 'FontSize', 14)
ylim([0 5.5])

figure(12) % Визуализация зависимости ошибки от ОСШ
subplot(1,2,1)
plot(SNR, Dfe_V(:, 2:end), 'LineWidth', 1);
grid on
title('Смещение выходного сигнала', 'FontSize', 14)
legend('\theta = -12\circ', '\theta = -11\circ', '\theta = -10\circ', '\theta = -9\circ', '\theta = -8\circ', '\theta = -7\circ', '\theta = -6\circ', '\theta = -5\circ', '\theta = -4\circ', '\theta = -3\circ', '\theta = -2\circ', '\theta = -1\circ', '\theta = 0\circ', '\theta = 1\circ', '\theta = 2\circ', '\theta = 3\circ', '\theta = 4\circ', '\theta = 5\circ', '\theta = 6\circ', '\theta = 7\circ', '\theta = 8\circ', '\theta = 9\circ', '\theta = 10\circ', '\theta = 11\circ', '\theta = 12\circ', 'Location', 'southoutside', 'Numcolumns', 5, 'FontSize', 14)
xlabel('SNR, дБ', 'FontSize', 14)
xlim([-30 10])
ylabel('U_{err}(SNR), В', 'FontSize', 14)
ylim([0 2.5])
subplot(1,2,2)
plot(SNR, Dfe_D(:, 2:end), 'LineWidth', 1);
grid on
title('Ошибка пеленгации', 'FontSize', 14)
xlabel('SNR, дБ', 'FontSize', 14)
xlim([-30 10])
ylabel('\theta_{err}(SNR), \circ', 'FontSize', 14)
ylim([0 5.5])

% 16. Расчет значений промаха
Fsi = 1;  % для диапазона углов
Lsi = 25; % -12..+12 град.
dl = zeros((Lsi - Fsi + 1), Dn);

for x = Fsi : Lsi 
    for y = 1 : Dn
        dl(x - Fsi + 1, y) = MissTgt(x, y, teta, di, Dfe_D);
    end
end

figure(13) % Визуализация значений промаха
plot(SNR, 1000*dl(:, :), 'LineWidth', 1);
grid on
title('Значения промахов', 'FontSize', 14)
legend('\theta = -12\circ', '\theta = -11\circ', '\theta = -10\circ', '\theta = -9\circ', '\theta = -8\circ', '\theta = -7\circ', '\theta = -6\circ', '\theta = -5\circ', '\theta = -4\circ', '\theta = -3\circ', '\theta = -2\circ', '\theta = -1\circ', '\theta = 0\circ', '\theta = 1\circ', '\theta = 2\circ', '\theta = 3\circ', '\theta = 4\circ', '\theta = 5\circ', '\theta = 6\circ', '\theta = 7\circ', '\theta = 8\circ', '\theta = 9\circ', '\theta = 10\circ', '\theta = 11\circ', '\theta = 12\circ', 'Location', 'eastoutside', 'FontSize', 14)
xlabel('SNR, дБ', 'FontSize', 14)
xlim([-30 10])
ylabel('\Deltal(SNR), м', 'FontSize', 14)
ylim([0 200])