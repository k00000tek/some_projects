function dl = MissTgt(ang_ind,snr_ind, rad, ind, Dfe_D)
% Вычисление промахов
ang = rad(ind(ang_ind));
if (ang < 0)
    err = ang - Dfe_D(ang_ind, (snr_ind + 1)) * (pi/180);
else
    err = ang + Dfe_D(ang_ind, (snr_ind + 1)) * (pi/180);
end
le = 2 * sin(err);
d = 2 * cos(err);
l = d * tan(ang);
dl = abs(le - l);
end