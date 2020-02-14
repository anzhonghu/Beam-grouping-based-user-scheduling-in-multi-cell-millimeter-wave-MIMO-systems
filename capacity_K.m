clear;
close all;
L = 3;
naz = 12;
nel = 100;
n_arr = naz * nel;
maz = 2;
mel = 2;
m_arr = maz * mel;
Rmin = 90;
Rmax = 100;
base = [0, 0; 0,  2 * Rmax; sqrt(3) * Rmax,  Rmax;];
Nr = 4;
P = 4;
K = 20;
variable_s = [10; 20; 30; 40; 50;];
Pvsigma2 = 1e18;
height = 10;
time_fre_resources = 20;
f = 80*10^9;%1G bandwidth
lambda = 3 * 10^8 / f;
miu = 0.5;
beta_m = 10.3 * pi / 180;
Nite = 1e4;
capacity = zeros(length(variable_s),  6);
jain = zeros(length(variable_s),  6);
angle = zeros(1, 2);
U = zeros(n_arr, n_arr);
for nx = 0 : naz-1
    for ny = 0 : nel-1
        angle(1, 2) = (-1+2*ny/nel);%el
        angle(1, 1) = (-1+2*nx/naz);%az
        n = ny * naz + nx + 1;
        for mx = 0 : naz-1
            for my = 0 : nel-1
                m = my * naz + 1 + mx;
                U(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(naz-1)) * angle(1, 1) + (my-0.5*(nel-1)) * angle(1, 2))) / sqrt(n_arr);
            end
        end
    end
end
Um = zeros(m_arr, m_arr);
for nx = 0 : maz-1
    for ny = 0 : mel-1
        angle(1, 2) = (-1+2*ny/mel);%el
        angle(1, 1) = (-1+2*nx/maz);%az
        n = ny * maz + nx + 1;
        for mx = 0 : maz-1
            for my = 0 : mel-1
                m = my * maz + 1 + mx;
                Um(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(maz-1)) * angle(1, 1) + (my-0.5*(mel-1)) * angle(1, 2))) / sqrt(m_arr);
            end
        end
    end
end
for variable_n = 1 : length(variable_s)
    K = variable_s(variable_n, 1);
    for ii = 1 : Nite
        H = zeros(m_arr, K*L*L*n_arr);
        p_store = zeros(K, L);
        beta_store = zeros(K, P*L);
        betab_store = zeros(K, P*L*L);
        thetab_store = zeros(K, P*L);
        thetabb_store = zeros(K, P*L*L);
        phib_store = zeros(K, P*L);
        phibb_store = zeros(K, P*L*L);
        thetam_store = zeros(K, P*L);
        phim_store = zeros(K, P*L);
        pos_store = zeros(K*L, 2);
        for l = 1 : L
            for k = 1 : K
                pos_temp = zeros(1, 2);
                while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
                    pos_temp(1, 1) = rand(1, 1) * Rmax;
                    pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
                end
                switch l
                    case 1
                        angle_temp = atan(pos_temp(1, 2) / pos_temp(1, 1)) + pi / 3;
                    case 2
                        angle_temp = atan(pos_temp(1, 2) / pos_temp(1, 1)) - pi / 3;
                    case 3
                        angle_temp = atan(pos_temp(1, 2) / pos_temp(1, 1)) + pi;
                    otherwise
                end
                d_temp = norm(pos_temp);
                pos_temp(1, 1) = d_temp * cos(angle_temp);
                pos_temp(1, 2) = d_temp * sin(angle_temp);
                pos_store((l-1)*K+k, :) = pos_temp;%relative to the serving BS
            end
        end
        pos = zeros(K, 3);
        theta = zeros(K, P);
        phi = zeros(K, P);
        beta = zeros(K, P);
        theta_m = zeros(K, P);
        phi_m = zeros(K, P);
        block_store = zeros(K, L*L);
        for l1 = 1  : L
            for l2 = 1 : L
                for k = 1 : K
                    for p = 1 : P
                        if 1==p
                            pos_temp = pos_store((l2-1)*K+k, :) + base(l2, :)-base(l1, :);%%relative position of MS in the l2-th cell to l1-th BS
                            pos(k, 3) = norm(pos_temp);
                            pos(k, 3) = norm([pos(k, 3), height]);%distance
                            if rand(1,1)>0.1
                                phi(k, p) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
                                switch l1
                                    case 1
                                        phi(k, p) = phi(k, p) - pi / 3;%az
                                    case 2
                                        phi(k, p) = phi(k, p) + pi / 3;%az
                                    case 3
                                        phi(k, p) = phi(k, p) - pi;%az
                                    otherwise
                                end
                                theta(k, p) = asin((pos_temp(1, 1) * sin(beta_m) - height * cos(beta_m)) / pos(k, 3));%el
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                phi_m(k, p) = pi - phi(k, p);
                                theta_m(k, p) = pi - theta(k, p);
                            else
                                block_store(k, (l1-1)*L+l2) = 1;
                                beta(k, p) = 0;
                                continue;
                            end
                        else
                            %take the reflector inside the cell
                            pos_temp = zeros(1, 2);
                            while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
                                pos_temp(1, 1) = rand(1, 1) * Rmax;
                                pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
                            end
                            switch l2
                                case 1
                                    angle_temp = atan(pos_temp(1, 2) / pos_temp(1, 1)) + pi / 3;
                                case 2
                                    angle_temp = atan(pos_temp(1, 2) / pos_temp(1, 1)) - pi / 3;
                                case 3
                                    angle_temp = atan(pos_temp(1, 2) / pos_temp(1, 1)) + pi;
                                otherwise
                            end
                            d_temp = norm(pos_temp);
                            pos_temp(1, 1) = d_temp * cos(angle_temp);
                            pos_temp(1, 2) = d_temp * sin(angle_temp);
                            pos_temp = pos_temp  + base(l2, :)-base(l1, :);
                            pos(k, 3) = norm(pos_temp);
                            pos(k, 3) = norm([pos(k, 3), height]);%distance
                            phi(k, p) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
                            switch l1
                                case 1
                                    phi(k, p) = phi(k, p) - pi / 3;%az
                                case 2
                                    phi(k, p) = phi(k, p) + pi / 3;%az
                                case 3
                                    phi(k, p) = phi(k, p) - pi;%az
                                otherwise
                            end
                            theta(k, p) = asin((pos_temp(1, 1) * sin(beta_m) - height * cos(beta_m)) / pos(k, 3));%el
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            phi_m(k, p) = pi - phi(k, p);
                            theta_m(k, p) = pi - theta(k, p);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Aaz = -min(12 * phi(k, p)^2 / (70/180*pi)^2, 25);
                        Ael = -min(12 * theta(k, p)^2 / (7/180*pi)^2, 20);
                        D0 = -min(-Aaz-Ael, 25);
                        D0 = 10^(D0*0.1);
                        if p < 2
                            if block_store(k, (l1-1)*L+l2) > 0
                            else
                                beta(k, p) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
                            end
                        else
                            beta(k, p) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi) * 10^((-rand(1,1) * 5 - 15)*0.05);%-15~-20dB loss
                        end
                        hb = zeros(n_arr, 1);
                        for nx = 0 : naz-1
                            for ny = 0 : nel-1
                                n = ny * naz+ 1 + nx;
                                hb(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz-1)) * cos(theta(k, p)) * sin(phi(k, p)) + (ny-0.5*(nel-1)) * sin(theta(k, p))));
                            end
                        end
                        hm = zeros(m_arr, 1);
                        for nx = 0 : maz-1
                            for ny = 0 : mel-1
                                n = ny * maz+ 1 + nx;
                                hm(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(maz-1)) * cos(theta_m(k, p)) * sin(phi_m(k, p)) + (ny-0.5*(mel-1)) * sin(theta_m(k, p))));
                            end
                        end
                        H(:,(l1-1)*K*L*n_arr+(l2-1)*K*n_arr+(k-1)*n_arr+1:(l1-1)*K*L*n_arr+(l2-1)*K*n_arr+k*n_arr) = ...
                            H(:,(l1-1)*K*L*n_arr+(l2-1)*K*n_arr+(k-1)*n_arr+1:(l1-1)*K*L*n_arr+(l2-1)*K*n_arr+k*n_arr) + beta(k, p) * hm * hb';
                        thetabb_store(k, (l1-1)*L*P+(l2-1)*P+p) = theta(k, p);
                        phibb_store(k, (l1-1)*L*P+(l2-1)*P+p) = phi(k, p);
                        betab_store(k, (l1-1)*L*P+(l2-1)*P+p) = abs(beta(k, p));
                        if l1==l2
                            thetab_store(k, (l1-1)*P+p) = theta(k, p);
                            phib_store(k, (l1-1)*P+p) = phi(k, p);
                            thetam_store(k, (l1-1)*P+p) = theta_m(k, p);
                            phim_store(k, (l1-1)*P+p) = phi_m(k, p);
                            beta_store(k, (l1-1)*P+p) = abs(beta(k, p));
                        else
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        greedy_appr;%1
        maxmi_appr;%2
        pf_appr;%3
        proposed;%4-6
        disp([variable_n, ii])
    end
end
capacity = capacity / Nite;
jain = jain / Nite;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1 = subplot(1,2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(variable_s, capacity(1:length(variable_s), 1), 'k--s','LineWidth',1,'MarkerSize',10)
hold on
plot(variable_s, capacity(1:length(variable_s), 2), 'k--*','LineWidth',1,'MarkerSize',10)
plot(variable_s, capacity(1:length(variable_s), 3), 'k--o','LineWidth',1,'MarkerSize',12)
plot(variable_s, capacity(1:length(variable_s), 4), 'k-^','LineWidth',1,'MarkerSize',10)
xlim([min(variable_s), max(variable_s)])
le = legend('Greedy','Max-min','PF','Proposed', 'Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',variable_s)
xlabel('Number of MSs','Fontname','Times')
ylabel('Sum rate (bps/Hz)','Fontname','Times')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
h2 = subplot(1,2,2);
plot(variable_s, jain(1:length(variable_s), 1), 'k--s','LineWidth',1,'MarkerSize',10)
hold on
plot(variable_s, jain(1:length(variable_s), 2), 'k--*','LineWidth',1,'MarkerSize',10)
plot(variable_s, jain(1:length(variable_s), 3), 'k--o','LineWidth',1,'MarkerSize',12)
plot(variable_s, jain(1:length(variable_s), 4), 'k-^','LineWidth',1,'MarkerSize',10)
xlim([min(variable_s), max(variable_s)])
set(le,'Fontname','Times')
set(gca,'XTick',variable_s)
xlabel('Number of MSs','Fontname','Times')
ylabel('Jain''s fairness index','Fontname','Times')
grid on