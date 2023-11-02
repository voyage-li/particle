close all; clear all; clc;
data = [0.1 1.0152 -4.75613E-15
        0.2 1.06398 -5.51074E-05
        0.3 1.15985 -0.0126204
        0.4 1.28506 -0.066128
        0.5 1.41566 -0.153359
        0.6 1.54571 -0.26411
        0.7 1.67387 -0.392401
        0.8 1.7999 -0.534552
        0.9 1.92387 -0.688109
        1.0 2.0459 -0.85133
        1.5 2.63233 -1.77571
        2 3.18914 -2.8272];
id = 7;
k = data(id, 1);
wr = data(id, 2); wi = data(id, 3);

% parameters
L = 2 * pi / k; dt = .02; nt = 3000; ntout = 200; ng = 32; np = 30000;
vb = 1.0; xp1 = 1.0e-2; vp1 = 0.0;
vt = 0.3; % note: the normalization sqrt(2) will be found in randn()
wp = 1; qm = -1;
q = wp ^ 2 / (qm * np / L); rho_back = -q * np / L; dx = L / ng;

% initial loading for the 2 Stream instability
xp = linspace(0, L, np)';
vp = vt * randn(np, 1) + (1 - 2 * mod([1:np]', 2)) .* vb;
% vp=vt*randn(np,1); % randn is {exp[-(x-mu)^2/(2*sgm^2)]}/[sgm*sqrt(2*pi)]

% Perturbation
vp = vp + vp1 * cos(k * xp);
xp = xp + xp1 * cos(k * xp);
p = 1:np; p = [p p];

% Main computational cycle
h = figure('Unit', 'Normalized', 'position', ...
    [0.02 0.3 0.6 0.6], 'DefaultAxesFontSize', 15);

for it = 1:nt
    % apply periodic bc on the particle positions
    xp = xp ./ L + 10.0; xp = L .* (xp - floor(xp));

    % diagnosing
    if (mod(it, nt / 4) == 1)
        subplot(2, 2, floor(4 * it / nt) + 1);
        plot(xp(1:2:end), vp(1:2:end), 'r.', xp(2:2:end), ...
            vp(2:2:end), 'b.', 'Markersize', 2);
        axis([0, L, -3 * (abs(vt) + abs(vb)), 3 * (abs(vt) + abs(vb))]);
        title(['Phase space plotting, t=', num2str((it - 1) * dt)]);
        xlabel('xp'); ylabel('vp'); pause(0.2);
        %         print(gcf, '-dpng', ['vp-x,t=',num2str(it*dt),'.png']);
    end

    % update xp
    xp = xp + vp * dt;

    % projection p->g
    g1 = floor(xp / dx - .5) + 1; g = [g1; g1 + 1];
    fraz1 = 1 - abs(xp / dx - g1 + .5);
    fraz = [fraz1; 1 - fraz1];

    % apply bc on the projection
    out = (g < 1); g(out) = g(out) + ng;
    out = (g > ng); g(out) = g(out) - ng;
    mat = sparse(p, g, fraz, np, ng);
    rho = full((q / dx) * sum(mat))' + rho_back;

    % computing fields, dE/dx
    Eg = zeros(ng, 1);

    for j = 1:ng - 1
        Eg(j + 1) = Eg(j) + (rho(j) + rho(j + 1)) * dx / 2;
    end

    Eg(1) = Eg(ng) + rho(ng) * dx;
    Eg = Eg - mean(Eg);

    % projection q->p and update of vp
    vp = vp + mat * qm * Eg * dt;

    EEk(it) = 0.5 * abs(q) * sum(vp .^ 2); % kinetic energy
    EEf(it) = 0.5 * sum(Eg .^ 2) * dx; % potential energy
    t(it) = it * dt;
end

%%
h = figure('Unit', 'Normalized', 'position', ...
    [0.02 0.4 0.6 0.3], 'DefaultAxesFontSize', 15);
subplot(121); plot(t, EEk, t, EEf, t, EEk + EEf, 'r:', 'LineWidth', 2);
% title(['(a) k=',num2str(k),', \omega_{theory}=',...
%     num2str(wr+1i*wi)],'fontsize',15);
xlabel('t'); ylabel('Energy'); legend('E_k', 'E_e', 'E_{tot}', 4);
legend('boxoff');

subplot(122);

% Find the corresponding indexes of the extreme max values
lndE = log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
it0 = floor(nt * 1/20); it1 = floor(nt * 17/20);
yy = lndE(it0:it1);
extrMaxIndex = find(diff(sign(diff(yy))) == -2) + 1;
t1 = t(it0 + extrMaxIndex(1)); t2 = t(it0 + extrMaxIndex(end));
y1 = yy(extrMaxIndex(1)); y2 = yy(extrMaxIndex(end));
plot(t, lndE, 'LineWidth', 2);
% plot(t,lndE,[t1,t2],[y1,y2],'r*--','LineWidth',2);
omega = pi / ((t2 - t1) / (length(extrMaxIndex) - 1));
gammas = (real(y2) - real(y1)) / (t2 - t1);
xlabel('t'); ylabel('ln({E_e^{1/2}})');
% title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas)]);
axis tight;
