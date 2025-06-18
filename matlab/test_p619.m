% This script tests the ongoing implementation of Recommendation ITU-R
% P.619-5.
% This implementation is not yet complete and is work under progress
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    24AUG24     Ivica Stevanovic, OFCOM         Initial version

close all
clear all

p = P619;

%tl = tl_p619_single(10, 100, 1)

%% TEST: Reproduce the results from Figure 1 in Recommendation ITU-R P.619-5

XPD = linspace(0,30,101);
[Ax, Ac] = p.co_and_cross_polar_attenuation(XPD);
figure
plot(XPD, Ax, 'b', 'LineWidth', 1.5);
legendstr(1) = "Cross-polar";
hold on
plot(XPD, Ac, 'r', 'LineWidth', 1.5);
legendstr(2) = "Co-polar";
xlabel('XPD (dB)')
ylabel('Loss (dB)');
title('ITU-R P.619-5, Figure 1')
grid on
set(gca, 'XLim', [0, 30], 'YLim', [0,30]);
title('ITU-R P.619-5, Figure 1')
legend(legendstr, 'Location', 'northwest')

%% TEST: Reproduce the results from Figure 1 in Rec. ITU-R P.531-15

fGHz = logspace(-2, 1, 101);


%% TEST: Reproduce Figure 8 on beam spreading loss in Recommendation ITU-R P.619-5
h = [0 0.5 1 2 3 5];
theta = linspace(0,6,100);
figure
hold on
count = 1;
clear legendstr;
for ih = 1:length(h)
    Abs = p.beam_spreading_loss(theta, h(ih));
    plot(theta, Abs, 'LineWidth', 1.5);
    legendstr(count) = sprintf("h = %g km", h(ih));
    count  = count + 1;
end
grid
xlabel('Free-space elevation angle (deg)');
ylabel('Beam spreading loss (dB)');
title('ITU-R P.619-5, Figure 8')
legend(legendstr);

%% TEST: Reproduce the results from Figure 9 in Recommendation ITU-R P.619-5

tp = logspace(-2, 1, 101);
theta = [5, 10, 20, 35];
Nwet = 42.5;
f = 30;
Ga = 0;
figure
count = 1;
clear legendstr
for ith = 1:length(theta)
    Ast_fades = p.tropospheric_scintillation(f, tp, Nwet, theta(ith), Ga);
    
    Ast_enh   = p.tropospheric_scintillation(f, 100-tp, Nwet, theta(ith), Ga);
    loglog(tp, abs(Ast_fades), ':', 'LineWidth', 1.5);
    legendstr(count) = sprintf("F: \\theta = %g^o", theta(ith));
    count  = count + 1;
    hold on
    loglog(tp, abs(Ast_enh), '-', 'LineWidth', 1.5);
    legendstr(count) = sprintf("E: \\theta = %g^o", theta(ith));
    count  = count + 1;
end
set(gca,'YLim', [0.1,10]);
legend(legendstr, 'Location', 'southwest', 'NumColumns', 2);
grid;
xlabel('Percentage time fades and enhancement exceeded')
ylabel('Enhancement, Fade (dB)')
title('ITU-R P.619-5, Figure 9')


%% TEST: Verify temperature and pressure profile for standard atmosphere in P.835

h = linspace(0,100,11);
T = zeros(size(h));
P = zeros(size(h));

[T, P, rho] = p.p835_reference_atmosphere(h);

Tref = 1e2*[
   2.881500000000000
   2.232520926479785
   2.166500000000000
   2.265090836113301
   2.503496461024211
   2.706500000000000
   2.470208847727968
   2.195848217750576
   1.986385762508688
   1.868673000000000
   1.950813443352469];

Pref = 1e3*[   
   1.013250000000000
   0.264998926632084
   0.055293585835330
   0.011970513284783
   0.002871516854551
   0.000797821781035
   0.000219595798590
   0.000052211125206
   0.000010525341342
   0.000001835996726
   0.000000320124364];

fail = 0;
pass = 0;
tol = 1e-8;

pass = sum(abs(T-Tref)<=tol);
fail = sum(abs(T-Tref)>tol);

fprintf(1,"Temp profile S1: %d out of %d tests passed\n", pass, pass + fail);

pass = sum(abs(P-Pref)<=tol);
fail = sum(abs(P-Pref)>tol);

fprintf(1,"Pressure profile S1: %d out of %d tests passed\n", pass, pass + fail);

%% TEST: Reproduce Figure 2 from Rec. ITU-R P.676-13
f = linspace(50, 70, 501);
h = [0, 5, 10, 15, 20];
figure

[T, P, rho] = p.p835_reference_atmosphere(h);
e = rho .* T / 216.7;          %(8) P.835
Pd = P - e; %dry air pressure
for i = 1:length(h)
    
    for ff = 1:length(f)

        [gammastd_0, gammastd_w] = p.p676d11_ga(f(ff), Pd(i), rho(i), T(i));

        Agas(ff) = gammastd_0 + gammastd_w;

    end

    semilogy(f, Agas);
    hold on
end

xlabel('Frequency (GHz)');
ylabel('Specific attenuation (dB/km)');
set(gca,'FontSize',14);
grid on
legend('0 km', '5 km', '10 km', '15 km', '20 km')
title('P.676, Figure 2')

%% TEST: Reproduce the results from Figure 4 in Recommendation ITU-R P.676-13

% First, compute the bounds for the layers from 0-100 km
h = p.p676_slant_path_geometry15();

% Second, compute the profiles for the standard atmosphere at the midpoint
% of each layer 1...N and thereform a vector of refractive indices
% vectorized form

hmid = 0.5 * ( h(1:end-1) + h(2:end) );

[T, P, rho, nstd] = p.p835_std_atm_profiles(hmid, 1);

ndry = p.p453_n(T, P, 0.*P);

e = rho .* T / 216.7;          %(8) P.835

% Compute the vector of path lengths in each layer for zenith (phi = 90)
phi = 90;
astd = p.p676_slant_path_geometry17(h, nstd, phi);
adry = p.p676_slant_path_geometry17(h, ndry, phi);

f = linspace(1,1000, 1001); % 1 GHz spacing
Agas_std = zeros(size(f));
Agas_dry = zeros(size(f));


for ff = 1:length(f)
    for i = 1:length(hmid)
        Pd = P(i)-e(i); %dry air pressure
        [gammastd_0, gammastd_w] = p.p676d11_ga(f(ff), P(i)-e(i), rho(i), T(i));
        [gammadry_0, gammadry_w] = p.p676d11_ga(f(ff), P(i),      0,      T(i));

        Agas_std(ff) = Agas_std(ff) + (gammastd_0 + gammastd_w)*  astd(i);
        Agas_dry(ff) = Agas_dry(ff) + (gammadry_0 + gammadry_w) * adry(i);
    
    end
end 

figure
semilogy(f, Agas_std, 'r');
hold on
semilogy(f, Agas_dry, 'b');
xlabel('Frequency (GHz)');
ylabel('Zenith attenuation (dB)');
set(gca,'FontSize',14);
legend('Standard','Dry')
grid on
title('P.676, Figure 4')


%% TEST: Tests for iterative method in solving equation (45)

atm_type = [1, 2, 31, 32, 41, 42];

Href = [26.078700121044676, 26.072252395625583, 26.074575727979209, 26.076356247471267, 26.072705633556325, 26.077733260472996];

fail = 0;
pass = 0;
tol = 1e-8;

for i = 1:length(atm_type)
    H(i) = p.solve45(30, -2, atm_type(i), 10);
end

pass = sum(abs(H-Href)<=tol);
fail = sum(abs(H-Href)>tol);

fprintf(1,"Equation (45): %d out of %d tests passed\n", pass, pass + fail);

%% TEST: Test the E2s path for zenith, it should give the same results and P.676
He = 0; % Altitude of the earth station antenna
Hs = 100; % Altitude of the space station antenna (it is beyond 100 km)
phi_e = 90; % elevation of the earth station antenna main beam - towards zenith)
phi_s = 90; % elevation of the space station antenna main beam - toward earth
Dphi_s = 2; % Beamwidth of the space station antenna - as they point to one another, this has no effect
std_atm = true; % we use the standard atmospheres from P.835
atm_type = 1;   % we use the standard atmosphere from Section 1.

for ff = 1:length(f)
    
   Ag_E2s_std(ff) = p.atm_attenuation_E2s( f(ff), He, Hs, phi_e, phi_s, Dphi_s, h, rho, T, P, nstd, std_atm, atm_type);
   Ag_E2s_dry(ff) = p.atm_attenuation_E2s( f(ff), He, Hs, phi_e, phi_s, Dphi_s, h, rho.*0, T, P, ndry, std_atm, atm_type);

end

figure
semilogy(f, Ag_E2s_std, 'r');
hold on
semilogy(f, Agas_std, 'b--');
xlabel('Frequency (GHz)');
ylabel('Zenith attenuation (dB)');
set(gca,'FontSize',14);
legend('P.619 E2s','P.676')
title('Standard atmosphere')
grid on

figure
semilogy(f, Ag_E2s_dry, 'r');
hold on
semilogy(f, Agas_dry, 'b--');
xlabel('Frequency (GHz)');
ylabel('Zenith attenuation (dB)');
set(gca,'FontSize',14);
legend('P.619 E2s','P.676')
title('Dry atmosphere')
grid on

%% TEST: Test the s2E path for zenith, it should give the same results and P.676
He = 0; % Altitude of the earth station antenna
Hs = 100; % Altitude of the space station antenna (it is beyond 100 km)
phi_e = -90; % elevation of the earth station antenna main beam - towards zenith)
phi_s = -90; % elevation of the space station antenna main beam - toward earth
Dphi_e = 2; % Beamwidth of the space station antenna - as they point to one another, this has no effect
std_atm = true; % we use the standard atmospheres from P.835
atm_type = 1;   % we use the standard atmosphere from Section 1.

for ff = 1:length(f)
    
   Ag_s2E_std(ff) = p.atm_attenuation_s2E( f(ff), He, Hs, phi_e, phi_s, Dphi_e, h, rho, T, P, nstd, std_atm);
   Ag_s2E_dry(ff) = p.atm_attenuation_s2E( f(ff), He, Hs, phi_e, phi_s, Dphi_e, h, rho.*0, T, P, ndry, std_atm);

end

figure
semilogy(f, Ag_s2E_std, 'r');
hold on
semilogy(f, Agas_std, 'b--');
xlabel('Frequency (GHz)');
ylabel('Zenith attenuation (dB)');
set(gca,'FontSize',14);
legend('P.619 s2E','P.676')
title('Standard atmosphere')
grid on

figure
semilogy(f, Ag_s2E_dry, 'r');
hold on
semilogy(f, Agas_dry, 'b--');
xlabel('Frequency (GHz)');
ylabel('Zenith attenuation (dB)');
set(gca,'FontSize',14);
legend('P.619 s2E','P.676')
title('Dry atmosphere')
grid on

%% TEST: test E2s with negative elevation angles phi_e
% TODO

%% TEST: Test Attachment D implementation with P453 digital maps
% TODO
f = 20;
tp = 2;
phi_e = 6;
phi_n = 46;
theta = 10;

Nwet = p.get_interp2_Nwet_Annual_time_location(tp, phi_e, phi_n);
Ga = 20;
Ast = p.tropospheric_scintillation(f,tp,Nwet,theta,Ga);


%% TEST: Test Attachment E implementation 
% TODO
Ht = 0.1;
theta1 = 5.0;  %degrees
[Dc, Hr] = p.beam_clearance(Ht, theta1);
figure
plot(Dc, Hr, 'r');
hold on
theta2 = 5.001;  %degrees
[Dc, Hr] = p.beam_clearance(Ht, theta2);
plot(Dc, Hr, 'b--');
xlabel('Distance (km)');
ylabel('Ray height (km)');
set(gca,'FontSize',14);
lstr1 = ['w/  Refraction: \theta = ' num2str(theta1)];
lstr2 = ['w/o Refraction: \theta = ' num2str(theta2)];
legend(lstr1,lstr2)
grid on


%% TEST: Pro-forma Test Attachment F implementation
% TODO
Pint = -100;
f = 20;
r1 = 100;
r2 = 200; 
Ptx = 100;
G = 10;
dtx = 10;
drx = 10;
[gamma_o, gamma_w] = p.p676d11_ga(f, 1015, 3, 300);
gamma_g = gamma_o + gamma_w;
cyl = 1;
theta = 10;
theta_path = 10;
R_rain = 1;
pol = 0; 

pr_scatt = p.test_precipitation_scatter(Pint, f, r1, r2, Ptx, G, dtx, drx, gamma_g, cyl, theta, theta_path, R_rain, pol);


