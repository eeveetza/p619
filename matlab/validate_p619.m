function results = validate_p619(vectorDir)
%VALIDATE_P619 Runs the P.619-6 validation vectors against P619.m
%
%   results = validate_p619()            uses ../validation/vectors relative to this file
%   results = validate_p619(vectorDir)
%
%   Each vector file is evaluated by an adapter (one local function per
%   file, below). To validate another implementation, only the adapters
%   need to be rewritten: each returns a struct whose field names are the
%   output columns of the CSV file, with one value per row (NaN = not
%   provided by the implementation, reported as "n/a").
%
%   Pass criterion per value:  |impl - ref| <= abs_tol + rel_tol*|ref|
%   with tolerances from tolerances.csv.
%
%   Works in MATLAB (R2019b+) and Octave (6+).

if nargin < 1
    here = fileparts(mfilename('fullpath'));
    vectorDir = fullfile(here, '.', 'validation', 'vectors');
end

p = P619;
tol = read_tolerances(fullfile(vectorDir, 'tolerances.csv'));

files = { ...
    '01_free_space.csv',            @run_free_space;
    '02_depolarisation.csv',        @run_depolarisation;
    '03_faraday.csv',               @run_faraday;
    '04_beam_spreading.csv',        @run_beam_spreading;
    '05_geometry_attA.csv',         @run_geometry;
    '06_attB_fs_to_apparent.csv',   @run_fs2app;
    '07_attB_apparent_to_fs.csv',   @run_app2fs;
    '08_specific_attenuation.csv',  @run_specific_attenuation;
    '09_atmosphere_magra.csv',      @run_magra;
    '10_ag_ascending.csv',          @run_ag_ascending;
    '11_ag_descending.csv',         @run_ag_descending;
    '12_hmin.csv',                  @run_hmin;
    '13_ag_negative_elevation.csv', @run_ag_negative;
    '14_scintillation_attD.csv',    @run_scintillation;
    '15_ray_height_attE.csv',       @run_ray_height;
    '16_ducting_ldtb.csv',          @run_ldtb;
    '17_end_to_end_eq14.csv',       @run_end_to_end;
    '18_precip_scatter_attF.csv',   @run_precip_scatter;
    '19_fresnel_knife_edge.csv',    @run_fresnel_knife_edge;
    '20_finite_width_screen.csv',   @run_finite_width_screen;
    '21_terrain_ldtb.csv',          @run_terrain_ldtb;
    '22_terrain_slant_p526.csv',    @run_terrain_slant;
    '23_bel_p2109.csv',             @run_bel;
    '24_multiple_entry_eq15.csv',   @run_multiple_entry;
    '25_clutter_p2108_2_draft.csv', @run_clutter_draft;
    '26_p836_rho.csv',              @run_p836;
    '27_gas_ag_p1.csv',             @run_gas_ag_p1;
    '28_nwet_p453.csv',             @run_nwet;
    '29_end_to_end_scint.csv',      @run_end_to_end_scint };

results = struct('file', {}, 'quantity', {}, 'n', {}, 'passed', {}, 'na', {}, 'max_abs_dev', {}, 'worst_id', {});
fprintf('\nP.619-6 validation of P619.m\n%s\n', repmat('-', 1, 92));
fprintf('%-30s %-16s %8s %8s %14s %8s\n', 'file', 'quantity', 'passed', 'n/a', 'max |dev|', 'worst id');

for k = 1:size(files, 1)
    [ref, txt] = read_vectors(fullfile(vectorDir, files{k, 1}));
    out = files{k, 2}(p, ref, txt);
    qnames = fieldnames(out);
    for q = 1:numel(qnames)
        qn = qnames{q};
        [at, rt] = lookup_tol(tol, files{k, 1}, qn);
        r = ref.(qn);
        v = out.(qn);
        valid = ~isnan(r);
        na = valid & isnan(v);
        dev = abs(v - r);
        ok = valid & ~na & (dev <= at + rt .* abs(r));
        dev(~valid | na) = NaN;
        [mx, iw] = max(dev);
        if isnan(mx)
            worst = '-';
        else
            worst = sprintf('%d', ref.id(iw));
        end
        res = struct('file', files{k, 1}, 'quantity', qn, 'n', nnz(valid), 'passed', nnz(ok), ...
                     'na', nnz(na), 'max_abs_dev', mx, 'worst_id', worst);
        results(end + 1) = res; %#ok<AGROW>
        flag = '';
        if nnz(ok) < nnz(valid & ~na)
            flag = '  <-- FAIL';
        end
        fprintf('%-30s %-16s %4d/%-3d %8d %14.3g %8s%s\n', files{k, 1}, qn, nnz(ok), nnz(valid), nnz(na), mx, worst, flag);
    end
end

nPass = sum([results.passed]);
nTot = sum([results.n]);
nNA = sum([results.na]);
fprintf('%s\nTotal: %d of %d values passed (%d not provided by the implementation)\n\n', repmat('-', 1, 92), nPass, nTot, nNA);

end

% =========================================================================
%                       Adapters for P619.m
% =========================================================================

function out = run_free_space(p, r, ~)
out.Lbfs_dB = arrayfun(@(f, d) p.tl_free_space(f, d), r.f_GHz, r.d_km);
end

function out = run_depolarisation(p, r, ~)
[out.Ax_dB, out.Ac_dB] = p.co_and_cross_polar_attenuation(r.XPD_dB);
out.Axq_dB = arrayfun(@(x) p.hydrometeor_depolarization_attenuation(x), r.XPD_dB);
end

function out = run_faraday(p, r, ~)
[out.AxF_dB, out.AcF_dB, out.thetaF_rad] = arrayfun(@(f, b, nt) p.faraday_co_and_cross_polar_attenuation(f, b, nt), ...
    r.f_GHz, r.Bav_T, r.NT_el_per_m2);
end

function out = run_beam_spreading(p, r, ~)
out.Abs_dB = arrayfun(@(t, h) p.beam_spreading_loss(t, h), r.theta0_deg, r.h_km);
end

function out = run_geometry(p, r, ~)
n = numel(r.id);
out.Dts_km = NaN(n, 1); out.theta0_deg = NaN(n, 1); out.psi_deg = NaN(n, 1);
for i = 1:n
    [out.Dts_km(i), out.theta0_deg(i), out.psi_deg(i)] = ...
        p.straight_line_e2s(r.Hs_km(i), r.Ht_km(i), r.lat_s_deg(i), r.lat_t_deg(i), r.delta_deg(i));
end
end

function out = run_fs2app(p, r, ~)
out.theta_deg = arrayfun(@(h, t) p.fs2apparent(t, h), r.Ht_km, r.theta0_deg);
end

function out = run_app2fs(p, r, ~)
out.theta0_deg = arrayfun(@(h, t) p.apparent2fs(t, h), r.Ht_km, r.theta_deg);
end

function out = run_specific_attenuation(p, r, ~)
n = numel(r.id);
out.gamma_o_dBkm = NaN(n, 1); out.gamma_w_dBkm = NaN(n, 1);
for i = 1:n
    [out.gamma_o_dBkm(i), out.gamma_w_dBkm(i)] = p.p676d11_ga(r.f_GHz(i), r.pdry_hPa(i), r.rho_gm3(i), r.T_K(i));
end
end

function out = run_magra(p, r, ~)
[T, P, rho, n] = p.p835_std_atm_profiles(r.h_km, 1);
out.T_K = T(:); out.P_hPa = P(:); out.rho_gm3 = rho(:);
out.e_hPa = out.rho_gm3 .* out.T_K / 216.7;
out.n = n(:);
end

function Ag = ag_case(p, f, He, Hs, phi_e, phi_s, descending)
hl = p.p676_slant_path_geometry16(max(He, 0), min(Hs, 100));
hm = 0.5 * (hl(1:end-1) + hl(2:end));
[T, P, rho, n] = p.p835_std_atm_profiles(hm, 1);
if descending
    Ag = p.atm_attenuation_s2E(f, He, Hs, phi_e, phi_s, 360, hl, rho, T, P, n, true);
else
    Ag = p.atm_attenuation_E2s(f, He, Hs, phi_e, 0, 360, hl, rho, T, P, n, true, 1);
end
end

function out = run_ag_ascending(p, r, ~)
out.Ag_dB = arrayfun(@(f, he, hs, ph) ag_case(p, f, he, hs, ph, 0, false), r.f_GHz, r.He_km, r.Hs_km, r.phi_e_deg);
end

function out = run_ag_descending(p, r, ~)
out.Ag_dB = arrayfun(@(f, he, hs, pe, ps) ag_case(p, f, he, hs, pe, ps, true), ...
    r.f_GHz, r.He_km, r.Hs_km, r.phi_e_arrival_deg, r.phi_s_deg);
end

function out = run_hmin(p, r, ~)
out.Hmin_km = arrayfun(@(he, ph) p.solve45(he, ph, 1), r.He_km, r.phi_e_deg);
end

function out = run_ag_negative(p, r, ~)
out.Hmin_km = arrayfun(@(he, ph) p.solve45(he, ph, 1), r.He_km, r.phi_e_deg);
out.Ag_dB = arrayfun(@(f, he, hs, ph) ag_case(p, f, he, hs, ph, 0, false), r.f_GHz, r.He_km, r.Hs_km, r.phi_e_deg);
end

function out = run_scintillation(p, r, ~)
[out.Ast_dB, out.sigma_dB] = arrayfun(@(f, pp, nw, th, ga) p.tropospheric_scintillation(f, pp, nw, th, ga), ...
    r.f_GHz, r.p_pct, r.Nwet, r.theta_deg, r.Ga_dBi);
end

function out = run_ray_height(p, r, ~)
n = numel(r.id);
out.Hr_km = NaN(n, 1);
for i = 1:n
    [Dc, Hr] = p.beam_clearance(r.Ht_km(i), r.theta_deg(i));
    k = find(abs(Dc - r.Dc_km(i)) < 1e-9, 1);
    if ~isempty(k)
        out.Hr_km(i) = Hr(k);
    end
end
end

function out = run_ldtb(p, r, ~)
n = numel(r.id);
out.Ldtb_dB = NaN(n, 1);
for i = 1:n
    [~, out.Ldtb_dB(i)] = evalc_quiet(@() p.tl_dtb(r.f_GHz(i), r.p_pct(i), r.lat_deg(i), ...
        r.theta_hoz_mrad(i), r.d_hoz_km(i), r.Ld_dB(i)));
end
end

function out = run_end_to_end(p, r, txt)
n = numel(r.id);
names = {'Dts_km', 'theta0_deg', 'theta_deg', 'Lbfs_dB', 'Ag_dB', 'Abs_dB', 'As_dB', 'Ldtb_dB', 'Lb_dB'};
for k = 1:numel(names)
    out.(names{k}) = NaN(n, 1);
end
for i = 1:n
    e2s = strcmp(txt.direction{i}, 'E2s');
    % antennas pointing along the path ([] beams), MAGRA (atm_type 1),
    % p1 = p2 = 50 %, Ga = 40 dBi (only used by scintillation, As ~ 0 at 50 %)
    L = p.tl_p619_single(r.f_GHz(i), [], e2s, 50, r.He_km(i), r.Hs_km(i), [], [], 360, 1, ...
        r.lat_s_deg(i), r.lat_e_deg(i), r.lon_s_deg(i), r.lon_e_deg(i), r.p2_pct(i), 40);
    out.Dts_km(i) = L.dkm;
    out.theta0_deg(i) = L.theta_0;
    out.theta_deg(i) = L.theta;
    out.Lbfs_dB(i) = L.Lbfs;
    out.Ag_dB(i) = L.Ag;
    out.Abs_dB(i) = L.Abs;
    out.As_dB(i) = L.Ast;
    out.Ldtb_dB(i) = 0;                          % no terrain profile given
    % P619.m has no default Axp (xpdMode 'none' gives Ax = 0); the vectors
    % use the §2.2 default Axp = 3 dB, so it is added here.
    out.Lb_dB(i) = L.Lb - L.Ax + r.Axp_dB(i);
end
end

function out = run_precip_scatter(p, r, ~)
n = numel(r.id);
names = {'L1_m', 'L2_m', 'S_dBWm2', 'Pin_dBW', 'gamma_r_dBkm', 'Pout_dBW', 'Pscat_dBW', ...
         'Peirps_dBW', 'Fnis_dB', 'Ptxs_dBW', 'full_calc'};
fields = {'L1', 'L2', 'S', 'Pin', 'gamma_r', 'Pout', 'Pscat', 'Peirps', 'Fnis', 'Ptxs'};
for k = 1:numel(names)
    out.(names{k}) = NaN(n, 1);
end
pol_of_tau = containers.Map({0, 90, 45}, {0, 1, 2});   % P619.m: 0 = h, 1 = v, 2 = circular
for i = 1:n
    [full, o] = p.test_precipitation_scatter(r.Pint_dBW(i), r.f_GHz(i), r.r1_m(i), r.r2_m(i), ...
        r.Ptx_dBW(i), r.G_dBi(i), r.dtx_km(i), r.drx_km(i), r.gamma_g_dBkm(i), r.cyl(i), ...
        r.theta_deg(i), r.theta_path_deg(i), r.R_mmh(i), pol_of_tau(r.tau_deg(i)));
    for k = 1:numel(fields)
        out.(names{k})(i) = o.(fields{k});
    end
    out.full_calc(i) = double(full);
end
end

function out = run_fresnel_knife_edge(p, r, ~)
[out.nu, out.R1_m] = p.fresnel_parameters(r.h_m, r.d1_km, r.f_GHz);
out.J_dB = p.knife_edge_loss(out.nu);
end

function out = run_finite_width_screen(p, r, ~)
[out.Jmin_dB, out.Jav_dB] = arrayfun(@(a, b, c) p.finite_width_screen_loss(a, b, c), r.nu_top, r.nu_left, r.nu_right);
end

function out = run_terrain_ldtb(p, r, ~)
n = numel(r.id);
names = {'nu_max', 'd_max_km', 'h_max_m', 'Ld_dB', 'theta_hoz_mrad', 'd_hoz_km', 'Ldtb_dB'};
fields = {'nu_max', 'd_max', 'h_max', 'Ld', 'th_hoz', 'd_hoz'};
for k = 1:numel(names)
    out.(names{k}) = NaN(n, 1);
end
for i = 1:n
    % synthetic ridge profile, see validation/README.md
    d = (0:round(r.dmax_km(i) / r.step_km(i)))' * r.step_km(i);
    h = r.h0_m(i) + r.A_m(i) * exp(-((d - r.dc_km(i)) / r.w_km(i)).^2);
    [Ldtb, o] = p.diffraction_ducting_loss(r.f_GHz(i), r.p_pct(i), r.lat_deg(i), r.Ht_km(i), ...
        r.theta_deg(i), d, h, r.ae_km(i), 'knife-edge');
    for k = 1:numel(fields)
        out.(names{k})(i) = o.(fields{k});
    end
    out.Ldtb_dB(i) = Ldtb;
end
end

function out = run_terrain_slant(p, r, ~)
n = numel(r.id);
names = {'nu', 'h_m', 'd1_km', 'J_dB', 'Le_dB', 'kb', 'Ld_dB', 'theta_hoz_mrad', 'd_hoz_km', 'Ldtb_dB'};
fields = {'nu_max', 'h_max', 'd_max', 'J', 'Le', 'kb', 'Ld', 'th_hoz', 'd_hoz'};
for k = 1:numel(names)
    out.(names{k}) = NaN(n, 1);
end
for i = 1:n
    % profile with two ridges, see validation/README.md
    d = (0:round(r.dmax_km(i) / r.step_km(i)))' * r.step_km(i);
    h = r.h0_m(i) + r.A1_m(i) * exp(-((d - r.dc1_km(i)) / r.w1_km(i)).^2) ...
                  + r.A2_m(i) * exp(-((d - r.dc2_km(i)) / r.w2_km(i)).^2);
    [Ldtb, o] = p.diffraction_ducting_loss(r.f_GHz(i), r.p_pct(i), r.lat_deg(i), r.Ht_km(i), ...
        r.theta_deg(i), d, h, r.ae_km(i), 'slant');
    for k = 1:numel(fields)
        out.(names{k})(i) = o.(fields{k});
    end
    out.Ldtb_dB(i) = Ldtb;
end
end

function out = run_bel(p, r, ~)
out.Lbe_dB = arrayfun(@(f, pp, b, th) p.bel_p2109(f, pp, b, th), r.f_GHz, r.p_pct, r.building, r.theta_deg);
end

function out = run_multiple_entry(p, r, txt)
n = numel(r.id);
names = {'Dts_km', 'theta0_deg', 'theta_deg', 'Lbfs_dB', 'Ag_dB', 'Abs_dB', 'Lc_dB', 'Lbe_dB', 'Ldtb_dB', 'Lb_dB'};
fields = {'dkm', 'theta_0', 'theta', 'Lbfs', 'Ag', 'Abs', 'Lc', 'Lbe', 'Ldtb', 'Lb'};
for k = 1:numel(names)
    out.(names{k}) = NaN(n, 1);
end
for i = 1:n
    args = {'Axp', r.Axp_dB(i), 'clutter', r.clutter(i) == 1, 'indoor', r.indoor(i) == 1};
    if r.clutter(i) == 1
        % default clutter model: draft P.2108-2 §3.3
        args = [args, {'pLc', r.pLc_pct(i), 'hs', r.hs_m(i), 'hm', r.hm_m(i)}];
    end
    if r.indoor(i) == 1
        args = [args, {'pLbe', r.pLbe_pct(i), 'buildingType', r.building(i)}];
    end
    L = p.tl_p619_multiple(r.f_GHz(i), [], true, r.p_pct(i), r.He_km(i), r.Hs_km(i), [], [], 360, 1, ...
        r.lat_s_deg(i), r.lat_e_deg(i), r.lon_s_deg(i), r.lon_e_deg(i), args{:});
    for k = 1:numel(fields)
        out.(names{k})(i) = L.(fields{k});
    end
end
end

function out = run_clutter_draft(p, r, ~)
out.Lc_dB = arrayfun(@(f, th, pp, h, hm) p.cl_p2108_3(f, th, pp, h, hm), r.f_GHz, r.theta_deg, r.p_pct, r.h_m, r.hm_m);
end

function out = run_p836(p, r, ~)
% Requires the P.836-6 maps (initiate_digital_maps.m); n/a otherwise
out.rho_gm3 = NaN(size(r.id));
try
    for i = 1:numel(r.id)
        out.rho_gm3(i) = p.p836_rho_surface(r.p_pct(i), r.lat_deg(i), r.lon_deg(i), r.alt_km(i));
    end
catch err
    fprintf(1, '   (26_p836_rho.csv: %s)\n', strtok(err.message, sprintf('\n')));
end
end

function out = run_gas_ag_p1(p, r, ~)
n = numel(r.id);
names = {'q_pct', 'rho_s_gm3', 'rho0_gm3', 'theta_deg', 'Ag_dB'};
fields = {'q_rho', 'rho_s', 'rho0', 'theta', 'Ag'};
for k = 1:numel(names)
    out.(names{k}) = NaN(n, 1);
end
for i = 1:n
    if r.mode(i) == 1
        args = {'rho0', r.rho0_in_gm3(i)};
    else
        args = {'waterVapour', 'P.836'};
    end
    try
        L = p.tl_p619_single(r.f_GHz(i), [], true, r.p1_pct(i), r.He_km(i), r.Hs_km(i), [], [], 360, 1, ...
            r.lat_s_deg(i), r.lat_e_deg(i), r.lon_s_deg(i), r.lon_e_deg(i), 50, 40, args{:});
    catch err      % P.836 maps not available
        continue
    end
    for k = 1:numel(fields)
        out.(names{k})(i) = L.(fields{k});
    end
end
end

function out = run_nwet(p, r, ~)
% Requires the P.453-14 Nwet maps (DigitalMaps_NWET_Annual_*.m or
% initiate_digital_maps.m); n/a otherwise
out.Nwet = NaN(size(r.id));
try
    for i = 1:numel(r.id)
        out.Nwet(i) = p.get_interp2_Nwet_Annual_time_location(r.p_pct(i), r.lon_deg(i), r.lat_deg(i));
    end
catch err
    fprintf(1, '   (28_nwet_p453.csv: %s)\n', strtok(err.message, sprintf('\n')));
end
end

function out = run_end_to_end_scint(p, r, ~)
n = numel(r.id);
names = {'theta0_deg', 'Nwet', 'As_dB', 'Lb_dB'};
for k = 1:numel(names)
    out.(names{k}) = NaN(n, 1);
end
for i = 1:n
    % antennas pointing along the path ([] beams), MAGRA (atm_type 1), p1 = 50 %
    try
        L = p.tl_p619_single(r.f_GHz(i), [], true, 50, r.He_km(i), r.Hs_km(i), [], [], 360, 1, ...
            r.lat_s_deg(i), r.lat_e_deg(i), r.lon_s_deg(i), r.lon_e_deg(i), r.p2_pct(i), r.Ga_dBi(i));
    catch err      % Nwet maps not available
        continue
    end
    out.theta0_deg(i) = L.theta_0;
    out.Nwet(i) = L.Nwet;
    out.As_dB(i) = L.Ast;
    % no default Axp in P619.m ('none' gives Ax = 0); the vectors use Axp = 3 dB
    out.Lb_dB(i) = L.Lb - L.Ax + r.Axp_dB(i);
end
end

% =========================================================================
%                              Helpers
% =========================================================================

function [txtout, val] = evalc_quiet(fh)
% Calls fh() while suppressing printed output (tl_dtb echoes dn).
val = [];
txtout = evalc('val = fh();');
end

function [num, txt] = read_vectors(fname)
% Reads a vector CSV: numeric columns into num.<name>, text columns
% (source, scenario, direction) into txt.<name> as cellstr.
fid = fopen(fname, 'r');
if fid < 0
    error('Cannot open %s', fname);
end
header = strsplit(strtrim(fgetl(fid)), ',');
lines = {};
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    if ~isempty(strtrim(l)), lines{end + 1} = l; end %#ok<AGROW>
end
fclose(fid);
nc = numel(header);
nr = numel(lines);
cells = cell(nr, nc);
for i = 1:nr
    cells(i, :) = split_csv_line(lines{i}, nc);
end
num = struct();
txt = struct();
for j = 1:nc
    col = cells(:, j);
    v = str2double(col);
    isnum = ~isnan(v) | strcmpi(col, 'nan');
    if all(isnum)
        num.(header{j}) = v;
    else
        txt.(header{j}) = col;
    end
end
end

function c = split_csv_line(l, nc)
% Minimal RFC-4180 split (quoted fields may contain commas).
c = cell(1, nc);
k = 1; buf = ''; inq = false;
for ch = l
    if ch == '"'
        inq = ~inq;
    elseif ch == ',' && ~inq
        c{k} = buf; k = k + 1; buf = '';
    else
        buf(end + 1) = ch; %#ok<AGROW>
    end
end
c{k} = buf;
end

function tol = read_tolerances(fname)
[~, t] = read_vectors(fname);
[n, ~] = read_vectors(fname);
tol.file = t.file;
tol.quantity = t.quantity;
tol.abs = n.abs_tol;
tol.rel = n.rel_tol;
end

function [at, rt] = lookup_tol(tol, file, q)
at = 0; rt = 0;
for i = 1:numel(tol.file)
    if strcmp(tol.file{i}, file) && any(strcmp(strsplit(tol.quantity{i}, ';'), q))
        at = tol.abs(i); rt = tol.rel(i);
        return
    end
end
error('No tolerance for %s / %s', file, q);
end
