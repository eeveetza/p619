function M = p453_nwet_maps()
%P453_NWET_MAPS Annual maps of Nwet of Rec. ITU-R P.453-14 (loaded once per session)
%
%   M.p     1 x 18 time percentages (0.1 ... 99 %)
%   M.nwet  241 x 481 x 18 Nwet (ppm); rows: latitude -90 ... 90 deg,
%           columns: longitude -180 ... 180 deg (0.75 deg grid)
%
%   Sources, in this order:
%   1) private/DigitalMaps_P453.mat, written by initiate_digital_maps.m from
%      the files NWET_Annual_<p>.TXT of Rec. ITU-R P.453-14
%   2) the functions private/DigitalMaps_NWET_Annual_<p>.m

persistent cache
if ~isempty(cache)
    M = cache;
    return
end

here = fileparts(mfilename('fullpath'));
ptags = {'01', '02', '03', '05', '1', '2', '3', '5', '10', '20', '30', '50', '60', '70', '80', '90', '95', '99'};
pvals = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 60 70 80 90 95 99];

matfile = fullfile(here, 'DigitalMaps_P453.mat');
if exist(matfile, 'file') == 2
    S = load(matfile);
    M = S.P453;
elseif exist(fullfile(here, 'DigitalMaps_NWET_Annual_50.m'), 'file') == 2
    M.p = pvals;
    M.nwet = zeros(241, 481, 18);
    for k = 1:18
        M.nwet(:, :, k) = feval(['DigitalMaps_NWET_Annual_' ptags{k}]);
    end
else
    error(['P.453-14 Nwet maps not found.\n' ...
        'Download the Rec. ITU-R P.453-14 digital maps (NWET_Annual_<p>.TXT) from\n' ...
        'https://www.itu.int/rec/R-REC-P.453-14-201908-I/en, place them in ./private/maps\n' ...
        'and run initiate_digital_maps.m.'], '');
end
cache = M;
end
