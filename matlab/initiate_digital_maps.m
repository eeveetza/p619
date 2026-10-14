% This Recommendation uses integral digital products. They form an integral
% part of the recommendation and may not be reproduced, by any means whatsoever,
% without written permission of ITU.
%
% The user should download the necessary digital maps that are used by this
% implementation directly from ITU-R web site and place the files in the
% folder ./private/maps. After that, the user should execute the script
% contained in this file "initiate_digital_maps.m". The script produces the
% necessary functions/files that contain the digital maps and are used for
% interpolations. All maps are optional; the corresponding functionality
% is only available once the maps have been processed.
%
% 1) Rec. ITU-R P.452-18 (median effective Earth radius, Section 2.6, when
%    'ae' is not given)
%    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.452-18-202310-I!!ZIP-E.zip
%       - DN50.TXT
%       --> private/DigitalMaps_DN50.m
%
% 2) Rec. ITU-R P.836-6 (surface water-vapour density, time-percentage
%    dependence of Ag(p1), 'waterVapour', 'P.836')
%    https://www.itu.int/rec/R-REC-P.836-6-201712-I/en
%       - the 18 annual maps of surface water-vapour density and of the
%         water-vapour scale height for p = 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5,
%         10, 20, 30, 50, 60, 70, 80, 90, 95 and 99 %, whose file names
%         start with RHO_<p> and VSCH_<p> (<p> = 01, 02, 03, 05, 1, 2, ... 99)
%       - TOPO_0DOT5 (topographic altitude, 0.5 deg grid)
%       The files may be placed in sub-folders of ./private/maps; file names
%       are matched case-insensitively.
%       --> private/DigitalMaps_P836.mat
%
% 3) Rec. ITU-R P.453-14 (wet term of the surface refractivity Nwet,
%    tropospheric scintillation, Attachment D)
%    https://www.itu.int/rec/R-REC-P.453-14-201908-I/en
%       - NWET_Annual_<p>.TXT, <p> = 01, 02, 03, 05, 1, 2, ... 99 (18 maps,
%         0.75 deg grid, 241 x 481)
%       --> private/DigitalMaps_P453.mat
%       Not needed if the functions private/DigitalMaps_NWET_Annual_<p>.m
%       are present (they are used when DigitalMaps_P453.mat is absent).

here = fileparts(mfilename('fullpath'));
if isempty(here)
    here = pwd;
end
mapdir = fullfile(here, 'private', 'maps');

%% 1) P.452-18 DN50 ------------------------------------------------------------

files = {'DN50.TXT'};

for i  = 1 : length(files)

    xx = files{i};
    filename = fullfile(mapdir, xx);
    functionname = ['DigitalMaps_' xx(1:end-4)];
    filename_m = fullfile(here, 'private', [functionname '.m']);

    if (exist(filename_m, 'file') == 2)
        continue
    end

    fid = fopen(filename, "r");

    if (fid == -1)
        warning(['%s not found in ./private/maps (from ITU-R P.452-18). ' ...
                 'It is only needed for Section 2.6 when ''ae'' is not given.'], xx);
        continue
    end

    fprintf(1, "Processing file: %s --> %s...\n", filename, filename_m);

    datastr = [];

    while (1)

        line = fgetl(fid);

        if (line == -1)
            break;
        end

        datastr = [datastr line '\n'];

    end

    fclose(fid);

    % write the data in the function.m

    fid = fopen(filename_m, 'w');

    datastr_before = ['function y = ' functionname '()\n' ...
        '%%%% ' functionname ' returns values from the table ' files{i} ,'\n' ...
        '\n' ...
        '\n' ...
        'y = [' ];

    datastr_after = [ '];\n' ...
        'return\n' ...
        'end\n'];

    fprintf(fid, datastr_before);
    fprintf(fid, datastr);
    fprintf(fid, datastr_after);
    fclose(fid);
end

%% 2) P.836-6 annual maps --------------------------------------------------------

matfile = fullfile(here, 'private', 'DigitalMaps_P836.mat');

if (exist(matfile, 'file') ~= 2)

    % list all files below ./private/maps (recursively)
    names = {};
    paths = {};
    queue = {mapdir};
    while ~isempty(queue)
        dd = queue{1};
        queue(1) = [];
        entries = dir(dd);
        for k = 1:numel(entries)
            e = entries(k);
            if any(strcmp(e.name, {'.', '..'}))
                continue
            end
            if e.isdir
                queue{end+1} = fullfile(dd, e.name); %#ok<SAGROW>
            else
                names{end+1} = e.name; %#ok<SAGROW>
                paths{end+1} = fullfile(dd, e.name); %#ok<SAGROW>
            end
        end
    end

    ptags = {'01', '02', '03', '05', '1', '2', '3', '5', '10', '20', '30', '50', '60', '70', '80', '90', '95', '99'};
    pvals = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 60 70 80 90 95 99];

    find_file = @(pattern) find(~cellfun(@isempty, regexpi(names, pattern, 'once')));

    itopo = find_file('^TOPO_0DOT5');
    missing = {};
    if isempty(itopo)
        missing{end+1} = 'TOPO_0DOT5';
    end
    irho = zeros(1, 18);
    ivsch = zeros(1, 18);
    for k = 1:18
        a = find_file(['^RHO_' ptags{k} '([_.]|$)']);
        b = find_file(['^VSCH_' ptags{k} '([_.]|$)']);
        if isempty(a), missing{end+1} = ['RHO_' ptags{k}]; else, irho(k) = a(1); end %#ok<SAGROW>
        if isempty(b), missing{end+1} = ['VSCH_' ptags{k}]; else, ivsch(k) = b(1); end %#ok<SAGROW>
    end

    if ~isempty(missing)
        warning(['P.836-6 maps not found in ./private/maps (missing: %s). ' ...
                 'They are only needed for ''waterVapour'', ''P.836''.'], strjoin(missing, ', '));
    else
        fprintf(1, 'Processing P.836-6 annual maps --> %s...\n', matfile);
        P836.p = pvals;
        P836.rho = zeros(161, 321, 18);
        P836.vsch = zeros(161, 321, 18);
        for k = 1:18
            A = load(paths{irho(k)});
            B = load(paths{ivsch(k)});
            if ~isequal(size(A), [161 321]) || ~isequal(size(B), [161 321])
                error('Unexpected size of %s or %s (expected 161 x 321, 1.125 deg grid).', names{irho(k)}, names{ivsch(k)});
            end
            fprintf(1, '   %s, %s\n', names{irho(k)}, names{ivsch(k)});
            P836.rho(:, :, k) = A;
            P836.vsch(:, :, k) = B;
        end
        T = load(paths{itopo(1)});
        if ~isequal(size(T), [363 723])
            error('Unexpected size of %s (expected 363 x 723, 0.5 deg grid with border).', names{itopo(1)});
        end
        fprintf(1, '   %s\n', names{itopo(1)});
        P836.topo = T;
        save(matfile, 'P836', '-v7');
    end
end

%% 3) P.453-14 annual Nwet maps ---------------------------------------------------

matfile = fullfile(here, 'private', 'DigitalMaps_P453.mat');

if (exist(matfile, 'file') ~= 2)

    names = {};
    paths = {};
    queue = {mapdir};
    while ~isempty(queue)
        dd = queue{1};
        queue(1) = [];
        entries = dir(dd);
        for k = 1:numel(entries)
            e = entries(k);
            if any(strcmp(e.name, {'.', '..'}))
                continue
            end
            if e.isdir
                queue{end+1} = fullfile(dd, e.name); %#ok<SAGROW>
            else
                names{end+1} = e.name; %#ok<SAGROW>
                paths{end+1} = fullfile(dd, e.name); %#ok<SAGROW>
            end
        end
    end

    ptags = {'01', '02', '03', '05', '1', '2', '3', '5', '10', '20', '30', '50', '60', '70', '80', '90', '95', '99'};
    pvals = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 60 70 80 90 95 99];

    inwet = zeros(1, 18);
    missing = {};
    for k = 1:18
        a = find(~cellfun(@isempty, regexpi(names, ['^NWET_Annual_' ptags{k} '\.TXT$'], 'once')));
        if isempty(a), missing{end+1} = ['NWET_Annual_' ptags{k}]; else, inwet(k) = a(1); end %#ok<SAGROW>
    end

    if ~isempty(missing)
        if (exist(fullfile(here, 'private', 'DigitalMaps_NWET_Annual_50.m'), 'file') ~= 2)
            warning(['P.453-14 Nwet maps not found in ./private/maps (missing: %s). ' ...
                     'They are needed for tropospheric scintillation unless ''Nwet'' is given.'], strjoin(missing, ', '));
        end
    else
        fprintf(1, 'Processing P.453-14 annual Nwet maps --> %s...\n', matfile);
        P453.p = pvals;
        P453.nwet = zeros(241, 481, 18);
        for k = 1:18
            A = load(paths{inwet(k)});
            if ~isequal(size(A), [241 481])
                error('Unexpected size of %s (expected 241 x 481, 0.75 deg grid).', names{inwet(k)});
            end
            fprintf(1, '   %s\n', names{inwet(k)});
            P453.nwet(:, :, k) = A;
        end
        save(matfile, 'P453', '-v7');
    end
end
