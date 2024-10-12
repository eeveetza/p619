%     Class implementing Recommendation ITU-R P.619-5
%     Recommendation ITU-R P.619 provides methods for predicting signal 
%     propagation losses for interfering signals between stations in space 
%     and stations on (or near to) the surface of the Earth in the overall 
%     frequency range of 100 MHz to 100 GHz, except for a few exceptions 
%     restricted to lower frequencies which will be specified where they 
%     are described. This Recommendation provides methods to predict 
%     the propagation losses not exceeded for 0.001%-50% of the time. 
%     Guidance is given for single entry as well as multiple entry 
%     propagation losses in analyses that determine interfering signals, 
%     where correlations of temporal variability and location variability 
%     may be influential.
%
%     This implementation is not yet complete and is work under progress
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    05JUL24     Ivica Stevanovic, OFCOM         Initial version

classdef P619

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                               Section 2.1                               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Lbfs = tl_free_space(obj, fGHz, dkm)
            %%tl_free_space Computes basic transmission loss in free space 
            %
            % Recommendation ITU-R P.619-5 §2.1
            %
            
            % Inputs
            % Variable    Unit     Type     Description
            % fGHz	      GHz	   float    Frequency
            % dkm         km       float    Path length
            %
            % Outputs:
            %
            % Lbfs        dB       float    Free space basic transmission loss


            Lbfs = 92.45 + 20*log10(fGHz*dkm);                        % (1)


        end


        function [g_0, g_w] = p676d11_ga(obj, f, p, rho, T)
            %p676d11_ga Specific attenuation due to dry air and water vapour
            % [g_0, g_w] = p676d11_ga(f, p, rho, T)
            % This function computes the specific attenuation due to dry air and water vapour,
            % at frequencies up to 1 000 GHz for different values of of pressure, temperature
            % and humidity by means of a summation of the individual resonance lines from
            % oxygen and water vapour according to ITU-R P.676-11
            %
            % Input parameters:
            % f       -   Frequency (GHz)
            % p       -   Dry air pressure (hPa)
            % rho     -   Water vapor density (g/m^3)
            % T       -   Temperature (K)
            %
            % Output parameters:
            % g_o, g_w   -   specific attenuation due to dry air and water vapour


            %% spectroscopic data for oxigen
            %             f0        a1    a2     a3   a4     a5     a6
            oxigen = [50.474214, 0.975, 9.651, 6.690, 0.0, 2.566, 6.850;
                50.987745, 2.529, 8.653, 7.170, 0.0, 2.246, 6.800;
                51.503360, 6.193, 7.709, 7.640, 0.0, 1.947, 6.729;
                52.021429, 14.320, 6.819, 8.110, 0.0, 1.667, 6.640;
                52.542418, 31.240, 5.983, 8.580, 0.0, 1.388, 6.526;
                53.066934, 64.290, 5.201, 9.060, 0.0, 1.349, 6.206;
                53.595775, 124.600, 4.474, 9.550, 0.0, 2.227, 5.085;
                54.130025, 227.300, 3.800, 9.960, 0.0, 3.170, 3.750;
                54.671180, 389.700, 3.182, 10.370, 0.0, 3.558, 2.654;
                55.221384, 627.100, 2.618, 10.890, 0.0, 2.560, 2.952;
                55.783815, 945.300, 2.109, 11.340, 0.0, -1.172, 6.135;
                56.264774, 543.400, 0.014, 17.030, 0.0, 3.525, -0.978;
                56.363399, 1331.800, 1.654, 11.890, 0.0, -2.378, 6.547;
                56.968211, 1746.600, 1.255, 12.230, 0.0, -3.545, 6.451;
                57.612486, 2120.100, 0.910, 12.620, 0.0, -5.416, 6.056;
                58.323877, 2363.700, 0.621, 12.950, 0.0, -1.932, 0.436;
                58.446588, 1442.100, 0.083, 14.910, 0.0, 6.768, -1.273;
                59.164204, 2379.900, 0.387, 13.530, 0.0, -6.561, 2.309;
                59.590983, 2090.700, 0.207, 14.080, 0.0, 6.957, -0.776;
                60.306056, 2103.400, 0.207, 14.150, 0.0, -6.395, 0.699;
                60.434778, 2438.000, 0.386, 13.390, 0.0, 6.342, -2.825;
                61.150562, 2479.500, 0.621, 12.920, 0.0, 1.014, -0.584;
                61.800158, 2275.900, 0.910, 12.630, 0.0, 5.014, -6.619;
                62.411220, 1915.400, 1.255, 12.170, 0.0, 3.029, -6.759;
                62.486253, 1503.000, 0.083, 15.130, 0.0, -4.499, 0.844;
                62.997984, 1490.200, 1.654, 11.740, 0.0, 1.856, -6.675;
                63.568526, 1078.000, 2.108, 11.340, 0.0, 0.658, -6.139;
                64.127775, 728.700, 2.617, 10.880, 0.0, -3.036, -2.895;
                64.678910, 461.300, 3.181, 10.380, 0.0, -3.968, -2.590;
                65.224078, 274.000, 3.800, 9.960, 0.0, -3.528, -3.680;
                65.764779, 153.000, 4.473, 9.550, 0.0, -2.548, -5.002;
                66.302096, 80.400, 5.200, 9.060, 0.0, -1.660, -6.091;
                66.836834, 39.800, 5.982, 8.580, 0.0, -1.680, -6.393;
                67.369601, 18.560, 6.818, 8.110, 0.0, -1.956, -6.475;
                67.900868, 8.172, 7.708, 7.640, 0.0, -2.216, -6.545;
                68.431006, 3.397, 8.652, 7.170, 0.0, -2.492, -6.600;
                68.960312, 1.334, 9.650, 6.690, 0.0, -2.773, -6.650;
                118.750334, 940.300, 0.010, 16.640, 0.0, -0.439, 0.079;
                368.498246, 67.400, 0.048, 16.400, 0.0, 0.000, 0.000;
                424.763020, 637.700, 0.044, 16.400, 0.0, 0.000, 0.000;
                487.249273, 237.400, 0.049, 16.000, 0.0, 0.000, 0.000;
                715.392902, 98.100, 0.145, 16.000, 0.0, 0.000, 0.000;
                773.839490, 572.300, 0.141, 16.200, 0.0, 0.000, 0.000;
                834.145546, 183.100, 0.145, 14.700, 0.0, 0.000, 0.000];

            %% spectroscopic data for water-vapor %Table 2, modified in version P.676-11
            %            f0       b1    b2    b3   b4   b5   b6
            vapor = [22.235080 .1079 2.144 26.38 .76 5.087 1.00;
                67.803960 .0011 8.732 28.58 .69 4.930 .82;
                119.995940 .0007 8.353 29.48 .70 4.780 .79;
                183.310087 2.273 .668 29.06 .77 5.022 .85;
                321.225630 .0470 6.179 24.04 .67 4.398 .54;
                325.152888 1.514 1.541 28.23 .64 4.893 .74;
                336.227764 .0010 9.825 26.93 .69 4.740 .61;
                380.197353 11.67 1.048 28.11 .54 5.063 .89;
                390.134508 .0045 7.347 21.52 .63 4.810 .55;
                437.346667 .0632 5.048 18.45 .60 4.230 .48;
                439.150807 .9098 3.595 20.07 .63 4.483 .52;
                443.018343 .1920 5.048 15.55 .60 5.083 .50;
                448.001085 10.41 1.405 25.64 .66 5.028 .67;
                470.888999 .3254 3.597 21.34 .66 4.506 .65;
                474.689092 1.260 2.379 23.20 .65 4.804 .64;
                488.490108 .2529 2.852 25.86 .69 5.201 .72;
                503.568532 .0372 6.731 16.12 .61 3.980 .43;
                504.482692 .0124 6.731 16.12 .61 4.010 .45;
                547.676440 .9785 .158 26.00 .70 4.500 1.00;
                552.020960 .1840 .158 26.00 .70 4.500 1.00;
                556.935985 497.0 .159 30.86 .69 4.552 1.00;
                620.700807 5.015 2.391 24.38 .71 4.856 .68;
                645.766085 .0067 8.633 18.00 .60 4.000 .50;
                658.005280 .2732 7.816 32.10 .69 4.140 1.00;
                752.033113 243.4 .396 30.86 .68 4.352 .84;
                841.051732 .0134 8.177 15.90 .33 5.760 .45;
                859.965698 .1325 8.055 30.60 .68 4.090 .84;
                899.303175 .0547 7.914 29.85 .68 4.530 .90;
                902.611085 .0386 8.429 28.65 .70 5.100 .95;
                906.205957 .1836 5.110 24.08 .70 4.700 .53;
                916.171582 8.400 1.441 26.73 .70 5.150 .78;
                923.112692 .0079 10.293 29.00 .70 5.000 .80;
                970.315022 9.009 1.919 25.50 .64 4.940 .67;
                987.926764 134.6 .257 29.85 .68 4.550 .90;
                1780.000000 17506. .952 196.3 2.00 24.15 5.00];

            a1 = oxigen(:,2);
            a2 = oxigen(:,3);
            a3 = oxigen(:,4);
            a4 = oxigen(:,5);
            a5 = oxigen(:,6);
            a6 = oxigen(:,7);

            b1 = vapor(:,2);
            b2 = vapor(:,3);
            b3 = vapor(:,4);
            b4 = vapor(:,5);
            b5 = vapor(:,6);
            b6 = vapor(:,7);



            theta = 300.0/T;

            e = rho * T / 216.7;        % equation (4)

            %% Oxigen computation
            fi = oxigen(:,1);

            Si = a1 .* 1e-7 * p * theta.^3 .*exp(a2 * (1.0 - theta));       % equation (3)

            df = a3 .* 1e-4 .* (p .* theta .^ (0.8-a4) + 1.1 * e * theta);  % equation (6a)

            % Doppler broadening

            df = sqrt( df.*df + 2.25e-6);                                   % equation (6b)

            delta = (a5 + a6 * theta) * 1e-4 * (p + e) .* theta.^(0.8);     % equation (7)

            Fi = f ./ fi .* (  (df - delta .* (fi - f))./( (fi - f).^2 + df.^2  ) + ...
                (df - delta .* (fi + f))./( (fi + f).^2 + df.^2  ));        % equation (5)

            d = 5.6e-4 * (p + e) * theta.^(0.8);                            % equation (9)

            Ndf = f * p * theta.^2 *( 6.14e-5/(d * (1 + (f/d).^2) ) + ...
                1.4e-12 * p * theta.^(1.5)/(1 + 1.9e-5 * f.^(1.5)) );       % equation (8)

            % specific attenuation due to dry air (oxygen, pressure induced nitrogen
            % and non-resonant Debye attenuation), equations (1-2)

            g_0 = 0.182 * f * (sum(Si .* Fi) + Ndf);
            

            %% vapor computation

            fi = vapor(:,1);

            Si = b1 .* 1e-1 .* e .* theta.^3.5 .*exp(b2 .* (1.0 - theta));      % equation (3)

            df = b3 .* 1e-4 .* (p .* theta .^ (b4) + b5 .* e .* theta.^b6);     % equation (6a)

            % Doppler broadening

            df = 0.535 .* df + sqrt( 0.217* df.*df + 2.1316e-12 * fi.*fi/theta); % equation (6b)

            delta = 0;                                                           % equation (7)

            Fi = f ./ fi .* (  (df - delta .* (fi - f))./( (fi - f).^2 + df.^2  ) + ...
                (df - delta.* (fi + f))./( (fi + f).^2 + df.^2  ));              % equation (5)

            % specific attenuation due to water vapour, equations (1-2)

            g_w = 0.182 * f * (sum(Si .* Fi) );

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        Attachment A to Annex 1                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Dts, theta_0, psi] = straight_line_e2s(obj, Hs, Ht, phi_s, phi_t, delta)
            %%straight_line_e2s Computes geometry of straigth-line Earth-space path 
            %
            % Recommendation ITU-R P.619-5 Attachment A to Annex 1
            %
            
            % Inputs
            % Variable    Unit     Type     Description
            % Hs	      km	   float    Altitude of space station, asl
            % Ht          km       float    Altitude of earth-based station, asl
            % phi_s       deg      float    Latitude of sub-satellite point (zero for geostationary satellite)
            % phi_t       deg      float    Latitude of earth-based station
            % delta       deg      float    Difference in longitude between the sub-satellite point and the earth
            %                               based station, limited to less than half a circle, positive when the space
            %                               station is to the east of the earth-based station 
            % Outputs:
            %
            % Dts         km       float    Straight-line distance between earth and space stations
            % theta_0     deg      float    Elevation angle above horizontal of the straight-line
            % psi         deg      float    Azimuthal bearing of the straight line eastwards from true NorthD


            % Step 1 Calculate the distances of the space station and the earth-based station
            % from the center of the Earth

            Re = 6371;                                               %(18c)
            Rs = Re + Hs;                                            %(18a)
            Rt = Re + Ht;                                            %(18b)

            % Step2: Calculate Cartesian coordinates of the space station where the axes origin is
            % at the center of the Earth, the Z axis is directed northwards
            % (such that the north pole is on the positive Z axis), and the 
            % X axis is in the meridian of the earth-based station:

            X1 = Rs*cosd(phi_s)*cosd(delta);                         %(19a)
            Y1 = Rs*cosd(phi_s)*sind(delta);                         %(19b)
            Z1 = Rs*sind(phi_s);                                     %(19c)

            % Step 3: Rotate the Cartesian axes around the Y axis such that
            % the Z axis passess through the earth-based station, and then
            % move the origin, without rotation, such that the origin coincides
            % with the earth-based station

            X2 = X1*sind(phi_t) - Z1*cosd(phi_t);                    %(20a)
            Y2 = Y1;                                                 %(20b)
            Z2 = Z1*sind(phi_t) + X1*cosd(phi_t) - Rt;               %(20c)

            % Step 4: Calculate the straight-line distance between the
            % earth-based station and the space station

            Dts = sqrt(X2.^2 + Y2.^2 + Z2.^2);                        %(21)

            % Step 5: Calculate the length of the line represented by Dts
            % projected into the X,Y plane

            Gts = sqrt(X2.^2 + Y2.^2);                                %(22)

            % Step 6: Calculate the elevation angle of the straight line
            %from the earth-based station to the space station

            theta_0 = atan2d(Gts, Z2);                                %(23)

            % Step 7: Initially calculate the azimuthal bearing of the
            % straight line from the earth-based station to the space
            % station relative to true South

            psi = atand2(X2, Y2);                                     %(24)

            % Step 8: Reassign psi to be eastwards from true North by
            % subtracting it from a half-circle. Depending on the
            % implementation of the atan2 function, the bearing may need to
            % be processed into the range (0-360 degrees) The bearing is
            % indeterminate if the elevation angle represents a vertical
            % path

            psi = 180 - psi;              % TODO: Make sure this is correct

            % Equation (23) gives the elevation angle of the ray at the
            % earth-based station which would exist in the absence of
            % tropospheric refraction, sometime referred to as the
            % free-space elevation angle. The apparent elevation angle
            % theta can be estimated from theta_0 using equation (25) in
            % Attachment B

            % theta = fs2apparent(obj, theta_0, Ht);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        Attachment B to Annex 1                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function theta = fs2apparent(obj, theta_0, Ht)
            %%fs2apparent Computes apparent elevation from free space elevation 
            %
            % Recommendation ITU-R P.619-5 Attachment B to Annex 1
            %
            
            % Inputs
            % Variable    Unit     Type     Description
            % theta_0     deg      float    Free-space elevation angle
            % Ht          km       float    Altitude of the earth-based station, asl
            %
            % Outputs:
            % theta       deg      float    Apparent or actual elevation
            
            % For an earth-based station at altitude Ht less than 3 km and
            % for theta_0 between -1 and 10 degrees (including) the
            % following estimation holds

            if (Ht > 3)
                warning ("Earth-based station altitude is higher than 3 km, approximation may not hold");
            end

            Tfs1 = 1.728   + 0.5411*theta_0  + 0.03723*theta_0.^2;   %(26a)
            Tfs2 = 0.1815  + 0.06272*theta_0 + 0.01380*theta_0.^2;   %(26b)
            Tfs3 = 0.01727 + 0.008288*theta_0;                       %(26c)

            tau_fs = 1.0/(Tfs1 + Tfs2 * Ht + Tfs3 * Ht.^2);           %(26)

            theta = theta_0 + tau_fs;                                 %(25)
        end

        function theta_0 = apparent2fs(obj, theta, Ht)
            %%apparent2fs Computes free space elevation from the apparent elevation 
            %
            % Recommendation ITU-R P.619-5 Attachment B to Annex 1
            %
            
            % Inputs
            % Variable    Unit     Type     Description
            % theta       deg      float    Apparent or actual elevation
            % Ht          km       float    Altitude of the earth-based station, asl
            %
            % Outputs:
            % theta_0     deg      float    Free-space elevation angle
            
            T1 = 1.314   + 0.6437*theta  + 0.02869*theta.^2;         %(28a)
            T2 = 0.2305  + 0.09428*theta + 0.01096*theta.^2;         %(28b)
            T3 = 0.008583;                                           %(28c)

            tau = 1.0/(T1 + T2 * Ht + T3 * Ht.^2);                    %(28)

            theta_0 = theta - tau;                                    %(27)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        Attachment C to Annex 1                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Ag = atm_attenuation_s2E(obj, f, He, Hs, phi_e, phi_s, Dphi_e, h, rho, T, P, n, std_atm)
            %%atm_attenuation_s2E Computes attenuation along space-Earth path (Descending ray)
            %
            % Recommendation ITU-R P.619-5 Attachment C.4 to Annex 1
            %

            % Inputs
            % Variable    Unit     Type     Description
            % f           GHz      float    Frequency
            % He	      km	   float    Height of Earth station (above mean sea level)
            % Hs          km       float    Height of space-based station (above mean sea level)
            % phi_e       deg      float    Elevation angle of the main beam of the earth based station antenna
            % phi_s       deg      float    Elevation angle of the main beam of the space-based station antenna
            % Dphi_e      deg      float    Half-power beamwidth of the earth based station antenna
            % h           km       float    Vector of radial heights of atmospheric layers (length N+1)
            % rho         g/m^3    float    Vector of atmospheric water-vapour density height profile (length N)
            % T           K        float    Vector of atmospheric temperature height profile (length N)
            % P           hPa      float    Vector of atmospheric pressure height profile (length N)
            % n           n.u.     float    Vector of atmospheric refractive index height profile (length N)
            % std_atm     n.u.     bool     true if standard atmosphere is used, false otherwise
            %
            % Outputs:
            %
            % Ag          dB       float    Attenuation along space-Earth path

            Ag = 1e20;
            Reff = 6371;
            re = He + Reff;
            rs = Hs + Reff;
            r = h + Reff;
            rn1 = r(2:end);
            rn  = r(1:end-1);

            ns = n(end);
            ne = n(1);

            N = length(P);

            % atmospheric water vapor partial pressure height profile

            e = rho .* T / 216.7;

            % dry atmospheric pressure profile

            Pd = P - e;

            % Step by step procedure

            % Only if a nadir path phi_s ~= -90 compute Steps 1-3

            if (phi_s ~= -90)

                % Step 1: Use Snell's law in polar coordinates to calculate the
                % incident elevation angle phi_ce, in degrees, at the earth
                % station antenna as given in equation (35)

                phi_ce = acosd((rs*ns)/(re*ne)*cosd(phi_s));

                % Step 2: Determine if the calculated elevation angle phi_ce is
                % within the beamwidth of the earth station antenna. If yes,
                % proceed to Step 3, otherwise stop

                if abs(phi_ce - phi_e) > Dphi_e/2.0

                    Ag = 1e20;
                   
                    return

                end

                % Step 3: Determine if the line-of-sight between the two
                % antennas is free from ducting. If a standard atmosphere
                % is being used, ducting does not occur

                if (~std_atm)

                    eq37 = ((rs * ns) ./ (r .* n) ) * cosd(phi_s);

                    kk = find( eq37 >= 1, 1);

                    if (~isempty( kk ))

                        Ag = 1e20;

                        return
                    end

                end

            else % skip Steps 1 to 3 and go to Step 4

                % Step 4: Calculate the length of the slant path lns within
                % each layer from equation (34)

                d2 = ( ns ./ n * rs * cosd(phi_s) ).^2;
                lns = sqrt( rn1.^2 - d2 ) - sqrt( rn.^2  - d2);       %(34)

                Ag = 0;

                for i = 1:N

                    % Step 5: Calculate the atmospheric specific attenuation
                    % gamma_n within each layer in terms of the atmospheric
                    % parameters within the layer from equation (1) of Annex 1
                    % of Recommendation ITU-R P.676

                    [gamma_0, gamma_w] = p676d11_ga(obj, f, Pd(i), rho(i), T(i));

                    % Step 6: Calculate the total gaseous attenuation along
                    % the space-Earth slant path from equation (33)

                    Ag = Ag + (gamma_0 + gamma_w) * lns(i);

                end

            end

            return

        end

        function Ag = atm_attenuation_E2s(obj, f, He, Hs, phi_e, phi_s, Dphi_s, h, rho, T, P, n, std_atm, atm_type)
            %%atm_attenuation_E2s Computes attenuation along Earth-space paths (Ascending ray)
            %
            % Recommendation ITU-R P.619-5 Attachment C.5 to Annex 1
            %

            % Inputs
            % Variable    Unit     Type     Description
            % f           GHz      float    Frequency
            % He	      km	   float    Height of Earth station (above mean sea level)
            % Hs          km       float    Height of space-based station (above mean sea level)
            % phi_e       deg      float    Elevation angle of the main beam of the earth based station antenna
            % phi_s       deg      float    Elevation angle of the main beam of the space-based station antenna
            % Dphi_s      deg      float    Half-power beamwidth of the space based station antenna
            % h           km       float    Vector of radial heights of atmospheric layers (length N+1)
            % rho         g/m^3    float    Vector of atmospheric water-vapour density height profile (length N)
            % T           K        float    Vector of atmospheric temperature height profile (length N)
            % P           hPa      float    Vector of atmospheric pressure height profile (length N)
            % n           n.u.     float    Vector of atmospheric refractive index height profile (length N)
            % std_atm     n.u.     bool     true if standard atmosphere is used, false otherwise
            % atm_type    n.u.     int      1  - ITU-R P.835 mean annual global reference atmosphere
            %                               2  - ITU-R P.835 low-latitude reference atmosphere
            %                               31 - ITU-R P.835 summer mid-latitude reference atmosphere
            %                               32 - ITU-R P.835 winter mid-latitude reference atmosphere
            %                               41 - ITU-R P.835 summer high-latitude reference atmosphere
            %                               42 - ITU-R P.835 winter high-latitude reference atmosphere
            %
            %
            % Outputs:
            %
            % Ag          dB       float    Attenuation along Earth-space path

            Ag = 1e20;
            
            Reff = 6371;
            re = He + Reff;
            rs = Hs + Reff;
            r = h + Reff;
            rn1 = r(2:end);
            rn  = r(1:end-1);

            ns = n(end);
            ne = n(1);

            N = length(P);

            % atmospheric water vapor partial pressure height profile

            e = rho .* T / 216.7;

            % dry atmospheric pressure profile

            Pd = P - e;

            if (phi_e >= 0)

                % Case 1

                % Step by step procedure

                % Only if a nadir path phi_s ~= 90 compute Steps 1-3

                if (phi_s ~= 90)

                    % Step 1: Use Snell's law in polar coordinates to calculate the
                    % incident elevation angle phi_cs, in degrees, at the soace
                    % station antenna as given in equation (41)

                    phi_cs = acosd((re*ne)/(rs*ns)*cosd(phi_e));

                    % Step 2: Determine if the calculated elevation angle phi_cs is
                    % within the beamwidth of the space station antenna. If yes,
                    % proceed to Step 3, otherwise stop

                    if abs(phi_cs - phi_s) > Dphi_s/2.0

                        Ag = 1e20;

                        return

                    end

                    % Step 3: Determine if the line-of-sight between the two
                    % antennas is free from ducting. If a standard atmosphere
                    % is being used, ducting does not occur

                    if (~std_atm)

                        eq43 = ((rs * ns) ./ (r .* n) ) * cosd(phi_s);

                        kk = find( eq43 >= 1, 1);

                        if (~isempty( kk ))

                            Ag = 1e20;

                            return
                        end

                    end

                else % skip Steps 1 to 3 and go to Step 4

                    % Step 4: Calculate the length of the slant path lns within
                    % each layer from equation (34)

                    d2 = ( ne ./ n * re * cosd(phi_e) ).^2;
                    lne = sqrt( rn1.^2 - d2 ) - sqrt( rn.^2  - d2);       %(40)

                    Ag = 0;

                    for i = 1:N

                        % Step 5: Calculate the atmospheric specific attenuation
                        % gamma_n within each layer in terms of the atmospheric
                        % parameters within the layer from equation (1) of Annex 1
                        % of Recommendation ITU-R P.676

                        [gamma_0, gamma_w] = p676d11_ga(obj, f, Pd(i), rho(i), T(i));

                        % Step 6: Calculate the total gaseous attenuation along
                        % the space-Earth slant path from equation (39)

                        Ag = Ag + (gamma_0 + gamma_w) * lne(i);

                    end

                    return

                end

                
            else % Case 2: phi_e < 0

                % In this case, the attenuatin in equation (29) can be
                % calculated as the sum of two paths, one from the height of
                % the eart station to a virtual termial a t the minimum
                % altitude height Hmin and one from the virtual terminal to
                % the space station

                % This subrutine assumes that the first order derivative of
                % the atmospheric refractive index w.r.t. altitude is not
                % given. For that reason equation (45) is solved using
                % iterative method as given in Recommendation ITU-R P.676.
                % Only standard atmospheres are used

                if (~std_atm)

                    error('Case 2 is implemented only for standard atmospheres');

                end

                % Step 1: Solve transcendental equation (45) using iterative method
                % following Recommendation ITU-R P.676

                Hmin = solve45(obj, He, phi_e, atm_type);

                % Step 2: Upon determining the value of Hmin, a virtual
                % terminal with zero elevation angle can be treated as if it
                % were at this altitude:
                % one propagation path extending from the Earth-based
                % transmitter to the virtual terminal to account for the
                % first term of equation (44)
                % one propagation path extending from the virtual terminal
                % to the space-based receiver to account for the second term
                % of equation (44)

                % TODO: Resolve the following issues:
                % 1 - the beamwidth of the virtual terminal is not defined in the Recommendation
                %     It is assumed to be 180 degrees to make sure the calculated elevation angle is always
                %     within the half-power beamwidth of the virtual
                %     terminal -> not needed when using P.676 defined paths 
                % 2 - Profiles rho, T, P, n need to be recalculated
                %     for h = Hmin, He and then Hmin, Hs?
                % 3 - P.676 defines the two paths as follows: 
                %     Path1 is the gaseous attenuation betwee a virtual earth station at height hG and 
                %     the actual earth station at a height of h1 at an apparent elevation antle of 0 degs
                %     Path2 is the gaseous attenuation between a virtual earth station at height hG and the maximum
                %     atmospheric height (typically 100 km) at an apparent elevation angle of 0 deg
                %     P.619 defines the two paths as follows (slightly different to P.676)
                %     Path1 is extending from the earth-based transmitter to the virtual terminal to account for the first
                %     term of equation (44). If Hmin = hG < h1 = He, we cannot use the method in Case 1 as the path is descending
                %     Path2 is extending from the virtual terminal to the
                %     space-based receiver to account for the second term
                %     of equation (44)
                %     Here we assume Hmin < He < Hs and follow P.676

                % find lower and upper boundaries for each layer in path 1

                h = p676_slant_path_geometry16(obj, Hmin, He);
                
                % compute the midpoint for each layer
                
                hmid = 0.5*( h(1:end-1) + h(2:end));

                % compute T, P, rho, n profiles at midpoint of each layer

                [T, P, rho, n] = std_atm_profiles(obj, hmid, atm_type);

                % compute the gaseous attenuation for path 1

                Ag1 = atm_attenuation_E2s(obj, f, Hmin, He, 0, phi_e, Dphi_e, h, rho, T, P, n, std_atm, atm_type);

                % find lower and upper boundaries for each layer in path 1

                h = p676_slant_path_geometry16(obj, Hmin, Hs);

                % compute the midpoint for each layer

                hmid = 0.5*( h(1:end-1) + h(2:end));

                % compute T, P, rho, n profiles at midpoint of each layer
                
                [T, P, rho, n] = std_atm_profiles(obj, h, atm_type);

                % compute the gaseous attenuation for path 2

                Ag2 = atm_attenuation_E2s(obj, f, Hmin, Hs, 0, phi_s, Dphi_s, h, rho, T, P, n, std_atm, atm_type);

                % compute the total gaseous attenuation
                
                Ag = Ag1 + Ag2;

                return

            end

            return

        end

        function [Hmin, err] = solve45(obj, He, phi_e, atm_type, nit)
            %%solve45 Solves iteratively equation (45) from P.619-5
            %
            % Recommendation ITU-R P.619-5, Equation (45)
            %
            % Inputs
            % Variable    Unit     Type     Description
            % He          km	   float    Geometric height of the earth station 
            % phi_e       deg      float    Elevation angle of the main beam of the earth based station antenna 
            % atm_type    n.u.     int      1  - ITU-R P.835 mean annual global reference atmosphere
            %                               2  - ITU-R P.835 low-latitude reference atmosphere
            %                               31 - ITU-R P.835 summer mid-latitude reference atmosphere
            %                               32 - ITU-R P.835 winter mid-latitude reference atmosphere
            %                               41 - ITU-R P.835 summer high-latitude reference atmosphere
            %                               42 - ITU-R P.835 winter high-latitude reference atmosphere
            % nit         n.u.     int      number of iterations
            % Outputs:
            % Hmin        km       float    Altitude of the virtual terminal where the radio beam is parallel to the Earth's surface
            
            Reff = 6371;

            [~, ~, ~, ne] = p835_std_atm_profiles(obj, He, atm_type);
           
            C = (Reff + He) * ne * cosd(phi_e);

            Hmin_old = He;
            Hmin_new = He;
            
            for i = 1:nit

                Hmin_old = Hmin_new;

                [~, ~, ~, nmin] = p835_std_atm_profiles(obj, Hmin_old, atm_type);

                Hmin_new = C / nmin - Reff;

                        
            end

            Hmin = Hmin_new;
            err = Hmin_new - Hmin_old;

            return

        end

            

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Recommendation ITU-R P.835-6: Reference standard atmospheres           %
        %  Annex 1                                                                %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [T, P, rho] = p835_s1_old(obj, h)
            %% p835_s1_old Computes T, P, rho height profiles for mean annual global reference atmosphere
            %
            % Recommendation ITU-R P.835-6 §1
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Geometric height 
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile


            if (0 <= h && h <= 86)
                % first height regime

                % Convert geometic height to geopotential height (1a)
                
                hp = 6356.766 * h / (6356.766 + h);

                if (0 <= hp && hp <= 11 )
                    
                    T = 288.15 - 6.5 * hp;                  %(2a)
                    P = 1013.25*(288.15/(288.15-6.5*hp)).^(-34.1632/6.5);   
                                                            %(3a)

                elseif ( hp <= 20 )

                    T = 216.65;                             %(2b)
                    P = 226.3226*exp(-34.1632*(hp-11)/216.65);
                                                            %(3b)

                elseif ( hp <= 32)

                    T = 216.65 + (hp - 20);                 %(2c)
                    P = 54.74980*(216.65/(216.65+hp-20)).^(34.1632);
                                                            %(3c)


                elseif (hp <= 47)

                    T = 228.65 + 2.8 * (hp - 32);           %(2d)
                    P = 8.680422*(228.65/(228.65+2.8*(hp-32))).^(34.1632/2.8);
                                                            %(3d)

                elseif (hp <= 51)

                    T = 270.65;                             %(2e)
                    P = 1.109106*exp(-34.1632*(hp-47)/270.65);  %(3e)

                elseif (hp <= 71)

                    T = 270.65 - 2.8 * (hp - 51);           %(2f)
                    P = 0.6694167*(270.65/(270.65-2.8*(hp-51))).^(-34.1632/2.8);
                                                            %(3f)

                else

                    T = 214.65 - 2.0 * (hp - 71);           %(2g)
                    P = 0.03956649*(214.65/(214.65-2.0*(hp-71))).^(-34.1632/2.0);
                                                            %(3g)

                end

            elseif (h <= 100)

                % Second height regime

                if (h <= 91)

                    T = 186.8673;                           %(4a)

                else

                    T = 263.1905 - 76.3232 * sqrt(1-((h-91)/19.9429).^2); 
                                                            %(4b)

                end

                a0 = 95.571899;
                a1 = -4.011801;
                a2 = 6.424731e-2;
                a3 = -4.789660e-4;
                a4 = 1.340543e-6;

                P = exp(a0 + a1*h + a2*h^2 + a3*h^3 + a4*h^4);   %(5)

            else
                error('Geometric height needs to be in the range 0, 100 km.')

            end

            % Section 1.2: Water-vapour pressure

            rho0 = 7.5;                   %(7)
            h0 = 2;              

            rho = rho0*exp(-h/h0);        %(6)


        end


        function [T, P, rho] = p835_s1(obj, h)
            %% p835_s1 Computes T, P, rho height profiles for mean annual global reference atmosphere
            %
            % Recommendation ITU-R P.835-6 §1
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Vector of geometric heights
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile


            if( min(h) < 0 || max(h) > 100 )
               
                error('Geometric heights must be within the range 0, 100 km.')

            end

            kk = (0 <= h & h <= 86);

            % first height regime

            % Convert geometic height to geopotential height (1a)

            hp = 6356.766 .* h(kk) ./ (6356.766 + h(kk));

            %if (0 <= hp && hp <= 11 )
            ll = (0 <= hp & hp <= 11);

            T(ll)  = 288.15 - 6.5 * hp(ll);                  %(2a)
            P(ll) = 1013.25.*(288.15./(288.15-6.5*hp(ll))).^(-34.1632/6.5);
            %(3a)

            ll = ( hp > 11 & hp <= 20 );

            T(ll) = 216.65;                             %(2b)
            P(ll) = 226.3226.*exp(-34.1632.*(hp(ll)-11)./216.65);
            %(3b)

            ll = (hp> 20 & hp <= 32);

            T(ll) = 216.65 + (hp(ll) - 20);                 %(2c)
            P(ll) = 54.74980.*(216.65./(216.65+hp(ll)-20)).^(34.1632);
            %(3c)


            ll = (hp > 32 & hp <= 47);

            T(ll) = 228.65 + 2.8 * (hp(ll) - 32);           %(2d)
            P(ll) = 8.680422.*(228.65./(228.65+2.8.*(hp(ll)-32))).^(34.1632/2.8);
            %(3d)

            ll = (hp > 47 & hp <= 51);

            T(ll) = 270.65;                             %(2e)
            P(ll) = 1.109106.*exp(-34.1632.*(hp(ll)-47)./270.65);  %(3e)

            ll = (hp > 51 & hp <= 71);

            T(ll) = 270.65 - 2.8 * (hp(ll) - 51);           %(2f)
            P(ll) = 0.6694167.*(270.65./(270.65-2.8.*(hp(ll)-51))).^(-34.1632/2.8);
            %(3f)

            ll = (hp > 71);

            T(ll) = 214.65 - 2.0 * (hp(ll) - 71);           %(2g)
            P(ll) = 0.03956649.*(214.65./(214.65-2.0.*(hp(ll)-71))).^(-34.1632/2.0);
            %(3g)


            % Second height regime

            ll = (h > 86 & h <= 91); 
                
            T(ll)  = 186.8673;                           %(4a)

            ll = (h > 91 & h <= 100); 
            
            T(ll) = 263.1905 - 76.3232 * sqrt(1-((h(ll)-91)/19.9429).^2);
            %(4b)

            a0 = 95.571899;
            a1 = -4.011801;
            a2 = 6.424731e-2;
            a3 = -4.789660e-4;
            a4 = 1.340543e-6;

            ll = (h > 86 & h <= 100);

            P(ll) = exp(a0 + a1.*h(ll) + a2.*h(ll).^2 + a3.*h(ll).^3 + a4.*h(ll).^4);   %(5)



            % Section 1.2: Water-vapour pressure

            rho0 = 7.5;                   %(7)
            h0 = 2;

            rho = rho0 * exp(-h./h0);        %(6)


            %Column vectors
            T = T.';
            P = P.';
            


        end


        function [T, P, rho] = p835_s2(obj, h)
            %% p835_s2 Computes T, P, rho height profiles for low-latitude reference atmosphere
            % For low latitudes (smaller than 22 deg) the seasonal
            % variations are not very important and a single annual profile
            % can be used
            %
            % Recommendation ITU-R P.835-6 §2
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Vector of geometric heights 
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile

            if (min(h) < 0 || max(h) >100)
                error ('Geometric heights must be within the range 0-100 km');
            end

            kk = (0 <= h & h < 17);

            T(kk) = 300.4222 - 6.3533*h(kk)  + 0.005886 * h(kk).^2;

            kk = (h>=17 & h < 47);

            T(kk) = 194 + (h(kk)-17)*2.533;

            kk = (h>= 47 & h < 52);

            T(kk) = 270;

            kk = (h >= 52 & h < 80);

            T(kk) = 270 - (h(kk)-52)*3.0714;

            kk = (h >= 80 & h <= 100);

            T(kk) = 184;



            kk = (0 <= h & h <= 10);

            P(kk) = 1012.0306 - 109.0338 * h(kk) + 3.6316 * h(kk).^2;

            kk = ( h > 10 & h <= 72);

            P10 = 1012.0306 - 109.0338 * 10 + 3.6316 * 100;
            P(kk) = P10.*exp(-0.147*(h(kk)-10));

            kk = (h > 72 & h <= 100);

            P10 = 1012.0306 - 109.0338 * 10 + 3.6316 * 100;
            P72 = P10*exp(-0.147*(72-10));
            P(kk) = P72.*exp(-0.165*(h(kk)-72));


            kk = ( 0  <= h & h <= 15);

            rho(kk) = 19.6542.*exp(-0.2313.*h(kk) - 0.1122.*h(kk).^2 + 0.01351.*h(kk).^3 - 0.0005923.*h(kk).^4);

            kk = ( h > 15 & h <= 100);
            rho(kk) = 0;


            %Column vectors
            T = T.';
            P = P.';
            rho = rho.';



        end

        function [T, P, rho] = p835_s2_old(obj, h)
            %% p835_s2 Computes T, P, rho height profiles for low-latitude reference atmosphere
            % For low latitudes (smaller than 22 deg) the seasonal
            % variations are not very important and a single annual profile
            % can be used
            %
            % Recommendation ITU-R P.835-6 §2
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Geometric height 
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile



           if (0 <= h && h < 17)

               T = 300.4222 - 6.3533*h  + 0.005886 * h.^2;

           elseif (h < 47)

               T = 194 + (h-17)*2.533;

           elseif (h < 52)

               T = 270;

           elseif (h < 80)

               T = 270 - (h-52)*3.0714;

           elseif (h <= 100)
               
               T = 184;

           else
                error ('Geometric height must be within the range 0-100 km');
           end 


           if (0 <= h && h <= 10)
               
               P = 1012.0306 - 109.0338 * h + 3.6316 * h.^2;

           elseif (h <= 72)

               P10 = 1012.0306 - 109.0338 * 10 + 3.6316 * 100;
               P = P10*exp(-0.147*(h-10));

           elseif (h <= 100)

               P10 = 1012.0306 - 109.0338 * 10 + 3.6316 * 100;
               P72 = P10*exp(-0.147*(72-10));
               P = P72*exp(-0.165*(h-72));

           else
                error ('Geometric height must be within the range 0-100 km');
           end
              
           if ( 0  <= h && h <= 15)

               rho = 19.6542*exp(-0.2313*h - 0.1122*h.^2 + 0.01351*h.^3 - 0.0005923*h.^4);

           else
               rho = 0;
           end


        end


        function [T, P, rho] = p835_s31(obj, h)
            %% p835_s31 Computes T, P, rho height profiles for summer mid-latitude reference atmosphere
            % Mid-latitudes (between 22 deg and 45 deg)
            %
            % Recommendation ITU-R P.835-6 §3.1
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Vector of geometric heights
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile

            if (min(h) < 0 || max(h) >100)
                error ('Geometric heights must be within the range 0-100 km');
            end

            kk = (0 <= h & h < 13);

            T(kk) = 294.9838 - 5.2159*h(kk)  - 0.07109 * h(kk).^2;

            kk = (h >=13 & h < 17);

            T(kk) = 215.15;

            kk = (h >=17 & h < 47);

            T(kk) = 215.15.*exp((h(kk)-17).*0.008128);

            kk = (h >= 47 & h < 53);

            T(kk) = 275;

            kk = (h >= 53 & h < 80);

            T(kk) = 275 + (1-exp((h(kk)-53)*0.06))*20;
            % NOTE: this function is discontinuous at h = 80 as
            % T(80) = 193.9382 \neq 175
            % Replacing 0.06 with 0.06636145 gives T(80) = 175.0000

            kk = (h >= 80 & h <= 100);

            T(kk) = 175;

            kk = (0 <= h & h <= 10);

            P(kk) = 1012.8186 - 111.5569 * h(kk) + 3.8646 * h(kk).^2;

            kk = (h > 10 & h <= 72);

            P10 = 1012.8186 - 111.5569 * 10 + 3.8646 * 100;
            P(kk) = P10*exp(-0.147*(h(kk)-10));

            kk = (h > 72 & h <= 100);

            P10 = 1012.8186 - 111.5569 * 10 + 3.8646 * 100;
            P72 = P10*exp(-0.147*(72-10));
            P(kk) = P72*exp(-0.165*(h(kk)-72));


            kk = ( 0  <= h & h <= 15);

            rho(kk) = 14.3542*exp(-0.4174*h(kk) - 0.02290*h(kk).^2 + 0.001007*h(kk).^3);

            kk = (h > 15 & h <=100);

            rho(kk) = 0;


            %Column vectors
            T = T.';
            P = P.';
            rho = rho.';



        end

        function [T, P, rho] = p835_s31_old(obj, h)
            %% p835_s31 Computes T, P, rho height profiles for summer mid-latitude reference atmosphere
            % Mid-latitudes (between 22 deg and 45 deg) 
            %
            % Recommendation ITU-R P.835-6 §3.1
            %
            
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Geometric height 
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile



           if (0 <= h && h < 13)

               T = 294.9838 - 5.2159*h  - 0.07109 * h.^2;

           elseif (h < 17)

               T = 215.15;

           elseif (h < 47)

               T = 215.15*exp((h-17)*0.008128);

           elseif (h < 53)

               T = 275;

           elseif (h < 80)
               
               T = 275 + (1-exp((h-53)*0.06))*20;
               % NOTE: this function is discontinuous at h = 80 as 
               % T(80) = 193.9382 \neq 175
               % Replacing 0.06 with 0.06636145 gives T(80) = 175.0000

           elseif (h <= 100)

               T = 175;

           else
                error ('Geometric height must be within the range 0-100 km');
           end 


           if (0 <= h && h <= 10)
               
               P = 1012.8186 - 111.5569 * h + 3.8646 * h.^2;

           elseif (h <= 72)

               P10 = 1012.8186 - 111.5569 * 10 + 3.8646 * 100;
               P = P10*exp(-0.147*(h-10));

           elseif (h <= 100)

               P10 = 1012.8186 - 111.5569 * 10 + 3.8646 * 100;
               P72 = P10*exp(-0.147*(72-10));
               P = P72*exp(-0.165*(h-72));

           else
                error ('Geometric height must be within the range 0-100 km');
           end
              
           if ( 0  <= h && h <= 15)

               rho = 14.3542*exp(-0.4174*h - 0.02290*h.^2 + 0.001007*h.^3);

           else
               rho = 0;
           end


        end


        function [T, P, rho] = p835_s32(obj, h)
            %% p835_s32 Computes T, P, rho height profiles for winter mid-latitude reference atmosphere
            % Mid-latitudes (between 22 deg and 45 deg)
            %
            % Recommendation ITU-R P.835-6 §3.2
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Vector of geometric heights
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile

            if (min(h) < 0 || max(h) >100)
                error ('Geometric heights must be within the range 0-100 km');
            end

            kk = (0 <= h & h < 10);

            T(kk) = 272.7241 - 3.6217*h(kk)  - 0.1759 * h(kk).^2;

            kk = (h >=10 & h < 33);

            T(kk) = 218;

            kk = ( h >= 33 &h < 47);

            T(kk) = 218 + (h(kk)-33)*3.3571;

            kk = (h >= 47 & h < 53);

            T(kk) = 265;

            kk =  (h >= 53 & h < 80);

            T(kk) = 265 - (h(kk)-53)*2.0370;

            kk =  (h >=80 & h <= 100);

            T(kk) = 210;



            kk = (0 <= h & h <= 10);

            P(kk) = 1018.8627 - 124.2954 * h(kk) + 4.8307 * h(kk).^2;

            kk = (h > 10 & h <= 72);

            P10 = 1018.8627 - 124.2954 * 10 + 4.8307 * 100;
            P(kk) = P10*exp(-0.147*(h(kk)-10));

            kk = (h > 72 & h <= 100);

            P10 = 1018.8627 - 124.2954 * 10 + 4.8307 * 100;
            P72 = P10*exp(-0.147*(72-10));
            P(kk) = P72*exp(-0.165*(h(kk)-72));


            kk = ( 0  <= h & h <= 15);

            rho(kk) = 3.4742*exp(-0.2697*h(kk) - 0.03604*h(kk).^2 + 0.0004489*h(kk).^3);

            kk = (h > 15 & h<=100);
            rho(kk) = 0;


            %Column vectors
            T = T.';
            P = P.';
            rho = rho.';

        end

        function [T, P, rho] = p835_s32_old(obj, h)
            %% p835_s32 Computes T, P, rho height profiles for winter mid-latitude reference atmosphere
            % Mid-latitudes (between 22 deg and 45 deg) 
            %
            % Recommendation ITU-R P.835-6 §3.2
            %
            
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Geometric height 
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile



           if (0 <= h && h < 10)

               T = 272.7241 - 3.6217*h  - 0.1759 * h.^2;

           elseif (h < 33)

               T = 218;

           elseif (h < 47)

               T = 218 + (h-33)*3.3571;

           elseif (h < 53)

               T = 265;

           elseif (h < 80)
               
               T = 265 - (h-53)*2.0370;

           elseif (h <= 100)

               T = 210;

           else
                error ('Geometric height must be within the range 0-100 km');
           end 


           if (0 <= h && h <= 10)
               
               P = 1018.8627 - 124.2954 * h + 4.8307 * h.^2;

           elseif (h <= 72)

               P10 = 1018.8627 - 124.2954 * 10 + 4.8307 * 100;
               P = P10*exp(-0.147*(h-10));

           elseif (h <= 100)

               P10 = 1018.8627 - 124.2954 * 10 + 4.8307 * 100;
               P72 = P10*exp(-0.147*(72-10));
               P = P72*exp(-0.165*(h-72));

           else
                error ('Geometric height must be within the range 0-100 km');
           end
              
           if ( 0  <= h && h <= 15)

               rho = 3.4742*exp(-0.2697*h - 0.03604*h.^2 + 0.0004489*h.^3);

           else
               rho = 0;
           end


        end

        function [T, P, rho] = p835_s41(obj, h)
            %% p835_s41 Computes T, P, rho height profiles for summer high-latitude reference atmosphere
            % High-latitudes (higher than 45 deg)
            %
            % Recommendation ITU-R P.835-6 §4.1
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Vector of geometric heights
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile

            if (min(h) < 0 || max(h) >100)
                error ('Geometric heights must be within the range 0-100 km');
            end

            kk = (0 <= h & h < 10);

            T(kk) = 286.8374 - 4.7805*h(kk)  - 0.1402 * h(kk).^2;

            kk = (h >=10 & h < 23);

            T(kk) = 225;

            kk = (h >= 23 & h < 48);

            T(kk) = 225*exp( (h(kk)-23)*0.008317 );

            kk = (h >=48 & h < 53);

            T(kk) = 277;

            kk = (h >= 53 & h < 79);

            T(kk) = 277 - (h(kk)-53)*4.0769;

            kk = (h >=79 & h <= 100);

            T(kk) = 171;



            kk = (0 <= h & h <= 10);

            P(kk) = 1018.0278 - 113.2494 * h(kk) + 3.9408 * h(kk).^2;

            kk = (h > 10 & h <= 72);

            P10 = 1018.0278 - 113.2494 * 10 + 3.9408 * 100;
            P(kk) = P10*exp(-0.140*(h(kk)-10));

            kk = (h > 72 & h <= 100);

            P10 = 1018.0278 - 113.2494 * 10 + 3.9408 * 100;
            P72 = P10*exp(-0.140*(72-10));
            P(kk) = P72*exp(-0.165*(h(kk)-72));


            kk = ( 0  <= h & h <= 15);

            rho(kk) = 8.988*exp(-0.3614*h(kk) - 0.005402*h(kk).^2 - 0.001955*h(kk).^3);

            kk = (h > 15 & h <= 100);

            rho(kk) = 0;

            %Column vectors
            T = T.';
            P = P.';
            rho = rho.';


        end

        function [T, P, rho] = p835_s41_old(obj, h)
            %% p835_s41 Computes T, P, rho height profiles for summer high-latitude reference atmosphere
            % High-latitudes (higher than 45 deg) 
            %
            % Recommendation ITU-R P.835-6 §4.1
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Geometric height 
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile



           if (0 <= h && h < 10)

               T = 286.8374 - 4.7805*h  - 0.1402 * h.^2;

           elseif (h < 23)

               T = 225;

           elseif (h < 48)

               T = 225*exp( (h-23)*0.008317 );

           elseif (h < 53)

               T = 277;

           elseif (h < 79)
               
               T = 277 - (h-53)*4.0769;

           elseif (h <= 100)

               T = 171;

           else
                error ('Geometric height must be within the range 0-100 km');
           end 


           if (0 <= h && h <= 10)
               
               P = 1018.0278 - 113.2494 * h + 3.9408 * h.^2;

           elseif (h <= 72)

               P10 = 1018.0278 - 113.2494 * 10 + 3.9408 * 100;
               P = P10*exp(-0.140*(h-10));

           elseif (h <= 100)

               P10 = 1018.0278 - 113.2494 * 10 + 3.9408 * 100;
               P72 = P10*exp(-0.140*(72-10));
               P = P72*exp(-0.165*(h-72));

           else
                error ('Geometric height must be within the range 0-100 km');
           end
              
           if ( 0  <= h && h <= 15)

               rho = 8.988*exp(-0.3614*h - 0.005402*h.^2 - 0.001955*h.^3);

           else
               rho = 0;
           end

        end

        function [T, P, rho] = p835_s42(obj, h)
            %% p835_s42 Computes T, P, rho height profiles for winter high-latitude reference atmosphere
            % High-latitudes (higher than 45 deg)
            %
            % Recommendation ITU-R P.835-6 §4.2
            %

            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Geometric height
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile

            
            if (min(h) < 0 || max(h) >100)
                error ('Geometric heights must be within the range 0-100 km');
            end

            kk = (0 <= h & h < 8.5);

            T(kk) = 257.4345 + 2.3474*h(kk) - 1.5479 * h(kk).^2 + 0.08473 * h(kk).^3;

            kk = (h >= 8.5 & h < 30);

            T(kk) = 217.5;

            kk = (h >= 30 & h < 50);

            T(kk) = 217.5 + ( (h(kk)-30)*2.125 );

            kk = (h >= 50 & h < 54);

            T(kk) = 260;

            kk = (h >= 54 & h <= 100);

            T(kk) = 260 - (h(kk)-54)*1.667;

            kk = (0 <= h & h <= 10);

            P(kk) = 1010.8828 - 122.2411 * h(kk)+ 4.554 * h(kk).^2;

            kk = (h > 10 & h <= 72);

            P10 = 1010.8828 - 122.2411 * 10 + 4.554 * 100;
            P(kk) = P10*exp(-0.147*(h(kk)-10));

            kk = (h > 72 & h <= 100);

            P10 = 1010.8828 - 122.2411 * 10 + 4.554 * 100;
            P72 = P10*exp(-0.147*(72-10));
            P(kk) = P72*exp(-0.150*(h(kk)-72));


            kk = ( 0  <= h & h <= 15);

            rho(kk) = 1.2319*exp(0.07481*h(kk) - 0.0981*h(kk).^2 + 0.00281*h(kk).^3);

            kk = (h > 15 & h <= 100);
            rho(kk) = 0;


            %Column vectors
            T = T.';
            P = P.';
            rho = rho.';


        end

        function [T, P, rho] = p835_s42_old(obj, h)
            %% p835_s42 Computes T, P, rho height profiles for winter high-latitude reference atmosphere
            % High-latitudes (higher than 45 deg) 
            %
            % Recommendation ITU-R P.835-6 §4.2
            %
            
            % Inputs
            % Variable    Unit     Type     Description
            % h	          km	   float    Geometric height 
            %
            % Outputs:
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Water vapour profile



           if (0 <= h && h < 8.5)

               T = 257.4345 + 2.3474*h  - 1.5479 * h.^2 + 0.08473 * h.^3;

           elseif (h < 30)

               T = 217.5;

           elseif (h < 50)

               T = 217.5 + ( (h-30)*2.125 );

           elseif (h < 54)

               T = 260;

           elseif (h <= 100)
               
               T = 260 - (h-54)*1.667;

           else
                error ('Geometric height must be within the range 0-100 km');
           end 


           if (0 <= h && h <= 10)
               
               P = 1010.8828 - 122.2411 * h + 4.554 * h.^2;

           elseif (h <= 72)

               P10 = 1010.8828 - 122.2411 * 10 + 4.554 * 100;
               P = P10*exp(-0.147*(h-10));

           elseif (h <= 100)

               P10 = 1010.8828 - 122.2411 * 10 + 4.554 * 100;
               P72 = P10*exp(-0.147*(72-10));
               P = P72*exp(-0.150*(h-72));

           else
                error ('Geometric height must be within the range 0-100 km');
           end
              
           if ( 0  <= h && h <= 15)

               rho = 1.2319*exp(0.07481*h - 0.0981*h.^2 + 0.00281*h.^3);

           else
               rho = 0;
           end

        end

        function [T, P, rho, n] = p835_std_atm_profiles(obj, h, atm_type)
            %% std_atm_profiles Computes T, P, rho, n profiles for a vector of heights h for std atmospheres
            %
            % Recommendation ITu-R P.835-6, Annex 1
            % Recommendation ITU-R P.453-14
            % Inputs
            % Variable    Unit     Type     Description
            % h           km       float    Vector of geometric heights 
            % atm_type    n.u.     int      1  - ITU-R P.835 mean annual global reference atmosphere
            %                               2  - ITU-R P.835 low-latitude reference atmosphere
            %                               31 - ITU-R P.835 summer mid-latitude reference atmosphere
            %                               32 - ITU-R P.835 winter mid-latitude reference atmosphere
            %                               41 - ITU-R P.835 summer high-latitude reference atmosphere
            %                               42 - ITU-R P.835 winter high-latitude reference atmosphere

            % Outputs:
            % T           K        float    Temperature profile 
            % P           hPa      float    Pressure profile 
            % rho         g/m^3    float    water-vapour density profile 
            % n           -        float    Athmospheric radio refractivindex profile 

            if (atm_type == 1)

                [T, P, rho] = p835_s1(obj, h);

            elseif (atm_type == 2)

                [T, P, rho] = p835_s2(obj, h);

            elseif (atm_type == 31)

                [T, P, rho] = p835_s31(obj, h);

            elseif (atm_type == 32)

                [T, P, rho] = p835_s32(obj, h);

            elseif (atm_type == 41)

                [T, P, rho] = p835_s41(obj, h);

            elseif (atm_type == 42)

                [T, P, rho] = p835_s42(obj, h);

            else

                error('This function is only defined for standard atmospheres from Rec. ITU-R P.835-6');

            end


            n = p453_n(obj, T, P, rho);

            
            return
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Recommendation ITU-R P.453-14: Formula for refractive index            %
        %  Annex 1, equations (1) and (6)                                         %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function n = p453_n(obj, T, P, rho)
            %%p453_n Computes the atmospheric radio refractive index 
            %
            % Recommendation ITU-R P.453-14 §1
            % Recommendation ITU-R P.835-6  §1.2
            
            % Inputs
            % Variable    Unit     Type     Description
            % T           K        float    Temperature profile
            % P           hPa      float    Pressure profile
            % rho         g/m^3    float    Standard ground-level water-vapour density
            % h           km       float    Vector of geometric heights
            %
            % Outputs:
            % n           -        float    Athmospheric radio refractive index

            % P.835, Section 1.2: Water-vapour pressure

            % Vapour pressure

            e = rho .* T / 216.7;          %(8) P.835
            
            % Water-vapour density decreases exponantially with increasing
            % altitude, up to an altitude where the mixing ration e(h)/P(h)
            % = 2e-6. Above this altitude, the mixing ration is assumed to
            % be constant

            kk = find(e < 2e-6*P);

            e(kk) = 2e-6*P(kk);          % P.835
            %TODO: this condition is mentioned in P.835 but neither in
            %P.619 nor in P.676

            N = 77.6*P./T - 5.6*e./T + 3.75e5 *e./(T.^2);   %(6) P.453

            n = 1 + N*1e-6;                            %(1) P.453

        end           

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Recommendation ITU-R P.676-13: Slant path geometry                     %
        %  Annex 1, §2.2.1, equation 15                                           %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function h = p676_slant_path_geometry15(obj)
            %%p676_slant_path_geometry15 Computes a vector of exponentially increasing
            %%layer heights for gaseous attenuation calculation
            %
            % Recommendation ITU-R P.676-14 §2.2.1 (14-15)
            
            
            % Inputs
            % Variable    Unit     Type     Description
            %
            % Outputs:
            % h           km       float    Vector of slant path heights

            n = 923;
    
            i = (1:n).';

            h = 0.0001*(exp(0.01*(i-1))-1)./(exp(0.01)-1);
            h(end) = 100;

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Recommendation ITU-R P.676-13: Slant path geometry                     %
        %  Annex 1, §2.2.1, equation 16                                           %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function h = p676_slant_path_geometry16(obj, hL, hU)
            %%p676_slant_path_geometry Computes a vector of exponentially increasing
            %%layer heights for gaseous attenuation calculation
            %
            % Recommendation ITU-R P.676-14 §2.2.1 (16)
            
            
            % Inputs
            % Variable    Unit     Type     Description
            % hL          km       float    Lower height of the slant path
            % hU          km       float    Upper height of the slant path
            %
            % Outputs:
            % h           km       float    Vector of slant path heights
            

            iL = floor(100*log(1e4*hL*(exp(0.01)-1)+1)+1);           %(16a)
            iU =  ceil(100*log(1e4*hU*(exp(0.01)-1)+1)+1);           %(16b)

            n = iU - iL + 1;

            h = zeros(n,1);
    
            m = (exp(0.02)-exp(0.01)) / (exp(iU/100.0)-exp(iL/100.0)) * ( hU - hL );
                                                                     %(16c)

            for i = iL : iU

                h(i - iL + 1) = hL + m * (exp(0.01*(i-1)) - exp(0.01*(iL-1)))/(exp(0.01)-1);

            end

            h = h.'; %column vector

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Recommendation ITU-R P.676-13: Slant path geometry                     %
        %  Annex 1, §2.2.1, equation 17                                           %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function a = p676_slant_path_geometry17(obj, h, n, phi)
            %%p676_slant_path_geometry Computes a vector of path lengths  
            %
            % Recommendation ITU-R P.676-14 §2.2.1 (17)
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % h           km       float    Vector of layer heights. Layer #i bounds are (h(i), h(i+1)), 
            %                               so that the size of this vector is N+1, where N is the number of layers
            % n           km       float    Vector of refractive indices at the middle of each layer 
            %                               at the height 0.5*(h(i+1) + h(i)) of size N
            % phi         deg      float    Local apparent elevation angle at height h(1)    
            %
            % Outputs:
            % a           km       float    Vector of path lengths for each layer, of size N 
            
            % Local zenith angle at or near the surface of the Earth
            % (complement of the apparent elevation angle)

            beta1 = 90 - phi;

            rE = 6371;

            r = rE + h(1:end-1);
            ri1 = rE + h(2:end);
            delta = ri1-r;

            beta  = asind(n(1)*r(1) ./(n.*r)  *sind(beta1));
            % alpha = asind(n(1)*r(1) ./(n.*ri1)*sind(beta1));
            
            a = -r.*cosd(beta) + sqrt(r.^2 .* (cosd(beta)).^2 + 2*r.*delta + delta.^2);
           
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        Attachment D to Annex 1                          %
        %                       Tropospheric scintillation                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Ast = tropospheric_scintillation(obj, f, p, phi_e, phi_n, theta, Ga)
            %% tropospheric_scintillation Computes the tropospheric scintillation attenuation
            %
            % Recommendation ITU-R P.619-5 Attachment D
            % Recommendation ITU-R P.618-12, Section 2.4.1 for scintillation intensity
            % Recommendation ITU-R P.453-14 for wet term of the surface refractivity
            %
            %
            % Inputs
            % Variable    Unit     Type     Description
            % f           GHz      float    Frequency (4 GHz <=f <= 20 GHz)
            % p           %        float    Time percentage (within the range 0 - 100)
            % pni_e       deg      float    Longitude of the earth-station site
            % phi_n       deg      float    Latitude of the earth-station site
            % theta       deg      float    Free-space elevation angle (theta >= 5)
            % Ga          dBi      float    Earth-based antenna gain in the direction of the path
            %
            % Outputs:
            % Ast         dB       float    Tropospheric scintillation not exceeded for p percent time
            %                               Negative values for p < 50 indicate an enhancement in singal level


            if (f < 4 || f > 20)
                warning('This function is defined for frequencies within the range 4 GHz to 20 GHz.');
                % The function computation continues with this warning
            end

            if (theta < 5)
                warning('This function is defined for free-space elevation angles larger than or equal to 5 degs.');   
                Ast = 0;
                return;
            end

            % Use method given in § 2.4.1 of Recommendation ITU-R P.618-12
            % to calculate the scintillation intensity.

            % In this function it is assumed that the wet term of the surface
            % refractivity exceeded for the average year is obtained from
            % the digital maps in Recommendation ITU-R P.453 and therefore
            % it starts with Step 3

            Nwet = get_interp2_Nwet_Annual_time_location(obj, p, phi_e, phi_n);

            % Step 3: Calculate the standard deviation of the  reference signal amplitude

            sigma_ref = 3.6e-3 + 1e-4 * Nwet;                         %(43) P.618

            % Step 4: Calculate the effective path length using 1000 m for
            % the height of the turbulent layer hL

            hL = 1000;
            L = 2*hL / (sqrt( sind(theta).^2 + 2.35e-4 ) + sind(theta)); %(44) P.618

            % Step 5 is modified and follows Rec ITU-R P.619-5

            % The effective aperture of the earth-based antenna can be
            % estimated from its antenna gain in the direction of the path
            % using equation (48)

            Deff = 0.3 .* 10.^(0.05*Ga) ./ (pi*f);                    %(48)

            % Step 6: Calculate the antenna averaging factor

            x = 1.22 * Deff.^2 * f / L;                              %(46b)  P.618

            if (x >= 7)
                Ast = 0;
                return
            end

            g = sqrt(3.86 * (x.^2 + 1).^(11.0/12.0) * sin( 11.0/6.0 * atan(1.0/x) ) - 7.08 * x.^(5.0/6.0));   %(46a) P.618

            % Step 7: Calculate the standard deviation of the signal for
            % the applicable period and propagation path

            sigma = sigma_ref * f.^(7.0/12.0) * g / ( sind(theta).^(1.2) );          %(47) P.618

            if p <= 50

                logp = log10(p);

                ast = 2.672 - 1.258 * logp - 0.0835 * logp.^2 - 0.0597 * logp.^3;    %(49a)

                Ast = - sigma * ast;                                  %(50)

            else  % p > 50

                q = 100 - p;

                logq = log10(q);

                ast = 3.0 - 1.71 * logq + 0.072 * logq.^2 - 0.061 * logq.^3;         %(49b)

                Ast  = sigma * ast;                                   %(50)

            end

            return

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        Attachment E to Annex 1                          %
        %           Beam clearance taking atmospheric refraction into account     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Dc, Hr] = beam_clearance(obj, Ht, theta)
            %% beam_clearance Computes the ray height profile 
            %
            % Recommendation ITU-R P.619-5 Attachment E
            %
            % Inputs
            % Variable    Unit     Type     Description
            % Ht          km       float    Altitude of earth-based station (above sea level)
            % theta       deg      float    Apparent elevation angle at the earth-based station 
            %
            % Outputs:
            % Dc          km       float    Array of horizontal distances of the ray over the curved Earth 
            % Hr          km       float    Array of ray altitudes at each horizontal distance Dc
            %
            %     Rev   Date        Author                          Description
            %     -------------------------------------------------------------------------------
            %     v0    20SEP24     Ivica Stevanovic, OFCOM         Initial version

           
            if (theta <= 5)

                % Initialize
                Nmax = 2000;
                Hr_tmp = zeros(Nmax);
                Dc_tmp = zeros(Nmax);

                Hr_tmp = [Ht];          % Ray altitude at the earth-based station (51)
                Dc_tmp = [0];           % Horizontal distance at the earth-based station (52)
                epsilon = theta * pi/180;    % Ray elevation angle above local horizontal (radians) (53)

                % Set the increment in horizontal distance over curved Earth

                delta_d = 1;                  %(54)

                Re = 6371.0;

                count = 2;

                while (1)
                    % Calculate the increment in ray elevation angle:

                    delta_eps = delta_d * (1.0/Re - 4.28715e-5 * exp(-Hr_tmp(count-1)/7.348));                %(55)

                    % Re-assign the ray height
                    Hr_tmp(count) = Hr_tmp(count - 1) + delta_d*epsilon;             %(56)

                    % Re-assign the ray elevation angle:
                    epsilon = epsilon + delta_eps;              %(57)

                    % Re-assign the horizontal distance over curved Earth
                    Dc_tmp(count) = Dc_tmp(count-1) + delta_d;           %(58)

                    count = count + 1;
                    if (count > Nmax)
                        break
                    end
                    if (Hr_tmp(count-1) >= 10)
                        break
                    end

                end

                Hr = Hr_tmp(1:count-1);
                Dc = Dc_tmp(1:count-1);
            else
                Re = 6371;
                Dcmax = sqrt(Re.^2 * tand(theta).^2 + 2 * Re * (10 - Ht) ) - Re*tand(theta);
                Dc = 0 : 1 : ceil(Dcmax);  % 1 km increments
                Hr = Ht + Dc.*tand(theta) + Dc.^2/(2*Re);
            end

            return
        end

        function Nw = get_interp2_Nwet_Annual_time_location(obj, p, phie, phin)
            %% get_interp2_Nwet_Annual_time_location Interpolates the value from NWET_Annual at a given phie,phin,p
            %
            %     Recommenation iTU-R P.453-14, Section 2.2
            %
            %     Input parameters:
            %     p       -   Average year percentage (%) within the range 0.1% to 99%
            %     phie    -   Longitude, positive to east (deg) within the range -180 to 180
            %     phin    -   Latitude, positive to north (deg) within the range -90 to 90
            %
            %     Output parameters:
            %     Nw      -    Interpolated value for Nwet from the radiometeorological maps at point (phie,phin,p)
            %
            %     Rev   Date        Author                          Description
            %     -------------------------------------------------------------------------------
            %     v0    23AUG24     Ivica Stevanovic, OFCOM         Initial version

            if (p < 0.1 || p > 99)
                error('This function is defined for the average year time probabilities within the range 0.1% to 99%.');
            end

            % a) Determine two probabilities pabove and pbelow above and
            % below the desired probability from the following set

            p_set = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 60 70 80 90 95 99];
            
            kk1 = find(p >= p_set);
            kk2 = find(p <= p_set);

            kbelow = kk1(end);
            kabove = kk2(1);

            if (kbelow == kabove)

                Nw = get_interp2_Nwet_Annual_location(obj, p, phie, phin);

            else

                % b) - c) determine wet term of the surface refractivity
                % Nwetabove and Nwetbelow at the probabilities pabove and
                % pbelow using 2D interpolation as described in
                % Recommendation ITU-R P.1144
                
                pabove = p_set(kabove);
                pbelow = p_set(kbelow);

                Nwetabove = get_interp2_Nwet_Annual_location(obj, pabove, phie, phin);
                Nwetbelow = get_interp2_Nwet_Annual_location(obj, pbelow, phie, phin);

                % d) determine the wet term of the surface refractivity
                % Nwet at the desired probability p interpolating Nwetabove
                % and Nwebelow vs pabove and pbelow to p on a linear Nwet
                % vs log(p) scale
                % Recommendation does not define the basis of log
                % Here we asume log10

                Nw = Nwetbelow + (Nwetabove - Nwetbelow)*(log10(p)-log10(pbelow))/(log10(pabove)-log10(pbelow));
            end

            return

        end


            function y = get_interp2_Nwet_Annual_location(obj, cs, phie, phin)
            %% get_interp2_Nwet_Annual_location Interpolates the value from NWET_Annual_cs at a given phie,phin
            %
            %     Input parameters:
            %     cs      -   float pointing to the radiometeorological map
            %                 allowed values of time percentages
            %                 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80, 90, 95, 99
            %                 (rows-latitude: -90 to 90, columns-longitude: -180 to 180)
            %     phie    -   Longitude, positive to east (deg)
            %     phin    -   Latitude, positive to north (deg)
            %
            %     Output parameters:
            %     y      -    Interpolated value from the radiometeorological map at point (phie,phin)
            %                 as defined in Rec. ITU.R P.1144
            %
            %     Rev   Date        Author                          Description
            %     -------------------------------------------------------------------------------
            %     v0    23AUG24     Ivica Stevanovic, OFCOM         Initial version

            nr = 241;
            nc = 481;
            spacing = 0.75;

            if (phin < -90 || phin > 90)
                error ('Latitude must be within the range -90 to 90 degrees');
            end

            if (phie < -180 || phie > 180)
                error('Longitude must be within the range -180 to 180 degrees');
            end

            switch cs

                case 0.1
                    map = DigitalMaps_NWET_Annual_01();

                case 0.2
                    map = DigitalMaps_NWET_Annual_02();

                case 0.3
                    map = DigitalMaps_NWET_Annual_03();

                case 0.5
                    map = DigitalMaps_NWET_Annual_05();

                case 1
                    map = DigitalMaps_NWET_Annual_1();

                case 2
                    map = DigitalMaps_NWET_Annual_2();

                case 3
                    map = DigitalMaps_NWET_Annual_3();

                case 5
                    map = DigitalMaps_NWET_Annual_5();

                case 10
                    map = DigitalMaps_NWET_Annual_10();

                case 20
                    map = DigitalMaps_NWET_Annual_20();

                case 30
                    map = DigitalMaps_NWET_Annual_30();

                case 50
                    map = DigitalMaps_NWET_Annual_50();

                case 60
                    map = DigitalMaps_NWET_Annual_60();

                case 70
                    map = DigitalMaps_NWET_Annual_70();

                case 80
                    map = DigitalMaps_NWET_Annual_80();

                case 90
                    map = DigitalMaps_NWET_Annual_90();

                case 95
                    map = DigitalMaps_NWET_Annual_95();

                case 99
                    map = DigitalMaps_NWET_Annual_99();

            end

            latitudeOffset = phin + 90;
            longitudeOffset = phie + 180;

            latitudeIndex  = floor(latitudeOffset / spacing)  + 1;
            longitudeIndex = floor(longitudeOffset / spacing) + 1;

            latitudeFraction  = (latitudeOffset / spacing)  - (latitudeIndex  - 1);
            longitudeFraction = (longitudeOffset / spacing) - (longitudeIndex - 1);

            val_ul = map(latitudeIndex, longitudeIndex);
            val_ur = map(latitudeIndex, min(longitudeIndex + 1, nc));
            val_ll = map(min(latitudeIndex + 1, nr), longitudeIndex);
            val_lr = map(min(latitudeIndex + 1, nr), min(longitudeIndex + 1, nc));

            y1 = longitudeFraction  * ( val_ur - val_ul ) + val_ul;
            y2 = longitudeFraction  * ( val_lr - val_ll ) + val_ll;
            y  = latitudeFraction * ( y2 - y1 ) + y1;

            return
        end

    end
end
