fGHz = 10;
dkm = 50;
e2sflag = true;
p1 = 50;
He = 0; % Altitude of the earth station antenna
Hs = 100; % Altitude of the space station antenna (it is beyond 100 km)
phi_e = 90; % elevation of the earth station antenna main beam - towards zenith)
phi_s = 90; % elevation of the space station antenna main beam - toward earth
Dphi = 2; % Beamwidth of the space station antenna - as they point to one another, this has no effect
std_atm = true; % we use the standard atmospheres from P.835
atm_type = 1;   % we use the standard atmosphere from Section 1.
% these are arbitrary and not connected to the parameters aboveChère 
lat_s = 0; %geostationary
lat_e = 45;
lon_s = 10;
lon_e = 9;
p2 = 50;
Ga = 20;


p = P619;

tl = p.tl_p619_single(fGHz, dkm, e2sflag, p1, He, Hs, phi_e, phi_s, Dphi, atm_type, lat_s, lat_e, lon_s, lon_e, p2, Ga)