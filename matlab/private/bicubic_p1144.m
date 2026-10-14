function v = bicubic_p1144(G, r, c)
%BICUBIC_P1144 Bicubic interpolation of Rec. ITU-R P.1144 Annex 1
%   v = bicubic_p1144(G, r, c) interpolates the grid G at the fractional
%   0-based row and column indices r, c using the 4 x 4 neighbouring points
%   and the kernel K(d) with a = -0.5. Indices outside the grid are clamped.

R = floor(r);
C = floor(c);
[nr, nc] = size(G);
v = 0;
for i = R-1:R+2
    ii = min(max(i, 0), nr-1) + 1;
    row = 0;
    for j = C-1:C+2
        jj = min(max(j, 0), nc-1) + 1;
        row = row + G(ii, jj) * kern(c - j);
    end
    v = v + row * kern(r - i);
end
end

function k = kern(d)
d = abs(d);
if (d <= 1)
    k = 1.5*d^3 - 2.5*d^2 + 1;
elseif (d <= 2)
    k = -0.5*d^3 + 2.5*d^2 - 4*d + 2;
else
    k = 0;
end
end
