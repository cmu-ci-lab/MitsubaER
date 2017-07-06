clc, clear, close all
rng(1)
% fileName = '/mnt/da64b98f-9fd9-4b2c-994e-ca7276846901/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/testVol.vol';
fileName = 'testVol.vol';

bbmin = [0 0 0];
bbmax = [1 1 1];

data = rand(20, 30, 23); % the convention is y X x X z
writeGridToVol(data, bbmin, bbmax, fileName);

[data1, bbmin1, bbmax1] = readGridToVol(fileName);

sum( abs(data - data1), 'all')
sum( abs(bbmin - bbmin1), 'all')
sum( abs(bbmax - bbmax1), 'all')

%Display to validate
% s = '';
% for z = 1:size(data, 3)
%     for y = 1:size(data, 2)
%         for x = 1:size(data, 1)
%             s = s + sprintf("%f, ", data(x, y, z));
%         end
%         s = s + sprintf("\n ");
%     end
%     s = s + sprintf("\n ");
% end
% s

% spline coefficients to validate
% addpath('/mnt/da64b98f-9fd9-4b2c-994e-ca7276846901/Dropbox/AccoustoOptics+InvRendering/CodeEtc/SkeletalRenderer/ercrdr_angletracingNEE/splinepractice');
% 
% S = Spline(3, bbmin, bbmax, size(data));
% S = S.build(data);
% s = '';
% for z = 1:size(data, 3)
%     for y = 1:size(data, 2)
%         for x = 1:size(data, 1)
%             s = s + sprintf("%f, ", S.coeff(x, y, z));
%         end
%         s = s + sprintf("\n ");
%     end
%     s = s + sprintf("\n ");
% end
% s