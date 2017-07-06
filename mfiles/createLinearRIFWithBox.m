
clc, clear, close all

%grow linearly from 1.3 to 2.0

nmax = 1.6;
nmin = 1.3;

bbmin = [-225, -225,  25];
bbmax = [ 225,  225, 125];
bbres = [ 226,  226,  51];

bbstride = (bbmax - bbmin)./(bbres - 1);

data = zeros(bbres);

for i=1:size(data,1)
    for j=1:size(data,2)
        for k=1:size(data,3)
            data(i, j, k) = nmin + (nmax - nmin)/(size(data,2) - 1)*(j-1);
        end
    end
end

writeGridToVol(data, bbmin, bbmax, '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxRIF.vol');

% cleanRIF('/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxRIFUnclean.vol', ...
%          '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxWN.vol', ...
%          '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxRIF.vol');
%          
