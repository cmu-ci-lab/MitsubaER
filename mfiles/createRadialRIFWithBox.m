
clc, clear, close all

%grow linearly from 1.3 to 2.0


bbmin = [-225, -225,  25];
bbmax = [ 225,  225, 125];
bbres = [ 226,  226,  51];

bbstride = (bbmax - bbmin)./(bbres - 1);

data = zeros(bbres);

center = (bbmax + bbmin)/2;
R  = max( norm(bbmax-center), norm(bbmin-center));
falloff = @(r, R) (2 - (r/R).^2);

for i=1:size(data,1)
    for j=1:size(data,2)
        for k=1:size(data,3)
            P = bbmin + bbstride.*([i j k] - 1);
            data(i, j, k) = falloff(norm(P - center), R);
        end
    end
end

% for k=1:size(data, 3)
%     imagesc(squeeze(data(:, :, k)));
%     caxis([min(data(:)) max(data(:))]);
%     pause(0.1);
% end

writeGridToVol(data, bbmin, bbmax, '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxRIF_Radial.vol');

% cleanRIF('/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxRIFUnclean.vol', ...
%          '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxWN.vol', ...
%          '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxRIF.vol');
%          
