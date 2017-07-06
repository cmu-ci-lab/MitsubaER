clc, clear, close all

% Creates RIF from Signed Distance Metric (SD)

% Strategy: 
% 1. n (SD <= 0) = 1 
% 2. highest(h) will have nmax. 
% 3. Inbetween will go as ^r. 
% 4. n = 1 + k*d^r; k = (n_max - 1)/h^r
%will have 

nmax = 1.50;
nmin = 1.10;
rs = [1 1.5 2 3 10]; %linear
visualizeRIF = 0;

%Box
% inFileName = '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxSD.vol';
% outFileNamePrefix = '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/meshes/BoxSDRIF';

%Lucy
% inFileName = '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/Lucy/lucySD.vol';
% outFileNamePrefix = '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/Lucy/lucySDRIF';
%Dragon
inFileName = '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/Dragon/dragonSD.vol';
outFileNamePrefix = '/home/igkiou/Dropbox/BlenderXMLFiles/Eikonal_ToG/Dragon/dragonSDRIF';

[data, bbmin, bbmax] = readVolToGrid(inFileName);

data = -data; %For dragon, SDF is flipped. 

data(data < 0) = 0;
h = max(data(:));
for r = rs
    k = (nmax - nmin)/h^r;
    RIF = nmin + k*data.^r;
    outFileName = strcat(outFileNamePrefix, '_', num2str(r), '.vol');
    writeGridToVol(RIF, bbmin, bbmax, outFileName);
    
    %Visualize RIF
    if(visualizeRIF)
        figure, 
        for i=1:size(RIF, 3) 
            imagesc(RIF(:,:, i)); 
            drawnow;
            pause(0.1);
        end
    end
    
end


