clc, clear, close all

M = zeros(100, 100, 100);

fid = fopen('~/Dropbox/BlenderXMLFiles/Eikonal_ToG/LaserScatteringMediumScene/delete.txt', 'r');

xmin = -2;
ymin = -2;
zmin = -2;

xmax =  3;
ymax =  3;
zmax =  3;

xres =  100;
yres =  100;
zres =  100;

x = linspace(xmin, xmax, xres);
y = linspace(ymin, ymax, yres);
z = linspace(zmin, zmax, zres);

[X, Y, Z] = meshgrid(x, y, z);

index= 1;
tline = fgetl(fid);
j = 1;
k = 1;
while(ischar(tline))
    if(~isempty(tline))
        C = strsplit(tline, ',');
        M(:, j, k) = str2double(C(1:100));
        j = j + 1;
    else
        j = 1;
        k = k+1;
    end
    tline = fgetl(fid);
%     s = fscanf(fid, '%f');
%     if(s ~= [])
%         M(index) = s;index = index + 1;
%     end
end

fclose(fid);

figure, scatter3(X(:), Y(:), Z(:), abs(M(:)))

J = [X(:) Y(:) Z(:) M(:)];

J(abs(J(:, 4)) < .5, :) = [];

figure, plot3(J(:, 1), J(:, 2), J(:, 3), 'r.');
