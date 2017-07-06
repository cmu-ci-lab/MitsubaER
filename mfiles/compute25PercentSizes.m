clc, clear, close all

originalMin = [-141.445 11.9041 -257.147];
originalMax = [153.547 866.509 205.833];

boundary = (originalMax - originalMin)*0.125;

bbMin = originalMin - boundary
bbMax = originalMax + boundary

(bbMax - bbMin)./boundary; % Should be 10.0

smallestRes = 100;
%bbRes heuristic: Get atleast 100 in smallest dimension
stride = min((bbMax - bbMin)/(smallestRes-1));
bbRes = round( (bbMax - bbMin)/stride + 1)
(bbMax - bbMin)./(bbRes-1); %Should be roughly equal
