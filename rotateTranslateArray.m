function [rotatedTranslatedArray]=rotateTranslateArray(arrayNo,arrayAB,arrayType,rotAngle,xTranslation,yTranslation)
%% Rotates and translates electrodes on array according to measurements
% Array No is assigned by experimenter (1-N) 
% Array AB can be ignored if not assigned by experimenter (0/1/2=Null/A/B)
% arrayXY is the array template
% rotAngle=angle between array base and intersect (see paper)
% xTranslate/yTranslate=x-/y-coord from reference point


%% MAKE TEMPLATE ELECTRODE ARRAYS FOR THE 3 TYPES USED IN THE EXPERIMENT
% (36-CHANNEL FMA, 18-CHANNEL FMA, 18-CHANNEL HD-FMA)
%Make a template array for each electrode array, containing the coordinates of each
%electrode within it. Values are with respect to a reference point at (0,0), and 0 rotation angle. 
%This is translated and rotated below, individually for each array
[FMA36,FMA18,HDFMA18]=makeArrayCoordinates();

if strcmp(arrayType,'FMA36')
    arrayTemplate=FMA36;
elseif  strcmp(arrayType,'FMA18')
    arrayTemplate=FMA18;
elseif  strcmp(arrayType,'HDFMA18')
    arrayTemplate=HDFMA18;
end


%% ROTATE AND TRANSLATE ARRAYS
% rotate by degrees clockwise around (0,0)
electrodeID=arrayTemplate(:,1);x=arrayTemplate(:,2);y=arrayTemplate(:,3);
radAngle = deg2rad(rotAngle);
xRot     = x*cos(radAngle) - y*sin(radAngle);
yRot     = x*sin(radAngle) + y*cos(radAngle);

% and translate by (xTranslation yTranslation)
xRotShift = xRot + xTranslation;
yRotShift = yRot + yTranslation;

rotatedTranslatedArray=[repmat(arrayNo,size(arrayTemplate,1),1),repmat(arrayAB,size(arrayTemplate,1),1),electrodeID,xRotShift,yRotShift];
end