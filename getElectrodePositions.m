function [electrodeMapping]=getElectrodePositions(electrodeArrayMeasurements)
%electrodeArrayMeasurements is an array that contains electrode array ID, FMA type, rotation angle
%and x- and y-distance from the posterior tip of Principal Sulcus (see paper for description)
%i.e. {'g12', 1, 'HDFMA18', 60.6, 1.24, 2.41}

%MAKE SURE ALL MEASUREMENTS FOLLOW THE PAPER TO A TEE, DO NOT TAKE MEASUREMENTS WITH
%ARRAYS THAT ARE 'UPSIDE DOWN' (WIRE = UP-SIDE)  


for arrayNo=1:length(electrodeArrayMeasurements)
    %extract necessary details
    arrayIDStr=electrodeArrayMeasurements{arrayNo,1};
    arrayID=electrodeArrayMeasurements{arrayNo,2};
    arrayType=electrodeArrayMeasurements{arrayNo,3};
    arrayRotAngle=electrodeArrayMeasurements{arrayNo,4};
    arrayXtranslate=electrodeArrayMeasurements{arrayNo,5};
    arrayYtranslate=electrodeArrayMeasurements{arrayNo,6};
    
    %rotate and translate
    [processedArray]=rotateTranslateArray(arrayID,0,arrayType,arrayRotAngle,arrayXtranslate,arrayYtranslate);
    %remove ground electrodes and re-label electrodes from 1:16 or 1:32
    [processedArray]=removeGroundElectrodes(processedArray);
    %add the array ID str 'g12' to every electrode for assigning neuron AP position later
    processedArray=[repmat({arrayIDStr},length(processedArray),1),processedArray];
    %concat into electrode mapping array, order of array does not matter
    electrodeMapping=[electrodeMapping;processedArray];
end
end