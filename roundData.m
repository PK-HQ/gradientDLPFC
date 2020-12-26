function roundedData=roundData(data,roundNearestVal)
switch roundNearestVal
    case {.1}
        roundedData=round(data,1);
    case {.25}
        %round up
        roundedData = ceil(data * 4) / 4;

        %round down
        %roundedData = floor(data * 4) / 4;
    case {.5}
        roundedData = round(data*2)/2;
        %{
        x = 16.625; 
        dist = mod(x, 0.5); 
        floorVal = x - dist; 
        newVal = floorVal; 
        if dist >= 0.25, newVal = newVal + 0.5; end
        %}
    case {1}
        roundedData=round(data,0);
    case {2}
        roundedData=2*floor(data/2);
    case {5}
        roundedData= round(data/ 5 ) * 5;
    case {10}
        roundedData= ceil(data/10)*10;
        
    otherwise
        roundedData=roundNearestVal*floor(data/roundNearestVal);
end

if ~isequal(size(roundedData),size(data))
    fprintf('Warning: size of output is different from size of data')
end
        
        
end