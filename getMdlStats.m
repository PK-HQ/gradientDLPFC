function [R2Adj_whole,R2Adj_perseg,RMSE_whole,RMSE_perseg]=getMdlStats(mdlA,mdlB,mdlC)
numModels=sum(~isempty(mdlA)+~isempty(mdlB)+~isempty(mdlC));
switch numModels
    case {1}
        %calc rmse
        RMSE_whole=sqrt((mdlA.SSE)/(mdlA.DFE));
        RMSE_perseg=[mdlA.RMSE];

        %calc adj r2
        sumN=mdlA.NumObservations;
        sumK=mdlA.NumCoefficients;
        sumSSE=mdlA.SSE;
        sumSST=mdlA.SST;
        R2Adj_whole=1-((sumN-1)*(sumSSE))/((sumN-sumK)*(sumSST));
        R2Adj_perseg=[mdlA.Rsquared.Adjusted];
        
    case {2}
        %calc rmse
        RMSE_whole=sqrt((mdlA.SSE + mdlB.SSE)/(mdlA.DFE + mdlB.DFE));
        RMSE_perseg=[mdlA.RMSE mdlB.RMSE];

        %calc adj r2
        sumN=mdlA.NumObservations+mdlB.NumObservations;
        sumK=mdlA.NumCoefficients+mdlB.NumCoefficients;
        sumSSE=mdlA.SSE+mdlB.SSE;
        sumSST=mdlA.SST+mdlB.SST;
        R2Adj_whole=1-((sumN-1)*(sumSSE))/((sumN-sumK)*(sumSST));
        R2Adj_perseg=[mdlA.Rsquared.Adjusted mdlB.Rsquared.Adjusted];
        
    case {3}
        %calc rmse
        RMSE_whole=sqrt((mdlA.SSE + mdlB.SSE + mdlC.SSE)/(mdlA.DFE + mdlB.DFE + mdlC.DFE));
        RMSE_perseg=[mdlA.RMSE mdlB.RMSE mdlC.RMSE];

        %calc adj r2
        sumN=mdlA.NumObservations+mdlB.NumObservations+mdlC.NumObservations;
        sumK=mdlA.NumCoefficients+mdlB.NumCoefficients+mdlC.NumCoefficients;
        sumSSE=mdlA.SSE+mdlB.SSE+mdlC.SSE;
        sumSST=mdlA.SST+mdlB.SST+mdlC.SST;
        R2Adj_whole=1-((sumN-1)*(sumSSE))/((sumN-sumK)*(sumSST));
        R2Adj_perseg=[mdlA.Rsquared.Adjusted mdlB.Rsquared.Adjusted mdlC.Rsquared.Adjusted];
        
end