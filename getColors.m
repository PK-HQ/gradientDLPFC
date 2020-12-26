function [colorGroup1,colorGroup2]=getColors(nColors,targetColors)
multihue=0;
scaling=1;

switch multihue
    case {0}
        greenRGB=[255 252 218
        247,252,185
        217,240,163
        173,221,142
        120,198,121
        65,171,93
        35,132,67
        0,104,55
        0,69,41]*scaling/255;

        magentaRGB=[252 238 236
        253,224,221
        252,197,192
        250,159,181
        247,104,161
        221,52,151
        174,1,126
        122,1,119
        73,0,106]*scaling/255;

        blueRGB=[240,248,255
        222,235,247
        198,219,239
        158,202,225
        107,174,214
        66,146,198
        33,113,181
        8,81,156
        8,48,107]*scaling/255;

        redRGB=[255,242,235
        254,224,210
        252,187,161
        252,146,114
        251,106,74
        239,59,44
        203,24,29
        165,15,21
        103,0,13]*scaling/255;

        [magentaOut]=upSaturation(magentaRGB,nColors);
        [greenOut]=upSaturation(greenRGB,nColors);
        [redOut]=upSaturation(redRGB,nColors);
        [blueOut]=upSaturation(blueRGB,nColors);
    
    case {1}
        %single hue, no yellow mixed in
        greenRGB=[247,252,253
        229,245,249
        204,236,230
        153,216,201
        102,194,164
        65,174,118
        35,139,69
        0,109,44
        0,68,27]*scaling/255;

        magentaRGB=[247,252,253
        224,236,244
        191,211,230
        158,188,218
        140,150,198
        140,107,177
        136,65,157
        129,15,124
        77,0,75]*scaling/255;
end


switch targetColors
    case {1}
        %target
        colorGroup1=redOut;
        colorGroup2=blueOut;
    case {0}
        %distractors
        colorGroup1=greenOut;
        colorGroup2=magentaOut;
end



%{

%oranges
colorScheme=    [    1.0000    0.9294    0.6275
    0.9961    0.8510    0.4627
    0.9961    0.6980    0.2980
    0.9922    0.5529    0.2353
    0.9882    0.3059    0.1647
    0.8902    0.1020    0.1098
    0.7412         0    0.1490
    0.5020         0    0.1490
    0.3922         0    0.0863];


green=[247,252,245
229,245,224
199,233,192
161,217,155
116,196,118
65,171,93
35,139,69
0,109,44
0,68,27
0,51,20
0,31,13]/255;

magenta=[252,251,253
239,237,245
218,218,235
188,189,220
158,154,200
128,125,186
106,81,163
84,39,143
63,0,125
50,0,101
38,0,75]/255;
%}
end