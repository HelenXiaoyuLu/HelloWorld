%% JEDI-2P
clc
% val1 = [0.753342516 0.687270608 0.739410633];
% val2 = [0.651811054 0.65775464];
% val3 = [0.588335662 0.549945044];
% [T1,p1] = ttest2(val1,val2) 
% [T2,p2] = ttest2(val2,val3) 
% [T3,p3] = ttest2(val1,val3) 
val1 = [423.2727341 346.4051769 397.7494332];
val2 = [578.372199 599.7080706];
val3 = [816.6399043 742.9418041];
[T1,p1] = ttest2(val1,val2) 
[T2,p2] = ttest2(val2,val3) 
[T3,p3] = ttest2(val1,val3) 

%% JEDI-1P
clc
% val1 = [0.474827498 0.481475942 0.478065992];
% val2 = [0.445763264 0.451398946];
% val3 = [0.436900936 0.447853278];
% [T1,p1] = ttest2(val1,val2) 
% [T2,p2] = ttest2(val2,val3) 
% [T3,p3] = ttest2(val1,val3) 
val1 = [301.4809999 304.0600725 332.1240396];
val2 = [484.8334172 583.005459];
val3 = [809.2105801 805.2165921];
[T1,p1] = ttest2(val1,val2) 
[T2,p2] = ttest2(val2,val3) 
[T3,p3] = ttest2(val1,val3) 
%% ASAP2s
clc
% val1 = [0.402289546 0.409220607];
% val2 = [0.404025655 0.409467765];
% val3 = [0.404009958 0.406910739];
% [T1,p1] = ttest2(val1,val2) 
% [T2,p2] = ttest2(val2,val3) 
% [T3,p3] = ttest2(val1,val3) 
val1 = [272.3143157 296.2886613];
val2 = [407.2029001 460.0403275];
val3 = [673.6645448 450.858299];
[T1,p1] = ttest2(val1,val2) 
[T2,p2] = ttest2(val2,val3) 
[T3,p3] = ttest2(val1,val3) 
%% dpOPT-CAAX
clc
% val1 = [0.34020364 0.353250514];
% val2 = [0.375831705 0.395630121];
% val3 = [0.406183665 0.41364425];
% [T1,p1] = ttest2(val1,val2) 
% [T2,p2] = ttest2(val2,val3) 
% [T3,p3] = ttest2(val1,val3) 
val1 = [2454.699588 2769.971792];
val2 = [3399.680189 3034.149106];
val3 = [3341.067573 2357.014158];
[T1,p1] = ttest2(val1,val2) 
[T2,p2] = ttest2(val2,val3) 
[T3,p3] = ttest2(val1,val3) 

%% cpOPT-CAAX
clc
% val1 = [0.896628012 0.886259381];
% val2 = [0.87557995 0.865727992];
% val3 = [0.847315947 0.841322633];
% [T1,p1] = ttest2(val1,val2) 
% [T2,p2] = ttest2(val2,val3) 
% [T3,p3] = ttest2(val1,val3) 
val1 = [1155.989599 1272.605666];
val2 = [1550.473744 2008.000166];
val3 = [1617.664131 1545.674986];
[T1,p1] = ttest2(val1,val2) 
[T2,p2] = ttest2(val2,val3) 
[T3,p3] = ttest2(val1,val3) 

%% lyn-EGFP
val1 = [5964.764204 5251.336944];
val2 = [7588.894513 6666.371524];
val3 = [5642.700138 5896.131325];
[T1,p1] = ttest2(val1,val2) 
[T2,p2] = ttest2(val2,val3) 
[T3,p3] = ttest2(val1,val3) 

%% EGFP-CAAX
val1 = [6583.218282 5685.538499];
val2 = [8178.471079 7074.967389];
val3 = [5867.094487 5917.439038];
[T1,p1] = ttest2(val1,val2) 
[T2,p2] = ttest2(val2,val3) 
[T3,p3] = ttest2(val1,val3) 