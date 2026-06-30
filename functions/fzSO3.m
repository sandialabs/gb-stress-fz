%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental zones in SO(3) parametrized by the Rodrigues vector
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Examples on how to create these FZ given at the end.
%%% Inputs:
% cs1,cs2   Point group symmetries (Oh,D8h,D6h,D4h,D2h,D3d,C2h,Ci)
%
%%% Outputs:
% V         Vertices (Rodrigues vector parametrization)
% F         Faces (vertex connectivity)
% thmax     Angle the FZ spans
% 

function [V,F,thmax] = fzSO3(cs1,cs2) 
    tanpi16 = tan((pi/8)/2);      % 0.198912
    tanpi12 = tan((pi/6)/2);      % 0.267949
    tanpi08 = tan((pi/4)/2);      % 0.414214
    mtanpi12 = 1-tanpi12;         % 0.732051
    m2tanpi8 = 1-2*tanpi08;       % 0.171573

    if strcmp(cs1,'Oh') && strcmp(cs2,'Oh')
        V = [       0        0        0             % O
              tanpi08        0        0             % A
              tanpi08  tanpi08        0             % B
                  1/3      1/3      1/3             % E
              tanpi08  tanpi08 m2tanpi8             % C
              tanpi08 1/2-tanpi08/2 1/2-tanpi08/2]; % D
        F = {[1  2  3];                             % OAB
             [1  3  5  4];                          % OBCE
             [1  2  6  4];                          % OADE
             [2  3  5  6];                          % ABCD
             [5  6  4]};                            % CDE
        thmax = pi/4;

    elseif (strcmp(cs1,'Oh') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'Oh'))
        V = [   -tanpi08   tanpi08  m2tanpi8
               -m2tanpi8   tanpi08   tanpi08
                -tanpi08  m2tanpi8   tanpi08
                 tanpi08  m2tanpi8   tanpi08
                 tanpi08   tanpi08  m2tanpi8
                m2tanpi8   tanpi08   tanpi08
                -tanpi08   tanpi08         0
                -tanpi08         0   tanpi08
                 tanpi08         0   tanpi08
                 tanpi08   tanpi08         0
                -tanpi08         0         0
                 tanpi08         0         0];
        F = {[5  10  12   9   4];
             [1   7  10   5   6   2];
             [2  6  4  9  8  3];
             [3   8  11   7   1];
             [4  6  5];
             [1  2  3]};
        thmax = pi;

    elseif (strcmp(cs1,'D8h') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'D8h'))
        V = [                -1                      0                      0
                             -1      0.198912367379658                      0
            -0.8477590650225735     0.5664544973505214                      0
            -0.5664544973505214     0.8477590650225735                      0
             -0.198912367379658                      1                      0
              0.198912367379658                      1                      0
             0.5664544973505214     0.8477590650225735                      0
             0.8477590650225735     0.5664544973505214                      0
                              1      0.198912367379658                      0
                              1                      0                      0
                             -1                      0      0.198912367379658
                             -1      0.198912367379658      0.198912367379658
            -0.8477590650225735     0.5664544973505214      0.198912367379658
            -0.5664544973505214     0.8477590650225735      0.198912367379658
             -0.198912367379658                      1      0.198912367379658
              0.198912367379658                      1      0.198912367379658
             0.5664544973505214     0.8477590650225735      0.198912367379658
             0.8477590650225735     0.5664544973505214      0.198912367379658
                              1      0.198912367379658      0.198912367379658
                              1                      0      0.198912367379658];
        F = {[1 2 12 11];
             [2 3 13 12];
             [3 4 14 13];
             [4 5 15 14];
             [5 6 16 15];
             [6 7 17 16];
             [7 8 18 17];
             [8 9 19 18];
             [9 10 20 19];
             [ 1  2  3  4  5  6  7  8  9 10];
             [11 12 13 14 15 16 17 18 19 20];
             [1 10 20 11]};
        thmax = pi;

    elseif (strcmp(cs1,'D6h') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'D6h'))
        V = [   tanpi12         1   tanpi12
                      1   tanpi12   tanpi12
                     -1   tanpi12   tanpi12
               mtanpi12  mtanpi12   tanpi12
              -mtanpi12  mtanpi12   tanpi12
               -tanpi12         1   tanpi12
               -tanpi12         1         0
                tanpi12         1         0
                      1         0   tanpi12
                      1   tanpi12         0
                     -1         0   tanpi12
                     -1  tanpi12          0
              -mtanpi12  mtanpi12         0
               mtanpi12  mtanpi12         0
                     -1         0         0
                      1         0         0];
        F = {[6   1   4   2   9  11   3   5];
             [10  16   9   2];
             [2   4  14  10];
             [4   1   8  14];
             [7  8  1  6];
             [6   5  13   7];
             [5   3  12  13];
             [3  11  15  12];
             [8   7  13  12  15  16  10  14];
             [16  15  11   9]};
        thmax = pi;

    elseif (strcmp(cs1,'D4h') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'D4h'))
        V = [ -tanpi08       1 tanpi08
               tanpi08       1 tanpi08
                     1 tanpi08 tanpi08
                    -1 tanpi08 tanpi08
                    -1 tanpi08       0
                     1       0 tanpi08
                    -1       0 tanpi08
              -tanpi08       1       0
               tanpi08       1       0
                     1 tanpi08       0
                    -1       0       0
                     1       0       0];
        F = {[2  3  6  7  4  1];
             [10  12   6   3];
             [3   2   9  10];
             [2  1  8  9];
             [5  8  1  4];
             [4   7  11   5];
             [5  11  12  10   9   8];
             [12  11   7   6]};
        thmax = pi;

    elseif strcmp(cs1,'D2h') && strcmp(cs2,'D2h')
        V = [  -1      1      1
                1      1      1
               -1      1      0
               -1      0      1
                1      1      0
                1      0      1
               -1      0      0
                1      0      0];
        F = {[5  8  6  2];
             [2  1  3  5];
             [2  6  4  1];
             [4  7  3  1];
             [5  3  7  8];
             [8  7  4  6]};
        thmax = pi;

    elseif (strcmp(cs1,'D3d') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'D3d'))
        V = [        -1   tanpi12   tanpi12
                     -1  -tanpi12   tanpi12
              -mtanpi12 -mtanpi12   tanpi12
               -tanpi12         1   tanpi12
               -tanpi12        -1   tanpi12
                tanpi12         1   tanpi12
                tanpi12        -1   tanpi12
               mtanpi12  mtanpi12   tanpi12
                      1   tanpi12   tanpi12
                      1  -tanpi12   tanpi12
              -mtanpi12  mtanpi12   tanpi12
               mtanpi12 -mtanpi12   tanpi12
                     -1  -tanpi12         0
                     -1   tanpi12         0
              -mtanpi12  mtanpi12         0
              -mtanpi12 -mtanpi12         0
               -tanpi12         1         0
               -tanpi12        -1         0
                tanpi12         1         0
                tanpi12        -1         0
               mtanpi12  mtanpi12         0
               mtanpi12 -mtanpi12         0
                      1  -tanpi12         0
                      1   tanpi12         0];
        F = {[12   7   5   3   2   1  11   4   6   8   9  10];
             [24  23  10   9];
             [23  22  12  10];
             [9   8  21  24];
             [7  12  22  20];
             [8   6  19  21];
             [20  18   5   7];
             [17  19   6   4];
             [3   5  18  16];
             [4  11  15  17];
             [2   3  16  13];
             [11   1  14  15];
             [13  14   1   2];
             [15  14  13  16  18  20  22  23  24  21  19  17]};
        thmax = 2*pi;

    elseif (strcmp(cs1,'C2h') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'C2h'))
        V = [  -1      1      1
                1      1      1
                1      1     -1
               -1      1     -1
               -1      0      1
               -1      0     -1
                1      0      1
                1      0     -1];
        F = {[2  3  8  7];
             [7  5  1  2];
             [4  6  8  3];
             [4  3  2  1];
             [1  5  6  4];
             [5  7  8  6]};
        thmax = pi;

    elseif (strcmp(cs1,'C2h*') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'C2h*'))
        V = [   1      1      1
                1     -1      1
                1     -1     -1
                1      1     -1
                0      1      1
                0      1     -1
                0     -1      1
                0     -1     -1];
        F = {[2  3  8  7];
             [7  5  1  2];
             [4  6  8  3];
             [4  3  2  1];
             [1  5  6  4];
             [5  7  8  6]};
        thmax = pi;

    elseif (strcmp(cs1,'Ci') && strcmp(cs2,'D2h')) || (strcmp(cs1,'D2h') && strcmp(cs2,'Ci'))
        V = [  -1 -1 -1
               -1 -1  1
               -1  1 -1
               -1  1  1
                1 -1 -1
                1 -1  1
                1  1 -1
                1  1  1];
        F = {[6  8  7  5];
             [3  7  8  4];
             [6  2  4  8];
             [3  1  5  7];
             [5  1  2  6];
             [4  2  1  3]};
        thmax = 2*pi;
    
    else
        error(['The fundamental zone ',cs1,'-',cs2,' is not included.'])
    end

end



%%% Scripts to generate the fundamental zones
%
% CASE 1: from MTEX; for example (D4h-D2h) as:
% cs1 = crystalSymmetry('D4h');
% cs2 = crystalSymmetry('D2h');
% oR = fundamentalRegion(cs1,cs2);
% figure; plot(oR,'Color','r','rodrigues');
% V = Rodrigues(oR.V)
% disp('F = {')
% for i = 1:length(oR.F)
%     disp(['[',num2str([oR.F{i}']),'];'])
% end
% disp('};')
%
%
% CASE 2: from polygon for (Dnh-D2h) fz; for example:
% n = 8;
% zmax = tan((pi/n)/2);
% Vpoly = nsidedpoly(2*n);
% Vpoly = Vpoly.Vertices/max(abs(Vpoly.Vertices),[],'all');
% Vpoly(Vpoly(:,2)<0,:) = [];
% Vpoly = [-1,0;...
%          Vpoly;...
%           1,0];
% Vpoly = [Vpoly;Vpoly];
% Vpoly = [Vpoly,[0*ones(size(Vpoly,1)/2,1);zmax*ones(size(Vpoly,1)/2,1)]];
% figure
% scatter3(Vpoly(:,1),Vpoly(:,2),Vpoly(:,3),10,'b','filled')
% axis equal
%


