function S = dlLoadStimuli

%--------------------------------------------------------------------------

% move to proper image directory
PWD = pwd;
cd('set_1');

    F(1,1) = dir('*_01.png');
    F(2,1) = dir('*_02.png');
    F(3,1) = dir('*_03.png');
    F(4,1) = dir('*_04.png');
    F(5,1) = dir('*_05.png');
    F(6,1) = dir('*_06.png');

    
    for i = 1:size(F,1)    
        S(i,1).id = [num2str(i,'%02.f')];
        S(i,1).imgfile = F(i,1).name;
    end
    
    % load images
    for i=1:size(F,1)
        S(i,1).img = imread(S(i,1).imgfile);
    end
    
    cd('../');

% end

% return to original directory
cd(PWD);