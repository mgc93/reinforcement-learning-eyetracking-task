function S = dlLoadStimuli2

%--------------------------------------------------------------------------

% move to proper image directory
PWD = pwd;
cd('set_2');

    F(1,1) = dir('*1.png');
    F(2,1) = dir('*2.png');
    F(3,1) = dir('*3.png');
    F(4,1) = dir('*4.png');
    F(5,1) = dir('*5.png');
    F(6,1) = dir('*6.png');

    
    for i = 1:size(F,1)    
        S(i,1).id = [num2str(i,'%.0f')];
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