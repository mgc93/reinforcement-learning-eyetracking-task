%% initialize psychtoolbox

% ctrl +c to stop the experiment
% sca for closing the window

%--------------------------------------------------------------------------
% Initial Setup  
%--------------------------------------------------------------------------


% Clear the workspace and the screen
close all;
clear;

% enter subject data
rng('shuffle');
sesInfo = inputdlg({'Subject ID:'},'Subject Details',1,{'1'});
sID     = str2num(sesInfo{1});
sCondition    = 1;

if mod(sID,2)~=0
    mode = 1;
elseif mod(sID,2)==0
    mode = 2;
end
  

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);
%screenNumber = 1;

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray1 = GrayIndex(screenNumber);
% different gray (lighter/silver)
gray = [192,192,192];

% screen and text color
backgroundColor = gray1;
textColor = white;
fixColor = white;
numberColor = white;    
frameColor = [1 0 0];

%black = gray1;

% Open an on screen window
%[window, windowRect] = PsychImaging('OpenWindow', screenNumber, backgroundColor, [0 0 800 600]);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, backgroundColor, []);




Screen('TextFont', window,'Arial');
Screen('TextSize', window,40);
Screen('TextStyle', window,1);



%--------------------------------------------------------------------------
% Welcome Screen
%--------------------------------------------------------------------------

HideCursor();

DrawFormattedText(window,['Welcome to the experiment! Press any key to start.'],'center','center',textColor);
Screen('Flip',window);
KbWait([], 2);
KbWait([], 1);

%% define coordinates

%--------------------------------------------------------------------------
% Define Rectangular Coordinates
%--------------------------------------------------------------------------

globalRect = windowRect;
[midX,midY] = RectCenter(windowRect);

% length of the small boxes where images are placed
% yLength = globalRect(4)/8;
% xLength = globalRect(3)/8;
% 
% xStart = [2*xLength, 5*xLength, midX-xLength/2, midX-xLength/2];
% yStart = [3*yLength,  3*yLength, 2*yLength, 4*yLength];

yLength = globalRect(4)/6;
xLength = globalRect(3)/6;

xStart = [xLength, 4*xLength, midX-xLength/2, midX-xLength/2];
yStart = [3*yLength-yLength/2, 3*yLength-yLength/2, yLength, 4*yLength];

% choice rectangles
rectC(1,:) = [xStart(1), yStart(1), xStart(1)+xLength, yStart(1)+yLength];
rectC(2,:) = [xStart(2), yStart(2), xStart(2)+xLength, yStart(2)+yLength];
% feedback rectangles
rectC(3,:) = [xStart(3), yStart(3), xStart(3)+xLength, yStart(3)+yLength];
rectC(4,:) = [xStart(4), yStart(4), xStart(4)+xLength, yStart(4)+yLength];


% fixation cross parameters
fixWidth    = 3.75;
fixHeight   = 30;

FixCross    = [midX-fixWidth,midY-fixHeight,midX+fixWidth,midY+fixHeight;
    midX-fixHeight,midY-fixWidth,midX+fixHeight,midY+fixWidth];

% placement of the number payoffs text
xposN = [midX, midX];
yposN = [yStart(3)+yLength/2-yLength/8, yStart(4)+yLength/2-yLength/8];

xText = [midX, xLength, xLength];
yText = [2*yLength, 2*yLength+yLength/2, 4*yLength+yLength/2];

%--------------------------------------------------------------------------
% Trial timing matrix
%--------------------------------------------------------------------------


% get timing information for 100 simulations  
[isi_t1_1, isi_t2_1, isi_t1_2, isi_t2_2, eff_1a_1, eff_2a_1,eff_all_1, eff_1a_2, eff_2a_2,eff_all_2] = createTimingMatrix();



%% part 1



%--------------------------------------------------------------------------
% Trial matrix
%--------------------------------------------------------------------------

% associate each square number with a color
% cyan        = [0.2 0.8 0.8];
% brown       = [0.2 0 0];
% orange      = [1 0.5 0];
% blue        = [0 0.5 1];
% green       = [0 0.6 0.3];
% red         = [1 0.2 0.2];
% purple      = [0.5 0 0.5];

% colorSquares = [0.5 0 0.5;
%     0 0.5 1;
%     0.2 0.8 0.8;
%     0 0.6 0.3;
%     1 0.5 0;
%     1 0.2 0.2];
% 
% colorNames = ['purple';
%     'blue';
%     'cyan';
%     'green';
%     'orange';
%     'red'];
% 
% % permute the rows for every subject
% [m, n] = size(colorSquares);
% rColor = randperm(m);
% colorSquares = colorSquares(rColor,:);
% colorNames = colorNames(rColor,:);


% magenta = [1 0 1];
% yellow = [1 0.8 0];
% turquoise = [0 0.4 0.4];
% brown = [0.2 0 0];
% light green = [0 1 0.2];
% pink = [0.8 0 0.4];

% colorSquares = [1 0 1;
%     1 0.8 0;
%     0 0.4 0.4;  
%     0.2 0 0;
%     0 1 0.2;
%     0.8 0 0.4];

% define payoffs for each square color

% randomize order of sets/payoffs

if(mod(sID,2) == 1)
    ord = 0;
else
    ord = 1;
end

if(sCondition == 1)
    if(ord==0)
        
        squaresPay1 = [11,7;
            6,10;
            9,5;
            4,8;
            7,3;
            2,6];
    else
        squaresPay1 = [7,11;
            10,6;
            5,9;
            8,4;
            3,7;
            6,2];
    end
    
else
    disp('NOT A VALID CONDITION');
%     if(ord<=0.5)
%         squaresPay1 = [11,7;
%             7,11;
%             9,5;
%             5,9;
%             7,3;
%             3,7];
%     else
%         squaresPay1 = [6,10;
%             10,6;
%             4,8;
%             8,4;
%             2,6;
%             6,2];
%     end
end


% permute payoffs for each subjects
[p1, q1] = size(squaresPay1);
squaresPay1 = squaresPay1(randperm(p1),:);

% trial matrix with the choice options
trialMat1 = createTrialMatrix();

% add noise payoffs
payTrueLeft1 = squaresPay1(trialMat1(:,1),:);
payTrueRight1 = squaresPay1(trialMat1(:,2),:);

payNoiseLeft1 = randi([0,3],size(trialMat1,1),2);
payNoiseRight1 = randi([0,3],size(trialMat1,1),2);

payTotalLeft1 = payTrueLeft1 + payNoiseLeft1;
payTotalRight1 = payTrueRight1 + payNoiseRight1;

payTotal1 = [payTotalLeft1, payTotalRight1];


% create trial matrix for the colors, one number/string for each color
% and for the payoffs
for trial = 1:size(trialMat1,1)
%     colorTrialMat(trial,:) = [colorSquares(trialMat(trial,1),:), colorSquares(trialMat(trial,2),:)];
%     colorNamesTrialMat(trial,:) =  [colorNames(trialMat(trial,1),:), colorNames(trialMat(trial,2),:)];
    payTrialMat1(trial,:) = [squaresPay1(trialMat1(trial,1),:), squaresPay1(trialMat1(trial,2),:)];
end

% number of trials = n
n = size(trialMat1,1);

% create trial number
trialNumber = (1:1:n)';

% create block variable
nStimuli = size(squaresPay1,1);
nBlockTrials = nchoosek(nStimuli,2)+nStimuli;
nBlocks = n/nBlockTrials;
blockNumber = repmat(1:1:nBlocks,[nBlockTrials,1]);

% part number
partNumber1 = ones(n,1);



%--------------------------------------------------------------------------
% Load Images
%--------------------------------------------------------------------------

%load sample images

S = dlLoadStimuli1;

for i = 1:size(S,1)
    tex{i,1} = Screen('MakeTexture',window,S(i,1).img);
end




% keyboard
% The avaliable keys to press
escapeKey = KbName('ESCAPE');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
secretKey = KbName('j');
secretKeyEnd = KbName('s');

%--------------------------------------------------------------------------
% Instructions
%--------------------------------------------------------------------------

    
% load images
images = dlLoadInstructions1;
nrImages = size(images,1);

 % create textutes
 for i = 1:size(images,1)
    textures{i,1} = Screen('MakeTexture',window,images(i,1).img);
end
    
    % draw the first texture
    ind = 1;
    secret = 0;
    while (secret==0)
        if(ind<=1)
            ind = 1;
            Screen('DrawTexture', window, textures{ind, 1});
            Screen('flip',window);
            ind = 2;
        end
            
        [secs, keyCode, deltaSecs] = KbWait([],2);

            if(keyCode(rightKey))
                
                keyCode(rightKey) = 0;
                
                ind = ind + 1;
                Screen('DrawTexture', window, textures{ind, 1});
                Screen('flip',window);
                if(ind==nrImages)
                    ind = ind - 1;
                end
                
            elseif(keyCode(leftKey))
                
                keyCode(leftKey) = 0;
                
                if(ind<=0)
                    Screen('DrawTexture', window, textures{1, 1});
                    Screen('flip',window);
                else
                    Screen('DrawTexture', window, textures{ind-1, 1});
                    Screen('flip',window);
                    
                end
                ind = ind - 1;
            elseif(keyCode(secretKey))
                secret = 1;
            else
                continue
            end
            
        
    end

%--------------------------------------------------------------------------
% Welcome Screen
%--------------------------------------------------------------------------


DrawFormattedText(window,['Press any key to begin the task.'],'center','center',textColor);
Screen('Flip',window);
KbWait([], 2);
KbWait([], 1);

rt1           = -1*ones(n,1);
ansBox1       = -1*ones(n,1);
ansSquare1    = -1*ones(n,1);
ansPay1       = -1*ones(n,1);
ansPayNow1    = -1*ones(n,1);
ansPayLater1  = -1*ones(n,1);
% ansColor1     = -1*ones(n,3);
f2_1 = zeros(1,n);
f1_1 = zeros(1,n);

 
  
for trial = 1:size(isi_t2_1,1)
    if(trial>2)
        
        %------------------------------------------------------------------
        % Get Response
        %------------------------------------------------------------------
        
        
        % Get the centre coordinate of the window
        [midX,midY] = RectCenter(windowRect);
        
        stTime = GetSecs;
        
        % Get a response
        while (GetSecs-stTime<=3)
            
            
            %--------------------------------------------------------------
            % Display choice screen
            %--------------------------------------------------------------
            
            
            % first 2 squares are the choice options
%            DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
            Screen('DrawTexture', window, tex{trialMat1(trial,1),:},[], rectC(1,:),0,[],1);
            Screen('DrawTexture', window, tex{trialMat1(trial,2),:},[], rectC(2,:),0,[],1);
%             Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%             Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
            
            % last 2 squares are the feedback options
            %             DrawFormattedText(window, 'chosen last round', xText(2), yText(2), black);
            %             DrawFormattedText(window, 'choosen before', xText(3), yText(3), black);
            %             Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(3,:));
            %             DrawFormattedText(window, ['+', num2str(squaresPay(ansSquare(trial-1,1),1))], 'center',yposN(1), numberColor);
            %             Screen('FillRect', window, colorSquares(ansSquare(trial-2,1),:), rectC(4,:));
            %             DrawFormattedText(window, ['+', num2str(squaresPay(ansSquare(trial-2,1),2))], 'center',yposN(2), numberColor);
            %
            
            
            %[~,x,y,whichButton] = GetClicks(window);
            Screen('Flip',window);
            
            [keyIsDown,secs, keyCode] = KbCheck;
            
            if(keyIsDown~=1)
                continue
            else

                if(keyCode(leftKey))
                    boxID = 1;
                    
%                    DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt1(trial,1)   = GetSecs-stTime;
                     pause(3-rt1);
                    
                elseif(keyCode(rightKey) && (GetSecs-stTime)<=3)
                    boxID = 2;
                    
%                    DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt1(trial,1)   = GetSecs-stTime;
                     pause(3-rt1);
                    
                else
                    continue
                end
                
                                % record answers
%                 rt(trial,1)   = GetSecs-stTime;
                ansBox1(trial,1) = boxID; % left or right box
                ansSquare1(trial,1) = trialMat1(trial,boxID); % square number
                if(boxID==1)
                    ansPay1(trial,1) = sum(payTotalLeft1(trial,:));
                    ansPayNow1(trial,1) = payTotalLeft1(trial,1);
                    ansPayLater1(trial,1) = payTotalLeft1(trial,2);
                else
                    ansPay1(trial,1) = sum(payTotalRight1(trial,:));
                    ansPayNow1(trial,1) = payTotalRight1(trial,1);
                    ansPayLater1(trial,1) = payTotalRight1(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat1(trial,boxID),:); % square color
                break;

            end
         
        end
        % if the timeout expired without an answer
        if(ansBox1(trial,1)==-1)
            % make a random choice
                    % put in nan
                boxID = randi([1,2],1);
                rt1(trial,1)   = NaN;
                ansSquare1(trial,1) = trialMat1(trial,boxID); % square number
                if(boxID==1)
                    ansPay1(trial,1) = sum(payTotalLeft1(trial,:));
                    ansPayNow1(trial,1) = payTotalLeft1(trial,1);
                    ansPayLater1(trial,1) = payTotalLeft1(trial,2);
                else
                    ansPay1(trial,1) = sum(payTotalRight1(trial,:));
                    ansPayNow1(trial,1) = payTotalRight1(trial,1);
                    ansPayLater1(trial,1) = payTotalRight1(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat1(trial,boxID),:); % square color
        end
        
        Screen('FillRect', window, fixColor, FixCross');
        Screen('Flip',window);

        %------------------------------------------------------------------
        % Display feedback screen
        %------------------------------------------------------------------
        if(~isnan(ansBox1(trial,1)))
            
        % wait time = feedback timing - stimulus timing - 3s choice
        % presentation time
        f1_1(trial) = isi_t2_1(trial,1)-isi_t1_1(trial,1)-3;
        WaitSecs(f1_1(trial));
            
        % last 2 squares are the feedback options
        % Screen('TextSize', window, 23  );
        
             if(~isnan(ansBox1(trial-1,1)))
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor); 
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare1(trial-1,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater1(trial-1,1))],'center',yposN(2), numberColor);
             elseif(isnan(ansBox1(trial-1,1)))
                trial_nonNaN = find(~isnan(ansBox1([1:(trial-1)],1)), 1, 'last');
                if(~isempty(trial_nonNaN))
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor); 
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare1(trial_nonNaN,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial_nonNaN,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater1(trial_nonNaN,1))],'center',yposN(2), numberColor);
                else
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1},[], rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
                end
             else
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor); 
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1}, [],rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
             end
        
        Screen('Flip',window);
        
%         % GetImage call. Alter the rect argument to change the location of the screen shot
%     imageArray = Screen('GetImage', window, []);
% 
%     % imwrite is a Matlab function, not a PTB-3 function
%     imwrite(imageArray, 'test.jpg')
    
        % wait 2 seconds at the feedback screen
        WaitSecs(2);

        end
        
        
    elseif(trial==1)
        %------------------------------------------------------------------
        % Get Response
        %------------------------------------------------------------------
        
        % Get the centre coordinate of the window
        [midX,midY] = RectCenter(windowRect);
        
        stTime = GetSecs;
        
        % Get a response
        while ((GetSecs-stTime)<=3)
            
            %--------------------------------------------------------------
            % Display choice screen
            %--------------------------------------------------------------
%             DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
            Screen('DrawTexture', window, tex{trialMat1(trial,1),1},[], rectC(1,:), 0, [], 1);
            Screen('DrawTexture', window, tex{trialMat1(trial,2),1},[], rectC(2,:), 0, [], 1);
%             Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%             Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
            
            
            %[~,x,y,whichButton] = GetClicks(window);
            Screen('Flip',window);
            
            [keyIsDown,secs, keyCode] = KbCheck;
            
            if(keyIsDown~=1)
                continue
            else
                if(keyCode(leftKey))
                    boxID = 1;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt1(trial,1)   = GetSecs-stTime;
                     pause(3-rt1);
                elseif(keyCode(rightKey))
                    boxID = 2;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt1(trial,1)   = GetSecs-stTime;
                     pause(3-rt1);
                else
                    continue
                end
                
                % record answers
%                 rt(trial,1)   = GetSecs-stTime;
                ansBox1(trial,1) = boxID; % left or right box
                ansSquare1(trial,1) = trialMat1(trial,boxID); % square number
                if(boxID==1)
                    ansPay1(trial,1) = sum(payTotalLeft1(trial,:));
                    ansPayNow1(trial,1) = payTotalLeft1(trial,1);
                    ansPayLater1(trial,1) = payTotalLeft1(trial,2);
                else
                    ansPay1(trial,1) = sum(payTotalRight1(trial,:));
                    ansPayNow1(trial,1) = payTotalRight1(trial,1);
                    ansPayLater1(trial,1) = payTotalRight1(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat1(trial,boxID),:); % square color
                break;
            end
            
        end
        
                % if the timeout expired without an answer
        if(ansBox1(trial,1)==-1)
            % make a random choice
                    % put in nan
                boxID = randi([1,2],1);
                rt1(trial,1)   = NaN;
                ansSquare1(trial,1) = trialMat1(trial,boxID); % square number
                if(boxID==1)
                    ansPay1(trial,1) = sum(payTotalLeft1(trial,:));
                    ansPayNow1(trial,1) = payTotalLeft1(trial,1);
                    ansPayLater1(trial,1) = payTotalLeft1(trial,2);
                else
                    ansPay1(trial,1) = sum(payTotalRight1(trial,:));
                    ansPayNow1(trial,1) = payTotalRight1(trial,1);
                    ansPayLater1(trial,1) = payTotalRight1(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat1(trial,boxID),:); % square color
        end
        Screen('FillRect', window, fixColor  , FixCross');
        Screen('Flip',window);
        
        %------------------------------------------------------------------
        % Display feedback screen
        %------------------------------------------------------------------
        if(~isnan(ansBox1(trial,1)))
        % wait time = feedback timing - stimulus timing - 3s choice
        % presentation time
        f1_1(trial) = isi_t2_1(trial,1)-isi_t1_1(trial,1)-3;
        WaitSecs(f1_1(trial));
        
        % last 1 squares are the feedback options
        % Screen('TextSize', window, 23);
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1}, [],rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
        Screen('Flip',window);
        
        WaitSecs(2);
        end
        
        
    else % if trial==2
        %------------------------------------------------------------------
        % Get Response
        %------------------------------------------------------------------
        
        % Get the centre coordinate of the window
        [midX,midY] = RectCenter(windowRect);
        
        stTime = GetSecs;
        
        % Get a response
        while ((GetSecs-stTime)<=3)
            
            
            %--------------------------------------------------------------
            % Display choice screen
            %--------------------------------------------------------------
%             DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
            Screen('DrawTexture', window, tex{trialMat1(trial,1),1},[], rectC(1,:), 0, [], 1);
            Screen('DrawTexture', window, tex{trialMat1(trial,2),1},[], rectC(2,:), 0, [], 1);
%             Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%             Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
            % and 1 feedback square from previous round
            %             Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(3,:));
            %             DrawFormattedText(window, ['+', num2str(squaresPay(ansSquare(trial-1,1),1))], 'center',yposN(1), numberColor);
            %             DrawFormattedText(window, 'chosen last round', xText(2), yText(2), black);
            %
            %[~,x,y,whichButton] = GetClicks(window);
            Screen('Flip',window);
            
            [keyIsDown,secs, keyCode] = KbCheck;
            
            if(keyIsDown~=1)
                continue
            else
                if(keyCode(leftKey))
                    boxID = 1;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt1(trial,1)   = GetSecs-stTime;
                     pause(3-rt1);
                elseif(keyCode(rightKey))
                    boxID = 2;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt1(trial,1)   = GetSecs-stTime;
                     pause(3-rt1);
                else
                    continue
                end
                
                % record answers
%                 rt(trial,1)   = GetSecs-stTime;
                ansBox1(trial,1) = boxID; % left or right box
                ansSquare1(trial,1) = trialMat1(trial,boxID); % square number
                if(boxID==1)
                    ansPay1(trial,1) = sum(payTotalLeft1(trial,:));
                    ansPayNow1(trial,1) = payTotalLeft1(trial,1);
                    ansPayLater1(trial,1) = payTotalLeft1(trial,2);
                else
                    ansPay1(trial,1) = sum(payTotalRight1(trial,:));
                    ansPayNow1(trial,1) = payTotalRight1(trial,1);
                    ansPayLater1(trial,1) = payTotalRight1(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat1(trial,boxID),:); % square color
                break;
            end
        end
        
                % if the timeout expired without an answer
        if(ansBox1(trial,1)==-1)
            % make a random choice
                    % put in nan
                boxID = randi([1,2],1);
                rt1(trial,1)   = NaN;
                ansSquare1(trial,1) = trialMat1(trial,boxID); % square number
                if(boxID==1)
                    ansPay1(trial,1) = sum(payTotalLeft1(trial,:));
                    ansPayNow1(trial,1) = payTotalLeft1(trial,1);
                    ansPayLater1(trial,1) = payTotalLeft1(trial,2);
                else
                    ansPay1(trial,1) = sum(payTotalRight1(trial,:));
                    ansPayNow1(trial,1) = payTotalRight1(trial,1);
                    ansPayLater1(trial,1) = payTotalRight1(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat1(trial,boxID),:); % square color
        end
        
        Screen('FillRect', window, fixColor  , FixCross');  
        Screen('Flip',window);
        
        %------------------------------------------------------------------
        % Display feedback screen
        %------------------------------------------------------------------
        if(~isnan(ansBox1(trial,1)))
        % wait time = feedback timing - stimulus timing - 3s choice
        % presentation time
        f1_1(trial) = isi_t2_1(trial,1)-isi_t1_1(trial,1)-3;
        WaitSecs(f1_1(trial));
        
        % last 2 squares are the feedback options
        % Screen('TextSize', window, 23  );

             if(~isnan(ansBox1(trial-1,1)))
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor); 
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare1(trial-1,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater1(trial-1,1))],'center',yposN(2), numberColor);
             elseif(isnan(ansBox1(trial-1,1)) && (trial-1)~=1)
                trial_nonNaN = find(~isnan(ansBox1([1:(trial-1)],1)), 1, 'last');
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor);  
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare1(trial_nonNaN,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial_nonNaN,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater1(trial_nonNaN,1))],'center',yposN(2), numberColor);
             else
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);  
        Screen('DrawTexture', window, tex{ansSquare1(trial,1),1}, [],rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
             end

        
        Screen('Flip',window);
        WaitSecs(2);
        end
    end
    %----------------------------------------------------------------------
    % Display fixation cross (iti)
    %----------------------------------------------------------------------
    
    Screen('FillRect', window, fixColor  , FixCross');
    Screen('Flip',window);

    % f2 = time of second stimulus - feedback time - 2s presentation of
    % feedback

    f2_1(trial) = isi_t1_1(trial+1,1)-isi_t2_1(trial,1)-2;
    WaitSecs(f2_1(trial));

end



%--------------------------------------------------------------------------
% Calculate points earned
%--------------------------------------------------------------------------


totalPay1 = sum(ansPay1(~isnan(ansPay1(:,1)),1)) - 5*sum(isnan(rt1));

totalMoney1 = 0.05*(totalPay1 - 1680);

%--------------------------------------------------------------------------
% Answer Matrix
%--------------------------------------------------------------------------


answerMat1 = [rt1, ansBox1, ansSquare1, ansPay1, ansPayNow1, ansPayLater1,...
    payTrueLeft1, payTrueRight1, payNoiseLeft1, payNoiseRight1, payTotalLeft1, payTotalRight1,...
    payTrialMat1, trialNumber, blockNumber(:), partNumber1];

isi_t1_1(106,:) = [];

timingMat1 = [isi_t1_1, isi_t2_1, f1_1', f2_1'];

%     'Your total points for this part are: ', num2str(totalPay1),'\n',...
%     'Your total earnings for this part are: ', num2str(totalMoney1),'\n',...


DrawFormattedText(window,['You can take a short break before part 2 starts','\n',...
    'When you are ready press any key to continue.'],'center','center',textColor);
Screen('Flip',window);
KbWait([], 2);
KbWait([], 1);


%% part 2



%--------------------------------------------------------------------------
% Trial matrix
%--------------------------------------------------------------------------

% associate each square number with a color
% cyan        = [0.2 0.8 0.8];
% brown       = [0.2 0 0];
% orange      = [1 0.5 0];
% blue        = [0 0.5 1];
% green       = [0 0.6 0.3];
% red         = [1 0.2 0.2];
% purple      = [0.5 0 0.5];

% colorSquares = [0.5 0 0.5;
%     0 0.5 1;
%     0.2 0.8 0.8;
%     0 0.6 0.3;
%     1 0.5 0;
%     1 0.2 0.2];
% 
% colorNames = ['purple';
%     'blue';
%     'cyan';
%     'green';
%     'orange';
%     'red'];
% 
% % permute the rows for every subject
% [m, n] = size(colorSquares);
% rColor = randperm(m);
% colorSquares = colorSquares(rColor,:);
% colorNames = colorNames(rColor,:);


% magenta = [1 0 1];
% yellow = [1 0.8 0];
% turquoise = [0 0.4 0.4];
% brown = [0.2 0 0];
% light green = [0 1 0.2];
% pink = [0.8 0 0.4];

% colorSquares = [1 0 1;
%     1 0.8 0;
%     0 0.4 0.4;  
%     0.2 0 0;
%     0 1 0.2;
%     0.8 0 0.4];

% define payoffs for each square color and random order for the 2 sets


if(sCondition == 1)
    if(ord == 0)
        
        squaresPay2 = [7,11;
            10,6;
            5,9;
            8,4;
            3,7;
            6,2];
    else
        squaresPay2 = [11,7;
            6,10;
            9,5;
            4,8;
            7,3;
            2,6];
    end
    
else
    disp('NOT A VALID CONDITION');
%     if(ord<=0.5)
%         squaresPay2 = [6,10;
%             10,6;
%             4,8;
%             8,4;
%             2,6;
%             6,2];
%     else
%         squaresPay2 = [11,7;
%             7,11;
%             9,5;
%             5,9;
%             7,3;
%             3,7];
%     end
end


% permute payoffs for each subjects
[p2, q2] = size(squaresPay2);
squaresPay2 = squaresPay2(randperm(p2),:);

% trial matrix with the choice options
trialMat2 = createTrialMatrix();

% add noise payoffs
payTrueLeft2 = squaresPay2(trialMat2(:,1),:);
payTrueRight2 = squaresPay2(trialMat2(:,2),:);

payNoiseLeft2 = randi([0,3],size(trialMat2,1),2);
payNoiseRight2 = randi([0,3],size(trialMat2,1),2);

payTotalLeft2 = payTrueLeft2 + payNoiseLeft2;
payTotalRight2 = payTrueRight2 + payNoiseRight2;

payTotal2 = [payTotalLeft2, payTotalRight2];


% create trial matrix for the colors, one number/string for each color
% and for the payoffs
for trial = 1:size(trialMat2,1)
%     colorTrialMat(trial,:) = [colorSquares(trialMat(trial,1),:), colorSquares(trialMat(trial,2),:)];
%     colorNamesTrialMat(trial,:) =  [colorNames(trialMat(trial,1),:), colorNames(trialMat(trial,2),:)];
    payTrialMat2(trial,:) = [squaresPay2(trialMat2(trial,1),:), squaresPay2(trialMat2(trial,2),:)];
end

% number of trials = n
n = size(trialMat2,1);

% create trial number
trialNumber = (1:1:n)';

% create block variable
nStimuli = size(squaresPay2,1);
nBlockTrials = nchoosek(nStimuli,2)+nStimuli;
nBlocks = n/nBlockTrials;
blockNumber = repmat(1:1:nBlocks,[nBlockTrials,1]);

% part number
partNumber2 = 2*ones(n,1);



%--------------------------------------------------------------------------
% Load Images
%--------------------------------------------------------------------------

%load sample images

S = dlLoadStimuli2;

for i = 1:size(S,1)
    tex{i,1} = Screen('MakeTexture',window,S(i,1).img);
end




% keyboard
% The avaliable keys to press
escapeKey = KbName('ESCAPE');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');

%--------------------------------------------------------------------------
% Instructions
%--------------------------------------------------------------------------

    
% load images
images = dlLoadInstructions2;
nrImages = size(images,1);

 % create textutes
 for i = 1:size(images,1)
    textures{i,1} = Screen('MakeTexture',window,images(i,1).img);
end
    
% draw the first texture
ind = 1;
%     Screen('DrawTexture', windowHandle, textures{ind, 1});
%     Screen('flip', windowHandle);
while (ind<=nrImages)
    if(ind<=1)
        ind = 1;
        Screen('DrawTexture', window, textures{ind, 1});
        Screen('flip',window);
        ind = 2;
    end
    
    [secs, keyCode, deltaSecs] = KbWait([],2);
    
    if(keyCode(rightKey))
        
        keyCode(rightKey) = 0;
        
        ind = ind + 1;
        if(ind==nrImages)
        	Screen('DrawTexture', window, textures{ind, 1});
            Screen('flip',window);
            WaitSecs(0.5);
            FlushEvents('keyDown');
            % Wait for button press
            while KbCheck
            end
        else
        
            Screen('DrawTexture', window, textures{ind, 1});
            Screen('flip',window);
        end

    elseif(keyCode(leftKey))
        
        keyCode(leftKey) = 0;
        
        if(ind<=0)
            Screen('DrawTexture', window, textures{1, 1});
            Screen('flip',window);
        else
            Screen('DrawTexture', window, textures{ind-1, 1});
            Screen('flip',window);
            
        end
        ind = ind - 1;
    else
        continue
    end
    
end


%--------------------------------------------------------------------------
% Welcome Screen
%--------------------------------------------------------------------------


DrawFormattedText(window,['Press any key to begin the task.'],'center','center',textColor);
Screen('Flip',window);
KbWait([], 2);
KbWait([], 1);

rt2           = -1*ones(n,1);
ansBox2       = -1*ones(n,1);
ansSquare2    = -1*ones(n,1);
ansPay2       = -1*ones(n,1);
ansPayNow2    = -1*ones(n,1);
ansPayLater2  = -1*ones(n,1);
% ansColor1     = -1*ones(n,3);
f2_2 = zeros(1,n);
f1_2 = zeros(1,n);

 
  
for trial = 1:size(isi_t2_2,1)
    if(trial>2)
        
        %------------------------------------------------------------------
        % Get Response
        %------------------------------------------------------------------
        
        
        % Get the centre coordinate of the window
        [midX,midY] = RectCenter(windowRect);
        
        stTime = GetSecs;
        
        % Get a response
        while (GetSecs-stTime<=3)
            
            
            %--------------------------------------------------------------
            % Display choice screen
            %--------------------------------------------------------------
            
            
            % first 2 squares are the choice options
%            DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
            Screen('DrawTexture', window, tex{trialMat2(trial,1),:},[], rectC(1,:),0,[],1);
            Screen('DrawTexture', window, tex{trialMat2(trial,2),:},[], rectC(2,:),0,[],1);
%             Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%             Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
            
            % last 2 squares are the feedback options
            %             DrawFormattedText(window, 'chosen last round', xText(2), yText(2), black);
            %             DrawFormattedText(window, 'choosen before', xText(3), yText(3), black);
            %             Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(3,:));
            %             DrawFormattedText(window, ['+', num2str(squaresPay(ansSquare(trial-1,1),1))], 'center',yposN(1), numberColor);
            %             Screen('FillRect', window, colorSquares(ansSquare(trial-2,1),:), rectC(4,:));
            %             DrawFormattedText(window, ['+', num2str(squaresPay(ansSquare(trial-2,1),2))], 'center',yposN(2), numberColor);
            %
            
            
            %[~,x,y,whichButton] = GetClicks(window);
            Screen('Flip',window);
            
            [keyIsDown,secs, keyCode] = KbCheck;
            
            if(keyIsDown~=1)
                continue
            else

                if(keyCode(leftKey))
                    boxID = 1;
                    
%                    DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt2(trial,1)   = GetSecs-stTime;
                     pause(3-rt2);
                    
                elseif(keyCode(rightKey) && (GetSecs-stTime)<=3)
                    boxID = 2;
                    
%                    DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt2(trial,1)   = GetSecs-stTime;
                     pause(3-rt2);
                    
                else
                    continue
                end
                
                                % record answers
%                 rt(trial,1)   = GetSecs-stTime;
                ansBox2(trial,1) = boxID; % left or right box
                ansSquare2(trial,1) = trialMat2(trial,boxID); % square number
                if(boxID==1)
                    ansPay2(trial,1) = sum(payTotalLeft2(trial,:));
                    ansPayNow2(trial,1) = payTotalLeft2(trial,1);
                    ansPayLater2(trial,1) = payTotalLeft2(trial,2);
                else
                    ansPay2(trial,1) = sum(payTotalRight2(trial,:));
                    ansPayNow2(trial,1) = payTotalRight2(trial,1);
                    ansPayLater2(trial,1) = payTotalRight2(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat2(trial,boxID),:); % square color
                break;

            end
         
        end
        % if the timeout expired without an answer
        if(ansBox2(trial,1)==-1)
            % make a random choice
                    % put in nan
                boxID = randi([1,2],1);
                rt2(trial,1)   = NaN;
                ansSquare2(trial,1) = trialMat2(trial,boxID); % square number
                if(boxID==1)
                    ansPay2(trial,1) = sum(payTotalLeft2(trial,:));
                    ansPayNow2(trial,1) = payTotalLeft2(trial,1);
                    ansPayLater2(trial,1) = payTotalLeft2(trial,2);
                else
                    ansPay2(trial,1) = sum(payTotalRight2(trial,:));
                    ansPayNow2(trial,1) = payTotalRight2(trial,1);
                    ansPayLater2(trial,1) = payTotalRight2(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat2(trial,boxID),:); % square color
        end
        
        Screen('FillRect', window, fixColor, FixCross');
        Screen('Flip',window);

        %------------------------------------------------------------------
        % Display feedback screen
        %------------------------------------------------------------------
        if(~isnan(ansBox2(trial,1)))
            
        % wait time = feedback timing - stimulus timing - 3s choice
        % presentation time
        f1_2(trial) = isi_t2_2(trial,1)-isi_t1_2(trial,1)-3;
        WaitSecs(f1_2(trial));
            
        % last 2 squares are the feedback options
        % Screen('TextSize', window, 23  );
        
             if(~isnan(ansBox2(trial-1,1)))
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor); 
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare2(trial-1,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater2(trial-1,1))],'center',yposN(2), numberColor);
             elseif(isnan(ansBox2(trial-1,1)))
                trial_nonNaN = find(~isnan(ansBox2([1:(trial-1)],1)), 1, 'last');
                if(~isempty(trial_nonNaN))
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor); 
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare2(trial_nonNaN,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial_nonNaN,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater2(trial_nonNaN,1))],'center',yposN(2), numberColor);
                else
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1},[], rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
                end
             else
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor); 
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1}, [],rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
             end
        
        Screen('Flip',window);
        
%         % GetImage call. Alter the rect argument to change the location of the screen shot
%     imageArray = Screen('GetImage', window, []);
% 
%     % imwrite is a Matlab function, not a PTB-3 function
%     imwrite(imageArray, 'test.jpg')
    
        % wait 2 seconds at the feedback screen
        WaitSecs(2);

        end
        
        
    elseif(trial==1)
        %------------------------------------------------------------------
        % Get Response
        %------------------------------------------------------------------
        
        % Get the centre coordinate of the window
        [midX,midY] = RectCenter(windowRect);
        
        stTime = GetSecs;
        
        % Get a response
        while ((GetSecs-stTime)<=3)
            
            %--------------------------------------------------------------
            % Display choice screen
            %--------------------------------------------------------------
%             DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
            Screen('DrawTexture', window, tex{trialMat2(trial,1),1},[], rectC(1,:), 0, [], 1);
            Screen('DrawTexture', window, tex{trialMat2(trial,2),1},[], rectC(2,:), 0, [], 1);
%             Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%             Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
            
            
            %[~,x,y,whichButton] = GetClicks(window);
            Screen('Flip',window);
            
            [keyIsDown,secs, keyCode] = KbCheck;
            
            if(keyIsDown~=1)
                continue
            else
                if(keyCode(leftKey))
                    boxID = 1;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt2(trial,1)   = GetSecs-stTime;
                     pause(3-rt2);
                elseif(keyCode(rightKey))
                    boxID = 2;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt2(trial,1)   = GetSecs-stTime;
                     pause(3-rt2);
                else
                    continue
                end
                
                % record answers
%                 rt(trial,1)   = GetSecs-stTime;
                ansBox2(trial,1) = boxID; % left or right box
                ansSquare2(trial,1) = trialMat2(trial,boxID); % square number
                if(boxID==1)
                    ansPay2(trial,1) = sum(payTotalLeft2(trial,:));
                    ansPayNow2(trial,1) = payTotalLeft2(trial,1);
                    ansPayLater2(trial,1) = payTotalLeft2(trial,2);
                else
                    ansPay2(trial,1) = sum(payTotalRight2(trial,:));
                    ansPayNow2(trial,1) = payTotalRight2(trial,1);
                    ansPayLater2(trial,1) = payTotalRight2(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat2(trial,boxID),:); % square color
                break;
            end
            
        end
        
                % if the timeout expired without an answer
        if(ansBox2(trial,1)==-1)
            % make a random choice
                    % put in nan
                boxID = randi([1,2],1);
                rt2(trial,1)   = NaN;
                ansSquare2(trial,1) = trialMat2(trial,boxID); % square number
                if(boxID==1)
                    ansPay2(trial,1) = sum(payTotalLeft2(trial,:));
                    ansPayNow2(trial,1) = payTotalLeft2(trial,1);
                    ansPayLater2(trial,1) = payTotalLeft2(trial,2);
                else
                    ansPay2(trial,1) = sum(payTotalRight2(trial,:));
                    ansPayNow2(trial,1) = payTotalRight2(trial,1);
                    ansPayLater2(trial,1) = payTotalRight2(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat2(trial,boxID),:); % square color
        end
        Screen('FillRect', window, fixColor  , FixCross');
        Screen('Flip',window);
        
        %------------------------------------------------------------------
        % Display feedback screen
        %------------------------------------------------------------------
        if(~isnan(ansBox2(trial,1)))
        % wait time = feedback timing - stimulus timing - 3s choice
        % presentation time
        f1_2(trial) = isi_t2_2(trial,1)-isi_t1_2(trial,1)-3;
        WaitSecs(f1_2(trial));
        
        % last 1 squares are the feedback options
        % Screen('TextSize', window, 23);
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1}, [],rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
        Screen('Flip',window);
        
        WaitSecs(2);
        end
        
        
    else % if trial==2
        %------------------------------------------------------------------
        % Get Response
        %------------------------------------------------------------------
        
        % Get the centre coordinate of the window
        [midX,midY] = RectCenter(windowRect);
        
        stTime = GetSecs;
        
        % Get a response
        while ((GetSecs-stTime)<=3)
            
            
            %--------------------------------------------------------------
            % Display choice screen
            %--------------------------------------------------------------
%             DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
            Screen('DrawTexture', window, tex{trialMat2(trial,1),1},[], rectC(1,:), 0, [], 1);
            Screen('DrawTexture', window, tex{trialMat2(trial,2),1},[], rectC(2,:), 0, [], 1);
%             Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%             Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
            % and 1 feedback square from previous round
            %             Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(3,:));
            %             DrawFormattedText(window, ['+', num2str(squaresPay(ansSquare(trial-1,1),1))], 'center',yposN(1), numberColor);
            %             DrawFormattedText(window, 'chosen last round', xText(2), yText(2), black);
            %
            %[~,x,y,whichButton] = GetClicks(window);
            Screen('Flip',window);
            
            [keyIsDown,secs, keyCode] = KbCheck;
            
            if(keyIsDown~=1)
                continue
            else
                if(keyCode(leftKey))
                    boxID = 1;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt2(trial,1)   = GetSecs-stTime;
                     pause(3-rt2);
                elseif(keyCode(rightKey))
                    boxID = 2;
%                     DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
%                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
%                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
%                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
%                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
%                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
%                     Screen('Flip',window);
                     rt2(trial,1)   = GetSecs-stTime;
                     pause(3-rt2);
                else
                    continue
                end
                
                % record answers
%                 rt(trial,1)   = GetSecs-stTime;
                ansBox2(trial,1) = boxID; % left or right box
                ansSquare2(trial,1) = trialMat2(trial,boxID); % square number
                if(boxID==1)
                    ansPay2(trial,1) = sum(payTotalLeft2(trial,:));
                    ansPayNow2(trial,1) = payTotalLeft2(trial,1);
                    ansPayLater2(trial,1) = payTotalLeft2(trial,2);
                else
                    ansPay2(trial,1) = sum(payTotalRight2(trial,:));
                    ansPayNow2(trial,1) = payTotalRight2(trial,1);
                    ansPayLater2(trial,1) = payTotalRight2(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat2(trial,boxID),:); % square color
                break;
            end
        end
        
                % if the timeout expired without an answer
        if(ansBox2(trial,1)==-1)
            % make a random choice
                    % put in nan
                boxID = randi([1,2],1);
                rt2(trial,1)   = NaN;
                ansSquare2(trial,1) = trialMat2(trial,boxID); % square number
                if(boxID==1)
                    ansPay2(trial,1) = sum(payTotalLeft2(trial,:));
                    ansPayNow2(trial,1) = payTotalLeft2(trial,1);
                    ansPayLater2(trial,1) = payTotalLeft2(trial,2);
                else
                    ansPay2(trial,1) = sum(payTotalRight2(trial,:));
                    ansPayNow2(trial,1) = payTotalRight2(trial,1);
                    ansPayLater2(trial,1) = payTotalRight2(trial,2);
                end
%                 ansColor1(trial,:) = colorSquares(trialMat2(trial,boxID),:); % square color
        end
        
        Screen('FillRect', window, fixColor  , FixCross');  
        Screen('Flip',window);
        
        %------------------------------------------------------------------
        % Display feedback screen
        %------------------------------------------------------------------
        if(~isnan(ansBox2(trial,1)))
        % wait time = feedback timing - stimulus timing - 3s choice
        % presentation time
        f1_2(trial) = isi_t2_2(trial,1)-isi_t1_2(trial,1)-3;
        WaitSecs(f1_2(trial));
        
        % last 2 squares are the feedback options
        % Screen('TextSize', window, 23  );

             if(~isnan(ansBox2(trial-1,1)))
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor); 
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare2(trial-1,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial-1,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater2(trial-1,1))],'center',yposN(2), numberColor);
             elseif(isnan(ansBox2(trial-1,1)) && (trial-1)~=1)
                trial_nonNaN = find(~isnan(ansBox2([1:(trial-1)],1)), 1, 'last');
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
%         DrawFormattedText(window, 'choosen before', xText(3), yText(3), textColor);  
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1},[], rectC(3,:), 0, [], 1);
        Screen('DrawTexture', window, tex{ansSquare2(trial_nonNaN,1),1},[], rectC(4,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
%         Screen('FillRect', window, colorSquares(ansSquare(trial_nonNaN,1),:), rectC(4,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
        DrawFormattedText(window, ['+', num2str(ansPayLater2(trial_nonNaN,1))],'center',yposN(2), numberColor);
             else
%         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);  
        Screen('DrawTexture', window, tex{ansSquare2(trial,1),1}, [],rectC(3,:), 0, [], 1);
%         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
        DrawFormattedText(window, ['+', num2str(ansPayNow2(trial,1))], 'center',yposN(1), numberColor);
             end

        
        Screen('Flip',window);
        WaitSecs(2);
        end
    end
    %----------------------------------------------------------------------
    % Display fixation cross (iti)
    %----------------------------------------------------------------------
    
    Screen('FillRect', window, fixColor  , FixCross');
    Screen('Flip',window);

    % f2 = time of second stimulus - feedback time - 2s presentation of
    % feedback

    f2_2(trial) = isi_t1_2(trial+1,1)-isi_t2_2(trial,1)-2;
    WaitSecs(f2_2(trial));

    
end



%--------------------------------------------------------------------------
% Calculate points earned
%--------------------------------------------------------------------------


totalPay2 = sum(ansPay2(~isnan(ansPay2(:,1)),1)) - 5*sum(isnan(rt2));

totalMoney2 = 0.05*(totalPay2 - 1680);

% total earnings



%--------------------------------------------------------------------------
% Answer Matrix
%--------------------------------------------------------------------------


answerMat2 = [rt2, ansBox2, ansSquare2, ansPay2, ansPayNow2, ansPayLater2,...
    payTrueLeft2, payTrueRight2, payNoiseLeft2, payNoiseRight2, payTotalLeft2, payTotalRight2,...
    payTrialMat2, trialNumber, blockNumber(:), partNumber2];

isi_t1_2(106,:) = [];

timingMat2 = [isi_t1_2, isi_t2_2, f1_2', f2_2'];

%     'Your total points for this part are: ', num2str(totalPay2),'\n',...
%     'Your total earnings ($) for this part are: ', num2str(totalMoney2),'\n',...


if(totalMoney1 + totalMoney2 < 0)
    totalMoney = 5;
else
    totalMoney = ceil(totalMoney1+totalMoney2) + 5;
end

    
    
secretEnd = 0;
while(secretEnd==0)
    
    DrawFormattedText(window,['Your total points for the experiment: ', num2str(totalPay1+totalPay2),'\n',...
        'Your total earnings ($) ','\n',...
        'for the experiment (including the show-up bonus): ', num2str(totalMoney)],'center','center',textColor);
    Screen('Flip',window);
    
    [secs, keyCode, deltaSecs] = KbWait([],2);
    
    if(~keyCode(secretKeyEnd))
        secretEnd = 0;
        DrawFormattedText(window,['Your total points for the experiment: ', num2str(totalPay1+totalPay2),'\n',...
            'Your total earnings ($) ','\n',...
            'for the experiment (including the show-up bonus): ', num2str(totalMoney)],'center','center',textColor);
        Screen('Flip',window);
    else
        secretEnd = 1;
    end
    
end





%% save data

answerMat = [answerMat1; answerMat2];
timingMat = [timingMat1; timingMat2];
trialMat  = [trialMat1; trialMat2];

if ~exist('completeData', 'dir')
       mkdir('completeData')
end
  
sVec = ones(size(answerMat,1),1).*sID;
cVec = ones(size(answerMat,1),1).*sCondition;
oVec = ones(size(answerMat,1),1).*ord;
dataFile = fullfile('completeData',[datestr(now,'mm-dd-yyyy'),'-S',num2str(sID),'-C',num2str(sCondition),'.mat']);
save(dataFile,'answerMat','trialMat','timingMat','sVec','cVec','oVec');
fclose('all');

  
Screen('CloseAll');

% totalPay = ceil(totalPay1+totalPay2);
% display(totalPay)

% totalMoney = ceil(totalMoney1+totalMoney2) + 5;
display(totalMoney);

