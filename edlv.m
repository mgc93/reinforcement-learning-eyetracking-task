%% initialize psychtoolbox

% ctrl +c to stop the experiment
% sca for closing the window


%--------------------------------------------------------------------------
% Initial Setup  
%--------------------------------------------------------------------------

try
% Clear the workspace and the screen
close all;
clear;

    %%%%%
    % if you want to do eye tracking, set this to one, otherwise, set it to
    % 0
    trackEye = 1;
    NTrial = 12; % number of trials - to change later
    %%Eye-Tracking Options
    fixDur = 1;
    fixSize = 5; %radius of fixation dot
    %%%%%


%%%%% check for the operating system
KbName('UnifyKeyNames');
if IsOSX==1
    TopPriority=0;
    sysTxtOff=0;  % no text adjustment needed
    ptbPipeMode=kPsychNeedFastBackingStore;  % enable imaging pipeline for osx
elseif IsWin==1
    TopPriority=1;
    sysTxtOff=1;  % windows draws text from the upper left corner instead of lower left.  to correct,
    % an adjustment factor of 1*letterheight is subtracted
    % from the y coordinate of the starting point
    ptbPipeMode=[];  % don't need to enable imaging pipeline
else
    ListenChar(0); clear screen; error('Operating system not OS X or Windows.')
end

    if IsWin==1
        skipSync = 0;
    else
        skipSync = 0; %if necessary to skip sync
    end %if ispc
%%%%%

% keyboard
    if IsOSX==1
        enterButton=[KbName('enter') KbName('return')];
    else
        enterButton=KbName('return');
    end %if ismac
    
% The avaliable keys to press
escapeKey = KbName('ESCAPE');
expKey = KbName('e');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
secretKey = KbName('j');
secretKeyEnd = KbName('s');


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
  

%--------------------------------------------------------------------------
% Set up Eyetracker
%--------------------------------------------------------------------------
  
    % make file for eyetracking data
    %%% eyetracking  data %%%
 
    rootDir = cd;
    expName='dlLearning';
    baseName=[num2str(sID),'_',expName];
    if trackEye == 1
        fileName = [baseName '_eyeTrackingData' '.csv'];
        eyeTrackDataFile = fopen(fileName, 'a');
        eyeLabel = '%s, %s, %s,%s,%s,%s,%s\n';
        fprintf(eyeTrackDataFile, eyeLabel, ...
            'Type','Part','sID', 'trial', 'look', 'fixStart', 'fixEnd');
    end
%%%%%  


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
greycol = 150;

% screen and text color
backgroundColor = gray1;
textColor = white;
fixColor = white;
numberColor = white;    
frameColor = [1 0 0];
mainSize = 40;
mainFont = 'Arial';

% screen visuals
fixThresh = 150; % can look 50 pixels away from fixation before it gets mad

% Open an on screen window
%[window, windowRect] = PsychImaging('OpenWindow', screenNumber, backgroundColor, [0 0 800 600]);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, backgroundColor, []);




Screen('TextFont', window,'Arial');
Screen('TextSize', window,40);
Screen('TextStyle', window,1);


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

ScreenX = globalRect(3);
ScreenY = globalRect(4);

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



%%%%% configure eyetracker

  %%% Transition Screen %%%
    Screen('FillRect', window, 0);
    if trackEye == 1
        line1 = '\n We will now calibrate the eye-tracker. Please alert the experimenter.';
    else
        line1 = '\n We will now move onto the next part of the study. Please alert the experimenter.';
    end
    
    
    %Draw all the text in one go
    Screen('TextSize', window, mainSize);
    Screen('TextFont', window, mainFont);
    DrawFormattedText(window, line1,'center', ScreenY*.35, white, [], [], [], 1.5);
    
    Screen('Flip', window);
    
    %Wait for button press
    while KbCheck
    end
    
    FlushEvents('keyDown');
    proceed = 0;
    while proceed == 0
        [keyIsDown, secs, keyCode] = KbCheck(-1);
        if keyIsDown
            if keyCode(escapeKey)
                ListenChar(0); clear screen; error('User terminated script with ESCAPE key.')
            elseif keyCode(expKey)
                proceed = 1;
            end
        end
    end %while proceed
    
    WaitSecs(.1);

if trackEye == 1
    
    % initiates the defaults for the eye tracker on the screen you
    % opened
    el = EyelinkInitDefaults(window);
    
    edfNTrials = 105;
    %create edf file
    edfFileName = ['ND_' num2str(sID)]; % CHANGE TO INCLUDE SUBJ. NUMBER, BUT LESS THAN 8 CHARACTERS. %['SelfControl' num2str(sID)];
    
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit() % Initializes Eyelink and Ethernet system. Returns: 0 if OK, -1 if error
        error('could not init connection to Eyelink')
    end
    
    % check the software version
    [~ , vs] = Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs);
    
    % open file to record data to (which is basically done
    % automatically)
    status = Eyelink('openfile', edfFileName);
    % if something goes wrong with creating the EDF, shut it down.
    if status~=0
        fprintf('Cannot create EDF file ''%s'' ', edfFileName);
        Eyelink('Shutdown');
        Screen('CloseAll');
        return;
    end
    
    
    % SET UP TRACKER CONFIGURATION
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    Eyelink('command', 'screen_pixel_coords = %ld %ld %ld %ld', 0, 0, ScreenX, ScreenY);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, ScreenX, ScreenY);
    % set calibration type. (This can also be done in the eye tracking
    % gui at the start of the study, but this makes sure that it's the
    % same for every participant)
    Eyelink('command', 'calibration_type = HV9');
    
    % set EDF file contents using the file_sample_data and
    % file-event_filter commands
    % set link data through link_sample_data and link_event_filter
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    
    % add "HTARGET" to record possible target data for EyeLink Remote
    % This is done if the version (vs) is >= 4. We don't actually have
    % to worry about this, because we use the headrest (this is for if
    % you use the target sticker and let the head move freely)
    if sscanf(vs(12:end),'%f') >=4
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERES,STATUS,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
    else
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    
    
    % make sure we're still connected.
    if Eyelink('IsConnected') == 0
        fprintf('not connected, clean up\n');
        Eyelink('ShutDown');
        Screen('CloseAll');
        return;
    end
    
    % Calibrate the eye tracker
    % setup the proper calibration foreground and background colors
    % Make sure your background color is different from the foreground
    % and msgfont colors, or else it will look like nothing is being
    % displayed
    el.backgroundcolour = 0;
    el.foregroundcolour = 255;
    el.msgfontcolour  = 255;
    el.calibrationtargetcolour = greycol;
    
    % parameters are in frequency, volume, and duration
    % set the second value in each line to 0 to turn off the sound
    el.cal_target_beep=[600, 0, 0.05];
    el.drift_correction_target_beep=[600, 0, 0.05];
    el.calibration_failed_beep=[400, 0, 0.25];
    el.calibration_success_beep=[800, 0, 0.25];
    el.drift_correction_failed_beep=[400, 0, 0.25];
    el.drift_correction_success_beep=[800, 0, 0.25];
    % you must call this function to apply the changes from above
    EyelinkUpdateDefaults(el);
    
    % Hide the mouse cursor
    HideCursor;
    EyelinkDoTrackerSetup(el); % Calibration
    %EyelinkDoDriftCorrection(el); - this is unneccessary to do
    %immediately after calibrating
    
    % Maximize Keyboard priority??
    %Priority(1);
   
end % of trackEye
    

%%%%%



%--------------------------------------------------------------------------
% Welcome Screen
%--------------------------------------------------------------------------

HideCursor();

DrawFormattedText(window,['Welcome to the experiment! Press any key to start.'],'center','center',textColor);
Screen('Flip',window);
KbWait([], 2);
KbWait([], 1);


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

% %%%%% 
%     % Check for drift (displays a dot at center, and corrects for any
%     % drift. I recommend putting this at the beginning of every section, or
%     % beginning of a block if you're doing blocks of trials.
%     EyelinkDoDriftCorrection(el);
% %%%%%

% NTrial = 12; % size(isi_t2_1,1);
eyeDataFeedback = [];
eyeDataChoice = [];
for trial = 1:7 %trialN

        
        %------------------------------------------------------------------
        % Get Response
        %------------------------------------------------------------------
        if trackEye == 1
            
            Eyelink('command', 'record_status_message "CHOICE TRIAL %d/%d"', trial, NTrial);
            
            % start recording eye position (preceded by a short pause so that
            % the tracker can finish the mode transition)
            % The paramerters for the 'StartRecording' call controls the
            % file_samples, file_events, link_samples, link_events availability
            Eyelink('Command', 'set_idle_mode');
            WaitSecs(0.05);
            Eyelink('StartRecording');
            % record a few samples before we actually start displaying
            % otherwise you may lose a few msec of data
            WaitSecs(0.1);
            eye_used = Eyelink('eyeavailable'); % get the eye that's tracked
            
            
            % used for syncing time, This will print "synctime_part1" into
            % the EDF at this timepoint. This is useful for keeping track
            % in the edf of when a trial starts, and what part of the study
            % it belongs to if you have several parts.
            Eyelink('Message', 'SYNCTIME_part1');
            % This is just telling us that the trial is valid (useful for
            % if you recallibrate during a trial, in which case this will
            % later be switched to 0 so you know to throw the trial out)
            Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 1);
            
            eye_used = Eyelink('eyeavailable'); % get eye that's tracked
            
            % Set up the two ROIs. The !V is asking for gaze/velocity
            % information, IAREA RECTANGLE is the type of roi (it's a
            % rectangle... but you can also change it to be IAREA ELLIPSE
            % for circles or IAREA FREEHAND for an irregular shape.)
            % The following "%d"s correspond to the following (for both
            % rectangle and ellipse):
            % 1: id #
            % 2 & 3: top left (x,y)
            % 4 & 5: bottom right (x,y)
            % For freehand, the first one is still id #, the following are
            % x,y coordinates for each outer portion of the ROI (only x,y pairs
            % have a comma between them).
            % The final %s corresponds to the label string you give your
            % ROI. The examples below are two ROIs for the left and right
            % side of the screen.
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, rectC(1,1),rectC(1,2),rectC(1,3),rectC(1,4),'left');
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, rectC(2,1),rectC(2,2),rectC(2,3),rectC(2,4),'right');

            
            stTime = GetSecs;
            
            % Get a response
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
            %             Screen('Flip',window);
            
            
            % Define when the trial starts here, would look something like:
            [~, trialOnset, ~, ~, ~] = Screen('Flip', window, 1);
            
            lookLeft = -1; % set default
            prevlookLeft = -1; % set default
            fixEnd = 0; % set default
            fixStart = 0;
            
            while (GetSecs-trialOnset<=3)
                
                
                % Check recording status
                err = Eyelink('CheckRecording');
                
                if err ~= 0
                    error('checkrecording problem, status: %d', err);
                end
                
                % Check for new sample
                status = Eyelink('NewFloatSampleAvailable');
                if status >= 0
                    evt = Eyelink('NewestFloatSample');
                    % if we do, get current gaze position from sample
                    eyeX = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                    eyeY = evt.gy(eye_used+1);
                    eyetime = GetSecs();
                    
                    % Make sure that the data you're recording is even
                    % valid
                    if eyeX ~= el.MISSING_DATA && eyeY ~= el.MISSING_DATA && evt.pa(eye_used+1)>0
                        
                        % See which ROI they're looking at - if they're
                        % outside of both roi's, it's counted as a -1
                        if eyeY >= rectC(1,2) && eyeY <= rectC(1,4)
                            if eyeX >= rectC(1,1) && eyeX <= rectC(1,3)
                                lookLeft = 1;
                            end
                        elseif eyeY >= rectC(2,2) && eyeY <= rectC(2,4)
                            if eyeX >= rectC(2,1) && eyeX <= rectC(2,3)
                                lookLeft = 0;
                            end
                        end
                        
                        % Making this independent of whether
                        % the eye was within an roi helps eliminate
                        % blinks
                        if lookLeft >= 0
                            if lookLeft ~= prevlookLeft
                                newFixStart = eyetime - trialOnset; % finds when the new dwell started
                                if prevlookLeft == -1
                                    fixStart = eyetime - trialOnset; %like a default
                                end
                            else
                                fixStart = newFixStart; % keeps the same fix start throughout the dwell
                                fixEnd = eyetime - trialOnset; %updates every time, until they change rois
                            end
                        end
                        
                        % If they have completed the dwell, aka changed
                        % rois, then write down the previous fixation
                        % information into your own file. This is so that
                        % you have all your information in two places (the
                        % EDF and your matlab output). Hopefully you'll
                        % only need one.
                        if lookLeft ~= prevlookLeft && fixEnd ~= 0
                            
                            % Write out the eye data to a file as it's
                            % happening, but only after a dwell has been
                            % completed  - also it says
                            % "prevlookTop" here because we want
                            % what they were just looking at, not
                            % the roi they switched to
                            fprintf(eyeTrackDataFile, eyeLabel, ...
                                'c','1',num2str(sID), num2str(trial), num2str(prevlookLeft), num2str(fixStart), num2str(fixEnd));
                            eyeDataChoice = [eyeDataChoice; 1, sID, trial, prevlookLeft, fixStart, fixEnd];
                            
                        end
                        
                        
                        % update the variables
                        prevlookLeft = lookLeft;
                        
                    end
                end
                
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
                        rt1(trial,1)   = GetSecs-trialOnset;
                        
                        
                        if trackEye == 1
                            Eyelink('Message', 'choiceMade_train1');
                            % If they made their choice while
                            % staring at an ROI, put down RT
                            % instead of "fixend" - also it says
                            % "prevLookLeft" here because we want
                            % what they were just looking at, not
                            % the roi they switched to
                            if status >= 0 % status is whether or not there's a new sample. Only do this if there's a new sample
                                fprintf(eyeTrackDataFile, eyeLabel, ...
                                    'c','1',num2str(sID), num2str(trial), num2str(prevlookLeft), num2str(fixStart), num2str(fixEnd));
                                eyeDataChoice = [eyeDataChoice; 1, sID, trial, prevlookLeft, fixStart, fixEnd];
                            end
                        end %if eyetracking
                        
                        
                        pause(3-rt1);
                        
                        
                        
                        
                    elseif(keyCode(rightKey) && (GetSecs-trialOnset)<=3)
                        boxID = 2;
                        
                        %                    DrawFormattedText(window, 'current options', 'center', yText(1), textColor);
                        %                     Screen('FillRect', window, colorSquares(trialMat(trial,1),:), rectC(1,:));
                        %                     Screen('FillRect', window, colorSquares(trialMat(trial,2),:), rectC(2,:));
                        %                     Screen('DrawTexture', window, tex{trialMat(trial,1),1},[], rectC(1,:), 0, [], 1);
                        %                     Screen('DrawTexture', window, tex{trialMat(trial,2),1},[], rectC(2,:), 0, [], 1);
                        %                     Screen('FrameRect',window,frameColor,rectC(boxID,:),6);
                        %                     Screen('Flip',window);
                        rt1(trial,1)   = GetSecs-trialOnset;
                        
                        if trackEye == 1
                            Eyelink('Message', 'choiceMade_train1');
                            % If they made their choice while
                            % staring at an ROI, put down RT
                            % instead of "fixend" - also it says
                            % "prevLookLeft" here because we want
                            % what they were just looking at, not
                            % the roi they switched to
                            if status >= 0 % status is whether or not there's a new sample. Only do this if there's a new sample
                                fprintf(eyeTrackDataFile, eyeLabel, ...
                                    'c','1',num2str(sID), num2str(trial), num2str(prevlookLeft), num2str(fixStart), num2str(fixEnd));
                                eyeDataChoice = [eyeDataChoice; 1, sID, trial, prevlookLeft, fixStart, fixEnd];
                            end
                        end %if eyetracking
                        
                        
                        
                        pause(3-rt1);
                    elseif(keyCode(escapeKey))
                        if trackEye == 1
                            Eyelink('CloseFile');
                            Eyelink('ReceiveFile');
                            Eyelink('ShutDown');
                        end
                        ListenChar(0); clear screen; error('User terminated script with ESCAPE key.')
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
                
            end % while
        end %trackEye
        
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
        
        
        %------------------------------------------------------------------
        % Display fixation
        %------------------------------------------------------------------
        
        % wait time = feedback timing - stimulus timing - 3s choice
        % presentation time
        
        %%%%%
        if trackEye == 1
            
            Eyelink('command', 'record_status_message "FIXCROSS1 TRIAL %d/%d"', trial, NTrial);
            
            % start recording eye position (preceded by a short pause so that
            % the tracker can finish the mode transition)
            % The paramerters for the 'StartRecording' call controls the
            % file_samples, file_events, link_samples, link_events availability
            Eyelink('Command', 'set_idle_mode');
            WaitSecs(0.05);
            Eyelink('StartRecording');
            % record a few samples before we actually start displaying
            % otherwise you may lose a few msec of data
            WaitSecs(0.1);
            eye_used = Eyelink('eyeavailable'); % get the eye that's tracked
            
            
            %             % used for syncing time, This will print "synctime_part1" into
            %             % the EDF at this timepoint. This is useful for keeping track
            %             % in the edf of when a trial starts, and what part of the study
            %             % it belongs to if you have several parts.
            %             Eyelink('Message', 'SYNCTIME_part1');
            %             % This is just telling us that the trial is valid (useful for
            %             % if you recallibrate during a trial, in which case this will
            %             % later be switched to 0 so you know to throw the trial out)
            %             Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 1);
            %
            %             eye_used = Eyelink('eyeavailable'); % get eye that's tracked
            %
            %             % Set up the two ROIs. The !V is asking for gaze/velocity
            %             % information, IAREA RECTANGLE is the type of roi (it's a
            %             % rectangle... but you can also change it to be IAREA ELLIPSE
            %             % for circles or IAREA FREEHAND for an irregular shape.)
            %             % The following "%d"s correspond to the following (for both
            %             % rectangle and ellipse):
            %             % 1: id #
            %             % 2 & 3: top left (x,y)
            %             % 4 & 5: bottom right (x,y)
            %             % For freehand, the first one is still id #, the following are
            %             % x,y coordinates for each outer portion of the ROI (only x,y pairs
            %             % have a comma between them).
            %             % The final %s corresponds to the label string you give your
            %             % ROI. The examples below are two ROIs for the left and right
            %             % side of the screen.
            %             Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, midX-50, midY-50, midX+50, midY+50, 'cross1');
            
            % PUT CODE FOR DRAWING THE FIXATION CROSS HERE
            Screen('FillRect', window, fixColor  , FixCross');
            %flip to screen
            [~, startTime, ~, ~]=Screen('Flip', window, 1);
            
            f1_1(trial) = isi_t2_1(trial,1)-isi_t1_1(trial,1)-3;
            fixDur = f1_1(trial);
            % if you want people to fixate the cross for a specific period
            % of time, here you go.
            while GetSecs - startTime <= fixDur % fixDur is how long you want them to fixate before moving on
                % Check if there's a new sample
                if Eyelink('NewFloatSampleAvailable') > 0
                    % get the sample in the form of an event structure
                    evt = Eyelink('NewestFloatSample');
                    if eye_used ~= -1 % do we know which eye to use yet?
                        % if we do, get current gaze position from sample
                        % x-coordinate
                        eyeX = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                        % y-coordinate
                        eyeY = evt.gy(eye_used+1);
                        % pupil area - not used in standard eye tracking
                        % really, but great if you're looking at pupil
                        % dilation or making sure they're not blinking
                        a = evt.pa(eye_used+1);
                        
                        % do we have valid data and is the pupil visible?
                        if eyeX ~= el.MISSING_DATA && eyeY ~= el.MISSING_DATA && a > 0
                            distFromFix = sqrt((eyeX - midX)^2 + (eyeY - midY)^2);
                        else
                            distFromFix = 99999; % if no eye is present,do not advance trial
                        end
                    end
                    
                    if distFromFix > fixThresh
                        startTime = GetSecs; % Do not advance the trial, and reset the counter
                    end
                    
                    % During this time, add in the option to recalibrate
                    FlushEvents('keydown');
                    [keyIsDown, ~, keyCode] = KbCheck(-1);
                    if keyIsDown
                        if keyCode(KbName('c'))
                            % Change this label to be 0, so you know the
                            % trial actually wasn't valid
                            Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 0);
                            Eyelink('StopRecording');
                            EyelinkDoTrackerSetup(el);
                            Eyelink('StartRecording');
                            Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 1);
                        end
                    end
                end
            end

        end
        %------------------------------------------------------------------
        % Display feedback screen
        %------------------------------------------------------------------
       
        if trackEye == 1
            
            Eyelink('command', 'record_status_message "FEEDBACK TRIAL %d/%d"', trial, NTrial);
            
            % start recording eye position (preceded by a short pause so that
            % the tracker can finish the mode transition)
            % The paramerters for the 'StartRecording' call controls the
            % file_samples, file_events, link_samples, link_events availability
            Eyelink('Command', 'set_idle_mode');
            WaitSecs(0.05);
            Eyelink('StartRecording');
            % record a few samples before we actually start displaying
            % otherwise you may lose a few msec of data
            WaitSecs(0.1);
            eye_used = Eyelink('eyeavailable'); % get the eye that's tracked
            
            
            % used for syncing time, This will print "synctime_part1" into
            % the EDF at this timepoint. This is useful for keeping track
            % in the edf of when a trial starts, and what part of the study
            % it belongs to if you have several parts.
            Eyelink('Message', 'SYNCTIME_part1');
            % This is just telling us that the trial is valid (useful for
            % if you recallibrate during a trial, in which case this will
            % later be switched to 0 so you know to throw the trial out)
            Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 1);
            
            eye_used = Eyelink('eyeavailable'); % get eye that's tracked
            
            % Set up the two ROIs. The !V is asking for gaze/velocity
            % information, IAREA RECTANGLE is the type of roi (it's a
            % rectangle... but you can also change it to be IAREA ELLIPSE
            % for circles or IAREA FREEHAND for an irregular shape.)
            % The following "%d"s correspond to the following (for both
            % rectangle and ellipse):
            % 1: id #
            % 2 & 3: top left (x,y)
            % 4 & 5: bottom right (x,y)
            % For freehand, the first one is still id #, the following are
            % x,y coordinates for each outer portion of the ROI (only x,y pairs
            % have a comma between them).
            % The final %s corresponds to the label string you give your
            % ROI. The examples below are two ROIs for the left and right
            % side of the screen.
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, rectC(3,1),rectC(3,2),rectC(3,3),rectC(3,4),'top');
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, rectC(4,1),rectC(4,2),rectC(4,3),rectC(4,4),'bottom');
            
          if(trial>2)  
            % draw textures for feedback
            if(~isnan(ansBox1(trial,1)))
                
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
            end
            
          elseif(trial==1)
              % last 1 squares are the feedback options
              % Screen('TextSize', window, 23);
              %         DrawFormattedText(window, 'chosen last round', xText(2), yText(2), textColor);
              Screen('DrawTexture', window, tex{ansSquare1(trial,1),1}, [],rectC(3,:), 0, [], 1);
              %         Screen('FillRect', window, colorSquares(ansSquare(trial,1),:), rectC(3,:));
              DrawFormattedText(window, ['+', num2str(ansPayNow1(trial,1))], 'center',yposN(1), numberColor);
          else % if trial == 2
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
          end
            % Define when the trial starts here, would look something like:
            [~, trialOnset, ~, ~, ~] = Screen('Flip', window, 1);
            
            lookTop = -1; % set default
            prevlookTop = -1; % set default
            fixEnd = 0; % set default
            
            while(GetSecs-trialOnset<2)
                
                % Check recording status
                err = Eyelink('CheckRecording');
                
                if err ~= 0
                    error('checkrecording problem, status: %d', err);
                end
                
                % Check for new sample
                status = Eyelink('NewFloatSampleAvailable');
                if status >= 0
                    evt = Eyelink('NewestFloatSample');
                    % if we do, get current gaze position from sample
                    eyeX = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                    eyeY = evt.gy(eye_used+1);
                    eyetime = GetSecs();
                    
                    % Make sure that the data you're recording is even
                    % valid
                    if eyeX ~= el.MISSING_DATA && eyeY ~= el.MISSING_DATA && evt.pa(eye_used+1)>0
                        
                        % See which ROI they're looking at - if they're
                        % outside of both roi's, it's counted as a -1
                        if eyeY >= rectC(3,2) && eyeY <= rectC(3,4)
                            if eyeX >= rectC(3,1) && eyeX <= rectC(3,3)
                                lookTop = 1;
                            end
                        elseif eyeY >= rectC(4,2) && eyeY <= rectC(4,4)
                            if eyeX >= rectC(4,1) && eyeX <= rectC(4,3)
                                lookTop = 0;
                            end
                        end
                        
                        % Making this independent of whether
                        % the eye was within an roi helps eliminate
                        % blinks
                        if lookTop >= 0
                            if lookTop ~= prevlookTop
                                newFixStart = eyetime - trialOnset; % finds when the new dwell started
                                if prevlookTop == -1
                                    fixStart = eyetime - trialOnset; %like a default
                                end
                            else
                                fixStart = newFixStart; % keeps the same fix start throughout the dwell
                                fixEnd = eyetime - trialOnset; %updates every time, until they change rois
                            end
                        end
                        
                        % If they have completed the dwell, aka changed
                        % rois, then write down the previous fixation
                        % information into your own file. This is so that
                        % you have all your information in two places (the
                        % EDF and your matlab output). Hopefully you'll
                        % only need one.
                        if lookTop ~= prevlookTop && fixEnd ~= 0
                            
                            % Write out the eye data to a file as it's
                            % happening, but only after a dwell has been
                            % completed  - also it says
                            % "prevlookTop" here because we want
                            % what they were just looking at, not
                            % the roi they switched to
                            fprintf(eyeTrackDataFile, eyeLabel, ...
                                'f','1',num2str(sID), num2str(trial), num2str(prevlookTop), num2str(fixStart), num2str(fixEnd));
                            eyeDataFeedback = [eyeDataFeedback; 1, sID, trial, prevlookTop, fixStart, fixEnd];
                            
                        end
                        
                        
                        % update the variables
                        prevlookTop = lookTop;
                        
                    end
                end
            end % while GetSecs-trialOnset<2
   
        end % eyeTrack
        
        %%%%%

    %----------------------------------------------------------------------
    % Display fixation cross (iti)
    %----------------------------------------------------------------------    
    
    if trackEye == 1
        
        Eyelink('command', 'record_status_message "FIXCROSS2 TRIAL %d/%d"', trial, NTrial);
        
        % start recording eye position (preceded by a short pause so that
        % the tracker can finish the mode transition)
        % The paramerters for the 'StartRecording' call controls the
        % file_samples, file_events, link_samples, link_events availability
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        % otherwise you may lose a few msec of data
        WaitSecs(0.1);
        eye_used = Eyelink('eyeavailable'); % get the eye that's tracked
        
        
%         % used for syncing time, This will print "synctime_part1" into
%         % the EDF at this timepoint. This is useful for keeping track
%         % in the edf of when a trial starts, and what part of the study
%         % it belongs to if you have several parts.
%         Eyelink('Message', 'SYNCTIME_part1');
%         % This is just telling us that the trial is valid (useful for
%         % if you recallibrate during a trial, in which case this will
%         % later be switched to 0 so you know to throw the trial out)
%         Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 1);
%         
%         eye_used = Eyelink('eyeavailable'); % get eye that's tracked
%         
%         % Set up the two ROIs. The !V is asking for gaze/velocity
%         % information, IAREA RECTANGLE is the type of roi (it's a
%         % rectangle... but you can also change it to be IAREA ELLIPSE
%         % for circles or IAREA FREEHAND for an irregular shape.)
%         % The following "%d"s correspond to the following (for both
%         % rectangle and ellipse):
%         % 1: id #
%         % 2 & 3: top left (x,y)
%         % 4 & 5: bottom right (x,y)
%         % For freehand, the first one is still id #, the following are
%         % x,y coordinates for each outer portion of the ROI (only x,y pairs
%         % have a comma between them).
%         % The final %s corresponds to the label string you give your
%         % ROI. The examples below are two ROIs for the left and right
%         % side of the screen.
%         Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, midX-50, midY-50, midX+50, midY+50, 'cross');
        
        % PUT CODE FOR DRAWING THE FIXATION CROSS HERE
        Screen('FillRect', window, fixColor  , FixCross');
        %flip to screen
        [~, startTime, ~, ~]=Screen('Flip', window, 1);
        
        f2_1(trial) = isi_t1_1(trial+1,1)-isi_t2_1(trial,1)-2;
        fixDur = f2_1(trial);
        % if you want people to fixate the cross for a specific period
        % of time, here you go.
        while GetSecs - startTime <= fixDur % fixDur is how long you want them to fixate before moving on
            % Check if there's a new sample
            if Eyelink('NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                evt = Eyelink('NewestFloatSample');
                if eye_used ~= -1 % do we know which eye to use yet?
                    % if we do, get current gaze position from sample
                    % x-coordinate
                    eyeX = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                    % y-coordinate
                    eyeY = evt.gy(eye_used+1);
                    % pupil area - not used in standard eye tracking
                    % really, but great if you're looking at pupil
                    % dilation or making sure they're not blinking
                    a = evt.pa(eye_used+1);
                    
                    % do we have valid data and is the pupil visible?
                    if eyeX ~= el.MISSING_DATA && eyeY ~= el.MISSING_DATA && a > 0
                        distFromFix = sqrt((eyeX - midX)^2 + (eyeY - midY)^2);
                    else
                        distFromFix = 99999; % if no eye is present,do not advance trial
                    end
                end
                
                if distFromFix > fixThresh
                    startTime = GetSecs; % Do not advance the trial, and reset the counter
                end
                
                % During this time, add in the option to recalibrate
                FlushEvents('keydown');
                [keyIsDown, ~, keyCode] = KbCheck(-1);
                if keyIsDown
                    if keyCode(KbName('c'))
                        % Change this label to be 0, so you know the
                        % trial actually wasn't valid
                        Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 0);
                        Eyelink('StopRecording');
                        EyelinkDoTrackerSetup(el);
                        Eyelink('StartRecording');
                        Eyelink('Message', '!V TRIAL_VAR VALID_TRIAL %d', 1);
                    end
                end
            end
        end
    end
    

end %trial loop

        if trackEye == 1
           Eyelink('StopRecording');
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


trialN = size(isi_t2_2,1);
  
for trial = 1:3 %trialN
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

%%%%% close eyetracker

    if trackEye == 1
        Eyelink('CloseFile');
        Eyelink('ReceiveFile');
        Eyelink('ShutDown');
    end

%%%%%

% totalPay = ceil(totalPay1+totalPay2);
% display(totalPay)

% totalMoney = ceil(totalMoney1+totalMoney2) + 5;
display(totalMoney);


catch err  %quit out if the code encounters an error
    disp('There is an ERROR!');
    matlabfile = ['matlab' num2str(sID) '.mat'];
    save(matlabfile);
    
    ListenChar(0); %restore normal keyboard use
    
    sca
    fclose('all');
    if trackEye == 1
        Eyelink('CloseFile');
        Eyelink('ReceiveFile');
        Eyelink('ShutDown');
    end
    ShowCursor;
    Screen('CloseAll')
    rethrow(err);
end
