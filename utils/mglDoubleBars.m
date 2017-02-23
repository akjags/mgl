% mglDoubleBars.m
%
%        usage: mglDoubleBars(numBlocks)
%           by: akshay jagadeesh
%         date: 12/15/2016
%      purpose: Displays a retinotopy stimulus with two flashing checkerboard bars.
%
%

function [myscreen,stimImage] = mglDoubleBars(blockNum, frameGrab)

global xLimit yLimit stimulus;

% Initialize the screen
if ~ieNotDefined('frameGrab')
  myscreen.background = 'black';
else
  myscreen.background = 'gray';
end
%myscreen.displayName = 'fMRIprojFlex'; % Set to stimscreen for framegrabbing
myscreen.displayName = 'stimscreen';
myscreen = initScreen(myscreen);
xLimit = myscreen.imageWidth / 2;
yLimit = myscreen.imageHeight / 2;

stimulus = [];
% Init Stimulus
myscreen = initStimulus('stimulus', myscreen);

%stimfile = getLastStimfile(myscreen);
%keyboard;

% Simulate ticks in run
%mglSimulateRun(5, 10, 10);

%% Frame Grab: set to 1 (and change screen above to stimscreen) to create 192x108 stimimage
if ~ieNotDefined('frameGrab')
  stimulus.frameGrab = 1;
else
  stimulus.frameGrab = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Task 1: Fixation Stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('framegrab is: %d', stimulus.frameGrab));
if stimulus.frameGrab == 1
  disp('frame grab is on');
  stimulusTaskNum = 1;
  myscreen.background = 'black';
  stimulus.stimImage = NaN(myscreen.screenWidth, myscreen.screenHeight, 3, 960);
else
  fixTaskNum = 1;
  stimulusTaskNum = 2;

  fixStimulus.diskSize = 0.5;
  fixStimulus.fixWidth = 1;
  fixStimulus.fixLineWidth = 3;
  fixStimulus.stimTime = 0.5;
  fixStimulus.stimTime = 0.5;
  fixStimulus.responseTime = 1;
  fixStimulus.pos = [0 0];

  % Set the first task to be the fixation staircase task
  [task{fixTaskNum} myscreen] = fixStairInitTask(myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the stimulus
stimulus.fixedRandom = 0;

% Set stimulus parameters
stimulus.numCycles = 10;
stimulus.initialHalfCycle = 0;
stimulus.stepsPerCycle = 24;
if stimulus.frameGrab == 1
  stimulus.stimulusPeriod = 5; 
else
  stimulus.stimulusPeriod = 48;
end
stimulus.xOffset = 0; stimulus.yOffset = 0;
stimulus.barWidth = 2;
stimulus.imageWidth = myscreen.imageWidth;
stimulus.imageHeight = myscreen.imageHeight;

%initialize bar rotation matrix
stimulus.barRotMatrix1 = [cos(0) sin(0); -cos(0) sin(0)];
stimulus.barRotMatrix2 = [cos(0) sin(0); -cos(0) sin(0)];

%counter of number of stimulus segments
stimulus.segNum = 0;

% Task 2 is the retinotopy double bars
task{stimulusTaskNum}{1}.waitForBacktick = 1;

%task{stimulusTaskNum}{1}.seglen(1:stimulus.stepsPerCycle) = 4*(stimulus.stimulusPeriod / stimulus.stepsPerCycle)/5;
%task{stimulusTaskNum}{1}.parameter.bar1Angle = 0:45:359;
%task{stimulusTaskNum}{1}.parameter.bar2Angle = 0:45:359;
task{stimulusTaskNum}{1}.parameter.contrast = 1:8;
%task{stimulusTaskNum}{1}.randVars.calculated.elementAngle = nan;
task{stimulusTaskNum}{1}.numSegs = stimulus.stepsPerCycle;
task{stimulusTaskNum}{1}.numTrials = stimulus.numCycles;
task{stimulusTaskNum}{1}.random = 1;
if ~ieNotDefined('synchEveryVol') && synchEveryVol == 1
  disp('Synching to Vol after every segment');
  task{stimulusTaskNum}{1}.seglen(1:stimulus.stepsPerCycle) = 4*(stimulus.stimulusPeriod / stimulus.stepsPerCycle)/5;
  task{stimulusTaskNum}{1}.synchToVol(1:stimulus.stepsPerCycle) = 1;
else
  disp('Only synching to vol after last segment');
  task{stimulusTaskNum}{1}.synchToVol(1:stimulus.stepsPerCycle) = 0;
  task{stimulusTaskNum}{1}.seglen(1:stimulus.stepsPerCycle) = (stimulus.stimulusPeriod / stimulus.stepsPerCycle);
  if stimulus.frameGrab == 0
    task{stimulusTaskNum}{1}.seglen(stimulus.stepsPerCycle) = 4*(stimulus.stimulusPeriod / stimulus.stepsPerCycle)/5;
    task{stimulusTaskNum}{1}.synchToVol(stimulus.stepsPerCycle) = 1;
  end
end

%%%%% Set up Conditions: 8 directions, 4 single bar contrasts, 3 double bar contrast combinations.
% Set up bar directions
dirs = 0:45:359; 
dirs1 = repmat(dirs, 4,1);
singleContrasts = [.0625 0; .125 0; .675 0; 1 0];
singleContrastsDirs = repmat(singleContrasts, 8,1);
singleBarContrastWithDirs = [singleContrastsDirs dirs1(:) repmat(-1, 32, 1)]; % size = (32, 4)

doubleContrasts = [.25 .75; .5 .5; .75 .25];
doubleContrastDirs = repmat(doubleContrasts, 28, 1); % length: 84
dirBar1 = [repmat(0, 21, 1); repmat(45, 18, 1); repmat(90, 15, 1); repmat(135, 12, 1); repmat(180, 9, 1); repmat(225, 6, 1); repmat(270, 3, 1)];
d = repmat(dirs(2:end), 3, 1); d = d(:); dirBar2 = [d; d(4:end); d(7:end); d(10:end); d(13:end); d(16:end); d(19:end)]; % length: 84
doubleBarContrastWithDirs = [doubleContrastDirs dirBar1 dirBar2]; % size = (84, 4)

%% 116 total conditions: first 32 in array are single bars (bar1 contrast = 0), last 84 are double bars
stimulus.conditions = [singleBarContrastWithDirs; doubleBarContrastWithDirs]; 
numConditions = size(stimulus.conditions,1);
task{stimulusTaskNum}{1}.parameter.conditionNum = 1:numConditions;

% Set up randomization
stimulus.blockNum = blockNum;

%Add traces
[task{stimulusTaskNum}{1} myscreen] = addTraces(task{stimulusTaskNum}{1},myscreen,'maskPhase');
[task{stimulusTaskNum}{1} myscreen] = addTraces(task{stimulusTaskNum}{1},myscreen,'blank');

% Init retinotopy stimulus
stimulus = initRetinotopyStimulus(stimulus,myscreen);
stimulus.cycleTime = mglGetSecs;

% Initialize Task
[task{stimulusTaskNum}{1} myscreen] = initTask(task{stimulusTaskNum}{1}, myscreen, @startSegmentCallback, @updateScreenCallback, [], @startTrialCallback);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phaseNum = 1;
while (phaseNum <= length(task{stimulusTaskNum})) && ~myscreen.userHitEsc
  % update the task
  [task{stimulusTaskNum} myscreen phaseNum] = updateTask(task{stimulusTaskNum},myscreen,phaseNum);
  if stimulus.frameGrab == 0
    [task{fixTaskNum} myscreen] = updateTask(task{fixTaskNum}, myscreen, 1);
  end
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
if isfield(stimulus, 'stimImage')
  stimImage = stimulus.stimImage;
else
  stimImage = -1;
end
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Task Structure Methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------
% startSegmentCallback
%     gets called at start of each segment
%--------------------------------------------
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

segNum = task.thistrial.thisseg;

%Draw stencils for bar position at start of each segment
mglStencilCreateBegin(1);
x = stimulus.xBar1(:,segNum);
y = stimulus.yBar1(:,segNum);
coords1 = stimulus.barRotMatrix1*[x(:,1)'; y(:,1)'];
mglQuad(coords1(1,:)', coords1(2,:)', [0; 0; 0]);
mglStencilCreateEnd;
mglClearScreen;

mglStencilCreateBegin(2);
x = stimulus.xBar2(:,segNum);
y = stimulus.yBar2(:,segNum);
coords2 = stimulus.barRotMatrix2*[x(:,1)';y(:,1)'];
mglQuad(coords2(1,:)', coords2(2,:)', [0;0;0]);
mglStencilCreateEnd;
mglClearScreen;

mglStencilCreateBegin(3);
bar1 = [coords1(1,:)' coords1(2,:)'];
bar2 = [coords2(1,:)' coords2(2,:)'];
% if 0 or 180, vertical
if abs(stimulus.bar1Angle - stimulus.bar2Angle) == 180 || stimulus.bar2Angle == -1
  % do nothing
  disp('no overlap drawn');
elseif stimulus.bar1Angle == 0 %|| stimulus.bar1Angle == 180)
  calcDrawOverlap(bar1, bar2, 1);
elseif stimulus.bar1Angle == 180
  bar1a = [bar1(3:4,:); bar1(1:2,:)];
  calcDrawOverlap(bar1a, bar2, 1);
elseif stimulus.bar2Angle == 0 %|| stimulus.bar2Angle == 180)
  calcDrawOverlap(bar1, bar2, 2);
elseif stimulus.bar2Angle == 180
  bar2a = [bar2(3:4,:); bar2(1:2,:)];
  calcDrawOverlap(bar1, bar2a, 2);
else
  calcDrawOverlap(bar1,bar2);
end
mglStencilCreateEnd;
mglClearScreen;

stimulus.segNum = stimulus.segNum+1;

%-------------------------------------------
% startTrialCallback
%     gets called at start of each trial
%-------------------------------------------
function [task myscreen] = startTrialCallback(task, myscreen)

global stimulus;

numConditions = size(stimulus.conditions,1);

% Get trial number
seed = 4084807311;
rng('default');
rng(seed);
conditionOrder = randsample(numConditions, numConditions, false);

% Set up bar directions at start of trial
task.thistrial.conditionNum = conditionOrder((stimulus.blockNum-1)*stimulus.numCycles + task.trialnum); 
disp(sprintf('----Condition Number: %d-----', task.thistrial.conditionNum));
trialCondition = stimulus.conditions(task.thistrial.conditionNum, :);
stimulus.contrast1 = trialCondition(1);
stimulus.contrast2 = trialCondition(2);
stimulus.bar1Angle = trialCondition(3);
stimulus.bar2Angle = trialCondition(4);

%task.thistrial.bar1Angle = dirs(bars(1), randsample(2,1));
%task.thistrial.bar2Angle = dirs(bars(2), randsample(2,1));
c = cos(pi*stimulus.bar1Angle/180);
s = sin(pi*stimulus.bar1Angle/180);
stimulus.barRotMatrix1 = [c s;-s c];

c = cos(pi*stimulus.bar2Angle/180);
s = sin(pi*stimulus.bar2Angle/180);
stimulus.barRotMatrix2 = [c s; -s c];

% Set up bar contrasts at start of trial
%contrastTable = [.1 .4; .2 .5; .25 .35; .5 .2; .5 .5; .1 .9; .9 .1; .4 .6];
%stimulus.contrast1 = contrastTable(task.thistrial.contrast, 1);
%stimulus.contrast2 = contrastTable(task.thistrial.contrast, 2);
disp(sprintf('Bar 1 Contrast: %.01f%%. Bar 2 Contrast: %d%%', 100*stimulus.contrast1, 100*stimulus.contrast2));
disp(sprintf('Bar 1 angle: %d. Bar 2 angle: %d', stimulus.bar1Angle, stimulus.bar2Angle));


%--------------------------------------------
% updateScreenCallback
%       gets called to draw stimulus each frame
%--------------------------------------------
function [task myscreen] = updateScreenCallback(task,myscreen)

global stimulus;
stimulus = updateRetinotopyStimulus(stimulus,myscreen);

%--------------------------------------------
% initRetinotopyStimulus
%     initializes bars paramters
%--------------------------------------------
function stimulus = initRetinotopyStimulus(stimulus, myscreen)

stimulus.barHeight = stimulus.imageWidth*2;
stimulus.barMaskWidth = stimulus.imageWidth*1.5;
stimulus.barSweepExtent = 'max';

barDirVec = abs(stimulus.barRotMatrix1*[1 0]');
sweepExtent = max(barDirVec(1)*stimulus.imageWidth, barDirVec(2)*stimulus.imageHeight);
%sweepExtent = min(stimulus.imageWidth, stimulus.imageHeight);
stepSize = sweepExtent / (stimulus.stepsPerCycle-1);
stimulus.barCenter1 = [];
if isodd(stimulus.stepsPerCycle)
  stimulus.barCenter1(:,1) = -stepSize*(stimulus.stepsPerCycle-1)/2:stepSize:stepSize*(stimulus.stepsPerCycle-1);
else
  stimulus.barCenter1(:,1) = -stepSize*stimulus.stepsPerCycle/2+stepSize/2:stepSize:stepSize*stimulus.stepsPerCycle/2-stepSize/2;
end
stimulus.barCenter1(:,2) = 0;
disp(sprintf('(mglRetinotopy) barCenter1: %s', num2str(stimulus.barCenter1(:,1)', '%0.02f ')));

barDirVec = abs(stimulus.barRotMatrix2*[1 0]');
sweepExtent = max(barDirVec(1)*stimulus.imageWidth, barDirVec(2)*stimulus.imageHeight);
stepSize = sweepExtent / (stimulus.stepsPerCycle-1);
stimulus.barCenter2 = [];
if isodd(stimulus.stepsPerCycle)
  stimulus.barCenter2(:,1) = -stepSize*(stimulus.stepsPerCycle-1)/2:stepSize:stepSize*(stimulus.stepsPerCycle-1);
else
  stimulus.barCenter2(:,1) = -stepSize*stimulus.stepsPerCycle/2+stepSize/2:stepSize:stepSize*stimulus.stepsPerCycle/2-stepSize/2;
end
stimulus.barCenter2(:,2) = 0;
disp(sprintf('mglRetinotopy) barCenter2: %s', num2str(stimulus.barCenter2(:,1)', '%0.02f ')));

% Compute the coordinates of the two bars
%
stimulus.xBar1 = []; stimulus.yBar1 = [];
for i = 1:size(stimulus.barCenter1,1)
  xC = stimulus.barCenter1(i,1);
  yC = stimulus.barCenter1(i,2);
  w = stimulus.barWidth/2; h = stimulus.barHeight/2;
  stimulus.xBar1(:, end+1) = [xC-w; xC-w; xC+w; xC+w];
  stimulus.yBar1(:, end+1) = [yC-h; yC+h; yC+h; yC-h];
end
stimulus.xBar2 = []; stimulus.yBar2 = [];
for i = 1:size(stimulus.barCenter2,1)
  xC = stimulus.barCenter2(i,1);
  yC = stimulus.barCenter2(i,1);
  w = stimulus.barWidth/2; h = stimulus.barHeight/2;
  stimulus.xBar2(:,end+1) = [xC-w; xC-w; xC+w; xC+w];
  stimulus.yBar2(:,end+1) = [yC-h; yC+h; yC+h; yC-h];
end
stimulus.xIntersect = []; stimulus.yIntersect = [];
for i = 1:size(stimulus.xBar2,2)
  x1 = stimulus.xBar1(:,i); y1 = stimulus.yBar1(:,i);
  x2 = stimulus.xBar2(:,i); y2 = stimulus.yBar2(:,i);

end
%stimulus.nRect = length(allPhases1);
stimulus.phaseNumRect = 1;
stimulus.xRectOffset = 0;
stimulus.phaseNum = 1;
stimulus.frameCount = 1;
stimulus.frameTime = 0;

%--------------------------------------------
% updateRetinotopyStimulus
%     draws retinotopy stimulus to screen
%--------------------------------------------
function stimulus = updateRetinotopyStimulus(stimulus, myscreen)
%stimulus.t = mglGetSecs;
if stimulus.frameGrab == 1

  xScrn = [-myscreen.imageWidth/2 -myscreen.imageWidth/2 myscreen.imageWidth/2 myscreen.imageWidth/2];
  yScrn = [-myscreen.imageHeight/2 myscreen.imageHeight/2 myscreen.imageHeight/2 -myscreen.imageHeight/2];

  mglClearScreen;
  mglStencilSelect(1);
  mglPolygon(xScrn, yScrn, stimulus.contrast1);
  mglStencilSelect(2);
  mglPolygon(xScrn, yScrn, stimulus.contrast2);
  mglStencilSelect(3);
  mglPolygon(xScrn, yScrn, stimulus.contrast1+stimulus.contrast2);
  mglStencilSelect(0);
  
  stimIm = mglFrameGrab;
  frameIndex = 4*(stimulus.segNum-1)+1;
  %disp(sprintf('Grabbing frames %d through %d', frameIndex, frameIndex+3));
  stimulus.stimImage(:,:,:,frameIndex) = stimIm;
  stimulus.stimImage(:,:,:,frameIndex+1) = stimIm;
  stimulus.stimImage(:,:,:,frameIndex+2) = stimIm;
  stimulus.stimImage(:,:,:,frameIndex+3) = stimIm;
end

if stimulus.frameGrab==0 %&& mglGetSecs(stimulus.frameTime) >= 2;
  i = stimulus.frameCount; 
  mglClearScreen;
  mglStencilSelect(1);
  drawCheckerboard(stimulus.contrast1, i);
  mglStencilSelect(2);
  drawCheckerboard(stimulus.contrast2, i);
  mglStencilSelect(3);
  drawCheckerboard(stimulus.contrast1+stimulus.contrast2, i);
  mglStencilSelect(0);

  if mglGetSecs(stimulus.frameTime)>=0.25
    stimulus.frameCount = stimulus.frameCount+1;
    stimulus.frameTime = mglGetSecs;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             Helper Methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------
% calcDrawOverlap
%     Calculates the coordinates of the area of overlap between two bars, 
%     and draws it using mglPolygon
% Arguments: 
%        bar1, bar2 --> coordinates arrays in the form [x1 y1; x2 y2; x3 y3; x4 y4]
%--------------------------------------------
function calcDrawOverlap(bar1, bar2, vert)

if ieNotDefined('vert')
  vert = 0;
end

if vert ~= 1
  % for each bar, first calculate the line (y=mx+b) between points 1 and 2, and points 3 and 4 (long sides)
  [rB1(1), rB1(2)] = calculateLine(bar1(1,:), bar1(2,:));
  [fB1(1), fB1(2)] = calculateLine(bar1(3,:), bar1(4,:));
end

if vert ~= 2
  [rB2(1), rB2(2)] = calculateLine(bar2(1,:), bar2(2,:));
  [fB2(1), fB2(2)] = calculateLine(bar2(3,:), bar2(4,:));
end

if vert == 1
  x1 = bar1(1,1); x2 = bar1(3,1);
  p1(1) = x1; p1(2) = rB2(1)*x1 + rB2(2);
  p2(1) = x2; p2(2) = rB2(1)*x2 + rB2(2);
  p3(1) = x2; p3(2) = fB2(1)*x2 + fB2(2);
  p4(1) = x1; p4(2) = fB2(1)*x1 + fB2(2);
elseif vert == 2
  x1 = bar2(1,1); x2 = bar2(3,1);
  p1(1) = x1; p1(2) = rB1(1)*x1 + rB1(2);
  p2(1) = x2; p2(2) = rB1(1)*x2 + rB1(2);
  p3(1) = x2; p3(2) = fB1(1)*x2 + fB1(2);
  p4(1) = x1; p4(2) = fB1(1)*x1 + fB1(2);
else
  [p1(1), p1(2)] = calculateIntersect(rB1, rB2);
  [p2(1), p2(2)] = calculateIntersect(rB1, fB2);
  [p3(1), p3(2)] = calculateIntersect(fB1, fB2);
  [p4(1), p4(2)] = calculateIntersect(fB1, rB2);
end

intersect = [p1; p2; p3; p4];
mglPolygon(intersect(:,1)', intersect(:,2)', 0);

%--------------------------------------------
% calculateIntersect
%        Calculates point of intersection between two lines
%   Arguments: 
%        l1, l2 --> lines in the form [slope, intersect]
%--------------------------------------------
function [intersectX,intersectY] = calculateIntersect(l1, l2)
intersectX = (l2(2) - l1(2)) / (l1(1) - l2(1));
intersectY = (l1(1)*l2(2) - l2(1)*l1(2)) / (l1(1) - l2(1));

%--------------------------------------------
% calculateLine
%     Calculates a line between two points
%   Arguments: 
%        p1, p2 --> points in the form [x,y]
%--------------------------------------------
function [slope,intercept] = calculateLine(p1, p2)
slope = (p2(2) - p1(2)) / (p2(1) - p1(1));
intercept = p2(2) - slope*p2(1);

%--------------------------------------------
% drawCheckerBoard
%     Draws a checkerboard to fill the screen using mglQuad
%     Starting in the top left, draws alternating black and white squares.
%     Then imposes a gray (0.50) semi-transparent window
%   Arguments: 
%        contrast --> 1-alpha : specifies transparency of window, i.e. contrast of checkerboard. 
%        index --> 0 or 1 specifies which color to start with (for flashing effect)
%--------------------------------------------
function drawCheckerboard(contrast, index, squaresize)

global myscreen xLimit yLimit;
if ieNotDefined('squaresize')
  squaresize = xLimit/20;
end

% Get array of x coordinates and y coordinates for squares
xc = -xLimit:squaresize:(xLimit+2*squaresize); 
yc = -yLimit:squaresize:(yLimit+2*squaresize);
if mod(length(xc),2) == 1
  xc = xc(1:end-1);
end
if mod(length(yc),2) == 1
  yc = yc(1:end-1);
end
numYsq = length(yc)-1;
numXsq = length(xc)-1;
yCoords1 = reshape(repmat(yc(1:end-1)', 1, numXsq)', 1, numYsq*numXsq);
yCoords2 = reshape(repmat(yc(2:end)', 1, numXsq)', 1, numYsq*numXsq);
xCoords1 = repmat(xc(1:end-1), 1, numYsq);
xCoords2 = repmat(xc(2:end), 1, numYsq);

% Get color array of squares
if mod(index,2) == 1
  colorMatrix = repmat([1 0; 1 0; 1 0], 1, numXsq*numYsq);
elseif mod(index,2) == 0
  colorMatrix = repmat([0 1; 0 1; 0 1], 1, numXsq*numYsq);
end
colorMatrix = colorMatrix(:,1:numXsq*numYsq);

mglQuad([xCoords1; xCoords1; xCoords2; xCoords2],[yCoords1; yCoords2; yCoords2; yCoords1], colorMatrix); 
mglPolygon([-xLimit, xLimit, xLimit, -xLimit], [yLimit, yLimit, -yLimit, -yLimit], [0.5 0.5 0.5 1-contrast]);
