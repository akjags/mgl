% mglDoubleBars.m
%
%        usage: mglDoubleBars(numBlocks)
%           by: akshay jagadeesh
%         date: 12/15/2016
%      purpose: Displays a retinotopy stimulus with two flashing checkerboard bars.
%
%

function [myscreen,stimImage] = mglDoubleBars(numBlocks)

global myscreen xLimit yLimit stimulus;

% Initialize the screen
myscreen = initScreen('stimscreen');
myscreen.background = 'gray';
xLimit = myscreen.imageWidth / 2;
yLimit = myscreen.imageHeight / 2;

% Init Stimulus
myscreen = initStimulus('stimulus', myscreen);

%% Frame Grab: set to 1 (and change screen above to stimscreen) to create 192x108 stimimage
stimulus.frameGrab = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Task 1: Fixation Stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('framegrab is: %d', stimulus.frameGrab));
if stimulus.frameGrab == 1
  disp('frame grab is on');
  stimulusTaskNum = 1;
  myscreen.background = 'black';
else
  mglClearScreen(0.5); mglFlush; mglClearScreen(0.5);
  myscreen.background = 'gray';
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
stimulus.numCycles = 20;
stimulus.initialHalfCycle = 0;
stimulus.stepsPerCycle = 24;
stimulus.stimulusPeriod = 24;
stimulus.xOffset = 0; stimulus.yOffset = 0;
stimulus.barWidth = 4;
stimulus.imageWidth = myscreen.imageWidth;
stimulus.imageHeight = myscreen.imageHeight;

%initialize bar rotation matrix
stimulus.barRotMatrix1 = [cos(0) sin(0); -cos(0) sin(0)];
stimulus.barRotMatrix2 = [cos(90) sin(90); -cos(90) sin(90)];

%counter of number of stimulus segments
stimulus.segNum = 0;

% Task 2 is the retinotopy double bars
task{stimulusTaskNum}{1}.waitForBacktick = 1;
task{stimulusTaskNum}{1}.seglen(1:stimulus.stepsPerCycle) = repmat(stimulus.stimulusPeriod / stimulus.stepsPerCycle, 1, stimulus.stepsPerCycle);
task{stimulusTaskNum}{1}.parameter.bar1Angle = 0:45:359;
task{stimulusTaskNum}{1}.parameter.bar2Angle = 0:45:359;
task{stimulusTaskNum}{1}.parameter.contrast = 1:8;
%task{stimulusTaskNum}{1}.randVars.calculated.elementAngle = nan;
task{stimulusTaskNum}{1}.numSegs = stimulus.stepsPerCycle;
task{stimulusTaskNum}{1}.numTrials = stimulus.numCycles;
task{stimulusTaskNum}{1}.random = 1;
task{stimulusTaskNum}{1}.synchToVol(1:stimulus.stepsPerCycle) = 0;

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
stimImage = stimulus.stimImage;
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
if(stimulus.bar1Angle == 0 || stimulus.bar1Angle == 180)
  calcDrawOverlap(bar1, bar2, 1);
elseif(stimulus.bar2Angle == 0 || stimulus.bar2Angle == 180)
  calcDrawOverlap(bar1, bar2, 2);
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

% Set up bar directions at start of trial
%dirs = [0 180; 45 225; 90 270; 135 315];
dirs = [45 225; 90 270; 135 315];
bars = randsample(size(dirs,1),2,false);
%bars = datasample(0:45:359, 2, 'Replace', false);
task.thistrial.bar1Angle = dirs(bars(1), randsample(2,1));
task.thistrial.bar2Angle = dirs(bars(2), randsample(2,1));
c = cos(pi*task.thistrial.bar1Angle/180);
s = sin(pi*task.thistrial.bar1Angle/180);
stimulus.bar1Angle = task.thistrial.bar1Angle;
stimulus.barRotMatrix1 = [c s;-s c];

c = cos(pi*task.thistrial.bar2Angle/180);
s = sin(pi*task.thistrial.bar2Angle/180);
stimulus.bar2Angle = task.thistrial.bar2Angle;
stimulus.barRotMatrix2 = [c s; -s c];

% Set up bar contrasts at start of trial
contrastTable = [.1 .4; .2 .5; .25 .35; .5 .2; .5 .5; .1 .9; .9 .1; .4 .6];
stimulus.contrast1 = contrastTable(task.thistrial.contrast, 1);
stimulus.contrast2 = contrastTable(task.thistrial.contrast, 2);
disp(sprintf('Bar 1 Contrast: %d%%. Bar 2 Contrast: %d%%', 100*stimulus.contrast1, 100*stimulus.contrast2));
disp(sprintf('Bar 1 angle: %d. Bar 2 angle: %d', task.thistrial.bar1Angle, task.thistrial.bar2Angle));


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
%sweepExtent = max(barDirVec(1)*stimulus.imageWidth, barDirVec(2)*stimulus.imageHeight);
sweepExtent = min(stimulus.imageWidth, stimulus.imageHeight);
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

%--------------------------------------------
% updateRetinotopyStimulus
%     draws retinotopy stimulus to screen
%--------------------------------------------
function stimulus = updateRetinotopyStimulus(stimulus, myscreen)

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

  stimulus.stimImage(:,:,:,stimulus.segNum) = mglFrameGrab;
end

if stimulus.frameGrab == 0
  i = 1 + mod(stimulus.frameCount,2);
  mglClearScreen;
  mglStencilSelect(1);
  drawCheckerboard(stimulus.contrast1, i);
  mglStencilSelect(2);
  drawCheckerboard(stimulus.contrast2, i);
  mglStencilSelect(3);
  drawCheckerboard(stimulus.contrast1+stimulus.contrast2, i);
  mglStencilSelect(0);

  stimulus.frameCount = stimulus.frameCount + 1;
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
  x1 = bar2(1,1); x2 = bar2(1,1);
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
