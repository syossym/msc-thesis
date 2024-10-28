function ipf(x,y)
% ipf(x,y) or ipf(DataMatrix)
% Keyboard-operated Interactive Peak Fitter for data vectors in 
% arguments x,y, or a single data matrix "DataMatrix", with 
% x values in row 1 and y values in row 2 (e.g. [x y]).
% This version uses keyboard commands only. Type "help ipf" for list.
% See http://www.wam.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
% T. C. O'Haver (toh@umd.edu).  Verion 4.1, April 18: Added Autozero
% ON/OFF notification at top, bug fixes, and "k" keypress. Version 4:
% Added keypress function "x" to refine the fit by performing 10 trial  
% fits with slightly different first guesses and take the one with the  
% lowest fittinig error. You can change the number of trials, "NumTrials",  
% in line 69 (the default is 10)
% Example 1: x=[0:.005:1];y=humps(x).^3;ipf(x,y)
% Example 2: x=[-10:.1:10];y=exp(-(x).^2);ipf(x,y)
% For large data sets, to view only a portion of the data over a specified x-
% axis range, you can type n1=val2ind(x,x1);n2=val2ind(x,x2);ipf(x(n1:n2),y(n1:n2))
% where x1 and x2 are the end limits of the x-axis values displayed.
%
% Keyboard Controls:
% Pan signal left and right:  Coarse pan: < and >   
%                             Fine pan: left and right cursor arrow keys
% Zoom in and out:            Coarse zoom: / and '   
%                             Fine zoom: up and down cursor arrow keys
% Select # of peaks:          Number keys 1-5
% Select peak shape:          g Gaussian
%                             l Lorentzian
%                             o Logistic
%                             p Pearson (use a,z keys to adjust shape)
%                             e exponentially-broadened Gaussian 
%                                 (use a,z keys to adjust broadening)
% Fit:                        f
% Toggle autozero off/on      t   (Newly added as of version 3.5)
% Baseline:                   b, then click left and right baseline
% Custom start position:      c, then click on each peak position
% Adjust 'extra' up or down:  a,z
% Print parameters & results: q
% Print fit results only:     r
% Print out data segment:     d  (Newly added as of version 3.3)
% Refine fit                  x  (Newly added as of version 4)
% Prints list of commands     k  (Newly added as of version 4.1)

global X
global Y
global xo
global dx
global NumPeaks
global NumTrials
global Shape
global MeanFitError
global PEAKHEIGHTS
global FitResults
global start
global extra
global delta
global AUTOZERO
global miny
close
format short g
format compact
warning off all

    % Assign arguments to internal global variables
if nargin==1,
    X=x(:,1);
    Y=x(:,2);
else
X=x;
Y=y;
end


% Adjust X and Y vector shape to 1 x n (rather than n x 1)
X=reshape(X,1,length(X));
Y=reshape(Y,1,length(Y));

% X=[1:length(Y)]; % Use this only to create an x vector if needed

% Remove excess offset from data
miny=min(Y);

% Initial values of parameters
NumPeaks=1; % Initial Number of peaks in model
NumTrials=10; % Number of repeat fits when X key is pressed
Shape=1; % Initial Shape of the peaks (1=Gaussian, 2=Lorentzian, etc)
extra=1;
delta=0; % Initial Random change in initial start guesses for each re-fit
xo=length(Y)/2; % Initial Pan setting
dx=length(Y)/4; % Initial Zoom setting
start=length(Y)/2;
PEAKHEIGHTS=zeros(NumPeaks,1);
AUTOZERO=1; % Sets autozero operation. Press T to toggle AUTOZERO off and on.

% Plot the signal and its fitted components and residuals
RedrawSignal(X,Y,xo,dx,NumPeaks);
h=figure(1);
h2=gca;
% Attaches KeyPress test function to the figure.
set(gcf,'KeyPressFcn',@ReadKey)
uicontrol('Style','text')
% end of outer function
% ----------------------------SUBFUNCTIONS--------------------------------
function ReadKey(obj,eventdata)
% Interprets key presses from the Figure window.
global X
global Y
global xx
global yy
global xo
global dx
global start
global FitResults
global NumPeaks
global NumTrials
global Shape
global delta
global PEAKHEIGHTS
global AUTOZERO
global miny
global extra
global MeanFitError
% When a key is pressed, interprets key and calls corresponding function.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
key=get(gcf,'CurrentCharacter');
if ischar(key),
  switch double(key),
    case 28
        % Pans one point down when left arrow pressed.
          xo=xo-1;
          if xo<1,xo=1;,end
          [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
         % h2=gca;
    case 29
        % Pans one point up when right arrow pressed.
        ly=length(Y);
        xo=xo+1;
        if xo>ly,xo=ly;,end
        [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
    case 44
        % Pans 2% down when < key pressed.
        ly=length(Y);
        xo=xo-ly/50;
        if xo<1,xo=1;,end
        [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
    case 46
        % Pans 2% up when > key pressed.
        ly=length(Y);
        xo=xo+ly/50;
        if xo>ly,xo=ly;,end
        [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
    case 31
        % Zooms one point up when up arrow pressed.
        dx=dx+2;
        [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
    case 30
        % Zooms one point down when down arrow pressed.
        dx=dx-2;
        [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
    case 47
        % Zooms 2% up when / pressed.
        ly=length(Y);
        dx=dx+ly/50;
        [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
    case 39
        % Zooms 2% down when ' pressed.
        ly=length(Y);
        dx=dx-ly/50;
        [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);      
    case 102
        % When 'f' key is pressed, tweaks start values, computes and plots fit.
        startnow=start;
        delta=(max(xx)-min(xx))/100;
          for k=1:2:2*NumPeaks,
            startnow(k)=start(k)+randn*delta;
          end
        [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,startnow,extra); 
    case 99
        % When 'c' key is pressed, user clicks graph to enter start positons, 
        % then fit is computed and graph re-drawn.q
        % Acquire first-guess peak positions from user mouse pointer
        figure(1);subplot(2,1,1);xlabel('Click on the estimated positions of each proposed component peak.')
        [clickX,clickY] = GINPUT(NumPeaks);
        % Create a "start" vector using these peak positions, with peak widths equal 
        % to 1/10 of the zoom region.
        n=max(xx)-min(xx);
        width=n/(5*NumPeaks);
        start=[];
        for k=1:NumPeaks,
            start=[start clickX(k) width];
        end
        [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);    
    case 98
        % When 'b' key is pressed, user clicks graph before and after peaks
        % to enter background points, then fit is computed and graph re-drawn.
          figure(1);subplot(2,1,1);xlabel('Click on the baseline to the LEFT the peak(s).')
          [X1,Y1] = GINPUT(1);
          figure(1);subplot(2,1,1);xlabel('Now click on the baseline to the RIGHT the peak(s).')
          [X2,Y2] = GINPUT(1);
          n=length(xx);
          %  Create "start" for this number of peaks
          yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
          [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
    case {49,50,51,52,53,54}
        % When a number key is pressed, sets the number of peaks to that number.  
        n=key-48;
        ShapeString=''; 
        if round(n)~=NumPeaks,
            NumPeaks=round(n);
            switch Shape
              case 1
                ShapeString='Gaussian';
              case 2
                ShapeString='Lorentzian';    
              case 3
                ShapeString='logistic';
              case 4
                ShapeString='Pearson7';
              case 5        
                ShapeString='ExpGaussian';
              otherwise
                ShapeString='';
            end % switch
          subplot(2,1,2)
          xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
        end % if
        n=max(xx)-min(xx);
        % Create a start value for this number of peaks
        start=[];
        startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
        for marker=1:NumPeaks,
            markx=startpos(marker);
            start=[start markx n/(5*NumPeaks)];
        end
        RedrawSignal(X,Y,xo,dx,NumPeaks);
          xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
    case {103,108,111,112,101}
         % Selects peak shape when G,L,O,P or E kay is pressed
        switch key
            case 103 % When 'g' key is pressed, peak shape is set to Gaussian. 
                n=1;
            case 108 % When 'l' key is pressed, peak shape is set to Lorentzian. 
                n=2;
            case 111 % When 'o' key is pressed, peak shape is set to Logistic. 
                n=3;
            case 112 % When 'p' key is pressed, peak shape is set to Pearson.
                n=4;
            case 101 % When 'e' key is pressed, peak shape is set to Exponentally-broadened gaussian.
                n=5;
        end % switch
        if round(n)~=Shape,
        Shape=round(n);
          switch Shape
            case 1
                ShapeString='Gaussian';
            case 2
                ShapeString='Lorentzian';    
            case 3
                ShapeString='logistic';
            case 4
                ShapeString='Pearson';
            case 5        
                ShapeString='ExpGaussian';
            otherwise
          end % switch
          figure(1);subplot(2,1,2)
          xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
       end % if
    case 97
        % When 'a' key is pressed, increases "extra" by 0.1
        extra=extra+.1;
        if extra==0, extra=.01;,end
        [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
    case 122
        % When 'z' key is pressed, decreases "extra" by 0.1
        extra=extra-.1;
        if extra==0, extra=.01;,end
        [FitResults,MeagfnFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
    case 114
        % When 'r' key is pressed, prints out table of fit results
       disp('Peak#  Position       Height        Width         Area') 
         for peak=1:NumPeaks,
               disp(num2str(FitResults(peak,:)))
         end % for
         start
    case 100
        % When 'd' key is pressed, prints out table of data segment
    disp('     x          y') 
         disp(num2str([xx(1:length(xx))' yy(1:length(xx))']))
          case 116
        % When 't' key is pressed, toggles AUTOZERO mode
       if AUTOZERO,
           AUTOZERO=0;
           % Y=Y-miny;
           subplot(2,1,1);
           title('Autozero OFF. Pan and Zoom to isolate peaks to be fit in upper window.')
       else
           AUTOZERO=1;
           % Y=Y+miny;
           subplot(2,1,1);
           title('Autozero ON. Pan and Zoom to isolate peaks to be fit in upper window.')
       end
    case 107
        % When 'k' key is pressed, prints out table of keyboar commands
        disp('Keyboard Controls:')
        disp(' Pan signal left and right:  Coarse pan: < and >')   
        disp('                             Fine pan: left and right cursor arrow keys')
        disp(' Zoom in and out:            Coarse zoom: / and "  ') 
        disp('                             Fine zoom: up and down cursor arrow keys')
        disp(' Select # of peaks:          Number keys 1-5')
        disp(' Select peak shape:          g Gaussian')
        disp('                             l Lorentzian')
        disp('                             o Logistic')
        disp('                             p Pearson (use a,z keys to adjust shape)')
        disp('                             e exponentially-broadened Gaussian ')
        disp('                                 (use a,z keys to adjust broadening)')
        disp(' Fit:                        f')
        disp(' Toggle autozero OFF/ON      t   (Newly added as of version 3.5)')
        disp(' Baseline:                   b, then click left and right baseline')
        disp(' Custom start position:      c, then click on each peak position')
        disp(' Adjust "extra" up or down:  a,z')
        disp(' Print parameters & results: q')
        disp(' Print fit results only:     r')
        disp(' Print out data segment:     d  (Newly added as of version 3.3)')
        disp(' Refine fit                  x  (Newly added as of version 4)')
   case 120 % When 'x' key is pressed, calls peakfit to take best of 'NumTrials' trial fits
       center=(max(xx)+min(xx))/2;
       window=max(xx)-min(xx);          
       % Prints out peakfit functions with arguments in command window
       % disp(['figure(2);[FitResults]=peakfit([x,y],' num2str(center) ',' num2str(window) ',' num2str(NumPeaks) ',' num2str(Shape) ',' num2str(extra) ',10,[' num2str(start) ']) '])
       [FitResults]=peakfit([X',Y'],center,window,NumPeaks,Shape,extra,NumTrials,start)
       figure(1)
   case 113
        % When 'q' key is pressed, prints out fitting parameters
        switch Shape
           case 1
               ShapeString='Gaussian';
           case 2
               ShapeString='Lorentzian';    
           case 3
               ShapeString='Logistic';
           case 4
               ShapeString='Pearson';
           case 5        
               ShapeString='Exponentially-broadened Gaussian';
           otherwise
         end % switch
       disp('--------------------------------------------------')
       disp(['Peak Shape = ' ShapeString])
       disp(['Number of peaks = ' num2str(NumPeaks)])
       if Shape>=4, disp(['Extra = ' num2str(extra)]), end
       disp(['Fitted range = ' num2str(min(xx)) ' - ' num2str(max(xx)) ' (' num2str(max(xx)-min(xx)) ')  (' num2str((max(xx)+min(xx))/2) ')  ' ])
       disp(['Percent Error = ' num2str(MeanFitError)])
       disp('Peak#  Position       Height        Width         Area') 
       for peak=1:NumPeaks,
           disp(num2str(FitResults(peak,:)))
       end % for
   otherwise  
       UnassignedKey=double(key)
       disp('Press k to print out list of keyboard commands')
   end % switch
end % if
% ----------------------------------------------------------------------    
function [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
% Plots the entire signal (X,Y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys.
% Top half of the figure shows original signal
global AUTOZERO
Startx=round(xo-(dx/2));
Endx=round(xo+(dx/2)-1);
if Endx>length(Y),
    Endx=length(Y);
end
if Startx<1,
     Startx=1;
end
PlotRange=[Startx:Endx];
if PlotRange<5, PlotRange=[xo:xo+5];,end
xx=X(PlotRange);
yy=Y(PlotRange); 
% Remove comments from the next 5 lines to enable auto-zero operation
if AUTOZERO==1,
  X1=min(xx);
  X2=max(xx);
  Y1=mean(yy(1:length(xx)/10));
  Y2=mean(yy((length(xx)-length(xx)/10):length(xx)));
  yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
end % if
figure(1);subplot(2,1,1);plot(xx,yy,'b.'); % Plot the original signal in blue
hold on
% Mark starting peak positions with vertical dashed lines
% Determine locations for peak (vertical line) markers
n=max(xx)-min(xx);
start=[];
startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
for marker=1:NumPeaks,
    markx=startpos(marker);
    start=[start markx n/20];
    figure(1);subplot(2,1,1);plot([markx markx],[min(yy) max(yy)],'m--')
end % for
hold off
if AUTOZERO,
    title('Autozero ON. Pan and Zoom to isolate peaks to be fit in upper window.')
    axis([X(Startx) X(Endx) min(yy) max(yy)]);
else
    title('Autozero OFF. Pan and Zoom to isolate peaks to be fit in upper window.')
    % axis([X(Startx) X(Endx) min(Y) max(yy)]);
    axis([X(Startx) X(Endx) 0 max(yy)]);
end
xlabel('Line up the estimated peak positions roughly with the vertical lines')

% Bottom half of the figure shows full signal
figure(1);subplot(2,1,2);plot(X,Y)
hold on
for marker=1:NumPeaks,
    markx=startpos(marker);
    figure(1);subplot(2,1,2);plot([markx markx],[min(Y) max(Y)],'m--')
end % for
hold off
xlabel([' # Peaks=1-5; Press F to fit; Press K to print out keyboard commands'])
% ----------------------------------------------------------------------
function [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra)
% Given isolated segment of signal (xx,yy), plots it in upper half, computes fit with
% "NumPeaks" component peaks of shape "Shape", starting with start values
% "start", then plots residuals in lower half. 
%  T. C. O'Haver (toh@umd.edu),  Version 2, March 2008, uses unconstrained
%  fit.
global PEAKHEIGHTS
global AUTOZERO
%global extra
%global MeanFitError
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
x0=min(xx);
h=figure(1);
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'TypicalX',x0,'Display','off' );
switch Shape
    case 1
        FitParameters=fminsearch(@fitgaussian,start,options,xx,yy);
        ShapeString='Gaussian';
    case 2
        FitParameters=fminsearch(@fitlorentzian,start,options,xx,yy);
        ShapeString='Lorentzian';    
    case 3
        FitParameters=fminsearch(@fitlogistic,start,options,xx,yy);
        ShapeString='Logistic';
    case 4
        FitParameters=fminsearch(@fitpearson,start,options,xx,yy,extra);
        ShapeString='Pearson';
    case 5
        FitParameters=fminsearch(@fitexpgaussian,start,options,xx,yy,-extra);
        ShapeString='ExpGaussian';
    otherwise
end % switch

% Construct model from fitted parameters
A=zeros(NumPeaks,n);
AA=zeros(NumPeaks,100);
xxx=linspace(min(xx),max(xx));
for m=1:NumPeaks,
   switch Shape
    case 1
        A(m,:)=gaussian(xx,FitParameters(2*m-1),FitParameters(2*m));
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        A(m,:)=lorentzian(xx,FitParameters(2*m-1),FitParameters(2*m));
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        A(m,:)=logistic(xx,FitParameters(2*m-1),FitParameters(2*m));
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        A(m,:)=pearson(xx,FitParameters(2*m-1),FitParameters(2*m),extra);
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        A(m,:)=expgaussian(xx,FitParameters(2*m-1),FitParameters(2*m),-extra)';
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra)';
    otherwise
  end % switch
end % for
model=PEAKHEIGHTS'*A;  % Multiplies each row by the corresponding amplitude and adds them up
mmodel=PEAKHEIGHTS'*AA;
% Top half of the figure shows original signal and the fitted model.
figure(1);subplot(2,1,1);plot(xx,yy,'b.'); % Plot the original signal in blue dots
hold on
for m=1:NumPeaks,
    figure(1);plot(xxx,PEAKHEIGHTS(m)*AA(m,:),'g')  % Plot the individual component peaks in green lines
    area(m)=trapz(xx,PEAKHEIGHTS(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end
% Mark starting peak positions with vertical dashed lines
for marker=1:NumPeaks,
    markx=start((2*marker)-1);
    figure(1);subplot(2,1,1);plot([markx markx],[0 max(yy)],'m--')
end % for
figure(1);plot(xxx,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
hold off;
if AUTOZERO,
    axis([min(xx) max(xx) min(yy) max(yy)]);
else
    axis([min(xx) max(xx) 0 max(yy)]);
end

title('Pan and Zoom to select peaks; Press k to print out list of keyboard commands')
xlabel('Vertical dotted lines indicate first guess peak positions');
% Bottom half of the figure shows the residuals and displays RMS error between original signal and model
residual=yy-model;
figure(1);subplot(2,1,2);plot(xx,residual,'b.')
MeanFitError=100*norm(residual)./(sqrt(n)*max(yy));
xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(MeanFitError) '     Extra = ' num2str(extra) ] )
axis([min(xx) max(xx) min(residual) max(residual)]);
% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
for m=1:NumPeaks,
    if m==1,
       FitResults=[round(m) FitParameters(2*m-1) PEAKHEIGHTS(m) abs(FitParameters(2*m)) area(m)];
    else
       FitResults=[FitResults ; [round(m) FitParameters(2*m-1) PEAKHEIGHTS(m) abs(FitParameters(2*m)) area(m)]];
    end % if
end % for

% Display Fit Results on upper graph
figure(1);subplot(2,1,1);
residualaxis=axis;
text(residualaxis(1), 0.8*residualaxis(4) ,num2str(FitResults));
% ----------------------------------------------------------------------
function [FitResults,LowestError,BestStart]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start)

global PEAKHEIGHTS
global AUTOZERO
X=signal(:,1);
Y=signal(:,2);
X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
Y=reshape(Y,1,length(Y));
% Y=Y-min(Y);  % Remove excess offset from data
% Isolate desired segment from data set for curve fitting
if nargin==1,center=(max(X)-min(X))/2;window=max(X)-min(X);,end
xoffset=center-window;
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
xx=X(n1:n2)-xoffset;
yy=Y(n1:n2);

% Remove linear baseline from data segment
if AUTOZERO==1,
  X1=min(xx);
  X2=max(xx);
  Y1=mean(yy(1:length(xx)/10));
  Y2=mean(yy((length(xx)-length(xx)/10):length(xx)));
  yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
end % if

% Define values of any missing arguments
switch nargin
    case 1
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
    case 4
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
    case 5
        extra=0;
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
    case 6
        NumTrials=1;
        start=calcstart(xx,NumPeaks,xoffset);
    case 7
        start=calcstart(xx,NumPeaks,xoffset);
    otherwise
end % switch nargin

PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;
switch NumPeaks
    case 1
        newstart(1)=start(1)-xoffset;
    case 2
        newstart(1)=start(1)-xoffset;
        newstart(3)=start(3)-xoffset;
    case 3
        newstart(1)=start(1)-xoffset;
        newstart(3)=start(3)-xoffset;
        newstart(5)=start(5)-xoffset;
    case 4
        newstart(1)=start(1)-xoffset;
        newstart(3)=start(3)-xoffset;
        newstart(5)=start(5)-xoffset;
        newstart(7)=start(7)-xoffset;
     case 5
        newstart(1)=start(1)-xoffset;
        newstart(3)=start(3)-xoffset;
        newstart(5)=start(5)-xoffset;
        newstart(7)=start(7)-xoffset;        
        newstart(9)=start(9)-xoffset;
    otherwise       
end % switch NumPeaks

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.001,'TypicalX',center,'Display','off' );
LowestError=1000;
FitParameters=zeros(NumPeaks.*2); 
BestStart=zeros(NumPeaks.*2); 
height=zeros(NumPeaks); 
bestmodel=zeros(length(yy));
for k=1:NumTrials,
    % disp(k) % Progress indicator during repeat trial fits
    switch NumPeaks
       case 1
          newstart=[newstart(1)*(1+randn/50) newstart(2)*(1+randn/10)]; 
       case 2
          newstart=[newstart(1)*(1+randn/50) newstart(2)*(1+randn/10) newstart(3)*(1+randn/50) newstart(4)*(1+randn/10)]; 
       case 3
          newstart=[newstart(1)*(1+randn/50) newstart(2)*(1+randn/10) newstart(3)*(1+randn/50) newstart(4)*(1+randn/10) newstart(5)*(1+randn/50) newstart(6)*(1+randn/10)]; 
       case 4
          newstart=[newstart(1)*(1+randn/50) newstart(2)*(1+randn/10) newstart(3)*(1+randn/50) newstart(4)*(1+randn/10) newstart(5)*(1+randn/50) newstart(6)*(1+randn/10)  newstart(7)*(1+randn/50) newstart(8)*(1+randn/10)]; 
       case 5
          newstart=[newstart(1)*(1+randn/50) newstart(2)*(1+randn/10) newstart(3)*(1+randn/50) newstart(4)*(1+randn/10) newstart(5)*(1+randn/50) newstart(6)*(1+randn/10)  newstart(7)*(1+randn/50) newstart(8)*(1+randn/10)  newstart(9)*(1+randn/50) newstart(10)*(1+randn/10)]; 
      otherwise       
    end % switch NumPeaks
  switch peakshape
    case 1
        TrialParameters=fminsearch(@fitgaussian,newstart,options,xx,yy);
        ShapeString='Gaussian';
    case 2
        TrialParameters=fminsearch(@fitlorentzian,newstart,options,xx,yy);
        ShapeString='Lorentzian';    
    case 3
        TrialParameters=fminsearch(@fitlogistic,newstart,options,xx,yy);
        ShapeString='Logistic';
    case 4
        TrialParameters=fminsearch(@fitpearson,newstart,options,xx,yy,extra);
        ShapeString='Pearson';
    case 5
        TrialParameters=fminsearch(@fitexpgaussian,newstart,options,xx,yy,-extra);
        ShapeString='ExpGaussian';
    otherwise
end % switch peakshape

% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
   switch peakshape
    case 1
        A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 2
        A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 3
        A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
    case 4
        A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
    case 5
        A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
    otherwise
  end % switch
end % for
% Multiplies each row by the corresponding amplitude and adds them up
model=PEAKHEIGHTS'*A;

% Compare trial model to data segment and computer fit error
  MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(PEAKHEIGHTS)>0,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=PEAKHEIGHTS; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min
  end % if MeanFitError
end % if k

% Construct model from best-fit parameters
AA=zeros(NumPeaks,100);
xxx=linspace(min(xx),max(xx));
for m=1:NumPeaks,
   switch peakshape
    case 1
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra)';
    otherwise
  end % switch
end % for

% Multiplies each row by the corresponding amplitude and adds them up
mmodel=height'*AA;

% Top half of the figure shows original signal and the fitted model.
figure(2);subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
hold on
for m=1:NumPeaks,
    figure(2);plot(xxx+xoffset,height(m)*AA(m,:),'g')  % Plot the individual component peaks in green lines
    area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
end

% Mark starting peak positions with vertical dashed lines
for marker=1:NumPeaks,
    markx=BestStart((2*marker)-1);
    figure(2);subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
end % for
figure(2);plot(xxx+xoffset,mmodel,'r');  % Plot the total model (sum of component peaks) in red lines
hold off;
figure(2);axis([min(xx)+xoffset max(xx)+xoffset min(yy) max(yy)]);
title(['Best of ' num2str(NumTrials) ' trial fits.' ])
xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Error = ' num2str(LowestError) '     Extra = ' num2str(extra) ] )

% Bottom half of the figure shows the residuals and displays RMS error
% between original signal and model
residual=yy-bestmodel;
figure(2);subplot(2,1,2);plot(xx+xoffset,residual,'b.')
figure(2);axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
xlabel('Residual Plot')

% Put results into a matrix, one row for each peak, showing peak index number,
% position, amplitude, and width.
for m=1:NumPeaks,
    if m==1,
       FitResults=[FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m) LowestError];
    else
       FitResults=[FitResults ; [FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m) LowestError]];
    end % if
end % for

% Display Fit Results on upper graph
subplot(2,1,1);
residualaxis=axis;
figure(2);text(residualaxis(1), 0.8*residualaxis(4) ,num2str(FitResults));
% ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker)+ xoffset;
      start=[start markx n/5];
  end % for marker
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band spectrum.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6006.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for single lorentzian, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function err = fitlogistic(lambda,t,y)
%	Fitting function for logistic, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlogistic assumes a logistic function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991 
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function err = fitlognormal(lambda,t,y)
%	Fitting function for lognormal, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlognormal assumes a lognormal function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lognormal(t,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991  
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band spectrum.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function. 
% g = pearson7(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
%  T. C. O'Haver, 1990  
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened Gaussian band spectrum.
%  T. C. O'Haver, October 23, 2006.
global PEAKHEIGHTS
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2006
g = exp(-((x-pos)./(0.6006.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) convolutes y by an exponential decay of time constant t
% by multiplying Fourier transforms and inverse transforming the result.
fy=fft(y);
a=exp(-[1:length(y)]./t);
fa=fft(a);
fy1=fy.*fa';           
yb=real(ifft(fy1))./sum(a);  
% ----------------------------------------------------------------------
function err = fitNewPeakShape(lambda,x,y)
% Replace "NewPeakShape" with the peak function of your choice.
%  T. C. O'Haver, 2006   
global PEAKHEIGHTS
A = zeros(length(x),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = NewPeakShape(x,lambda(2*j-1),lambda(2*j))';
end
PEAKHEIGHTS = A\y';
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function g = NewPeakShape(x,a,b)
% Replace this with the peak function of your choice.
% a and b are the two non-linear parameters.
g=sin(a.*x+b);
% ----------------------------------------------------------------------
% Note: Add this to your Matlab path if you don't already have it.
%%  function [index,closestval]=val2ind(x,val)
%%  Returns the index and the value of the element of vector x that is closest to val
%%  If more than one element is equally close, returns vectors of indicies and values
%%  Tom O'Haver (toh@umd.edu) October 2006
%%  Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
%%  [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
% dif=abs(x-val);
% index=find((dif-min(dif))==0);
% closestval=x(index);