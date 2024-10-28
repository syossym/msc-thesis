function [ML,Chi2,MLRT,SimML,P_SimML]=max_likelihood_exp(Lambda, Events, DoF, Nsim, Psig, NFactor,AddOffset,InterpMethod);
%----------------------------------------------------------------------------
% max_likelihood_exp function   Assuming exponential distribution in the range
%                            [0:Inf] and list
%                            of "events". Calculates the maximum likelihood
%                            for the events.
%                            Calculate the ML ratio test.
%                            Using Monte-Carlo simulation of the parent
%                            numerical distribution, calculate the ML
%                            probability distribution.
% Input  : - Lambda of the exponential distribution
%            (e.g., Lambda.*exp(-Lambda.*X)
%          - List of event [X].
%          - Number of Degree's of Freedom. Default is 1
%            (Free parameters of the model).
%          - Number of Monte-Carlo simulations. Default is 1000.
%          - Column vector of Probabilities for the ML ratio-test.
%            Return the MLRT for each probability.
%            Default is [0.6827; 0.9545; 0.9973].
%            Usefull when ML is calculated as function of a free parameters.
%          - Optional value: If (ParNumDist) is single-column, then
%            generate ParNumDist so that each bin will have (Nfactor) points.
%            Default is 100.
%          - Optional offset to add to the objects and to the
%            Monte-Carlo objects.
%            In each MC simulation, add an offset to the simulated events
%            before calculating the ML.
%            Default is 0 - no offset.
%          - Interpolation method in calculating the probability from
%            the numerical function. Default is 'linear'.
% Output : - log Likelihood.
%          - Chi2 per DoF.
%          - ML Ratio-test.
%          - Vector of simulated MLs.
%          - Given the ML and Simulated ML, return the probability to
%            get the Events ML from the Parent Simulated ML distribution:
%            P(SimML>ML).
% Tested : Matlab 5.3
%     By : Eran O. Ofek          September 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------------
NsimDef    = 1000;   % default number of simulations
PsigDef    = [0.6827; 0.9545; 0.9973];
NFactorDef = 100;
InterpMethodDef = 'linear';
ProbType   = 'd';
NewSeed    = 'n';

if (nargin==2),
   DoF     = 1;
   Nsim    = NsimDef;
   Psig    = PsigDef;
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==3),
   Nsim    = NsimDef;
   Psig    = PsigDef;
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==4),
   Psig    = PsigDef;
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==5),
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==6),
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==7),
   InterpMethod = InterpMethodDef;
elseif (nargin==8),
   % do nothing
else
   error('Illegal number of input arguments');
end


Ne      = length(Events);

% sorted events vector
SEvents = sort(Events);


% calculate the probability density per event
SEvents = SEvents + AddOffset;
PE = Lambda.*exp(-Lambda.*SEvents);   % integral [0:Inf] is normalize to 1

%K       = find(isnan(PE)==1);
%PE(K)   = MinProb;
%K       = find(PE==0);
%PE(K)   = MinProb;
% calculate the ML
ML      = sum(log(PE));


% ML ratio test
Chi2 = chi2inv(Psig,DoF);
MLRT = -0.5.*Chi2;


% Monte-Carlo simulation
if (nargout>3),
   SimML = zeros(Nsim,1);
   for I=1:1:Nsim,

      %SimEvents  = randgen([ParNumDist(:,1), ParNumDist(:,2)./max(ParNumDist(:,2))],Ne,ProbType,InterpMethod,NewSeed);
      %SimEvents  = randgen([ParNumDist(:,1), ParNumDist(:,2)],Ne,ProbType,InterpMethod,NewSeed);

      SimEvents = exprnd((1./Lambda),Ne,1);

      SimSEvents = sort(SimEvents);

      % optional offset
      SimSEvents = SimSEvents + AddOffset;

      % calculate the probability density per simulated event
      %PSE      = interp1(ParNumDist(:,1), ParNumDist(:,2), SimSEvents,'linear');
      PSE = Lambda.*exp(-Lambda.*SimSEvents);


      %K        = find(isnan(PSE)==1);
      %PSE(K)   = MinProb;
      %K        = find(PSE==0);
      %PSE(K)   = MinProb;
      SimML(I) = sum(log(PSE));

   end

   % The probability to get the ML from the Monte-Carlo.
   P_SimML = length(find(SimML>ML))./Nsim;

end

