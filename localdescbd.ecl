
:- local store(recsteps).
:- setval(total,0).

% benefit(K,T,N,Bd,Max,Benefit)
% K starting cost ; T target cost; N neighbour count ; 
% Bd bound on cost difference between neighbours ; Max worst cost
% Benefit: Difference between blind search and local blind search
% 
% Calculate the number of steps for local descent 
% starting at Steps_{T+1} = steps(T+1, T,N,Bd,Max)
% storing the value of Steps_{T+1} which is used when
% computing Steps_{i} for all i>T+1,
% then calculating Steps_{i} for i from T+1 up to K
% Finally the Benefit is Steps_{k} - blind(T)
benefit(K,T,N,Bd,Max,Benefit) :-
        NT is T+1,
        (for(I,NT,K), param(T,N,Bd,Max)
        do 
           Steps is steps(I,T,N,Bd,Max),
           store_set(recsteps,I,Steps)
        ),
        store_get(recsteps,K,AllSteps),
        Benefit is blind(Max,T) - AllSteps.

% Call `benefit' for each value of K from T+1 up, 
% until the returned benefit (blind search steps - local search steps)
% is negative
bentestk(T,N,B,Max,BenK) :-
        int(T+1, Max // 2, K),
        Ben is benefit(K,T,N,B,Max),
        Ben<0-1,
        BenK is K-1.


% succN returns the probability SuccP there is a better neighbour
% and the expected number of step SuccCt needed to find it
steps(K,T,N,B,Max,Steps) :-
        succN(K,N,B,Max,SuccP,SuccCt), 
        nextsteps(K,T,B,Max,NSteps),
        blind(Max,T, Blind),
        Steps is SuccP*(SuccCt+NSteps)
               + (1-SuccP)*(N+ Blind). 

% The probability of improving, SuccP, 
% is the probability of failing (Pf=1-Pn)
% N times: Succp = 1-Pf^N
% Count steps returns the expected number of steps to succeed
succN(K,N,B,Max,SuccP,Ct) :-
        pn(K,B,Max,Pn),
        Pf is 1-Pn,
        SuccP is 1-Pf^N,
        count_steps(Pn,N,Ct).

% The probability of improving is the probability ct 
% for each possible cost better than K
% divided by the probability cost of all possible neighbours
pn(K,B,Max,Prob) :-
        KP is min(Max,K+B),
        KN is max(0,K-B),
        KN1 is K-1,
        sumnci(Max,KN,KN1,Better), 
        sumnci(Max,K,KP,Worse),
        Prob is Better / (Better+Worse).  

% Sum the probability cost Ct (= nci(Max,K)) of all costs in the range S..E
sumnci(Max,S,E,Prob) :-
        (fromto(0,This,Next,Prob),for(I,S,E),param(Max)
        do  Next is This+nci(Max,I) ).


% If the starting cost is the target plus 1, then 
% an improving neighbour must reach the target, 
% so no more steps will be needed
nextsteps(K,T,_,_,SBen) :-
        K =:= T+1, !, SBen=0.
% Otherwise, count all the improving neighbours: TotalCt
% Now for each improving cost above the target cost
% count the points with this cost (CtJ)
% The probability a neighbour has this cost is CtJ/TotalCt
% (by uniform Nweight below Bd, and normal neighbourhood)
% The remaining steps from starting cost J is JSteps
% which is incurred with probability PJ
% Add all these weighted probabilities to get the expected number of steps.
nextsteps(K,T,Bd,Max,Steps) :-
        KN is K-1,
        TP is T+1,
        S is max(K-Bd,0),
        sumnci(Max,S,KN,TotalCt),
        (for(J,TP,KN),fromto(0,This,Next,Steps),
         param(Max,TotalCt)
        do
            nci(Max,J,CtJ),
            PJ is CtJ/TotalCt,
            store_get(recsteps,J,JSteps),
            Next is This + PJ * JSteps 
        ).

% Sum up, for I between 0 and Max, 
% the number of ways of choosing I from Max
total(Max,Total) :-
    (for(I,0,Max),fromto(0,T,N,Total),param(Max) 
        do N is T+nci(Max,I)
    ). 
% The inverse of the probability, Prob, 
% of choosing a point with cost This or better
% = 1/Prob
% Calculated by summing up, for I from 0 to This, 
%  the number of ways of choosing I from Max divided by Total 
%  = Prob,  and taking the inverse 1/Prob
blind(Max,This,Blind) :-
       getval(total,T), 
       (T=:=0 -> total(Max,RMax), setval(total,RMax) ; RMax=T),
       (for(I,0,This),fromto(0,T,N,RThis),param(Max) 
        do true, %writeln(I),
           N is T+nci(Max,I)),
       Blind is RMax/RThis.

% Mathematical expression 
% equivalent to "rsucc" below
count_steps(Pn,N,Ct) :-
        Pf is 1-Pn,
        Ct is (1/Pn) - Pf^N * ((1/Pn) + N). %Same result as rsucc

% Just used to check on the previous 
rsucc(K,N,B,Max,Ct) :-
        pn(K,B,Max,Pn),
        recs(1,N,Pn,Ct).

% The step count for the Mth trial is
% M times the probability of success on the Mth trial
% which is the probability of success times 
% the prob of failure to the (M-1)th power

% If M is the limit, then this is the final count, Ct
% Otherwise keep counting the expected number of steps 
% for more than M trials, and add the result to Ct
recs(M,N,Pn,Ct) :-
        Ct1 is M*Pn*((1-Pn)^(M-1)),
        (N=M -> Ct is Ct1 ;
                M2 is M+1, recs(M2,N,Pn,Ct2),
                Ct is Ct1+Ct2
        ).


%%%%%%%% Unused  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


enimp(K,B,Max, Enimp) :-
        KP is min(Max,K+B),
        KN is max(0,K-B),
        KN1 is K-1,
        sumnci(Max,KN,KP,NeighCt),
        (for(J,KN,KN1), 
         fromto(0,This,Next,Enimp),param(Max,K,NeighCt)
        do
            I is K-J,
            nci(Max,J,CtJ),
            Next is This + (CtJ*I/NeighCt)
        ).

eimp(K,Max,Eimp) :-
       KN1 is K-1,
       total(Max,Total), 
       (for(J,0,KN1), 
         fromto(0,This,Next,Eimp),param(Max,K,Total)
        do
            I is K-J,
            nci(Max,J,CtJ),
            Next is This + (CtJ*I/Total)
       ).       

%test
eimptestk(T,B,Max,ImpK) :-
        int(T+1,Max // 2,K),
        enimp(K,B,Max, Enimp),
        eimp(K,Max,Eimp),
        Eimp>Enimp+1,
        ImpK is K-1.

binit(K,T,N,Bd,Max) :-
        K=100, T=40, N=50, Bd=20, Max=200,
        total(Max,Total), setval(total,Total).

disp(L2) :-
List = [1-23,3-37,5-28,7-23,9-22,11-23,13-24,15-26,17-28,19-30],
    (foreach(Bd-K,List),foreach(EnI-EI,L2) do
    enimp(K,Bd,200,EnI), eimp(K,200,EI)).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Useful generic predicates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Error checking
check(Goal,_Error) :-
        call(Goal), !.
check(_,Error) :- writeln(Error), abort.

% Return an integer I between Min and Max
int(Min,Max,I) :-
        X is eval(Min), Y is eval(Max),
        int1(X,Y,I).
int1(M,_,M).
int1(Min,Max,X) :- Next is Min+1, Next=<Max, int1(Next,Max,X).

% Factorial: fac(N) = R
fac(N,R) :-
        N>1, !,
        M is N-1,
        fac(M,RM),
        R is N*RM.
fac(N,1) :- N>=0, !.
fac(Z,_) :- writeln(error(fac(Z,_))), abort.

% fac(N)/fac(I)
fac(N,I,R) :- N=I, !, R=1.
fac(N,I,R) :- N>I, !, M is N-1, fac(M,I,RM), R is N*RM.
fac(N,I,_) :- writeln(error(fac(N,I,_))), abort.

%Exponential
exp(X,Y,R) :- Y>0, !, Y2 is Y-1, exp(X,Y2,R2), R is X * R2.
exp(_X,0,1).

% Number of ways of choosing I/4 from N/4 
nci(N,I,R) :-
        N4 is (N div 4)+1,
        I4 is I div 4,
        M is I mod 4,
        (M =:= 0 -> nchoosei(N4,I4,R) ; ncim(N4,I4,M,R)).

% Number of ways of choosing I from N  (nCi) 
nchoosei(N,I,R) :-
        NI is N-I,
        R1 is fac(N,I), R2 is fac(NI), 
        R is R1 div R2.

% M/4th of the distance along the line nCi to nC(i+1)
ncim(N,I,M,R) :-
        nchoosei(N,I,RL),
        I1 is I+1,
        nchoosei(N,I1,RH),
        R is RL + (RH-RL) * M / 4.


