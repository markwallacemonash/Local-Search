% 2-SAT
% 100 clauses each with two variables
% 50 distinct variables
% Each variable occurs in 4 clauses
% Minimise violated clauses

% Assume X/100 clauses are false
% We find a neighbour by flipping a variable
% Flipping a variable in an unsatisfied clause always makes it
% satisfied.
% Flipping a variable in a satisfied 2-clause makes it unsatisfiable
% with probability 1/3

% The change after flipping is 0 iff
% (a) all four clauses in which the flipped variable appears are satisfied
% and  the clause stays satisfied with probability 2/3
% or (b) one of the clauses is satisfied and flipping the variable makes the
% clause false, and one of the clauses is unsatisfied and flipping the
% variable makes it true (necessarily), and two of the
% clauses are satisfied but flipping the variable does not make it
% false (there are 12 ways of this occurring in four clauses)
% or (c) two of the clauses are unsatisfied (and necessarily become
% true) and two are satisfied and flipping the variable makes
% both false (there are 6 ways of this occurring in four clauses)
?- lib(lists).


% pn(X,Y,P): 
% if X is the number of false clauses, 
% and Y is the number of false clauses after flipping one variable
% P is the probability of the transition X-Y by flipping one variable
pn(_X,Y,Pn) :- (Y<0;Y>100), !, Pn is 0.
pn(XE,YE,Pn) :- Diff is eval(XE)-eval(YE), abs(Diff)>4, !, Pn=0.
% If flipping a variable makes a true clause false call it +1
% If flipping a variable leaves a true clause true call it 0
% If flipping a variable makes a false clause true call it -1
% For each length 4 sequence that adds to the change (XE-YE) compute its
% probability.
% e.g. if XE-YE = 1 then the possible sequences are:
% <0,0,0,1>, <0,0,1,0>, <0,1,0,0>, <1,0,0,0>, <0,1,1,-1>, <1,0,1,-1>, 
% <1,1,0,-1>, <0,1,-1,1>, <0,-1,1,1>, <-1,0,1,1>, <1,0,-1,1>, 
% <1,-1,0,1>, <-1,1,0,1>, <1,1,-1,0>, <1,-1,1,0>, <-1,1,1,0>
pn(XE,YE,Pn) :- X is eval(XE), Y is eval(YE), Diff is X-Y,
        genall(Diff,ABCSet),
        ( foreach([A,B,C],ABCSet),foreach(PSum,ProbList), param(X)
        do
            findall(P,pcalc(100,X,0,A,B,C,P),PList), sum(PList,PSum)
        ),
        sum(ProbList,Pn).

% Find all combinations of [+1,0,-1] of size 4 whose members add to Diff
% eg genall(1,List) => List= [[1,3,0],[2,1,1]]
genall(Diff,ABCSet) :-
        genall(Diff,0,0,0,ABCList),
        setof(X,member(X,ABCList),ABCSet).
genall(_Diff,NA,NB,NC,[[NA,NB,NC]]) :- 
        NA+NB+NC =:= 4, !.
genall(Diff,NA,NB,NC,ABCList) :-
        NextC is NC+1, 
        ((NA-NextC) + (4 - (NA+NB+NextC))  >= Diff ->
                genall(Diff,NA,NB,NextC,List1) ;
                List1=[]
        ),
        NextA is NA+1, 
        ((NextA - NC) - (4- (NextA+NB+NC)) =< Diff ->
                genall(Diff,NextA,NB,NC,List2) ;
                List2 = []
        ),
        NextB is NB+1,
        ((NA-NC) + (4 - (NA+NextB+NC))  >= Diff,
         (NA - NC) - (4- (NA+NextB+NC)) =< Diff ->
                genall(Diff,NA,NextB,NC,List3) ;
                List3 = []
        ),
        append([List1,List2,List3],ABCList).
        

% pcalc(4,F,V,A,B,C,Prob) finds one possible sequence of 4 clauses in which
% the flipped variable occurs, and multiplies probabilities of each clause
% pcalc is recursive, and calls itself with reduced values of its
% arguments until the first argument is 0.
% pcalc is also non-deterministic, and finds more sequences, and their
% probabilities, on backtracking.

% Ct is the total number of remaining clauses in the problem
% F is the number of false clauses remaining
% V is the number of clauses treated so far
% A is the number of -1 clauses: false clauses (where a variable flip always makes it true)
% B is the number of  0 clauses: true clauses where a variable flip leaves it true
% C is the number of +1 clauses: true clauses where a variable flip makes it false
% Prob is the probability of these clauses being encountered in a certain order
% The "findall" above ensures the probability of each possible order is calculated

% If 4 clauses already treated, then stop 
% without changing the probability (hence Prob=1)
pcalc(_,_,4,_,_,_,Prob) :- !, Prob=1.
% Otherwise choose an A (-1), a B (0) or a C (+1) for the next clause
% calculate its probability,
% call pcalc again for the remaining clauses returning their probability
% multiply the probabilities
pcalc(Ct,F,V,A,B,C,Prob) :-
     % Variable flip makes a true clause false
     (C>0, P1 is ((Ct-F)/Ct) * (1/3), F2=F, C2 is C-1, B2=B,A2=A ;
     % Variable flip makes a true clause remain true
      B>0, P1 is ((Ct-F)/Ct)*(2/3),   F2=F,      B2 is B-1, A2=A,C2=C ;
     % Variable flip makes a false clause true
      A>0, P1 is (F/Ct),              F2 is F-1, A2 is A-1, C2=C,B2=B
     ),
     P1>0,  % If the probability P1=0, don't compute further
            % probabilities for this sequence (i.e. fail this sequence
            % immediately)
     V2 is V+1,  % number of clauses treated increased by 1
     Ct2 is Ct-1, % Number of unseen clauses remaining in the system
                  % decreased by 1
     pcalc(Ct2,F2,V2,A2,B2,C2,P2),
     Prob is P1*P2.


        
% The probability a 2-clause is true is 3/4 (each clause has two
% distinct variables, and they can be (TT, TX or XT).  The probability
% it is false is 1/4 (FF).
% The probability that X out of 100 clauses are true is
% (3/4)^X*(1/4)^(100-X) * 100!/(X!*(100-X)!)

% p(X,P) 
% P is the probability that precisely X clauses (out of 100) are false
p(X,P) :- (X<0;X>100), !, P=0.
p(X,P) :- NX is X, PX is 100-NX,   
        P is 0.75^PX * 0.25^NX * (fac(100)/(fac(NX)*fac(PX))).

% This clause writes to a file (called sat2p) for each number of false
% clauses (I) from 1 to 99, its probability (P)
genp :-
        open('sat2p.txt',write,S),
        (for(I,1,99),param(S) 
        do
            p(I,P),approx(P,1000000,AP),
            write(S,I),write(S,'     '),write(S,AP),nl(S)
        ),
        close(S).

approx(W,M,AW) :-
            WW is eval(W),
            I is WW*M, II is fix(round(I)),
            AW is II/M.

% r(K,Delta,R) holds if r(K,Delta) = R 
% This is computed as:  (pn(X,X+N)+pn(X,X-N))/(p(X+N)+p(X-N))
r(X,0,R) :- !, pn(X,X,Pn),p(X,P), R is Pn/P.
r(X,N,R) :- pn(X,X+N,PPn), 
            pn(X,X-N,NPn),
            p(X+N,PM1), 
            p(X-N,PP1), 
            (PM1+PP1 =:= 0 -> true ; R is (PPn+NPn)/(PM1+PP1)).      

% posr(K,Delta) = pn(K,K-Delta) - r(K,Delta) * p(K-Delta)
posr(K,D,P) :- pn(K,K-D,P1), p(K-D,P2), r(K,D,R),
           (var(R) -> P=0 ; P is P1-(P2*R)).

% The average value of R in the range 1..(100-X) is as coded here,
% because if I>4 R(X,I)=0
avr(X,Av) :- 
        r(X,1,R1), r(X,2,R2), r(X,3,R3),r(X,4,R4),
        (var(R1) -> Av=1 ;
                   desc([R1,R2,R3,R4]), Av is (R1+R2+R3+R4)/X
        ).
% If R(X,D) is undetermined (since p(X,D) is 0), 
% set it so that R is descending.
desc([_]) :- !.
desc([H1,H2|T]) :- var(H2), !, H2=H1, desc([H2|T]).
desc([H1,H2|T]) :- desc([H2|T]).

% The modal number of true clauses is the number, N, with the highest
% probability p(N)
modal(Mode) :- ( for(J,1,99),fromto((0,0),(T,TJ),(N,NJ),(Max,Mode))
                do
                    P is p(J), (P>T -> NJ=J,N=P ; (N,NJ)=(T,TJ))
                ),
                writeln(maxprob(Max)).

% Since no flip can change the number of unsatisfied clauses by more
% than 4, the "good enough" level below which the NSC properties must
% hold is the modal fitness - 4.
% The flip neighbourhood is unskewed if posr(K,Delta) >= 0 for all
% values of K below Mode-4, and Delta=<4
unskewed :-
     modal(Mode),Max is Mode-4,
     (multifor([C,D],[1,1],[Max,4]) do
         posr(C,D) >= 0
%          pn(C,C-D)/p(C-D) >= pn(C,C+D)/p(C+D)
     ).

% The largest number in the range 1..Mode-4 for which local
% improvement is beneficial is the first number N for which 
% either r(N,0)=<1 or the average avr(N)>=1
min_ben(Max) :- %modal(Mode), 
              First is 99, %Mode-4
              rec_ben(First,0,Max).

rec_ben(Min,Min,_) :- writeln("Local search not beneficial"), !.
%% rec_ben(X,Min,Max) :-
%%         r(X,0,R0),avr(X,Av),
%%         (R0=<1 -> Max=X, writeln("R0=<1") ;
%%          Av>=1 -> Max=X, writeln("Average R >= 1") ;
%%          true  -> NX is X-1, rec_ben(NX,Min,Max)
%%         ).
rec_ben(X,Min,Max) :-
        PX is X-1,
        (for(I,0,PX),foreach(P,PList) do p(I,P)),
        (for(I,0,PX),foreach(Pn,PnList),param(X) do pn(X,I,Pn)),
        (sum(PnList)>sum(PList) -> Max=X, 
                Ratio is sum(PnList)/sum(PList), writeln(ratio(Ratio))
        ; rec_ben(PX,Min,Max)
        ).

% agt(K,P)
% P is the difference in probability between blind search returning a
% point with cost between K+1 and 2K, and neighbourhood search starting 
% at a point with cost K returning a neighbour with cost in this range
agt(K,P) :-
        pgt(K,P1),pngt(K,P2),
        P is P1-P2.

% agtgt(K,P)
% P is the difference in probability between blind search returning a
% point with cost worse than 2K+1, and neighbourhood search starting 
% at a point with cost K returning a neighbour with cost in this range
agtgt(K,P) :-
        pgtgt(K,P1), pngtgt(K,P2),
        P is P1-P2.

% plt(K,Prob)
% Prob is the probability of blind search returning a point with cost
% better than K
plt(K,Prob) :-
        (for(I,0,K-1),foreach(P,PList) do p(I,P)), Prob is sum(PList).
% pnlt(K,Prob)
% Prob is the probability of neighbourhood search starting at a point
% with cost K, returning a neighbour with cost better than K
pnlt(K,Prob) :-
        (for(I,0,K-1),foreach(Pn,PnList),param(K) do pn(K,I,Pn)), 
        Prob is sum(PnList).

pgt(K,Prob) :-
        (for(I,K+1,2*K),foreach(P,PList) do p(I,P)), Prob is sum(PList).        
pngt(K,Prob) :-
        (for(I,K+1,2*K),foreach(Pn,PnList),param(K) do pn(K,I,Pn)), 
        Prob is sum(PnList).

pgtgt(K,Prob) :-
        (for(I,(2*K)+1,100),foreach(P,PList) do p(I,P)), 
        Prob is sum(PList).
pngtgt(K,Prob) :-
        (for(I,(2*K)+1,100),foreach(Pn,PnList),param(K) do pn(K,I,Pn)), 
        Prob is sum(PnList). 

%% % This should write 1.0-1.0 for every value of K
%% test :- for(K,1,25) do
%%             (P is p(K)+pgt(K)+pgtgt(K)+plt(K),
%%              Pn is pn(K,K)+pngt(K)+pngtgt(K)+pnlt(K),
%%              write(K),write(': '),
%%              writeln(P-Pn)).

% For the paper table 2satbeneficial
t(K,P) :- P is pn(K,K)-p(K).
a(K,P) :- P is agt(K)+agtgt(K).
table(satbeneficial) :-
        for(K,2,26,3) do t(K,P1),a(K,P2),write(K),write(': '),writeln(tk:P1-ak:P2).


fac(N,Fac) :-
      for(J,1,N),fromto(1,T,N,Fac) do N is T*J.

max_nsf(Max) :- 
        rec_nsf(2,100,Max).

rec_nsf(Mode,Mode,_) :- !, writeln('NSF property holds everywhere').
rec_nsf(X,Mode,Max) :-
        r(X,1,R1), 
        r(X,2,R2), R2=<R1,
        r(X,3,R3), R3=<R2,
        r(X,4,R4), R4=<R3, !,
        Y is X+1, 
        rec_nsf(Y,Mode,Max).
rec_nsf(X,_,Y) :- Y is X-1.

intin(X,X,_).
intin(K,Min,Max) :- Min<Max, Next is Min+1, intin(K,Next,Max).
