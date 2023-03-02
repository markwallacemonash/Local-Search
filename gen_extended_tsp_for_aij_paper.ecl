/*
?- tsp(Cities,ValueCt,Range,   Total).
*/

% Input
% Cities: The number of cities
% ValueCt: The number of distinct edge cost values (the values are
% 1..ValueCt)
% Count: the number of problem instances
% Range: The amount of discretisation of the objective

% Description:
% The procedure creates a list of (edges and their) costs.
% Values in 1..ValueCt are randomly assigned to the edges.
% Finally for each total cost C between Cities and ValueCt*Cities, 
% the system finds all sets of edges whose total cost is C, and counts them.

% Output
% The number of solutions at each cost
% The probability p(K) a solution will have a given cost K
% The probability pn(K,K-d) that a neighbour of a solution with cost K
% has cost K-d
% The r(K,d) which is the ratio pn(K,K-d)/p(K-d)
% For K in the top quarter of actual costs, the average (over d) of r(K,d)

% I ran beneficial_tsp(10,25,100) on Feb 25 2023 and it succeeded
beneficial_tsp(Cities,ValueCt,Count) :-
    create_arrays(Cities,ValueCt,Size),
    seed(1),  % Make the solution repeatable
    setval(exct,0),
    repeat,
        dist_matrix(Cities,ValueCt,Dist),
        all(route(Cities,Dist)),
        count_all_solutions(Size),
        min_max_mode(Size),
        record_p_pn(Size),
    beneficial,
    incval(exct),
    getval(exct,Count).

check(Goal) :- Goal, !.
check(Goal) :- getval(exct,This), writeln(This-Goal-failed), abort.

beneficial :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), param(Total,Min)
    do 
        getval(neighct(K,0),KNeigh),
        calc_r(K,Min,Total,KNeigh,R0),
        (nonvar(R0), R0>=1 -> true ;  % don't check further if r(k,k-k_opt) >=1
           check(avr_gt_1(K,Total,Min)), % check average NWeight avr(k) >= 1
           check(pnlow_gt_avrxp(K,Total,Min))  % Check pn<(k) > avr(k) * p<(k)
        )
    ).

calc_r(K,Min,Total,KNeigh,R0) :-
         D is K-Min,
         ct_points_at_dist(K,D,Ct),
         ct_neigh_at_dist(K,D,NCt),
         P is Ct/Total,
        ( KNeigh>0 -> 
           PN is NCt/KNeigh,
           R0 is PN/P 
          ; R0=1
        ).

avr_gt_1(K,Total,Min) :-
    calc_r_list(K,Total,Min,RList),
    AvR is sum(RList)/length(RList), 
    AvR>=1.

pnlow_gt_avrxp(K,Total,Min) :-
        calc_r_list(K,Total,Min,RList),
        AvR is sum(RList)/length(RList),
        KNeigh is getval(neighct(K,0)),
        MaxDiff is (K-Min),
        (for(D,1,MaxDiff), foreach(P,PList), foreach(Pn,PnList),
         param(Total,K,KNeigh)
        do 
            LowFit is K-D,
            getval(dynct(LowFit),Ct),
            P is Ct/Total,
            NCt is getval(neighct(K,LowFit)),
            Pn is NCt/KNeigh            
        ),
        PnLow is sum(PnList),
        PAvR is sum(PList) * AvR,
        PnLow >= PAvR.              

pbr_gt_avrxp(K,Total,Min) :-
        calc_r_list(K,Total,Min,RList),
        AvR is sum(RList)/length(RList), 
        MaxDiff is (K-Min),
        (for(D,1,MaxDiff), foreach(R,RList), foreach(PR,PRList), foreach(P,PList),
         param(Total,K)
        do 
            LowFit is K-D,
            getval(dynct(LowFit),Ct),
            P is Ct/Total,
            PR is P*R
        ),
        PBR is sum(PRList),
        PAvR is sum(PList) * AvR,
        PBR >= PAvR.      

pnlow_gt_pbrlow(K,Total,Min) :-
        calc_r_list(K,Total,Min,RList),
        KNeigh is getval(neighct(K,0)),
        MaxDiff is (K-Min),
        (for(D,1,MaxDiff), foreach(R,RList), foreach(PR,PRList), foreach(Pn,PnList),
         param(Total,K,KNeigh)
        do 
            LowFit is K-D,
            getval(dynct(LowFit),Ct),
            P is Ct/Total,
            PR is P*R,
            NCt is getval(neighct(K,LowFit)),
            Pn is NCt/KNeigh
        ),
        PBR is sum(PRList),
        PnLow is sum(PnList),
        PnLow >= PBR.     

% I ran tsp(10,25,1)
tsp(Cities,ValueCt,Count) :-
    create_arrays(Cities,ValueCt,Size),
    seed(1),  % Make the solution repeatable
    setval(exct,0),
    repeat,
    dist_matrix(Cities,ValueCt,Dist),
    all(route(Cities,Dist)),
    incval(exct),
    getval(exct,Count), !,
    count_all_solutions(Size),
    min_max_mode(Size),
    record_p_pn(Size),
    show_p(Cities,ValueCt,Count),
    show_p_pn(Cities,ValueCt,Count),
    show_r_kge(Cities,ValueCt,Count),
    show_pn_kge(Cities,ValueCt,Count),
    show_p_r_kge(Cities,ValueCt,Count),
    show_r_av(Cities,ValueCt,Count),
    show_low_pn(Cities,ValueCt,Count),
    show_low_pbr(Cities,ValueCt,Count),
    show_low_pavr(Cities,ValueCt,Count),
    show_high_pbr(Cities,ValueCt,Count),
    show_high_pavr(Cities,ValueCt,Count).


% Assign costs to only the edges [I,J] with I<J.
% For edge [I,I] the cost is 0
% Assuming a symmetric tsp, the costs are the same in both directions
dist_matrix(Cities,ValueCt,Dist) :-
        dim(Dist,[Cities,Cities]),
        (multifor([I,J],[1,1],[Cities,Cities]), param(Dist,ValueCt)
        do
             (I=J -> 0 is Dist[I,J] ;
              I<J -> 
                     E is 1+nrand(ValueCt), %(random mod ValueCt)+1,
                     E is Dist[I,J],
                     E is Dist[J,I] ;
              I>J -> true
             )
        ).

% The array dynct will hold the number of solutions for each cost
% The array neighct will hold the number of neighbours linking each pair
% of costs
% Thus neighct(I,J) records the number of pairs of neighbouring
% solutions S1 and S2 where S1 has cost I and S2 has cost J  
create_arrays(Cities,ValueCt,Size) :- 
    Size is Cities*ValueCt,
    Dim is Size+1,
    (current_array(dynct(_Dim),_) -> erase_array(dynct/1) ; true),
    (local array(dynct(Dim))),
    (for(J,0,Dim-1) do setval(dynct(J),0)),
    (current_array(neighct(_Dim,_Dim),_) -> erase_array(neighct/2) ; true),
    (local array(neighct(Dim,Dim))),
    (current_array(r_kge(_Dim),_) -> erase_array(r_kge/1) ; true),
    (local array(r_kge(Dim))),
    (for(J,0,Dim-1) do setval(r_kge(J),0)),
    (multifor([K,J],[0,0],[Dim-1,Dim-1]) do setval(neighct(K,J),0)),
    (current_array(tspp(_Dim),_) -> erase_array(tspp/1) ; true),
    (local array(tspp(Dim),float)),
    (current_array(tsppn(_Dim),_) -> erase_array(tsppn/1) ; true),
    (local array(tsppn(Dim),float)).

% Route returns, one at a time on backtracking, all possible
% Hamiltonian routes through Cities points.  
% It record in dynct the number of routes of each length
% It also uses all_neigh to record the number of neighbours which are
% recorded in neighct
route(Cities,Dist) :-
    (for(C,2,Cities),foreach(C,CList) do true),
    (fromto(CList,This,Next,[]), fromto(1,T,N,Last), 
     foreach(T,TRoute),foreach(D,Dists),
     param(Dist)
    do
        remove(C,This,Next), D is Dist[T,C], N=C
    ),
    append(TRoute,[Last],Route),
    DLast is Dist[Last,1],
    Cost is sum(Dists)+DLast,
    Index is Cost,
    incval(dynct(Index)),
    all_neigh(Dist,Route,Cost).

% For a given route, two_swap finds all its neighbours by swapping two
% points and reversing the routes for all the points between them.
% Assuming the distnace matrix is symmetrical, it computes the change
% in length (D3+D4)-(D1+D2), and from that the new cost.
% It increments the relevant neighct(I,J) and also increments the
% total number of neighbours neighct(I,0) of routes with length I+Cities-1
all_neigh(Dist,Route,Cost) :-  
        two_swap(Route,C1,C2,Prev,Succ),
        D1 is Dist[Prev,C1], D2 is Dist[C2,Succ], D3 is Dist[Prev,C2],
        D4 is Dist[C1,Succ],
        NewCost is Cost + D3+D4 - (D1+D2), 
        Index1 is Cost,
        Index2 is NewCost,
        incval(neighct(Index1,Index2)),
        incval(neighct(Index1,0)),
        true,
        fail.
all_neigh(_,_,_ ). 

% Swap any pair of distinct cities, C1 and C2, but assuming the chain
% C1->C2 is reversed, C2 cannot be the predecessor of C1 or even the
% predecessor of the predecessor of C1 because in either case the new
% route will have the same length (assuming symmetric distances)
two_swap(Route,C1,C2,Prev,Succ) :-          
        Route = [C1,Prev|_Tail], append(_,[Succ,C2],Route),
        Prev \= Succ.
two_swap(Route,C1,C2,Prev,Succ) :-          
        Route = [C1|Tail], append(_,[Prev],Route),
        append(_, [C2,Succ,_|_],Tail). 
two_swap(Route,C1,C2,Prev,Succ) :-           
        append(_,[_,Prev,C1|Tail],Route), 
        append(_,[C2],Tail), Route = [Succ|_].
two_swap(Route,C1,C2,Prev,Succ) :-          
        append(_,[Prev,C1|Tail],Route), 
        append(_, [C2,Succ|_],Tail).

record_p_pn(Size) :- 
    writeln('recording...'),
    getval(min_max_mode,(Min,_Max,Mode)),
%    GE is (Min+Mode)//2,
    (for(D,1,Size),foreach(Ct,Counts) do getval(dynct(D),Ct)),
    BigSum is sum(Counts),  
    SMin is Min+1,
    ( for(I,SMin,Mode), param(BigSum,Min)
    do  
        PrevI is I-1,
        ( for(F,Min,PrevI),foreach(Ct,Cts),foreach(NCt,NCts),param(I)
        do
            getval(dynct(F),Ct),
            getval(neighct(I,F),NCt)
        ),
        PI is sum(Cts)/BigSum,
        getval(neighct(I,0),INCt),
        (INCt=0 ->
            PNI=0.0 ;
            PNI is sum(NCts)/INCt
        ),
        setval(tspp(I),PI),
        setval(tsppn(I),PNI)
    ).

intin(L,_,L).
intin(L,Max,K) :-
        L<Max, Next is L+1,
        intin(Next,Max,K).

show_pn_kge(Cities,ValueCt,Count) :- 
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Max is GE-Min,
    getval(neighct(GE,0),GENeighs),
    (for(D,1,Max), foreach(Pn,PnList), 
     param(GE,GENeighs)
    do 
        IPN is GE-D,
        getval(neighct(GE,IPN),Ct),
        Pn is Ct/GENeighs
    ),
    create_filename(pn_kge,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,1,Max),foreach(Pn,PnList),param(S)
    do
        write(S,D),write(S,'    '), 
        approx(Pn,1000,APn),
        write(S,APn),nl(S)
    ),
    close(S).

show_p_r_kge(Cities,ValueCt,Count) :- 
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Max is GE-Min,
    getval(dynct(0),Total),
    (for(I,1,Max), foreach(PXR,PXRList), 
     param(GE,Total)
    do 
        IP is GE-I,
        getval(dynct(IP),Ct),
        getval(r_kge(I),R),
        PXR is Ct*R/Total
    ),
    create_filename(p_r_kge,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,1,Max),foreach(PXR,PXRList),param(S)
    do
        write(S,D),write(S,'    '), 
        approx(PXR,1000,APXR),
        write(S,APXR),nl(S)
    ),
    close(S).
    

approx_list(WList,M,AWList) :-
        foreach(W,WList),foreach(AW,AWList), param(M)
        do  approx(W,M,AW).

approx(W,M,AW) :-
            WW is eval(W),
            I is WW*M, II is fix(round(I)),
            AW is II/M.

show_p_pn(Cities,ValueCt,Count) :-
    create_filename(p_pn,Cities,ValueCt,Count,FileName),
    getval(min_max_mode,(Min,_Max,Mode)),
    open(FileName,write,S),
    write(S,'xx    yy    zz'),nl(S),
    (for(D,Min,Mode), 
     param(S) %,Size)
    do 
        getval(tspp(D),PD),
        getval(tsppn(D),PND),
        write(S,D), write(S,'    '), 
        approx(PD,1000,APD),
        write(S,APD), write(S,'    '), 
        approx(PND,1000,APND),
        write(S,APND), nl(S)
    ),
    close(S).
  

show_p(Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,Max,_Mode)),
    (for(D,Min,Max), foreach(P,PList), 
     param(Total)
    do 
        getval(dynct(D),Ct),
        P is Ct/Total
    ),
    create_filename(p,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,Min,Max),foreach(P,PList),param(S)
    do
        write(S,D),write(S,'    '), 
        approx(P,1000,AP),
        write(S,AP),nl(S)
    ),
    close(S).

% To show decreasing r values
% Find modal fitness value (with the most solutions) Mode
% Find "good enough" value GE
% For each Delta in 1..GE, 
%   Find ct(GE+Delta)+ct(GE-Delta)/Bigsum = p(GE+-Delta)
%   Find neighct(GE,GE+Delta)+neighct(GE,GE-Delta) / neighct(GE,0)
%        = pn(GE,GE+-Delta)
%   Calculate r(GE,Delta) = pn(GE,GE+-Delta)/p(GE+-Delta)
% Write Delta and r(GE,Delta) to FileName
show_r_kge(Cities,ValueCt,Count) :-
    record_r_kge(MaxDiff),
    create_filename(r_kge,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,1,MaxDiff),param(S)
    do
        getval(r_kge(D),R),
        write(S,D),write(S,'    '),
        approx(R,1000,AR),
        write(S,AR),nl(S)
    ),
    close(S).

record_r_kge(MaxDiff) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    MaxDiff is GE-Min, %ValueCt*2,
    GENeigh is getval(neighct(GE,0)),
    (for(D,1,MaxDiff), foreach(R,RList), 
     param(Total,GE,GENeigh)
    do 
        ct_points_at_dist(GE,D,Ct),
        ct_neigh_at_dist(GE,D,NCt),
        P is Ct/Total,
        PN is NCt/GENeigh,
        (Ct>0 -> R is PN/P ; var(R))
    ),
    set_vars(RList),
    (for(I,1,MaxDiff),foreach(R,RList) do setval(r_kge(I),R)).


show_low_pn(Cities,ValueCt,Count) :-
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(PnLow,PnLowList), param(Min)
    do
        KNeigh is getval(neighct(K,0)),
        MaxLow is K-1,
        (for(I,Min,MaxLow), foreach(NCt,NCtList), param(K)
        do 
            getval(neighct(K,I),NCt)
        ),
        PnLow is sum(NCtList)/KNeigh
    ),
    create_filename(pn_low,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(K,Low,GE),foreach(PnLow,PnLowList),param(S)
    do
        write(S,K),write(S,'    '),
        approx(PnLow,1000,APnLow),
        write(S,APnLow),nl(S)
    ),
    close(S).

show_low_pbr(Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(PBRK,PBRList), param(Total,Min)
    do
        calc_r_list(K,Total,Min,RList),
        MaxDiff is K-Min,
        (for(D,1,MaxDiff), foreach(R,RList), foreach(PR,PRKList),
         param(Total,K)
        do 
            LowFit is K-D,
            getval(dynct(LowFit),Ct),
            P is Ct/Total,
            PR is P*R
        ),
        PBRK is sum(PRKList)
    ),
    create_filename(pbr_low,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(K,Low,GE),foreach(PBR,PBRList),param(S)
    do
        write(S,K),write(S,'    '),
        approx(PBR,1000,APBR),
        write(S,APBR),nl(S)
    ),
    close(S).

show_high_pbr(Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(PBRK,PBRList), param(Total,Max)
    do
        KNeigh is getval(neighct(K,0)),
        MaxDiff is Max-K,
        (for(D,1,MaxDiff), foreach(R,RList), 
             param(Total,K,KNeigh)
        do 
            ct_points_at_dist(K,D,Ct),
            ct_neigh_at_dist(K,D,NCt),
            P is Ct/Total,
            PN is NCt/KNeigh,
            (Ct>0 -> R is PN/P ; var(R))
        ),
        set_vars(RList),
        (for(D,1,MaxDiff), foreach(R,RList), foreach(PR,PRKList),
         param(Total,K)
        do 
            HighFit is K+D,
            getval(dynct(HighFit),Ct),
            P is Ct/Total,
            PR is P*R
        ),
        PBRK is sum(PRKList)
    ),
    create_filename(pbr_high,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,Low,GE),foreach(PBR,PBRList),param(S)
    do
        write(S,D),write(S,'    '),
        approx(PBR,1000,APBR),
        write(S,APBR),nl(S)
    ),
    close(S).

show_r_av(Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(AvR,AvRList), param(Total,Min)
    do  calc_r_list(K,Total,Min,RList),
        AvR is sum(RList)/length(RList)
    ),
    create_filename(avr,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,Low,GE),foreach(AvR,AvRList),param(S)
    do
        write(S,D),write(S,'    '),
        approx(AvR,1000,AAvR),
        write(S,AAvR),nl(S)
    ),
    close(S).

show_low_pavr(Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(PAvR,PAvRList), param(Total,Min)
    do   calc_r_list(K,Total,Min,RList),
         AvR is sum(RList)/length(RList),
         MaxDiff is K-Min,
        (for(D,1,MaxDiff), foreach(PKAvR,PAvRKList),
         param(Total,K,AvR)
        do 
            LowFit is K-D,
            getval(dynct(LowFit),Ct),
            P is Ct/Total,
            PKAvR is P*AvR
        ),
        PAvR is sum(PAvRKList)
    ),
    create_filename(pavr_low,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,Low,GE),foreach(PAvR,PAvRList),param(S)
    do
        write(S,D),write(S,'    '),
        approx(PAvR,1000,APAvR),
        write(S,APAvR),nl(S)
    ),
    close(S).

show_high_pavr(Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(PAvR,PAvRList), param(Total,Max,Min)
    do calc_r_list(K,Total,Min,RList),
       AvR is sum(RList)/length(RList),
        MaxDiff is Max-K,
        (for(D,1,MaxDiff), foreach(PKAvR,PAvRKList),
         param(Total,K,AvR)
        do 
            HighFit is K+D,
            getval(dynct(HighFit),Ct),
            P is Ct/Total,
            PKAvR is P*AvR
        ),
        PAvR is sum(PAvRKList)
    ),
    create_filename(pavr_high,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,Low,GE),foreach(PAvR,PAvRList),param(S)
    do
        write(S,D),write(S,'    '),
        approx(PAvR,1000,APAvR),
        write(S,APAvR),nl(S)
    ),
    close(S).

calc_r_avr(K,RList,AvR) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,_Mode)),
    calc_r_list(K,Total,Min,RList),
    AvR is sum(RList)/length(RList).    

calc_r_list(K,Total,Min,RList) :-
         KNeigh is getval(neighct(K,0)),
         KNeigh > 0, !,
         MaxDiff is K-Min,
        (for(D,1,MaxDiff), foreach(R,RList), 
         param(Total,K,KNeigh)
        do 
            ct_points_at_dist(K,D,Ct),
            ct_neigh_at_dist(K,D,NCt),
            P is Ct/Total,
            PN is NCt/KNeigh,
            (Ct>0 -> R is PN/P ; var(R))
        ),
        set_vars(RList).   
calc_r_list(K,_Total,Min,RList) :-
        % KNeigh=0
        writeln(K-has_no_neighbours),
        MaxDiff is K-Min,
        (for(_D,1,MaxDiff), foreach(1,RList) do true).

% If the list has some variables,
% e.g. [10,X,8,Y,Z,5,4,V,W]
% instantiate the variables with intermediate values
% i.e. [10,9,8,7,6,5,4,4,4]

set_vars([V|Vars]) :- var(V), !,
       set_first_vars([V],Vars).
set_vars([Val|Vars]) :-
        set_mid_vars(Val,[],Vars).

%Handle the case where the first elemenet of the list is a variable
% In this case find the first non-variable element and set the earlier
% variables to this value
set_first_vars(List,[V|Vars]) :- var(V), !,
        set_first_vars([V|List],Vars).
set_first_vars(List,[Val|Vars]) :-
        (foreach(V,List),param(Val) do V=Val),
        set_mid_vars(Val,[],Vars).
set_first_vars(_,[]) :- writeln('Error - all r values variable'), abort.

% Val is the last value previously in the list
% Continue collecting variables until the first actual value
set_mid_vars(Val,List,[V|Vars]) :- var(V), !,
        append(List,[V],NList),
        set_mid_vars(Val,NList,Vars).
% Once we have two values, the variables in between are equally spaced
% between the two values
set_mid_vars(V1,List,[V2|Vars]) :-
        length(List,N),
        Diff is (V1-V2)/(1+N),
        (foreach(V,List),fromto(V1,T,N,_),param(Diff)
        do 
            N is T+Diff, V=N
        ),
        set_mid_vars(V2,[],Vars).
% Deals with the case where the last elements are variables
% They are set to equal the last value in the list
set_mid_vars(Val,List,[]) :-
        (foreach(V,List),param(Val) do V=Val).

count_all_solutions(Size) :-
    (for(D,1,Size),foreach(Ct,Counts) do getval(dynct(D),Ct)),
    Total is sum(Counts),
    setval(dynct(0),Total),
    writeln(total_count:Total).

% Find the longest (Max) and shortest (Min) route recorded in dynct
% and the most frequently occurring route length (Mode)
min_max_mode(Size) :-
        (for(D,1,Size),
         fromto((0,0),(TCt,TMode),(NCt,NMode),(_MaxCt,Mode)),
         fromto(Size,Tmin,Nmin,Min),
         fromto(1,Tmax,Nmax,Max)
        do
            Ct is getval(dynct(D)),
            ( Ct >0 ->
                Nmin is min(Tmin,D),
                Nmax is max(Tmax,D),
                (Ct>TCt -> (NCt,NMode)=(Ct,D) 
                ;          (NCt,NMode)=(TCt,TMode)
                )
                
            ;
                Nmin=Tmin, Nmax=Tmax,(NCt,NMode)=(TCt,TMode)
            )
        ),
        GE is (Min+Mode)//2,
        writeln((min(Min),mode(Mode),max(Max),ge(GE))),
        setval(min_max_mode,(Min,Max,Mode)).

% In case GE-D is less than or equal to zero, there are no points 
% with such a cost, so set the count to 0
ct_points_at_dist(GE,D,Ct) :-
                Low is GE-D,  
                ( Low=<0 -> Ct=0 ;
                     getval(dynct(Low),LowDCt),
                     High is GE+D, getval(dynct(High),HighDCt),
                     Ct is LowDCt+HighDCt
                ).
% More complicated to handle Range>1 which may not be needed
% Counts all the neighbours of GE in the two ranges
% GE+D..GE+D+Range and GE-(D+Range)..GE-D
% In case GE-(D+Range) is less than zero, there are no neighbours 
% with such a cost, so set the neighbour count to 0
ct_neigh_at_dist(GE,D,NCt) :-
                Low is GE-D,  
                (Low=<0 -> NCt=0 ;
                    getval(neighct(GE,Low),LowDNCt),
                    High is GE+D, 
                    getval(neighct(GE,High),HighDNCt),
                    NCt is LowDNCt+HighDNCt
                ).

% File name is <Cities><Range>tspp_r.dat
create_filename(Type,Cities,ValueCt,Count,FileName) :-
    number_string(Cities,CA),
    number_string(ValueCt,VA),
    number_string(Count,CtA),
    append_strings(CA,VA,CVA),
    append_strings(CVA,CtA,CCtA),
    append_strings(CCtA,".dat",Tail),
    type_string(Type,String),
    append_strings(String,Tail,FileString),
    atom_string(FileName,FileString).

type_string(pavr_high,"tsp_multi_pavr_high").
type_string(pbr_high,"tsp_multi_pbr_high").
type_string(pavr_low,"tsp_multi_pavr_low").
type_string(pn_low,"tsp_multi_pn_low").
type_string(pbr_low,"tsp_multi_pbr_low").
type_string(avr,"tsp_multi_av_r").
type_string(r_kge,"tsp_multi_r_kge").
type_string(pn_kge,"tsp_multi_pn_kge").
type_string(p_r_kge,"tsp_multi_p_r_kge").
type_string(p,"tsp_multi_p").
type_string(p_pn,"tsp_multi_ppn").

remove(H,[H|T],T).
remove(X,[H|T],[H|R]) :- remove(X,T,R).

all(Goal) :- call(Goal), fail.
all(_).

genint(Min,_,Min).
genint(Min,Max,Val) :-
       Min<Max,N is Min+1, genint(N,Max,Val).

% Generate random values in the range 0..Max
% so that the probability a value occurs is normally distributed
% around Max/2
nrand(Max,Res) :-
        init_nd(ND),
        frandom(R0), (R0>=0.5 -> R1=R0 ; R1=1-R0),
        lookup(R1,ND,0,Row,Col),
        Val is ((Row*10+Col)*(Max-1))/(360*2),
        Mid is (Max-1)/2,
        (R0>=0.5 -> Res is fix(round(Mid-Val)) ; Res is fix(round(Mid+Val))).

lookup(R,ND,N,N,Col) :- (N=36 ; M is N+1, R=<ND[M,1]), !, row(R,ND,N,0,Col).
lookup(R,ND,N,Row,Col) :- M is N+1, lookup(R,ND,M,Row,Col).
row(R,ND,N,C,C) :- (C=9 ;D is C+1, R=<ND[N,D]), !.
row(R,ND,N,C,Row) :- D is C+1, row(R,ND,N,D,Row).
    

init_nd(ND) :- ND= [](
[](0.5000, 0.5040, 0.5080, 0.5120, 0.5160, 0.5199, 0.5239, 0.5279, 0.5319, 0.5359 ),
[](0.5398, 0.5438, 0.5478, 0.5517, 0.5557, 0.5596, 0.5636, 0.5675, 0.5714, 0.5753 ),
[](0.5793, 0.5832, 0.5871, 0.5910, 0.5948, 0.5987, 0.6026, 0.6064, 0.6103, 0.6141 ),
[](0.6179, 0.6217, 0.6255, 0.6293, 0.6331, 0.6368, 0.6406, 0.6443, 0.6480, 0.6517 ),
[](0.6554, 0.6591, 0.6628, 0.6664, 0.6700, 0.6736, 0.6772, 0.6808, 0.6844, 0.6879 ),
[](0.6915, 0.6950, 0.6985, 0.7019, 0.7054, 0.7088, 0.7123, 0.7157, 0.7190, 0.7224 ),
[](0.7257, 0.7291, 0.7324, 0.7357, 0.7389, 0.7422, 0.7454, 0.7486, 0.7517, 0.7549 ),
[](0.7580, 0.7611, 0.7642, 0.7673, 0.7704, 0.7734, 0.7764, 0.7794, 0.7823, 0.7852 ),
[](0.7881, 0.7910, 0.7939, 0.7967, 0.7995, 0.8023, 0.8051, 0.8078, 0.8106, 0.8133 ),
[](0.8159, 0.8186, 0.8212, 0.8238, 0.8264, 0.8289, 0.8315, 0.8340, 0.8365, 0.8389 ),
[](0.8413, 0.8438, 0.8461, 0.8485, 0.8508, 0.8531, 0.8554, 0.8577, 0.8599, 0.8621 ),
[](0.8643, 0.8665, 0.8686, 0.8708, 0.8729, 0.8749, 0.8770, 0.8790, 0.8810, 0.8830 ),
[](0.8849, 0.8869, 0.8888, 0.8907, 0.8925, 0.8944, 0.8962, 0.8980, 0.8997, 0.9015 ),
[](0.9032, 0.9049, 0.9066, 0.9082, 0.9099, 0.9115, 0.9131, 0.9147, 0.9162, 0.9177 ),
[](0.9192, 0.9207, 0.9222, 0.9236, 0.9251, 0.9265, 0.9279, 0.9292, 0.9306, 0.9319 ),
[](0.9332, 0.9345, 0.9357, 0.9370, 0.9382, 0.9394, 0.9406, 0.9418, 0.9429, 0.9441 ),
[](0.9452, 0.9463, 0.9474, 0.9484, 0.9495, 0.9505, 0.9515, 0.9525, 0.9535, 0.9545 ),
[](0.9554, 0.9564, 0.9573, 0.9582, 0.9591, 0.9599, 0.9608, 0.9616, 0.9625, 0.9633 ),
[](0.9641, 0.9649, 0.9656, 0.9664, 0.9671, 0.9678, 0.9686, 0.9693, 0.9699, 0.9706 ),
[](0.9713, 0.9719, 0.9726, 0.9732, 0.9738, 0.9744, 0.9750, 0.9756, 0.9761, 0.9767 ),
[](0.9772, 0.9778, 0.9783, 0.9788, 0.9793, 0.9798, 0.9803, 0.9808, 0.9812, 0.9817 ),
[](0.9821, 0.9826, 0.9830, 0.9834, 0.9838, 0.9842, 0.9846, 0.9850, 0.9854, 0.9857 ),
[](0.9861, 0.9864, 0.9868, 0.9871, 0.9875, 0.9878, 0.9881, 0.9884, 0.9887, 0.9890 ),
[](0.9893, 0.9896, 0.9898, 0.9901, 0.9904, 0.9906, 0.9909, 0.9911, 0.9913, 0.9916 ),
[](0.9918, 0.9920, 0.9922, 0.9925, 0.9927, 0.9929, 0.9931, 0.9932, 0.9934, 0.9936 ),
[](0.9938, 0.9940, 0.9941, 0.9943, 0.9945, 0.9946, 0.9948, 0.9949, 0.9951, 0.9952 ),
[](0.9953, 0.9955, 0.9956, 0.9957, 0.9959, 0.9960, 0.9961, 0.9962, 0.9963, 0.9964 ),
[](0.9965, 0.9966, 0.9967, 0.9968, 0.9969, 0.9970, 0.9971, 0.9972, 0.9973, 0.9974 ),
[](0.9974, 0.9975, 0.9976, 0.9977, 0.9977, 0.9978, 0.9979, 0.9979, 0.9980, 0.9981 ),
[](0.9981, 0.9982, 0.9982, 0.9983, 0.9984, 0.9984, 0.9985, 0.9985, 0.9986, 0.9986 ),
[](0.9987, 0.9987, 0.9987, 0.9988, 0.9988, 0.9989, 0.9989, 0.9989, 0.9990, 0.9990 ),
[](0.9990, 0.9991, 0.9991, 0.9991, 0.9992, 0.9992, 0.9992, 0.9992, 0.9993, 0.9993 ),
[](0.9993, 0.9993, 0.9994, 0.9994, 0.9994, 0.9994, 0.9994, 0.9995, 0.9995, 0.9995 ),
[](0.9995, 0.9995, 0.9995, 0.9996, 0.9996, 0.9996, 0.9996, 0.9996, 0.9996, 0.9997 ),
[](0.9997, 0.9997, 0.9997, 0.9997, 0.9997, 0.9997, 0.9997, 0.9997, 0.9997, 0.9998 ),
[](0.9998, 0.9998, 0.9998, 0.9998, 0.9998, 0.9998, 0.9998, 0.9998, 0.9998, 0.9998 ),
[](0.9998, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 1.0)).

