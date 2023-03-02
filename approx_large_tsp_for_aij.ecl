
% Generate a TSP, estimate p(j) for all j, find K_{ge}, 
% estimate p(Kge,i) for all i
% and thus calculate r(Kge,i) for all i
% I ran sample_tsp(80,25,1,400000)
sample_tsp(Cities,ValueCt,ProbCt,SampleSize) :-
      create_arrays(Cities,ValueCt,Dim),
      seed(1),
      first_prob(Cities,SampleSize,ValueCt,Dim),
      setval(probct,ProbCt),
      ( ProbCt=1 ->
          test_count(Dim,BigSum),
          show_res(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum) ;
        ProbCt<1 ->
          writeln(mustbegtzero(ProbCt)), abort ;
        ProbCt>1 ->                      
          repeat,
              dist_matrix(Cities,ValueCt,Dist),
              sample_counts(SampleSize,Cities,Dist),
              decval(probct),
              getval(probct,1), !,
          test_count(Dim,BigSum),
          show_res(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum)
       ).

test_count(Dim,BigSum) :-
    find_min_max(Dim,Min,Max,Mode),
    (for(D,Min,Max),foreach(Ct,Counts) do getval(dynct(D),Ct)),
    BigSum is sum(Counts),
    (BigSum = 0 -> writeln(BigSum=0), abort ; 
                   setval(dynct(0),BigSum),
                   writeln(total_count:BigSum)
    ),
    calc_rlists(BigSum,Min,Max,Mode).  

calc_rlists(Total,Min,Max,Mode) :-
    Low is Min+1,
    GE is (Min+Mode)//2,
    ( for(K,Low,GE),param(Total,Min,Max)
    do
        MaxLowDiff is K-Min,
        calc_rlist(K,Total,MaxLowDiff,RListLow),
        setval(rlistlow(K),RListLow),
        MaxHighDiff is Max-K,
        calc_rlist(K,Total,MaxHighDiff,RListHigh),
        setval(rlisthigh(K),RListHigh)
    ).

calc_rlist(K,Total,MaxDiff,RList) :-
        KNeigh is getval(neighct(K,0)),
        (for(D,1,MaxDiff), foreach(R,RList), 
         param(Total,K,KNeigh)
        do 
            (KNeigh=0 -> var(R) 
            ;
                ct_points_at_dist(K,D,Ct),
                ct_neigh_at_dist(K,D,NCt),
                P is Ct/Total,
                PN is NCt/KNeigh,
                (Ct>0 -> R is PN/P ; var(R))
            )
        ),
        set_vars(RList).   

show_res(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum) :-
          show_p(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum),
          show_impp_imppn(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum),
          show_wk_cond(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum),
          show_r_kge(Cities,ValueCt,ProbCt),
          show_pn_kge(Cities,ValueCt,ProbCt),
          show_p_r_kge(Cities,ValueCt,ProbCt),
          show_r_av(Cities,ValueCt,ProbCt),
          show_pbr(low,Cities,ValueCt,ProbCt),
          show_pavr(low,Cities,ValueCt,ProbCt),
          show_pbr(high,Cities,ValueCt,ProbCt),
          show_pavr(high,Cities,ValueCt,ProbCt).

% The array dynct will hold the number of solutions for each cost
% The array neighct will hold the number of neighbours linking each pair
% of costs
% Thus neighct(I,J) records the number of pairs of neighbouring
% solutions S1 and S2 where S1 has cost I and S2 has cost J  
create_arrays(Cities,ValueCt,Size) :- 
    Size is Cities*ValueCt,
    Dim is Size+1,
    (current_array(dynct(_Dim),_) -> erase_array(dynct/1) ; true),
    (local array(dynct(Dim),integer)),
    (for(J,0,Dim-1) do setval(dynct(J),0)),
    (current_array(neighct(_Dim,_Dim),_) -> erase_array(neighct/2) ; true),
    (local array(neighct(Dim,Dim))),
    (multifor([K,J],[0,0],[Dim-1,Dim-1]) do setval(neighct(K,J),0)),
    (current_array(r_kge(_Dim),_) -> erase_array(r_kge/1) ; true),
    (local array(r_kge(Dim))),
    (for(J,0,Dim-1) do setval(r_kge(J),0)),
    (current_array(rlistlow(_Dim),_) -> erase_array(rlistlow/1) ; true),
    (local array(rlistlow(Dim))),
    (for(J,0,Dim-1) do setval(rlistlow(J),[])),
    (current_array(rlisthigh(_Dim),_) -> erase_array(rlisthigh/1) ; true),
    (local array(rlisthigh(Dim))),
    (for(J,0,Dim-1) do setval(rlisthigh(J),[])),    
    (current_array(tspp(_Dim),_) -> erase_array(tspp/1) ; true),
    (local array(tspp(Dim),float)),
    (current_array(tsppn(_Dim),_) -> erase_array(tsppn/1) ; true),
    (local array(tsppn(Dim),float)).

first_prob(Cities,SampleSize,ValueCt,Dim) :-
     dist_matrix(Cities,ValueCt,Dist),
     sample_counts(SampleSize,Cities,Dist),
     find_min_max(Dim,Min,Max,Mode),
     GE is (Mode+Min)//2,
     writeln((min(Min),max(Max),mode(Mode),ge(GE))).

% Dist matrix records the (random) distance (in 1..ValueCt) between each
% pair of cities
% Assign costs to only the edges [I,J] with I<J.
% For edge [I,I] the cost is 0
% Assuming a symmetric tsp, the costs are the same in both directions
dist_matrix(Cities,ValueCt,Dist) :-
        %seed(1),  % Make the solution repeatable
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

% Generate SampleSize random routes, and count the routes with each
% level of fitness
%% sample_counts(SampleSize,Cities,Dist) :- 
%%       (for(_,1,SampleSize),fromto([],List,NList,_),param(Cities,Dist)
%%       do
%%           random_route(Cities,Dist,Cost,Route),
%%           ( add_route(Route,List,NList) ->
%%             all_neigh(Route,Cost,Dist) 
%%           ; 
%%             NList=List
%%           ),
%%           incval(dynct(Cost))
%%       ).
sample_counts(SampleSize,Cities,Dist) :- 
      setval(ssize,0),
      repeat,
      random_route(Cities,Dist,Cost,Route),
      all_neigh(Route,Cost,Dist),
      incval(dynct(Cost)),
      incval(ssize),
      getval(ssize,SampleSize).


% Find the lowest and highest fitness with a non-zero count.
% Assume the modal fitness is halfway between the above Min and Max
% and set Kge to halway between the Min and the modal fitness.
find_min_max(Dim,Min,Max,Mode) :-     
      (for(I,1,Dim),
       fromto(Dim,TMin,NMin,Min),
       fromto(0,TMax,NMax,Max),
       fromto(0-0,TCt-TMode,NCt-NMode,_MaxCt-Mode)
      do
          getval(dynct(I),Val),
          (Val < 20 -> NMin=TMin, NMax=TMax ;
                       NMin is min(TMin,I),
                       NMax is max(TMax,I)
          ),
          (Val>TCt -> (NCt-NMode) = (Val-I) ;
                       (NCt-NMode) = (TCt-TMode)
          )
      ),
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
                    High is GE+D, getval(neighct(GE,High),HighDNCt),
                    NCt is LowDNCt+HighDNCt
                ).


% Finding all neighbours means that around the best randomly generated
% solution, additional, even better, neighbours can be elicited
record_p_pn(Dim,ValueCt,BigSum) :- 
    find_min_max(Dim,Min,_Max,Mode),
    ( for(I,Min,Mode), foreach(PNI,PNIList),param(BigSum,ValueCt)
    do  
        (getval(dynct(I),0) -> true
        ;
            PrevI is I-1,
            NMin is I-2*ValueCt,
            ( for(F,NMin,PrevI),foreach(Ct,Cts),foreach(NCt,NCts),param(I)
            do
              getval(dynct(F),Ct),
              getval(neighct(I,F),NCt)
            ),
            PI is sum(Cts)/BigSum,
            getval(neighct(I,0),INCt),
            (INCt=0 ->
              true ;  %PNI=0.0 ;
              PNI is sum(NCts)/INCt
            ),
            setval(tspp(I),PI) 
        )
    ),
    set_vars(PNIList),
    (foreach(PNI,PNIList),for(I,Min,Mode) do setval(tsppn(I),PNI)).



approx(W,M,AW) :-
        approx(W,20,M,AW).
approx(W,N,M,AW) :-
        W>10^N, ! , 
        Val is (M*W)/(10^N),
        IVal is fix(round(Val)),
        AW is (10^N * IVal)/M.
approx(W,N,M,AW) :-
        N > -4, !,
        PrevN is N-1,
        approx(W,PrevN,M,AW).
approx(_,_,_,0.0).
        
calc_avrpg_rpg_avrpl_rpl(ValueCt,Dim,BigSum,AvRs,PnLs,AvRPLs,PnHs,AvRPHs) :-      
    DMax is ValueCt*2,
    find_min_max(Dim,Min,Max,Mode),
    NMin is max(0,Min-DMax),
    NMax is min(Dim-1,Max+DMax),
    ( for(K,Min,Mode), 
      foreach(AvR,AvRs),
      foreach(AvRPL,AvRPLs),foreach(PnL,PnLs),
      foreach(AvRPH,AvRPHs),foreach(PnH,PnHs),
      param(BigSum,NMin,NMax)
    do
      (getval(dynct(K),0) -> AvRPL=0,PnL=0,AvRPH=0,PnH=0
      ;
        KNeigh is getval(neighct(K,0)),
        KDeltaLowMax is K-NMin,
        (for(D,1,KDeltaLowMax),
         foreach(R,RList),foreach(P,PLows),foreach(Pn,PnLows),
         param(K,KNeigh,BigSum)
        do
          (KNeigh=0 -> var(R)
          ;
              Pt is K-D,
              getval(dynct(Pt),Ct),
              getval(neighct(K,Pt),NCt),
              P is Ct/BigSum,
              Pn is NCt/KNeigh,
              (Ct>0 -> R is Pn/P ; var(R))
          )
        ),
      set_vars(RList),
      AvR is sum(RList)/K, 
      AvRPL  is AvR*sum(PLows),
      PnL is sum(PnLows),
      KDeltaHighMax is NMax-K, 
       (for(D,1,KDeltaHighMax),
         foreach(P,PHighs),foreach(Pn,PnHighs),
         param(K,KNeigh,BigSum)
        do
          (KNeigh=0 -> true
          ;          
              Pt is K+D,
              getval(dynct(Pt),Ct),
              getval(neighct(K,Pt),NCt),
              P is Ct/BigSum,
              Pn is NCt/KNeigh
          )
        ), 
      AvRPH is AvR*sum(PHighs),
      PnH is sum(PnHighs)
      )
      ).

% This should be the same as pavr_high, pavr_low, pbr_high and pbr_low
%  It actually doesn't seem to be the same - needs checking!

show_wk_cond(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum) :- 
    calc_avrpg_rpg_avrpl_rpl(ValueCt,Dim,BigSum,AvRs,PnLs,AvRPLs,PnHs,AvRPHs),
    show_wk_low_cond(SampleSize,ProbCt,Cities,Dim,AvRs,PnLs,AvRPLs),
    show_wk_high_cond(SampleSize,ProbCt,Cities,Dim,PnHs,AvRPHs).

show_wk_low_cond(SampleSize,ProbCt,Cities,Dim,AvRs,PnLs,AvRPLs) :-
    create_p_filename(weak_low_p_pn,Cities,SampleSize,ProbCt,FileName),
    find_min_max(Dim,Min,_Max,Mode),
    open(FileName,write,S),
    write(S,'ww     xx        yy      zz'),nl(S),
    (for(K,Min,Mode),foreach(AvR,AvRs),foreach(PnL,PnLs),foreach(AvRPL,AvRPLs), 
     param(S) %,Size)
    do 
        (getval(dynct(K),0) -> true 
        ;
            write(S,K), write(S,'    '), 
            approx(AvR,1000,AAvR),
            write(S,AAvR), write(S,'    '),
            approx(PnL,1000,APnL),
            write(S,APnL), write(S,'    '), 
            approx(AvRPL,1000,AAvRPL),
            write(S,AAvRPL), nl(S)
        )
    ),
    close(S). 

show_wk_high_cond(SampleSize,ProbCt,Cities,Dim,PnHs,AvRPHs) :-
    create_p_filename(weak_high_p_pn,Cities,SampleSize,ProbCt,FileName),
    find_min_max(Dim,Min,_Max,Mode),
    open(FileName,write,S),
    write(S,'xx    yy    zz'),nl(S),
    (for(K,Min,Mode),foreach(PnH,PnHs),foreach(AvRPH,AvRPHs), 
     param(S) %,Size)
    do 
        (getval(dynct(K),0) -> true 
        ;
            write(S,K), write(S,'    '), 
            approx(PnH,1000,APnH),
            write(S,APnH), write(S,'    '), 
            approx(AvRPH,1000,AAvRPH),
            write(S,AAvRPH), nl(S)
        )
    ),
    close(S). 

    
show_p(SampleSize,ProbCt,Cities,_ValueCt,Dim,BigSum) :-
     find_min_max(Dim,Min,Max,_Mode),
    ( for(I,Min,Max),foreach(P,PList),param(BigSum)
    do 
        getval(dynct(I),Ct),
        P is Ct/BigSum
    ),
    create_p_filename(p,Cities,SampleSize,ProbCt,FileName),
    open(FileName,write,S),
    ( for(I,Min,Max),foreach(P,PList),param(S)
    do
        write(S,I),write(S,'    '), 
        approx(P,1000,AP),
        write(S,AP),nl(S)
    ),
    close(S).

show_impp_imppn(SampleSize,ProbCt,Cities,ValueCt,Dim,BigSum) :-
    record_p_pn(Dim,ValueCt,BigSum),
    create_p_filename(p_pn,Cities,SampleSize,ProbCt,FileName),
    find_min_max(Dim,Min,_Max,Mode),
    open(FileName,write,S),
    write(S,'xx    yy    zz'),nl(S),
    (for(D,Min,Mode), 
     param(S) %,Size)
    do 
        (getval(dynct(D),0) -> true 
        ;
            getval(tspp(D),PD),
            getval(tsppn(D),PND),
            write(S,D), write(S,'    '), 
            approx(PD,1000,APD),
            write(S,APD), write(S,'    '), 
            approx(PND,1000,APND),
            write(S,APND), nl(S)
        )
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
    record_r_kge(ValueCt,Max),
    create_p_filename(r_kge,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,1,Max),param(S)
    do
        getval(r_kge(D),R),
        write(S,D),write(S,'    '),
        approx(R,1000,AR),
        write(S,AR),nl(S)
    ),
    close(S).

show_pn_kge(Cities,ValueCt,Count) :- 
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Max is min(2*ValueCt,GE-Min),
    getval(neighct(GE,0),GENeighs),  % Assume it is non-zero
    (GENeighs=0 -> writeln("No neighbours for Kge!"), abort ; true),
    (for(D,1,Max), foreach(Pn,PnList), 
     param(GE,GENeighs)
    do 
        IPN is GE-D,
        getval(neighct(GE,IPN),Ct),
        Pn is Ct/GENeighs
    ),
    create_p_filename(pn_kge,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,1,Max),foreach(Pn,PnList),param(S)
    do
        write(S,D),write(S,'    '), 
        approx(Pn,1000,APn),
        write(S,APn),nl(S)
    ),
    close(S).

% p(k) * r(Kge,k) for each k in Min..Kge-1
show_p_r_kge(Cities,ValueCt,Count) :- 
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Max is min(2*ValueCt,GE-Min),
    getval(dynct(0),Total),
    (for(I,1,Max), foreach(PXR,PXRList), 
     param(GE,Total)
    do 
        IP is GE-I,
        getval(dynct(IP),Ct),
        getval(r_kge(I),R),
        PXR is Ct*R/Total
    ),
    create_p_filename(p_r_kge,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,1,Max),foreach(PXR,PXRList),param(S)
    do
        write(S,D),write(S,'    '), 
        approx(PXR,1000,APXR),
        write(S,APXR),nl(S)
    ),
    close(S).
    
record_r_kge(ValueCt,Max) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Max is min(2*ValueCt,GE-Min),
    GENeigh is getval(neighct(GE,0)),
    (GENeigh=0 -> writeln("No neighbours for Kge!"), abort ; true),    
    (for(D,1,Max), foreach(R,RList), 
     param(Total,GE,GENeigh)
    do 
        ct_points_at_dist(GE,D,Ct),
        ct_neigh_at_dist(GE,D,NCt),
        P is Ct/Total,
        PN is NCt/GENeigh,
        (Ct>0 -> R is PN/P ; var(R))
    ),
    set_vars(RList),
    (for(I,1,Max),foreach(R,RList) do setval(r_kge(I),R)).

show_pbr(LH,Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(PBRK,PBRList), param(LH,Total,Min,Max)
    do
        (LH=low ->  MaxDiff is K-Min, getval(rlistlow(K),RList) ;
         LH=high -> MaxDiff is Max-K, getval(rlisthigh(K),RList)
        ),
        (for(D,1,MaxDiff), foreach(R,RList), foreach(PR,PRKList),
         param(LH,Total,K)
        do 
            (LH= low -> Fit is K-D ;
             LH=high -> Fit is K+D),           
            getval(dynct(Fit),Ct),
            getval(dynct(K),ThisCt),
            (ThisCt=0 -> var(PR) ;
                P is Ct/Total,
                PR is P*R
            )
        ),
        set_vars(PRKList),
        (getval(dynct(K),0) -> var(PBRK) ; PBRK is sum(PRKList))
    ),
    set_vars(PBRList),
    (LH=low -> Name = pbr_low ; 
     LH=high ->Name = pbr_high),
    write_to_file(Name,Cities,ValueCt,Count,Low,GE,PBRList).


show_r_av(Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,_Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(AvR,AvRList), param(Total,Min)
    do
         KNeigh is getval(neighct(K,0)),
         MaxDiff is K-Min,
        (for(D,1,MaxDiff), foreach(R,RList), 
         param(Total,K,KNeigh)
        do 
            (KNeigh=0 -> var(R) 
            ;
                ct_points_at_dist(K,D,Ct),
                ct_neigh_at_dist(K,D,NCt),
                P is Ct/Total,
                PN is NCt/KNeigh,
                (Ct>0 -> R is PN/P ; var(R))
            )
        ),
        set_vars(RList),
        AvR is sum(RList)/length(RList)
    ),
    create_p_filename(avr,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,Low,GE),foreach(AvR,AvRList),param(S)
    do
        write(S,D),write(S,'    '),
        approx(AvR,1000,AAvR),
        write(S,AAvR),nl(S)
    ),
    close(S).

show_pavr(LH,Cities,ValueCt,Count) :-
    getval(dynct(0),Total),
    getval(min_max_mode,(Min,Max,Mode)),
    GE is (Min+Mode)//2,
    Low is Min+1,
    ( for(K,Low,GE), foreach(PAvR,PAvRList), param(LH,Total,Min,Max)
    do
        (LH=low ->  MaxDiff is K-Min, getval(rlistlow(K),RList) ;
         LH=high -> MaxDiff is Max-K, getval(rlisthigh(K),RList)
        ),
        AvR is sum(RList)/length(RList),
        (for(D,1,MaxDiff), foreach(PKAvR,PAvRKList),
         param(LH,K)
        do 
            (LH = low -> Fit is K-D ;
             LH = high -> Fit is K+D ),
            getval(dynct(Fit),Ct), 
            ( Ct=0 -> var(PKAvR) ; PKAvR=Ct)
         ),
         set_vars(PAvRKList),
         (AvR=0.0 -> var(PAvR) ; PAvR is (sum(PAvRKList)*AvR)/Total)
    ),
    set_vars(PAvRList),
    (LH=low -> Name = pavr_low ; 
     LH=high ->Name = pavr_high),
    write_to_file(Name,Cities,ValueCt,Count,Low,GE,PAvRList).

write_to_file(Name,Cities,ValueCt,Count,Min,Max,List) :-
    create_p_filename(Name,Cities,ValueCt,Count,FileName),
    open(FileName,write,S),
    ( for(D,Min,Max),foreach(P,List),param(S)
    do
        write(S,D),write(S,'    '),
        approx(P,1000,AP),
        write(S,AP),nl(S)
    ),
    close(S).        

% Make a list of distinct cities, starting with city 1
% The next city is selected randomly from ones not yet chosen,
% and the cost is the distance between the previous last city and the
% new one.  The total cost is the sum of the distances
random_route(Cities,Dist,Cost,Route) :-
    Len is Cities-1, 
    ( for(C,2,Cities),foreach(C,CList) do true ),
    ( fromto(CList,This,Next,[]), fromto(1,T,N,Last),
      for(L,Len,1,-1),
      foreach(T,TRoute),foreach(D,Dists),
      param(Dist)
    do
        remove_random(C,L,This,Next), D is Dist[T,C], N=C
    ),
    append(TRoute,[Last],Route),
    DLast is Dist[Last,1],
    Cost is sum(Dists)+DLast.

% R is a randon number between 1 and L
% Remove the Rth member of the list
remove_random(C,L,This,Next) :-
        R is (random mod L)+1,
        remove_nth(R,C,T-T,This,Next).

remove_nth(1,C,Prev-Tail,[C|This],Prev) :- !, Tail=This.
remove_nth(N,C,Old-Tail,[H|List],Next) :-
        N>1, M is N-1, Tail=[H|More], 
        remove_nth(M,C,Old-More,List,Next).

% Keeps a forest recording all the routes found so far
% A new route is added as a branch in the tree,
% branching from the first point it differs from any stored route
% If the route is already in the tree, it fails so preventing the same
% route being processed more than once

% [H] is the last city in a route
% If it isn't already in the firest add it as a leaf
add_route([H],[],[H]) :- !.
% If it is already there, fail
add_route([H],[H|_],_) :- !, fail.
% If there are still leaves,and H>G, try inserting it further on
add_route([H],[G|T],Last) :- H>G, !, add_route([H],T,Last).
% Otherwise insert it here
add_route([H],[J|T],[H,J|T]) :- !.

% [H|T] is the route, or the rest of the route not yet matched
% If there is nothing at this depth in the forest, add a new subtree
% holding the remainder of this route
add_route([H|T],[],[H-Tree]) :-
        add_route(T,[],Tree).
% If there's a match, insert the remainder of the route into the subtree
add_route([H|T],[H-Tree|Tail],[H-New|Tail]) :- !,
        add_route(T,Tree,New).
% If the first city in the route is larger than the current node in
% the forest, insert the route later
add_route([H|T],[G-Tree|More],[G-Tree|New]) :- H>G, !,
        add_route([H|T],More,New).
% otherwise insert the route here
add_route([H|T],[G-Tree|More],[H-NewTree,G-Tree|More]) :-
        add_route(T,[],NewTree).


% For a given route, two_swap finds all its neighbours by swapping two
% points and reversing the routes for all the points between them.
% Assuming the distance matrix is symmetrical, it computes the change
% in length (D3+D4)-(D1+D2), and from that the new cost.
% It increments the relevant neighct(I,J) and also increments the
% total number of neighbours neighct(I,0) of routes with length I+Cities-1

all_neigh(Route,Cost,Dist) :-  
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
all_neigh(_,_,_). 

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
set_first_vars(List,[]) :- %writeln('Error - all r values variable'), abort.
        (foreach(V,List) do V=0).
% Val is the last value previously in the list
% Continue collecting variables until the first actual value
set_mid_vars(Val,List,[V|Vars]) :- var(V), !,
        append(List,[V],NList),
        set_mid_vars(Val,NList,Vars).
% Once we have two values, the variables in between are equally spaced
% between the two values
set_mid_vars(_V1,[],[V2|Vars]) :- !, 
        set_mid_vars(V2,[],Vars).
set_mid_vars(V1,List,[V2|Vars]) :-
        length(List,Len),
        Diff is (V2-V1)/(1+Len),
        (foreach(V,List),fromto(V1,T,N,_),param(Diff)
        do 
            N is T+Diff, V=N
        ),
        set_mid_vars(V2,[],Vars).
% Deals with the case where the last elements are variables
% They are set to equal the last value in the list
set_mid_vars(Val,List,[]) :-
        (foreach(V,List),param(Val) do V=Val).


% Counts all the points in 
% GE+D and GE-D
% In case GE-D is less than zero, there are no points 
% with such a cost, so set the count to 0
ct_points(GE,D,DCt) :-
       Low is GE-D,  
       ( Low=<0 -> DCt=0 ;
           getval(dynct(Low),LowDCt),
           High is GE+D, getval(dynct(High),HighDCt),
           DCt is LowDCt+HighDCt
       ).

% Counts all the neighbours of GE in
% GE+D and GE-D
% In case GE-D is less than zero, there are no neighbours 
% with such a cost, so set the neighbour count to 0
ct_neighs(GE,D,DNCt) :-
       Low is GE-D,  
       (Low=<0 -> DNCt=0 ;
           getval(neighct(GE,Low),LowDNCt),
           High is GE+D, getval(neighct(GE,High),HighDNCt),
           DNCt is LowDNCt+HighDNCt
       ).


% File name is <Cities><Range>tspp_r.dat
create_big_filename(Cities,Size,ProbCt,FileName) :-
    number_string(Cities,CA),
    append_strings(CA,"_",CCA),
     number_string(Size,SA),
     append_strings(CCA,SA,CSA),
    append_strings(CSA,"_",CCSA),
    append_strings(CCSA,"_",CCSRA),
    number_string(ProbCt,TA),
    append_strings(CCSRA,TA,CTSRA),
    append_strings(CTSRA,".dat",Tail),
    append_strings("tsp_approx_multi_r",Tail,FileString),
    atom_string(FileName,FileString).

create_p_filename(Name,Cities,Size,ProbCt,FileName) :-
    number_string(Cities,CA),
    append_strings(CA,"_",CCA),
     number_string(Size,SA),
     append_strings(CCA,SA,CSA),
    append_strings(CSA,"_",CCSA),
    number_string(ProbCt,TA),
    append_strings(CCSA,TA,CTSRA),
    append_strings(CTSRA,".dat",Tail),
    atom_string(Name,StringName),
    append_strings("tsp_approx_",StringName,FullStringName),
    append_strings(FullStringName,Tail,FileString),
    atom_string(FileName,FileString).


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

