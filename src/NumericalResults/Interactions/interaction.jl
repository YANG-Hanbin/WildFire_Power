gurobiResultList = Dict{Any, Any}()

indexSets = load("data/testData_RTS/demand/indexSets.jld2")["indexSets"]
paramOPF = load("data/testData_RTS/demand/paramOPF.jld2")["paramOPF"]
paramDemand = load("data/testData_RTS/demand/paramDemand.jld2")["paramDemand"]
Ω_rv = load("data/testData_RTS/Ω_rvList.jld2")["Ω_rvList"][200, 13]
indexSets.Ω = [1:length(keys(Ω_rv))...];

pList = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.985, 1.0];



for p in pList
    prob = Dict{Int64, Float64}();
    for ω in indexSets.Ω 
        prob[ω] = p/200;
    end

    prob[33] = (1 - p) + prob[length(keys(Ω_rv))] #empty scenario
    gurobiResultList[p, "end"] = gurobiOptimizeTest!(indexSets, 
                                                        paramDemand, 
                                                        paramOPF, 
                                                        Ω_end,
                                                        prob; 
                                                        mipGap = 1e-1, timelimit = 3000); 
    gurobiResultList[p, "exg"] = gurobiOptimizeTest!(indexSets, 
                                                        paramDemand, 
                                                        paramOPF, 
                                                        Ω_exg,
                                                        prob; 
                                                        mipGap = 1e-1, timelimit = 3000); 
    save("src/NumericalResults/Interactions/gurobiResultList.jld2", "gurobiResultList", gurobiResultList)
end

save("src/NumericalResults/Interactions/gurobiResultList.jld2", "gurobiResultList", gurobiResultList)

Ω_end = Dict{Int64, RandomVariables}()
Ω_exg = Dict{Int64, RandomVariables}()

ub = Dict{Int64, Int64}()
vb = Dict{Int64, Int64}()
Ibb = Dict{Int64, Vector{Int64}}()
Ibg = Dict{Int64, Vector{Int64}}()
Ibl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
for i in keys(Ω_rv[i].ub)
    ub[i] = 0
    vb[i] = 0
    Ibb[i] = []
    Ibg[i] = []
    Ibl[i] = []
end

ug = Dict{Int64, Int64}()
vg = Dict{Int64, Int64}()
Igb = Dict{Int64, Vector{Int64}}()
Igg = Dict{Int64, Vector{Int64}}()
Igl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
for i in keys(Ω_rv[i].ug)
    ug[i] = 0
    vg[i] = 0
    Igb[i] = []
    Igg[i] = []
    Igl[i] = []
end

ul = Dict{Tuple{Int64, Int64}, Int64}()
vl = Dict{Tuple{Int64, Int64}, Int64}()
Ilb = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
Ilg = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
Ill = Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}()
for i in keys(Ω_rv[i].ul)
    ul[i] = 0
    vl[i] = 0
    Ilb[i] = []
    Ilg[i] = []
    Ill[i] = []
end

for i in 1:200
    rv_end = RandomVariables(Ω_rv[i].τ, Ω_rv[i].ub, Ω_rv[i].ug, Ω_rv[i].ul, vb, vg, vl, Ω_rv[i].Ibb, Ω_rv[i].Ibg, Ω_rv[i].Ibl, Ω_rv[i].Igb, Ω_rv[i].Igg, Ω_rv[i].Igl, Ω_rv[i].Ilb, Ω_rv[i].Ilg, Ω_rv[i].Ill)
    rv_exg = RandomVariables(Ω_rv[i].τ, ub, ug, ul, Ω_rv[i].vb, Ω_rv[i].vg, Ω_rv[i].vl, Ibb, Ibg, Ibl, Igb, Igg, Igl, Ilb, Ilg, Ill)

    Ω_end[i] = rv_end
    Ω_exg[i] = rv_exg
end


## ----------------------------------------------- test --------------------------------------------- ##
totalCost = Dict{Any, Float64}()
costDistribution = Dict()

indexSets = load("data/testData_RTS/demand/indexSets.jld2")["indexSets"]
paramOPF = load("data/testData_RTS/demand/paramOPF.jld2")["paramOPF"]
paramDemand = load("data/testData_RTS/demand/paramDemand.jld2")["paramDemand"]
Ω_rv = load("data/testData_RTS/Ω_rv5000.jld2")["Ω_rv"]
indexSets.Ω = [1:length(keys(Ω_rv))...]
prob = Dict{Int64, Float64}();
for ω in 1:5000 
    prob[ω] = 1/5000;
end

k = 1


Ω_end = Dict{Int64, RandomVariables}()
Ω_exg = Dict{Int64, RandomVariables}()

ub = Dict{Int64, Int64}()
vb = Dict{Int64, Int64}()
Ibb = Dict{Int64, Vector{Int64}}()
Ibg = Dict{Int64, Vector{Int64}}()
Ibl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
for i in keys(Ω_rv[i].ub)
    ub[i] = 0
    vb[i] = 0
    Ibb[i] = []
    Ibg[i] = []
    Ibl[i] = []
end

ug = Dict{Int64, Int64}()
vg = Dict{Int64, Int64}()
Igb = Dict{Int64, Vector{Int64}}()
Igg = Dict{Int64, Vector{Int64}}()
Igl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
for i in keys(Ω_rv[i].ug)
    ug[i] = 0
    vg[i] = 0
    Igb[i] = []
    Igg[i] = []
    Igl[i] = []
end

ul = Dict{Tuple{Int64, Int64}, Int64}()
vl = Dict{Tuple{Int64, Int64}, Int64}()
Ilb = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
Ilg = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
Ill = Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}()
for i in keys(Ω_rv[i].ul)
    ul[i] = 0
    vl[i] = 0
    Ilb[i] = []
    Ilg[i] = []
    Ill[i] = []
end

for i in 1:200
    rv_end = RandomVariables(Ω_rv[i].τ, Ω_rv[i].ub, Ω_rv[i].ug, Ω_rv[i].ul, vb, vg, vl, Ω_rv[i].Ibb, Ω_rv[i].Ibg, Ω_rv[i].Ibl, Ω_rv[i].Igb, Ω_rv[i].Igg, Ω_rv[i].Igl, Ω_rv[i].Ilb, Ω_rv[i].Ilg, Ω_rv[i].Ill)
    rv_exg = RandomVariables(Ω_rv[i].τ, ub, ug, ul, Ω_rv[i].vb, Ω_rv[i].vg, Ω_rv[i].vl, Ibb, Ibg, Ibl, Igb, Igg, Igl, Ilb, Ilg, Ill)

    Ω_end[i] = rv_end
    Ω_exg[i] = rv_exg
end



for sol in ["mix", "end", "exg"]
    for testSet in [Ω_rv, Ω_exg, Ω_end]
        for p in pList 
            state_variable = gurobiResultList[p, sol].first_state_variable;

            ## -------------------------------- solve the first stage model -------------------------------- ##
            (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω);
            (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L);

            Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                                "OutputFlag" => 1, 
                                                "Threads" =>0, 
                                                "MIPGap" => 5e-3, 
                                                "TimeLimit" => 600) 
                                                );
                                                
            @variable(Q, θ_angle[B, 1:T]);      ## phase angle of the bus i
            @variable(Q, P[L, 1:T]);       ## real power flow on line l; elements in L is Tuple (i, j)
            @variable(Q, s[G, 1:T]);       ## real power generation at generator g
            @variable(Q, 0 <= x[D, 1:T] <= 1);  ## load shedding
            # @variable(Q, slack_variable[B, 1:T])

            # constraint 3b 3c
            for l in L
                i = l[1];
                j = l[2];
                @constraint(Q, [t in 1:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - state_variable[:zl][l, t] ) ) );
                @constraint(Q, [t in 1:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - state_variable[:zl][l, t] ) ) );
            end

            ## constraint 1d
            @constraint(Q, [l in L, t in 1:T], P[l, t] >= - paramOPF.W[l] * state_variable[:zl][l, t] );
            @constraint(Q, [l in L, t in 1:T], P[l, t] <= paramOPF.W[l] * state_variable[:zl][l, t] );

            ## constraint 1e
            @constraint(Q, [i in B, t in 1:T], sum(s[g, t] for g in Gᵢ[i]) + 
                                                sum(P[(i, j), t] for j in out_L[i]) - 
                                                sum(P[(j, i), t] for j in in_L[i]) 
                                                .>=  sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) );
            
            ## constraint 1f
            @constraint(Q, [g in G, t in 1:T], s[g, t] >= paramOPF.smin[g] * state_variable[:zg][g, t] );
            @constraint(Q, [g in G, t in 1:T], s[g, t] <= paramOPF.smax[g] * state_variable[:zg][g, t] );

            # constraint 1g h i j
            for i in B 
                @constraint(Q, [t in 1:T, d in Dᵢ[i]], state_variable[:zb][i, t] >= x[d, t] );
                @constraint(Q, [t in 1:T, g in Gᵢ[i]], state_variable[:zb][i, t] >= state_variable[:zg][g, t]);
                @constraint(Q, [t in 1:T, j in out_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(i, j), t] );
                @constraint(Q, [t in 1:T, j in in_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(j, i), t] );
            end


            ## objective function
            @objective(Q, Min, sum( prob[ω] * 
                                            ( sum( sum(paramDemand.w[d] * (1 - x[d, t]) for d in D ) for t in 1:testSet[ω].τ - 1 ) )
                                                                                                                                            for ω in Ω )  
                    );
            optimize!(Q);
            state_value  = JuMP.objective_value(Q); 
                                                                                                                                            
            ## stage 2
            c = 0.0;
            for ω in indexSets.Ω
                @info "$p, $ω"
                randomVariables = testSet[ω];
                forward2Info = forward_stage2_model!(indexSets, 
                                                        paramDemand,
                                                        paramOPF,
                                                        randomVariables              
                                                        );

                ẑ = Dict(   :zg => state_variable[:zg][:, randomVariables.τ - 1], 
                            :zb => state_variable[:zb][:, randomVariables.τ - 1], 
                            :zl => state_variable[:zl][:, randomVariables.τ - 1]
                            );
                ## modify the constraints according to the first stage state variables
                forward_stage2_modify_constraints!(indexSets, 
                                                    forward2Info,
                                                    1,
                                                    ẑ,
                                                    randomVariables
                                                    );

                ####################################################### solve the model and display the result ###########################################################
                optimize!(forward2Info.model);
                damageCost = sum(paramDemand.cb[i] * value(forward2Info.νb[i]) for i in indexSets.B) + 
                                        sum(paramDemand.cg[g] * value(forward2Info.νg[g]) for g in indexSets.G) + 
                                                        sum(paramDemand.cl[l] * value(forward2Info.νl[l]) for l in indexSets.L);
                loadShedding = sum( sum(paramDemand.w[d] * (1 - value(x[d, t])) for d in D ) for t in 1:testSet[ω].τ - 1 ) + (JuMP.objective_value(forward2Info.model) - damageCost);
                c = c + prob[ω] * JuMP.objective_value(forward2Info.model);
                costDistribution[sol, testSet, p, ω] = (damageCost = damageCost, loadShedding = loadShedding)
            end

            # sum(costDistribution[p, ω].damageCost for ω in indexSets.Ω)/5000
            # sum(costDistribution[p, ω].loadShedding for ω in indexSets.Ω)/5000


            u = state_value + c;
            totalCost[sol, testSet, p] = u;
            save("src/NumericalResults/Interactions/totalCost.jld2", "totalCost", totalCost)

        end
    end
end

save("src/NumericalResults/Interactions/totalCost.jld2", "totalCost", totalCost)

totalCost = load("src/NumericalResults/Interactions/totalCost.jld2")["totalCost"]

save("src/NumericalResults/Interactions/costDistribution.jld2", "costDistribution", costDistribution)
costDistribution = load("src/NumericalResults/Interactions/costDistribution.jld2")["costDistribution"]


for p in pList
    @info "$p, $(costDistribution[p, 1090])"
end

p = 0.3
round((sum(costDistribution[p, ω].damageCost for ω in 1:5000) - costDistribution[p, 1090].damageCost * 103)/(5000 - 103), digits = 1)
round((sum(costDistribution[p, ω].loadShedding for ω in 1:5000) - costDistribution[p, 1090].loadShedding * 103)/(5000 - 103), digits = 1)

