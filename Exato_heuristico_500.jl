include("Dados500.jl")
#include("Dados500.jl")
#include("Dados500.jl")
CPUtime = 1200

#------------------------------------------------------------------------------------------
#Functions que controlam as heurísticas
function Plantio_Colheita_Primeiro_Corte(V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime,Xa1,Xe1,Ta1,Te1)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit",CPUtime)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    set_optimizer_attribute(modelo, "DisplayInterval", 20)

    @variable(modelo, xa[i in Va, j in Ja, h in H],Bin)
    @variable(modelo, ta[i in Va, j in Ja, c in C, m in M, d in D],Bin)
    @variable(modelo, xe[i in Ve, j in Je, h in H],Bin)
    @variable(modelo, te[i in Ve, j in Je, c in C, m in M, d in D],Bin)
    @variable(modelo, y[j in J],Bin)

    if size(Xa1,1) > 0
        for i in Va, j in Ja, h in H
            set_start_value(xa[i,j,h],Xa1[i,j,h])
        end
        for i in Ve, j in Je, h in H
            set_start_value(xe[i,j,h],Xe1[i,j,h])
        end
        for i in Va, j in Ja, m in M, d in D
            set_start_value(ta[i,j,1,m,d],Ta1[i,j,1,m,d])
        end
        for i in Ve, j in Je, m in M, d in D
            set_start_value(te[i,j,1,m,d],Te1[i,j,1,m,d])
        end
    end

    @objective(modelo,Max,
    a1*(
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D)
    ) +
    a2*(
    sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D))
    )

    @constraint(modelo,[j in Ja],
    sum(xa[i,j,h] for i in Va for h in H) == 1
    )

    @constraint(modelo,[j in Je],
    sum(xe[i,j,h] for i in Ve for h in H) == 1
    )

    @constraint(modelo,[i in Va],
    sum(xa[i,j,h] for j in Ja for h in H) <= 0.3*ka
    )

    @constraint(modelo,[i in Ve],
    sum(xe[i,j,h] for j in Je for h in H) <= 0.3*ke
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Ja],
    6*y[j] <= sum(h*xa[i,j,h] for i in Va for h in H)
    )

    @constraint(modelo,[j in Je],
    6*y[j] <= sum(h*xe[i,j,h] for i in Ve for h in H)
    )

    @constraint(modelo,[j in Ja,c in 1],
    sum(i*xa[i,j,h] for i in Va for h in H) == sum(i*ta[i,j,c,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je,c in 1],
    sum(i*xe[i,j,h] for i in Ve for h in H) == sum(i*te[i,j,c,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[j in Ja,c in 1],
    sum(ta[i,j,c,m,d] for i in Va for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Je,c in 1],
    sum(te[i,j,c,m,d] for i in Ve for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) == sum((m-d)*ta[i,j,1,m,d] for i in Va for m in M for d in D) - 6*(1-y[j])
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) == sum((m-d)*te[i,j,1,m,d] for i in Ve for m in M for d in D) - 6*(1-y[j])
    )


    @constraint(modelo,[m in M, c in 1],
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DS[m,c]
    )

    @constraint(modelo,[m in M, c in 1],
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DF[m,c]
    )

    @constraint(modelo,[m in M, c in 1],
    sum(((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) <= Cap[m,c]
    )

    status=optimize!(modelo)


    println("-------------------------------------------------------------")
    println("Dimensionamento do plantio e colheita em cada lote j")
    println("-------------------------------------------------------------")
    println(" Plantio Colheita")
    println(" ------------ -------------------")
    println("   j    i    h          m    d    ciclo")
    println("-------------------------------------------------------------")
    for j in Ja
        @printf("%4.1d %4.1d %4.1d %10.1d %4.1d %6.1d\n",
        j,
        sum(i*value(xa[i,j,h]) for i in Va for h in H),
        sum(h*value(xa[i,j,h]) for i in Va for h in H),
        sum(m*value(ta[i,j,1,m,d]) for i in Va for m in M for d in D),
        sum(d*value(ta[i,j,1,m,d]) for i in Va for m in M for d in D),
        value(y[j])
        )
    end
    for j in Je
        @printf("%4.1d %4.1d %4.1d %10.1d %4.1d %6.1d\n",
        j,
        sum(i*value(xe[i,j,h]) for i in Ve for h in H),
        sum(h*value(xe[i,j,h]) for i in Ve for h in H),
        sum(m*value(te[i,j,1,m,d]) for i in Ve for m in M for d in D),
        sum(d*value(te[i,j,1,m,d]) for i in Ve for m in M for d in D),
        value(y[j])
        )
    end
    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DS")
    println("-------------------------------------------------------------")
    for c in 1, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(SPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(SPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DS[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DF")
    println("-------------------------------------------------------------")
    for c in 1, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(FPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(FPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DF[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c    Moagem      Cap.")
    println("-------------------------------------------------------------")
    for c in 1, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(((1-αa)^(c-1))*PSC[i]*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(((1+αe)^(c-1))*PEC[i]*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D) ,
        Cap[m,c])
    end
    println("-------------------------------------------------------------")

    Ta = value.(ta)
    Xa = value.(xa)
    Te = value.(te)
    Xe = value.(xe)
    Y = value.(y)
    tempo = MOI.get(modelo, MOI.SolveTime())
    return (Xa,Xe,Ta,Te,Y,tempo)
end

function Colheita_anos_seguintes(Ta,Te,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit",CPUtime)
    set_optimizer_attribute(modelo, "MIPGap", 0.005)
    set_optimizer_attribute(modelo, "DisplayInterval", 20)

    @variable(modelo, ta[i in Va, j in Ja, c in C,m in M, d in D],Bin)
    @variable(modelo, te[i in Ve, j in Je, c in C,m in M, d in D],Bin)

    @objective(modelo,Max,
    a1*(
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D)
    ) +
    a2*(
    sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D)
    )
    )

    @constraint(modelo,[i in Va, j in Ja,m in M, d in D],
    ta[i,j,Colheita-1,m,d] == Ta[i,j,Colheita-1,m,d]
    )

    @constraint(modelo,[i in Ve, j in Je,m in M, d in D],
    te[i,j,Colheita-1,m,d] == Te[i,j,Colheita-1,m,d]
    )

    @constraint(modelo,[j in Ja],
    sum(ta[i,j,Colheita,m,d] for i in Va for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Je],
    sum(te[i,j,Colheita,m,d] for i in Ve for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Ja],
    sum(i*ta[i,j,Colheita,m,d] for i in Va for m in M for d in D) == sum(i*Ta[i,j,Colheita-1,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je],
    sum(i*te[i,j,Colheita,m,d] for i in Ve for m in M for d in D) == sum(i*Te[i,j,Colheita-1,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[j in Ja],
    sum((m-d)*ta[i,j,Colheita,m,d] for i in Va for m in M for d in D) == sum(m*Ta[i,j,Colheita-1,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je],
    sum((m-d)*te[i,j,Colheita,m,d] for i in Ve for m in M for d in D) == sum(m*Te[i,j,Colheita-1,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[m in M, c in Colheita],
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DS[m,c]
    )

    @constraint(modelo,[m in M, c in Colheita],
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DF[m,c]
    )

    @constraint(modelo,[m in M, c in Colheita],
    sum(((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) <= Cap[m,c]
    )

    status=optimize!(modelo)

    println("-------------------------------------------------------------")
    println("Dimensionamento do plantio e colheita em cada lote j")
    println("-------------------------------------------------------------")
    println(" Plantio Colheita")
    println(" ------------ -------------------")
    println("   j    i    c    m    d          c+1   m    d ")
    println("-------------------------------------------------------------")
    for j in Ja
        @printf("%4.1d %4.1d %4.1d %4.1d %4.1d  %10.1d %4.1d %4.1d\n",
        j,
        sum(i*value(ta[i,j,Colheita-1,m,d]) for i in Va for m in M for d in D),
        Colheita-1,
        sum(m*value(ta[i,j,Colheita-1,m,d]) for i in Va for m in M for d in D),
        sum(d*value(ta[i,j,Colheita-1,m,d]) for i in Va for m in M for d in D),
        Colheita,
        sum(m*value(ta[i,j,Colheita,m,d]) for i in Va for m in M for d in D),
        sum(d*value(ta[i,j,Colheita,m,d]) for i in Va for m in M for d in D)
        )
    end
    for j in Je
        @printf("%4.1d %4.1d %4.1d %4.1d %4.1d  %10.1d %4.1d %4.1d\n",
        j,
        sum(i*value(te[i,j,Colheita-1,m,d]) for i in Ve for m in M for d in D),
        Colheita-1,
        sum(m*value(te[i,j,Colheita-1,m,d]) for i in Ve for m in M for d in D),
        sum(d*value(te[i,j,Colheita-1,m,d]) for i in Ve for m in M for d in D),
        Colheita,
        sum(m*value(te[i,j,Colheita,m,d]) for i in Ve for m in M for d in D),
        sum(d*value(te[i,j,Colheita,m,d]) for i in Ve for m in M for d in D)
        )
    end

    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DS")
    println("-------------------------------------------------------------")
    for c in Colheita, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(SPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(SPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DS[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DF")
    println("-------------------------------------------------------------")
    for c in Colheita, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(FPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(FPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DF[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c    Moagem      Cap.")
    println("-------------------------------------------------------------")
    for c in Colheita, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(((1-αa)^(c-1))*PSC[i]*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(((1+αe)^(c-1))*PEC[i]*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D) ,
        Cap[m,c])
    end
    println("-------------------------------------------------------------")
    Ta = value.(ta)
    Te = value.(te)
    tempo = MOI.get(modelo, MOI.SolveTime())
    return (Ta,Te,tempo)
end
#------------------------------------------------------------------------------------------
#Functions que controlam as heurísticas na solução compromisso
function Plantio_Colheita_Primeiro_Corte_Compromisso(z1max,z1min,z2max,z2min,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

    Ideal = [z1max,z2min]

    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit",CPUtime)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    set_optimizer_attribute(modelo, "DisplayInterval", 20)

    @variable(modelo, xa[i in Va, j in Ja, h in H],Bin)
    @variable(modelo, ta[i in Va, j in Ja, c in C, m in M, d in D],Bin)
    @variable(modelo, xe[i in Ve, j in Je, h in H],Bin)
    @variable(modelo, te[i in Ve, j in Je, c in C, m in M, d in D],Bin)
    @variable(modelo, y[j in J],Bin)
    @variable(modelo, u >=0)

    @objective(modelo,Min, u +
    0.001*((Ideal[1] - (sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D)))/(z1max - z1min) +
    (-Ideal[2] + sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D))/(z2max-z2min))
    )

    @constraint(modelo,
    (Ideal[1] -
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) -
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D) -
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) -
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D))/(z1max-z1min) <=u
    )

    @constraint(modelo,
    (sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D) - Ideal[2])/(z2max-z2min) <=u
    )

    @constraint(modelo,[j in Ja],
    sum(xa[i,j,h] for i in Va for h in H) == 1
    )

    @constraint(modelo,[j in Je],
    sum(xe[i,j,h] for i in Ve for h in H) == 1
    )

    @constraint(modelo,[i in Va],
    sum(xa[i,j,h] for j in Ja for h in H) <= 0.3*ka
    )

    @constraint(modelo,[i in Ve],
    sum(xe[i,j,h] for j in Je for h in H) <= 0.3*ke
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Ja],
    6*y[j] <= sum(h*xa[i,j,h] for i in Va for h in H)
    )

    @constraint(modelo,[j in Je],
    6*y[j] <= sum(h*xe[i,j,h] for i in Ve for h in H)
    )

    @constraint(modelo,[j in Ja,c in 1],
    sum(i*xa[i,j,h] for i in Va for h in H) == sum(i*ta[i,j,c,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je,c in 1],
    sum(i*xe[i,j,h] for i in Ve for h in H) == sum(i*te[i,j,c,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[j in Ja,c in 1],
    sum(ta[i,j,c,m,d] for i in Va for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Je,c in 1],
    sum(te[i,j,c,m,d] for i in Ve for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) == sum((m-d)*ta[i,j,1,m,d] for i in Va for m in M for d in D) - 6*(1-y[j])
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) == sum((m-d)*te[i,j,1,m,d] for i in Ve for m in M for d in D) - 6*(1-y[j])
    )


    @constraint(modelo,[m in M, c in 1],
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DS[m,c]
    )

    @constraint(modelo,[m in M, c in 1],
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DF[m,c]
    )

    @constraint(modelo,[m in M, c in 1],
    sum(((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) <= Cap[m,c]
    )

    status=optimize!(modelo)


    println("-------------------------------------------------------------")
    println("Dimensionamento do plantio e colheita em cada lote j")
    println("-------------------------------------------------------------")
    println(" Plantio Colheita")
    println(" ------------ -------------------")
    println("   j    i    h          m    d    ciclo")
    println("-------------------------------------------------------------")
    for j in Ja
        @printf("%4.1d %4.1d %4.1d %10.1d %4.1d %6.1d\n",
        j,
        sum(i*value(xa[i,j,h]) for i in Va for h in H),
        sum(h*value(xa[i,j,h]) for i in Va for h in H),
        sum(m*value(ta[i,j,1,m,d]) for i in Va for m in M for d in D),
        sum(d*value(ta[i,j,1,m,d]) for i in Va for m in M for d in D),
        value(y[j])
        )
    end
    for j in Je
        @printf("%4.1d %4.1d %4.1d %10.1d %4.1d %6.1d\n",
        j,
        sum(i*value(xe[i,j,h]) for i in Ve for h in H),
        sum(h*value(xe[i,j,h]) for i in Ve for h in H),
        sum(m*value(te[i,j,1,m,d]) for i in Ve for m in M for d in D),
        sum(d*value(te[i,j,1,m,d]) for i in Ve for m in M for d in D),
        value(y[j])
        )
    end
    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DS")
    println("-------------------------------------------------------------")
    for c in 1, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(SPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(SPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DS[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DF")
    println("-------------------------------------------------------------")
    for c in 1, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(FPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(FPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DF[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c    Moagem      Cap.")
    println("-------------------------------------------------------------")
    for c in 1, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(((1-αa)^(c-1))*PSC[i]*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(((1+αe)^(c-1))*PEC[i]*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D) ,
        Cap[m,c])
    end
    println("-------------------------------------------------------------")

    Ta = value.(ta)
    Xa = value.(xa)
    Te = value.(te)
    Xe = value.(xe)
    Y = value.(y)
    tempo = MOI.get(modelo, MOI.SolveTime())
    return (Xa,Xe,Ta,Te,Y,tempo)
end

function Colheita_anos_seguintes_Compromisso(z1max,z1min,z2max,z2min,Ta,Te,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit",CPUtime)
    set_optimizer_attribute(modelo, "MIPGap", 0.005)
    set_optimizer_attribute(modelo, "DisplayInterval", 20)

    @variable(modelo, ta[i in Va, j in Ja, c in C,m in M, d in D],Bin)
    @variable(modelo, te[i in Ve, j in Je, c in C,m in M, d in D],Bin)
    @variable(modelo, u>=0)

    Ideal = [z1max,z2min]

    @objective(modelo,Min, u +
    0.001*((Ideal[1] - (sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D)))/(z1max - z1min) +
    (-Ideal[2] + sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D))/(z2max-z2min))
    )

    @constraint(modelo,
    (Ideal[1] -
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) -
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D) -
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) -
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D))/(z1max-z1min) <=u
    )


    @constraint(modelo,[i in Va, j in Ja,m in M, d in D],
    ta[i,j,Colheita-1,m,d] == Ta[i,j,Colheita-1,m,d]
    )

    @constraint(modelo,[i in Ve, j in Je,m in M, d in D],
    te[i,j,Colheita-1,m,d] == Te[i,j,Colheita-1,m,d]
    )

    @constraint(modelo,[j in Ja],
    sum(ta[i,j,Colheita,m,d] for i in Va for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Je],
    sum(te[i,j,Colheita,m,d] for i in Ve for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Ja],
    sum(i*ta[i,j,Colheita,m,d] for i in Va for m in M for d in D) == sum(i*Ta[i,j,Colheita-1,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je],
    sum(i*te[i,j,Colheita,m,d] for i in Ve for m in M for d in D) == sum(i*Te[i,j,Colheita-1,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[j in Ja],
    sum((m-d)*ta[i,j,Colheita,m,d] for i in Va for m in M for d in D) == sum(m*Ta[i,j,Colheita-1,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je],
    sum((m-d)*te[i,j,Colheita,m,d] for i in Ve for m in M for d in D) == sum(m*Te[i,j,Colheita-1,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[m in M, c in Colheita],
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DS[m,c]
    )

    @constraint(modelo,[m in M, c in Colheita],
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DF[m,c]
    )

    @constraint(modelo,[m in M, c in Colheita],
    sum(((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) <= Cap[m,c]
    )

    status=optimize!(modelo)

    println("-------------------------------------------------------------")
    println("Dimensionamento do plantio e colheita em cada lote j")
    println("-------------------------------------------------------------")
    println(" Plantio Colheita")
    println(" ------------ -------------------")
    println("   j    i    c    m    d          c+1   m    d ")
    println("-------------------------------------------------------------")
    for j in Ja
        @printf("%4.1d %4.1d %4.1d %4.1d %4.1d  %10.1d %4.1d %4.1d\n",
        j,
        sum(i*value(ta[i,j,Colheita-1,m,d]) for i in Va for m in M for d in D),
        Colheita-1,
        sum(m*value(ta[i,j,Colheita-1,m,d]) for i in Va for m in M for d in D),
        sum(d*value(ta[i,j,Colheita-1,m,d]) for i in Va for m in M for d in D),
        Colheita,
        sum(m*value(ta[i,j,Colheita,m,d]) for i in Va for m in M for d in D),
        sum(d*value(ta[i,j,Colheita,m,d]) for i in Va for m in M for d in D)
        )
    end
    for j in Je
        @printf("%4.1d %4.1d %4.1d %4.1d %4.1d  %10.1d %4.1d %4.1d\n",
        j,
        sum(i*value(te[i,j,Colheita-1,m,d]) for i in Ve for m in M for d in D),
        Colheita-1,
        sum(m*value(te[i,j,Colheita-1,m,d]) for i in Ve for m in M for d in D),
        sum(d*value(te[i,j,Colheita-1,m,d]) for i in Ve for m in M for d in D),
        Colheita,
        sum(m*value(te[i,j,Colheita,m,d]) for i in Ve for m in M for d in D),
        sum(d*value(te[i,j,Colheita,m,d]) for i in Ve for m in M for d in D)
        )
    end

    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DS")
    println("-------------------------------------------------------------")
    for c in Colheita, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(SPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(SPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DS[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c     Prod.       DF")
    println("-------------------------------------------------------------")
    for c in Colheita, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(FPSH(i,d,c)*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(FPEH(i,d,c)*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D),
        DF[m,c])
    end

    println("-------------------------------------------------------------")
    println("   m    c    Moagem      Cap.")
    println("-------------------------------------------------------------")
    for c in Colheita, m in M
        @printf("%4.1d %4.1d %10.2f %10.2f\n",m,c,
        sum(((1-αa)^(c-1))*PSC[i]*L[j]*value(ta[i,j,c,m,d]) for j in Ja for i in Va for d in D) +
        sum(((1+αe)^(c-1))*PEC[i]*L[j]*value(te[i,j,c,m,d]) for j in Je for i in Ve for d in D) ,
        Cap[m,c])
    end
    println("-------------------------------------------------------------")
    Ta = value.(ta)
    Te = value.(te)
    tempo = MOI.get(modelo, MOI.SolveTime())
    return (Ta,Te,tempo)
end

#------------------------------------------------------------------------------------------
#Functions que controlam o método exato
function Exato(Xa1,Xe1,Ta1,Te1,Ta2,Te2,Ta3,Te3,Ta4,Te4,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit",CPUtime)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    set_optimizer_attribute(modelo, "MIPFocus", 1)
    set_optimizer_attribute(modelo, "DisplayInterval", 20)

    @variable(modelo, xa[i in Va, j in Ja, h in H],Bin)
    @variable(modelo, ta[i in Va, j in Ja, c in C,m in M, d in D],Bin)
    @variable(modelo, xe[i in Ve, j in Je, h in H],Bin)
    @variable(modelo, te[i in Ve, j in Je, c in C,m in M, d in D],Bin)
    @variable(modelo, y[j in J],Bin)

    for i in Va, j in Ja, h in H
        set_start_value(xa[i,j,h],Xa1[i,j,h])
    end
    for i in Ve, j in Je, h in H
        set_start_value(xe[i,j,h],Xe1[i,j,h])
    end
    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,1,m,d],Ta1[i,j,1,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,1,m,d],Te1[i,j,1,m,d])
    end

    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,2,m,d],Ta2[i,j,2,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,2,m,d],Te2[i,j,2,m,d])
    end

    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,3,m,d],Ta3[i,j,3,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,3,m,d],Te3[i,j,3,m,d])
    end

    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,4,m,d],Ta4[i,j,4,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,4,m,d],Te4[i,j,4,m,d])
    end

    @objective(modelo,Max,
    a1*(
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D)
    ) +
    a2*(
    sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D)
    )
    )

    @constraint(modelo,[j in Ja],
    sum(xa[i,j,h] for i in Va for h in H) == 1
    )

    @constraint(modelo,[j in Je],
    sum(xe[i,j,h] for i in Ve for h in H) == 1
    )

    @constraint(modelo,[i in Va],
    sum(xa[i,j,h] for j in Ja for h in H) <= 0.3*ka
    )

    @constraint(modelo,[i in Ve],
    sum(xe[i,j,h] for j in Je for h in H) <= 0.3*ke
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Ja],
    6*y[j] <= sum(h*xa[i,j,h] for i in Va for h in H)
    )

    @constraint(modelo,[j in Je],
    6*y[j] <= sum(h*xe[i,j,h] for i in Ve for h in H)
    )

    @constraint(modelo,[j in Ja,c in C],
    sum(i*xa[i,j,h] for i in Va for h in H) == sum(i*ta[i,j,c,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je,c in C],
    sum(i*xe[i,j,h] for i in Ve for h in H) == sum(i*te[i,j,c,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[j in Ja,c in C],
    sum(ta[i,j,c,m,d] for i in Va for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Je,c in C],
    sum(te[i,j,c,m,d] for i in Ve for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) == sum((m-d)*ta[i,j,1,m,d] for i in Va for m in M for d in D) - 6*(1-y[j])
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) == sum((m-d)*te[i,j,1,m,d] for i in Ve for m in M for d in D) - 6*(1-y[j])
    )

    @constraint(modelo,[j in Ja, c in [2,3,4]],
    sum((m-d)*ta[i,j,c,m,d] for i in Va for m in M for d in D) == sum(m*ta[i,j,c-1,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je, c in [2,3,4]],
    sum((m-d)*te[i,j,c,m,d] for i in Ve for m in M for d in D) == sum(m*te[i,j,c-1,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[m in M, c in C],
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DS[m,c]
    )

    @constraint(modelo,[m in M, c in C],
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DF[m,c]
    )

    @constraint(modelo,[m in M, c in C],
    sum(((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) <= Cap[m,c]
    )

    status=optimize!(modelo)
    Ta = value.(ta)
    Te = value.(te)
    Xa = value.(xa)
    Xe = value.(xe)
    Y=value.(y)
    tempo = MOI.get(modelo, MOI.SolveTime())
    Gap = MOI.get(modelo, MOI.RelativeGap())
    return (Xa,Xe,Ta,Te,Y,tempo,Gap)
end

function Exato_Compromisso(z1max,z1min,z2max,z2min,Xa1,Xe1,Ta1,Te1,Ta2,Te2,Ta3,Te3,Ta4,Te4,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit",CPUtime)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    set_optimizer_attribute(modelo, "MIPFocus", 1)
    set_optimizer_attribute(modelo, "DisplayInterval", 20)

    @variable(modelo, xa[i in Va, j in Ja, h in H],Bin)
    @variable(modelo, ta[i in Va, j in Ja, c in C,m in M, d in D],Bin)
    @variable(modelo, xe[i in Ve, j in Je, h in H],Bin)
    @variable(modelo, te[i in Ve, j in Je, c in C,m in M, d in D],Bin)
    @variable(modelo, y[j in J],Bin)
    @variable(modelo, u>=0)

    for i in Va, j in Ja, h in H
        set_start_value(xa[i,j,h],Xa1[i,j,h])
    end
    for i in Ve, j in Je, h in H
        set_start_value(xe[i,j,h],Xe1[i,j,h])
    end
    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,1,m,d],Ta1[i,j,1,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,1,m,d],Te1[i,j,1,m,d])
    end

    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,2,m,d],Ta2[i,j,2,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,2,m,d],Te2[i,j,2,m,d])
    end

    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,3,m,d],Ta3[i,j,3,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,3,m,d],Te3[i,j,3,m,d])
    end

    for i in Va, j in Ja, m in M, d in D
        set_start_value(ta[i,j,4,m,d],Ta4[i,j,4,m,d])
    end
    for i in Ve, j in Je, m in M, d in D
        set_start_value(te[i,j,4,m,d],Te4[i,j,4,m,d])
    end

    Ideal = [z1max,z2min]

    @objective(modelo,Min, u +
    0.001*((Ideal[1] - (sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D)))/(z1max - z1min) +
    (-Ideal[2] + sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D))/(z2max-z2min))
    )

    @constraint(modelo,
    (Ideal[1] -
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) -
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D) -
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) -
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D))
    /(z1max-z1min) <=u
    )

    @constraint(modelo,
    (sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D) - Ideal[2])/(z2max-z2min) <=u
    )

    @constraint(modelo,[j in Ja],
    sum(xa[i,j,h] for i in Va for h in H) == 1
    )

    @constraint(modelo,[j in Je],
    sum(xe[i,j,h] for i in Ve for h in H) == 1
    )

    @constraint(modelo,[i in Va],
    sum(xa[i,j,h] for j in Ja for h in H) <= 0.3*ka
    )

    @constraint(modelo,[i in Ve],
    sum(xe[i,j,h] for j in Je for h in H) <= 0.3*ke
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) - 4 <= 6*y[j]
    )

    @constraint(modelo,[j in Ja],
    6*y[j] <= sum(h*xa[i,j,h] for i in Va for h in H)
    )

    @constraint(modelo,[j in Je],
    6*y[j] <= sum(h*xe[i,j,h] for i in Ve for h in H)
    )

    @constraint(modelo,[j in Ja,c in C],
    sum(i*xa[i,j,h] for i in Va for h in H) == sum(i*ta[i,j,c,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je,c in C],
    sum(i*xe[i,j,h] for i in Ve for h in H) == sum(i*te[i,j,c,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[j in Ja,c in C],
    sum(ta[i,j,c,m,d] for i in Va for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Je,c in C],
    sum(te[i,j,c,m,d] for i in Ve for m in M for d in D) == 1
    )

    @constraint(modelo,[j in Ja],
    sum(h*xa[i,j,h] for i in Va for h in H) == sum((m-d)*ta[i,j,1,m,d] for i in Va for m in M for d in D) - 6*(1-y[j])
    )

    @constraint(modelo,[j in Je],
    sum(h*xe[i,j,h] for i in Ve for h in H) == sum((m-d)*te[i,j,1,m,d] for i in Ve for m in M for d in D) - 6*(1-y[j])
    )

    @constraint(modelo,[j in Ja, c in [2,3,4]],
    sum((m-d)*ta[i,j,c,m,d] for i in Va for m in M for d in D) == sum(m*ta[i,j,c-1,m,d] for i in Va for m in M for d in D)
    )

    @constraint(modelo,[j in Je, c in [2,3,4]],
    sum((m-d)*te[i,j,c,m,d] for i in Ve for m in M for d in D) == sum(m*te[i,j,c-1,m,d] for i in Ve for m in M for d in D)
    )

    @constraint(modelo,[m in M, c in C],
    sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DS[m,c]
    )

    @constraint(modelo,[m in M, c in C],
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) >= DF[m,c]
    )

    @constraint(modelo,[m in M, c in C],
    sum(((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for j in Ja for i in Va for d in D) + sum(((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for j in Je for i in Ve for d in D) <= Cap[m,c]
    )

    status=optimize!(modelo)
    Ta = value.(ta)
    Te = value.(te)
    Xa = value.(xa)
    Xe = value.(xe)
    Y = value.(y)
    U = value.(u)
    tempo = MOI.get(modelo, MOI.SolveTime())
    Gap = MOI.get(modelo, MOI.RelativeGap())
    return (Xa,Xe,Ta,Te,Y,tempo,Gap)
end

#------------------------------------------------------------------------------------------
#Function que calcula os valores objetivos
function Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,Ta_exato,Te_exato,Ta1,Te1,Ta2,Te2,Ta3,Te3,Ta4,Te4,Va,Ja,M,D,Ve,Je,C)

    if size(Ta1,1) > 0
        prod_col_1 = sum(SPSH(i,d,1)*L[j]*Ta1[i,j,1,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(SPEH(i,d,1)*L[j]*Te1[i,j,1,m,d] for i in Ve for j in Je for m in M for d in D)+
        sum(FPSH(i,d,1)*L[j]*Ta1[i,j,1,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(FPEH(i,d,1)*L[j]*Te1[i,j,1,m,d] for i in Ve for j in Je for m in M for d in D)

        prod_col_2 = sum(SPSH(i,d,2)*L[j]*Ta2[i,j,2,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(SPEH(i,d,2)*L[j]*Te2[i,j,2,m,d] for i in Ve for j in Je for m in M for d in D)+
        sum(FPSH(i,d,2)*L[j]*Ta2[i,j,2,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(FPEH(i,d,2)*L[j]*Te2[i,j,2,m,d] for i in Ve for j in Je for m in M for d in D)

        prod_col_3 = sum(SPSH(i,d,3)*L[j]*Ta3[i,j,3,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(SPEH(i,d,3)*L[j]*Te3[i,j,3,m,d] for i in Ve for j in Je for m in M for d in D)+
        sum(FPSH(i,d,3)*L[j]*Ta3[i,j,3,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(FPEH(i,d,3)*L[j]*Te3[i,j,3,m,d] for i in Ve for j in Je for m in M for d in D)

        prod_col_4 = sum(SPSH(i,d,4)*L[j]*Ta4[i,j,4,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(SPEH(i,d,4)*L[j]*Te4[i,j,4,m,d] for i in Ve for j in Je for m in M for d in D)+
        sum(FPSH(i,d,4)*L[j]*Ta4[i,j,4,m,d] for i in Va for j in Ja for m in M for d in D)+
        sum(FPEH(i,d,4)*L[j]*Te4[i,j,4,m,d] for i in Ve for j in Je for m in M for d in D)

        prod_total = [prod_col_1;prod_col_2;prod_col_3;prod_col_4]

        custo_col_1 = sum(CA[j,1]*((1-αa)^(1-1))*PSC[i]*L[j]*Ta1[i,j,1,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,1]*((1+αe)^(1-1))*PEC[i]*L[j]*Te1[i,j,1,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_col_2 = sum(CA[j,2]*((1-αa)^(2-1))*PSC[i]*L[j]*Ta2[i,j,2,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,2]*((1+αe)^(2-1))*PEC[i]*L[j]*Te2[i,j,2,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_col_3 = sum(CA[j,3]*((1-αa)^(3-1))*PSC[i]*L[j]*Ta3[i,j,3,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,3]*((1+αe)^(3-1))*PEC[i]*L[j]*Te3[i,j,3,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_col_4 = sum(CA[j,4]*((1-αa)^(4-1))*PSC[i]*L[j]*Ta4[i,j,4,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,4]*((1+αe)^(4-1))*PEC[i]*L[j]*Te4[i,j,4,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_total = [custo_col_1;custo_col_2;custo_col_3;custo_col_4]

    end
    if size(Ta_exato,1) > 0
        prod_col_1 = sum(SPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D)+
        sum(SPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D)+
        sum(FPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D)+
        sum(FPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D)

        prod_col_2 = sum(SPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 2 for m in M for d in D)+
        sum(SPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 2 for m in M for d in D)+
        sum(FPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 2 for m in M for d in D)+
        sum(FPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 2 for m in M for d in D)

        prod_col_3 = sum(SPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 3 for m in M for d in D)+
        sum(SPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 3 for m in M for d in D)+
        sum(FPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 3 for m in M for d in D)+
        sum(FPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 3 for m in M for d in D)

        prod_col_4 = sum(SPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 4 for m in M for d in D)+
        sum(SPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 4 for m in M for d in D)+
        sum(FPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 4 for m in M for d in D)+
        sum(FPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 4 for m in M for d in D)

        prod_total = [prod_col_1;prod_col_2;prod_col_3;prod_col_4]

        custo_col_1 = sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
        sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D)

        custo_col_2 = sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 2 for m in M for d in D) +
        sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 2 for m in M for d in D)

        custo_col_3 = sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 3 for m in M for d in D) +
        sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 3 for m in M for d in D)

        custo_col_4 = sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in 4 for m in M for d in D) +
        sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in 4 for m in M for d in D)

        custo_total = [custo_col_1;custo_col_2;custo_col_3;custo_col_4]
    end
    return (prod_total,custo_total)
end

#------------------------------------------------------------------------------------------
#Resolução Heurística - Máxima produção
#a1=1 e a2=0,para maximizar produção
a1 = 1
a2 = 0

(Xa1_MP,Xe1_MP,Ta1_MP,Te1_MP,Y1_MP,Tempo1_MP) = Plantio_Colheita_Primeiro_Corte(V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,1.1*DS,1.0*DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime,[],[],[],[])

Colheita = 2
(Ta2_MP,Te2_MP,Tempo2_MP) = Colheita_anos_seguintes(Ta1_MP,Te1_MP,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 3
(Ta3_MP,Te3_MP,Tempo3_MP) = Colheita_anos_seguintes(Ta2_MP,Te2_MP,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 4
(Ta4_MP,Te4_MP,Tempo4_MP) = Colheita_anos_seguintes(Ta3_MP,Te3_MP,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Tempo_MP = Tempo1_MP + Tempo2_MP + Tempo3_MP + Tempo4_MP

FileIO.save("Xa_heuristico_1_MP500.jld2","RH",Xa1_MP)
FileIO.save("Xe_heuristico_1_MP500.jld2","RH",Xe1_MP)
FileIO.save("Y_heuristico_MP500.jld2","RH",Y1_MP)
FileIO.save("Ta_heuristico_1_MP500.jld2","RH",Ta1_MP)
FileIO.save("Ta_heuristico_2_MP500.jld2","RH",Ta2_MP)
FileIO.save("Ta_heuristico_3_MP500.jld2","RH",Ta3_MP)
FileIO.save("Ta_heuristico_4_MP500.jld2","RH",Ta4_MP)
FileIO.save("Te_heuristico_1_MP500.jld2","RH",Te1_MP)
FileIO.save("Te_heuristico_2_MP500.jld2","RH",Te2_MP)
FileIO.save("Te_heuristico_3_MP500.jld2","RH",Te3_MP)
FileIO.save("Te_heuristico_4_MP500.jld2","RH",Te4_MP)
FileIO.save("Tempo_heuristicoMP500.jld2","RH",Tempo_MP)
#------------------------------------------------------------------------------------------
#Resolução Heurística - Mínimo Custo
#a1=0 e a2=-1,para mínimo custo
a1 = 0
a2 = -1

(Xa1_MC,Xe1_MC,Ta1_MC,Te1_MC,Y1_MC,Tempo1_MC) = Plantio_Colheita_Primeiro_Corte(V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,1.1*DS,1.0*DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime,Xa1_MP,Xe1_MP,Ta1_MP,Te1_MP)

Colheita = 2
(Ta2_MC,Te2_MC,Tempo2_MC) = Colheita_anos_seguintes(Ta1_MC,Te1_MC,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 3
(Ta3_MC,Te3_MC,Tempo3_MC) = Colheita_anos_seguintes(Ta2_MC,Te2_MC,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 4
(Ta4_MC,Te4_MC,Tempo4_MC) = Colheita_anos_seguintes(Ta3_MC,Te3_MC,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Tempo_MC = Tempo1_MC + Tempo2_MC + Tempo3_MC + Tempo4_MC

FileIO.save("Xa_heuristico_1_MC500.jld2","RH",Xa1_MC)
FileIO.save("Xe_heuristico_1_MC500.jld2","RH",Xe1_MC)
FileIO.save("Y_heuristico_MC500.jld2","RH",Y1_MC)
FileIO.save("Ta_heuristico_1_MC500.jld2","RH",Ta1_MC)
FileIO.save("Ta_heuristico_2_MC500.jld2","RH",Ta2_MC)
FileIO.save("Ta_heuristico_3_MC500.jld2","RH",Ta3_MC)
FileIO.save("Ta_heuristico_4_MC500.jld2","RH",Ta4_MC)
FileIO.save("Te_heuristico_1_MC500.jld2","RH",Te1_MC)
FileIO.save("Te_heuristico_2_MC500.jld2","RH",Te2_MC)
FileIO.save("Te_heuristico_3_MC500.jld2","RH",Te3_MC)
FileIO.save("Te_heuristico_4_MC500.jld2","RH",Te4_MC)
FileIO.save("Tempo_heuristicoMC500.jld2","RH",Tempo_MC)

#------------------------------------------------------------------------------------------
#Rodamos o exato para determinar as soluções lexicográficas exatas partindo das soluções heurísticas
#Máxima Produção
a1 = 1
a2 = 0
(Xa_exato_MP,Xe_exato_MP,Ta_exato_MP,Te_exato_MP,Y_exato_MP,Tempo_exato_MP,Gap_MP) = Exato(Xa1_MP,Xe1_MP,Ta1_MP,Te1_MP,Ta2_MP,Te2_MP,Ta3_MP,Te3_MP,Ta4_MP,Te4_MP,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

FileIO.save("Xa_exato_MP500.jld2","RH",Xa_exato_MP)
FileIO.save("Xe_exato_MP500.jld2","RH",Xe_exato_MP)
FileIO.save("Y_exato_MP500.jld2","RH",Y_exato_MP)
FileIO.save("Ta_exato_MP500.jld2","RH",Ta_exato_MP)
FileIO.save("Te_exato_MP500.jld2","RH",Te_exato_MP)
FileIO.save("Tempo_exatoMP500.jld2","RH",Tempo_exato_MP)
FileIO.save("Gap_exatoMP500.jld2","RH",Gap_MP)

#Mínimo Custo
a1 = 0
a2 = -1
(Xa_exato_MC,Xe_exato_MC,Ta_exato_MC,Te_exato_MC,Y_exato_MC,Tempo_exato_MC,Gap_MC) = Exato(Xa1_MC,Xe1_MC,Ta1_MC,Te1_MC,Ta2_MC,Te2_MC,Ta3_MC,Te3_MC,Ta4_MC,Te4_MC,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

FileIO.save("Xa_exato_MC500.jld2","RH",Xa_exato_MC)
FileIO.save("Xe_exato_MC500.jld2","RH",Xe_exato_MC)
FileIO.save("Y_exato_MC500.jld2","RH",Y_exato_MC)
FileIO.save("Ta_exato_MC500.jld2","RH",Ta_exato_MC)
FileIO.save("Te_exato_MC500.jld2","RH",Te_exato_MC)
FileIO.save("Tempo_exatoMC500.jld2","RH",Tempo_exato_MC)
FileIO.save("Gap_exatoMC500.jld2","RH",Gap_MC)

#Calculamos as coordenadas do vetor ideal exato
(prod_exato_MP,custo_exato_MP) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,Ta_exato_MP,Te_exato_MP,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C)
(prod_exato_MC,custo_exato_MC) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,Ta_exato_MC,Te_exato_MC,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C)

z1max1 = prod_exato_MP[1]
z1max2 = prod_exato_MP[2]
z1max3 = prod_exato_MP[3]
z1max4 = prod_exato_MP[4]

z1min1 = prod_exato_MC[1]
z1min2 = prod_exato_MC[2]
z1min3 = prod_exato_MC[3]
z1min4 = prod_exato_MC[4]

z2max1 = custo_exato_MP[1]
z2max2 = custo_exato_MP[2]
z2max3 = custo_exato_MP[3]
z2max4 = custo_exato_MP[4]

z2min1 = custo_exato_MC[1]
z2min2 = custo_exato_MC[2]
z2min3 = custo_exato_MC[3]
z2min4 = custo_exato_MC[4]

z1maxT = sum([z1max1,z1max2,z1max3,z1max4])
z1minT = sum([z1min1,z1min2,z1min3,z1min4])

z2maxT = sum([z2max1,z2max2,z2max3,z2max4])
z2minT = sum([z2min1,z2min2,z2min3,z2min4])

#------------------------------------------------------------------------------------------
#Determinamos a solução compromisso aproximada tendo como referência as coordenadas do Ideal Exato
(XaC1,XeC1,TaC1,TeC1,YC1,TempoC1) = Plantio_Colheita_Primeiro_Corte_Compromisso(z1max1,z1min1,z2max1,z2min1,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 2
(TaC2,TeC2,TempoC2) = Colheita_anos_seguintes_Compromisso(z1max2,z1min2,z2max2,z2min2,TaC1,TeC1,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 3
(TaC3,TeC3,TempoC3) = Colheita_anos_seguintes_Compromisso(z1max3,z1min3,z2max3,z2min3,TaC2,TeC2,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 4
(TaC4,TeC4,TempoC4) = Colheita_anos_seguintes_Compromisso(z1max4,z1min4,z2max4,z2min4,TaC3,TeC3,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Tempo_C = TempoC1 + TempoC2 + TempoC3 + TempoC4

FileIO.save("Xa_heuristico_1_CO500.jld2","RH",XaC1)
FileIO.save("Xe_heuristico_1_CO500.jld2","RH",XeC1)
FileIO.save("Y_heuristico_CO500.jld2","RH",YC1)
FileIO.save("Ta_heuristico_1_CO500.jld2","RH",TaC1)
FileIO.save("Ta_heuristico_2_CO500.jld2","RH",TaC2)
FileIO.save("Ta_heuristico_3_CO500.jld2","RH",TaC3)
FileIO.save("Ta_heuristico_4_CO500.jld2","RH",TaC4)
FileIO.save("Te_heuristico_1_CO500.jld2","RH",TeC1)
FileIO.save("Te_heuristico_2_CO500.jld2","RH",TeC2)
FileIO.save("Te_heuristico_3_CO500.jld2","RH",TeC3)
FileIO.save("Te_heuristico_4_CO500.jld2","RH",TeC4)
FileIO.save("Tempo_heuristicoCO500.jld2","RH",Tempo_C)

#------------------------------------------------------------------------------------------
#Determinamos a solução compromisso exata tendo como referência as coordenadas do Ideal Exato e usando a solução compromisso aproximada
(Xa_exato_compromisso,Xe_exato_compromisso,Ta_exato_compromisso,Te_exato_compromisso,Y_exato_compromisso,Tempo_exato_compromisso,Gap_exatoC) = Exato_Compromisso(z1maxT,z1minT,z2maxT,z2minT,XaC1,XeC1,TaC1,TeC1,TaC2,TeC2,TaC3,TeC3,TaC4,TeC4,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

FileIO.save("Xa_exato_CO500.jld2","RH",Xa_exato_compromisso)
FileIO.save("Xe_exato_CO500.jld2","RH",Xe_exato_compromisso)
FileIO.save("Y_exato_CO500.jld2","RH",Y_exato_compromisso)
FileIO.save("Ta_exato_CO500.jld2","RH",Ta_exato_compromisso)
FileIO.save("Te_exato_CO500.jld2","RH",Te_exato_compromisso)
FileIO.save("Tempo_exatoCO500.jld2","RH",Tempo_exato_compromisso)
FileIO.save("Gap_exatoCO500.jld2","RH",Gap_exatoC)

#------------------------------------------------------------------------------------------
