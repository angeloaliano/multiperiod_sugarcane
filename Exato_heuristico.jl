include("Dados30.jl")
#include("Dados65.jl")
#include("Dados150.jl")
CPUtime = 600

#------------------------------------------------------------------------------------------
#Functions que controlam as heurísticas
function Plantio_Colheita_Primeiro_Corte(V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)
    modelo = Model(Gurobi.Optimizer)
    set_optimizer_attribute(modelo, "TimeLimit",CPUtime)
    set_optimizer_attribute(modelo, "MIPGap", 0.01)
    set_optimizer_attribute(modelo, "DisplayInterval", 20)

    @variable(modelo, xa[i in Va, j in Ja, h in H],Bin)
    @variable(modelo, ta[i in Va, j in Ja, c in C, m in M, d in D],Bin)
    @variable(modelo, xe[i in Ve, j in Je, h in H],Bin)
    @variable(modelo, te[i in Ve, j in Je, c in C, m in M, d in D],Bin)
    @variable(modelo, y[j in J],Bin)

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

    Ideal = [z1max,z2min]/4

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
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D)))/(z1max/4 - z1min/4) +
    (-Ideal[2] + sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D))/(z2max/4-z2min/4))
    )

    @constraint(modelo,
    (Ideal[1] - sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D))/(z1max/4-z1min/4) <=u
    )

    @constraint(modelo,
    (sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in 1 for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in 1 for m in M for d in D) - Ideal[2])/(z2max/4-z2min/4) <=u
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

    Ideal = [z1max,z2min]/4

    @objective(modelo,Min, u +
    0.001*((Ideal[1] - (sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D)))/(z1max/4 - z1min/4) +
    (-Ideal[2] + sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D))/(z2max/4-z2min/4))
    )

    @constraint(modelo,
    (Ideal[1] - sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D))/(z1max/4-z1min/4) <=u
    )

    @constraint(modelo,
    (sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in Colheita for m in M for d in D) +
    sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in Colheita for m in M for d in D) - Ideal[2])/(z2max/4-z2min/4) <=u
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
    tempo = MOI.get(modelo, MOI.SolveTime())
    return (Xa,Xe,Ta,Te,tempo)
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
    (Ideal[1] - sum(SPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(SPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D) +
    sum(FPSH(i,d,c)*L[j]*ta[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
    sum(FPEH(i,d,c)*L[j]*te[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D))/(z1max-z1min) <=u
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
    tempo = MOI.get(modelo, MOI.SolveTime())
    return (Xa,Xe,Ta,Te,tempo)
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

        prod_total = prod_col_1 + prod_col_2 + prod_col_3 + prod_col_4

        custo_col_1 = sum(CA[j,1]*((1-αa)^(1-1))*PSC[i]*L[j]*Ta1[i,j,1,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,1]*((1+αe)^(1-1))*PEC[i]*L[j]*Te1[i,j,1,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_col_2 = sum(CA[j,2]*((1-αa)^(2-1))*PSC[i]*L[j]*Ta2[i,j,2,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,2]*((1+αe)^(2-1))*PEC[i]*L[j]*Te2[i,j,2,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_col_3 = sum(CA[j,3]*((1-αa)^(3-1))*PSC[i]*L[j]*Ta3[i,j,3,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,3]*((1+αe)^(3-1))*PEC[i]*L[j]*Te3[i,j,3,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_col_4 = sum(CA[j,4]*((1-αa)^(4-1))*PSC[i]*L[j]*Ta4[i,j,4,m,d] for i in Va for j in Ja for m in M for d in D) +
        sum(CE[j,4]*((1+αe)^(4-1))*PEC[i]*L[j]*Te4[i,j,4,m,d] for i in Ve for j in Je for m in M for d in D)

        custo_total = custo_col_1 + custo_col_2 + custo_col_3 + custo_col_4
    end
    if size(Ta_exato,1) > 0
        prod_total = sum(SPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D)+
        sum(SPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D)+
        sum(FPSH(i,d,c)*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D)+
        sum(FPEH(i,d,c)*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D)

        custo_total = sum(CA[j,c]*((1-αa)^(c-1))*PSC[i]*L[j]*Ta_exato[i,j,c,m,d] for i in Va for j in Ja for c in C for m in M for d in D) +
        sum(CE[j,c]*((1+αe)^(c-1))*PEC[i]*L[j]*Te_exato[i,j,c,m,d] for i in Ve for j in Je for c in C for m in M for d in D)
    end
    return (prod_total,custo_total)
end

#Function que calcula os desvios
function Desvios(Ta_exato,Te_exato,Ta1,Te1,Ta2,Te2,Ta3,Te3,Ta4,Te4,Va,Ja,M,D,Ve,Je,C)

    if size(Ta1,1) > 0
        Des1=zeros(length(Ja)+length(Je))
        for j in Ja
            Des1[j] = sum(abs(d)*Ta1[i,j,1,m,d] for i in Va for m in M for d in D)
        end
        for j in Je
            Des1[j] = sum(abs(d)*Te1[i,j,1,m,d] for i in Ve for m in M for d in D)
        end
        Des2=zeros(length(Ja)+length(Je))
        for j in Ja
            Des2[j] = sum(abs(d)*Ta2[i,j,2,m,d] for i in Va for m in M for d in D)
        end
        for j in Je
            Des2[j] = sum(abs(d)*Te2[i,j,2,m,d] for i in Ve for m in M for d in D)
        end
        Des3=zeros(length(Ja)+length(Je))
        for j in Ja
            Des3[j] = sum(abs(d)*Ta3[i,j,3,m,d] for i in Va for m in M for d in D)
        end
        for j in Je
            Des3[j] = sum(abs(d)*Te3[i,j,3,m,d] for i in Ve for m in M for d in D)
        end
        Des4=zeros(length(Ja)+length(Je))
        for j in Ja
            Des4[j] = sum(abs(d)*Ta4[i,j,4,m,d] for i in Va for m in M for d in D)
        end
        for j in Je
            Des4[j] = sum(abs(d)*Te4[i,j,4,m,d] for i in Ve for m in M for d in D)
        end

        Des_total = Des1 + Des2 + Des3 + Des4
    end

    if size(Ta_exato,1) > 0
        Des=zeros(length(Ja)+length(Je))
        for j in Ja
            Des[j] = sum(abs(d)*Ta_exato[i,j,c,m,d] for i in Va for c in C for m in M for d in D)
        end
        for j in Je
            Des[j] = sum(abs(d)*Te_exato[i,j,c,m,d] for i in Ve for c in C for m in M for d in D)
        end

        Des_total = Des
    end

    return (Des_total)
end

#Função que mede a concentração das variedades mais produtivas na plantação
function Concentracao(L,PSC,PEC,Xa_exato,Xe_exato,Xa,Xe,Ja,H,Je)

    Var_A_mais_prod = sortperm(PSC,rev=true)[1:2]
    Var_E_mais_prod = 21

    if size(Xa,1) > 0
        PropAMais = 100*sum(L[j]*Xa[i,j,h] for i in Var_A_mais_prod for j in Ja for h in H)/sum(L[j] for j in Ja)
        PropEMais = 100*sum(L[j]*Xe[i,j,h] for i in Var_E_mais_prod for j in Je for h in H)/sum(L[j] for j in Je)
    end

    if size(Xa_exato,1) > 0
        PropAMais = 100*sum(L[j]*Xa_exato[i,j,h] for i in Var_A_mais_prod for j in Ja for h in H)/sum(L[j] for j in Ja)
        PropEMais = 100*sum(L[j]*Xe_exato[i,j,h] for i in Var_E_mais_prod for j in Je for h in H)/sum(L[j] for j in Je)
    end

    Var_A_menos_prod = sortperm(PSC[findall(PSC .> 0)])[1:2]
    Var_E_menos_prod = 25

    if size(Xa,1) > 0
        PropAMenos = 100*sum(L[j]*Xa[i,j,h] for i in Var_A_menos_prod for j in Ja for h in H)/sum(L[j] for j in Ja)
        PropEMenos = 100*sum(L[j]*Xe[i,j,h] for i in Var_E_menos_prod for j in Je for h in H)/sum(L[j] for j in Ja)
    end

    if size(Xa_exato,1) > 0
        PropAMenos = 100*sum(L[j]*Xa_exato[i,j,h] for i in Var_A_menos_prod for j in Ja for h in H)/sum(L[j] for j in Ja)
        PropEMenos = 100*sum(L[j]*Xe_exato[i,j,h] for i in Var_E_menos_prod for j in Je for h in H)/sum(L[j] for j in Je)
    end

    return (PropAMais,PropEMais,PropAMenos,PropEMenos)
end

#------------------------------------------------------------------------------------------
#Resolução Heurística - Máxima produção
#a1=1 e a2=0,para maximizar produção
a1 = 1
a2 = 0

(Xa1_MP,Xe1_MP,Ta1_MP,Te1_MP,Y1_MP,Tempo1_MP) = Plantio_Colheita_Primeiro_Corte(V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,1.0*DS,1.0*DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 2
(Ta2_MP,Te2_MP,Tempo2_MP) = Colheita_anos_seguintes(Ta1_MP,Te1_MP,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 3
(Ta3_MP,Te3_MP,Tempo3_MP) = Colheita_anos_seguintes(Ta2_MP,Te2_MP,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 4
(Ta4_MP,Te4_MP,Tempo4_MP) = Colheita_anos_seguintes(Ta3_MP,Te3_MP,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

FileIO.save("Xa_heuristico_1_MP.jld2","RH",Xa1_MP)
FileIO.save("Xe_heuristico_1_MP.jld2","RH",Xe1_MP)
FileIO.save("Y_heuristico_1_MP.jld2","RH",Y1_MP)
FileIO.save("Ta_heuristico_1_MP.jld2","RH",Ta1_MP)
FileIO.save("Ta_heuristico_2_MP.jld2","RH",Ta2_MP)
FileIO.save("Ta_heuristico_3_MP.jld2","RH",Ta3_MP)
FileIO.save("Ta_heuristico_4_MP.jld2","RH",Ta4_MP)
FileIO.save("Te_heuristico_1_MP.jld2","RH",Te1_MP)
FileIO.save("Te_heuristico_2_MP.jld2","RH",Te2_MP)
FileIO.save("Te_heuristico_3_MP.jld2","RH",Te3_MP)
FileIO.save("Te_heuristico_4_MP.jld2","RH",Te4_MP)
#------------------------------------------------------------------------------------------
#Resolução Heurística - Mínimo Custo
#a1=0 e a2=-1,para mínimo custo
a1 = 0
a2 = -1

(Xa1_MC,Xe1_MC,Ta1_MC,Te1_MC,Y1_MC,Tempo1_MC) = Plantio_Colheita_Primeiro_Corte(V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,1.1*DS,1.1*DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 2
(Ta2_MC,Te2_MC,Tempo2_MC) = Colheita_anos_seguintes(Ta1_MC,Te1_MC,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 3
(Ta3_MC,Te3_MC,Tempo3_MC) = Colheita_anos_seguintes(Ta2_MC,Te2_MC,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 4
(Ta4_MC,Te4_MC,Tempo4_MC) = Colheita_anos_seguintes(Ta3_MC,Te3_MC,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

FileIO.save("Xa_heuristico_1_MC.jld2","RH",Xa1_MC)
FileIO.save("Xe_heuristico_1_MC.jld2","RH",Xe1_MC)
FileIO.save("Y_heuristico_1_MC.jld2","RH",Y1_MC)
FileIO.save("Ta_heuristico_1_MC.jld2","RH",Ta1_MC)
FileIO.save("Ta_heuristico_2_MC.jld2","RH",Ta2_MC)
FileIO.save("Ta_heuristico_3_MC.jld2","RH",Ta3_MC)
FileIO.save("Ta_heuristico_4_MC.jld2","RH",Ta4_MC)
FileIO.save("Te_heuristico_1_MC.jld2","RH",Te1_MC)
FileIO.save("Te_heuristico_2_MC.jld2","RH",Te2_MC)
FileIO.save("Te_heuristico_3_MC.jld2","RH",Te3_MC)
FileIO.save("Te_heuristico_4_MC.jld2","RH",Te4_MC)
#------------------------------------------------------------------------------------------
#Rodamos o exato para determinar as soluções lexicográficas exatas partindo das soluções heurísticas
#Máxima Produção
a1 = 1
a2 = 0
(Xa_exato_MP,Xe_exato_MP,Ta_exato_MP,Te_exato_MP,Tempo_exato_MP) = Exato(Xa1_MP,Xe1_MP,Ta1_MP,Te1_MP,Ta2_MP,Te2_MP,Ta3_MP,Te3_MP,Ta4_MP,Te4_MP,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

#Mínimo Custo
a1 = 0
a2 = -1
(Xa_exato_MC,Xe_exato_MC,Ta_exato_MC,Te_exato_MC,Tempo_exato_MC) = Exato(Xa1_MC,Xe1_MC,Ta1_MC,Te1_MC,Ta2_MC,Te2_MC,Ta3_MC,Te3_MC,Ta4_MC,Te4_MC,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

#Calculamos as coordenadas do vetor ideal exato
(prod_exato_MP,custo_exato_MP) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,Ta_exato_MP,Te_exato_MP,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C)
(prod_exato_MC,custo_exato_MC) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,Ta_exato_MC,Te_exato_MC,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C)

z1max = prod_exato_MP
z1min = prod_exato_MC
z2max = custo_exato_MP
z2min = custo_exato_MC
#------------------------------------------------------------------------------------------
#Determinamos a solução compromisso aproximada tendo como referência as coordenadas do Ideal Exato
(XaC1,XeC1,TaC1,TeC1,YC1,TempoC1) = Plantio_Colheita_Primeiro_Corte_Compromisso(z1max,z1min,z2max,z2min,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 2
(TaC2,TeC2,TempoC2) = Colheita_anos_seguintes_Compromisso(z1max,z1min,z2max,z2min,TaC1,TeC1,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 3
(TaC3,TeC3,TempoC3) = Colheita_anos_seguintes_Compromisso(z1max,z1min,z2max,z2min,TaC2,TeC2,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

Colheita = 4
(TaC4,TeC4,TempoC4) = Colheita_anos_seguintes_Compromisso(z1max,z1min,z2max,z2min,TaC3,TeC3,Colheita,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)

#------------------------------------------------------------------------------------------
#Determinamos a solução compromisso exata tendo como referência as coordenadas do Ideal Exato e usando a solução compromisso aproximada
(Xa_exato_compromisso,Xe_exato_compromisso,Ta_exato_compromisso,Te_exato_compromisso,Tempo_exato_compromisso) = Exato_Compromisso(z1max,z1min,z2max,z2min,XaC1,XeC1,TaC1,TeC1,TaC2,TeC2,TaC3,TeC3,TaC4,TeC4,V,Va,Ve,J,Ja,Je,H,C,M,D,SPSH,SPEH,FPSH,FPEH,ka,ke,DS,DF,Cap,L,PSC,PEC,αa,αe,a1,a2,CPUtime)
#------------------------------------------------------------------------------------------
#Calculamos os valores objetivos das soluções heuristicas
(prod_heuristico_MP,custo_heuristico_MP) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,[],[],Ta1_MP,Te1_MP,Ta2_MP,Te2_MP,Ta3_MP,Te3_MP,Ta4_MP,Te4_MP,Va,Ja,M,D,Ve,Je,C)
(prod_heuristico_MC,custo_heuristico_MC) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,[],[],Ta1_MC,Te1_MC,Ta2_MC,Te2_MC,Ta3_MC,Te3_MC,Ta4_MC,Te4_MC,Va,Ja,M,D,Ve,Je,C)
(prod_heuristico_compromisso,custo_heuristico_compromisso) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,[],[],TaC1,TeC1,TaC2,TeC2,TaC3,TeC3,TaC4,TeC4,Va,Ja,M,D,Ve,Je,C)
#Calculamos os valores objetivos da solução exata
(prod_exato_compromisso,custo_exato_compromisso) = Valores_obj(SPSH,SPEH,FPSH,FPEH,CA,αa,PSC,L,CE,αe,PEC,Ta_exato_compromisso,Te_exato_compromisso,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C)

AV_heuristico_MP = mean(Desvios([],[],Ta1_MP,Te1_MP,Ta2_MP,Te2_MP,Ta3_MP,Te3_MP,Ta4_MP,Te4_MP,Va,Ja,M,D,Ve,Je,C))
AV_heuristico_MC = mean(Desvios([],[],Ta1_MC,Te1_MC,Ta2_MC,Te2_MC,Ta3_MC,Te3_MC,Ta4_MC,Te4_MC,Va,Ja,M,D,Ve,Je,C))
AV_heuristico_compromisso = mean(Desvios([],[],TaC1,TeC1,TaC2,TeC2,TaC3,TeC3,TaC4,TeC4,Va,Ja,M,D,Ve,Je,C))

AV_exato_MP = mean(Desvios(Ta_exato_MP,Te_exato_MP,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C))
AV_exato_MC = mean(Desvios(Ta_exato_MC,Te_exato_MC,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C))
AV_exato_compromisso = mean(Desvios(Ta_exato_compromisso,Te_exato_compromisso,[],[],[],[],[],[],[],[],Va,Ja,M,D,Ve,Je,C))


(PropAMais_heuristico_MP,PropEMais_heuristico_MP,PropAMenos_heuristico_MP,PropEMenos_heuristico_MP) = Concentracao(L,PSC,PEC,[],[],Xa1_MP,Xe1_MP,Ja,H,Je)
(PropAMais_heuristico_MC,PropEMais_heuristico_MC,PropAMenos_heuristico_MC,PropEMenos_heuristico_MC) = Concentracao(L,PSC,PEC,[],[],Xa1_MC,Xe1_MC,Ja,H,Je)
(PropAMais_heuristico_C,PropEMais_heuristico_C,PropAMenos_heuristico_C,PropEMenos_heuristico_C) = Concentracao(L,PSC,PEC,[],[],XaC1,XeC1,Ja,H,Je)

(PropAMais_exato_MP,PropEMais_exato_MP,PropAMenos_exato_MP,PropEMenos_exato_MP) = Concentracao(L,PSC,PEC,Xa_exato_MP,Xe_exato_MP,[],[],Ja,H,Je)
(PropAMais_exato_MC,PropEMais_exato_MC,PropAMenos_exato_MC,PropEMenos_exato_MC) = Concentracao(L,PSC,PEC,Xa_exato_MP,Xe_exato_MP,[],[],Ja,H,Je)
(PropAMais_exato_C,PropEMais_exato_C,PropAMenos_exato_C,PropEMenos_exato_C) = Concentracao(L,PSC,PEC,Xa_exato_compromisso,Xe_exato_compromisso,[],[],Ja,H,Je)
