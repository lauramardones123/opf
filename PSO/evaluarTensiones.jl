function evaluarTensiones(datosLinea::DataFrame, datosGenerador::DataFrame, datosNodo::DataFrame,
    nNodos::Int64, nLineas::Int64, bMVA::Float64, generadores_activos::Vector{Float64}, 
    potencias_P::Vector{Float64}, potencias_Q::Vector{Float64})
    
    println("Construir la matriz de admitancias")
    Y = zeros(Complex{Float64}, nNodos, nNodos)    

    # Llenar la matriz Y con los datos de las líneas
    for i in 1:nLineas
        from_bus = datosLinea.F_BUS[i]
        to_bus = datosLinea.T_BUS[i]

        # Usar los nombres correctos de las columnas
        r = datosLinea.R[i]  # Cambio de BR_R a R
        x = datosLinea.X[i]  # Cambio de BR_X a X
        b = datosLinea.BSh[i]  # Cambio de BR_B a BSh

        y = 1/(r + im*x)

        Y[from_bus,from_bus] += y + im*b/2
        Y[to_bus,to_bus] += y + im*b/2

        Y[from_bus,to_bus] -= y
        Y[to_bus,from_bus] -= y
    end
    println("matriz Y rellenada")

    # Inicializar vector de tensiones (flat start)
    V = ones(Complex{Float64}, nNodos)

    # Newton-Raphson para calcular tensiones
    max_iter = 100
    tol = 1e-6

    for iter in 1:max_iter
        # Calcular potencias inyectadas
        S_calc = zeros(Complex{Float64}, nNodos)
        for i in 1:nNodos
            S_calc[i] = V[i] * conj(sum(Y[i,j] * V[j] for j in 1:nNodos))
        end

        # Calcular desbalances
        ΔP = zeros(Float64, nNodos)
        ΔQ = zeros(Float64, nNodos)

        for i in 1:nNodos
            # Potencia generada (si hay generador activo en el nodo)
            Pgen = 0.0
            Qgen = 0.0
            for (idx, gen) in enumerate(datosGenerador.BUS)  # Cambio de GEN_BUS a BUS
                if gen == i && generadores_activos[idx] == 1
                    Pgen = potencias_P[idx]
                    Qgen = potencias_Q[idx]
                end
            end

            # Potencia demandada
            Pdem = datosNodo.PD[i]
            Qdem = datosNodo.QD[i]

            # Desbalances considerando P y Q
            ΔP[i] = (Pgen - Pdem)/bMVA - real(S_calc[i])
            ΔQ[i] = (Qgen - Qdem)/bMVA - imag(S_calc[i])
        end

        # Verificar convergencia
        if maximum(abs.(vcat(ΔP, ΔQ))) < tol
            break
        end

        # Actualizar tensiones usando Newton-Raphson
        for i in 1:nNodos
            V[i] += 0.1 * (ΔP[i] + im*ΔQ[i])
            V[i] = V[i]/abs(V[i])  # Normalizar magnitud
        end
    end

    # Calcular violaciones de límites de tensión
    violaciones = 0.0
    for i in 1:nNodos
        Vmin = datosNodo.Vmin[i]  # Cambio de VMIN a Vmin
        Vmax = datosNodo.Vmax[i]  # Cambio de VMAX a Vmax
        Vmag = abs(V[i])

        if Vmag < Vmin
            violaciones += (Vmin - Vmag)^2
        elseif Vmag > Vmax
            violaciones += (Vmag - Vmax)^2
        end
    end
    println("violaciones: ", violaciones)
    println("V: ", V)

    return abs.(V), violaciones
end