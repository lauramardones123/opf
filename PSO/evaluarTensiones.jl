# Función que coge las potencias activas y reactivas y los estados de los generadores del PSO
# y calcula las tensiones y las violaciones de tensiones

using SparseArrays
using LinearAlgebra

function evaluarTensiones(datosLinea::DataFrame, datosGenerador::DataFrame, datosNodo::DataFrame,
    nNodos::Int64, nLineas::Int64, bMVA::Float64, potencias_P::Vector{Float64}, 
    potencias_Q::Vector{Float64}, estados_u::Vector{Float64}, 
    Y::Matrix{Complex{Float64}}, 
    log_file::Union{IOStream, Nothing}, log_enabled::Bool)
    
    log_to_file(log_file, "\nEvaluando tensiones:", log_enabled)
    log_to_file(log_file, "Potencias P: $potencias_P", log_enabled) # Potencias activas obtenidas del PSO
    log_to_file(log_file, "Potencias Q: $potencias_Q", log_enabled) # Potencias reactivas obtenidas del PSO
    
    # Ajustar potencias según estados (solo para el cálculo, no modifica los valores originales)
    potencias_P_calc = copy(potencias_P)
    potencias_Q_calc = copy(potencias_Q)
    
    for i in 1:length(potencias_P)
        if estados_u[i] < 0.5  # Si el generador está apagado las potencias deben ser 0 para 
                               # cumplir con las restricción Pmin*x <= P <= Pmax*x
                               # si está encendido, se mantienen las potencias obtenidas del PSO
            potencias_P_calc[i] = 0.0 
            potencias_Q_calc[i] = 0.0
        end
    end

    # Inicializar vector de tensiones (flat start)
    V = ones(Complex{Float64}, nNodos)

    # Newton-Raphson para calcular tensiones
    max_iter = 100
    tol = 1e-6

    for iter in 1:max_iter
        # Calcular potencias inyectadas
        S_calc = zeros(Complex{Float64}, nNodos)
        for i in 1:nNodos
            # Restricción de la relación entre los nodos: S_Gen - S_Demand = V * conj(Y * V)
            S_calc[i] = V[i] * conj(sum(Y[i,j] * V[j] for j in 1:nNodos))
        end

        # Calcular desbalances
        ΔP = zeros(Float64, nNodos)
        ΔQ = zeros(Float64, nNodos)

        for i in 1:nNodos
            # Potencia generada (si hay generador activo en el nodo)
            Pgen = 0.0
            Qgen = 0.0
            for (idx, gen) in enumerate(datosGenerador.BUS)
                if gen == i
                    Pgen = potencias_P_calc[idx]  # Usar las potencias calculadas
                    Qgen = potencias_Q_calc[idx]
                end
            end

            # Potencia demandada
            Pdem = datosNodo.PD[i]
            Qdem = datosNodo.QD[i]

            # Desbalances considerando P y Q (La potencia inyectada debe ser igual a la generada
            # menos la demandada)
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
        Vmin = datosNodo.Vmin[i]
        Vmax = datosNodo.Vmax[i]
        Vmag = abs(V[i]) # Módulo de la tensión
        Vang = angle(V[i]) * (180/π)  # Convertir a grados

        println("Vmin: ", Vmin, " Vmax: ", Vmax, " Vmag: ", Vmag)
        log_to_file(log_file, "\nNodo $i:", log_enabled)
        log_to_file(log_file, "Magnitud de tensión: $(round(Vmag, digits=4)) p.u.", log_enabled)
        log_to_file(log_file, "Ángulo de tensión: $(round(Vang, digits=4))°", log_enabled)

        if Vmag < Vmin
            #violaciones += (Vmin - Vmag)^2
            violaciones = Inf
        elseif Vmag > Vmax
            # violaciones += (Vmag - Vmax)^2
            violaciones = Inf
        end
    end
    println("violaciones: ", violaciones)
    println("V: ", V)

    log_to_file(log_file, "\nViolaciones totales: $violaciones", log_enabled)
    log_to_file(log_file, "\nTensiones completas (p.u.):", log_enabled)
    for i in 1:nNodos
        log_to_file(log_file, "V[$i] = $(round(abs(V[i]), digits=4))∠$(round(angle(V[i])*(180/π), digits=4))°", log_enabled)
    end

    # Mostrar información de potencias y límites
    log_to_file(log_file, "\nInformación de Generadores:", log_enabled)
    log_to_file(log_file, "----------------------------", log_enabled)
    for i in 1:length(potencias_P)
        log_to_file(log_file, "\nGenerador $i:", log_enabled)
        log_to_file(log_file, "Estado u: $(round(estados_u[i], digits=3))", log_enabled)
        log_to_file(log_file, "P inicial: $(round(potencias_P[i], digits=2)) MW", log_enabled)
        log_to_file(log_file, "Q inicial: $(round(potencias_Q[i], digits=2)) MVAr", log_enabled)
        log_to_file(log_file, "Límites P: [$(round(datosGenerador.P_MIN[i], digits=2)), $(round(datosGenerador.P_MAX[i], digits=2))] MW", log_enabled)
        log_to_file(log_file, "Límites Q: [$(round(datosGenerador.Q_MIN[i], digits=2)), $(round(datosGenerador.Q_MAX[i], digits=2))] MVAr", log_enabled)
        log_to_file(log_file, "P ajustada: $(round(potencias_P_calc[i], digits=2)) MW", log_enabled)
        log_to_file(log_file, "Q ajustada: $(round(potencias_Q_calc[i], digits=2)) MVAr", log_enabled)
        log_to_file(log_file, "Estado: $(estados_u[i] >= 0.5 ? "Encendido" : "Apagado")", log_enabled)
    end
    println("\n")

    return V, violaciones
end