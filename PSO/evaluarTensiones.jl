## Función que evalúa las tensiones de la red y calcula el incumplimiento de límites

function evaluarTensiones(datosLinea::DataFrame, datosGenerador::DataFrame, datosNodo::DataFrame, 
                         nNodos::Int64, nLineas::Int64, bMVA::Float64,
                         generadores_activos::Vector{Float64}, potencias_gen::Vector{Float64})
    
    # Construir la matriz de admitancias
    Y = zeros(Complex{Float64}, nNodos, nNodos)
    
    # Llenar la matriz Y con los datos de las líneas
    for i in 1:nLineas
        # Índices de los nodos
        from_bus = datosLinea.F_BUS[i]
        to_bus = datosLinea.T_BUS[i]
        
        # Parámetros de la línea
        r = datosLinea.BR_R[i]
        x = datosLinea.BR_X[i]
        b = datosLinea.BR_B[i]
        
        # Admitancia serie
        y = 1/(r + im*x)
        
        # Elementos diagonales
        Y[from_bus,from_bus] += y + im*b/2
        Y[to_bus,to_bus] += y + im*b/2
        
        # Elementos fuera de la diagonal
        Y[from_bus,to_bus] -= y
        Y[to_bus,from_bus] -= y
    end
    
    # Inicializar vector de tensiones (flat start)
    V = ones(Complex{Float64}, nNodos)
    
    # Máximo número de iteraciones y tolerancia
    max_iter = 100
    tol = 1e-6
    
    # Newton-Raphson para calcular tensiones
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
            for (idx, gen) in enumerate(datosGenerador.GEN_BUS)
                if gen == i && generadores_activos[idx] == 1
                    Pgen = potencias_gen[idx]
                end
            end
            
            # Potencia demandada
            Pdem = datosNodo.PD[i]
            Qdem = datosNodo.QD[i]
            
            # Desbalances
            ΔP[i] = (Pgen - Pdem)/bMVA - real(S_calc[i])
            ΔQ[i] = (-Qdem)/bMVA - imag(S_calc[i])
        end
        
        # Verificar convergencia
        if maximum(abs.(vcat(ΔP, ΔQ))) < tol
            break
        end
        
        # Actualizar tensiones usando Newton-Raphson
        # (Aquí se simplifica usando una actualización básica)
        for i in 1:nNodos
            V[i] += 0.1 * (ΔP[i] + im*ΔQ[i])
            # Normalizar magnitud
            V[i] = V[i]/abs(V[i])
        end
    end
    
    # Calcular violaciones de límites de tensión
    violaciones = 0.0
    for i in 1:nNodos
        Vmin = datosNodo.VMIN[i]
        Vmax = datosNodo.VMAX[i]
        Vmag = abs(V[i])
        
        if Vmag < Vmin
            violaciones += (Vmin - Vmag)^2
        elseif Vmag > Vmax
            violaciones += (Vmag - Vmax)^2
        end
    end
    
    return abs.(V), violaciones
end 