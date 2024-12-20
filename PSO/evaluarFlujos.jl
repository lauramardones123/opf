function evaluarFlujos(datosLinea::DataFrame, Y::Matrix{Complex{Float64}}, 
                      V::Vector{Complex{Float64}}, bMVA::Float64,
                      y_series::Vector{Complex{Float64}},
                      y_shunts::Vector{Complex{Float64}},
                      log_file::Union{IOStream, Nothing}=nothing, log_enabled::Bool=false)
    
    nLineas = size(datosLinea, 1)
    violaciones_flujo = 0.0
    
    log_to_file(log_file, "\nEvaluando flujos de potencia:", log_enabled)
    
    # Calcular Y_0 y Y_Sh para cada línea
    for k in 1:nLineas
        i = datosLinea.F_BUS[k]
        j = datosLinea.T_BUS[k]
        
        # Usar los valores precalculados
        y_serie = y_series[k]
        y_sh = y_shunts[k]
        
        # Calcular corrientes usando la formulación dada
        I_ij = V[i] * y_sh + (V[i] - V[j]) * y_serie
        I_ji = V[j] * y_sh + (V[j] - V[i]) * y_serie
        
        # Calcular flujos de potencia aparente
        S_ij = V[i] * conj(I_ij)
        S_ji = V[j] * conj(I_ji)
        
        # Calcular cuadrados de las magnitudes de los flujos
        S_ij_mag_sq = real(S_ij)^2 + imag(S_ij)^2
        S_ji_mag_sq = real(S_ji)^2 + imag(S_ji)^2
        
        # Límite de potencia aparente normalizado al cuadrado
        S_max_sq = (datosLinea.L_SMAX[k] / bMVA)^2
        
        # Calcular violaciones (similar a la restricción dada)
        if S_ij_mag_sq > S_max_sq
            violaciones_flujo += Inf
        end
        if S_ij_mag_sq < -S_max_sq
            violaciones_flujo += Inf
        end
        
        # Logging
        if log_enabled
            log_to_file(log_file, "\nLínea $k ($i->$j):", log_enabled)
            log_to_file(log_file, "Flujo S_ij: $(round(sqrt(S_ij_mag_sq) * bMVA, digits=2)) MVA", log_enabled)
            log_to_file(log_file, "Flujo S_ji: $(round(sqrt(S_ji_mag_sq) * bMVA, digits=2)) MVA", log_enabled)
            log_to_file(log_file, "Límite: $(round(sqrt(S_max_sq) * bMVA, digits=2)) MVA", log_enabled)
        end
    end
    
    return violaciones_flujo
end