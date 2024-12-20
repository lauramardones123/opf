using SparseArrays
using LinearAlgebra

function calcularAdmitancias(datosLinea::DataFrame, nNodos::Int64, nLineas::Int64)
    # Crear matriz de incidencia usando SparseArrays
    A = sparse(datosLinea.F_BUS, 1:nLineas, 1, nNodos, nLineas) + 
        sparse(datosLinea.T_BUS, 1:nLineas, -1, nNodos, nLineas)
    
    # Calcular impedancias y admitancias serie
    Z = datosLinea.R .+ im .* datosLinea.X
    Y_serie = A * spdiagm(1 ./ Z) * A'
    
    # Calcular admitancias shunt
    y_sh = 0.5 .* (im .* datosLinea.BSh)
    Y_sh = spdiagm(diag(A * spdiagm(y_sh) * A'))
    
    # Matriz de admitancias total
    Y = Y_serie + Y_sh

    # Vector de admitancias serie y shunt por línea
    y_series = 1 ./ Z
    y_shunts = 0.5 .* (im .* datosLinea.BSh)

    return Y, y_series, y_shunts
end 