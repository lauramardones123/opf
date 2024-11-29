#### Una vez que se selecciona el estudio con la función "selectEstudio.jl", esta función guarda los datos .csv del
#### caso seleccionado en DataFrames*, calcula el número de nodos del sistema con estos datos y entraga la potencia
#### base del sistema y la ruta .m que necesita para optimizar posteriormente en PowerModels (comprueba el código, no se si lo haré yo)

#### *Un DataFrame almacena los valores del .csv en tablas

function extraerDatos(c::String)
    # Datos de las lineas
    datosLinea = CSV.read("Casos/$c/datosLineas.csv", DataFrame)

    # Datos de los generadores
    datosGenerador = CSV.read("Casos/$c/datosGeneradores.csv", DataFrame)

    # Datos de la demanda
    datosNodo = CSV.read("Casos/$c/datosNodos.csv", DataFrame)

    # Número de nodos
    nNodos = maximum([datosLinea.F_BUS; datosLinea.T_BUS])

    # Número de líneas
    nLineas = size(datosLinea, 1)

    # Potencia base
    bMVA = 100

    # Ruta al archivo .m
    rutaArchivoM = "Casos/$c/$c.m"

    if isfile(rutaArchivoM)
        ruta = rutaArchivoM
    else
        ruta = "None"
    end

    # Devuelve todos los DataFrames y variables generadas
    return(datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA, ruta)
end