### Función que genera valores aleatorios de una misma instancia (red con la misma topología)

function generar_instancias_aleatorias!(datosNodo, datosGenerador, datosLinea, n_instancias) #La exclamación indica que se van a modificar los datos del argumento de la función
    
    #Límites para los parámetros de nodos

    min_pd = 0.0   # Límite inferior de la demanda de activa
    max_pd = 150.0  # Límite superior de la demanda de activa

    min_qd = 0.0   # Límite inferior de la demanda de reactiva
    max_qd = 1.0  # Límite superior de la demanda de reactiva

    #Límites para los parámetros de generadores

    min_pmax = 10.0   # Límite inferior de la potencia máxima
    max_pmax = 1500.0  # Límite superior de la potencia máxima

    min_pmin = 0.0   # Límite inferior de la potencia mínima
    max_pmin = 100.0  # Límite superior de la potencia mínima

    min_qmax = 15.0   # Límite inferior de la potencia máxima
    max_qmax = 30.0  # Límite superior de la potencia máxima

    min_qmin = -100.0   # Límite inferior de la potencia máxima
    max_qmin = -10.0  # Límite superior de la potencia máxima

    min_coste0 = 0.0 # Límite inferior del coste 0
    max_coste0 = 10.0 #Límite superior del coste 0

    min_coste1 = 20.0 # Límite inferior del coste 1
    max_coste1 = 70.0 #Límite superior del coste 1

    min_coste2 = 0.0 # Límite inferior del coste 2
    max_coste2 = 1.0 #Límite superior del coste 2

    #Límites para los parámetros de las líneas

    min_resistencia = 0 #Límte inferior de la resistencia
    max_resistencia = 0.5 #Límite superior de la resistencia 

    min_inductancia = 0 #Límte inferior de la resistencia
    max_inductancia = 1 #Límite superior de la resistencia 


    instancias = []

    for i in 1:n_instancias
        # Copia los datos para no modificar los originales
        copia_datosNodo = copy(datosNodo)
        copia_datosGenerador = copy(datosGenerador)
        copia_datosLinea = copy(datosLinea)

        # Modificar demanda de forma aleatoria
        copia_datosNodo.PDEMANDA .= min_pd .+ (max_pd - min_pd) .* rand(size(copia_datosNodo, 1))
        copia_datosNodo.QDEMANDA .= min_qd .+ (max_qd - min_qd) .* rand(size(copia_datosNodo, 1))

        #Modificar valor de impedancia de forma aleatoria. 
        copia_datosLinea.RESISTENCIA .= min_resistencia .+ (max_resistencia - min_resistencia) .* rand(size(copia_datosLinea, 1))
        copia_datosLinea.INDUCTANCIA .= min_inductancia .+ (max_inductancia - min_inductancia) .* rand(size(copia_datosLinea, 1))


        # Modificar capacidad de generadores de forma aleatoria
        copia_datosGenerador.PMAX .= min_pmax .+ (max_pmax - min_pmax) .* rand(size(copia_datosGenerador, 1))
        copia_datosGenerador.PMIN .= min_pmin .+ (max_pmin - min_pmin) .* rand(size(copia_datosGenerador, 1))
        copia_datosGenerador.QMAX .= min_qmax .+ (max_qmax - min_qmax) .* rand(size(copia_datosGenerador, 1))
        copia_datosGenerador.QMIN .= min_qmin .+ (max_qmin - min_qmin) .* rand(size(copia_datosGenerador, 1))
    
        #Modificar coste generación de forma aleatoria
        copia_datosGenerador.COSTE0 .= min_coste0 .+ (max_coste0 - min_coste0) .* rand(size(copia_datosGenerador, 1))
        copia_datosGenerador.COSTE1 .= min_coste1 .+ (max_coste1 - min_coste1) .* rand(size(copia_datosGenerador, 1))
        copia_datosGenerador.COSTE2 .= min_coste2 .+ (max_coste2 - min_coste2) .* rand(size(copia_datosGenerador, 1))
        
        push!(instancias, (copia_datosNodo, copia_datosGenerador, copia_datosLinea))
    end

    return instancias
end