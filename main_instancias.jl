####Función principal que se usa para generar un numero de instancias con datos aleatorios

# Se cargan todas las librerías
include("./Funciones/cargarLibrerias.jl")

# Se cargan las funciones
include("./Funciones/cargarFunciones.jl")

Logging.disable_logging(Logging.Error)

# Se inicializa el programa con diferentes test
# principalmente para cargar los solvers y resolver con mayor rapidez el caso pedido por el usuario
boot()

# Variable para salir del bucle
finPrograma = false
# En caso de que no sea fin de programa
while !finPrograma

    # Limpiza del terminal
    limpiarTerminal()

    # Se entra en un bucle para que el usuario seleccione el caso que se quiere estudiar
    casoEstudio, opfTipo, s = selectEstudio()

    # Limpiza del terminal
    limpiarTerminal()

    # Se extrae los datos del caso de estudio
    # Donde:
    #   datos[1] = datos de las líneas
    #   datos[2] = datos de los generadores
    #   datos[3] = datos de la demanda
    #   datos[4] = número de nodos
    #   datos[5] = número de líneas
    #   datos[6] = potencia base
    #   datos[7] = ruta al archivo .m del caso
    println("\nExtrayendo datos...")
    datos = extraerDatos(casoEstudio)
    println("Datos extraídos.")

    # Se generan n instancias aleatorias
    n_instancias = 50
    instancias = generar_instancias_aleatorias!(datos[3], datos[2], datos[1], n_instancias)

    # Almaceno resultados
    resultados = []

    # Almaceno variable para luego calcular el tiempo de resolución
    tiempo_resolucion = Float64[] 

    for i in 1:n_instancias
        println("\nGenerando y resolviendo instancia $i de $n_instancias...")

        # Obtengo los datos modificados para esta instancia
        datosNodo_mod, datosGenerador_mod, datosLinea_mod = instancias[i]

        # Sustituyo las columnas antiguas por las nuevas
        datosNodo_mod.PD .= datosNodo_mod.PDEMANDA
        datosNodo_mod.QD .= datosNodo_mod.QDEMANDA
        datosGenerador_mod.P_MAX .= datosGenerador_mod.PMAX
        datosGenerador_mod.P_MIN .= datosGenerador_mod.PMIN
        datosGenerador_mod.Q_MAX .= datosGenerador_mod.QMAX
        datosGenerador_mod.Q_MIN .= datosGenerador_mod.QMIN
        datosGenerador_mod.P_COSTE0 .= datosGenerador_mod.COSTE0
        datosGenerador_mod.P_COSTE1 .= datosGenerador_mod.COSTE1
        datosGenerador_mod.P_COSTE2 .= datosGenerador_mod.COSTE2
        datosLinea_mod.R .= datosLinea_mod.RESISTENCIA
        datosLinea_mod.X .= datosLinea_mod.INDUCTANCIA

        # Elimino la nueva columna que se forma en el Dataframe
        datosNodo_mod = select!(datosNodo_mod, Not(:PDEMANDA))
        datosNodo_mod = select!(datosNodo_mod, Not(:QDEMANDA))
        datosGenerador_mod = select!(datosGenerador_mod, Not(:PMAX))
        datosGenerador_mod = select!(datosGenerador_mod, Not(:PMIN))
        datosGenerador_mod = select!(datosGenerador_mod, Not(:QMAX))
        datosGenerador_mod = select!(datosGenerador_mod, Not(:QMIN))
        datosGenerador_mod = select!(datosGenerador_mod, Not(:COSTE0))
        datosGenerador_mod = select!(datosGenerador_mod, Not(:COSTE1))
        datosGenerador_mod = select!(datosGenerador_mod, Not(:COSTE2))
        datosLinea_mod = select!(datosLinea_mod, Not(:RESISTENCIA))
        datosLinea_mod = select!(datosLinea_mod, Not(:INDUCTANCIA))

        # Guardo el tiempo en el que se inicia la resolución del problema
        tiempo_inicio = Dates.now()

        # Resolver el problema para esta instancia según el tipo de OPF
        if opfTipo == "LP-OPF"
            m, solGen, solFlujos, solAngulos = LP_OPF(datosLinea_mod, datosGenerador_mod, datosNodo_mod, datos[4], datos[5], datos[6], s)
            push!(resultados, (m, solGen, solFlujos, solAngulos, solBinaria, datos[7], opfTipo, s, coste_total))

        elseif opfTipo == "AC-OPF"
            m, solGen, solFlujos, solAngulos, solBinaria, coste_total = AC_OPF(datosLinea_mod, datosGenerador_mod, datosNodo_mod, datos[4], datos[5], datos[6], s)
            push!(resultados, (m, solGen, solFlujos, solAngulos, solBinaria, datos[7], opfTipo, s, coste_total))
        else
            println("ERROR: Fallo en cargar el tipo de OPF")
            break
        end

        #Guardo el momento en el que finaliza la resolución del problema
        tiempo_fin = Dates.now()

        #Ahora calculo la duración de la resolución haciendo una resta 
        duracion = Dates.value(tiempo_fin - tiempo_inicio)/1000 # Función que realiza una resta entre variables de tipo "time" en milisegundos

        #Almaceno el tiempo de resolución
        push!(tiempo_resolucion, duracion)

        # Limpieza del terminal
        limpiarTerminal()

        #Imprimo los datos modificados de forma aleatoria
        println("Datos de nodos modificados (instancia $i):")
        show(datosNodo_mod, allrows = true)
            
        println()

        println("Datos de generadores modificados (instancia $i):")
        show(datosGenerador_mod, allrows = true)  

        println()

        println("Datos de líneas modificados (instancia $i):")
        show(datosLinea_mod, allrows = true)  

        println()

        #Genero los resultados 
        println("Generando resultados para instancia $i...")
        gestorResultados(resultados[i]...)
    end

    limpiarTerminal()

    #Cálculo de la desviación típica del tiempo de resolución de una instancia tipo
    sd = std(tiempo_resolucion)

    #Cálculo del tiempo medio de resolución de una instancia tipo
    media = mean(tiempo_resolucion)

    #Imprimo la desviación y la media
    println("\nTiempo medio de resolución de las instancias: $media segundos")
    println("Desviación típica del tiempo de resolución de las instancias: $sd segundos")

    # Generar los resultados de optimización para todas las instancias
    println("Instancias resueltas")
    
    # Preguntar al usuario si quiere continuar en el bucle para estudiar otro caso
    println("\nPulsa la tecla ENTER para continuar o cualquier otra entrada para salir.")
    if readline() == ""
        # Se mantiene la variable en falso para continuar en el bucle
        finPrograma = false
    else
        # Actualización de la variable para salir del bucle
        finPrograma = true
        exit()
    end
end