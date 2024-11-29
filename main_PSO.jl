## Función principal, se ejecuta en la terminal. Realiza el PSO y, después ejecuta el AC OPF en caso de 
## que el usuario quiera resolverlo 

# Se cargan todas las librerías
include("./PSO/cargarLibrerias_PSO.jl")

# Se cargan todas las funciones relacionadas con el PSO
include("./PSO/cargarFunciones_PSO.jl")

# Se desabilitan los mensajes de warning o info del logging
Logging.disable_logging(Logging.Error)

"""
# Se inicializa el programa con diferentes test
# principalmente para cargar los solvers y resolver con mayor rapidez el caso pedido por el usuario
boot()
"""

# Variable para salir del bucle
finPrograma = false

# En caso de que no sea fin de programa
while !finPrograma
    
    # Vamos a generar el PSO seleccionando el caso de estudio y las características 
    println("\nSeleccionando características del PSO...")

    nDim, nParticle, nInter, casoEstudio, s = selectCaracteristicas()
    
    println("\nExtrayendo datos...")
    datos = extraerDatos_PSO(casoEstudio)
    println("\nDatos extraídos.")

    println("\nGenerando PSO...")

    # Vector que almacena la mejor combinación de apertura/cierre de generadores
    xs = Array{Float}(undef, 1) 

    # Variable que almacena el coste que supone la mejor combinación
    ys = zeros(1) 

    # Resolución del PSO
    xs, ys = runPSO(fitFunc, nDim, nParticle, nInter, casoEstudio)

    limpiarTerminal() 

    # Imprimir la solución en la terminal
    println("Mejor solución global:")
    println("Posición: ", xs)
    println("Coste: ", ys)

    # Preguntar al usuario si quiere resolver el AC OPF 
    println("\nPulsa la tecla ENTER para realizar el AC OPF con esta solución y cualquier otra para salir")
    
    if readline() == ""
        # Ejecución del AC OPF
        limpiarTerminal()
        println("\nGenerando OPF...")
        m, solGen, solFlujos, solAngulos, solBinaria, coste_total = AC_OPF_PSO(datos[1], datos[2], datos[3], datos[4], datos[5], datos[6], s, xs)
        limpiarTerminal()
        println("\nProblema resuelto.")
        gestorResultados_PSO(m, solGen, solFlujos, solAngulos, solBinaria, datos[7], s, coste_total)    
    else
        # Actualización de la variable para salir del bucle
        finPrograma = true
        exit()
    end

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