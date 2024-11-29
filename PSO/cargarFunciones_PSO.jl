## Se cargan todas las funciones que necesita el main_PSO 

include("./pso_binarias.jl")                    # Script que recoge todas las funciones necesarias para resolver el PSO
include("./selectCaracteristicas.jl")           # Selección del número de soluciones a analizar en el PSO 
include("./runPSO.jl")                          # Ejecución del PSO 

include("../Funciones/elegirOpcion.jl")         # Elegir el caso a estudiar 
include("./selectEstudio_PSO.jl")               # Selección del archivo que se quiere estudiar
include("./extraerDatos_PSO.jl")                # Extracción de los datos del archivo seleccionado

include("../OPF/AC_OPF/AC_OPF_PSO.jl")          # Resolución del AC OPF con las binarias devueltas por el PSO
include("./gestorResultados_PSO.jl")            # Imprimir los resultados del AC OPF 

include("../Funciones/limpiarTerminal.jl")      # Limpiar la terminal

include("./pso_hibrido.jl")                     # Cargar todas las funciones relacionadas con el PSO híbrido
include("./evaluarTensiones.jl")                # Cargar funciones de evaluación de tensiones




