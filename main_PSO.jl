# Carga de librerías y funciones
include("./PSO/cargarLibrerias_PSO.jl")
include("./PSO/evaluarTensiones.jl")
include("./PSO/cargarFunciones_PSO.jl")
Logging.disable_logging(Logging.Error)

"""
Función principal que ejecuta el PSO y opcionalmente el AC OPF
"""
function ejecutar_optimizacion(caso_estudio::String, parametros::Dict)
    # Extraer parámetros
    n_particulas = parametros["n_particulas"]
    n_iteraciones = parametros["n_iteraciones"]
    
    println("\nExtrayendo datos...")
    datos = extraerDatos_PSO(caso_estudio)
    println("Datos extraídos.")

    if parametros["tipo_pso"] == "binario"
        println("\nGenerando PSO binario...")    
        n_dim = size(datos[2], 1)
        println("datos[2]: ", datos[2])
        mejor_solucion, mejor_coste = runPSO(fitFunc, n_dim, n_particulas, n_iteraciones, caso_estudio)
    else            
        println("\nGenerando PSO híbrido...")
        n_dim = size(datos[2], 1)
        println("n_dim: ", n_dim)
        mejor_solucion, potencias, mejor_coste = runPSOHibrido(datos, n_particulas, n_iteraciones)
    end

    # if parametros["ejecutar_ac_opf"] == 1
    #     # println("\nGenerando OPF...")
    #     m, solGen, solFlujos, solAngulos, solBinaria, coste_total = AC_OPF_PSO(
    #         datos[1], datos[2], datos[3], datos[4], datos[5], datos[6], 
    #         caso_estudio, mejor_solucion
    #     )
    #     # println("\nProblema resuelto.")
    #     gestorResultados_PSO(m, solGen, solFlujos, solAngulos, solBinaria, datos[7], caso_estudio, coste_total)
    # end

    return mejor_solucion, mejor_coste
end

# Configuración y ejecución principal
if abspath(PROGRAM_FILE) == @__FILE__
    println("Ejecutando main_PSO.jl")

    # Parámetros configurables
    parametros = Dict(
        "caso_estudio" => "EjemploTwitter_kyrib",  # Caso de estudio a resolver
        "tipo_pso" => "hibrido",                   # "binario" o "hibrido"
        "n_particulas" => 20,                      # Número de partículas
        "n_iteraciones" => 3,                   # Número de iteraciones
        "ejecutar_ac_opf" => 0                     # 0 para false, 1 para true
    )

    # Casos de estudio disponibles
    casos_disponibles = [
        "EjemploTwitter_kyrib",
        "pglib_opf_case3_lmbd",
        "pglib_opf_case5_pjm",
        "pglib_opf_case14_ieee",
        "pglib_opf_case30_ieee",
        "pglib_opf_case118_ieee",
        "pglib_opf_case300_ieee",
        "pglib_opf_case1354_pegase"
    ]

     # try
        mejor_solucion, mejor_coste = ejecutar_optimizacion(
            parametros["caso_estudio"],  # TO DO: simplificar y que el argumento sea solo parametros
            parametros
        )
        
        # Guardar resultados si se desea
        # if parametros["guardar_resultados"]
        #     # Aquí se puede añadir código para guardar los resultados
        #     println("\nResultados guardados.")
        # end

    # catch e
        # println("\nError durante la ejecución:")
        # println(e)
    # end
end