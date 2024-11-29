## Función que ejecuta el PSO 

include("./pso_binarias.jl")

function runPSO(fitFunc::Function, nDim::Int, nParticle:: Int, nInter:: Int, c::String)

    ## fitFunc: función de coste de generación
    ## nDim: número de generadores
    ## nParticle: número de soluciones exploradas por las partículas 
    ## nInter: número de iteraciones del PSO
    ## c: archivo de estudio

    # Crear el enjambre 
    s = Swarm(fitFunc, nDim, nParticle = nParticle, nInter = nInter, c = c) 

    # Inicializar los valores del fit y actualizar las soluciones 
    initSwarm(s, c) 

    # Actualizar las soluciones del enjambre nInter veces. 
    updateSwarm(s, nInter, c) 
    
    # Devuelve la solución global 
    return s.gBest, s.fitgBest
end