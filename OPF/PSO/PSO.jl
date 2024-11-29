function RunPSO(fitFunc::Function, nDim::Int, nParticle:: Int, nInter:: Int)

    s = Swarm(fitFunc, nDim, nParticle = nParticle, nInter = nInter) # Creamos el enjambre 
    initSwarm(s)  # Inicializar los valores del fit y actualizar las soluciones 

    for i in s.nInter
        updateSwarm(s) # Vamos a actualizar las soluciones del enjambre n interacciones. 
    end

    # Imprimir la solución global 
    s.gBest
end