const Float = Float64
using CSV
using DataFrames

mutable struct ParticleHibrida
    nGeneradores::Int
    
    # Posiciones binarias (on/off de generadores)
    position_bin::Array{Float, 1}
    velocity_bin::Array{Float, 1}
    pBest_bin::Array{Float, 1}
    lBest_bin::Array{Float, 1}
    
    # Posiciones continuas (potencias de generadores)
    position_pot::Array{Float, 2}  # Matriz nx2 (P y Q para cada generador)
    velocity_pot::Array{Float, 2}
    pBest_pot::Array{Float, 2}
    lBest_pot::Array{Float, 2}
    
    # Valores de fitness
    fitValue::Float
    fitpBest::Float
    fitlBest::Float
    
    nFitEval::Int
    
    function ParticleHibrida(nGeneradores::Int, datosGenerador::DataFrame) 
        # Inicialización binaria
        position_bin = rand(0:1.0, nGeneradores)
        velocity_bin = rand(nGeneradores) - position_bin
        pBest_bin = copy(position_bin)
        lBest_bin = copy(position_bin)
        
        # Inicialización de potencias
        position_pot = zeros(Float, nGeneradores, 2)
        velocity_pot = zeros(Float, nGeneradores, 2)
        
        # Inicializar potencias dentro de límites
        for i in 1:nGeneradores
            # Potencia activa
            Pmin = datosGenerador.PMIN[i]
            Pmax = datosGenerador.PMAX[i]
            position_pot[i,1] = Pmin + rand()*(Pmax - Pmin)
            
            # Potencia reactiva
            Qmin = datosGenerador.QMIN[i]
            Qmax = datosGenerador.QMAX[i]
            position_pot[i,2] = Qmin + rand()*(Qmax - Qmin)
        end
        
        pBest_pot = copy(position_pot)
        lBest_pot = copy(position_pot)
        
        fitValue = Inf
        fitpBest = Inf
        fitlBest = Inf
        
        nFitEval = 0
        
        new(nGeneradores, position_bin, velocity_bin, pBest_bin, lBest_bin,
            position_pot, velocity_pot, pBest_pot, lBest_pot,
            fitValue, fitpBest, fitlBest, nFitEval)
    end       
end

mutable struct SwarmHibrido
    fitFunc::Function
    nGeneradores::Int
    datos::Tuple
    
    nParticle::Int
    nNeibor::Int
    nInter::Int
    
    c1::Float
    c2::Float
    
    wMax::Float
    wMin::Float
    w::Float
    
    gBest_bin::Array{Float, 1}    
    gBest_pot::Array{Float, 2}
    fitgBest::Float
    
    particles::Array{ParticleHibrida, 1}
    
    nFitEvals::Int
    
    function SwarmHibrido(fitFunc::Function, nGeneradores::Int, datos::Tuple;
            nParticle::Int=3, nNeibor::Int=3, nInter::Int=2000,
            c1::Float=2.0, c2::Float=2.0,
            wMax::Float=0.9, wMin::Float=0.4)
        
        if nNeibor > nParticle
            error("El número de partículas en un grupo local no debe exceder el número total de partículas")
        end
        
        w = wMax
        datosGenerador = datos[2]
        
        # Inicializar partículas
        particles = [ParticleHibrida(nGeneradores, datosGenerador) for i in 1:nParticle]
        
        # Inicializar mejores globales
        gBest_bin = rand(0:1.0, nGeneradores)
        gBest_pot = zeros(Float, nGeneradores, 2)
        fitgBest = Inf
        
        nFitEvals = 0
        
        new(fitFunc, nGeneradores, datos, nParticle, nNeibor, nInter,
            c1, c2, wMax, wMin, w, gBest_bin, gBest_pot, fitgBest,
            particles, nFitEvals)
    end
end

function evaluarParticula(p::ParticleHibrida, datos::Tuple)
    datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA = datos
    
    # Extraer potencias activas y reactivas
    potencias_P = p.position_pot[:,1]
    potencias_Q = p.position_pot[:,2]
    
    # Evaluar tensiones y violaciones
    V_mag, violaciones = evaluarTensiones(datosLinea, datosGenerador, datosNodo, 
                                        nNodos, nLineas, bMVA,
                                        p.position_bin, potencias_P)
    
    # Calcular coste de generación
    coste = 0.0
    for i in 1:p.nGeneradores
        if p.position_bin[i] == 1
            coste += datosGenerador.P_COSTE1[i] * potencias_P[i]
        end
    end
    
    # Penalización por violaciones de tensión
    coste_total = coste + 1000 * violaciones
    
    return coste_total, V_mag
end

function updateParticle!(p::ParticleHibrida, w::Float, c1::Float, c2::Float, datos::Tuple)
    # Actualizar velocidades
    p.velocity_bin = w * p.velocity_bin + 
                    c1 * rand() * (p.pBest_bin - p.position_bin) + 
                    c2 * rand() * (p.lBest_bin - p.position_bin)
    
    p.velocity_pot = w * p.velocity_pot + 
                    c1 * rand() * (p.pBest_pot - p.position_pot) + 
                    c2 * rand() * (p.lBest_pot - p.position_pot)
    
    # Actualizar posiciones
    p.position_bin += p.velocity_bin
    p.position_pot += p.velocity_pot
    
    # Discretizar binarias
    for i in 1:p.nGeneradores
        p.position_bin[i] = p.position_bin[i] >= 0.5 ? 1.0 : 0.0
    end
    
    # Limitar potencias
    datosGenerador = datos[2]
    for i in 1:p.nGeneradores
        p.position_pot[i,1] = clamp(p.position_pot[i,1], 
                                  datosGenerador.PMIN[i], 
                                  datosGenerador.PMAX[i])
        p.position_pot[i,2] = clamp(p.position_pot[i,2], 
                                  datosGenerador.QMIN[i], 
                                  datosGenerador.QMAX[i])
    end
    
    # Evaluar nueva posición
    p.fitValue, _ = evaluarParticula(p, datos)
    p.nFitEval += 1
    
    # Actualizar mejor personal si corresponde
    if p.fitValue < p.fitpBest
        p.fitpBest = p.fitValue
        p.pBest_bin = copy(p.position_bin)
        p.pBest_pot = copy(p.position_pot)
    end
end

function initFitValue!(fitFunc::Function, p::ParticleHibrida, datos::Tuple)
    p.fitValue, _ = fitFunc(p, datos)
    p.nFitEval += 1
    nothing
end

function updatepBestAndFitpBest!(p::ParticleHibrida)
    if p.fitValue < p.fitpBest
        p.fitpBest = p.fitValue
        p.pBest_bin = copy(p.position_bin)
        p.pBest_pot = copy(p.position_pot)
    end
    nothing
end

function updatelBestAndFitlBest!(s::SwarmHibrido)
    for i in 1:s.nParticle
        neiborIds = neiborIndices(i, s.nNeibor, s.nParticle)
        neiborFits = [s.particles[Id].fitValue for Id in neiborIds]
        fitlBest, index = findmin(neiborFits)
        
        if fitlBest < s.particles[i].fitlBest
            s.particles[i].lBest_bin = copy(s.particles[neiborIds[index]].position_bin)
            s.particles[i].lBest_pot = copy(s.particles[neiborIds[index]].position_pot)
            s.particles[i].fitlBest = fitlBest
        end
    end
    nothing
end

function updategBestAndFitgBest!(s::SwarmHibrido)
    gFits = [particle.fitValue for particle in s.particles]
    fitgBest, index = findmin(gFits)
    
    if fitgBest < s.fitgBest
        s.gBest_bin = copy(s.particles[index].position_bin)
        s.gBest_pot = copy(s.particles[index].position_pot)
        s.fitgBest = fitgBest
    end
    nothing
end

function neiborIndices(i::Int, nNeibor::Int, nParticle::Int)
    nNeibor = max(3, nNeibor)
    nLeft = (nNeibor - 1) ÷ 2
    startIndex = (i - nLeft)
    endIndex = startIndex + nNeibor - 1
    indices = collect(startIndex:endIndex)
    
    for i in 1:nNeibor
        if indices[i] < 1
            indices[i] += nParticle
        elseif indices[i] > nParticle
            indices[i] -= nParticle
        end
    end
    
    indices
end

function initSwarm(s::SwarmHibrido)
    for particle in s.particles
        initFitValue!(s.fitFunc, particle, s.datos)
        updatepBestAndFitpBest!(particle)
    end
    
    updatelBestAndFitlBest!(s)
    updategBestAndFitgBest!(s)
    
    nothing
end

function updateInertia!(s::SwarmHibrido)
    dw = (s.wMax - s.wMin)/s.nInter
    s.w -= dw
    nothing
end

function updateSwarm(s::SwarmHibrido)
    for i in 1:s.nInter
        # Actualizar cada partícula
        for p in s.particles
            updateParticle!(p, s.w, s.c1, s.c2, s.datos)
        end
        
        # Actualizar mejores locales y globales
        updatelBestAndFitlBest!(s)
        updategBestAndFitgBest!(s)
        updateInertia!(s)
        
        # Opcional: mostrar progreso
        if i % 100 == 0
            println("Iteración $i: Mejor fitness = $(s.fitgBest)")
        end
    end
    nothing
end

# Modificar runPSOHibrido para usar las nuevas funciones
function runPSOHibrido(datos::Tuple, nParticle::Int, nInter::Int)
    nGeneradores = size(datos[2], 1)
    
    # Crear enjambre
    swarm = SwarmHibrido(evaluarParticula, nGeneradores, datos, 
                        nParticle=nParticle, nInter=nInter)
    
    # Inicializar el enjambre
    initSwarm(swarm)
    
    # Ejecutar el algoritmo
    updateSwarm(swarm)
    
    return swarm.gBest_bin, swarm.gBest_pot, swarm.fitgBest
end 