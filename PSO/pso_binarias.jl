const Float = Float64
using CSV
using DataFrames


##    Creación de la partícula

## `np::Int`: número de partículas que se van a usar  

## `position::Array{Float, 1}`: posición actual de la partícula

## `velocity::Array{Float, 1}`: velocidad actual de la partícula

## `pBest::Array{Float, 1}`: mejor posición local -> mejor posición de la partícula en todo su historial de posiciones

## `lBest::Array{Float, 1}`: mejor posición global -> mejor posición del historial de posiciones del grupo o enjambre de partículas 

## `fitValue::Float`: valor que indica la calidad de la solución encontrada 

## `fitpBest::Float`: calidad de la solución local

## `fitlBest::Float`: calidad de la solución global

## `nFitEval::Int`: número de veces que se ha evaluado la calidad de la solución durante la ejecución 
## (solo se contabilizan las soluciones que estén dentro de los límites definidos del espacio de partículas)


# Definimos la estructura que va a tener la partícula. Mutable implica que después de definir la partícula, estos datos puede modificarse.
# Incluye una función: asignación de valores para caracterizar la partícula 
mutable struct Particle
    nParticle::Int

    position::Array{Float, 1}
    velocity::Array{Float, 1}
    pBest::Array{Float, 1}
    lBest::Array{Float, 1}
        
    fitValue::Float
    fitpBest::Float
    fitlBest::Float
        
    nFitEval::Int
        
    # inicializamos la partícula 
    # con lBest = pBest = position = números aleatorios
    # y fitlBest = fitpBest = fitValue = Infinito
    # Esta configuración se actualizará al inicializar la función swarm
    function Particle(nParticle::Int) #Indico en la función el número de partículas (tamaño del enjambre)
        position = rand(0:1.0, nParticle) # La posición debe estar en 0 o en 1 
        velocity = rand(nParticle) - position
        pBest = position
        lBest = position
            
        fitValue = Inf
        fitpBest = fitValue
        fitlBest = fitValue
            
        nFitEval = 0
            
        new(nParticle, position, velocity, pBest, lBest, 
                fitValue, fitpBest, fitlBest, nFitEval) # Función que crea una nueva partícula con la estructura definida previamente 
    end       
end

    # test: creo una instancia con 4 partículas donde las posiciones son 0 o 1, las valoraciones son infinito
    # el número de evaluaciones es 0 y la velocidad un número aleatorio. 
    # p1 = Particle(4) 

 ##    function fitFunc(p::Particle)

    ## Defino la función que necesito para calcular la mejor combinación de apertura de generadores. Hay que tener en 
    ## cuenta una serie de limitaciones:
    ##        - La generación debe satisfacer la demanda 
    ##        - No se debe transportar más energía de la que es capaz de soportar la línea
    ##        - La generación debe llegar a todos los nodos, por lo que hay que tener en cuenta la conexión gen-bus
    
function fitFunc(p::Particle, c::String)
    function generacion_total(p::Particle)
        datosGeneradores = CSV.read("Casos/$c/datosGeneradores.csv", DataFrame)
        p_tot = 0.0
        for i in 1: nrow(datosGeneradores)
            if p.position[i] == 1
                    p_tot += datosGeneradores.P_MAX[i]
            end
        end
        return p_tot
    end
    
    function coste_total(p::Particle)
        cost = 0.0
        datosGeneradores = CSV.read("Casos/$c/datosGeneradores.csv", DataFrame)
        for i in 1:nrow(datosGeneradores)
            # Cuando un generador está abierto, se calcula es coste suponiendo que genera su máxima potencia
            if p.position[i] == 1 
                cost_unit = datosGeneradores.P_COSTE1[i]
                p_max = datosGeneradores.P_MAX[i]
                cost += cost_unit*p_max
            end
        end
        return cost
    end
    
    function demanda_total() # Función que penaliza si la generación no cubre la demanda
        datosNodos = CSV.read("Casos/$c/datosNodos.csv", DataFrame)
        demandap = 0.0
        for i in 1:nrow(datosNodos)
            demandap = demandap + datosNodos.PD[i]
        end
    
        return demandap
    end
    
    p_tot = generacion_total(p)
    demp_tot = demanda_total()
    if p_tot < demp_tot 
        coste_tot = Inf
    else 
        coste_tot = coste_total(p)
    end
    
    return coste_tot
end

##   initFitValue!(fitFunc::Function, p::Particle)

##Se crea una función para inicializar el valor de la variable "fitValue" desde la partícula "p" utilizando la función "fitFunc".

function initFitValue!(fitFunc::Function, p::Particle, c::String)
    p.fitValue = fitFunc(p, c) 
        
    # Actualizamos el valor nFitEval
    p.nFitEval += 1
    nothing
end


    # test
    # p2 = Particle(3) El número de partículas es el número de generadores que hay
    # fitFunc(p2)
    # println(p2)
    # initFitValue!(fitFunc, p2)
    # println(p2) El vaalor actual de fitValue es el resultado de fitFunc y se suma uno al número de evaluaciones 



##    updatePositionAndFitValue!(fitFunc::Function, nDim::Int, p::Particle)

## Función que actualiza la posición y el fitValue de la partícula p usando la función "fitFunc" que se defina por el usuario
## con "nParticle" parámetros.

function updatePositionAndFitValue!(fitFunc::Function, p::Particle, c::String)
    p.position += p.velocity #según las ecuaciones del PSO "position[k+1] = position[k] + velocity[k+1]"

    # Cambiar las posiciones a binarias
    for i in 1:p.nParticle
        if p.position[i] >= 0.5 # Si la posición es mayor que 0.5, el generador se  abre 
            p.position[i] = 1.0
        else 
            p.position[i] = 0.0
        end
    end    
    # si la posición está fuera de los límites del espacio definidos, fijamos fitValue = infinito
    for x in p.position
        if (x < 0 || x > 1)
            p.fitValue = Inf
            return
        end
    end
    # actualizar nFitEval
    p.nFitEval += 1
    p.fitValue = fitFunc(p, c)
        
    nothing
end

    # test
    # p2 = Particle(3)
    # fitFunc(p2)
    # initFitValue!(fitFunc, p2)
    # println(p2)
    # updatePositionAndFitValue!(fitFunc, p2) Se calcula la nueva posición y la calidad del nuevo resultado (FitValue)
    # println(p2)



##    updatepBestAndFitpBest!(p::Particle)

## Función que actualiza la mejor solución local y la valoración de cada partícula p.

function updatepBestAndFitpBest!(p::Particle)
    if p.fitValue < p.fitpBest # Si el coste es menor que la mejor solución que había hasta ahora
                               # se actualiza la pBest
        p.fitpBest  = p.fitValue
        p.pBest = p.position
    end
    nothing
end

    # test
    # p = Particle (3) # 3 generadores con diferentes posiciones cada una 
    # updatePositionAndFitValue!(fitFunc, p)
    # p #Se actualiza el valor de fitValue, que ahora es el coste y el número de evaluaciones, que ahora es 1 
    # updatepBestAndFitpBest!(p)
    # p #Se actualiza la mejor solución de la partícula y la mejor solución local, así como sus fit





##    updateVelocity!(p::Particle, w::Float, c1::Float=c1, c2::Float=c2)

## Función que actualiza la velocidad de la partícula p (se siguen las ecuaciones del PSO)

function updateVelocity!(p::Particle, w::Float, c1::Float, c2::Float)
    p.velocity = w * p.velocity + 
    c1 * rand() * (p.pBest - p.position) + 
    c2 * rand() * (p.lBest - p.position)
        
    nothing
end

    # test
    # updateVelocity!(p2, 0.8, 2.0, 2.0)



##    neiborIndices(i::Int, nNeighbor::Int, nParticle::Int)

## Devuelve el índice de los vecinos de la partícula i. Esto sirve para establecer las conexiones y buscar la solcuión global.

function neiborIndices(i::Int, nNeibor::Int, nParticle::Int)
        
    # número de vecinos debe ser superior a 3, la función devuelve el máximo entre 3 y nNeibor
    nNeibor = max(3, nNeibor)
        
    # número de vecinos por la izquierda de la partícula i
    nLeft = (nNeibor - 1) ÷ 2
        
    # índice de la partícula inicial en el grupo de solciones locales 
    startIndex = (i - nLeft)
        
    # índice de la útlima partícula del grupo local
    endIndex = startIndex + nNeibor -1
        
    # guardamos los índices desde el inicial al final
    indices = collect(startIndex:endIndex)
        
    # ajuste de índices para que estén dentro del rango(1:nParticle)
    for i in 1:nNeibor
        if indices[i] < 1
            indices[i] += nParticle
        elseif indices[i] > nParticle
            indices[i] -= nParticle
        end
    end
        
    indices
end

    # test
    # neiborIndices(1, 3, 40)  
    # start index = 0, end index = 2 => 0+40 (porque 0<1), 1, 2



##    Creación del enjambre

## `fitFunc::Function`: función que da la valoración de la solución respecto a la óptima. 

## `nDim::Int`: número de generadores que tenemos en la red 

## `nParticle::Int`: número de posibles soluciones que vamos a explorar 
        
## `nNeibor::Int`: número de vecinos (partículas) en un grupo local
        
## `nInter::Int`: número de iteraciones que hay que realizar para que cada partícula avance hacia una mejor solución

## `c1::Float`: constante cognitiva 

## `c2::Float`: constante social

## `wMax::Float`: el valor máximo del peso de inercia (una componente de la velocidad es la inercia de la partícula)

## `wMin::Float`: valor mínimo de la inercia

## `w::Float`: valor actual de la inercia
        
## `gBest::Array{Float, 1}`: posición en la que el enjambre tiene la mejor solución de todo su historial
        
## `fitgBest::Float`: valoración de la solución global
        
## `particles::Array{Particle, 1}`: partículas del enjambre (no el número, sino su caracterización) t

## `nFitEvals::Int`: número de evaluaciones que hay que realizar para la valoración de la solución, no incluye los que salen fuera de
##   los límites del espacio definidos

# estructura cuyos valores pueden modificarse posteriormente. 
# Inlcuye una función: asignar valores al enjambre para caracterizarlo 
mutable struct Swarm
    fitFunc::Function
    nDim::Int
        
    nParticle::Int
    nNeibor::Int
    nInter::Int
        
    c1::Float
    c2::Float
        
    wMax::Float
    wMin::Float
    w::Float
        
    gBest::Array{Float, 1}    
    fitgBest::Float
        
    particles::Array{Particle, 1}
        
    nFitEvals::Int

    c::String
        
    # inicializar el enjambre
    function Swarm(fitFunc::Function, nDim::Int; 
            nParticle::Int=3, 
            nNeibor::Int=3, nInter::Int=2000,
            c1::Float=2.0, c2::Float=2.0,
        wMax::Float=0.9, wMin::Float=0.4, c::String)
            
        if nNeibor > nParticle
            error("El número de partículas en un grupo local no debe exceder 
                el número total de partículas en el enjambre")
        end    
            
        w = wMax
            
        gBest = rand(nDim) # Damos la mejor solución global a una posición aleatoria
        fitgBest = Inf # con una valoración infinita (esto es malo porque queremos que la diferencia entre el óptimo y la solución sean 
                       # mínimas)
            
        # Inicializar el enjambre con nParticle
        particles = [Particle(nDim) for i in 1:nParticle] 
            
            
        nFitEvals = 0 # Definimos el valor incial de las evaluaciones realizas
                    
        new(fitFunc, nDim, nParticle, nNeibor, nInter, 
            c1, c2, wMax, wMin, w, gBest, 
            fitgBest, particles, nFitEvals, c)        
    end       
end

    # test
    # nDim = 2
    # s = Swarm(fitFunc, nDim)



##    updatelBestAndFitlBest!(s::Swarm)

## Actualiza el valor de la mejor solución local (lBest) y su valoración en el enjambre s.
        
function updatelBestAndFitlBest!(s::Swarm)
    for i in 1:s.nParticle
        neiborIds = neiborIndices(i, s.nNeibor, s.nParticle) # Vemos los índices de las partículas que están en el mismo grupo
        neiborFits = [s.particles[Id].fitValue for Id in neiborIds] # Vemos las valoraciones de las soluciones obtenidas por cada partícula
        fitlBest, index = findmin(neiborFits) # Función que devuelve la solución con menor fit y el índice que la contiene
            
        # Si el fit encontrado es menor que todos los que hay en el enjambre, se actualiza:
                # la mejor solución local
                # la mejor solución local en el enjambre 
                # la mejor valoración local
        if fitlBest < s.particles[i].fitlBest 
            lBest = s.particles[neiborIds[index]].position
            s.particles[i].lBest = lBest
            s.particles[i].fitlBest = fitlBest
        end
    end
    nothing
end

    # test
    # nDim = 2
    # s = Swarm(fitFunc, nDim)  
    # updatelBestAndFitlBest!(s)
    # s.particles[1]




##    updategBestAndFitgBest!(s::Swarm)

## Actualiza la mejor solución global y su valoración para el enjambre s.
      
function updategBestAndFitgBest!(s::Swarm)
        
    gFits = [particle.fitValue for particle in s.particles] # Almaceno las valoraciones globales de todas las partículas del enjambre
    fitgBest, index = findmin(gFits) # Guardo la mínima valoración global (mejor solución) y el índice de la partícula que la tiene
        
    # Si la valoración global que he guardado es menor a la que tenía antes en el enjambre, actualizo la posición de la mejor solución 
    # global y actualizo la valoración de la solución global. 
    if fitgBest < s.fitgBest
        s.gBest = s.particles[index].position   
        s.fitgBest = fitgBest
    end
    nothing
end

    # test
    # fitFunc(x) = x[1]^2 + x[2]^2
    # nDim = 2
    # s = Swarm(fitFunc, nDim)  
    # updategBestAndFitgBest!(s)
    # s.gBest



##    initSwarm(s::Swarm)

## Inicialización del enjambre s.

function initSwarm(s::Swarm, c::String)
        
    # inicializar el fitValue para cada partícula 
    for particle in s.particles
        initFitValue!(s.fitFunc, particle, c)
        updatepBestAndFitpBest!(particle) #Busco las mejores soluciones de las partículas 
    end
        
    # actualizar la mejor solución local y su valoración para el enjambre 
    updatelBestAndFitlBest!(s)
        
    # actualizar la mejor solución global y su valoración para el enjambre 
    updategBestAndFitgBest!(s)
        
    nothing
end

    # test
    # nDim = 2
    # s = Swarm(fitFunc, nDim, nParticle=4)
    # initSwarm(s)
    # s




##    updateInertia!(s::Swarm)

## Actualizar la inercia después de cada iteración. 

function updateInertia!(s::Swarm)
    dw = (s.wMax - s.wMin)/s.nInter
    s.w -= dw
        
    nothing
end

    # test
    # fitFunc(x) = x[1]^2 + x[2]^2
    # nDim = 2
    # s = Swarm(fitFunc, nDim, nParticle=4)
    # updateInertia!(s)
    # s.w



##    updateVelocity!(s::Swarm)

## Actualizar la velocidad de cada partícula del enjambre. 

function updateVelocity!(s::Swarm)
    for particle in s.particles
        updateVelocity!(particle, s.w, s.c1, s.c2)
    end        
    nothing
end

    # test
    # fitFunc(x) = x[1]^2 + x[2]^2
    # nDim = 2
    # s = Swarm(fitFunc, nDim, nParticle=4)
    # initSwarm(s)
    # println(s.particles[1])
    # updateVelocity!(s)
    # println(s.particles[1])



##     updatePositionAndFitValue!(s::Swarm)

## Actualizar la posición y la valoración de la solución de cada partícula del enjambre. 

function updatePositionAndFitValue!(s::Swarm, c::String)
    for particle in s.particles
        updatePositionAndFitValue!(s.fitFunc, particle, c)
        updatepBestAndFitpBest!(particle)
    end        
    nothing
end

    # test
    # fitFunc(x) = x[1]^2 + x[2]^2
    # nDim = 2
    # s = Swarm(fitFunc, nDim, nParticle=4)
    # initSwarm(s)
    # println(s.particles[1])
    # updatePositionAndFitValue!(s)
    # println(s.particles[1])



##    updateSwarm(s::Swarm)
## Función que junta todas las funciones previas. 

function updateSwarm(s::Swarm, nInter:: Int64, c::String)
    for i in 1:nInter
        # actualiza la velocidad de cada partícula en el enjambre 
        updateVelocity!(s::Swarm)
            
        # actualiza la posición y la valoración de la solución de cada partícula del enjambre 
        updatePositionAndFitValue!(s::Swarm, c::String)
            
        # actualiza la mejor solución local y su valoración para cada partícula del enjambre
        updatelBestAndFitlBest!(s::Swarm)
            
        # uactualiza la mejor solución global y su valoración para cada partícula del enjambre
        updategBestAndFitgBest!(s::Swarm) 
            
        # actualiza la inercia de cada partícula del enjambre
        updateInertia!(s::Swarm)

        println(s)
    end

    nothing 
end

    # test
    # nDim = 3
    # s = Swarm(fitFunc, nDim, nParticle=4) #Vamos a explorar 4 posibles soluciones de tres generadores 
    # initSwarm(s)
    # s Me  inicializa todas las partículas y cambia el valor de la solución local con fitFunc, también guarda la mejor
        # solución local ( que será la misma para todas las partículas)
    # println(s.particles[1]) Me imprime la posición, tengo dos valores porque la dimensión es 2
    # updateSwarm(s, 2) Repito la actualización 2 veces para obtener la solución final
    # println(s.particles[1]) Puedo ver la solución de la partícula con índice i (en este caso índice 1)
    # updateSwarm(s)
    # s Se actualizan las posiciones y las velocidades, obteniendo una nueva solución. También se guarda la mejor 
    # de cada partícula y la mejor solución local. 
    # println(s.particles[1])
