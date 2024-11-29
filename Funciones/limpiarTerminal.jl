#### Funcón llamada por varias funciones del código para limpiar la terminal

function limpiarTerminal()

    # En caso de que el terminal sea de Windows
    if Sys.iswindows()
        Base.run(`cmd /c cls`)

    # En caso de otros terminales basados en Unix
    else
        Base.run(`clear`)
    end
    
end