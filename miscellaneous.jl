function jwstop()
    println("---------------------------")
    println("---------------------------")
    println("---------------------------")
    println("--------No Problem!--------")
    println("---------------------------")
    println("---------------------------")
    println("---------------------------")
    error("stop here")
 end
 
 function jwdisplay(x,str)
    println("---------------------------")
    println("$str is ")
    display(x)
 end
 
 function jwchk(x::Bool)
    if !x
        error("sanity check failed")
    end
 end