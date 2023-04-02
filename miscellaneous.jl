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
        error("\n
        ================================
        sanity check failed
        ================================\n")
    end
 end

 function getmaxelm(Γ::ITensor)
   max1 = maximum(abs.(array(Γ)))
  #  println("maximal absolute matrix element is $(max1)")
   return max1
 end
 function getmaxelm(Γ::Matrix)
   max1 = maximum(abs.(Γ))
  #  println("maximal absolute matrix element is $(max1)")
   return max1
 end
 function getmaxelm(Γ::Array)
   max1 = maximum(abs.(Γ))
  #  println("maximal absolute matrix element is $(max1)")
   return max1
 end