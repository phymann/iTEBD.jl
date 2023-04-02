let
  a = 1
  while true
    a += 1
    @show a
    if a > 10
      break
    end
  end
end