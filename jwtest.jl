function f1(x)
  return x, x^2
end

function f2(x)
  return f1(x)
end

f2(3)