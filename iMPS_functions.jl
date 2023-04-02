"""
 normalize a centrally canonical iMPS
"""
function normalizeiMPS(Γ, λ, η=1.0)
  λnrm = norm(λ)
  λ /= λnrm
  Γ *= λnrm/√η
  return Γ, λ
end