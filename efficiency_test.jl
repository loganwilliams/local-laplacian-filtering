A = rand(1000,1000)

function remap(A)
  for i = 1:1000
    for j = 1:1000
      if (A[i,j] > 0.5)
        A[i,j] = A[i,j]^2
      else
        A[i,j] = A[i,j]/2
      end
    end
  end
end

function remap2(subimage, sigma_r, g_0, f_d, f_e)
  for i = 1:size(subimage)[1]
    for j = 1:size(subimage)[2]
      pv = subimage[i,j]

      if abs(pv - g_0) < sigma_r
        pv = g_0 + sign(pv - g_0)*sigma_r*f_d(abs(pv-g_0)/sigma_r)
      else
        pv = g_0 + sign(pv - g_0)*(f_e(abs(pv-g_0) - sigma_r) + sigma_r)
      end

      subimage[i,j] = pv
    end
  end

  return subimage
end

function remap3(subimage, sigma_r, g_0, alpha, beta)
  for i = 1:size(subimage)[1]
    for j = 1:size(subimage)[2]
      pv = subimage[i,j]

      if abs(pv - g_0) < sigma_r
        pv = g_0 + sign(pv - g_0)*sigma_r*(abs(pv-g_0)/sigma_r)^alpha
      else
        pv = g_0 + sign(pv - g_0)*(beta*(abs(pv-g_0) - sigma_r) + sigma_r)
      end

      subimage[i,j] = pv
    end
  end

  return subimage
end

function generate_powercurve(alpha)
  function powercurve(i)
    return i^alpha
  end

  return powercurve
end

function generate_tonemapping(beta)
  function tonemapping(i)
    return i*beta
  end

  return tonemapping
end


function remap5(A)
  for i = 1:1000
    for j = 1:1000
        A[i,j] = A[i,j]^2
    end
  end
  return A
end

# @time remap(A)
# @time remap(A)
# @time remap2(A, 0.1, 0.5, generate_powercurve(0.1), generate_tonemapping(1))
# @time remap2(A, 0.1, 0.5, generate_powercurve(0.1), generate_tonemapping(1))
# @time remap3(A, 0.1, 0.5, 0.1, 1)
# @time remap3(A, 0.1, 0.5, 0.1, 1)
# @time remap4(A, 0.1, 0.5, 0.1, 1)
# @time remap4(A, 0.1, 0.5, 0.1, 1)

@time remap5(A)
@time remap5(A)

@time A.^2
@time A.^2