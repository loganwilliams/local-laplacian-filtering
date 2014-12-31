# generates a lower level of the Gaussian pyramid decomposition of an image. Returns both a downsampled image and an original size image
function impyramid(A; sizes="downsampled")
  a = 0.375
  kernel = [1/4 - a/2, 1/4, a, 1/4, 1/4 - a/2]
  B = imfilter(imfilter(A,kernel), kernel') # slowest part of impyramid, by far

  if (sizes == "downsampled")
    return (B[1:2:end, 1:2:end])
  elseif (sizes == "both")
    return (B, B[1:2:end,1:2:end])
  end
end

function gaussianpyramid(A, levels)
  pyramid = ()
  Ns = A

  for n = 1:levels
    pyramid = tuple(pyramid..., Ns)
    Ns = impyramid(Ns)
  end

  return pyramid
end

# generates a laplacian pyramid decomposition of an image A, with nlevels + 1 residual level
function lappyramid(A, nlevels)
  residual = A
  pyramid = ()

  for n = 1:nlevels
    (N, Ns) = impyramid(residual, sizes="both")
    L = residual - N
    residual = Ns
    pyramid = tuple(pyramid..., L)
  end

  return tuple(pyramid..., residual)
end

# reconstruct an image (approximate, for the time being), from its pyramidal decomposition
function reconstruct(L)
  im = last(L)

  for n = (length(L)-1):-1:1
     # upscale
    im2 = zeros(2*size(im)[1], 2*size(im)[2])
    im2[1:2:end,1:2:end, :] = im

    # gaussian interpolation
    kernel = [1/18, 1/2, 16/18, 1/2, 1/18];
    im2 = imfilter(imfilter(im2, kernel, "reflect"), kernel', "reflect")

    im = im2
    im += L[n]
  end

  return im
end

function remap(subimage, sigma_r, g_0, alpha, beta)
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

function locallaplacianfilter(I; sigma_r=0.1, levels=5, alpha=0.1, beta=1)
  G = gaussianpyramid(I, levels)
  Loutput = lappyramid(I, levels)

  for n = 1:levels
    for i = 1:size(G[n])[1]
      for j = 1:size(G[n])[2]
        print("Current index: ") # printing is kind of slow
        println((n,i,j))
        g_0 = G[n][i,j]

        # determine subbregion (TBD)
        R = copy(I)

        # remap subregion
        Rtilde = remap(R, sigma_r, g_0, alpha, beta)

        # generate laplacian pyramid of remapped region
        Ltilde = lappyramid(Rtilde, n)

        #update output laplacian pyramid
        Loutput[n][i,j] = Ltilde[n][i,j]
      end
    end
  end

  return reconstruct(Loutput)
end

function fastlocallaplacianfilter(I; sigma_r=0.1, levels=5, alpha=0.1, beta=1)
  G = gaussianpyramid(I, levels)
  Loutput = lappyramid(I, levels)
  pyramid_stack = ()

  # this should be set by looking at the spectrum of the remapping function. this is good enough, though.
  spacing = 0.01

  println("Pregenerating pyramid stack")
  for g = 0:spacing:1
    pyramid_stack = tuple(pyramid_stack..., lappyramid(remap(copy(I), sigma_r, g, alpha, beta), levels))
  end

  println("Starting pyramid rebuild")
  for n = 1:levels
    for i = 1:size(G[n])[1]
      for j = 1:size(G[n])[2]
        g_0 = G[n][i,j]

        #update output laplacian pyramid
        i1 = floor(g_0/spacing) + 1
        i2 = i1 + 1
        if (i2 >= (1/spacing + 1))
          i1 -= 1
          i2 -= 1
        end

        alpha = (g_0 - (i1-1)*spacing)/spacing

        Loutput[n][i,j] = pyramid_stack[i1][n][i,j] * alpha + pyramid_stack[i2][n][i,j] * (1 - alpha)
      end
    end
  end

  return reconstruct(Loutput)
end
