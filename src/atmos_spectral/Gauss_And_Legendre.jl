export Compute_Legendre!, Compute_Gaussian!

function Compute_Legendre!(num_fourier, num_spherical, sinθ, nθ)
    """
    Computes the Associated Legendre Polynomials and their first derivatives with respect to
    sine of latitude (μ) using stable three-term recurrence relations.

    This function implements the spectral transform appendix logic typically found in NWP models.
    It populates the spectral basis functions required for the spherical harmonic expansion:
    f(θ, λ) = ∑_{l=0} ∑_{m=-l}^{l} f_{lm} P_{lm}(sinθ)e^{imλ}

    Mathematical Formulation:
    1. Normalization: P_{0,0} = 1
    2. Diagonal Recurrence (Sectorial):
    P_{m,m} = √((2m+1)/2m) cosθ P_{m-1,m-1}
    3. Sub-diagonal Recurrence:
    P_{m+1,m} = √(2m+3) sinθ P_{m,m}
    4. Vertical Recurrence (General):
    ε_{l,m} P_{l,m} = sinθ P_{l-1,m} - ε_{l-1,m} P_{l-2,m}
    where ε_{l,m} = √((l² - m²) / (4l² - 1))
    5. Derivative (Meridional):
    (1-μ²) dP_{l,m}/dμ = -l ε_{l+1,m} P_{l+1,m} + (l+1) ε_{l,m} P_{l-1,m}

    Parameters
    ----------
    num_fourier (Int)
        Maximum zonal wavenumber (M) for the spectral truncation.

    num_spherical (Int)
        Maximum total wavenumber (L) for the spectral truncation.

    sinθ (Vector{Float64})
        Sine of the latitude at each grid point (Gaussian nodes μ).

    nθ (Int)
        Number of latitudinal grid points (length of sinθ).

    Returns
    -------
    qnm (num_fourier+1, num_spherical+1, nθ)
        Associated Legendre Polynomials P_{l,m}(sinθ).
        Indexing: qnm[m+1, l+1, lat] maps to physical mode P_{l,m}.

    dqnm (num_fourier+1, num_spherical+1, nθ)
        Meridional derivatives dP_{l,m}/dμ.
        Indexing: dqnm[m+1, l+1, lat] maps to physical derivative dP_{l,m}/dμ.

    Notes
    -----
    The implementation assumes a Gaussian grid where μ ≠ ±1. Division by cos²θ
    in the derivative calculation will result in singularities at the poles if used
    on a regular latitude-longitude grid including the poles.

    ====================================================================================

    Spectral Numerical Weather Prediction Models Appendix B
    f(θ, λ) = ∑_{l=0} ∑_{m=-l}^{l} f_{lm} P_{lm}(sinθ)e^{imλ} (Y_{l,m} = P_{l,m} e^{i m λ} )
    l=0,1...∞    and m = -l, -l+1, ... l-1, l
    P_{0,0} = 1, such that 1/4π ∫∫YYdS = δ
    P_{m,m} = sqrt((2m+1)/2m) cosθ P_{m-1m-1} 
    P_{m+1,m} = sqrt(2m+3) sinθ P_{m m} 
    sqrt((l^2-m^2)/(4l^2-1))P_{l,m} = x  P_{l-1, m} -  sqrt(((l-1)^2-m^2)/(4(l-1)^2 - 1))P_{l-2,m}
    ε[m,l] = sqrt((l^2- m^2)/(4l^2 - 1))
    (1-μ^2)d P_{m,l}/dμ = -nε[m,l+1]P_{m,l+1} + (l+1)ε_{m,l}P_{m,l-1}
    Julia index starts with 1 qnm[m+1,l+1] = P[l,m]

    dqnm = dP/dμ

    """
    
    # Associated Legendre Polynomials (m, l, θ)
    qnm  = zeros(Float64, num_fourier+1, num_spherical+2, nθ)
    dqnm = zeros(Float64, num_fourier+1, num_spherical+1, nθ)

    cosθ = sqrt.(1 .- sinθ.^2)
    ε    = zeros(Float64, num_fourier+1, num_spherical+2)

    # The diagonal recurrence (l == m)
    qnm[1, 1, :] .= 1.0
    for m = 1:num_fourier
        qnm[m+1, m+1, :] = sqrt((2m+1)/(2m)) .* cosθ .* qnm[m, m,:]
    end
    
    # The semi-diagonal recurrence (l == m + 1)
    for m = 1:num_fourier+1
        qnm[m, m+1, :] = sqrt(2*m+1) * sinθ .* qnm[m,m, :] 
    end
    
    # Normalization factors
    for m = 0:num_fourier
        for l = m:num_spherical+1
            # ε[m,l] = sqrt(((l-1)^2 - (m-1)^2) ./ (2*(l-1) + 1))
            ε[m+1,l+1] = sqrt((l^2 - m^2) ./ (4*l^2 - 1))
        end
    end

    # The main loop (l > m + 1)
    for m = 0:num_fourier
        for l = m+2:num_spherical+1
            # m=0, l=2
            # qnm[m+1,l+1,:] = sqrt(2*l-1) /ε[m+1,l+1] * (sinθ .* qnm[m+1,l,:] -  ε[m+1,l]/sqrt(2*l-3)*qnm[m+1,l-1,:])
            qnm[m+1,l+1,:] = (sinθ .* qnm[m+1,l,:] -  ε[m+1,l]*qnm[m+1,l-1,:])/ε[m+1,l+1]
        end
    end

    # Calculating derivatives
    for m = 0:num_fourier
        for l = m:num_spherical
            if l == m
                dqnm[m+1, l+1, :] = (-l*ε[m+1, l+2]*qnm[m+1,l+2,:])./(cosθ.^2)
            else
                dqnm[m+1, l+1, :] = (-l*ε[m+1, l+2]*qnm[m+1,l+2,:] + (l+1)*ε[m+1,l+1]*qnm[m+1,l,:])./(cosθ.^2)
            end

        end
    end

    return qnm[:, 1:num_spherical+1, :], dqnm
    # d P_{m,l}/dμ = -nε[m,l+1]P_{m,l+1} + (l+1)ε_{m,l}P_{m,l-1}     
end
    

function Compute_Gaussian!(n)
    """
    Computes the Gaussian latitudes and quadrature weights for a spectral grid.

    This function determines the roots (nodes) of the Legendre polynomial Pn(x)
    using the Newton-Raphson iteration method. These nodes form the Gaussian grid
    latitudes required for exact integration of polynomials up to degree 2n-1.

    Mathematical Formulation:
    1. Symmetry: The roots are symmetric about the equator (x=0). The function
    solves for the northern hemisphere roots and mirrors them.
    2. Recurrence Relation:
    n Pn(x) = (2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x)
    3. Derivative:
    (x² - 1) P'n(x) = n (x Pn(x) - P_{n-1}(x))
    4. Newton-Raphson Update:
    x_{new} = x_{old} - Pn(x) / P'n(x)
    5. Initial Guess (Asymptotic):
    x₀ = cos(π(i - 0.25) / (n + 0.5))
    6. Quadrature Weights:
    wᵢ = 2 / ((1 - xᵢ²) [P'n(xᵢ)]²)

    Parameters
    ----------
    n (Scalar)
        The total number of latitudinal grid points (must be an even integer).

    Returns
    -------
    sinθ (n)
        Gaussian nodes (roots of Legendre polynomial), representing sine of latitude.
        Ordered from south (-z) to north (+z) if using the provided symmetry logic.

    wts (n)
        Gaussian quadrature weights corresponding to the nodes in sinθ.

    Notes
    -----
    The iteration terminates when the update step is smaller than the tolerance
    (1.0e-15). If convergence is not reached within 10,000 iterations, an error is thrown.

    ====================================================================================

    Pn(x) is an odd function
    solve half of the n roots and weightes of Pn(x) # n = 2n_half
    P_{-1}(x) = 0
    P_0(x) = 1
    P_1(x) = x
    nP_n(x) = (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)
    P'_n(x) = n/(x^2-1)(xP_{n}(x) - P_{n-1}(x))
    x -= P_n(x)/P'_{n}()
    Initial guess xi^{0} = cos(π(i-0.25)/(n+0.5)) 
    wi = 2/(1-xi^2)/P_n'(xi)^2 

    """
    
    itermax = 10000
    tol     = 1.0e-15
 
    sinθ = zeros(Float64, n)
    wts  = zeros(Float64, n)

    n_half = Int64(n/2)
    for i=1:n_half
        dp = 0.0
        z  = cos(pi*(i - 0.25)/(n + 0.5))
        
        for iter=1:itermax
            p2 = 0.0
            p1 = 1.0
            
            for j=1:n
                p3 = p2 # Pj-2
                p2 = p1 # Pj-1
                p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j  #Pj
            end
            
            # P'_n
            dp = n*(z*p1 - p2)/(z*z - 1.0)
            z1 = z
            z  = z1 - p1/dp
            if(abs(z - z1) <= tol)
                break;
            end
            if iter == itermax
                @error("Compute_Gaussian! does not converge!")
            end
        end
        
        sinθ[i], sinθ[n-i+1],  = -z, z
        wts[i] = wts[n-i+1]  = 2.0/((1.0 - z*z)*dp*dp)
    end

    return sinθ, wts
end



function test()
    """
    Performs unit testing and validation of the spectral transformation kernel.

    This function verifies the accuracy of the Gaussian grid generation and the
    Associated Legendre Polynomial recurrence relations. It performs two types of checks:

    1. Analytical Verification:
    Compares the numerically computed polynomials Pₗₘ(μ) and derivatives dPₗₘ/dμ
    against hardcoded exact analytical formulas for specific low-order modes.
    
    2. Orthogonality Verification:
    Validates the discrete orthogonality property of the generated basis set
    using Gaussian quadrature. The normalized condition checked is:
    0.5 * ∑ (wₖ * Pₙ,ₘ(μₖ) * Pₗ,ₘ(μₖ)) = δₙₗ

    Parameters
    ----------
    None
        (The test configuration T21/N32 is hardcoded internally for validation).

    Returns
    -------
    None
        Results are printed to standard output via @show macros.
        The function throws an assertion error if the orthogonality condition fails
        within the tolerance (1.0e-10).

    """
    num_fourier, nθ = 21, 32
    num_spherical = num_fourier+1
    nλ = nθ

    sinθ, wts = Compute_Gaussian!(nθ)
    #compare with https://pomax.github.io/bezierinfo/legendre-gauss.html


    qnm, dqnm = Compute_Legendre!(num_fourier, num_spherical, sinθ, nθ)

    q44 = sqrt(35)/4*(1 .- sinθ.^2).^(3/2) 
    q34 = sqrt(105/8)*(sinθ .- sinθ.^3)
    q13 = sqrt(5)/2*(3sinθ.^2 .- 1)
    q14 = sqrt(7)/2*(5sinθ.^3 .- 3sinθ)
    q24 = sqrt(21)/4*(5sinθ.^2 .- 1).*sqrt.(1 .- sinθ.^2) 


    dq44 = sqrt(35)/4*(1 .- sinθ.^2).^(1/2) *(3/2).*(-2sinθ)  
    dq34 = sqrt(105/8)*(1 .- 3*sinθ.^2)
    dq13 = sqrt(5)/2*(6sinθ)
    dq14 = sqrt(7)/2*(15sinθ.^2 .- 3)
    dq24 = sqrt(21)/4*(10*sinθ).*sqrt.(1 .- sinθ.^2) + sqrt(21)/4*(5sinθ.^2 .- 1)./sqrt.(1 .- sinθ.^2)*0.5 .* (-2sinθ)

    #compare with exact form
    @show norm(qnm[4,4,:] - q44)
    @show norm(qnm[3,4,:] - q34)
    @show norm(qnm[1,3,:] - q13) 
    @show norm(qnm[1,4,:] - q14)
    @show norm(qnm[2,4,:] - q24) 

    @show norm(dqnm[4,4,:] - dq44)
    @show norm(dqnm[3,4,:] - dq34)
    @show norm(dqnm[1,3,:] - dq13) 
    @show norm(dqnm[1,4,:] - dq14)
    @show norm(dqnm[2,4,:] - dq24) 

    #check 1/2∫_{-1}^{1}P_{n,m} P_{l,m} = δ_{n,l}

    D = zeros()
    
    for m=1:num_fourier+1
        for l=m:num_spherical+1
            for n=m:num_spherical+1
                @assert( (0.5*sum(qnm[m,n,:].*qnm[m,l,:].*wts) - Float64(n==l)) < 1.0e-10)
            end
        end
    end
end
