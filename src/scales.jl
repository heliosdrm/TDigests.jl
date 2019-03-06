"""
    approx_asin(x)
    
Fast approximation to the `asin` function by Abramowitz & Stegun . (4.4.45, p. 81)
(M. Abramowitz and I. A. Stegun, *Handbook of mathematical functions*, 10th printing.
National Bureau of Standards, url: http://people.math.sfu.ca/~cbm/aands/)
"""
function approx_asin(x)
    if x ≥ 0
        return π/2 - √(1-x)*(1.5707288 - 0.2121144x + 0.0742610(x^2) - 0.0187293(x^3))
    else
        return -π/2 + √(1+x)*(1.5707288 + 0.2121144x + 0.0742610(x^2) + 0.0187293(x^3))
    end
end

"""
    approx_sin(x)
    
Fast approximation to the `sin` function by Abramowitz & Stegun. (4.3.98, p. 76)
(M. Abramowitz and I. A. Stegun, *Handbook of mathematical functions*, 10th printing.
National Bureau of Standards, url: http://people.math.sfu.ca/~cbm/aands/)
"""
function approx_sin(x)
    if x ≥ 0
        y = (0.5π - x)^2
        return 1.0 - 0.49670(y) + 0.03705(y*y)
    else
        y = (0.5π + x)^2
        return -1.0 + 0.49670(y) - 0.03705(y*y)
    end
end

"""
    k0(q, δ)
    k1(q, δ)

Scale functions for the *t*-digest.

`k0` is a linear scale function equal to (δ*q)/2; `k1` is the arcsine function,
defined between ±(δ/4) for `q` between [0,1] (see equation 3 in [1]).

---

    k0_inv(k, δ)
    k1_inv(k, δ)

Inverse scale functions, such that `k0_inv(k0(q, δ), δ) == q`, etc.

## References:

[1] Dunning, T. & Ertl, O. "Computing extremely accurate quantiles using t-digests",
arXiv: 1902.04023v1, 2019.

"""
k0(q, δ) = 0.5*(δ*q)
k1(q, δ) = δ/(2π)*approx_asin(2q - 1)

k0_inv(k, δ) = 2.0*k/δ
k1_inv(k, δ) = 0.5*approx_sin(2.0*π*k/δ) + 0.5

