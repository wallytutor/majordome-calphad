---@meta
---@diagnostic disable: undefined-global

------------------------------------------------------------------------------
-- N. C. Schieltz and M. R. Soliman
-- 10.1346/ccmn.1964.0130139
------------------------------------------------------------------------------

local QUARTZ_ALPHA = Substance {
  name = "QUARTZ_ALPHA",
  ranges = {
    Range(298.15, 6000.00, GibbsPolynomial {
        a =  8.58911225e+05, -- T^0
        b =  2.88907450e+02, -- T^1
        c = -4.69444800e+01, -- T ln T
        d = -1.71544000e-02, -- T^2
        e =  0.00000000e+00, -- T^3
        f =  0.00000000e+00, -- T^4
        g =  5.64840000e+05  -- T^-1
    })
  },
  elements = { Si = 1, O = 2 },
  aggregation_type = "Solid"
}

return { QUARTZ_ALPHA = QUARTZ_ALPHA }
