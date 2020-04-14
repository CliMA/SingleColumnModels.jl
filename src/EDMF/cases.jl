#### Cases

export Case
export BOMEX

"""
    Case

An abstract case, which encompasses
all EDMF benchmark cases.
"""
abstract type Case end

"""
    BOMEX

The Barbados Oceanographic Meteorological
Experiment (BOMEX).

Reference:
  Kuettner, Joachim P., and Joshua Holland.
  "The BOMEX project." Bulletin of the American
  Meteorological Society 50.6 (1969): 394-403.
"""
struct BOMEX <: Case end

"""
    Soares

Reference:
  Soares, P.M.M., Miranda, P.M.A., Siebesma, A.P. and Teixeira, J. (2004),
  An eddy‐diffusivity/mass‐flux parametrization for dry and shallow cumulus convection.
  Q.J.R. Meteorol. Soc., 130: 3365-3383. doi:10.1256/qj.03.223
"""
struct Soares <: Case end

