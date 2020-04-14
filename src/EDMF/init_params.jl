#### InitParams

export Params

"""
    Params(::Case)

Initialize stand-alone input parameters to
solve the EDMF equations for a given case.
"""
function Params end

function Params(param_set, ::BOMEX)
  FT = Float64
  IT = Int
  params = Dict()

  #####
  ##### TODO: Parameters that need to be added to CLIMAParameters
  #####
  params[:k_Karman] = 0.4 # "Von Karman constant (unit-less)"

  #####
  ##### Filter parameters
  #####
  params[:param_set] = param_set
  params[:N_subdomains] = 3

  #####
  ##### IO parameters
  #####

  params[:export_data] = false
  params[:plot_single_fields] = true
  params[:export_frequency] = 2000


  #####
  ##### Filter parameters
  #####

  params[:a_bounds] = [1e-3, 1-1e-3]                                  # filter for a
  params[:w_bounds] = [0.0, 10000.0]                                  # filter for w
  params[:q_bounds] = [0.0, 1.0]                                      # filter for q

  #####
  ##### Space and time discretization
  #####

  params[:UniformGridParams] = UniformGridParams{FT,IT}(;z_min=0.0,
                                                         z_max=3000.0,
                                                         n_elems=75)
  grid_params = params[:UniformGridParams]
  params[:Δz] = (grid_params.z_max-grid_params.z_min)/grid_params.n_elems

  params[:Δt]     = FT(20.0)
  params[:Δt_min] = FT(20.0)
  params[:t_end]  = FT(21600.0)
  params[:CFL]    = FT(0.8)

  #####
  ##### External conditions
  #####

  # Looks good
  params[:ForcingType] = NoForcing()

  # Looks okay
  # params[:ForcingType] = StandardForcing(apply_subsidence=true,
  #                                        apply_coriolis=true,
  #                                        coriolis_param=FT(0.376e-4))
  #####
  ##### Physical models
  #####

  # params[:EntrDetrModel]          = BOverW2{FT}(1, 1)
  params[:EntrDetrModel]          = RH_Diff{FT}(0.1, 0.5, 2)
  # Looks okay
  params[:MixingLengthModel]      = ConstantMixingLength{FT}(100)
  # Looks okay
  # params[:MixingLengthModel]      = SCAMPyMixingLength{FT}(StabilityDependentParam{FT}(2.7,-100.0),
  #                                                          StabilityDependentParam{FT}(-1.0,-0.2))

  # Getting NaNs for TKE and other fields. Something needs to be fixed
  # params[:MixingLengthModel]      = IgnaciosMixingLength(StabilityDependentParam{FT}(2.7,-100.0),
  #                                                        StabilityDependentParam{FT}(-1.0,-0.2),
  #                                                        0.1, 0.12, 0.4, 40/13, 0.74)

  params[:EddyDiffusivityModel]   = SCAMPyEddyDiffusivity{FT}(0.1)

  params[:PressureModel]          = SCAMPyPressure{FT}(;buoy_coeff=FT(1.0/3.0),
                                                        drag_coeff=FT(0.375),
                                                        plume_spacing=FT(500.0))
  params[:prandtl_number]         = FT(1.0)
  params[:tke_diss_coeff]         = FT(2.0)

  params[:f_coriolis] = 0.376e-4                                      # coriolis force coefficient
  params[:f_c] = 0                                                    # buoyancy gradient factor

  params[:zrough] = 1.0e-4                                            # surface roughness

  params[:SurfaceModel]           = SurfaceFixedFlux{FT}(param_set;
                                                         T=FT(300.4),
                                                         P=FT(1.015e5),
                                                         q_tot=FT(0.02245),
                                                         ustar=FT(0.28),
                                                         windspeed_min=FT(0.0),
                                                         tke_tol=FT(0.01),
                                                         area=FT(0.1))

  params[:inversion_height] = [1.0 for i in 1:params[:N_subdomains]]  # inversion height
  params[:Ri_bulk_crit] = 0.0                                         # inversion height parameters
  _molmass_ratio::FT = FT(molmass_ratio(param_set))
  _grav::FT = FT(grav(param_set))
  params[:bflux] = (_grav * ((8.0e-3 + (_molmass_ratio-1)*(299.1 * 5.2e-5  + 22.45e-3 * 8.0e-3)) /(299.1 * (1.0 + (_molmass_ratio-1) * 22.45e-3))))
  params[:cq] = 0.001133                                              # Some surface parameter in SCAMPy
  params[:ch] = 0.001094                                              # Some surface parameter in SCAMPy
  params[:cm] = 0.001229                                              # Some surface parameter in SCAMPy
  params[:grid_adjust] = (log(20.0/params[:zrough])/log(params[:Δz]/2/params[:zrough]))^2 # Some surface parameter in SCAMPy
  params[:cq] *= params[:grid_adjust]                                 # Some surface parameter in SCAMPy
  params[:ch] *= params[:grid_adjust]                                 # Some surface parameter in SCAMPy
  params[:cm] *= params[:grid_adjust]                                 # Some surface parameter in SCAMPy

  return params
end
