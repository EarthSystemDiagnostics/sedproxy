# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 09:32:35 2016
Last modified, Tue Oct 10 20:28:00 CEST 2017

@author: Didier M. Roche a.k.a. dmr
"""

__version__ = "1.0"

def delta_c(Tc, delta_w):
    # Inputs: Tc, temperature in C, delta_w d18O water in per mil
    # The equation to be computed is the solution
    # to the Kim & O'Neil 1997 equation
    # Written as: T = 16.1 - 4.64*(delta_c-delta_w) + 0.09*(delta_c-delta_w)**2
    # Noting in the following ukn = delta_c-delta_w

    a = 0.09
    b = -4.64
    c = 16.1-Tc
    delta = b**2 - 4*a*c

    # The only likely solution is ukn_2, the other one is non-physical
    # ukn_1 = (b*-1.0 + delta**(0.5))/(2*a)

    ukn_2 = (b*-1.0 - delta**(0.5))/(2*a)

    return ukn_2+delta_w-0.27
#enddef delta_c

# from http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_closest(A,target) :
    import numpy as np
    return (np.abs(A-target)).argmin()
#enddef find_closest

def make_cf_compliant(nc_var_nm):

    import sys

    secondsINyear = 31104000.0
    fac = secondsINyear * 10.0

    standard_names = {
                  "t2m":"air_temperature",
                  "pp":"precipitation_flux",
                  "swrs":"surface_downwelling_shortwave_flux_in_air",
                  "lon":"longitude",
                  "lat":"latitude",
                  "longitude":"longitude",
                  "latitude":"latitude",
                  "wetdays":"number_of_days_with_lwe_thickness_of_precipitation_amount_above_threshold",
                  "ts":"surface_temperature"
    }
    standard_units = {
                  "t2m":"K",
                  "pp":"kg m-2 s-1",
                  "swrs":"W m-2",
                  "lon":"degrees_east",
                  "lat":"degrees_north",
                  "longitude":"degrees_east",
                  "latitude":"degrees_north",
                  "wetdays":"days",
                  "ts":"K"
    }
    standard_facta = {
                  "t2m":"273.15",
                  "pp":"0.0",
                  "swrs":"0.0",
                  "lon":"0.0",
                  "lat":"0.0",
                  "longitude":"0.0",
                  "latitude":"0.0",
                  "wetdays":"0.0",
                  "ts":"273.15"
    }
    standard_factm = {"t2m":"1.0",
                      "pp":""+str(fac),
                      "swrs":"1.0",
                      "lon":"1.0",
                      "lat":"1.0",
                      "longitude":"1.0",
                      "latitude":"1.0",
                      "wetdays":"1.0",
                      "ts":"1.0"
    }

    try:
       return ( standard_names[nc_var_nm], standard_units[nc_var_nm],
                standard_facta[nc_var_nm], standard_factm[nc_var_nm] )

    except KeyError:
#      if nc_var_nm in ["bounds_latitude","bounds_longitude",]:
      # raise ValueError('Unkown variable ... '+str(nc_var_nm))
      return ( nc_var_nm, "", "0.0", "1.0" )


       # sys.exit()

#end def make_cf_compliant

def famed(woa_oxyg_dh,woa_oxyg_dh_m,temp_cl,temp_cl_m,depth):

    import numpy as np
    from numpy import ma

# Main FAME computation ...
# =========================

#	INPUTS FOR THE FAME CALCULATIONS
#		woa_oxyg_dh   == d18sw       in seawater <time mean> shape is [102, 180, 360]    e.g. depth, lat, lon
#		woa_oxyg_dh_m == d18sw       in seawater monthly     shape is [12,57, 180, 360]  e.g. depth, lat, lon
#		temp_cl       == temperature in seawater <time mean> shape is [1, 102, 180, 360] e.g. time (degenerate), depth, lat, lon
#		temp_cl_m     == temperature in seawater monthly     shape is [12, 57, 180, 360] e.g. time             , depth, lat, lon
#		depth         == depth of the levels in meters assumed positive here

    # FAME is coded in Kelvins ...
    temp_kl = temp_cl + 273.15
    temp_kl_m = temp_cl_m + 273.15

    # Computation of equilibrium calcite from WOA fields ...
    delt_dh_init = delta_c(temp_cl[0,...],woa_oxyg_dh)
    delt_dh_init_m = delta_c(temp_cl_m,woa_oxyg_dh_m)

    # i.e. auld method, d18Oca averaged over 50 meters ... == "Lukas Jonkers" methodology
    depth_50m = find_closest(depth,50.0)
    depth_00m = find_closest(depth, 0.0)

    # NOTA: all *_lj variables have a dimension without time and depth
    #       e.g. [180,360], lat, lon

    d18Osw_lj = ma.mean(woa_oxyg_dh[depth_00m:depth_50m,...],axis=0)
    tempcl_lj = ma.mean(temp_cl[0,depth_00m:depth_50m,...],axis=0)
    d18Oca_lj = delta_c(tempcl_lj,d18Osw_lj)

    d18Osw_ol = woa_oxyg_dh[depth_00m,...]
    tempcl_ol = temp_cl[0,depth_00m,...]
    d18Oca_ol = delta_c(tempcl_ol,d18Osw_ol)

    import forams_prod_l09 as fpl

    # Maximum shape of result: nb_forams, lat, lon
    max_shape_final = (len(fpl.l09_cnsts_dic),)+d18Osw_ol.shape

    # Create a placeholder for the foram result
    delt_forams = ma.zeros(max_shape_final,np.float32)

    for foram_specie in fpl.l09_cnsts_dic :

        # Rate of growth from the Lombard et al., 2009 methodology
        foram_growth = fpl.growth_rate_l09_array(foram_specie,temp_kl[0,...])
        foram_growth_m = fpl.growth_rate_l09_array(foram_specie,temp_kl_m)

        # Get the depth of the STD living foram in FAME
        f_dept = fpl.get_living_depth(foram_specie)

        # Find this depth as an index in the array water column
        indx_dfm = find_closest(depth,abs(float(f_dept[0])))

        # Shrink the FAME arrays to the foram living depth
        foram_growth = foram_growth[depth_00m:indx_dfm,...] # shape is depth, lat, lon
        foram_growth_m = foram_growth_m[:,depth_00m:indx_dfm,...] # shape is time, depth, lat, lon

        # Do the same for the equilibrium calcite from WOA
        delt_dh = delt_dh_init[depth_00m:indx_dfm,...] # idem

        delt_dh_m = delt_dh_init_m[:,depth_00m:indx_dfm,...] # idem

        # Get the location where there is SOME growth, based on a certain epsilon
        epsilon_growth = 0.1*fpl.l09_maxgrowth_dic[foram_specie][0] # or 0.032

        # Mask out the regions where the foram_growth is less than the epsilon
        masked_f_growth = ma.masked_less_equal(foram_growth,epsilon_growth)

        #~ nb_points_growth = (masked_f_growth * 0.0 + 1.0).filled(0.0)
        #~ if monthly is True:
           #~ nb_points_growth_m = (masked_f_growth_m * 0.0 + 1.0).filled(0.0)

        #~ nb_points_growth = ma.where(foram_growth > epsilon_growth,1,0) # 0.000001
        #~ if monthly is True:
           #~ nb_points_growth_m = ma.where(foram_growth_m > epsilon_growth,1,0) # 0.000001

        # Now sum the growth over the depth ...
        f_growth = ma.sum(masked_f_growth,axis=0)
        #~ n_growth = ma.where(ma.sum(nb_points_growth,axis=0)>0,1,0) # axis 0 = depth

        masked_f_growth_m = ma.masked_less_equal(foram_growth_m,epsilon_growth)
        f_growth_m = ma.sum(ma.sum(masked_f_growth_m,axis=1),axis=0)
        location_max_foramprod = ma.argmax(masked_f_growth_m,axis=1)
        location_max_foramprod = ma.masked_array(location_max_foramprod,mask=masked_f_growth_m[:,depth_00m:depth_50m,...].mask.all(axis=1))
        # location_max_foramprod = ma.masked_array(location_max_foramprod,mask=masked_f_growth_m[:,0,...].mask)

        # Computing the weighted sum for d18Ocalcite using growth over depth
        delt_fp = ma.sum(delt_dh*masked_f_growth,axis=0)
        delt_fp_m = ma.sum(ma.sum(delt_dh_m*masked_f_growth_m,axis=1),axis=0)

        # Mask out the points where no growth occur at all, in order to avoid NaNs ...
        delt_fp = delt_fp / ma.masked_less_equal(f_growth,0.0)
        delt_fp_m = delt_fp_m / ma.masked_less_equal(f_growth_m,0.0)

        # Result of FAME
        Z_om_fm = delt_fp
        Z_om_fm_m = ma.masked_array(delt_fp_m,mask=ma.max(location_max_foramprod[:,...],axis=0).mask)

        if foram_specie == "pachy_s":
            Z_om_fm = Z_om_fm + 0.1 # in per mil
            Z_om_fm_m = Z_om_fm_m + 0.1 # in per mil

        index_for = list(fpl.l09_cnsts_dic.keys()).index(foram_specie)
        delt_forams[index_for,...] = Z_om_fm_m

    #endfor on foram_specie

    # For comparison with Lukas Jonkers: old method on first 50 meters
    Z_om_lj = d18Oca_lj

    # For comparison with previous figures: old method on first 00 meters
    Z_om_ol = d18Oca_ol

    return delt_forams, Z_om_ol

#end def famed


def write_fame(nc_in, equi_calc,resultats_fame, nc_out="out-test.nc", latvar="lat",lonvar="lon",depthvar="depth"):

    import numpy as np
    import netCDF4
    import forams_prod_l09 as fpl

    # Write out the variables thus created

    dst = netCDF4.Dataset(nc_out,'w',format='NETCDF4_CLASSIC')
    dst.set_fill_on()

    # dmr --- use a modified solution from:
    #         CITE: http://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one
    with netCDF4.Dataset(nc_in,'r') as src:
        for nm_dim, dimension in src.dimensions.items():

            if nm_dim in [latvar,lonvar]:
               dst.createDimension(nm_dim,
                   len(dimension) if not dimension.isunlimited() else None)
               try:
                   v_dim = src.variables[nm_dim]
               except:
                   print("No data associated with variable", nm_dim)

               t_var = dst.createVariable(nm_dim,v_dim.datatype,v_dim.dimensions)

               stdnm, stdunit, stdadd, stdmul = make_cf_compliant(nm_dim)

               t_var.units = stdunit
               t_var.standard_name = stdnm

               if ( stdmul != "1.0" ):
                   temp = v_dim[:]*float(stdmul)
               else:
                   temp = v_dim[:]
               #endif
               if ( stdadd != "1.0" ):
                   t_var[:] = temp[:]+float(stdadd)
               else:
                   t_var[:] = temp[:]
               #endif
            #endif on nm_dim
        #endfor on nm_dim
    #endwith

    import time
    dst.history = 'Created ' + time.ctime(time.time()) \
                      + ' by FAME (v'+__version__+') (dmr,jyp,cw)'

    var_dimout = ()

    for key in ('lat','lon'):
        var_dimout += (dst.dimensions[key].name,)

    fill_value = netCDF4.default_fillvals["f4"]

    local_var = dst.createVariable(
                "d18Oc_std",np.float32().dtype,
                var_dimout, fill_value=fill_value
                                  )
    local_var.units = "per mil versus PDB"
    local_var.standard_name = "delta oxygen-18 of inorganic calcite"
    local_var[:] = equi_calc.filled(fill_value)


    for foram_species in fpl.l09_cnsts_dic :
    # for foram_species in ["dutertrei",] :
      local_var = dst.createVariable(
                  "d18Oc_"+foram_species,np.float32().dtype,
                  var_dimout, fill_value=fill_value
                                    )
      local_var.units = "per mil versus PDB"
      local_var.standard_name = "delta oxygen-18 of calcite in "+ foram_species

      indx_f = list(fpl.l09_cnsts_dic.keys()).index(foram_species)
      local_var[:] = resultats_fame[indx_f,...].filled(fill_value)

      #~ n_var = dst.createVariable(
                  #~ "nmg_"+foram_species,v_list[0].datatype,
                  #~ var_dimout, fill_value=fill_value
                                    #~ )
      #~ n_var.units = "month"
      #~ n_var.standard_name = "number months of non-zero growth "+ foram_species

      #~ n_var[:] = ngrowth_forams[indx_f,...].filled(fill_value)

      dst.sync()

    dst.close()

    return

#end def write_fame

# End of main FAME computation ...
# ================================

# The End of All Things (op. cit.)
