allometry.o : $(ED_UTILS)/allometry.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

an_header.o: $(ED_IO)/an_header.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

average_utils.o : $(ED_IO)/average_utils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_LOWO_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

bdf2_solver.o : $(ED_DYNAMICS)/bdf2_solver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

budget_utils.o : $(ED_UTILS)/budget_utils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

canopy_air_coms.o : $(ED_MEMORY)/canopy_air_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

canopy_layer_coms.o : $(ED_MEMORY)/canopy_layer_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

canopy_radiation_coms.o : $(ED_MEMORY)/canopy_radiation_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

canopy_struct_dynamics.o : $(ED_DYNAMICS)/canopy_struct_dynamics.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

c34constants.o : $(ED_MEMORY)/c34constants.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

charutils.o: $(ED_UTILS)/charutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

consts_coms.o : $(ED_MEMORY)/consts_coms.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

dateutils.o: $(ED_UTILS)/dateutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

decomp_coms.o : $(ED_MEMORY)/decomp_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

detailed_coms.o : $(ED_MEMORY)/detailed_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

disturbance.o : $(ED_DYNAMICS)/disturbance.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

disturb_coms.o : $(ED_MEMORY)/disturb_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

edio.o : $(ED_IO)/edio.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_1st.o : $(ED_DRIVER)/ed_1st.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_bigleaf_init.o : $(ED_INIT)/ed_bigleaf_init.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_driver.o : $(ED_DRIVER)/ed_driver.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_filelist.o : $(ED_UTILS)/ed_filelist.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2) 

ed_grid.o : $(ED_UTILS)/ed_grid.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1) 

ed_init.o : $(ED_INIT)/ed_init.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_init_atm.o : $(ED_INIT)/ed_init_atm.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2) 

ed_init_full_history.o : $(ED_IO)/ed_init_full_history.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_LOWO_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_load_namelist.o : $(ED_IO)/ed_load_namelist.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_max_dims.o : $(ED_MEMORY)/ed_max_dims.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_mem_alloc.o : $(ED_MEMORY)/ed_mem_alloc.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1) 

ed_mem_grid_dim_defs.o : $(ED_MEMORY)/ed_mem_grid_dim_defs.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_met_driver.o : $(ED_DRIVER)/ed_met_driver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(HDF5_INCS) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_misc_coms.o : $(ED_MEMORY)/ed_misc_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_model.o : $(ED_DRIVER)/ed_model.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_mpass_init.o : $(ED_MPI)/ed_mpass_init.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_nbg_init.o : $(ED_INIT)/ed_nbg_init.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_node_coms.o : $(ED_MPI)/ed_node_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1) 

ed_opspec.o : $(ED_IO)/ed_opspec.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_params.o : $(ED_INIT)/ed_params.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_para_coms.o : $(ED_MPI)/ed_para_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1) 

ed_para_init.o : $(ED_MPI)/ed_para_init.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_print.o : $(ED_IO)/ed_print.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_read_ed10_20_history.o : $(ED_IO)/ed_read_ed10_20_history.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_read_ed21_history.o : $(ED_IO)/ed_read_ed21_history.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_state_vars.o : $(ED_MEMORY)/ed_state_vars.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_LOWO_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

ed_therm_lib.o : $(ED_UTILS)/ed_therm_lib.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_type_init.o : $(ED_INIT)/ed_type_init.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_var_tables.o : $(ED_MEMORY)/ed_var_tables.f90 
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_work_vars.o : $(ED_MEMORY)/ed_work_vars.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ed_xml_config.o : $(ED_IO)/ed_xml_config.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

ename_coms.o : $(ED_MEMORY)/ename_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

euler_driver.o : $(ED_DYNAMICS)/euler_driver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

events.o : $(ED_DYNAMICS)/events.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

farq_leuning.o : $(ED_DYNAMICS)/farq_leuning.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

fatal_error.o : $(ED_UTILS)/fatal_error.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

fire.o : $(ED_DYNAMICS)/fire.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

forestry.o : $(ED_DYNAMICS)/forestry.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

fusion_fission_coms.o : $(ED_MEMORY)/fusion_fission_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

fuse_fiss_utils.o : $(ED_UTILS)/fuse_fiss_utils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

great_circle.o : $(ED_UTILS)/great_circle.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1) 

grid_coms.o : $(ED_MEMORY)/grid_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1) 

growth_balive.o : $(ED_DYNAMICS)/growth_balive.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

h5_output.o : $(ED_IO)/h5_output.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

hdf5_coms.o : $(ED_MEMORY)/hdf5_coms.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

hdf5_utils.o : $(ED_UTILS)/hdf5_utils.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(HDF5_INCS) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

heun_driver.o: $(ED_DYNAMICS)/heun_driver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

hybrid_driver.o : $(ED_DYNAMICS)/hybrid_driver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

hydrology_coms.o: $(ED_MEMORY)/hydrology_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

hydrology_constants.o: $(ED_MEMORY)/hydrology_constants.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

init_hydro_sites.o : $(ED_INIT)/init_hydro_sites.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90) 
	$(RM_COMMAND_1)

invmondays.o : $(ED_UTILS)/invmondays.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

isotopes.o : $(ED_MEMORY)/isotopes.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

iso_alloc.o : $(ED_DYNAMICS)/iso_alloc.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

iso_utils.o : $(ED_UTILS)/iso_utils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

landuse_init.o : $(ED_INIT)/landuse_init.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

lapse.o : $(ED_UTILS)/lapse.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

leaf_database.o : $(ED_IO)/leaf_database.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

libxml2f90.f90_pp.o : $(ED_UTILS)/libxml2f90.f90_pp.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

lsm_hyd.o : $(ED_DYNAMICS)/lsm_hyd.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

mem_polygons.o : $(ED_MEMORY)/mem_polygons.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

met_driver_coms.o : $(ED_MEMORY)/met_driver_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

mortality.o : $(ED_DYNAMICS)/mortality.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

multiple_scatter.o : $(ED_DYNAMICS)/multiple_scatter.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

numutils.o: $(ED_UTILS)/numutils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

old_twostream_rad.o : $(ED_DYNAMICS)/old_twostream_rad.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

optimiz_coms.o : $(ED_MEMORY)/optimiz_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1) 

phenology_aux.o : $(ED_DYNAMICS)/phenology_aux.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

phenology_coms.o : $(ED_MEMORY)/phenology_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

phenology_driv.o : $(ED_DYNAMICS)/phenology_driv.f90
	cp -f $< $(<F:.f90=.f90) 
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

phenology_startup.o : $(ED_INIT)/phenology_startup.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

photosyn_driv.o : $(ED_DYNAMICS)/photosyn_driv.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

physiology_coms.o : $(ED_MEMORY)/physiology_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

pft_coms.o : $(ED_MEMORY)/pft_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

radiate_driver.o : $(ED_DYNAMICS)/radiate_driver.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

radiate_utils.o : $(ED_UTILS)/radiate_utils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

reproduction.o : $(ED_DYNAMICS)/reproduction.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

rk4_coms.o : $(ED_MEMORY)/rk4_coms.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

rk4_derivs.o : $(ED_DYNAMICS)/rk4_derivs.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

rk4_driver.o : $(ED_DYNAMICS)/rk4_driver.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

rk4_integ_utils.o : $(ED_DYNAMICS)/rk4_integ_utils.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

rk4_misc.o : $(ED_DYNAMICS)/rk4_misc.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

rk4_stepper.o : $(ED_DYNAMICS)/rk4_stepper.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

rsys.o: $(ED_UTILS)/rsys.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

soil_coms.o : $(ED_MEMORY)/soil_coms.F90
	cp -f $< $(<F:.F90=.F90)
	$(FPP_COMMAND) $(<F:.F90=.F90)
	$(RM_COMMAND_2)

soil_respiration.o : $(ED_DYNAMICS)/soil_respiration.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

stable_cohorts.o : $(ED_UTILS)/stable_cohorts.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

structural_growth.o : $(ED_DYNAMICS)/structural_growth.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

therm_lib.o: $(ED_UTILS)/therm_lib.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

therm_lib8.o: $(ED_UTILS)/therm_lib8.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

twostream_rad.o : $(ED_DYNAMICS)/twostream_rad.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

update_derived_props.o : $(ED_UTILS)/update_derived_props.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

utils_c.o: $(ED_UTILS)/utils_c.c
	cp -f $< $(<F:.c=.c)
	$(CXX_COMMAND) $<
	rm -f $(<F:.c=.c)

utils_f.o: $(ED_UTILS)/utils_f.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

vegetation_dynamics.o : $(ED_DYNAMICS)/vegetation_dynamics.f90
	cp -f $< $(<F:.f90=.f90)
	$(F90_COMMAND) $(<F:.f90=.f90)
	$(RM_COMMAND_1)

include dependency.mk