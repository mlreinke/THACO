device, static_gray=1
device,true_color=24,retain=2,decomposed=0  
loadct,12,/silent
!p.background=0
!path='/usr/local/cmod/codes/efit/idl/:'+!path
!path='/home/labombard/edge/analysis: '+!path	;used for Brian's awesome splined flux surface tools
!path='/home/labombard/edge/geometry: '+!path
!path='/home/labombard/edge/modeling: '+!path
!path='/home/labombard/idl_lib: '+!path

.compile /usr/local/cmod/idl/GENIE/mlr_functions.pro
.compile /usr/local/cmod/idl/GENIE/genie_help.pro
.compile /usr/local/cmod/idl/GENIE/bspline_fs.pro
.compile /usr/local/cmod/idl/GENIE/makesym.pro
.compile /usr/local/cmod/idl/GENIE/Trapow.pro
.compile /usr/local/cmod/idl/GENIE/Trapex.pro
.compile /usr/local/cmod/idl/GENIE/Trapez.pro
.compile /usr/local/cmod/idl/GENIE/trap_int.pro
.compile /usr/local/cmod/idl/GENIE/oploterror.pro
.compile /usr/local/cmod/idl/GENIE/mpfit.pro
.compile /usr/local/cmod/idl/GENIE/mpfitfun.pro
.compile /usr/local/cmod/idl/GENIE/tvimage.pro
.compile /usr/local/cmod/idl/GENIE/genie_line.pro
.compile /usr/local/cmod/idl/GENIE/read_atomicmass_tables.pro
.compile /usr/local/cmod/idl/GENIE/read_xray_data.pro
.compile /usr/local/cmod/idl/GENIE/GENPOS/genpos.pro
.compile /usr/local/cmod/idl/GENIE/GENPOS/genpos_sphere.pro
.compile /usr/local/cmod/idl/GENIE/GENPOS/genpos_utility.pro
.compile /usr/local/cmod/idl/GENIE/GENSPEC/genspec.pro
.compile /usr/local/cmod/idl/GENIE/GENSPEC/genspec_sphere.pro
.compile /usr/local/cmod/idl/GENIE/GENRAD/read_ark_table.pro

.compile /usr/local/cmod/idl/HIREXSR/hirexsr_fit_ellipse.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_load_data.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_fit_spectra.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_bin_spectra.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_calc_moments.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_calc_profiles.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_shot_analysis_tools.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_tree_utilities
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_zoom_kill.pro
.compile /usr/local/cmod/idl/HIREXSR/wbin_zoom_kill.pro
.compile /usr/local/cmod/idl/HIREXSR/hirexsr_bsfit.pro

.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_calib.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_det_align.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_binning.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_he_moments.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_moments.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_profiles.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_compare.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_tifit.pro
.compile /usr/local/cmod/idl/HIREXSR/w_hirexsr_omfit.pro
.compile /usr/local/cmod/idl/HIREXSR/w_thaco.pro

makesym,10
W_THACO,shot=1120224033,tht=0,line=2
