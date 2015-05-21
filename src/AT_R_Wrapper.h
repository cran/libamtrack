#ifndef AT_R_WRAPPER_H_
#define AT_R_WRAPPER_H_
// Automatically created header file

#include <stdlib.h>
#include <stdbool.h>

#include "AT_Algorithms_CPP.h"
#include "AT_Algorithms_GSM.h"
#include "AT_Algorithms_IGK.h"
#include "AT_DataMaterial.h"
#include "AT_DataParticle.h"
#include "AT_DataRange.h"
#include "AT_ElectronRange.h"
#include "AT_EnergyLoss.h"
#include "AT_PhysicsRoutines.h"
#include "AT_RDD.h"
#include "AT_SPC.h"
#include "AT_StoppingPower.h"
#include "AT_SuccessiveConvolutions.h"


void AT_run_CPPSC_method_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2_or_dose_Gy,
		const int*		material_no,
		const int*		stopping_power_source_no,
		const int*		rdd_model,
		const float*	rdd_parameters,
		const int*		er_model,
		const int*		gamma_model,
		const float*	gamma_parameters,
		int*			N2,
		const float*	fluence_factor,
		const int*		write_output,
		const int*		shrink_tails,
		const float*	shrink_tails_under,
		const int*		adjust_N2,
		const int*		lethal_events_mode,
		float*			relative_efficiency,
		float*			d_check,
		float*			S_HCP,
		float*			S_gamma,
		float*			mean_number_of_tracks_contrib,
		float*			start_number_of_tracks_contrib,
		int*			n_convolutions,
		float*			lower_Jensen_bound,
		float*			upper_Jensen_bound
);


void AT_run_GSM_method_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2_or_dose_Gy,
		const int*		material_no,
		const int*		stopping_power_source_no,
		const int*		rdd_model,
		const float*	rdd_parameters,
		const int*		er_model,
		const int*		gamma_model,
		const float*	gamma_parameters,
		const int*		N_runs,
		const int*		write_output,
		const int*		nX,
		const float*	voxel_size_m,
		const int*		lethal_events_mode,
		float*			relative_efficiency,
		float*			d_check,
		float*			S_HCP,
		float*			S_gamma,
		float*			n_particles,
		float*			sd_relative_efficiency,
		float*			sd_d_check,
		float*			sd_S_HCP,
		float*			sd_S_gamma,
		float*			sd_n_particles
);


void AT_run_IGK_method_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2_or_dose_Gy,
		const int*		material_no,
		const int*		stopping_power_source_no,
		const int*		rdd_model,
		const float*	rdd_parameters,
		const int*		er_model,
		const int*		gamma_model,
		const float*	gamma_parameters,
		const float*	saturation_cross_section_factor,
		const int*		write_output,
		float*			relative_efficiency,
		float*			S_HCP,
		float*			S_gamma,
		float*			sI_cm2,
		float*			gamma_dose_Gy,
		float*			P_I,
		float*			P_g
);


void AT_set_user_material_from_composition_R( const int*		n,
		const float*	density_g_cm3,
		const int*		A,
		const int*		Z,
		const float*	weight_fraction,
		int*			status
);


void AT_set_user_material_R( const float*	density_g_cm3,
		const float*	I_eV,
		const float*	average_A,
		const float*	average_Z,
		int*			status
);


void AT_I_eV_from_composition_R( const int*		n,
		const int*		Z,
		const int*		A,
		const float*	weight_fraction,
		float*			I_eV
);


void AT_effective_Z_from_composition_R( const int*		n,
		const int*		Z,
		const float*	weight_fraction,
		const float*	electron_densities_cm3,
		const float*	exponent,
		float*			effective_Z
);


void AT_average_Z_from_composition_R( const int*		n,
		const int*		Z,
		const float*	weight_fraction,
		float*			average_Z
);


void AT_average_A_from_composition_R( const int*		n,
		const int*		A,
		const float*	weight_fraction,
		float*			average_A
);


void AT_electron_density_m3_from_composition_R( const int*		n,
		const float*	density_g_cm3,
		const int*		Z,
		const int*		A,
		const float*	weight_fraction,
		float*			electron_density_m3
);


void AT_electron_density_m3_multi_R( const int*		n,
		const float*	density_g_cm3,
		const float*	average_Z,
		const float*	average_A,
		float*			electron_density_m3
);


void AT_electron_density_m3_from_material_no_multi_R( const int*		n,
		const int*		material_no,
		float*			electron_density_m3
);


void AT_get_materials_data_R( const int*		number_of_materials,
		const int*		material_no,
		float*			density_g_cm3,
		float*			I_eV,
		float*			alpha_g_cm2_MeV,
		float*			p_MeV,
		float*			m_g_cm2,
		float*			average_A,
		float*			average_Z
);


void AT_nuclear_spin_from_particle_no_multi_R( const int*		n,
		const int*		particle_no,
		float*			I,
		int*			returnValue
);


void AT_Z_from_particle_no_R( const int*		n,
		const int*		particle_no,
		int*			Z,
		int*			returnValue
);


void AT_atomic_weight_from_Z_R( const int*		n,
		const int*		Z,
		float*			atomic_weight,
		int*			returnValue
);


void AT_A_from_particle_no_R( const int*		n,
		const int*		particle_no,
		int*			A,
		int*			returnValue
);


void AT_particle_no_from_Z_and_A_R( const int*		n,
		const int*		Z,
		const int*		A,
		int*			particle_no,
		int*			returnValue
);


void AT_WEPL_Bethe_multi_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_m,
		float*			WEPL
);


void AT_CSDA_energy_after_slab_E_MeV_u_multi_R( const int*		n,
		const float*	E_initial_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_m,
		float*			E_final_MeV_u
);


void AT_CSDA_range_Bethe_g_cm2_multi_R( const int*		n,
		const float*	E_initial_MeV_u,
		const float*	E_final_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			CSDA_range_cm2_g
);


void AT_max_electron_ranges_m_R( const int*		number_of_particles,
		const float*	E_MeV_u,
		const int*		material_no,
		const int*		er_model,
		float*			max_electron_range_m
);


#ifdef HAVE_CERNLIB
void AT_Vavilov_PDF_R( const int*		n,
		const float*	lambda_V,
		const float*	kappa,
		const float*	beta,
		float*			density
);
#endif


#ifdef HAVE_CERNLIB
void AT_Landau_PDF_R( const int*		n,
		const float*	lambda,
		float*			density
);
#endif


void AT_Rutherford_SDCS_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const int*		n,
		const float*	T_MeV,
		float*			dsdT_m2_MeV,
		int*			returnValue
);


void AT_Bethe_mean_energy_loss_MeV_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			returnValue
);


void AT_kappa_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			returnValue
);


#ifdef HAVE_CERNLIB
void AT_lambda_from_energy_loss_R( const int*		n,
		const float*	energy_loss_keV,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			lambda_V
);
#endif


#ifdef HAVE_CERNLIB
void AT_Vavilov_energy_loss_distribution_R( const int*		n,
		const float*	energy_loss_keV,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			fDdD
);
#endif


#ifdef HAVE_CERNLIB
void AT_energy_loss_distribution_R( const int*		n,
		const float*	energy_loss_keV,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			fDdD
);
#endif


#ifdef HAVE_CERNLIB
void AT_energy_loss_mode_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			returnValue
);
#endif


#ifdef HAVE_CERNLIB
void AT_energy_loss_FWHM_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			returnValue
);
#endif


void AT_beta_from_E_R( const int*		n,
		const float*	E_MeV_u,
		float*			beta,
		int*			returnValue
);


void AT_E_from_beta_R( const int*		n,
		const float*	beta,
		float*			E_MeV_u,
		int*			returnValue
);


void AT_E_MeV_u_from_momentum_MeV_c_u_R( const int*		n,
		const float*	momentum_MeV_c_u,
		float*			E_MeV_u,
		int*			returnValue
);


void AT_gamma_from_E_R( const int*		n,
		const float*	E_MeV_u,
		float*			gamma,
		int*			returnValue
);


void AT_energy_straggling_MeV2_cm2_g_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			dsE2dz_MeV2_cm2_g
);


void AT_energy_straggling_after_slab_E_MeV_u_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_m,
		const float*	initial_sigma_E_MeV_u,
		float*			sigma_E_MeV_u
);


void AT_effective_charge_from_E_MeV_u_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		float*			effective_charge,
		int*			returnValue
);


void AT_max_E_transfer_MeV_R( const int*		n,
		const float*	E_MeV_u,
		float*			max_E_transfer_MeV,
		int*			returnValue
);


void AT_momentum_MeV_c_u_from_E_MeV_u_R( const int*		n,
		const float*	E_MeV_u,
		float*			momentum_MeV_c,
		int*			returnValue
);


void AT_dose_Gy_from_fluence_cm2_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			dose_Gy
);


void AT_fluence_cm2_from_dose_Gy_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	D_Gy,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			fluence_cm2
);


void AT_beam_par_physical_to_technical_R( const int*		n,
		const float*	fluence_cm2,
		const float*	sigma_cm,
		float*			N,
		float*			FWHM_mm
);


void AT_beam_par_technical_to_physical_R( const int*		n,
		const float*	N,
		const float*	FWHM_mm,
		float*			fluence_cm2,
		float*			sigma_cm
);


void AT_total_D_Gy_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
);


void AT_total_fluence_cm2_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	D_Gy,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
);


void AT_fluence_weighted_E_MeV_u_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const float*	fluence_cm2,
		float*			returnValue
);


void AT_dose_weighted_E_MeV_u_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
);


void AT_fluence_weighted_LET_MeV_cm2_g_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
);


void AT_dose_weighted_LET_MeV_cm2_g_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
);


void AT_stopping_power_ratio_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		reference_material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
);


void AT_mean_number_of_tracks_contrib_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		er_model,
		const int*		stopping_power_source_no,
		float*			returnValue
);


void AT_r_RDD_m_R( const int*		n,
		const float*	D_RDD_Gy,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const int*		rdd_model,
		const float*	rdd_parameter,
		const int*		er_model,
		const int*		stopping_power_source_no,
		float*			r_RDD_m,
		int*			returnValue
);


void AT_D_RDD_Gy_R( const int*		n,
		const float*	r_m,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const int*		rdd_model,
		const float*	rdd_parameter,
		const int*		er_model,
		const int*		stopping_power_source_no,
		float*			D_RDD_Gy,
		int*			returnValue
);


void AT_SPC_read_data_from_filename_fast_R( const char**	filename,
		int*			n,
		int*			depth_step,
		float*			depth_g_cm2,
		float*			E_MeV_u,
		float*			DE_MeV_u,
		int*			particle_no,
		float*			fluence_cm2,
		int*			returnValue
);


void AT_SPC_read_header_from_filename_fast_R( const char**	filename,
		float*			E_MeV_u,
		float*			peak_position_g_cm2,
		int*			particle_no,
		int*			material_no,
		float*			normalisation,
		int*			depth_steps_no,
		int*			returnValue
);


void AT_SPC_get_number_of_bins_from_filename_fast_R( const char**	filename,
		int*			returnValue
);


void AT_Mass_Stopping_Power_R( const char**	stopping_power_source,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_MeV_cm2_g,
		int*			returnValue
);


void AT_Stopping_Power_R( const char**	stopping_power_source,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_keV_um,
		int*			returnValue
);


void AT_Mass_Stopping_Power_with_no_R( const int*		stopping_power_source_no,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_MeV_cm2_g,
		int*			returnValue
);


void AT_Stopping_Power_with_no_R( const int*		stopping_power_source_no,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_keV_um,
		int*			returnValue
);


void AT_translate_dose_into_DSB_distribution_R( const int*		n_bins_f,
		const float*	f_d_Gy,
		const float*	f_dd_Gy,
		const float*	f,
		const float*	enhancement_factor,
		const float*	DSB_per_Gy_per_domain,
		const int*		domains_per_nucleus,
		const int*		write_output,
		float*			total_pDSBs,
		float*			total_nDSBs,
		float*			number_of_iDSBs,
		float*			number_of_cDSBs,
		float*			avg_number_of_DSBs_in_cDSBs
);


#endif

