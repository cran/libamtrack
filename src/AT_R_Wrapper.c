// Automatically created header and body file

#include "AT_R_Wrapper.h"

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
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long gamma_model_long = (long)(*gamma_model);
  long N2_long = (long)(*N2);
  const double fluence_factor_double = (double)(*fluence_factor);
  const bool write_output_bool = (bool)(*write_output);
  const bool shrink_tails_bool = (bool)(*shrink_tails);
  const double shrink_tails_under_double = (double)(*shrink_tails_under);
  const bool adjust_N2_bool = (bool)(*adjust_N2);
  const bool lethal_events_mode_bool = (bool)(*lethal_events_mode);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_or_dose_Gy_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_or_dose_Gy_double[i] = (double)fluence_cm2_or_dose_Gy[i];
  }

//Allocate space for the input parameter.
  double* rdd_parameters_double = (double*)calloc(4,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 4; i++){
	rdd_parameters_double[i] = (double)rdd_parameters[i];
  }

//Allocate space for the input parameter.
  double* gamma_parameters_double = (double*)calloc(9,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 9; i++){
	gamma_parameters_double[i] = (double)gamma_parameters[i];
  }

//Define type-casted output variables
	double relative_efficiency_double = 0;
	double d_check_double = 0;
	double S_HCP_double = 0;
	double S_gamma_double = 0;
	double mean_number_of_tracks_contrib_double = 0;
	double start_number_of_tracks_contrib_double = 0;
	long n_convolutions_long = 0;
	double lower_Jensen_bound_double = 0;
	double upper_Jensen_bound_double = 0;

  AT_run_CPPSC_method( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_or_dose_Gy_double,
	material_no_long,
	stopping_power_source_no_long,
	rdd_model_long,
	rdd_parameters_double,
	er_model_long,
	gamma_model_long,
	gamma_parameters_double,
	N2_long,
	fluence_factor_double,
	write_output_bool,
	shrink_tails_bool,
	shrink_tails_under_double,
	adjust_N2_bool,
	lethal_events_mode_bool,
	&relative_efficiency_double,
	&d_check_double,
	&S_HCP_double,
	&S_gamma_double,
	&mean_number_of_tracks_contrib_double,
	&start_number_of_tracks_contrib_double,
	&n_convolutions_long,
	&lower_Jensen_bound_double,
	&upper_Jensen_bound_double);

//Results:
  *N2 = (int)N2_long;

  *relative_efficiency = (float)relative_efficiency_double;

  *d_check = (float)d_check_double;

  *S_HCP = (float)S_HCP_double;

  *S_gamma = (float)S_gamma_double;

  *mean_number_of_tracks_contrib = (float)mean_number_of_tracks_contrib_double;

  *start_number_of_tracks_contrib = (float)start_number_of_tracks_contrib_double;

  *n_convolutions = (int)n_convolutions_long;

  *lower_Jensen_bound = (float)lower_Jensen_bound_double;

  *upper_Jensen_bound = (float)upper_Jensen_bound_double;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_or_dose_Gy_double);
  free(rdd_parameters_double);
  free(gamma_parameters_double);
}



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
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long gamma_model_long = (long)(*gamma_model);
  const long N_runs_long = (long)(*N_runs);
  const bool write_output_bool = (bool)(*write_output);
  const long nX_long = (long)(*nX);
  const double voxel_size_m_double = (double)(*voxel_size_m);
  const bool lethal_events_mode_bool = (bool)(*lethal_events_mode);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_or_dose_Gy_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_or_dose_Gy_double[i] = (double)fluence_cm2_or_dose_Gy[i];
  }

//Allocate space for the input parameter.
  double* rdd_parameters_double = (double*)calloc(4,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 4; i++){
	rdd_parameters_double[i] = (double)rdd_parameters[i];
  }

//Allocate space for the input parameter.
  double* gamma_parameters_double = (double*)calloc(9,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 9; i++){
	gamma_parameters_double[i] = (double)gamma_parameters[i];
  }

//Define type-casted output variables
	double relative_efficiency_double = 0;
	double d_check_double = 0;
	double S_HCP_double = 0;
	double S_gamma_double = 0;
	double n_particles_double = 0;
	double sd_relative_efficiency_double = 0;
	double sd_d_check_double = 0;
	double sd_S_HCP_double = 0;
	double sd_S_gamma_double = 0;
	double sd_n_particles_double = 0;

  AT_run_GSM_method( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_or_dose_Gy_double,
	material_no_long,
	stopping_power_source_no_long,
	rdd_model_long,
	rdd_parameters_double,
	er_model_long,
	gamma_model_long,
	gamma_parameters_double,
	N_runs_long,
	write_output_bool,
	nX_long,
	voxel_size_m_double,
	lethal_events_mode_bool,
	&relative_efficiency_double,
	&d_check_double,
	&S_HCP_double,
	&S_gamma_double,
	&n_particles_double,
	&sd_relative_efficiency_double,
	&sd_d_check_double,
	&sd_S_HCP_double,
	&sd_S_gamma_double,
	&sd_n_particles_double);

//Results:
  *relative_efficiency = (float)relative_efficiency_double;

  *d_check = (float)d_check_double;

  *S_HCP = (float)S_HCP_double;

  *S_gamma = (float)S_gamma_double;

  *n_particles = (float)n_particles_double;

  *sd_relative_efficiency = (float)sd_relative_efficiency_double;

  *sd_d_check = (float)sd_d_check_double;

  *sd_S_HCP = (float)sd_S_HCP_double;

  *sd_S_gamma = (float)sd_S_gamma_double;

  *sd_n_particles = (float)sd_n_particles_double;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_or_dose_Gy_double);
  free(rdd_parameters_double);
  free(gamma_parameters_double);
}



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
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long gamma_model_long = (long)(*gamma_model);
  const double saturation_cross_section_factor_double = (double)(*saturation_cross_section_factor);
  const bool write_output_bool = (bool)(*write_output);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_or_dose_Gy_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_or_dose_Gy_double[i] = (double)fluence_cm2_or_dose_Gy[i];
  }

//Allocate space for the input parameter.
  double* rdd_parameters_double = (double*)calloc(4,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 4; i++){
	rdd_parameters_double[i] = (double)rdd_parameters[i];
  }

//Allocate space for the input parameter.
  double* gamma_parameters_double = (double*)calloc(9,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 9; i++){
	gamma_parameters_double[i] = (double)gamma_parameters[i];
  }

//Define type-casted output variables
	double relative_efficiency_double = 0;
	double S_HCP_double = 0;
	double S_gamma_double = 0;
	double sI_cm2_double = 0;
	double gamma_dose_Gy_double = 0;
	double P_I_double = 0;
	double P_g_double = 0;

  AT_run_IGK_method( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_or_dose_Gy_double,
	material_no_long,
	stopping_power_source_no_long,
	rdd_model_long,
	rdd_parameters_double,
	er_model_long,
	gamma_model_long,
	gamma_parameters_double,
	saturation_cross_section_factor_double,
	write_output_bool,
	&relative_efficiency_double,
	&S_HCP_double,
	&S_gamma_double,
	&sI_cm2_double,
	&gamma_dose_Gy_double,
	&P_I_double,
	&P_g_double);

//Results:
  *relative_efficiency = (float)relative_efficiency_double;

  *S_HCP = (float)S_HCP_double;

  *S_gamma = (float)S_gamma_double;

  *sI_cm2 = (float)sI_cm2_double;

  *gamma_dose_Gy = (float)gamma_dose_Gy_double;

  *P_I = (float)P_I_double;

  *P_g = (float)P_g_double;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_or_dose_Gy_double);
  free(rdd_parameters_double);
  free(gamma_parameters_double);
}



void AT_set_user_material_from_composition_R( const int*		n,
		const float*	density_g_cm3,
		const int*		A,
		const int*		Z,
		const float*	weight_fraction,
		int*			status
){
  long i;
  const long n_long = (long)(*n);
  const double density_g_cm3_double = (double)(*density_g_cm3);

//Allocate space for the input parameter.
  long* A_long = (long*)calloc(n_long,sizeof(long));
  long* Z_long = (long*)calloc(n_long,sizeof(long));
  double* weight_fraction_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	A_long[i] = (long)A[i];
	Z_long[i] = (long)Z[i];
	weight_fraction_double[i] = (double)weight_fraction[i];
  }

//Define type-casted output variables
	long status_long = 0;

  AT_set_user_material_from_composition( n_long,
	density_g_cm3_double,
	A_long,
	Z_long,
	weight_fraction_double,
	&status_long);

//Results:
  *status = (int)status_long;


//Free allocated space
  free(A_long);
  free(Z_long);
  free(weight_fraction_double);
}



void AT_set_user_material_R( const float*	density_g_cm3,
		const float*	I_eV,
		const float*	average_A,
		const float*	average_Z,
		int*			status
){
  const double density_g_cm3_double = (double)(*density_g_cm3);
  const double I_eV_double = (double)(*I_eV);
  const double average_A_double = (double)(*average_A);
  const double average_Z_double = (double)(*average_Z);

//Define type-casted output variables
	long status_long = 0;

  AT_set_user_material( density_g_cm3_double,
	I_eV_double,
	average_A_double,
	average_Z_double,
	&status_long);

//Results:
  *status = (int)status_long;


//Free allocated space
}



void AT_I_eV_from_composition_R( const int*		n,
		const int*		Z,
		const int*		A,
		const float*	weight_fraction,
		float*			I_eV
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* Z_long = (long*)calloc(n_long,sizeof(long));
  long* A_long = (long*)calloc(n_long,sizeof(long));
  double* weight_fraction_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	Z_long[i] = (long)Z[i];
	A_long[i] = (long)A[i];
	weight_fraction_double[i] = (double)weight_fraction[i];
  }

//Define type-casted output variables
	double I_eV_double = 0;

  AT_I_eV_from_composition( n_long,
	Z_long,
	A_long,
	weight_fraction_double,
	&I_eV_double);

//Results:
  *I_eV = (float)I_eV_double;


//Free allocated space
  free(Z_long);
  free(A_long);
  free(weight_fraction_double);
}



void AT_effective_Z_from_composition_R( const int*		n,
		const int*		Z,
		const float*	weight_fraction,
		const float*	electron_densities_cm3,
		const float*	exponent,
		float*			effective_Z
){
  long i;
  const long n_long = (long)(*n);
  const double exponent_double = (double)(*exponent);

//Allocate space for the input parameter.
  long* Z_long = (long*)calloc(n_long,sizeof(long));
  double* weight_fraction_double = (double*)calloc(n_long,sizeof(double));
  double* electron_densities_cm3_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	Z_long[i] = (long)Z[i];
	weight_fraction_double[i] = (double)weight_fraction[i];
	electron_densities_cm3_double[i] = (double)electron_densities_cm3[i];
  }

//Define type-casted output variables
	double effective_Z_double = 0;

  AT_effective_Z_from_composition( n_long,
	Z_long,
	weight_fraction_double,
	electron_densities_cm3_double,
	exponent_double,
	&effective_Z_double);

//Results:
  *effective_Z = (float)effective_Z_double;


//Free allocated space
  free(Z_long);
  free(weight_fraction_double);
  free(electron_densities_cm3_double);
}



void AT_average_Z_from_composition_R( const int*		n,
		const int*		Z,
		const float*	weight_fraction,
		float*			average_Z
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* Z_long = (long*)calloc(n_long,sizeof(long));
  double* weight_fraction_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	Z_long[i] = (long)Z[i];
	weight_fraction_double[i] = (double)weight_fraction[i];
  }

//Define type-casted output variables
	double average_Z_double = 0;

  AT_average_Z_from_composition( n_long,
	Z_long,
	weight_fraction_double,
	&average_Z_double);

//Results:
  *average_Z = (float)average_Z_double;


//Free allocated space
  free(Z_long);
  free(weight_fraction_double);
}



void AT_average_A_from_composition_R( const int*		n,
		const int*		A,
		const float*	weight_fraction,
		float*			average_A
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* A_long = (long*)calloc(n_long,sizeof(long));
  double* weight_fraction_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	A_long[i] = (long)A[i];
	weight_fraction_double[i] = (double)weight_fraction[i];
  }

//Define type-casted output variables
	double average_A_double = 0;

  AT_average_A_from_composition( n_long,
	A_long,
	weight_fraction_double,
	&average_A_double);

//Results:
  *average_A = (float)average_A_double;


//Free allocated space
  free(A_long);
  free(weight_fraction_double);
}



void AT_electron_density_m3_from_composition_R( const int*		n,
		const float*	density_g_cm3,
		const int*		Z,
		const int*		A,
		const float*	weight_fraction,
		float*			electron_density_m3
){
  long i;
  const long n_long = (long)(*n);
  const double density_g_cm3_double = (double)(*density_g_cm3);

//Allocate space for the input parameter.
  long* Z_long = (long*)calloc(n_long,sizeof(long));
  long* A_long = (long*)calloc(n_long,sizeof(long));
  double* weight_fraction_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	Z_long[i] = (long)Z[i];
	A_long[i] = (long)A[i];
	weight_fraction_double[i] = (double)weight_fraction[i];
  }

//Define type-casted output variables
	double electron_density_m3_double = 0;

  AT_electron_density_m3_from_composition( n_long,
	density_g_cm3_double,
	Z_long,
	A_long,
	weight_fraction_double,
	&electron_density_m3_double);

//Results:
  *electron_density_m3 = (float)electron_density_m3_double;


//Free allocated space
  free(Z_long);
  free(A_long);
  free(weight_fraction_double);
}



void AT_electron_density_m3_multi_R( const int*		n,
		const float*	density_g_cm3,
		const float*	average_Z,
		const float*	average_A,
		float*			electron_density_m3
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* density_g_cm3_double = (double*)calloc(n_long,sizeof(double));
  double* average_Z_double = (double*)calloc(n_long,sizeof(double));
  double* average_A_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	density_g_cm3_double[i] = (double)density_g_cm3[i];
	average_Z_double[i] = (double)average_Z[i];
	average_A_double[i] = (double)average_A[i];
  }

//Allocate space for the results.
  double* electron_density_m3_double = (double*)calloc(n_long,sizeof(double));

  AT_electron_density_m3_multi( n_long,
	density_g_cm3_double,
	average_Z_double,
	average_A_double,
	electron_density_m3_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	electron_density_m3[i] = (float)electron_density_m3_double[i];
  }

//Free allocated space
  free(density_g_cm3_double);
  free(average_Z_double);
  free(average_A_double);
  free(electron_density_m3_double);
}



void AT_electron_density_m3_from_material_no_multi_R( const int*		n,
		const int*		material_no,
		float*			electron_density_m3
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* material_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	material_no_long[i] = (long)material_no[i];
  }

//Allocate space for the results.
  double* electron_density_m3_double = (double*)calloc(n_long,sizeof(double));

  AT_electron_density_m3_from_material_no_multi( n_long,
	material_no_long,
	electron_density_m3_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	electron_density_m3[i] = (float)electron_density_m3_double[i];
  }

//Free allocated space
  free(material_no_long);
  free(electron_density_m3_double);
}



void AT_get_materials_data_R( const int*		number_of_materials,
		const int*		material_no,
		float*			density_g_cm3,
		float*			I_eV,
		float*			alpha_g_cm2_MeV,
		float*			p_MeV,
		float*			m_g_cm2,
		float*			average_A,
		float*			average_Z
){
  long i;
  const long number_of_materials_long = (long)(*number_of_materials);

//Allocate space for the input parameter.
  long* material_no_long = (long*)calloc(number_of_materials_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_materials_long; i++){
	material_no_long[i] = (long)material_no[i];
  }

//Allocate space for the results.
  double* density_g_cm3_double = (double*)calloc(number_of_materials_long,sizeof(double));
  double* I_eV_double = (double*)calloc(number_of_materials_long,sizeof(double));
  double* alpha_g_cm2_MeV_double = (double*)calloc(number_of_materials_long,sizeof(double));
  double* p_MeV_double = (double*)calloc(number_of_materials_long,sizeof(double));
  double* m_g_cm2_double = (double*)calloc(number_of_materials_long,sizeof(double));
  double* average_A_double = (double*)calloc(number_of_materials_long,sizeof(double));
  double* average_Z_double = (double*)calloc(number_of_materials_long,sizeof(double));

  AT_get_materials_data( number_of_materials_long,
	material_no_long,
	density_g_cm3_double,
	I_eV_double,
	alpha_g_cm2_MeV_double,
	p_MeV_double,
	m_g_cm2_double,
	average_A_double,
	average_Z_double);

//Results:
  for(i = 0 ; i < number_of_materials_long; i++){
	density_g_cm3[i] = (float)density_g_cm3_double[i];
	I_eV[i] = (float)I_eV_double[i];
	alpha_g_cm2_MeV[i] = (float)alpha_g_cm2_MeV_double[i];
	p_MeV[i] = (float)p_MeV_double[i];
	m_g_cm2[i] = (float)m_g_cm2_double[i];
	average_A[i] = (float)average_A_double[i];
	average_Z[i] = (float)average_Z_double[i];
  }

//Free allocated space
  free(material_no_long);
  free(density_g_cm3_double);
  free(I_eV_double);
  free(alpha_g_cm2_MeV_double);
  free(p_MeV_double);
  free(m_g_cm2_double);
  free(average_A_double);
  free(average_Z_double);
}



void AT_nuclear_spin_from_particle_no_multi_R( const int*		n,
		const int*		particle_no,
		float*			I,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* I_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_nuclear_spin_from_particle_no_multi( n_long,
	particle_no_long,
	I_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	I[i] = (float)I_double[i];
  }

//Free allocated space
  free(particle_no_long);
  free(I_double);
}



void AT_Z_from_particle_no_R( const int*		n,
		const int*		particle_no,
		int*			Z,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  long* Z_long = (long*)calloc(n_long,sizeof(long));

  int returnValue_internal = 	AT_Z_from_particle_no( n_long,
	particle_no_long,
	Z_long);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	Z[i] = (int)Z_long[i];
  }

//Free allocated space
  free(particle_no_long);
  free(Z_long);
}



void AT_atomic_weight_from_Z_R( const int*		n,
		const int*		Z,
		float*			atomic_weight,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* Z_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	Z_long[i] = (long)Z[i];
  }

//Allocate space for the results.
  double* atomic_weight_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_atomic_weight_from_Z( n_long,
	Z_long,
	atomic_weight_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	atomic_weight[i] = (float)atomic_weight_double[i];
  }

//Free allocated space
  free(Z_long);
  free(atomic_weight_double);
}



void AT_A_from_particle_no_R( const int*		n,
		const int*		particle_no,
		int*			A,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  long* A_long = (long*)calloc(n_long,sizeof(long));

  int returnValue_internal = 	AT_A_from_particle_no( n_long,
	particle_no_long,
	A_long);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	A[i] = (int)A_long[i];
  }

//Free allocated space
  free(particle_no_long);
  free(A_long);
}



void AT_particle_no_from_Z_and_A_R( const int*		n,
		const int*		Z,
		const int*		A,
		int*			particle_no,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  long* Z_long = (long*)calloc(n_long,sizeof(long));
  long* A_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	Z_long[i] = (long)Z[i];
	A_long[i] = (long)A[i];
  }

//Allocate space for the results.
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));

  int returnValue_internal = 	AT_particle_no_from_Z_and_A( n_long,
	Z_long,
	A_long,
	particle_no_long);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	particle_no[i] = (int)particle_no_long[i];
  }

//Free allocated space
  free(Z_long);
  free(A_long);
  free(particle_no_long);
}



void AT_CSDA_range_g_cm2_multi_R( const int*		n,
		const float*	E_initial_MeV_u,
		const float*	E_final_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			CSDA_range_cm2_g
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* E_initial_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  double* E_final_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_initial_MeV_u_double[i] = (double)E_initial_MeV_u[i];
	E_final_MeV_u_double[i] = (double)E_final_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* CSDA_range_cm2_g_double = (double*)calloc(n_long,sizeof(double));

  AT_CSDA_range_g_cm2_multi( n_long,
	E_initial_MeV_u_double,
	E_final_MeV_u_double,
	particle_no_long,
	material_no_long,
	CSDA_range_cm2_g_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	CSDA_range_cm2_g[i] = (float)CSDA_range_cm2_g_double[i];
  }

//Free allocated space
  free(E_initial_MeV_u_double);
  free(E_final_MeV_u_double);
  free(particle_no_long);
  free(CSDA_range_cm2_g_double);
}



void AT_max_electron_ranges_m_R( const int*		number_of_particles,
		const float*	E_MeV_u,
		const int*		material_no,
		const int*		er_model,
		float*			max_electron_range_m
){
  long i;
  const long number_of_particles_long = (long)(*number_of_particles);
  const int material_no_long = (int)(*material_no);
  const int er_model_long = (int)(*er_model);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_particles_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_particles_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

//Allocate space for the results.
  double* max_electron_range_m_double = (double*)calloc(number_of_particles_long,sizeof(double));

  AT_max_electron_ranges_m( number_of_particles_long,
	E_MeV_u_double,
	material_no_long,
	er_model_long,
	max_electron_range_m_double);

//Results:
  for(i = 0 ; i < number_of_particles_long; i++){
	max_electron_range_m[i] = (float)max_electron_range_m_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(max_electron_range_m_double);
}



void AT_mean_energy_loss_keV_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			returnValue
){
  const double E_MeV_u_double = (double)(*E_MeV_u);
  const long particle_no_long = (long)(*particle_no);
  const long material_no_long = (long)(*material_no);
  const double slab_thickness_um_double = (double)(*slab_thickness_um);

  double returnValue_internal = 	AT_mean_energy_loss_keV( E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
}



void AT_xi_keV_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			returnValue
){
  const double E_MeV_u_double = (double)(*E_MeV_u);
  const long particle_no_long = (long)(*particle_no);
  const long material_no_long = (long)(*material_no);
  const double slab_thickness_um_double = (double)(*slab_thickness_um);

  double returnValue_internal = 	AT_xi_keV( E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
}



void AT_kappa_multi_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			kappa
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* slab_thickness_um_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	slab_thickness_um_double[i] = (double)slab_thickness_um[i];
  }

//Allocate space for the results.
  double* kappa_double = (double*)calloc(n_long,sizeof(double));

  AT_kappa_multi( n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	kappa_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	kappa[i] = (float)kappa_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(slab_thickness_um_double);
  free(kappa_double);
}



void AT_Landau_PDF_R( const int*		n,
		const float*	lambda_landau,
		float*			density
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* lambda_landau_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	lambda_landau_double[i] = (double)lambda_landau[i];
  }

//Allocate space for the results.
  double* density_double = (double*)calloc(n_long,sizeof(double));

  AT_Landau_PDF( n_long,
	lambda_landau_double,
	density_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	density[i] = (float)density_double[i];
  }

//Free allocated space
  free(lambda_landau_double);
  free(density_double);
}



void AT_Landau_IDF_R( const int*		n,
		const float*	rnd,
		float*			lambda_landau
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* rnd_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	rnd_double[i] = (double)rnd[i];
  }

//Allocate space for the results.
  double* lambda_landau_double = (double*)calloc(n_long,sizeof(double));

  AT_Landau_IDF( n_long,
	rnd_double,
	lambda_landau_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	lambda_landau[i] = (float)lambda_landau_double[i];
  }

//Free allocated space
  free(rnd_double);
  free(lambda_landau_double);
}



void AT_lambda_landau_from_energy_loss_multi_R( const int*		n,
		const float*	energy_loss_keV,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			lambda_landau
){
  long i;
  const long n_long = (long)(*n);
  const double E_MeV_u_double = (double)(*E_MeV_u);
  const long particle_no_long = (long)(*particle_no);
  const long material_no_long = (long)(*material_no);
  const double slab_thickness_um_double = (double)(*slab_thickness_um);

//Allocate space for the input parameter.
  double* energy_loss_keV_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	energy_loss_keV_double[i] = (double)energy_loss_keV[i];
  }

//Allocate space for the results.
  double* lambda_landau_double = (double*)calloc(n_long,sizeof(double));

  AT_lambda_landau_from_energy_loss_multi( n_long,
	energy_loss_keV_double,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	lambda_landau_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	lambda_landau[i] = (float)lambda_landau_double[i];
  }

//Free allocated space
  free(energy_loss_keV_double);
  free(lambda_landau_double);
}



void AT_lambda_mean_multi_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			lambda_mean
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* slab_thickness_um_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	slab_thickness_um_double[i] = (double)slab_thickness_um[i];
  }

//Allocate space for the results.
  double* lambda_mean_double = (double*)calloc(n_long,sizeof(double));

  AT_lambda_mean_multi( n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	lambda_mean_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	lambda_mean[i] = (float)lambda_mean_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(slab_thickness_um_double);
  free(lambda_mean_double);
}



void AT_lambda_max_multi_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			lambda_max
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* slab_thickness_um_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	slab_thickness_um_double[i] = (double)slab_thickness_um[i];
  }

//Allocate space for the results.
  double* lambda_max_double = (double*)calloc(n_long,sizeof(double));

  AT_lambda_max_multi( n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	lambda_max_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	lambda_max[i] = (float)lambda_max_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(slab_thickness_um_double);
  free(lambda_max_double);
}



void AT_energy_loss_from_lambda_landau_multi_R( const int*		n,
		const float*	lambda_landau,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			energy_loss_keV
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* lambda_landau_double = (double*)calloc(n_long,sizeof(double));
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* slab_thickness_um_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	lambda_landau_double[i] = (double)lambda_landau[i];
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	slab_thickness_um_double[i] = (double)slab_thickness_um[i];
  }

//Allocate space for the results.
  double* energy_loss_keV_double = (double*)calloc(n_long,sizeof(double));

  AT_energy_loss_from_lambda_landau_multi( n_long,
	lambda_landau_double,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	energy_loss_keV_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	energy_loss_keV[i] = (float)energy_loss_keV_double[i];
  }

//Free allocated space
  free(lambda_landau_double);
  free(E_MeV_u_double);
  free(particle_no_long);
  free(slab_thickness_um_double);
  free(energy_loss_keV_double);
}



void AT_Vavilov_PDF_R( const int*		n,
		const float*	lambda_vavilov,
		const float*	kappa,
		const float*	beta,
		float*			density
){
  long i;
  const long n_long = (long)(*n);
  const double kappa_double = (double)(*kappa);
  const double beta_double = (double)(*beta);

//Allocate space for the input parameter.
  double* lambda_vavilov_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	lambda_vavilov_double[i] = (double)lambda_vavilov[i];
  }

//Allocate space for the results.
  double* density_double = (double*)calloc(n_long,sizeof(double));

  AT_Vavilov_PDF( n_long,
	lambda_vavilov_double,
	kappa_double,
	beta_double,
	density_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	density[i] = (float)density_double[i];
  }

//Free allocated space
  free(lambda_vavilov_double);
  free(density_double);
}



void AT_Vavilov_IDF_R( const int*		n,
		const float*	rnd,
		const float*	kappa,
		const float*	beta,
		float*			lambda_vavilov
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* rnd_double = (double*)calloc(n_long,sizeof(double));
  double* kappa_double = (double*)calloc(n_long,sizeof(double));
  double* beta_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	rnd_double[i] = (double)rnd[i];
	kappa_double[i] = (double)kappa[i];
	beta_double[i] = (double)beta[i];
  }

//Allocate space for the results.
  double* lambda_vavilov_double = (double*)calloc(n_long,sizeof(double));

  AT_Vavilov_IDF( n_long,
	rnd_double,
	kappa_double,
	beta_double,
	lambda_vavilov_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	lambda_vavilov[i] = (float)lambda_vavilov_double[i];
  }

//Free allocated space
  free(rnd_double);
  free(kappa_double);
  free(beta_double);
  free(lambda_vavilov_double);
}



void AT_lambda_vavilov_from_energy_loss_multi_R( const int*		n,
		const float*	energy_loss_keV,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			lambda_vavilov
){
  long i;
  const long n_long = (long)(*n);
  const double E_MeV_u_double = (double)(*E_MeV_u);
  const long particle_no_long = (long)(*particle_no);
  const long material_no_long = (long)(*material_no);
  const double slab_thickness_um_double = (double)(*slab_thickness_um);

//Allocate space for the input parameter.
  double* energy_loss_keV_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	energy_loss_keV_double[i] = (double)energy_loss_keV[i];
  }

//Allocate space for the results.
  double* lambda_vavilov_double = (double*)calloc(n_long,sizeof(double));

  AT_lambda_vavilov_from_energy_loss_multi( n_long,
	energy_loss_keV_double,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	lambda_vavilov_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	lambda_vavilov[i] = (float)lambda_vavilov_double[i];
  }

//Free allocated space
  free(energy_loss_keV_double);
  free(lambda_vavilov_double);
}



void AT_energy_loss_from_lambda_vavilov_multi_R( const int*		n,
		const float*	lambda_vavilov,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			energy_loss_keV
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* lambda_vavilov_double = (double*)calloc(n_long,sizeof(double));
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* slab_thickness_um_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	lambda_vavilov_double[i] = (double)lambda_vavilov[i];
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	slab_thickness_um_double[i] = (double)slab_thickness_um[i];
  }

//Allocate space for the results.
  double* energy_loss_keV_double = (double*)calloc(n_long,sizeof(double));

  AT_energy_loss_from_lambda_vavilov_multi( n_long,
	lambda_vavilov_double,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	energy_loss_keV_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	energy_loss_keV[i] = (float)energy_loss_keV_double[i];
  }

//Free allocated space
  free(lambda_vavilov_double);
  free(E_MeV_u_double);
  free(particle_no_long);
  free(slab_thickness_um_double);
  free(energy_loss_keV_double);
}



void AT_Gauss_PDF_R( const int*		n,
		const float*	lambda_gauss,
		float*			density
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* lambda_gauss_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	lambda_gauss_double[i] = (double)lambda_gauss[i];
  }

//Allocate space for the results.
  double* density_double = (double*)calloc(n_long,sizeof(double));

  AT_Gauss_PDF( n_long,
	lambda_gauss_double,
	density_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	density[i] = (float)density_double[i];
  }

//Free allocated space
  free(lambda_gauss_double);
  free(density_double);
}



void AT_Gauss_IDF_R( const int*		n,
		const float*	rnd,
		float*			lambda_gauss
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* rnd_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	rnd_double[i] = (double)rnd[i];
  }

//Allocate space for the results.
  double* lambda_gauss_double = (double*)calloc(n_long,sizeof(double));

  AT_Gauss_IDF( n_long,
	rnd_double,
	lambda_gauss_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	lambda_gauss[i] = (float)lambda_gauss_double[i];
  }

//Free allocated space
  free(rnd_double);
  free(lambda_gauss_double);
}



void AT_energy_loss_from_lambda_gauss_multi_R( const int*		n,
		const float*	lambda_gauss,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_um,
		float*			energy_loss_keV
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* lambda_gauss_double = (double*)calloc(n_long,sizeof(double));
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* slab_thickness_um_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	lambda_gauss_double[i] = (double)lambda_gauss[i];
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	slab_thickness_um_double[i] = (double)slab_thickness_um[i];
  }

//Allocate space for the results.
  double* energy_loss_keV_double = (double*)calloc(n_long,sizeof(double));

  AT_energy_loss_from_lambda_gauss_multi( n_long,
	lambda_gauss_double,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_um_double,
	energy_loss_keV_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	energy_loss_keV[i] = (float)energy_loss_keV_double[i];
  }

//Free allocated space
  free(lambda_gauss_double);
  free(E_MeV_u_double);
  free(particle_no_long);
  free(slab_thickness_um_double);
  free(energy_loss_keV_double);
}



void AT_beta_from_E_R( const int*		n,
		const float*	E_MeV_u,
		float*			beta,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

//Allocate space for the results.
  double* beta_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_beta_from_E( n_long,
	E_MeV_u_double,
	beta_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	beta[i] = (float)beta_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(beta_double);
}



void AT_E_from_beta_R( const int*		n,
		const float*	beta,
		float*			E_MeV_u,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* beta_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	beta_double[i] = (double)beta[i];
  }

//Allocate space for the results.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_E_from_beta( n_long,
	beta_double,
	E_MeV_u_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	E_MeV_u[i] = (float)E_MeV_u_double[i];
  }

//Free allocated space
  free(beta_double);
  free(E_MeV_u_double);
}



void AT_E_MeV_u_from_momentum_MeV_c_u_R( const int*		n,
		const float*	momentum_MeV_c_u,
		float*			E_MeV_u,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* momentum_MeV_c_u_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	momentum_MeV_c_u_double[i] = (double)momentum_MeV_c_u[i];
  }

//Allocate space for the results.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_E_MeV_u_from_momentum_MeV_c_u( n_long,
	momentum_MeV_c_u_double,
	E_MeV_u_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	E_MeV_u[i] = (float)E_MeV_u_double[i];
  }

//Free allocated space
  free(momentum_MeV_c_u_double);
  free(E_MeV_u_double);
}



void AT_gamma_from_E_R( const int*		n,
		const float*	E_MeV_u,
		float*			gamma,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

//Allocate space for the results.
  double* gamma_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_gamma_from_E( n_long,
	E_MeV_u_double,
	gamma_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	gamma[i] = (float)gamma_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(gamma_double);
}



void AT_energy_straggling_MeV2_cm2_g_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			dsE2dz_MeV2_cm2_g
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* dsE2dz_MeV2_cm2_g_double = (double*)calloc(n_long,sizeof(double));

  AT_energy_straggling_MeV2_cm2_g( n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	dsE2dz_MeV2_cm2_g_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	dsE2dz_MeV2_cm2_g[i] = (float)dsE2dz_MeV2_cm2_g_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(dsE2dz_MeV2_cm2_g_double);
}



void AT_energy_straggling_after_slab_E_MeV_u_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const float*	slab_thickness_m,
		const float*	initial_sigma_E_MeV_u,
		float*			sigma_E_MeV_u
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const double slab_thickness_m_double = (double)(*slab_thickness_m);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* initial_sigma_E_MeV_u_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	initial_sigma_E_MeV_u_double[i] = (double)initial_sigma_E_MeV_u[i];
  }

//Allocate space for the results.
  double* sigma_E_MeV_u_double = (double*)calloc(n_long,sizeof(double));

  AT_energy_straggling_after_slab_E_MeV_u( n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	slab_thickness_m_double,
	initial_sigma_E_MeV_u_double,
	sigma_E_MeV_u_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	sigma_E_MeV_u[i] = (float)sigma_E_MeV_u_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(initial_sigma_E_MeV_u_double);
  free(sigma_E_MeV_u_double);
}



void AT_effective_charge_from_E_MeV_u_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		float*			effective_charge,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* effective_charge_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_effective_charge_from_E_MeV_u( n_long,
	E_MeV_u_double,
	particle_no_long,
	effective_charge_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	effective_charge[i] = (float)effective_charge_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(effective_charge_double);
}



void AT_max_E_transfer_MeV_new_R( const int*		n,
		const float*	E_MeV_u,
		const int*		A,
		float*			max_E_transfer_MeV,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* A_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	A_long[i] = (long)A[i];
  }

//Allocate space for the results.
  double* max_E_transfer_MeV_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_max_E_transfer_MeV_new( n_long,
	E_MeV_u_double,
	A_long,
	max_E_transfer_MeV_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	max_E_transfer_MeV[i] = (float)max_E_transfer_MeV_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(A_long);
  free(max_E_transfer_MeV_double);
}



void AT_max_E_transfer_MeV_R( const int*		n,
		const float*	E_MeV_u,
		float*			max_E_transfer_MeV,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

//Allocate space for the results.
  double* max_E_transfer_MeV_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_max_E_transfer_MeV( n_long,
	E_MeV_u_double,
	max_E_transfer_MeV_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	max_E_transfer_MeV[i] = (float)max_E_transfer_MeV_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(max_E_transfer_MeV_double);
}



void AT_momentum_MeV_c_u_from_E_MeV_u_R( const int*		n,
		const float*	E_MeV_u,
		float*			momentum_MeV_c,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

//Allocate space for the results.
  double* momentum_MeV_c_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_momentum_MeV_c_u_from_E_MeV_u( n_long,
	E_MeV_u_double,
	momentum_MeV_c_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	momentum_MeV_c[i] = (float)momentum_MeV_c_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(momentum_MeV_c_double);
}



void AT_dose_Gy_from_fluence_cm2_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			dose_Gy
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

//Allocate space for the results.
  double* dose_Gy_double = (double*)calloc(n_long,sizeof(double));

  AT_dose_Gy_from_fluence_cm2( n_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_double,
	material_no_long,
	stopping_power_source_no_long,
	dose_Gy_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	dose_Gy[i] = (float)dose_Gy_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
  free(dose_Gy_double);
}



void AT_fluence_cm2_from_dose_Gy_R( const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	D_Gy,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			fluence_cm2
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* D_Gy_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	D_Gy_double[i] = (double)D_Gy[i];
  }

//Allocate space for the results.
  double* fluence_cm2_double = (double*)calloc(n_long,sizeof(double));

  AT_fluence_cm2_from_dose_Gy( n_long,
	E_MeV_u_double,
	particle_no_long,
	D_Gy_double,
	material_no_long,
	stopping_power_source_no_long,
	fluence_cm2_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	fluence_cm2[i] = (float)fluence_cm2_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(D_Gy_double);
  free(fluence_cm2_double);
}



void AT_beam_par_physical_to_technical_R( const int*		n,
		const float*	fluence_cm2,
		const float*	sigma_cm,
		float*			N,
		float*			FWHM_mm
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* fluence_cm2_double = (double*)calloc(n_long,sizeof(double));
  double* sigma_cm_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	fluence_cm2_double[i] = (double)fluence_cm2[i];
	sigma_cm_double[i] = (double)sigma_cm[i];
  }

//Allocate space for the results.
  double* N_double = (double*)calloc(n_long,sizeof(double));
  double* FWHM_mm_double = (double*)calloc(n_long,sizeof(double));

  AT_beam_par_physical_to_technical( n_long,
	fluence_cm2_double,
	sigma_cm_double,
	N_double,
	FWHM_mm_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	N[i] = (float)N_double[i];
	FWHM_mm[i] = (float)FWHM_mm_double[i];
  }

//Free allocated space
  free(fluence_cm2_double);
  free(sigma_cm_double);
  free(N_double);
  free(FWHM_mm_double);
}



void AT_beam_par_technical_to_physical_R( const int*		n,
		const float*	N,
		const float*	FWHM_mm,
		float*			fluence_cm2,
		float*			sigma_cm
){
  long i;
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* N_double = (double*)calloc(n_long,sizeof(double));
  double* FWHM_mm_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	N_double[i] = (double)N[i];
	FWHM_mm_double[i] = (double)FWHM_mm[i];
  }

//Allocate space for the results.
  double* fluence_cm2_double = (double*)calloc(n_long,sizeof(double));
  double* sigma_cm_double = (double*)calloc(n_long,sizeof(double));

  AT_beam_par_technical_to_physical( n_long,
	N_double,
	FWHM_mm_double,
	fluence_cm2_double,
	sigma_cm_double);

//Results:
  for(i = 0 ; i < n_long; i++){
	fluence_cm2[i] = (float)fluence_cm2_double[i];
	sigma_cm[i] = (float)sigma_cm_double[i];
  }

//Free allocated space
  free(N_double);
  free(FWHM_mm_double);
  free(fluence_cm2_double);
  free(sigma_cm_double);
}



void AT_total_D_Gy_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

  double returnValue_internal = 	AT_total_D_Gy( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_double,
	material_no_long,
	stopping_power_source_no_long);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
}



void AT_total_fluence_cm2_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	D_Gy,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* D_Gy_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	D_Gy_double[i] = (double)D_Gy[i];
  }

  double returnValue_internal = 	AT_total_fluence_cm2( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	D_Gy_double,
	material_no_long,
	stopping_power_source_no_long);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(D_Gy_double);
}



void AT_fluence_weighted_E_MeV_u_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const float*	fluence_cm2,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  double* fluence_cm2_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

  double returnValue_internal = 	AT_fluence_weighted_E_MeV_u( number_of_field_components_long,
	E_MeV_u_double,
	fluence_cm2_double);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(fluence_cm2_double);
}



void AT_dose_weighted_E_MeV_u_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

  double returnValue_internal = 	AT_dose_weighted_E_MeV_u( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_double,
	material_no_long,
	stopping_power_source_no_long);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
}



void AT_fluence_weighted_LET_MeV_cm2_g_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

  double returnValue_internal = 	AT_fluence_weighted_LET_MeV_cm2_g( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_double,
	material_no_long,
	stopping_power_source_no_long);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
}



void AT_dose_weighted_LET_MeV_cm2_g_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

  double returnValue_internal = 	AT_dose_weighted_LET_MeV_cm2_g( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_double,
	material_no_long,
	stopping_power_source_no_long);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
}



void AT_stopping_power_ratio_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		reference_material_no,
		const int*		stopping_power_source_no,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long reference_material_no_long = (long)(*reference_material_no);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

  double returnValue_internal = 	AT_stopping_power_ratio( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_double,
	material_no_long,
	reference_material_no_long,
	stopping_power_source_no_long);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
}



void AT_mean_number_of_tracks_contrib_R( const int*		number_of_field_components,
		const float*	E_MeV_u,
		const int*		particle_no,
		const float*	fluence_cm2,
		const int*		material_no,
		const int*		er_model,
		const int*		stopping_power_source_no,
		float*			returnValue
){
  long i;
  const long number_of_field_components_long = (long)(*number_of_field_components);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(number_of_field_components_long,sizeof(double));
  long* particle_no_long = (long*)calloc(number_of_field_components_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(number_of_field_components_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < number_of_field_components_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
	fluence_cm2_double[i] = (double)fluence_cm2[i];
  }

  double returnValue_internal = 	AT_mean_number_of_tracks_contrib( number_of_field_components_long,
	E_MeV_u_double,
	particle_no_long,
	fluence_cm2_double,
	material_no_long,
	er_model_long,
	stopping_power_source_no_long);

//Results:

	*returnValue = ( float )returnValue_internal;


//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
}



void AT_Rutherford_SDCS_R( const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		const int*		n,
		const float*	T_MeV,
		float*			dsdT_m2_MeV,
		int*			returnValue
){
  long i;
  const double E_MeV_u_double = (double)(*E_MeV_u);
  const long particle_no_long = (long)(*particle_no);
  const long material_no_long = (long)(*material_no);
  const long n_long = (long)(*n);

//Allocate space for the input parameter.
  double* T_MeV_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	T_MeV_double[i] = (double)T_MeV[i];
  }

//Allocate space for the results.
  double* dsdT_m2_MeV_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_Rutherford_SDCS( E_MeV_u_double,
	particle_no_long,
	material_no_long,
	n_long,
	T_MeV_double,
	dsdT_m2_MeV_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	dsdT_m2_MeV[i] = (float)dsdT_m2_MeV_double[i];
  }

//Free allocated space
  free(T_MeV_double);
  free(dsdT_m2_MeV_double);
}



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
){
  long i;
  const long n_long = (long)(*n);
  const double E_MeV_u_double = (double)(*E_MeV_u);
  const long particle_no_long = (long)(*particle_no);
  const long material_no_long = (long)(*material_no);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* D_RDD_Gy_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	D_RDD_Gy_double[i] = (double)D_RDD_Gy[i];
  }

//Allocate space for the results.
  double* r_RDD_m_double = (double*)calloc(n_long,sizeof(double));

//Allocate space for the input parameter.
  double* rdd_parameter_double = (double*)calloc(4,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 4; i++){
	rdd_parameter_double[i] = (double)rdd_parameter[i];
  }

  int returnValue_internal = 	AT_r_RDD_m( n_long,
	D_RDD_Gy_double,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	rdd_model_long,
	rdd_parameter_double,
	er_model_long,
	stopping_power_source_no_long,
	r_RDD_m_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	r_RDD_m[i] = (float)r_RDD_m_double[i];
  }

//Free allocated space
  free(D_RDD_Gy_double);
  free(rdd_parameter_double);
  free(r_RDD_m_double);
}



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
){
  long i;
  const long n_long = (long)(*n);
  const double E_MeV_u_double = (double)(*E_MeV_u);
  const long particle_no_long = (long)(*particle_no);
  const long material_no_long = (long)(*material_no);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);

//Allocate space for the input parameter.
  double* r_m_double = (double*)calloc(n_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	r_m_double[i] = (double)r_m[i];
  }

//Allocate space for the results.
  double* D_RDD_Gy_double = (double*)calloc(n_long,sizeof(double));

//Allocate space for the input parameter.
  double* rdd_parameter_double = (double*)calloc(4,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < 4; i++){
	rdd_parameter_double[i] = (double)rdd_parameter[i];
  }

  int returnValue_internal = 	AT_D_RDD_Gy( n_long,
	r_m_double,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	rdd_model_long,
	rdd_parameter_double,
	er_model_long,
	stopping_power_source_no_long,
	D_RDD_Gy_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	D_RDD_Gy[i] = (float)D_RDD_Gy_double[i];
  }

//Free allocated space
  free(r_m_double);
  free(rdd_parameter_double);
  free(D_RDD_Gy_double);
}



void AT_SPC_read_data_from_filename_fast_R( const char**	filename,
		int*			n,
		int*			depth_step,
		float*			depth_g_cm2,
		float*			E_MeV_u,
		float*			DE_MeV_u,
		int*			particle_no,
		float*			fluence_cm2,
		int*			returnValue
){
  long i;
  int n_long = (int)(*n);
  const char* filename_char = (char*)(*filename);

//Allocate space for the results.
  int* depth_step_long = (int*)calloc(n_long,sizeof(int));
  double* depth_g_cm2_double = (double*)calloc(n_long,sizeof(double));
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  double* DE_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));
  double* fluence_cm2_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_SPC_read_data_from_filename_fast( filename_char,
	n_long,
	depth_step_long,
	depth_g_cm2_double,
	E_MeV_u_double,
	DE_MeV_u_double,
	particle_no_long,
	fluence_cm2_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	depth_step[i] = (int)depth_step_long[i];
	depth_g_cm2[i] = (float)depth_g_cm2_double[i];
	E_MeV_u[i] = (float)E_MeV_u_double[i];
	DE_MeV_u[i] = (float)DE_MeV_u_double[i];
	particle_no[i] = (int)particle_no_long[i];
	fluence_cm2[i] = (float)fluence_cm2_double[i];
  }

//Free allocated space
  free(depth_step_long);
  free(depth_g_cm2_double);
  free(E_MeV_u_double);
  free(DE_MeV_u_double);
  free(particle_no_long);
  free(fluence_cm2_double);
}



void AT_SPC_read_header_from_filename_fast_R( const char**	filename,
		float*			E_MeV_u,
		float*			peak_position_g_cm2,
		int*			particle_no,
		int*			material_no,
		float*			normalisation,
		int*			depth_steps_no,
		int*			returnValue
){
  const char* filename_char = (char*)(*filename);

//Define type-casted output variables
	double E_MeV_u_double = 0;
	double peak_position_g_cm2_double = 0;
	long particle_no_long = 0;
	int material_no_long = 0;
	double normalisation_double = 0;
	int depth_steps_no_long = 0;

  int returnValue_internal = 	AT_SPC_read_header_from_filename_fast( filename_char,
	&E_MeV_u_double,
	&peak_position_g_cm2_double,
	&particle_no_long,
	&material_no_long,
	&normalisation_double,
	&depth_steps_no_long);

//Results:

	*returnValue = ( int )returnValue_internal;

  *E_MeV_u = (float)E_MeV_u_double;

  *peak_position_g_cm2 = (float)peak_position_g_cm2_double;

  *particle_no = (int)particle_no_long;

  *material_no = (int)material_no_long;

  *normalisation = (float)normalisation_double;

  *depth_steps_no = (int)depth_steps_no_long;


//Free allocated space
}



void AT_SPC_get_number_of_bins_from_filename_fast_R( const char**	filename,
		int*			returnValue
){
  const char* filename_char = (char*)(*filename);

  long returnValue_internal = 	AT_SPC_get_number_of_bins_from_filename_fast( filename_char);

//Results:

	*returnValue = ( int )returnValue_internal;


//Free allocated space
}



void AT_Mass_Stopping_Power_R( const char**	stopping_power_source,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_MeV_cm2_g,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const char* stopping_power_source_char = (char*)(*stopping_power_source);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* stopping_power_MeV_cm2_g_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_Mass_Stopping_Power( stopping_power_source_char,
	n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	stopping_power_MeV_cm2_g_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	stopping_power_MeV_cm2_g[i] = (float)stopping_power_MeV_cm2_g_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(stopping_power_MeV_cm2_g_double);
}



void AT_Stopping_Power_R( const char**	stopping_power_source,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_keV_um,
		int*			returnValue
){
  long i;
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const char* stopping_power_source_char = (char*)(*stopping_power_source);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* stopping_power_keV_um_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_Stopping_Power( stopping_power_source_char,
	n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	stopping_power_keV_um_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	stopping_power_keV_um[i] = (float)stopping_power_keV_um_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(stopping_power_keV_um_double);
}



void AT_Mass_Stopping_Power_with_no_R( const int*		stopping_power_source_no,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_MeV_cm2_g,
		int*			returnValue
){
  long i;
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* stopping_power_MeV_cm2_g_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_Mass_Stopping_Power_with_no( stopping_power_source_no_long,
	n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	stopping_power_MeV_cm2_g_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	stopping_power_MeV_cm2_g[i] = (float)stopping_power_MeV_cm2_g_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(stopping_power_MeV_cm2_g_double);
}



void AT_Stopping_Power_with_no_R( const int*		stopping_power_source_no,
		const int*		n,
		const float*	E_MeV_u,
		const int*		particle_no,
		const int*		material_no,
		float*			stopping_power_keV_um,
		int*			returnValue
){
  long i;
  const long stopping_power_source_no_long = (long)(*stopping_power_source_no);
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);

//Allocate space for the input parameter.
  double* E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));


//Fill in the input parameter.
  for(i = 0 ; i < n_long; i++){
	E_MeV_u_double[i] = (double)E_MeV_u[i];
	particle_no_long[i] = (long)particle_no[i];
  }

//Allocate space for the results.
  double* stopping_power_keV_um_double = (double*)calloc(n_long,sizeof(double));

  int returnValue_internal = 	AT_Stopping_Power_with_no( stopping_power_source_no_long,
	n_long,
	E_MeV_u_double,
	particle_no_long,
	material_no_long,
	stopping_power_keV_um_double);

//Results:

	*returnValue = ( int )returnValue_internal;

  for(i = 0 ; i < n_long; i++){
	stopping_power_keV_um[i] = (float)stopping_power_keV_um_double[i];
  }

//Free allocated space
  free(E_MeV_u_double);
  free(particle_no_long);
  free(stopping_power_keV_um_double);
}



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
){
  long i;
  const long n_bins_f_long = (long)(*n_bins_f);
  const double DSB_per_Gy_per_domain_double = (double)(*DSB_per_Gy_per_domain);
  const long domains_per_nucleus_long = (long)(*domains_per_nucleus);
  const bool write_output_bool = (bool)(*write_output);

//Allocate space for the input parameter.
  double* f_d_Gy_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double* f_dd_Gy_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double* f_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double* enhancement_factor_double = (double*)calloc(n_bins_f_long,sizeof(double));


//Fill in the input parameter.
  for(i = 0 ; i < n_bins_f_long; i++){
	f_d_Gy_double[i] = (double)f_d_Gy[i];
	f_dd_Gy_double[i] = (double)f_dd_Gy[i];
	f_double[i] = (double)f[i];
	enhancement_factor_double[i] = (double)enhancement_factor[i];
  }

//Define type-casted output variables
	double total_pDSBs_double = 0;
	double total_nDSBs_double = 0;
	double number_of_iDSBs_double = 0;
	double number_of_cDSBs_double = 0;
	double avg_number_of_DSBs_in_cDSBs_double = 0;

  AT_translate_dose_into_DSB_distribution( n_bins_f_long,
	f_d_Gy_double,
	f_dd_Gy_double,
	f_double,
	enhancement_factor_double,
	DSB_per_Gy_per_domain_double,
	domains_per_nucleus_long,
	write_output_bool,
	&total_pDSBs_double,
	&total_nDSBs_double,
	&number_of_iDSBs_double,
	&number_of_cDSBs_double,
	&avg_number_of_DSBs_in_cDSBs_double);

//Results:
  *total_pDSBs = (float)total_pDSBs_double;

  *total_nDSBs = (float)total_nDSBs_double;

  *number_of_iDSBs = (float)number_of_iDSBs_double;

  *number_of_cDSBs = (float)number_of_cDSBs_double;

  *avg_number_of_DSBs_in_cDSBs = (float)avg_number_of_DSBs_in_cDSBs_double;


//Free allocated space
  free(f_d_Gy_double);
  free(f_dd_Gy_double);
  free(f_double);
  free(enhancement_factor_double);
}



