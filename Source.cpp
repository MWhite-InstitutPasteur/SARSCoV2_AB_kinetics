#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include "randlib.h"
#include <omp.h>
#include <vector>
#include <algorithm>

using namespace std;


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//          //                                        //
//   ####   //  Setting up objects,                   //
//  ##  ##  //  declaring array sizes                 //
//  ##  ##  //                                        //
//  ##  ##  //                                        //
//   ####   //                                        //
//          //                                        //
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
// 0.1 Declare global variables of sizes of arrays

#define N_mcmc 30000000         // Number of MCMC iterations: indexed by mc

#define N_adapt 600000         // Number of MCMC iterations in adaptive phase
#define N_tune_start 1000      // Number of MCMC iterations in adaptive phase
#define N_tune_end 500000

#define N_loc_par 7            // Number of individual-level parameters

#define N_cohort_max 10        // Maximum number of cohorts

#define LARGE 1e12                      // large number needed for priors
#define log2 0.69314718055995
#define sqrt2 1.414214


/////////////////////////////////////////////////////////////////	
// 0.2 Create structure to hold data for participant n
//     and local parameter estimates

struct part_n
{
	//////////////////////////////////////
	// Covariate information

	int status;               // 0 = negative; 1 = positive

	int cohort;               // 0 = Pasteur; 1 = Iyer; 2 = Wang; 3 = Tang


	//////////////////////////////////////
	// Antibody data

	int N_sam;                // Number of samples of antibody data

	vector<double> tt;        // vector of times
	vector<double> AB;        // Vector of antibody levels  

	vector<double> lAB;       // vector of antibody levels (log scale)

	vector<double> y_AB;      // difference between model and data (log scale)


	//////////////////////////////////////
	// Assay range factors

	double AB_min;        // Minimum antibody level for this individual
	double lAB_min;

	double AB_max;        // Maximum antibody level for this individual
	double lAB_max;

	double AB_assay_min;  // Minimum antibody level for the assay
	double lAB_assay_min;

	double AB_assay_max;  // Maximum antibody level for the assay
	double lAB_assay_max;


	//////////////////////////////////////
	// Individual-level parameters

	double AB_bg;         // background antibody level
	double beta;          // boost in ASCs - vaccine dose 1
	double tau;           // delay in boosting of antibody responses
	double t_short;       // half-life of short-lived ASCs
	double t_long;        // half-life of long-lived ASCs
	double t_IgG;         // half-life of IgG molecules
	double rho;           // proportion of short-lived ASCs

	double lAB_bg;        // log(boost in antibody levels)
	double lbeta;         // log(boost in ASCs) 
	double ltau;          // delay in boosting of antibody responses
	double lt_short;      // half-life of short-lived ASCs
	double lt_long;       // half-life of long-lived ASCs
	double lt_IgG;        // half-life of IgG molecules
	double logitrho;      // logit proportion of short-lived ASCs - vaccine dose 1

	double r_short;       // drug decay rate
	double r_long;        // drug decay rate
	double r_IgG;         // drug decay rate


	//////////////////////////////////////
	// Likelihood

	double data_like;     // data likelihood
	double mix_like;      // mixed-effects likelihood

	double lhood;         // individual-level likelihood
};



/////////////////////////////////////////////////////////////////
// 0.3 Create structure for global parameters to be estimated

struct params
{
	/////////////////////////////////////////////////////////////
	// Number of cohorts

	int N_cohort;

	int N_glob_par;


	/////////////////////////////////////////////////////////////
	// Population-level parameters describing mixed effects

	vector<double> mu_par;
	vector<double> tau_par;


	/////////////////////////////////////////////////////////////
	// Count of numbers of samples per cohort

	vector<int> N_sample;


	/////////////////////////////////////////////////////////////
	// Parameter for observational error
	// Separate parameter for each cohort

	vector<double> sig_obs;        // precision of observational error

	vector<double> tau_obs;        // precision of observational error

	vector<double> log_sig_obs;


	/////////////////////////////////////////////////////////////
	// Log likelihood and prior

	double loglike;
	double prior;


	/////////////////////////////////////////////////////////////
	// Prior distributions - global parameters

	vector<double> glob_prior_mean;
	vector<double> glob_prior_mean_CV;
	vector<double> glob_prior_SD;
	vector<double> glob_prior_SD_CV;

	vector<double> glob_prior_logN_mu;
	vector<double> glob_prior_logN_mu_CV;
	vector<double> glob_prior_logN_sig;
	vector<double> glob_prior_logN_sig_CV;

	vector<double> glob_prior_mu;
	vector<double> glob_prior_tau;
	vector<double> glob_prior_k;
	vector<double> glob_prior_theta;


	/////////////////////////////////////////////////////////////
	// Individual-level parameter book-keeping - global parameters

	vector<double> Y_par;
	vector<double> Ymu2_par;


	/////////////////////////////////////////////////////////////
	// Prior distributions - observational parameters

	vector<double> obs_prior_SIG;
	vector<double> obs_prior_SIG_CV;

	vector<double> obs_prior_k;
	vector<double> obs_prior_theta;


	/////////////////////////////////////////////////////////////
	// Individual-level parameter book-keeping - global parameters

	vector<double> Y_obs2;
};


/////////////////////////////////////////////////////////////////
// 0.4 Individual-level structure for MCMC tuning

struct part_n_MCMC
{
	float par_vec[N_loc_par];                             // Parameter vector (in float format for setgmn) (lAB_0, rr)

	float par_vec_test[N_loc_par];                        // Test parameter vector for MCMC update (in float format for setgmn)
	float work[N_loc_par];                                // Scratch vector for setgmn

	double par_S1[N_loc_par];                             // Sum of parameters
	double par_S2[N_loc_par][N_loc_par];                  // Sum of product of pairs

	float COV_MAT[N_loc_par][N_loc_par];                  // covariance matrix (in float format for setgmn)
	float COV_MAT_dummy[N_loc_par][N_loc_par];            // dummy covariance matrix: setgmn gives back sqrt(COV_MAT) or similar so we feed it a dummy

	float GMN_parm[(N_loc_par)*(N_loc_par + 3) / 2 + 1];  // array for setgmn output

	int denom;                                            // denominator for tracking SD calculations

	double step_scale;                                    // scalar for tuning acceptance rate

	int accept;                                           // number of accepted steps
};


////////////////////////////////////////////////////
// 0.5 Initialise functions

double data_like_n(part_n* p, params* theta);
double mix_like_n(part_n* p, params* theta);
double global_prior(params* priors);
double local_prior(part_n* p);
double rm_scale(double step_scale, int step, int N_step_adapt, double log_prob);
double gammln(const double xx);


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//        //                                            //
//   ##   //  Initialise main object, read in data and  //
//  ###   //  fill out objects                          //
//   ##   //                                            //
//   ##   //                                            // 
//  ####  //                                            //
//        //                                            //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// 1.1 Initialise main object 

int main(int argc, char** argv)
{

	//////////////////////////////////////////////////////
	// 1.1.1 Read in file names

	char* glob_prior_File     = argv[1];

	int N_cohort = atoi(argv[2]);

	vector<int> N_t(N_cohort);
	vector<int> N_part(N_cohort);
	vector<int> N_neg(N_cohort);
	char* AB_pos_input_File[N_cohort_max];
	char* AB_neg_input_File[N_cohort_max];
	char* local_output_File[N_cohort_max];

	for (int c = 0; c < N_cohort; c++)
	{
		N_t[c]               = atoi(argv[3 + 6*c]);    // Maximum number of time points
		N_part[c]            = atoi(argv[4 + 6*c]);    // Number of positive participants
		N_neg[c]             = atoi(argv[5 + 6*c]);    // Number of negative controls
		AB_pos_input_File[c] =      argv[6 + 6*c];     // File name for positive samples  
		AB_neg_input_File[c] =      argv[7 + 6*c];     // File name for negative samples  
		local_output_File[c] =      argv[8 + 6*c];     // File name for writing local output
	}

	char* global_output_File = argv[8 + 6*(N_cohort-1) + 1];  // Check this number


	// do we have the correct command line?

	if ( argc != 4 + 6*N_cohort )
	{
		std::cout << "We should be handling " << N_cohort << " cohorts." << endl;
		std::cout << "Please check the command line" << endl;

		system("PAUSE");
	}

	cout << endl;
	cout << "Fitting model to data from " << N_cohort << " cohorts." << endl;

	for (int c = 0; c < N_cohort; c++)
	{
		cout << "Cohort " << c << " with " << N_part[c] << " participants and " << N_neg[c] << " neg controls." << endl;
	}
	cout << endl;


	//////////////////////////////////////////////////////
	// 1.1.2 Count up the number of participants

	int N_part_all = 0;    // Total number of positive participants
	int N_neg_all = 0;     // Total number of negative participants

	for (int c = 0; c < N_cohort; c++)
	{
		N_part_all = N_part_all + N_part[c];
		N_neg_all  = N_neg_all  + N_neg[c];
	}


	//////////////////////////////////////////////////////
	// 1.2 Declare seed, buffer for writing to and clock

	setall(time(NULL), 7);

	int cl = clock();


	///////////////////////////////////////////////////////
	// 1.3 Read in global_priors

	int N_glob_par = 2*N_cohort + 5;  // Number of global parameters


	std::ifstream glob_prior_Stream(glob_prior_File);

	if (glob_prior_Stream.fail())
	{
		std::cout << "Failure reading in data." << endl;
	}

	vector<double> glob_prior_read;
	string discard;

	glob_prior_read.resize(4*N_glob_par + 2*N_cohort);

	for (int i = 0; i < (4*N_glob_par + 2*N_cohort); i++)
	{
		glob_prior_Stream >> glob_prior_read[i] >> discard;
	}

	glob_prior_Stream.close();


	///////////////////////////////////////////////////////
	// 1.4 Read in antibody data (infected individuals)

	vector<vector<vector<double>>> AB_pos_data_read;

	AB_pos_data_read.resize(N_cohort);


	for (int c = 0; c < N_cohort; c++)
	{
		std::ifstream AB_pos_Stream(AB_pos_input_File[c]);

		if (AB_pos_Stream.fail())
		{
			std::cout << "Failure reading in data." << endl;
		}

		AB_pos_data_read[c].resize(N_part[c]);
		for (int i = 0; i < N_part[c]; i++)
		{
			AB_pos_data_read[c][i].resize(2 * N_t[c]);
		}

		for (int i = 0; i < N_part[c]; i++)
		{
			for (int j = 0; j < (2 * N_t[c]); j++)
			{
				AB_pos_Stream >> AB_pos_data_read[c][i][j];
			}
		}

		AB_pos_Stream.close();
	}


	///////////////////////////////////////////////////////
	// 1.5 Read in antibody data (negative controls)

	vector<vector<double>> AB_neg_data_read;

	AB_neg_data_read.resize(N_cohort);

	for (int c = 0; c < N_cohort; c++)
	{
		if (N_neg[c] > 0)
		{
			std::ifstream AB_neg_Stream(AB_neg_input_File[c]);

			if (AB_neg_Stream.fail())
			{
				std::cout << "Failure reading in data." << endl;
			}

			AB_neg_data_read[c].resize(N_neg[c]);

			for (int i = 0; i < N_neg[c]; i++)
			{
				AB_neg_Stream >> AB_neg_data_read[c][i];
			}

			AB_neg_Stream.close();
		}
	}


	///////////////////////////////////////////////////////
	// 1.6 Apply assay ranges and minimum antibody levels
	//     (Population level values)

	///////////////////////////////////////////////////////
	// Minimum antibody level

	vector<double> AB_assay_min(N_cohort);   // Minimum of the assay 

	for (int c = 0; c < N_cohort; c++)
	{
		AB_assay_min[c] = 1e10;

		if (N_neg[c] > 0)
		{
			for (int i = 0; i < N_neg[c]; i++)
			{
				if (AB_neg_data_read[c][i] < AB_assay_min[c])
				{
					AB_assay_min[c] = AB_neg_data_read[c][i];
				}
			}
		}
		if (N_neg[c] == 0)
		{
			AB_assay_min[c] = 0.001;
		}
	}


	for (int c = 0; c < N_cohort; c++)
	{
		for (int i = 0; i < N_part[c]; i++)
		{
			for (int j = N_t[c]; j < (2 * N_t[c]); j++)
			{
				if (AB_pos_data_read[c][i][j] > -50.0)
				{
					if (AB_pos_data_read[c][i][j] < AB_assay_min[c])
					{
						AB_assay_min[c] = AB_pos_data_read[c][i][j];
					}
				}
			}
		}
	}


	///////////////////////////////////////////////////////
	// Maximum antibody level

	vector<double> AB_assay_max(N_cohort);   // Maximum of the assay 

	for (int c = 0; c < N_cohort; c++)
	{
		AB_assay_max[c] = 1e-10;

		for (int i = 0; i < N_part[c]; i++)
		{
			for (int j = N_t[c]; j < (2 * N_t[c]); j++)
			{
				if (AB_pos_data_read[c][i][j] > AB_assay_max[c])
				{
					AB_assay_max[c] = AB_pos_data_read[c][i][j];
				}
			}
		}
	}


	//////////////////////////////////////////////////////
	// 1.7 Create parameter objects

	params theta, theta_p1;

	theta.N_cohort = N_cohort;

	theta.N_glob_par = 2 * theta.N_cohort + 5;


	theta.N_sample.resize(N_cohort);

	theta.sig_obs.resize(N_cohort);
	theta.tau_obs.resize(N_cohort);
	theta.log_sig_obs.resize(N_cohort);

	theta.obs_prior_SIG.resize(N_cohort);
	theta.obs_prior_SIG_CV.resize(N_cohort);
	theta.obs_prior_k.resize(N_cohort);
	theta.obs_prior_theta.resize(N_cohort);
	theta.Y_obs2.resize(N_cohort);


	theta.mu_par.resize(N_glob_par);
	theta.tau_par.resize(N_glob_par);

	theta.glob_prior_mean.resize(N_glob_par);
	theta.glob_prior_mean_CV.resize(N_glob_par);
	theta.glob_prior_SD.resize(N_glob_par);
	theta.glob_prior_SD_CV.resize(N_glob_par);

	theta.glob_prior_logN_mu.resize(N_glob_par);
	theta.glob_prior_logN_mu_CV.resize(N_glob_par);
	theta.glob_prior_logN_sig.resize(N_glob_par);
	theta.glob_prior_logN_sig_CV.resize(N_glob_par);

	theta.glob_prior_mu.resize(N_glob_par);
	theta.glob_prior_tau.resize(N_glob_par);
	theta.glob_prior_k.resize(N_glob_par);
	theta.glob_prior_theta.resize(N_glob_par);

	theta.Y_par.resize(N_glob_par);
	theta.Ymu2_par.resize(N_glob_par);


	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// 1.8.1 Global/Populations parameters

	//////////////////////////////////////////////////////
	// 1.8.1.1 Priors on global parameter

	for (int p = 0; p < N_glob_par; p++)
	{
		theta.glob_prior_mean[p]    = glob_prior_read[4*p];           // Mean value of the parameter in the population
		theta.glob_prior_mean_CV[p] = glob_prior_read[4*p + 1];       // Coefficient of variation in estimate of population - level parameter

		theta.glob_prior_SD[p]    = glob_prior_read[4*p + 2];       // Between Person standard deviation
		theta.glob_prior_SD_CV[p] = glob_prior_read[4*p + 3];       // Coefficient of variation in Between Person sd
	}


	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// 1.8.1.2 Transformation of priors 

	//////////////////////////////////////////////////////
	// Calculation of parameters on log scale based on linear 
	// Note that we make the assumption that CV is invariant
	// under the log transformation - not necessarily true 
	// and needs to be investigated further.

	for (int p = 0; p < (N_glob_par-1); p++)
	{
		theta.glob_prior_logN_mu[p]     = log( theta.glob_prior_mean[p] / sqrt(1.0 + pow(theta.glob_prior_SD[p] / theta.glob_prior_mean[p], 2.0)) );
		theta.glob_prior_logN_mu_CV[p]  = theta.glob_prior_mean_CV[p];

		theta.glob_prior_logN_sig[p]    = sqrt(log( 1.0 + pow(theta.glob_prior_SD[p] / theta.glob_prior_mean[p], 2.0) ));
		theta.glob_prior_logN_sig_CV[p] = theta.glob_prior_SD_CV[p];
	}


	//////////////////////////////////////////////////////
	// Calculation of parameters on logit scale based on linear
	// These parameters are based on R script for priors 

	theta.glob_prior_logN_mu[N_glob_par-1]    = 3.306613;
	theta.glob_prior_logN_mu_CV[N_glob_par-1] = 2.5*theta.glob_prior_mean_CV[N_glob_par-1];

	theta.glob_prior_logN_sig[N_glob_par-1]    = 0.9374749;
	theta.glob_prior_logN_sig_CV[N_glob_par-1] = 2.5*theta.glob_prior_SD_CV[N_glob_par-1];


	//////////////////////////////////////////////////////
	// Calculation of parameters for conjugate prior

	for (int p = 0; p < N_glob_par; p++)
	{
		theta.glob_prior_mu[p]    = theta.glob_prior_logN_mu[p];
		theta.glob_prior_tau[p]   = 1.0 / pow(theta.glob_prior_logN_mu[p] * theta.glob_prior_logN_mu_CV[p], 2.0);

		theta.glob_prior_k[p]     = 1.0 / pow(2.0*theta.glob_prior_logN_sig_CV[p], 2.0);
		theta.glob_prior_theta[p] = pow(2.0*theta.glob_prior_logN_sig_CV[p] / theta.glob_prior_logN_sig[p], 2.0);
	}


	//////////////////////////////////////////////////////
	// Random initialization of parameters

	for (int p = 0; p < N_glob_par; p++)
	{
		theta.mu_par[p]  = genunf(0.9, 1.1)*theta.glob_prior_mu[p];
		theta.tau_par[p] = 0.1*genunf(0.9, 1.1)*theta.glob_prior_k[p] * theta.glob_prior_theta[p];
	}


	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// 1.9.2 Observational parameters

	//////////////////////////////////////////////////////
	// 1.9.2.1 Priors on observational parameters

	for (int c = 0; c < N_cohort; c++)
	{
		theta.obs_prior_SIG[c]    = glob_prior_read[4*N_glob_par + 2*c];        // Between Person standard deviation
		theta.obs_prior_SIG_CV[c] = glob_prior_read[4*N_glob_par + 2*c + 1];    // Coefficient of variation in Between Person sd
	}


	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// 1.9.2.2 Transformation of priors 

	//////////////////////////////////////////////////////
	// Calculation of parameters for conjugate prior

	for (int c = 0; c < N_cohort; c++)
	{
		theta.obs_prior_k[c]     = 1.0 / pow(2.0*theta.obs_prior_SIG_CV[c], 2.0);
		theta.obs_prior_theta[c] = pow(2.0*theta.obs_prior_SIG_CV[c] / theta.obs_prior_SIG[c], 2.0);
	}


	//////////////////////////////////////////////////////
	// Random initialization of parameters

	for (int c = 0; c < N_cohort; c++)
	{
		theta.tau_obs[c] = 0.1*genunf(0.9, 1.1)*theta.obs_prior_k[c] * theta.obs_prior_theta[c];

		theta.sig_obs[c] = 1.0 / sqrt(theta.tau_obs[c]);
		theta.log_sig_obs[c] = log(theta.sig_obs[c]);
	}


	//////////////////////////////////////////////////////////
	// 1.10 Create individual-level objects for participant n

	part_n* part;
	part = new part_n[N_part_all + N_neg_all];

	for (int c = 0; c < N_cohort; c++)
	{
		theta.N_sample[c] = 0;
	}


	//////////////////////////////////////////////////////////
	// 1.10.1 Read in positives

	int n_count = 0;
	double beta_target;


	for (int c = 0; c < N_cohort; c++)
	{

		for (int n = 0; n < N_part[c]; n++)
		{
			part[n_count].cohort = c;

			part[n_count].status = 1;


			/////////////////////////////////////////////////
			// Fill antibody data

			part[n_count].N_sam = 0;

			for (int j = 0; j < N_t[c]; j++)
			{
				if (AB_pos_data_read[c][n][j] > -50.0)
				{
					part[n_count].tt.push_back(AB_pos_data_read[c][n][j]);
					part[n_count].AB.push_back(AB_pos_data_read[c][n][N_t[c] + j]);

					part[n_count].lAB.push_back(log(AB_pos_data_read[c][n][N_t[c] + j]));

					part[n_count].y_AB.push_back(100.0);

					part[n_count].N_sam = part[n_count].N_sam + 1;
				}
			}

			theta.N_sample[c] = theta.N_sample[c] + part[n_count].N_sam;


			/////////////////////////////////////////////////
			// Minimum and maximum for this individual

			part[n_count].AB_min = 1e10;
			part[n_count].AB_max = 1e-10;

			for (int j = 0; j < part[n_count].N_sam; j++)
			{
				if (part[n_count].AB[j] < part[n_count].AB_min)
				{
					part[n_count].AB_min = part[n_count].AB[j];
				}
			}


			for (int j = 0; j < part[n_count].N_sam; j++)
			{
				if (part[n_count].AB[j] > part[n_count].AB_max)
				{
					part[n_count].AB_max = part[n_count].AB[j];
				}
			}

			part[n_count].lAB_min = log(part[n_count].AB_min);
			part[n_count].lAB_max = log(part[n_count].AB_max);


			/////////////////////////////////////////////////
			// Set minimum and maximum antibody level of assay

			part[n_count].AB_assay_min = AB_assay_min[c];
			part[n_count].lAB_assay_min = log(part[n_count].AB_assay_min);

			part[n_count].AB_assay_max = AB_assay_max[c];
			part[n_count].lAB_assay_max = log(part[n_count].AB_assay_max);


			/////////////////////////////////////////////////
			// Randomly assign individual-level parameters
			// (except for beta)

			part[n_count].AB_bg   = genunf( 0.9*part[n_count].AB_assay_min, part[n_count].AB_min );
			part[n_count].tau     = genunf(0.1, 5.0);
			part[n_count].t_short = genunf(5.0, 10.0);
			part[n_count].t_long  = genunf(500.0, 1500.0);
			part[n_count].t_IgG   = genunf(10.0, 30.0);
			part[n_count].rho     = genunf(0.8, 0.9);

			part[n_count].r_short = log2 / part[n_count].t_short;
			part[n_count].r_long  = log2 / part[n_count].t_long;
			part[n_count].r_IgG   = log2 / part[n_count].t_IgG;


			/////////////////////////////////////////////////
			// Randomly assign beta within target range

			beta_target = ( (part[n_count].rho / (part[n_count].r_short - part[n_count].r_IgG))*(exp(-part[n_count].r_IgG * (20.0 - part[n_count].tau)) - exp(-part[n_count].r_short * (20.0 - part[n_count].tau))) +
			                ((1.0 - part[n_count].rho) / (part[n_count].r_long - part[n_count].r_IgG))*(exp(-part[n_count].r_IgG * (20.0 - part[n_count].tau)) - exp(-part[n_count].r_long *  (20.0 - part[n_count].tau))));

			part[n_count].beta = genunf( 0.1*part[n_count].AB_max/beta_target, 10.0*part[n_count].AB_max/beta_target );  // KEEP AN EYE ON THIS..............

			part[n_count].lAB_bg   = log(part[n_count].AB_bg);
			part[n_count].lbeta    = log(part[n_count].beta);
			part[n_count].ltau     = log(part[n_count].tau);
			part[n_count].lt_short = log(part[n_count].t_short);
			part[n_count].lt_long  = log(part[n_count].t_long);
			part[n_count].lt_IgG   = log(part[n_count].t_IgG);

			part[n_count].logitrho = log(part[n_count].rho / (1.0 - part[n_count].rho));

			
			/////////////////////////////////////////////////
			// Calculate individual-level likelihood

			part[n_count].data_like = data_like_n(&part[n_count], &theta);

			part[n_count].mix_like = mix_like_n(&part[n_count], &theta);

			part[n_count].lhood = part[n_count].data_like + part[n_count].mix_like;

			n_count = n_count + 1;

		}

	}

	AB_pos_data_read.clear();


	//////////////////////////////////////////////////////////
	// 1.10.2 Read in negatives

	for (int c = 0; c < N_cohort; c++)
	{

		for (int n = 0; n < N_neg[c]; n++)
		{
			part[n_count].cohort = c;

			part[n_count].status = 0;


			/////////////////////////////////////////////////
			// Fill antibody data

			part[n_count].N_sam = 1;

			part[n_count].tt.push_back(0.0);
			part[n_count].AB.push_back(AB_neg_data_read[c][n]);

			part[n_count].lAB.push_back(log(AB_neg_data_read[c][n]));

			part[n_count].y_AB.push_back(100.0);

			theta.N_sample[c] = theta.N_sample[c] + part[n_count].N_sam;


			/////////////////////////////////////////////////
			// Background antibody level set equal to single measurement

			part[n_count].AB_bg  = part[n_count].AB[0];
			part[n_count].lAB_bg = log(part[n_count].AB_bg);


			/////////////////////////////////////////////////
			// Calculate individual-level likelihood

			part[n_count].data_like = data_like_n(&part[n_count], &theta);

			part[n_count].mix_like = mix_like_n(&part[n_count], &theta);

			part[n_count].lhood = part[n_count].data_like + part[n_count].mix_like;

			n_count = n_count + 1;

		}
	}

	AB_neg_data_read.clear();


	//////////////////////////////////////////////////////
	// 1.11 Initialise adaptive MCMC object for individual-level parameters
	//      One object for each participant.

	part_n_MCMC* part_MCMC;
	part_MCMC = new part_n_MCMC[N_part_all + N_neg_all];

	for (int n = 0; n<(N_part_all + N_neg_all); n++)
	{
		///////////////////////////////////
		// Parameter vector for MVN update

		part_MCMC[n].par_vec[0] = part[n].lAB_bg;
		part_MCMC[n].par_vec[1] = part[n].lbeta;
		part_MCMC[n].par_vec[2] = part[n].ltau;
		part_MCMC[n].par_vec[3] = part[n].lt_short;
		part_MCMC[n].par_vec[4] = part[n].lt_long;
		part_MCMC[n].par_vec[5] = part[n].lt_IgG;
		part_MCMC[n].par_vec[6] = part[n].logitrho;


		/////////////////////////////
		// Initialise diagonal covariance matrix

		for (int p = 0; p<N_loc_par; p++)
		{
			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].COV_MAT[p][q] = 0.0;
			}
		}

		part_MCMC[n].COV_MAT[0][0] = 0.2*0.2;
		part_MCMC[n].COV_MAT[1][1] = 0.2*0.2;
		part_MCMC[n].COV_MAT[2][2] = 0.2*0.2;
		part_MCMC[n].COV_MAT[3][3] = 0.2*0.2;
		part_MCMC[n].COV_MAT[4][4] = 0.2*0.2;
		part_MCMC[n].COV_MAT[5][5] = 0.5*0.5;
		part_MCMC[n].COV_MAT[6][6] = 0.2*0.2;


		/////////////////////////////
		// Counting moments

		for (int p = 0; p<N_loc_par; p++)
		{
			part_MCMC[n].par_S1[p] = part_MCMC[n].par_vec[p];

			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].par_S2[p][q] = part_MCMC[n].par_vec[p] * part_MCMC[n].par_vec[q];
			}
		}

		part_MCMC[n].denom = 1;


		/////////////////////////////
		// Set up dummy covariance matrix including
		// step-size scaling

		part_MCMC[n].step_scale = 0.1;

		for (int p = 0; p<N_loc_par; p++)
		{
			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].COV_MAT_dummy[p][q] = part_MCMC[n].step_scale*part_MCMC[n].COV_MAT[p][q];
			}
		}

		part_MCMC[n].accept = 0.0;
	}


	////////////////////////////////////////////////////////
	// 1.12 Book-keeping

	for (int p = 0; p < N_glob_par; p++)
	{
		theta.Y_par[p]    = 0.0;
		theta.Ymu2_par[p] = 0.0;

		for (int n = 0; n < (N_part_all + N_neg_all); n++)
		{
			theta.Y_par[p]    = theta.Y_par[p] + part_MCMC[n].par_vec[p];
			theta.Ymu2_par[p] = theta.Ymu2_par[p] + (part_MCMC[n].par_vec[p] - theta.mu_par[p])*(part_MCMC[n].par_vec[p] - theta.mu_par[p]);
		}
	}

	theta_p1 = theta;


	////////////////////////////////////////////////////////
	// 1.13 Create objects for updating local parameters

	part_n* part_p1;
	part_p1 = new part_n[N_part_all + N_neg_all];

	for (int n = 0; n<(N_part_all + N_neg_all); n++)
	{
		part_p1[n] = part[n];
	}


	///////////////////////////////////////////////////////////////////////////
	// 1.11 Test output of likelihood (and print to screen)

	cout << "Some test output from individuals." << endl;

	//for (int n = 0; n < N_part_all; n++)
	for (int n = 0; n < 10; n++)
	{
		std::cout << n << "\t" << "AB_bg: " << part[n].lAB_bg << "\t" << "beta: " << part[n].lbeta << "\t" << "tau: " << part[n].ltau << "\t" 
			                   << "t_short: " << part[n].lt_short << "\t" << "t_long: " << part[n].lt_long << "\t" << "t_IgG: " << part[n].lt_IgG << "\t"
			                   << "rho: " << part[n].logitrho << "\t" << "logL: " << part[n].lhood << endl;
	}



	//for (int n = N_part_all; n < (N_part_all + N_neg_all); n++)
	for (int n = N_part_all; n < (N_part_all + 10); n++)
	{
		std::cout << n << "\t" << "AB_bg: " << part[n].lAB_bg << "\t" <<  "logL: " << part[n].lhood << endl;
	}

	cout << endl;


	///////////////////////////////////////////////////////////////////////////
	// 1.12 Initialise parameters for MCMC likelihood, Robbins-Munro 
	//      acceptance and output

	double loglike = global_prior(&theta);
	for (int n = 0; n<(N_part_all + N_neg_all); n++)
	{
		loglike = loglike + part[n].lhood;
	}


	double log_prob, loglike_p1;

	vector<double> log_loc_prob(N_part_all + N_neg_all);


	int screen_out = max(2, (int)((int)N_mcmc) / 1000);
	int glob_out = max(2, (int)((int)N_mcmc) / 10000);
	int loc_out = max(2, (int)((int)N_mcmc) / 1000);


	vector<double> randomU(N_part_all + N_neg_all);

	vector<double> loglike_vec_p1(N_part_all + N_neg_all);


	///////////////////////////////////////////////////////////////////////////
	// 1.13 Open file for output and write first line

	///////////////////////////////////////////////////////////////////////////
	// 1.13.1 Global output to screen

	std::cout << "START MCMC" << endl;

	std::cout << 0 << "\t";
	for (int p = 0; p < N_glob_par; p++)
	{
		std::cout << theta.mu_par[p] << "\t";
	}
	for (int p = 0; p < N_glob_par; p++)
	{
		std::cout << theta.tau_par[p] << "\t";
	}
	for (int c = 0; c < N_cohort; c++)
	{
		std::cout << theta.sig_obs[c] << "\t";
	}
	std::cout << loglike << "\t" << global_prior(&theta) << endl;


	///////////////////////////////////////////////////////////////////////////
	// 1.13.2 Global output to file

	std::ofstream global_MCMC_Stream(global_output_File);
	
	for (int p = 0; p < N_glob_par; p++)
	{
		global_MCMC_Stream << theta.mu_par[p] << "\t";
	}
	for (int p = 0; p < N_glob_par; p++)
	{
		global_MCMC_Stream << theta.tau_par[p] << "\t";
	}
	for (int c = 0; c < N_cohort; c++)
	{
		global_MCMC_Stream << theta.sig_obs[c] << "\t";
	}
	global_MCMC_Stream << "\t" << loglike << "\t" << global_prior(&theta) << endl;


	///////////////////////////////////////////////////////////////////////////
	// 1.13.3 Local output to file

	for (int c = 0; c < N_cohort; c++)
	{
		ofstream local_MCMC_Stream(local_output_File[c]);

		for (int n = 0; n < N_part_all; n++)
		{
			if (part[n].cohort == c)
			{
				local_MCMC_Stream << part[n].AB_bg << "\t" << part[n].beta << "\t" << part[n].tau << "\t" 
					              << part[n].t_short << "\t" << part[n].t_long << "\t" << part[n].t_IgG << "\t" << part[n].rho << "\t" << part[n].lhood << "\t";
			}
		}
		local_MCMC_Stream << endl;

		local_MCMC_Stream.close();
	}


	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	//          //                                         //
	//   ####   //  Begin MCMC fitting procedure           //
	//  ##  ##  //                                         //
	//     ##   //                                         //
	//    ##    //                                         //
	//   #####  //                                         //
	//          //                                         //
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////

	for (int mc = 1; mc<N_mcmc; mc++)
	{
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////
		//       //                                    //
		//  2.1  //  UPDATE STAGE 1: INDIVIDUAL-LEVEL  //
		//       //  Metropolis-Hastings sampler       //
		//       //  (Positive individuals only)       //
		//       //                                    //
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////

		/////////////////////////////////////////////
		// 2.1.1. Proposal step

		for (int n = 0; n < (N_part_all); n++)
		{
			////////////////////////////////////////////////
			// Update COV_MAT_dummay

			for (int p = 0; p < N_loc_par; p++)
			{
				for (int q = 0; q < N_loc_par; q++)
				{
					part_MCMC[n].COV_MAT_dummy[p][q] = part_MCMC[n].step_scale*part_MCMC[n].COV_MAT[p][q];
				}
			}


			///////////////////////////////////////////////
			// Multi-variate Normal proposal step

			setgmn(part_MCMC[n].par_vec, *part_MCMC[n].COV_MAT_dummy, N_loc_par, part_MCMC[n].GMN_parm);

			genmn(part_MCMC[n].GMN_parm, part_MCMC[n].par_vec_test, part_MCMC[n].work);

			part_p1[n].lAB_bg   = part_MCMC[n].par_vec_test[0];
			part_p1[n].lbeta    = part_MCMC[n].par_vec_test[1];
			part_p1[n].ltau     = part_MCMC[n].par_vec_test[2];
			part_p1[n].lt_short = part_MCMC[n].par_vec_test[3];
			part_p1[n].lt_long  = part_MCMC[n].par_vec_test[4];
			part_p1[n].lt_IgG   = part_MCMC[n].par_vec_test[5];
			part_p1[n].logitrho = part_MCMC[n].par_vec_test[6];

			randomU[n] = genunf(0.0, 1.0);
		}


		/////////////////////////////////////////////
		// 2.1.2. Update step
		
		//#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < (N_part_all); n++)
		{
			////////////////////////////////////////////////////////
			// 2.1.2.1. Only proceed if allowable parameters proposed

			if (local_prior(&part_p1[n]) > -0.5*LARGE)
			{
				part_p1[n].AB_bg = exp(part_p1[n].lAB_bg);

				part_p1[n].beta = exp(part_p1[n].lbeta);

				part_p1[n].tau = exp(part_p1[n].ltau);

				part_p1[n].t_short = exp(part_p1[n].lt_short);
				part_p1[n].r_short = log2 / part_p1[n].t_short;

				part_p1[n].t_long = exp(part_p1[n].lt_long);
				part_p1[n].r_long = log2 / part_p1[n].t_long;

				part_p1[n].t_IgG = exp(part_p1[n].lt_IgG);
				part_p1[n].r_IgG = log2 / part_p1[n].t_IgG;

				part_p1[n].rho = exp(part_p1[n].logitrho) / (1.0 + exp(part_p1[n].logitrho));
				

				part_p1[n].data_like = data_like_n(&part_p1[n], &theta);

				part_p1[n].mix_like = mix_like_n(&part_p1[n], &theta);

				part_p1[n].lhood = part_p1[n].data_like + part_p1[n].mix_like;


				double log_prob_n = part_p1[n].lhood - part[n].lhood;

				log_loc_prob[n] = _finite(log_prob_n) ? std::min(log_prob_n, 0.0) : -LARGE;


				////////////////////////////////////////
				// 2.1.2.2. Update if necessary

				if (log(randomU[n]) < log_loc_prob[n])
				{
					part[n] = part_p1[n];

					for (int p = 0; p < N_loc_par; p++)
					{
						part_MCMC[n].par_vec[p] = part_MCMC[n].par_vec_test[p];
					}

					part_MCMC[n].accept = part_MCMC[n].accept + 1;
				}


				////////////////////////////////////////
				// 2.1.3. Adjust step-size with Robbins-Monro
				//        Only do this for a local step within allowed range

				if (mc < N_adapt)
				{
					part_MCMC[n].step_scale = rm_scale(part_MCMC[n].step_scale, mc, N_adapt, log_loc_prob[n]);
				}
			}


			////////////////////////////////////////////////////////////
			// Running account of sums and sums of squares

			if (mc < N_tune_end)
			{
				for (int p = 0; p < N_loc_par; p++)
				{
					part_MCMC[n].par_S1[p] = part_MCMC[n].par_S1[p] + part_MCMC[n].par_vec[p];
					
					for (int q = 0; q < N_loc_par; q++)
					{
						part_MCMC[n].par_S2[p][q] = part_MCMC[n].par_S2[p][q] + part_MCMC[n].par_vec[p] * part_MCMC[n].par_vec[q];
					}
				}

				part_MCMC[n].denom = part_MCMC[n].denom + 1;
			}


			////////////////////////////////////////////////////////////
			// 2.1.4. TUNING STAGE 1

			/////////////////////////////////
			// Update covariance matrix

			if ((mc >= N_tune_start) && (mc < N_tune_end))
			{
				for (int p = 0; p < N_loc_par; p++)
				{
					for (int q = 0; q < N_loc_par; q++)
					{
						if (part_MCMC[n].accept / part_MCMC[n].denom > 0.01)
						{
							part_MCMC[n].COV_MAT[p][q] = part_MCMC[n].par_S2[p][q] / (part_MCMC[n].denom) - part_MCMC[n].par_S1[p] * part_MCMC[n].par_S1[q] / (part_MCMC[n].denom*part_MCMC[n].denom);
						}
					}
				}
			}
		}

		//////////////////////////////////////////////////////////////
		// 2.1.5. Update the total likelihood given the local updates

		loglike = global_prior(&theta);
		
		for (int n = 0; n < (N_part_all + N_neg_all); n++)
		{
			loglike = loglike + part[n].lhood;
		}


		///////////////////////////////////////////////////
		///////////////////////////////////////////////////
		//       //                                      //
		//  2.2  //  UPDATE STAGE 2 - POPULATION-LEVEL   //
		//       //  Gibbs sampler                       // 
		//       //                                      //
		///////////////////////////////////////////////////
		///////////////////////////////////////////////////

		///////////////////////////////////////////////////
		// AB_bg: all cohorts 

		for( int c = 0; c < N_cohort; c++ )
		{
			theta.Y_par[c] = 0.0;

			for (int n = 0; n < (N_part_all + N_neg_all); n++)
			{
				if (part[n].cohort == c)
				{
					theta.Y_par[c] = theta.Y_par[c] + part_MCMC[n].par_vec[0];
				}
			}

			theta.mu_par[c] = gennor( (theta.glob_prior_tau[c] * theta.glob_prior_mu[c] + theta.tau_par[c] * theta.Y_par[c]) / (theta.glob_prior_tau[c] + (N_part[c] + N_neg[c])* theta.tau_par[c]),
				                      1.0 / sqrt(theta.glob_prior_tau[c] + (N_part[c] + N_neg[c])* theta.tau_par[c]));

			theta.Ymu2_par[c] = 0.0;

			for (int n = 0; n < (N_part_all + N_neg_all); n++)
			{
				if (part[n].cohort == c)
				{
					theta.Ymu2_par[c] = theta.Ymu2_par[c] + (part_MCMC[n].par_vec[0] - theta.mu_par[c])*(part_MCMC[n].par_vec[0] - theta.mu_par[c]);
				}
			}

			theta.tau_par[c] = gengam( 1.0 / theta.glob_prior_theta[c] + 0.5*theta.Ymu2_par[c],
				                       0.5*(N_part[c] + N_neg[c]) + theta.glob_prior_k[c]);
		}


		///////////////////////////////////////////////////
		// beta: all cohorts 

		for (int c = 0; c < N_cohort; c++)
		{
			theta.Y_par[N_cohort + c] = 0.0;

			for (int n = 0; n < N_part_all; n++)
			{
				if (part[n].cohort == c)
				{
					theta.Y_par[N_cohort + c] = theta.Y_par[N_cohort + c] + part_MCMC[n].par_vec[1];
				}
			}

			theta.mu_par[N_cohort + c] = gennor( (theta.glob_prior_tau[N_cohort + c] * theta.glob_prior_mu[N_cohort + c] + theta.tau_par[N_cohort + c] * theta.Y_par[N_cohort + c]) / (theta.glob_prior_tau[N_cohort + c] + N_part[c] * theta.tau_par[N_cohort + c]),
				                                  1.0 / sqrt(theta.glob_prior_tau[N_cohort + c] + N_part[c] * theta.tau_par[N_cohort + c]));

			theta.Ymu2_par[N_cohort + c] = 0.0;

			for (int n = 0; n < N_part_all; n++)
			{
				if (part[n].cohort == c)
				{
					theta.Ymu2_par[N_cohort + c] = theta.Ymu2_par[N_cohort + c] + (part_MCMC[n].par_vec[1] - theta.mu_par[N_cohort + c])*(part_MCMC[n].par_vec[1] - theta.mu_par[N_cohort + c]);
				}
			}

			theta.tau_par[N_cohort + c] = gengam( 1.0 / theta.glob_prior_theta[N_cohort + c] + 0.5*theta.Ymu2_par[N_cohort + c],
				                                  0.5*N_part[c] + theta.glob_prior_k[N_cohort + c]);
		}


		///////////////////////////////////////////////////
		// tau 

		theta.Y_par[2*N_cohort] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Y_par[2 * N_cohort] = theta.Y_par[2 * N_cohort] + part_MCMC[n].par_vec[2];
		}

		theta.mu_par[2 * N_cohort] = gennor( (theta.glob_prior_tau[2 * N_cohort] * theta.glob_prior_mu[2 * N_cohort] + theta.tau_par[2 * N_cohort] * theta.Y_par[2 * N_cohort]) / (theta.glob_prior_tau[2 * N_cohort] + N_part_all * theta.tau_par[2 * N_cohort]),
			                                 1.0 / sqrt(theta.glob_prior_tau[2 * N_cohort] + N_part_all * theta.tau_par[2 * N_cohort]) );

		theta.Ymu2_par[2 * N_cohort] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Ymu2_par[2 * N_cohort] = theta.Ymu2_par[2 * N_cohort] + (part_MCMC[n].par_vec[2] - theta.mu_par[2 * N_cohort])*(part_MCMC[n].par_vec[2] - theta.mu_par[2 * N_cohort]);
		}

		theta.tau_par[2 * N_cohort] = gengam( 1.0 / theta.glob_prior_theta[2 * N_cohort] + 0.5*theta.Ymu2_par[2 * N_cohort],
			                                  0.5*N_part_all + theta.glob_prior_k[2 * N_cohort]);


		///////////////////////////////////////////////////
		// t_short

		theta.Y_par[2 * N_cohort + 1] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Y_par[2 * N_cohort + 1] = theta.Y_par[2 * N_cohort + 1] + part_MCMC[n].par_vec[3];
		}

		theta.mu_par[2 * N_cohort + 1] = gennor( (theta.glob_prior_tau[2 * N_cohort + 1] * theta.glob_prior_mu[2 * N_cohort + 1] + theta.tau_par[2 * N_cohort + 1] * theta.Y_par[2 * N_cohort + 1]) / (theta.glob_prior_tau[2 * N_cohort + 1] + N_part_all *theta.tau_par[2 * N_cohort + 1]),
				                                 1.0 / sqrt(theta.glob_prior_tau[2 * N_cohort + 1] + N_part_all *theta.tau_par[2 * N_cohort + 1]) );

		theta.Ymu2_par[2 * N_cohort + 1] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Ymu2_par[2 * N_cohort + 1] = theta.Ymu2_par[2 * N_cohort + 1] + (part_MCMC[n].par_vec[3] - theta.mu_par[2 * N_cohort + 1])*(part_MCMC[n].par_vec[3] - theta.mu_par[2 * N_cohort + 1]);
		}

		theta.tau_par[2 * N_cohort + 1] = gengam( 1.0/theta.glob_prior_theta[2 * N_cohort + 1] + 0.5*theta.Ymu2_par[2 * N_cohort + 1],
			                                      0.5*N_part_all + theta.glob_prior_k[2 * N_cohort + 1]);
		

		///////////////////////////////////////////////////
		// t_long

		theta.Y_par[2 * N_cohort + 2] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Y_par[2 * N_cohort + 2] = theta.Y_par[2 * N_cohort + 2] + part_MCMC[n].par_vec[4];
		}

		theta.mu_par[2 * N_cohort + 2] = gennor( (theta.glob_prior_tau[2 * N_cohort + 2] * theta.glob_prior_mu[2 * N_cohort + 2] + theta.tau_par[2 * N_cohort + 2] * theta.Y_par[2 * N_cohort + 2]) / (theta.glob_prior_tau[2 * N_cohort + 2] + N_part_all *theta.tau_par[2 * N_cohort + 2]),
			                                      1.0 / sqrt(theta.glob_prior_tau[2 * N_cohort + 2] + N_part_all *theta.tau_par[2 * N_cohort + 2]) );

		theta.Ymu2_par[2 * N_cohort + 2] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Ymu2_par[2 * N_cohort + 2] = theta.Ymu2_par[2 * N_cohort + 2] + (part_MCMC[n].par_vec[4] - theta.mu_par[2 * N_cohort + 2])*(part_MCMC[n].par_vec[4] - theta.mu_par[2 * N_cohort + 2]);
		}

		theta.tau_par[2 * N_cohort + 2] = gengam( 1.0 / theta.glob_prior_theta[2 * N_cohort + 2] + 0.5*theta.Ymu2_par[2 * N_cohort + 2],
			                                      0.5*N_part_all + theta.glob_prior_k[2 * N_cohort + 2] );


		///////////////////////////////////////////////////
		// t_IgG

		theta.Y_par[2 * N_cohort + 3] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Y_par[2 * N_cohort + 3] = theta.Y_par[2 * N_cohort + 3] + part_MCMC[n].par_vec[5];
		}

		theta.mu_par[2 * N_cohort + 3] = gennor( (theta.glob_prior_tau[2 * N_cohort + 3] * theta.glob_prior_mu[2 * N_cohort + 3] + theta.tau_par[2 * N_cohort + 3] * theta.Y_par[2 * N_cohort + 3]) / (theta.glob_prior_tau[2 * N_cohort + 3] + N_part_all *theta.tau_par[2 * N_cohort + 3]),
			                                     1.0 / sqrt(theta.glob_prior_tau[2 * N_cohort + 3] + N_part_all *theta.tau_par[2 * N_cohort + 3]));

		theta.Ymu2_par[2 * N_cohort + 3] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Ymu2_par[2 * N_cohort + 3] = theta.Ymu2_par[2 * N_cohort + 3] + (part_MCMC[n].par_vec[5] - theta.mu_par[2 * N_cohort + 3])*(part_MCMC[n].par_vec[5] - theta.mu_par[2 * N_cohort + 3]);
		}

		theta.tau_par[2 * N_cohort + 3] = gengam( 1.0 / theta.glob_prior_theta[2 * N_cohort + 3] + 0.5*theta.Ymu2_par[2 * N_cohort + 3],
			                                      0.5*N_part_all + theta.glob_prior_k[2 * N_cohort + 3]);


		///////////////////////////////////////////////////
		// rho

		theta.Y_par[2 * N_cohort + 4] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Y_par[2 * N_cohort + 4] = theta.Y_par[2 * N_cohort + 4] + part_MCMC[n].par_vec[6];
		}

		theta.mu_par[2 * N_cohort + 4] = gennor( (theta.glob_prior_tau[2 * N_cohort + 4] * theta.glob_prior_mu[2 * N_cohort + 4] + theta.tau_par[2 * N_cohort + 4] * theta.Y_par[2 * N_cohort + 4]) / (theta.glob_prior_tau[2 * N_cohort + 4] + N_part_all *theta.tau_par[2 * N_cohort + 4]),
			                                     1.0 / sqrt(theta.glob_prior_tau[2 * N_cohort + 4] + N_part_all *theta.tau_par[2 * N_cohort + 4]) );

		theta.Ymu2_par[2 * N_cohort + 4] = 0.0;

		for (int n = 0; n < N_part_all; n++)
		{
			theta.Ymu2_par[2 * N_cohort + 4] = theta.Ymu2_par[2 * N_cohort + 4] + (part_MCMC[n].par_vec[6] - theta.mu_par[2 * N_cohort + 4])*(part_MCMC[n].par_vec[6] - theta.mu_par[2 * N_cohort + 4]);
		}

		theta.tau_par[2 * N_cohort + 4] = gengam( 1.0 / theta.glob_prior_theta[2 * N_cohort + 4] + 0.5*theta.Ymu2_par[2 * N_cohort + 4],
			                                      0.5*N_part_all + theta.glob_prior_k[2 * N_cohort + 4] );

		theta_p1 = theta;


		//////////////////////////////////////////////////////////////
		// 2.2.1. Update the total likelihood given the local updates

		loglike = global_prior(&theta);

		#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < (N_part_all + N_neg_all); n++)
		{
			//part[n].data_like = data_like_n(&part[n], &theta);

			part[n].mix_like = mix_like_n(&part[n], &theta);

			part[n].lhood = part[n].data_like + part[n].mix_like;

			loglike = loglike + part[n].lhood;
		}


		///////////////////////////////////////////////////
		///////////////////////////////////////////////////
		//       //                                      //
		//  2.3  //  UPDATE STAGE 3 - OBS-VARIANCE       //
		//       //  Gibbs sampler                       // 
		//       //                                      //
		///////////////////////////////////////////////////
		///////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////
		// 2.3.1. Gibbs sampler update

		///////////////////////////////////////////////////
		// Observational variance: 

		for (int c = 0; c < N_cohort; c++)
		{
			theta.Y_obs2[c] = 0.0;

			for (int n = 0; n < (N_part_all + N_neg_all); n++)
			{
				if (part[n].cohort == c)
				{
					for (int j = 0; j < part[n].N_sam; j++)
					{
						theta.Y_obs2[c] = theta.Y_obs2[c] + part[n].y_AB[j] * part[n].y_AB[j];
					}
				}
			}

			theta.tau_obs[c] = gengam( 1.0 / theta.obs_prior_theta[c] + 0.5*theta.Y_obs2[c],
				                       0.5*theta.N_sample[c] + theta.obs_prior_k[c]);

			theta.sig_obs[c]     = 1.0 / sqrt(theta.tau_obs[c]);
			theta.log_sig_obs[c] = log(theta.sig_obs[c]);	
		}


		//////////////////////////////////////////////////////////////
		// 2.3.2. Update the total likelihood given the local updates

		loglike = global_prior(&theta);

		#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < (N_part_all + N_neg_all); n++)
		{
			part[n].data_like = data_like_n(&part[n], &theta);

			//part[n].mix_like = mix_like_n(&part[n], &theta);

			part[n].lhood = part[n].data_like + part[n].mix_like;

			loglike = loglike + part[n].lhood;
		}


		//////////////////////////////////////////////////////////////
		// 2.3.3. Update the total likelihood given the local updates

		loglike = global_prior(&theta);

		for (int n = 0; n < (N_part_all + N_neg_all); n++)
		{
			loglike = loglike + part[n].lhood;
		}


		/////////////////////////////////////////
		/////////////////////////////////////////
		//       //                            //
		//  2.4  //  Output results to file    //
		//       //                            //
		/////////////////////////////////////////
		/////////////////////////////////////////

		//////////////////////////////////////////////////////////////
		// 2.4.1 Global output to screen

		if (mc%screen_out == 0)
		{
			std::cout << 100 * ((double)mc) / ((double)N_mcmc) << "% complete." << "\t";
			for (int p = 0; p < N_glob_par; p++)
			{
				std::cout << theta.mu_par[p] << "\t";
			}
			for (int p = 0; p < N_glob_par; p++)
			{
				std::cout << theta.tau_par[p] << "\t";
			}
			for (int c = 0; c < N_cohort; c++)
			{
				std::cout << theta.sig_obs[c] << "\t";
			}
			std::cout << loglike << "\t" << global_prior(&theta) << endl;
		}


		//////////////////////////////////////////////////////////////
		// 2.4.2 Global output to file

		if (mc%glob_out == 0)
		{
			for (int p = 0; p < N_glob_par; p++)
			{
				global_MCMC_Stream << theta.mu_par[p] << "\t";
			}
			for (int p = 0; p < N_glob_par; p++)
			{
				global_MCMC_Stream << theta.tau_par[p] << "\t";
			}
			for (int c = 0; c < N_cohort; c++)
			{
				global_MCMC_Stream << theta.sig_obs[c] << "\t";
			}
			global_MCMC_Stream << loglike << "\t" << global_prior(&theta) << endl;
		}


		//////////////////////////////////////////////////////////////
		// 2.4.3 Local output to file

		if (mc%loc_out == 0)
		{
			for (int c = 0; c < N_cohort; c++)
			{
				ofstream local_MCMC_Stream(local_output_File[c], ios_base::app);

				for (int n = 0; n < N_part_all; n++)
				{
					if (part[n].cohort == c)
					{
						local_MCMC_Stream << part[n].AB_bg << "\t" << part[n].beta << "\t" <<
							                 part[n].tau << "\t" << part[n].t_short << "\t" << part[n].t_long << "\t" << part[n].t_IgG << "\t" <<
							                 part[n].rho << "\t" << part[n].lhood << "\t";
					}
				}
				local_MCMC_Stream << endl;

				local_MCMC_Stream.close();
			}
		}

	}


	//////////////////////////////////////
	// 2.5. Calculate and output acceptance rates

	std::cout << "local acceptance:  " << part_MCMC[0].accept << "\t" << (double(part_MCMC[0].accept)) / (double(N_mcmc)) << endl;

	std::cout << "Time taken: " << cl << endl;

	return 0;
}



///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//          //                                       //
//   ####   //  Functions                            //
//  ##  ##  //                                       //
//     ##   //                                       //
//  ##  ##  //                                       //
//   ####   //                                       //
//          //                                       // 
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// 3.1 Individual-level likelihood

double data_like_n(part_n* p, params* theta)
{

	///////////////////////////////////////
	// 3.1.1 Model-predicted antibody titre

	vector<double> model_AB(p->N_sam), model_lAB(p->N_sam);

	// Negative controls

	if (p->status == 0)
	{
		model_AB[0] = p->AB_bg;
	}


	// Individuals with longitudinal samples

	if ( p->status == 1 )
	{

		for (int j = 0; j < p->N_sam; j++)
		{
			model_AB[j] = p->AB_bg;  // (background does not decay)   **exp(-p->r_long*p->tt[j])*/;

			if (p->tt[j] > p->tau)
			{
				model_AB[j] = model_AB[j] +
						      p->beta*( (p->rho / (p->r_short - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->tau)) - exp(-p->r_short * (p->tt[j] - p->tau))) +
					                    ((1.0 - p->rho) / (p->r_long - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->tau)) - exp(-p->r_long *  (p->tt[j] - p->tau))));
			}
		}
	}



	for (int j = 0; j < p->N_sam; j++)
	{
		model_lAB[j] = log(model_AB[j]);

		p->y_AB[j] = fabs(model_lAB[j] - p->lAB[j]);
	}


	///////////////////////////////////////
	// 3.1.2 Likelihood
	//       log-Normally distributed error 

	///////////////////////////////////////
	// Sum up likelihood across all samples

	double loglike = 0.0;

	for (int j = 0; j < p->N_sam; j++)
	{
		/////////////////////////////////////////////////////////////
		// CASE 1: Data lies between assay minimum and maximum

		loglike = loglike - 0.9189385 - theta->log_sig_obs[p->cohort] - p->lAB[j] - 0.5*(p->lAB[j] - model_lAB[j])*(p->lAB[j] - model_lAB[j]) / (theta->sig_obs[p->cohort] * theta->sig_obs[p->cohort]);
	}

	return loglike;
}


///////////////////////////////////////////////////////
// 3.2 Mixed-effects likelihood

double mix_like_n(part_n* p, params* theta)
{
	double mixlike = 0.0;

	//////////////////////////////////////////
	// Mixed-effects component of likelihood

	if (p->status == 0)
	{
		mixlike = - 0.9189385 + 0.5*log(theta->tau_par[p->cohort])
			                  - 0.5*theta->tau_par[p->cohort] * (p->lAB_bg - theta->mu_par[p->cohort]) * (p->lAB_bg - theta->mu_par[p->cohort]);
	}

	if (p->status == 1)
	{
		mixlike = -6.43257 + 0.5*log( theta->tau_par[p->cohort] * theta->tau_par[theta->N_cohort + p->cohort] *
			                            theta->tau_par[2*theta->N_cohort] * theta->tau_par[2*theta->N_cohort + 1] * theta->tau_par[2*theta->N_cohort + 2] * 
			                            theta->tau_par[2*theta->N_cohort + 3] * theta->tau_par[2*theta->N_cohort + 4] )
			                 - 0.5*theta->tau_par[p->cohort]                   * (p->lAB_bg   - theta->mu_par[p->cohort])                   * (p->lAB_bg - theta->mu_par[p->cohort])
			                 - 0.5*theta->tau_par[theta->N_cohort + p->cohort] * (p->lbeta    - theta->mu_par[theta->N_cohort + p->cohort]) * (p->lbeta - theta->mu_par[theta->N_cohort + p->cohort])
			                 - 0.5*theta->tau_par[2*theta->N_cohort]           * (p->ltau     - theta->mu_par[2*theta->N_cohort])           * (p->ltau     - theta->mu_par[2*theta->N_cohort])
			                 - 0.5*theta->tau_par[2*theta->N_cohort + 1]       * (p->lt_short - theta->mu_par[2*theta->N_cohort + 1])       * (p->lt_short - theta->mu_par[2*theta->N_cohort + 1])
			                 - 0.5*theta->tau_par[2*theta->N_cohort + 2]       * (p->lt_long  - theta->mu_par[2*theta->N_cohort + 2])       * (p->lt_long  - theta->mu_par[2*theta->N_cohort + 2])
			                 - 0.5*theta->tau_par[2*theta->N_cohort + 3]       * (p->lt_IgG   - theta->mu_par[2*theta->N_cohort + 3])       * (p->lt_IgG   - theta->mu_par[2*theta->N_cohort + 3])
			                 - 0.5*theta->tau_par[2*theta->N_cohort + 4]       * (p->logitrho - theta->mu_par[2*theta->N_cohort + 4])       * (p->logitrho - theta->mu_par[2*theta->N_cohort + 4]);
	}
	

	return mixlike;
}



///////////////////////////////////////////////////////
// 3.3 Individual-level prior likelihood: excludes ineligible steps

double local_prior(part_n* p)
{
	if (p->lAB_bg < p->lAB_assay_min - 2.302585) { return -LARGE; }
	if (p->lAB_bg > p->lAB_min) { return -LARGE; }
	

	if( p->status == 1  )
	{
		if (p->ltau > 3.912023) { return -LARGE; }  // log(50) = 3.912023

		if (p->lt_short < -0.6931472) { return -LARGE; }    // log(0.5) = -0.6931472

		if (p->lt_short > p->lt_long) { return -LARGE; }

		if (p->lt_IgG < 0) { return -LARGE; }         // log(1) = 0
		if (p->lt_IgG > 4.094345) { return -LARGE; }  // log(60) = 4.094345

		if (p->rho < 0.0) { return -LARGE; }
		if (p->rho > 1.0) { return -LARGE; }
	}

	return 0.0;
}


///////////////////////////////////////////////////////////////////
// 3.4 Population-level prior likelihood: excludes ineligible steps

double global_prior(params* priors)
{
	double logprior = 0.0;


	///////////////////////////////////////////////////////
	// 3.4.1 Population/Global priors

	for (int p = 0; p < priors->N_glob_par; p++)
	{
		///////////////////////////
		// Normal prior on mean_par[p] 

		logprior = logprior - 0.9189385 + 0.5*log(priors->glob_prior_tau[p]) - 0.5*priors->glob_prior_tau[p] * (priors->mu_par[p] - priors->glob_prior_mu[p])*(priors->mu_par[p] - priors->glob_prior_mu[p]);


		///////////////////////////
		// Gamma prior on tau_par[p]

		logprior = logprior + (priors->glob_prior_k[p] - 1.0)*log(priors->tau_par[p]) - priors->tau_par[p] / priors->glob_prior_theta[p] - priors->glob_prior_k[p] *log(priors->glob_prior_theta[p]) - gammln(priors->glob_prior_k[p]);
	}


	///////////////////////////////////////////////////////
	// 3.4.2 Observational priors

	for (int c = 0; c < priors->N_cohort; c++)
	{
		///////////////////////////
		// Gamma prior on tau_par[p]

		logprior = logprior + (priors->obs_prior_k[c] - 1.0)*log(priors->tau_obs[c]) - priors->tau_obs[c] / priors->obs_prior_theta[c] - priors->obs_prior_k[c] * log(priors->obs_prior_theta[c]) - gammln(priors->obs_prior_k[c]);
	}


	return logprior;
}


/////////////////////////////////////////////////////////
// 3.5 Robbins-Monro algorithm for setting acceptance rate.
//
// Robbins-Munroe stochastic approximation adjustment
//	return adjustment
//	i iteration number
//	a scaling factor for adjustment size
//	m number of iterations when adjustment halves from value at i=0
//	d outcome at iteration i (e.g. calculated acceptance probability (min of MH ratio and 1),
//								or acceptance or no acceptance)
//	p desired mean for d (e.g. desired acceptance probability)

double rm_scale(double step_scale, int step, int N_step_adapt, double log_prob)
{
	double dd = exp(log_prob);
	if (dd < -30) { dd = 0.0; }
	dd = std::min(dd, 1.0);

	double rm_temp = (dd - 0.23) / ((double(step) + 1) / (0.01*(double(N_step_adapt)) + 1));

	double out = step_scale*exp(rm_temp);

	out = std::max(out, 0.02);
	out = std::min(out, 5.0);

	return out;
}


/////////////////////////////////////////////////////////
// 3.6 Log gamma function, based on gamma.h from NRC3

double gammln(const double xx) {
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = { 57.1562356658629235, -59.5979603554754912,
		14.1360979747417471, -0.491913816097620199, 0.339946499848118887e-4,
		0.465236289270485756e-4, -0.983744753048795646e-4, 0.158088703224912494e-3,
		-0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3,
		0.844182239838527433e-4, -0.261908384015814087e-4, 0.368991826595316234e-5 };
	if (xx <= 0) throw("bad arg in gammln");
	y = x = xx;
	tmp = x + 5.242187500000000;
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j<14; j++) ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005*ser / x);
}

