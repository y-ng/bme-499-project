 /****************************************************
  *  wpparc PPA 2 - clusters log Leyton.c            *
  *                                                  *
  *  WEAVER++/ARC application to the three clusters  *
  *  of the logopenic variant of primary progressive *
  *  aphasia observed by Leyton et al. (2015)        *
  *                                                  *
  *  Simulation of profiles of logopenic clusters    *
  *                                                  *
  *  Conducted in 2020-2021                          *
  *                                                  *
  *  Ardi Roelofs, DCC, Radboud University           *
  *                                                  *
  ****************************************************/

  /*

  Simulations reported in:
  Roelofs, A. (2022). A neurocognitive computational account of
  word production, comprehension, and repetition in primary progressive
  aphasia. Brain and Language, 227, 105094.

  Three logopenic clusters reported in:
  Leyton, C. E., Hodges, J. R., McLean, C. A., Kril, J. J., Piguet, O., 
  & Ballard, K. J. (2015). Is the logopenic-variant of primary progressive
  aphasia a unitary disorder? Cortex, 67, 122�133.

*/


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>


#define STEP_SIZE 25   /* duration time step in ms */
#define N_STEPs 80     /* 2000 ms in total */
#define N_CONCEPTs 5   
#define N_LEMMAs 5     
#define N_MORPHEMEs 5  
#define N_PHONEMEs 10   
#define N_SYLLABLEs 5  

#define N_lesion_values 100 /* for 100 for weight lesion, 66 (!) for decay lesion */

#define N_GROUPs 4 /* Normal, Logopenic Clusters 1 to 3 */
#define NORMAL 0
#define CLUSTER_1 1
#define CLUSTER_2 2
#define CLUSTER_3 3

#define N_TASKs 3 /* Naming, Comprehension, Repetition */
#define NAMING 0
#define COMPREHENSION 1
#define REPETITION 2

#define N_ASSESSMENTs 1
#define LEYTON 0

#define Y 1.0     /* connection present */
#define N 0.0     /* connection absent */

 /* labeling network nodes */
#define CAT 0
#define DOG 1
#define MAT 2
#define FOG 3
#define FISH 4

#define pK 0
#define pE 1
#define pT 2
#define pD 3
#define pO 4
#define pG 5
#define pM 6
#define pF 7
#define pI 8
#define pS 9

#define Cat 0
#define Dog 1
#define Mat 2
#define Fog 3
#define Fish 4


 /* connections conceptual stratum */
 double  CC_con[N_CONCEPTs][N_CONCEPTs] =  {
  	          /* CAT   DOG  MAT  FOG  FISH  */
 /* CAT  */   {   N,    Y,   N,   N,    Y },
 /* DOG  */   {   Y,    N,   N,   N,    Y },
 /* MAT  */   {   N,    N,   N,   N,    N },
 /* FOG  */   {   N,    N,   N,   N,    N },
 /* FISH */   {   Y,    Y,   N,   N,    N }
 };


 /* connections between concept and lemma nodes */
 double  CL_con[N_CONCEPTs][N_LEMMAs] = {
    { Y,  N,  N,  N,  N },
    { N,  Y,  N,  N,  N },
    { N,  N,  Y,  N,  N },
    { N,  N,  N,  Y,  N },
    { N,  N,  N,  N,  Y }
 };


 /* connections between lemma nodes and morpheme nodes */
 double  LM_con[N_LEMMAs][N_MORPHEMEs] = {
    { Y,  N,  N,  N,  N },
    { N,  Y,  N,  N,  N },
    { N,  N,  Y,  N,  N },
    { N,  N,  N,  Y,  N },
    { N,  N,  N,  N,  Y }
 };

 /* connections between morpheme nodes and output phoneme nodes */
 double MP_con[N_MORPHEMEs][N_PHONEMEs] = {
                  /* K  E  T  D  O  G  M  F  I  S */   
 /* <cat>   */   {   Y, Y, Y, N, N, N, N, N, N, N },
 /* <dog>   */   {   N, N, N, Y, Y, Y, N, N, N, N },
 /* <mat>   */   {   N, Y, Y, N, N, N, Y, N, N, N },
 /* <fog>   */   {   N, N, N, N, Y, Y, N, Y, N, N },
 /* <fish>  */   {   N, N, N, N, N, N, N, Y, Y, Y }

 };


 /* connections between output phoneme nodes and syllable program nodes */
 double PS_con[N_PHONEMEs][N_SYLLABLEs] = {
	        /* Cat Dog  Mat Fog  Fish */
	 /* K */ { Y,   N,   N,  N,   N },
	 /* E */ { Y,   N,   Y,  N,   N },
	 /* T */ { Y,   N,   Y,  N,   N },
	 /* D */ { N,   Y,   N,  N,   N },
	 /* O */ { N,   Y,   N,  Y,   N },
	 /* G */ { N,   Y,   N,  Y,   N },
	 /* M */ { N,   N,   Y,  N,   N },
	 /* F */ { N,   N,   N,  Y,   Y },
	 /* I */ { N,   N,   N,  N,   Y },
	 /* S */ { N,   N,   N,  N,   Y }
 };


 /* connections between input and output phoneme nodes */
 double PP_con[N_PHONEMEs][N_PHONEMEs] = {
	        /* K    E    T   D    O    G   M   F   I  S */
	 /* K */ { Y,   N,   N,  N,   N,   N,  N,  N,  N, N },
	 /* E */ { N,   Y,   N,  N,   N,   N,  N,  N,  N, N  },
	 /* T */ { N,   N,   Y,  N,   N,   N,  N,  N,  N, N  },
	 /* D */ { N,   N,   N,  Y,   N,   N,  N,  N,  N, N  },
	 /* O */ { N,   N,   N,  N,   Y,   N,  N,  N,  N, N  },
	 /* G */ { N,   N,   N,  N,   N,   Y,  N,  N,  N, N  },
	 /* M */ { N,   N,   N,  N,   N,   N,  Y,  N,  N, N  },
	 /* F */ { N,   N,   N,  N,   N,   N,  N,  Y,  N, N  },
	 /* I */ { N,   N,   N,  N,   N,   N,  N,  N,  Y, N  },
	 /* S */ { N,   N,   N,  N,   N,   N,  N,  N,  N, Y  }
 };


 /* connections between input phoneme nodes and input morpheme nodes */
 double PiM_con[N_PHONEMEs][N_MORPHEMEs] = {
	        /* Cat Dog  Mat Fog  Fish */
	 /* K */ { Y,   N,   N,  N,   N },
	 /* E */ { Y,   N,   Y,  N,   N },
	 /* T */ { Y,   N,   Y,  N,   N },
	 /* D */ { N,   Y,   N,  N,   N },
	 /* O */ { N,   Y,   N,  Y,   N },
	 /* G */ { N,   Y,   N,  Y,   N },
	 /* M */ { N,   N,   Y,  N,   N },
	 /* F */ { N,   N,   N,  Y,   Y },
	 /* I */ { N,   N,   N,  N,   Y },
	 /* S */ { N,   N,   N,  N,   Y }
 };



 /* connections between input morpheme and output morpheme nodes */
 double  iMM_con[N_MORPHEMEs][N_MORPHEMEs] = {
    { Y,  N,  N,  N,  N },
    { N,  Y,  N,  N,  N },
    { N,  N,  Y,  N,  N },
    { N,  N,  N,  Y,  N },
    { N,  N,  N,  N,  Y }
 };

 /* connections between input morpheme and lemma nodes */
 double  iML_con[N_MORPHEMEs][N_LEMMAs] = {
    { Y,  N,  N,  N,  N },
    { N,  Y,  N,  N,  N },
    { N,  N,  Y,  N,  N },
    { N,  N,  N,  Y,  N },
    { N,  N,  N,  N,  Y }
 };


/* Real data Leyton et al. (2015) */
double REAL_DATA_LEYTON[N_GROUPs][N_TASKs] = {
                /* Naming  Comprehension Repetition */
/* Control       */ { 92.7,      97.3,      94.0 },
/* Cluster 1     */ { 67.3,      93.0,      93.3 },
/* Cluster 2	 */ { 29.0,      77.7,      93.3 },
/* Cluster 3     */ { 32.0,      85.3,      52.7 }
 };



double REAL_DATA[N_GROUPs][N_TASKs];
double SIM_DATA[N_GROUPs][N_TASKs];
double GOODNESS_OF_FIT[N_lesion_values];

double WEIGHT_value[N_lesion_values];
double DECAY_value[N_lesion_values];


 /* concept and lemma */
 double C_node_act[N_CONCEPTs], L_node_act[N_LEMMAs];
 /* output lexeme */
 double M_node_act[N_MORPHEMEs], oP_node_act[N_PHONEMEs], S_node_act[N_SYLLABLEs];
 /* input phonemes */
 double iM_node_act[N_MORPHEMEs], iP_node_act[N_PHONEMEs];


 /* input buffer */
 double input_C[N_CONCEPTs];
 double input_L[N_LEMMAs];
 double input_M[N_MORPHEMEs];
 double input_iM[N_MORPHEMEs];
 double input_iP[N_PHONEMEs];
 double input_oP[N_PHONEMEs];
 double input_S[N_SYLLABLEs];


 int T, step;     /* time in ms, step */
 int assessment;
 int group, task, lesion_value;

 /* parameter values */
 int    CYCLE_TIME = 25;                 /* ms per link */
 double SEM_rate = 0.0101 * STEP_SIZE;   /* prop per step_size ms */
 double LEM_rate = 0.0074 * STEP_SIZE;   /* prop per step_size ms */
 double LEX_rate = 0.0120 * STEP_SIZE;   /* prop per step_size ms */
 double DECAY_rate = 0.0240 * STEP_SIZE; /* prop per step_size ms */
 double EXTIN = 0.1965 * STEP_SIZE;      /* act_units per step_size ms */
 double LEMLEXFRAC = 0.3; 
  /* fraction of LEX_rate spread between lemmas and morphemes */
  /* implementing weak cascading of activation, see Roelofs (2008, JEP:LMC) */ 

 double FR = 0.10;  /* fraction of connection weight for input phoneme 
					   to morpheme and lemma, see Roelofs (1997, Cognition) */
 int    SEGMENT_DURATION = 125;  /* ms */
 int    PICTURE_DURATION = 125;  /* 125 ms */


 int WEIGHT_LESION = 1;
 int DECAY_LESION = 0;
 
 int SHOW_RESULTS_ALL_VALUES = 0;

/* Aphasia parameters */

/* weight lesion */

 double CONNECTION_DECREASE_CLUSTER_1;
 /* connections to and from lexical output forms,
   and between input and output phonemes */
 
 double CONNECTION_DECREASE_CLUSTER_2;
 /* connections to and from lexical output forms,
   and between lemmas and concepts */
 
 double CONNECTION_DECREASE_CLUSTER_3;
 /* connections to and from lexical output forms,
	and between input and output phonemes */


/* decay lesion */

 double DECAY_INCREASE_CLUSTER_1;
 /* lexical output forms and output phonemes */

 double DECAY_INCREASE_CLUSTER_2;
 /* lexical output forms and lemmas */

 double DECAY_INCREASE_CLUSTER_3;
 /* lexical output forms and output phonemes */



double ACT_C[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs];
double ACT_S[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs];

/* Activation target concept */
double ACT_CT[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation conceptual alternative */
double ACT_CR[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation target lemma */
double ACT_LT[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation lemma semantic alternative */
double ACT_LR[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation target syllable */
double ACT_ST[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation syllabic alternative */
double ACT_SR[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 


double TOTAL_ACT_C[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_C[N_lesion_values][N_GROUPs][N_TASKs];
double TOTAL_ACT_S[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_S[N_lesion_values][N_GROUPs][N_TASKs];

/* T = target, R = relative */
double TOTAL_ACT_CT[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_CT[N_lesion_values][N_GROUPs][N_TASKs];
double TOTAL_ACT_CR[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_CR[N_lesion_values][N_GROUPs][N_TASKs];

double TOTAL_ACT_LT[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_LT[N_lesion_values][N_GROUPs][N_TASKs];
double TOTAL_ACT_LR[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_LR[N_lesion_values][N_GROUPs][N_TASKs];

double TOTAL_ACT_ST[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_ST[N_lesion_values][N_GROUPs][N_TASKs];
double TOTAL_ACT_SR[N_lesion_values][N_GROUPs][N_TASKs];
double MEAN_ACT_SR[N_lesion_values][N_GROUPs][N_TASKs];

void set_real_data_matrix();
void reset_network();
void update_network();
void set_spreading_rates();
void set_input_to_zero();
void get_external_input();
void get_internal_input();
void update_activation_of_nodes();
void print_heading();
void print_parameters();
void compute_fits_and_print_results_on_screen();
void set_aphasic_parameters();
void compute_activation_results();
void determine_activation_critical_nodes();


/*****************
 * MAIN ROUTINES *
 *****************/

int main()
 {
	double ls;

    print_heading();

	print_parameters();

	set_spreading_rates();

	if (WEIGHT_LESION)
		for (lesion_value = 0, ls = 0.0; lesion_value < N_lesion_values; lesion_value++, ls += 0.01)
			WEIGHT_value[lesion_value] = ls;

	if (DECAY_LESION)
		for (lesion_value = 0, ls = 1.01; lesion_value < N_lesion_values; lesion_value++, ls += 0.01)
			DECAY_value[lesion_value] = ls;

	for (assessment = 0; assessment < N_ASSESSMENTs; assessment++) {

		set_real_data_matrix();

		for (group = 0; group < N_GROUPs; group++) {

			for (task = 0; task < N_TASKs; task++) {

				for (lesion_value = 0; lesion_value < N_lesion_values; lesion_value++) {

					reset_network();

					set_aphasic_parameters();

					for (T = 0, step = 0; T < (N_STEPs * STEP_SIZE); T += STEP_SIZE, step++) {

						update_network();
						determine_activation_critical_nodes();

					}
				}
			}

		}

			compute_activation_results();


			compute_fits_and_print_results_on_screen();
		
			
		getchar();
		
	}

	return 0;
 }



void set_real_data_matrix()
{
	
	for (group = 0; group < N_GROUPs; group++)
		for (task = 0; task < N_TASKs; task++)
			if (assessment == LEYTON)
				REAL_DATA[group][task] = REAL_DATA_LEYTON[group][task];
	
}

 void set_spreading_rates()
 {
   int i,j;
   
for(lesion_value=0; lesion_value < N_lesion_values; lesion_value++) 
 for(group=0; group < N_GROUPs; group++)
  for(task=0; task < N_TASKs; task++) 
	 for(step=0; step < N_STEPs; step++) {
	  ACT_C[lesion_value][step][group][task] = 0.0;
	  ACT_S[lesion_value][step][group][task] = 0.0;
      ACT_CT[lesion_value][step][group][task] = 0.0;
      ACT_CR[lesion_value][step][group][task] = 0.0;
      ACT_LT[lesion_value][step][group][task] = 0.0;
      ACT_LR[lesion_value][step][group][task] = 0.0;
      ACT_ST[lesion_value][step][group][task] = 0.0;
      ACT_SR[lesion_value][step][group][task] = 0.0;
	 }
 

  for(i=0;i<N_CONCEPTs;i++)
     for(j=0;j<N_CONCEPTs;j++) 
	CC_con[i][j]*=(SEM_rate);

  for(i=0;i<N_CONCEPTs;i++)
     for(j=0;j<N_LEMMAs;j++) 
	CL_con[i][j]*=LEM_rate;

	 
   for(i=0;i<N_LEMMAs;i++)
     for(j=0;j<N_MORPHEMEs;j++) 
	LM_con[i][j]*=LEX_rate;
	 

  for(i=0;i<N_MORPHEMEs;i++)
     for(j=0;j<N_PHONEMEs;j++) 
	MP_con[i][j]*=LEX_rate;

  for(i=0;i<N_PHONEMEs;i++)
     for(j=0;j<N_SYLLABLEs;j++) 
	PS_con[i][j]*=LEX_rate;


  /* connections for input phonemes to output phonemes, input morphemes, 
     input morphemes to lemmas */
  for(i=0;i<N_PHONEMEs;i++)
     for(j=0;j<N_PHONEMEs;j++) 
	PP_con[i][j]*=LEX_rate;

  for(i=0;i<N_PHONEMEs;i++)
     for(j=0;j<N_MORPHEMEs;j++) 
	PiM_con[i][j]*= (FR * LEX_rate);


  for(i=0;i<N_MORPHEMEs;i++)
     for(j=0;j<N_MORPHEMEs;j++) 
	iMM_con[i][j]*=LEX_rate;

  for(i=0;i<N_MORPHEMEs;i++)
     for(j=0;j<N_LEMMAs;j++) 
	iML_con[i][j]*=LEX_rate;

  }




 void reset_network()
 {
   int i;


   for(i=0;i<N_CONCEPTs;i++) 
     C_node_act[i]=0.0;

   for(i=0;i<N_LEMMAs;i++) 
     L_node_act[i]=0.0;

   for(i=0;i<N_MORPHEMEs;i++)
     M_node_act[i]=0.0;

   for(i=0;i<N_MORPHEMEs;i++) 
     iM_node_act[i]=0.0;
  
   for(i=0;i<N_PHONEMEs;i++)
     iP_node_act[i]=0.0;

   for(i=0;i<N_PHONEMEs;i++)
     oP_node_act[i]=0.0;
   
   for(i=0;i<N_SYLLABLEs;i++) 
     S_node_act[i]=0.0;

 }


void set_aphasic_parameters()
{

	double WEIGHT_FACTOR;
	double DECAY_FACTOR;

  if(WEIGHT_LESION)
     WEIGHT_FACTOR = WEIGHT_value[lesion_value];  
  else
     WEIGHT_FACTOR = 1.0;

  if (DECAY_LESION)
	  DECAY_FACTOR = DECAY_value[lesion_value];   
  else
	 DECAY_FACTOR = 1.0;


  /* setting of weight parameters */

  if(group == CLUSTER_1) 
     CONNECTION_DECREASE_CLUSTER_1 = WEIGHT_FACTOR;
  else 
     CONNECTION_DECREASE_CLUSTER_1 = 1.0; /* normal */
  
  if (group == CLUSTER_2)
	  CONNECTION_DECREASE_CLUSTER_2 = WEIGHT_FACTOR;
  else 
	  CONNECTION_DECREASE_CLUSTER_2 = 1.0; /* normal */
  
  if (group == CLUSTER_3)
	  CONNECTION_DECREASE_CLUSTER_3 = WEIGHT_FACTOR;
  else 
	  CONNECTION_DECREASE_CLUSTER_3 = 1.0; /* normal */


	/* setting of decay parameters */

  if (group == CLUSTER_1)
	  DECAY_INCREASE_CLUSTER_1 = DECAY_FACTOR;
  else
	  DECAY_INCREASE_CLUSTER_1 = 1.0; /* normal */

  if (group == CLUSTER_2)
	  DECAY_INCREASE_CLUSTER_2 = DECAY_FACTOR;
  else
	  DECAY_INCREASE_CLUSTER_2 = 1.0; /* normal */

  if (group == CLUSTER_3)
	  DECAY_INCREASE_CLUSTER_3 = DECAY_FACTOR;
  else
	  DECAY_INCREASE_CLUSTER_3 = 1.0; /* normal */

}



/*****************************
 * NETWORK UPDATING ROUTINES *
 *****************************/

 void update_network()
 {
   set_input_to_zero();
   get_external_input();
   get_internal_input();
   update_activation_of_nodes();
 }



 void set_input_to_zero()
 {
   int i;

   for(i=0;i<N_CONCEPTs;i++) 
       input_C[i]=0.0;

   for(i=0;i<N_LEMMAs;i++)
       input_L[i]=0.0;

   for(i=0;i<N_MORPHEMEs;i++)
       input_M[i]=0.0;

   for(i=0;i<N_MORPHEMEs;i++) 
       input_iM[i]=0.0;
   
   for(i=0;i<N_PHONEMEs;i++) 
       input_iP[i]=0.0;

   for(i=0;i<N_PHONEMEs;i++) 
       input_oP[i]=0.0;

   for(i=0;i<N_SYLLABLEs;i++) 
       input_S[i]=0.0;
 

  }


 void get_external_input()
 {


   if(task == NAMING) {

    /* picture input */
	   if (T >= 0 && T < PICTURE_DURATION) 
		   input_C[CAT] += EXTIN;

    /* enhancement */
  	  if( (T >= (0 + 1 * CYCLE_TIME)) 
		         &&  (T < (1 * CYCLE_TIME + PICTURE_DURATION))  ) 
	     input_C[CAT]+= EXTIN;
	}
   
 
    if(task == COMPREHENSION || task == REPETITION) {
 
    /* spoken word input */
      if((0 <= T) && (T < SEGMENT_DURATION)) {
	         input_iP[pK]+=EXTIN ;
             }

      if((SEGMENT_DURATION <= T) 
          && (T < (2 * SEGMENT_DURATION))) {
	         input_iP[pE]+=EXTIN ;
             }


     if((2 * SEGMENT_DURATION <= T) 
         && (T < (3 * SEGMENT_DURATION))) {
	         input_iP[pT]+=EXTIN ;
             }
 
      }

 }


 void get_internal_input()
 {
   int i,j;


 /* input activation for concept nodes */

  for(i=0;i<N_CONCEPTs;i++)
	  for(j=0;j<N_CONCEPTs;j++) 
       input_C[i]+=(C_node_act[j] 
	                * (CC_con[j][i] ) );

  for(i=0;i<N_CONCEPTs;i++)
     for(j=0;j<N_LEMMAs;j++)
       input_C[i]+=(L_node_act[j] 
	              * CL_con[j][i] * CONNECTION_DECREASE_CLUSTER_2 );


 /* input activation for lemma nodes */
  for(i=0;i<N_LEMMAs;i++)
	  for(j=0;j<N_CONCEPTs;j++) 
        input_L[i]+=( C_node_act[j] * CL_con[j][i] 
	                  *  CONNECTION_DECREASE_CLUSTER_2 );

  for(i=0;i<N_LEMMAs;i++)
     for(j=0;j<N_MORPHEMEs;j++)                      
        input_L[i]+=( iM_node_act[j] * iML_con[j][i]);


 /* input activation for output morpheme nodes */
   for(i=0;i<N_MORPHEMEs;i++)
     for(j=0;j<N_LEMMAs;j++) 
	    input_M[i]+=( L_node_act[j] * LEMLEXFRAC * LM_con[j][i] 
			* CONNECTION_DECREASE_CLUSTER_1 * CONNECTION_DECREASE_CLUSTER_2 * CONNECTION_DECREASE_CLUSTER_3);


  for(i=0;i<N_MORPHEMEs;i++)
     for(j=0;j<N_MORPHEMEs;j++)
        input_M[i]+=( iM_node_act[j] * iMM_con[j][i] 
			* CONNECTION_DECREASE_CLUSTER_1 * CONNECTION_DECREASE_CLUSTER_2 * CONNECTION_DECREASE_CLUSTER_3);



 /* input activation for output phoneme nodes */
  for(i=0;i<N_PHONEMEs;i++)
     for(j=0;j<N_MORPHEMEs;j++) 
	  input_oP[i]+=( M_node_act[j] * MP_con[j][i] 
		  * CONNECTION_DECREASE_CLUSTER_1 * CONNECTION_DECREASE_CLUSTER_2 * CONNECTION_DECREASE_CLUSTER_3);

 for (i = 0; i < N_PHONEMEs; i++)
	  for (j = 0; j < N_PHONEMEs; j++)
		input_oP[i] += (iP_node_act[j] * PP_con[j][i] * CONNECTION_DECREASE_CLUSTER_1 * CONNECTION_DECREASE_CLUSTER_3);


 /* input activation for syllable program nodes */
   for(i=0;i<N_SYLLABLEs;i++)
     for(j=0;j<N_PHONEMEs;j++)  
	input_S[i]+=( oP_node_act[j] * PS_con[j][i] );



 /* input activation for input phoneme nodes */
  for(i=0;i<N_PHONEMEs;i++)
     for(j=0;j<N_PHONEMEs;j++) 
	input_iP[i]+=( oP_node_act[j] * PP_con[j][i] * CONNECTION_DECREASE_CLUSTER_1 * CONNECTION_DECREASE_CLUSTER_3);
 
 /* input activation for input morpheme nodes */
	  for(i=0;i<N_MORPHEMEs;i++)
        for(j=0;j<N_PHONEMEs;j++) 
	       input_iM[i]+=( iP_node_act[j] * PiM_con[j][i]);
 
 }



 void update_activation_of_nodes()
 {
   int i;

   for(i=0;i<N_CONCEPTs;i++)
     C_node_act[i]=((C_node_act[i] 
                     * (1.0 - DECAY_rate ) ) + input_C[i]);

   for(i=0;i<N_LEMMAs;i++) 
     L_node_act[i]=((L_node_act[i] * (1.0 - (DECAY_rate * DECAY_INCREASE_CLUSTER_2) )) + input_L[i]);

    
   for(i=0;i<N_MORPHEMEs;i++)
       M_node_act[i]=((M_node_act[i] * (1.0 - (DECAY_rate * DECAY_INCREASE_CLUSTER_1 
		              * DECAY_INCREASE_CLUSTER_2 * DECAY_INCREASE_CLUSTER_3) ) ) + input_M[i]);

     for(i=0;i<N_PHONEMEs;i++) 
       oP_node_act[i]=((oP_node_act[i] * (1.0 - (DECAY_rate * DECAY_INCREASE_CLUSTER_1
		   * DECAY_INCREASE_CLUSTER_3) ) ) + input_oP[i]);
       
     for(i=0;i<N_PHONEMEs;i++) 
       iP_node_act[i]=((iP_node_act[i] * (1.0 - DECAY_rate)) + input_iP[i]);

   for(i=0;i<N_MORPHEMEs;i++)
       iM_node_act[i]=((iM_node_act[i] * (1.0 - DECAY_rate)) + input_iM[i]);
   
   for(i=0;i<N_SYLLABLEs;i++) 
       S_node_act[i]=((S_node_act[i] * (1.0 - DECAY_rate)) + input_S[i]);


 }


void determine_activation_critical_nodes()
{
	ACT_C[lesion_value][step][group][task] = C_node_act[CAT];
	ACT_S[lesion_value][step][group][task] = S_node_act[CAT];
	ACT_CT[lesion_value][step][group][task] = C_node_act[CAT];
	ACT_CR[lesion_value][step][group][task] = C_node_act[DOG];
	ACT_LT[lesion_value][step][group][task] = L_node_act[CAT];
	ACT_LR[lesion_value][step][group][task] = L_node_act[DOG];
	ACT_ST[lesion_value][step][group][task] = S_node_act[CAT];
	ACT_SR[lesion_value][step][group][task] = S_node_act[MAT];

}


void compute_activation_results()
{

  int i;

  for(lesion_value=0; lesion_value < N_lesion_values; lesion_value++)
	for(group=0; group < N_GROUPs; group++) 
		for(task=0; task < N_TASKs; task++) {
	  TOTAL_ACT_C[lesion_value][group][task] = 0.0;
	  TOTAL_ACT_S[lesion_value][group][task] = 0.0;
	  MEAN_ACT_C[lesion_value][group][task] = 0.0;
	  MEAN_ACT_S[lesion_value][group][task] = 0.0;

	  TOTAL_ACT_CT[lesion_value][group][task] = 0.0;
	  TOTAL_ACT_CR[lesion_value][group][task] = 0.0;
	  MEAN_ACT_CT[lesion_value][group][task] = 0.0;
	  MEAN_ACT_CR[lesion_value][group][task] = 0.0;

	  TOTAL_ACT_LT[lesion_value][group][task] = 0.0;
	  TOTAL_ACT_LR[lesion_value][group][task] = 0.0;
	  MEAN_ACT_LT[lesion_value][group][task] = 0.0;
	  MEAN_ACT_LR[lesion_value][group][task] = 0.0;

	  TOTAL_ACT_ST[lesion_value][group][task] = 0.0;
	  TOTAL_ACT_SR[lesion_value][group][task] = 0.0;
	  MEAN_ACT_ST[lesion_value][group][task] = 0.0;
	  MEAN_ACT_SR[lesion_value][group][task] = 0.0;

	}
  
 for(lesion_value=0; lesion_value < N_lesion_values; lesion_value++)
   for(group=0; group < N_GROUPs; group++) 
	   for(task=0; task < N_TASKs; task++) {
         for(i = 0; i < N_STEPs; i++){ 
 	       TOTAL_ACT_C[lesion_value][group][task] 
		     += ACT_C[lesion_value][i][group][task];

	       TOTAL_ACT_S[lesion_value][group][task] 
		     += ACT_S[lesion_value][i][group][task];
	       TOTAL_ACT_CT[lesion_value][group][task] 
		     += ACT_CT[lesion_value][i][group][task];
	       TOTAL_ACT_CR[lesion_value][group][task] 
		     += ACT_CR[lesion_value][i][group][task];
	       TOTAL_ACT_LT[lesion_value][group][task] 
		     += ACT_LT[lesion_value][i][group][task];
	       TOTAL_ACT_LR[lesion_value][group][task] 
		     += ACT_LR[lesion_value][i][group][task];
	       TOTAL_ACT_ST[lesion_value][group][task] 
		     += ACT_ST[lesion_value][i][group][task];
	       TOTAL_ACT_SR[lesion_value][group][task] 
		     += ACT_SR[lesion_value][i][group][task];
	     }

	     MEAN_ACT_C[lesion_value][group][task] 
		   = (TOTAL_ACT_C[lesion_value][group][task] / N_STEPs);


	     MEAN_ACT_S[lesion_value][group][task] 
		   = (TOTAL_ACT_S[lesion_value][group][task] / N_STEPs);
	     MEAN_ACT_CT[lesion_value][group][task] 
		   = (TOTAL_ACT_CT[lesion_value][group][task] / N_STEPs);
	     MEAN_ACT_CR[lesion_value][group][task] 
		   = (TOTAL_ACT_CR[lesion_value][group][task] / N_STEPs);
	     MEAN_ACT_LT[lesion_value][group][task] 
		   = (TOTAL_ACT_LT[lesion_value][group][task] / N_STEPs);
	     MEAN_ACT_LR[lesion_value][group][task] 
		   = (TOTAL_ACT_LR[lesion_value][group][task] / N_STEPs);
	     MEAN_ACT_ST[lesion_value][group][task] 
		   = (TOTAL_ACT_ST[lesion_value][group][task] / N_STEPs);
	     MEAN_ACT_SR[lesion_value][group][task] 
		   = (TOTAL_ACT_SR[lesion_value][group][task] / N_STEPs);
  
	   }
}




/*********************
 * FITS AND PRINTING *
 *********************/

 void print_heading()
 {
	 printf("\n");
	 printf("WEAVER++/ARC model simulation of primary progressive aphasia (c) Ardi Roelofs\n");
	 printf("Simulating clusters of logopenic variant\n");
	 printf("Empirical data on clusters from Leyton et al. (2015)\n");

 }


 void print_parameters()
 {
	 

     printf("Parameter values:\n");
     printf("cycle time : %6d [ms]\n", CYCLE_TIME);
     printf("sem_rate   : %.4f [prop/ms]\n",SEM_rate/STEP_SIZE);
     printf("lem_rate   : %.4f [prop/ms]\n",LEM_rate/STEP_SIZE);
     printf("exin       : %.4f [act_units/ms]\n",EXTIN/STEP_SIZE);
     printf("d          : %.4f [prop/ms]\n",DECAY_rate/STEP_SIZE);


    printf("press any key to continue ");

    getchar();

 }


 void compute_fits_and_print_results_on_screen()
 {

	 double LV; /* lesion value */
	 int i, j, a;
  
	 for (i = 0; i < N_GROUPs; i++)
		 for (j = 0; j < N_TASKs; j++)
			 SIM_DATA[i][j] = 0.0;

	 for (i = 0; i < N_lesion_values; i++)
		 GOODNESS_OF_FIT[i] = 0.0;

	 if (assessment == LEYTON)
		 printf("\nAssessment is Leyton et al. (2015) \n");


   for(group=0; group <  N_GROUPs; group++) {
        printf(" \n");

        if(group == NORMAL) 
        printf("NORMAL \n");
		else if (group == CLUSTER_1)
			printf("CLUSTER 1 \n");
		else if (group == CLUSTER_2)
			printf("CLUSTER 2  \n");
		else if (group == CLUSTER_3)
			printf("CLUSTER 3  \n");


		printf("        Naming   Comprehension  Repetition \n");
		printf("Real:   %5.2f         %5.2f        %5.2f \n", 
		  REAL_DATA[group][NAMING], REAL_DATA[group][COMPREHENSION], REAL_DATA[group][REPETITION]);
		printf("Lesion:                                    MAE\n");


		for (lesion_value = 0; lesion_value < N_lesion_values; lesion_value++) {

			SIM_DATA[group][NAMING] = (MEAN_ACT_ST[lesion_value][group][NAMING]
				- MEAN_ACT_SR[lesion_value][group][NAMING])
				/ (MEAN_ACT_ST[lesion_value][NORMAL][NAMING]
					- MEAN_ACT_SR[lesion_value][NORMAL][NAMING]) * 100.0;

			SIM_DATA[group][COMPREHENSION] = (MEAN_ACT_CT[lesion_value][group][COMPREHENSION]
				- MEAN_ACT_CR[lesion_value][group][COMPREHENSION])
				/ (MEAN_ACT_CT[lesion_value][NORMAL][COMPREHENSION]
					- MEAN_ACT_CR[lesion_value][NORMAL][COMPREHENSION]) * 100.0;

			SIM_DATA[group][REPETITION] = (MEAN_ACT_ST[lesion_value][group][REPETITION]
				- MEAN_ACT_SR[lesion_value][group][REPETITION])
				/ (MEAN_ACT_ST[lesion_value][NORMAL][REPETITION]
					- MEAN_ACT_SR[lesion_value][NORMAL][REPETITION]) * 100.0;


			if (group == NORMAL)
				LV = 1.0;
			else if (WEIGHT_LESION)
				LV = WEIGHT_value[lesion_value];
			else if (DECAY_LESION)
				LV = DECAY_value[lesion_value];

			GOODNESS_OF_FIT[lesion_value] = (fabs(REAL_DATA[group][NAMING] - SIM_DATA[group][NAMING])
				+ fabs(REAL_DATA[group][COMPREHENSION] - SIM_DATA[group][COMPREHENSION])
				+ fabs(REAL_DATA[group][REPETITION] - SIM_DATA[group][REPETITION])) / 3.0;

			if (SHOW_RESULTS_ALL_VALUES)
			printf("%5.2f   %5.2f        %5.2f        %5.2f     %5.2f\n",
				LV, SIM_DATA[group][NAMING], SIM_DATA[group][COMPREHENSION], SIM_DATA[group][REPETITION],
				
				(fabs( REAL_DATA[group][NAMING] - SIM_DATA[group][NAMING])
					+ fabs( REAL_DATA[group][COMPREHENSION] - SIM_DATA[group][COMPREHENSION])
					+ fabs( REAL_DATA[group][REPETITION] - SIM_DATA[group][REPETITION])) / 3.0);

		}

		for (a = 0, i = 0; i < N_lesion_values; i++)
			if (GOODNESS_OF_FIT[a] > GOODNESS_OF_FIT[i])
				a = i;

		if (WEIGHT_LESION)
			printf("Best fit weight value = %.2f   MAE = %.2f\n", WEIGHT_value[a], GOODNESS_OF_FIT[a]);
		if (DECAY_LESION)
			printf("Best fit decay value = %.2f   MAE = %.2f\n", DECAY_value[a], GOODNESS_OF_FIT[a]);

		printf("Sim:   %5.2f         %5.2f        %5.2f \n",

			(MEAN_ACT_ST[a][group][NAMING]
				- MEAN_ACT_SR[a][group][NAMING])
			/ (MEAN_ACT_ST[a][NORMAL][NAMING]
				- MEAN_ACT_SR[a][NORMAL][NAMING]) * 100.0,

				(MEAN_ACT_CT[a][group][COMPREHENSION]
					- MEAN_ACT_CR[a][group][COMPREHENSION])
			/ (MEAN_ACT_CT[a][NORMAL][COMPREHENSION]
				- MEAN_ACT_CR[a][NORMAL][COMPREHENSION]) * 100.0,

				(MEAN_ACT_ST[a][group][REPETITION]
					- MEAN_ACT_SR[a][group][REPETITION])
			/ (MEAN_ACT_ST[a][NORMAL][REPETITION]
				- MEAN_ACT_SR[a][NORMAL][REPETITION]) * 100.0);
	 
   }
 }




 