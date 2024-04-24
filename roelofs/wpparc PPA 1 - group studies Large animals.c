/****************************************************
 *  wpparc PPA 1 - group studies Large animals.c    *
 *                                                  *
 *  WEAVER++/ARC application to the three variants  *
 *  of primary progressive aphasia (PPA)            *
 *                                                  *
 *  Basic profiles and progression of disease       *
 *                                                  *
 *  Simulation of group studies                     *
 *  Testing the effect of a larger network          *
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

English data for PPA group performance on single-word tasks:
Savage, S., Hsieh, S., Leslie, F., Foxe, D., Piguet, O., 
& Hodges, J.R. (2013). Distinguishing subtypes in primary progressive
aphasia: Application of the Sydney Language Battery. Dementia and
Geriatric Cognitive Disorders, 35, 208–218.

Dutch data for PPA group performance on single-word tasks:
Janssen, N., Roelofs, A., Van den Berg, E., Holleman, M. A., In de Braek, 
D. M. J. M., Piguet, O., Piai, V., & Kessels, R. P. C. (2022). 
The diagnostic value of language screening in primary progressive aphasia: 
Validation and application of the Sydney Language Battery. Journal
of Speech, Language, and Hearing Research, 65(1), 200–214.

Data for PPA disease progression, group performance on single-word tasks:

Brambati, S. M., Amici, S., Racine, C. A., Neuhaus, J., Miller, Z., Ogar, J., 
Dronkers, N., Miller, B. L., Rosen, H., & Gorno-Tempini, M. L. (2015). 
Longitudinal gray matter contraction in three variants of primary progressive
aphasia: A tenser-based morphometry study. NeuroImage: Clinical, 8, 345–355.

Rohrer, J. D., Caso, F., Mahoney, C., Henry, M., Rosen, H. J., Rabinovici, G., 
Rossor, M. N., Miller, B., Warren, J. D., Fox, N. C., Ridgway, G. R., & 
Gorno-Tempini, M. L. (2013). Patterns of longitudinal brain atrophy in the 
logopenic variant of primary progressive aphasia. Brain and Language, 
127(2), 121–126.

Mandelli, M. L., Vilaplana, E., Brown, J. A., Hubbard, H. I., Binney, R. J., 
Attygalle, S., Santos-Santos, M. A., Miller, Z. A., Pakvasa, M., Henry, M. L., 
Rosen, H. J., Henry, R. G., Rabinovici, G. D., Miller, B. L., Seeley, W. W., 
& Gorno-Tempini, M. L. (2016). Healthy brain connectivity predicts atrophy progression 
in non-fluent variant of primary progressive aphasia. Brain, 139(10), 2778–2791.

*/


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>



#define STEP_SIZE 25   /* duration time step in ms */
#define N_STEPs 80     /* 2000 ms in total */
#define N_CONCEPTs 12   
#define N_LEMMAs 12     
#define N_MORPHEMEs 12  
#define N_PHONEMEs 22   
#define N_SYLLABLEs 28  


#define N_lesion_values 100 /* for 100 for weight lesion, 66 (!) for decay lesion */

#define N_GROUPs 4 /* Normal, Nonfluent_agrammatic, Semantic_dementia, Logopenic */
#define NORMAL 0
#define NONFLUENT_AGRAMMATIC 1
#define SEMANTIC_DEMENTIA 2
#define LOGOPENIC 3

#define N_TASKs 3 /* Naming, Comprehension, Repetition */
#define NAMING 0
#define COMPREHENSION 1
#define REPETITION 2

#define N_ASSESSMENTs 6
#define ENGLISH 0
#define DUTCH 1
#define BRAMBATI_T1 2 /* baseline */
#define BRAMBATI_T2 3 /* follow up */
#define ROHRERMANDELLI_T1 4 /* baseline */
#define ROHRERMANDELLI_T2 5 /* follow up */

#define Y 1.0     /* connection present */
#define N 0.0     /* connection absent */




/* labeling network nodes */
/* target stethoscope ste θə skōp */
#define CAT 7 /* target */
#define DOG 8 /* relative */
#define pK 5 
#define pA 15
#define pT 14

#define sKAT 11 /* target */
#define sMAT 14 /* relative  */



 /* connections conceptual stratum */
double  CC_con[N_CONCEPTs][N_CONCEPTs] = {
	/*0		1	2	3	4	5	6	7	8	9	10	11 */
	/* 0	butterfly */	{ N,	Y,	Y,	Y,	Y,	Y,	Y,	Y,	Y,	N,	N,	Y	},
	/* 1	elephant */		{ Y,	N,	Y,	Y,	Y,	Y,	Y,	Y,	Y,	N,	N,	Y	},
	/* 2	caterpillar */	{ Y,	Y,	N,	Y,	Y,	Y,	Y,	Y,	Y,	N,	N,	Y	},
	/* 3	dinosaur */		{ Y,	Y,	Y,	N,	Y,	Y,	Y,	Y,	Y,	N,	N,	Y	},
	/* 4	rhinoceros */	{ Y,	Y,	Y,	Y,	N,	Y,	Y,	Y,	Y,	N,	N,	Y	},
	/* 5	hippopotamus */	{ Y,	Y,	Y,	Y,	Y,	N,	Y,	Y,	Y,	N,	N,	Y	},
	/* 6	orangutan */ 	{ Y,	Y,	Y,	Y,	Y,	Y,	N,	Y,	Y,	N,	N,	Y	},
	/* 7	cat */			{ Y,	Y,	Y,	Y,	Y,	Y,	Y,	N,	Y,	N,	N,	Y	},
	/* 8	dog */			{ Y,	Y,	Y,	Y,	Y,	Y,	Y,	Y,	Y,	N,	N,	Y	},
	/* 9	mat */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 10	fog */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 11	fish */			{ Y,	Y,	Y,	Y,	Y,	Y,	Y,	Y,	Y,	N,	N,	N	}
};


/* connections between concept and lemma nodes */
double  CL_con[N_CONCEPTs][N_LEMMAs] = {
	/*0		1	2	3	4	5	6	7	8	9	10	11 */
	/* 0	butterfly */	{ Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 1	elephant */		{ N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 2	caterpillar */	{ N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 3	dinosaur */		{ N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 4	rhinoceros */	{ N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/* 5	hippopotamus */	{ N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N	},
	/* 6	orangutan */ 	{ N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N	},
	/* 7	cat */			{ N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N	},
	/* 8	dog */			{ N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N	},
	/* 9	mat */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N	},
	/* 10	fog */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N	},
	/* 11	fish */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	}
};


/* connections between lemma nodes and morpheme nodes */
double  LM_con[N_LEMMAs][N_MORPHEMEs] = {
	/*0		1	2	3	4	5	6	7	8	9	10	11 */
	/* 0	butterfly */	{ Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 1	elephant */		{ N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 2	caterpillar */	{ N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 3	dinosaur */		{ N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 4	rhinoceros */	{ N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/* 5	hippopotamus */	{ N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N	},
	/* 6	orangutan */ 	{ N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N	},
	/* 7	cat */			{ N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N	},
	/* 8	dog */			{ N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N	},
	/* 9	mat */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N	},
	/* 10	fog */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N	},
	/* 11	fish */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	}
};



/* connections between morpheme nodes and output phoneme nodes */
double MP_con[N_MORPHEMEs][N_PHONEMEs] = {
	/*	b	d	f	g	h	k	l	m	n	ŋ	p	r	s	ʃ	t	a	ä	ə	e  	ī	i	ȯ	*/
	/* 0	bətərflī */		{	Y,	N,	Y,		N,	N,	Y,	N,	N,	N,	N,	Y,	N,	N,	Y,	N,	N,	Y,	N,	Y,	N,	N	},
	/* 1	eləfənt */		{	N,	N,	Y,	N,	N,	N,	Y,	N,	Y,	N,	N,	N,	N,	N,	Y,	N,	N,	Y,	Y,	N,	N,	N	},
	/* 2	katərpilər */	{	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	Y,	Y,	N,	N,	Y,	Y,	N,	Y,	N,	N,	Y,	N	},
	/* 3	dīnəsȯr */		{	N,	Y,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	Y,	Y,	N,	N,	N,	N,	Y,	N,	N,	Y,	Y	},
	/* 4	rīnäsrəs */		{	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	Y,	Y,	N,	N,	N,	Y,	Y,	N,	Y,	N,	N	},
	/* 5	hipəpätəməs */	{	N,	N,	N,	N,	Y,	N,	N,	Y,	N,	N,	Y,	N,	Y,	N,	Y,	N,	Y,	Y,	N,	N,	Y,	N	},
	/* 6	əraŋətaŋ */		{	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	Y,	N,	N,	Y,	Y,	N,	Y,	N,	N,	N,	N	},
	/* 7	kat */			{	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N	},
	/* 8	dȯg */			{	N,	Y,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	},
	/* 9	mat */			{	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N	},
	/* 10	fȯg */			{	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	},
	/* 11	fiʃ */			{	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	Y,	N	}
};



/* connections between output phoneme nodes and syllable program nodes */
double PS_con[N_PHONEMEs][N_SYLLABLEs] = {
	/*	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	*/
	/*	bə	dī  dȯg	ə	e  	fənt flī fȯg	fiʃ	hi	ka	kat	lə	lər	mat	məs	näs	nə	pä	pə	pi	raŋ	rəs	rī	sȯr	taŋ	tə	tər	*/
	/*	0	b	*/	{	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	1	d	*/	{	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	2	f	*/	{	N,	N,	N,	N,	N,	Y,	Y,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	3	g	*/	{	N,	N,	Y,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	4	h	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	5	k	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	6	l	*/	{	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	7	m	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	8	n	*/	{	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	9	ŋ	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	Y,	N,	N	},
	/*	10	p	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/*	11	r	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	Y,	Y,	N,	N,	Y	},
	/*	12	s	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	Y,	N,	Y,	N,	N,	N	},
	/*	13	ʃ	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	14	t	*/	{	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	Y,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	Y	},
	/*	15	a	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	Y,	N,	N	},
	/*	16	ä	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	17	ə	*/	{	N,	N,	N,	Y,	N,	Y,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	Y,	N,	Y,	N,	Y,	N,	N,	Y,	N,	N,	N,	Y,	Y	},
	/*	18	e  	*/	{	Y,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	19	ī	*/	{	N,	Y,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N	},
	/*	20	i	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/*	21	ȯ	*/	{	N,	N,	Y,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N	}
};




/* connections between input and output phoneme nodes */
double PP_con[N_PHONEMEs][N_PHONEMEs] = {
	/*	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	*/
	/*	0	b	*/	{	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	1	d	*/	{	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	2	f	*/	{	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	3	g	*/	{	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	4	h	*/	{	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	5	k	*/	{	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	6	l	*/	{	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	7	m	*/	{	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	8	n	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	9	ŋ	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	10	p	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	11	r	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	12	s	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	13	ʃ	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	14	t	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/*	15	a	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N	},
	/*	16	ä	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N	},
	/*	17	ə	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N	},
	/*	18	e  	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N	},
	/*	19	ī	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N	},
	/*	20	i	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N	},
	/*	21	ȯ	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	}
};



/* connections between input phoneme nodes and input morpheme nodes */
double PiM_con[N_PHONEMEs][N_MORPHEMEs] = {
	/*	0	1	2	3	4	5	6	7	8	9	10	11	*/
	/*	b	*/	{	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	d	*/	{	N,	N,	N,	Y,	N,	N,	N,	N,	Y,	N,	N,	N	},
	/*	f	*/	{	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	Y	},
	/*	g	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	Y,	N	},
	/*	h	*/	{	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N	},
	/*	k	*/	{	N,	N,	Y,	N,	N,	N,	N,	Y,	N,	N,	N,	N	},
	/*	l	*/	{	Y,	Y,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	m	*/	{	N,	N,	N,	N,	N,	Y,	N,	N,	N,	Y,	N,	N	},
	/*	n	*/	{	N,	Y,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/*	ŋ	*/	{	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N	},
	/*	p	*/	{	N,	N,	Y,	N,	N,	Y,	N,	N,	N,	N,	N,	N	},
	/*	r	*/	{	Y,	N,	Y,	Y,	Y,	N,	Y,	N,	N,	N,	N,	N	},
	/*	s	*/	{	N,	N,	N,	Y,	Y,	Y,	N,	N,	N,	N,	N,	N	},
	/*	ʃ	*/	{	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	},
	/*	t	*/	{	Y,	Y,	Y,	N,	N,	Y,	Y,	Y,	N,	Y,	N,	N	},
	/*	a	*/	{	N,	N,	Y,	N,	N,	N,	Y,	Y,	N,	Y,	N,	N	},
	/*	ä	*/	{	N,	N,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N	},
	/*	ə	*/	{	Y,	Y,	Y,	Y,	Y,	Y,	Y,	N,	N,	N,	N,	N	},
	/*	e  	*/	{	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/*	ī	*/	{	Y,	N,	N,	Y,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/*	i	*/	{	N,	N,	Y,	N,	N,	Y,	N,	N,	N,	N,	N,	Y	},
	/*	ȯ	*/	{	N,	N,	N,	Y,	N,	N,	N,	N,	Y,	N,	Y,	N	}
};



/* connections between input morpheme and output morpheme nodes */
double  iMM_con[N_MORPHEMEs][N_MORPHEMEs] = {
	/*0		1	2	3	4	5	6	7	8	9	10	11 */
	/* 0	butterfly */	{ Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 1	elephant */		{ N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 2	caterpillar */	{ N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 3	dinosaur */		{ N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 4	rhinoceros */	{ N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/* 5	hippopotamus */	{ N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N	},
	/* 6	orangutan */ 	{ N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N	},
	/* 7	cat */			{ N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N	},
	/* 8	dog */			{ N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N	},
	/* 9	mat */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N	},
	/* 10	fog */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N	},
	/* 11	fish */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	}
};

/* connections between input morpheme and lemma nodes */
double  iML_con[N_MORPHEMEs][N_LEMMAs] = {
	/*0		1	2	3	4	5	6	7	8	9	10	11 */
	/* 0	butterfly */	{ Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 1	elephant */		{ N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 2	caterpillar */	{ N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 3	dinosaur */		{ N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N,	N	},
	/* 4	rhinoceros */	{ N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N,	N	},
	/* 5	hippopotamus */	{ N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N,	N	},
	/* 6	orangutan */ 	{ N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N,	N	},
	/* 7	cat */			{ N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N,	N	},
	/* 8	dog */			{ N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N,	N	},
	/* 9	mat */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N,	N	},
	/* 10	fog */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y,	N	},
	/* 11	fish */			{ N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	N,	Y	}
};






/* English data on PPA for single word tasks: Savage et al. (2013) */
double REAL_DATA_ENGLISH[N_GROUPs][N_TASKs] = {
                        /* Naming  Comprehension Repetition */
/* Control              */ { 88.7,      97.0,      99.7 },
/* Nonfluent/agrammatic */ { 78.3,      94.3,      79.7 },
/* Semantic	            */ { 22.7,      63.3,      95.3 },
/* Logopenic            */ { 41.3,      84.7,      84.7 }
 };


/* Dutch data on PPA for single word tasks: Janssen et al. (2021) */
double REAL_DATA_DUTCH[N_GROUPs][N_TASKs] = {
	                          /* Naming  Comprehension Repetition */
	/* Control              */ { 90.3,      96.3,      96.7 },
	/* Nonfluent/agrammatic */ { 77.3,      97.7,      89.3 },
	/* Semantic	            */ { 29.0,      78.0,      96.3 },
	/* Logopenic            */ { 66.3,      93.7,      91.3 }
};


/* Brambati T1 data on PPA for single word tasks: Brambati et al. (2015) */
double REAL_DATA_BRAMBATI_T1[N_GROUPs][N_TASKs] = {
	                         /* Naming  Comprehension Repetition */
	/* Control              */ { 90.3,      96.3,      96.7 }, /* dummy, from Savage */
	/* Nonfluent/agrammatic */ { 85.3,      99.7,      83.7 },
	/* Semantic	            */ { 26.7,      88.0,      90.6 },
	/* Logopenic            */ { 69.3,      95.0,      69.0 }
};

/* Brambati T2 data on PPA for single word tasks: Brambati et al. (2015) */
double REAL_DATA_BRAMBATI_T2[N_GROUPs][N_TASKs] = {
	/* Naming  Comprehension Repetition */
	/* Control              */ { 90.3,      96.3,      96.7 }, /* dummy, from Savage */
	/* Nonfluent/agrammatic */ { 83.3,      94.8,      68.0 },
	/* Semantic	            */ { 19.3,      66.7,      82.3 },
	/* Logopenic            */ { 52.7,      95.0,      58.8 }
};


/*
Rohrer et al. (2013), logopenic patients (N=21), T1 baseline and T2 one year later
Mandelli et al. (2016), nonfluent/agrammatic patients (N=34), T1 baseline and T2 one year later
*/

/* RohrerMandelli T1 data on PPA: Rohrer et al. (2013), Mandelli et al. (2016) */
double REAL_DATA_ROHRERMANDELLI_T1[N_GROUPs][N_TASKs] = {
	                         /* Naming  Comprehension Repetition */
	/* Control              */ { 90.3,      96.3,      96.7 }, /* dummy, from Savage */
	/* Nonfluent/agrammatic */ { 76.7,      99.0,      81.5 },
	/* Semantic	            */ { 26.7,      88.0,      90.6 }, /* dummy, from Brambati */
	/* Logopenic            */ { 61.0,      94.0,      94.0 }
};

/* RohrerMandelli T2 data on PPA: Rohrer et al. (2013), Mandelli et al. (2016) */
double REAL_DATA_ROHRERMANDELLI_T2[N_GROUPs][N_TASKs] = {
	                         /* Naming  Comprehension Repetition */
	/* Control              */ { 90.3,      96.3,      96.7 }, /* dummy, from Savage */
	/* Nonfluent/agrammatic */ { 66.0,      90.0,      65.5 },
	/* Semantic	            */ { 26.7,      88.0,      90.6 }, /* dummy, from Brambati */
	/* Logopenic            */ { 43.0,      85.0,      77.0 }
};



double REAL_DATA[N_GROUPs][N_TASKs];
double SIM_DATA[N_GROUPs][N_TASKs];
double GOODNESS_OF_FIT[N_lesion_values];


double WEIGHT_value[N_lesion_values];
double DECAY_value[N_lesion_values];


/* concept and lemma */
 double C_node_act[N_CONCEPTs], L_node_act[N_LEMMAs];
 /* output form */
 double M_node_act[N_MORPHEMEs], oP_node_act[N_PHONEMEs], S_node_act[N_SYLLABLEs];
 /* input form */
 double iM_node_act[N_MORPHEMEs], iP_node_act[N_PHONEMEs];


 /* input buffer for nodes */
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
 double SEM_rate = 0.2 * 0.0101 * STEP_SIZE;   /* prop per step_size ms */
                                    /* adjusted for larger network: 0.2 */
 double LEM_rate = 0.0074 * STEP_SIZE;   /* prop per step_size ms */
 double LEX_rate = 0.0120 * STEP_SIZE;   /* prop per step_size ms */
 double DECAY_rate = 0.0240 * STEP_SIZE; /* prop per step_size ms */
 double EXTIN = 0.1965 * STEP_SIZE;      /* act_units per step_size ms */
 double LEMLEXFRAC = 0.3; 
  /* fraction of LEX_rate spread between lemmas and output morphemes */
  /* implementing weak cascading of activation, see Roelofs (2008, JEP:LMC) */ 

 double FR = 0.10;  /* fraction of connection weight for input phoneme 
					   to input morpheme, cf. Roelofs (1997, Cognition) */
 int    SEGMENT_DURATION = 125;  /* ms */
 int    PICTURE_DURATION = 125;  /* ms */

 /* set here to simulate weight or decay lesion and what to print */
 int WEIGHT_LESION = 1;
 int DECAY_LESION = 0;

 int SHOW_RESULTS_ALL_VALUES = 0; /* set here whether to print all values */
 

/* Aphasia parameters */

/* weight lesion */

 double CONNECTION_DECREASE_NONFLUENT_AGRAMMATIC;
 /* connections to and from output phonemes */

 double CONNECTION_DECREASE_SEMANTIC_DEMENTIA;
 /* connections to, within, and from conceptual network*/
 
 double CONNECTION_DECREASE_LOGOPENIC;
 /* connections to and from lexical output forms, 
    and between input and output phonemes */

/* decay lesion */

 double DECAY_INCREASE_NONFLUENT_AGRAMMATIC;
 /* output phonemes */

 double DECAY_INCREASE_SEMANTIC_DEMENTIA;
 /* concepts */

 double DECAY_INCREASE_LOGOPENIC;
 /* lexical output forms */



double ACT_C[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs];
double ACT_S[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs];

/* Activation of target concept, cat */
double ACT_CT[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation of conceptual relative, dog */
double ACT_CR[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs];

/* Activation of target lemma, cat */
double ACT_LT[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation of lemma relative, i.e., semantically related, dog */
double ACT_LR[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs];

/* Activation of target syllable, cat */
double ACT_ST[N_lesion_values][N_STEPs][N_GROUPs][N_TASKs]; 
/* Activation of syllabic relative, mat */
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

	double ls; /* exact lesion value */

    print_heading();

	print_parameters();

	set_spreading_rates();

	
	if (WEIGHT_LESION)
	for (lesion_value = 0, ls = 0.0; lesion_value < N_lesion_values; lesion_value++, ls += 0.01)
		WEIGHT_value[lesion_value] = ls; /* values between maximally damaged, 0.0, 
										    and minimally damaged, 0.99 */

	if (DECAY_LESION)
 	for (lesion_value = 0, ls = 1.01; lesion_value < N_lesion_values; lesion_value++, ls += 0.01)
		DECAY_value[lesion_value] = ls; /* values between minimally damaged, 1.01, 
										   and maximally damaged, i.e., full decay, 1.66 */


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
			if (assessment == ENGLISH)
				REAL_DATA[group][task] = REAL_DATA_ENGLISH[group][task];
			else if (assessment == DUTCH)
				REAL_DATA[group][task] = REAL_DATA_DUTCH[group][task];
			else if (assessment == BRAMBATI_T1)
				REAL_DATA[group][task] = REAL_DATA_BRAMBATI_T1[group][task];
			else if (assessment == BRAMBATI_T2)
				REAL_DATA[group][task] = REAL_DATA_BRAMBATI_T2[group][task];
			else if (assessment == ROHRERMANDELLI_T1)
				REAL_DATA[group][task] = REAL_DATA_ROHRERMANDELLI_T1[group][task];
			else if (assessment == ROHRERMANDELLI_T2)
				REAL_DATA[group][task] = REAL_DATA_ROHRERMANDELLI_T2[group][task];
	
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
	CC_con[i][j]*=SEM_rate;

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

  if (group == NONFLUENT_AGRAMMATIC)
	  CONNECTION_DECREASE_NONFLUENT_AGRAMMATIC = WEIGHT_FACTOR;
  else
	  CONNECTION_DECREASE_NONFLUENT_AGRAMMATIC = 1.0; /* normal */

  if(group == SEMANTIC_DEMENTIA) 
     CONNECTION_DECREASE_SEMANTIC_DEMENTIA = WEIGHT_FACTOR; 
  else 
     CONNECTION_DECREASE_SEMANTIC_DEMENTIA = 1.0; /* normal */
  
  if (group == LOGOPENIC) 
	  CONNECTION_DECREASE_LOGOPENIC = WEIGHT_FACTOR;
  else 
	  CONNECTION_DECREASE_LOGOPENIC = 1.0; /* normal */


	/* setting of decay parameters */

  if (group == NONFLUENT_AGRAMMATIC)
	  DECAY_INCREASE_NONFLUENT_AGRAMMATIC = DECAY_FACTOR;
  else
	  DECAY_INCREASE_NONFLUENT_AGRAMMATIC = 1.0; /* normal */

  if (group == SEMANTIC_DEMENTIA)
	  DECAY_INCREASE_SEMANTIC_DEMENTIA = DECAY_FACTOR;
  else
	  DECAY_INCREASE_SEMANTIC_DEMENTIA = 1.0; /* normal */

  if (group == LOGOPENIC)
	  DECAY_INCREASE_LOGOPENIC = DECAY_FACTOR;
  else
	  DECAY_INCREASE_LOGOPENIC = 1.0; /* normal */

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
		   input_C[CAT] += CONNECTION_DECREASE_SEMANTIC_DEMENTIA * EXTIN;

    /* enhancement */
  	  if( (T >= (0 + 1 * CYCLE_TIME)) 
		         &&  (T < (1 * CYCLE_TIME + PICTURE_DURATION))  ) 
	     input_C[CAT]+= EXTIN;
	}
   


    if( task == COMPREHENSION || task == REPETITION) {
 
    /* spoken word input */
      if((0 <= T) && (T < SEGMENT_DURATION)) {
	         input_iP[pK]+=EXTIN;
             }

      if((SEGMENT_DURATION <= T) 
          && (T < (2 * SEGMENT_DURATION))) {
	         input_iP[pA]+=EXTIN;
             }


     if((2 * SEGMENT_DURATION <= T) 
         && (T < (3 * SEGMENT_DURATION))) {
	         input_iP[pT]+=EXTIN;
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
	                * (CC_con[j][i] * CONNECTION_DECREASE_SEMANTIC_DEMENTIA ) );

  for(i=0;i<N_CONCEPTs;i++)
     for(j=0;j<N_LEMMAs;j++)
       input_C[i]+=(L_node_act[j] 
	              * CL_con[j][i] * CONNECTION_DECREASE_SEMANTIC_DEMENTIA );


 /* input activation for lemma nodes */
  for(i=0;i<N_LEMMAs;i++)
	  for(j=0;j<N_CONCEPTs;j++) 
        input_L[i]+=( C_node_act[j] * CL_con[j][i] 
	                  *  CONNECTION_DECREASE_SEMANTIC_DEMENTIA );

  for(i=0;i<N_LEMMAs;i++)
     for(j=0;j<N_MORPHEMEs;j++)                      
        input_L[i]+=( iM_node_act[j] * iML_con[j][i]);


 /* input activation for output morpheme nodes */
   for(i=0;i<N_MORPHEMEs;i++)
     for(j=0;j<N_LEMMAs;j++) 
	    input_M[i]+=( L_node_act[j] * LEMLEXFRAC * LM_con[j][i] * CONNECTION_DECREASE_LOGOPENIC); 


  for(i=0;i<N_MORPHEMEs;i++)
     for(j=0;j<N_MORPHEMEs;j++)
        input_M[i]+=( iM_node_act[j] * iMM_con[j][i] * CONNECTION_DECREASE_LOGOPENIC);



 /* input activation for output phoneme nodes */
  for(i=0;i<N_PHONEMEs;i++)
     for(j=0;j<N_MORPHEMEs;j++) 
	  input_oP[i]+=( M_node_act[j] * MP_con[j][i] 
		  * CONNECTION_DECREASE_NONFLUENT_AGRAMMATIC * CONNECTION_DECREASE_LOGOPENIC);

 for (i = 0; i < N_PHONEMEs; i++)
	  for (j = 0; j < N_PHONEMEs; j++)
		input_oP[i] += (iP_node_act[j] * PP_con[j][i] * CONNECTION_DECREASE_NONFLUENT_AGRAMMATIC 
			                  * CONNECTION_DECREASE_LOGOPENIC);


 /* input activation for syllable program nodes */
   for(i=0;i<N_SYLLABLEs;i++)
     for(j=0;j<N_PHONEMEs;j++)  
	input_S[i]+=( oP_node_act[j] * PS_con[j][i] * CONNECTION_DECREASE_NONFLUENT_AGRAMMATIC);



 /* input activation for input phoneme nodes */
  for(i=0;i<N_PHONEMEs;i++)
     for(j=0;j<N_PHONEMEs;j++) 
	input_iP[i]+=( oP_node_act[j] * PP_con[j][i] * CONNECTION_DECREASE_NONFLUENT_AGRAMMATIC 
		              * CONNECTION_DECREASE_LOGOPENIC);
 
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
                     * (1.0 - (DECAY_rate * DECAY_INCREASE_SEMANTIC_DEMENTIA)) ) + input_C[i]);

   for(i=0;i<N_LEMMAs;i++) 
     L_node_act[i]=((L_node_act[i] * (1.0 - DECAY_rate)) + input_L[i]);

    
   for(i=0;i<N_MORPHEMEs;i++)
       M_node_act[i]=((M_node_act[i] * (1.0 - (DECAY_rate * DECAY_INCREASE_LOGOPENIC) )) + input_M[i]);

     for(i=0;i<N_PHONEMEs;i++) 
       oP_node_act[i]=((oP_node_act[i] * (1.0 - (DECAY_rate * DECAY_INCREASE_NONFLUENT_AGRAMMATIC))) + input_oP[i]);
       
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
	ACT_ST[lesion_value][step][group][task] = S_node_act[sKAT];
	ACT_SR[lesion_value][step][group][task] = S_node_act[sMAT];

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
	 printf("Simulation of group studies \n");
	 printf("Testing the effect of a larger network: Large, all SYDBAT animals \n");

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

	 if (assessment == ENGLISH)
		 printf("\nAssessment is Savage et al. (2013), English\n");
	 if (assessment == DUTCH)
		 printf("\nAssessment is Janssen et al. (2022), Dutch\n");
	 if (assessment == BRAMBATI_T1)
		 printf("\nAssessment is Brambati et al. (2015), baseline T1\n");
	 if (assessment == BRAMBATI_T2)
		 printf("\nAssessment is Brambati et al. (2015), follow up T2\n");
	 if (assessment == ROHRERMANDELLI_T1)
		 printf("\nAssessment is Rohrer et al. (2013) and Mandelli et al. (2016), baseline T1\n");
	 if (assessment == ROHRERMANDELLI_T2)
		 printf("\nAssessment is Rohrer et al. (2013) and Mandelli et al. (2016), follow up T2\n");


   for(group=0; group <  N_GROUPs; group++) {
        printf(" \n");

        if(group == NORMAL) 
            printf("NORMAL \n");
		else if (group == NONFLUENT_AGRAMMATIC)
			printf("NONFLUENT/AGRAMMATIC \n");
		else if (group == SEMANTIC_DEMENTIA)
			printf("SEMANTIC DEMENTIA  \n");
		else if (group == LOGOPENIC)
			printf("LOGOPENIC  \n");

		if (assessment == ENGLISH)
			printf("\nSavage et al. (2013), English\n");
		if (assessment == DUTCH)
			printf("\nAssessment is Janssen et al. (2022), Dutch\n");
		if (assessment == BRAMBATI_T1)
			printf("\nAssessment is Brambati et al. (2015), baseline T1\n");
		if (assessment == BRAMBATI_T2)
			printf("\nAssessment is Brambati et al. (2015), follow up T2\n");
		if (assessment == ROHRERMANDELLI_T1)
			printf("\nAssessment is Rohrer et al. (2013) and Mandelli et al. (2016), baseline T1\n");
		if (assessment == ROHRERMANDELLI_T2)
			printf("\nAssessment is Rohrer et al. (2013) and Mandelli et al. (2016), follow up T2\n");

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
				+ fabs(REAL_DATA[group][REPETITION] - SIM_DATA[group][REPETITION]) ) / 3.0;

			if(SHOW_RESULTS_ALL_VALUES) /* toggle for printing the results for all lesion values */
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




 