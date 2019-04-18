////
/////////////////////////////////////////////////////////
// Spatial Operating Model built by Daniel Goethel (NMFS SEFSC)  
// edited by Katelyn Bosley and others
//////////////////////////////////////////////////////////

//need to add all_natal calcs for all values that might want to compare across an area, because regular non-natal
//homing calcs do not account for different weights across natal populations so if if any part of calc involves weight
//it needs to use the biomass_all_natal value not just biomass_AM

GLOBALS_SECTION
  #include "admodel.h"
  #include "qfclib.h"
  #define EOUT(var) cout <<#var<<" "<<var<<endl;

TOP_OF_MAIN_SECTION
  arrmblsize=500000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);

DATA_SECTION
////////////////////////////////////////////////////////////////////////////////////
/////MODEL STRUCTURE INPUTS/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
 init_number OM_structure //
  // ==0 the OM is panmictic
  // ==1 the OM is metamictic
  // ==2 the OM is metapop
  // ==3 the OM is natal homing - not there yet...

  init_int nages //number of ages
  init_int nyrs //number of years for simulation
  init_int npops //number of populations
///////////// Indices for ragged arrays (can't be type init so need to establish as integer)
  !! int np=npops;
  !! int ny=nyrs;
  !! int na=nages;
//////////////////////////////////////////////////////
  init_ivector nregions(1,np) //number of regions within a population - for metamictic regions = areas, populations = 1
  !! ivector nreg=nregions;
  init_ivector nfleets(1,np) //number of fleets in each region by each population
  !! ivector nf=nfleets;
  init_ivector nfleets_survey(1,np) //number of fleets in each region by each population
  !! ivector nfs=nfleets_survey;

 // INDICES FOR RNG...use these to input the MAX values across ALL runs/scenarios for a simulation experiment so that always have consistent RNG values
  init_int max_pops
  init_int max_regs
  init_int max_ages
  init_int max_yrs
  init_int max_flts
  init_int max_surv_flts
  init_int max_tag_yrs
////////////////////////////////////////////////////////////////////////////////////
//////////////SWITCHES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
  init_number model_type_switch
  ///// Changes the type of harvest model
  //==0 use to simulate catch based on input F (use if searching for MSY)
  //==1 use input TAC to set F (use if want to input a desired harvest level and use NR to determine associated F)
  //==2 use input harvest rate
  init_number do_tag
  //==0 do not calculate tagging data
  //==1 calculate tagging data
 //##############################################################################
 //##############################################################################
 //### used when model type>0 ###################################################
 //###########################################################################
  init_number parse_TAC
  //==0 do not alter the input TAC or harvest rate
  //==1 use observed data source to parse TAC or harvest rate (used when allocating harvest but pop structure unknown)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// USING calc_TAC_from_uMSY==1 IS PROBABLY ONLY CORRECT WAY TO PARSE u /////////////////////////////////////////////////
/////////////// IF PARSE input_u directly, then sum(u) unlikely to equal input_u because of differences in population sizes/////////
/////////////// ie applying less than the full uMSY to each area/region ////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_number calc_TAC_from_uMSY
  //MUST HAVE MODEL_TYPE==1, even though uses input u (harvest rate)
  //==0 don't use
  //==1 calculate a yearly TAC from uMSY(input)*biomass_population(j,y) to get a yearly TAC that obtains MSY in equil without crashing the stock
  //==2 calculate a yearly TAC from uMSY(input)*biomass_region(j,r,y) to get a yearly TAC that obtains MSY in equil without crashing the stock
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_number parse_TAC_source
  //==0 recruitment index_BM, assume index occurs at tspawn so always have 1 year timelag in using rec index to parse TAC
  //==1 recruitment index_AM, assume index occurs at tspawn so always have 1 year timelag in using rec index to parse TAC
  //==2 survey biomass, if tsurvey==0 can use survey in current year or assume timelag
  //==3 equal apportionment across all fleets in a given area
  init_number TAC_survey_parse_timelag_switch
  //==0 no timelag, use survey apportionment in current year (if tsurvey==0) or previous year (if tsurvey>0)
  //==1 use timelag, use survey apportionment from y-TAC_survey_parse_timelag, assume equal apportionment of TAC among fleets in first year
  init_number TAC_survey_parse_timelag
  //whole number value to implement a time lag in year that use survey apportionment from
  //##############################################################################
  //##############################################################################
 
  init_matrix tsurvey(1,np,1,nreg) //time of survey in proportion of year (0,1)
  init_number larval_move_switch
  ///// Changes the type of larval movement pattern (sets age class 1 movements)
  //==(-1) larval (age-1) movement is determined in same way as other ages (e.g., based on move_switch) 
  //==0 no movement
  //==1 input movement
  //==2 movement within population only based on residency (symmetric)
  //==3 symmetric movement but only allow movement within a population (ie regions within a population) not across populations
  //==4 symmetric movement across all populations and regions
  //==5 allow movement across all regions and populations, based on population/region specific residency (symmetric off-diag)

  init_number move_switch
  ///// Sets the type of adult movement pattern (sets age class>1 movements)
  //==0 no movement
  //==1 input movement
  //==2 movement within population only based on input residency (symmetric off diagnol)
  //==3 symmetric movement but only allow movement within a population (ie regions within a population) not across populations
  //==4 symmetric movement across all populations and regions
  //==5 allow movement across all regions and populations, based on population/region specific residency (symmetric off diagnol)
  //==6 natal return, based on age of return and return probability (certain fraction of fish make return migration to natal population eg ontogenetic migration)
  //==7 larvae stay in the population that they move to (i.e., for overlap, do not return to natal population if adult movement==0...otherwise with natal
  //    homing would return to natal population because natal residency is 100% and use natal movement rates (not current population movement rates like with metapopulation/random movement))
  //==8 density dependent movement based on relative biomass among potential destination area/regions, partitions (1-input_residency) based on a logistic function of biomass in current area/region and 'suitability' of destination area/regions
  //// uses use_input_Bstar switch
  //// DD MOVEMENT CAN BE AGE BASED OR CONSTANT ACROSS AGES...FOR AGE BASED MAKE SURE DD_move_age_switch==1, FOR AGE-INVARIANT DD_move_age_switch==0
  //==9 use input T_year to allow T to vary by year
  //==10 use T_FUL_Input to input time and age varying T

  init_number rand_move
  ///adjust DD movement to include random variation based on input sigma and lognormal random variable
  //IF move_switch==0 or larval_move_switch==0, then NO RANDOMNESS IS EMPLOYED (I.E., RESIDENCY STAYS AT 100%)
  //==0 no randomness, just use movement from T calcs
  //==1 add randomness to T (bounded so Fract move cannot exceed 1 or go below 0)
  
  init_number DD_move_age_switch
  /////// Allow age-based movement when using DD movement (Y/N) (1/0)
  //#==0 no age-based DD movement (assumes movement based on total relative biomass after movement, before mortality, in previous year)   ///can't use current year because biomass matrix not filled out at time of DD calcs, bec need T for calc of abund_AM which occurs before moving to next age loop
  //#==1 DD movement is age-based (assumes movement based on age-specfic relative biomass at beginning of current year)

  init_number use_input_Bstar
  //==0 set Bstar for DD movement equal to SSB_zero*SSB_zero_appor (if nreg==1, Bstar=SSB0), **NOTE** for Rec_type!=2 (not BH), defaults to using input_Bstar since no SSB0 calcs are done
  //==1 use input_Bstar

////// Population Structure switches
  init_number natal_homing_switch
  //==0 no natal homing (SSB is sum of SSB in population regardless of natal origin; weight/mat/fecund/ are based on current population not natal population) - Metapopulation/metamictic
  //==1 do natal homing (a fish only adds to SSB if it is in its natal population at spawning time; weight/mat/fecund/ are based on natal population) - Natal homing
  //natal homing  assumes genetic based life history and contribution to SSB (i.e., natal homing and no demographic mixing), natal homing==0 assumes demographic mixing (e.g. metapopulations where life history is more location based)

  init_number spawn_return_switch
   //==0 if natal_homing_switch==1 then only fish that are in natal population add to SSB
   //==1 natal_homing_switch==1 a fraction of fish return to natal population to spawn (inpsantaneous migration to natal population and back at time of spawning) based spawn_return_prob; weight/mat/fecund/ are based on natal population)
//////////////////////////////////////////////////////

  init_number select_switch
  //==0 input selectivity
  //==1 logistic selectivity based on input sel_beta1 and sel_beta2
  //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4
/////////////////////////////////////////////////////

  init_number select_switch_survey
  //==0 input selectivity
  //==1 logistic selectivity based on input sel_beta1 and sel_beta2
  //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4
/////////////////////////////////////////////////////

  init_number F_switch
  //==1 input F
  //==2 input single overall F (FMSY)
  //==3 Estimate F that minimizes difference in input and estimated total harvest rate
  //==4 overall F (FMSY) is split evenly among populations (each fleet uses population F)
  //==5 overall F (FMSY) is is split evenly among all regions (each fleet uses region F)
  //==6 overall F (FMSY) is split evenly among fleets
  //==7 F devs about input F based on sigma_F
  //==8 random walk in F
  
  init_number recruit_devs_switch
  //==0 use stock-recruit relationphip directly
  //==1 allow lognormal error around SR curve (i.e., include randomness based on input sigma_recruit)

  init_number recruit_randwalk_switch
  //==0 no random walk recruitment deviations
  //==1 have random walk lognormal recruitment deviations (requirs recruit_devs_switch==1)....NEEDS WORK!!!

  init_number init_abund_switch
  //==0 input init abundance
  //==1 decay from R_ave

 //determine how to estimate R0 when there are multiple regions within a population that have different vital rates
  init_number maturity_switch_equil
  //==0 for equal by area or average
  //==1 weighted average using equil_ssb_apportion to determine proportional contribution to equil vital rates by region
  //SSB0 must be calculated to determine stock-recruit function (if only know steepness and R0 for the population)
  //Use equilibrium SPR calcs to get SSB0, but to do so requires vital rates (maturity, weight), which are typically constant across a population
  //With multiple regions within a pop each with different vitals, must make assumption regarding the proportional contribution of each region's demograhics to equil SSB
  //When ==1 just assume equal (average) contributions, when ==1 input a proportional contribution (ie assume one region has higher carrying capacity and contributes more to equil SSB)
  
  init_number SSB_type
  //==1 fecundity based SSB
  //==2 weight based SSB
  
  init_number apportionment_type
  //==-1 no recruitment apportionment to regions within a population (each region within a population gets full amount of recruits from SR curve)
  //==0 apportionment to each region is based on relative SSB in region compared to population SSB
  //==1 input apportionment
  //==2 recruits are apportioned equally to each region within a population
  //==3 recruits are apportioned in a completely random manner with uniform equilibrium distribution
  //==4 recruits are apportioned stochastically with normal error surrounding the input proportions...uses the multivariate logistic random variables (Cox and Krunland 2008, FIsheries Research)
  //==5 recruits are approtioned based on theoretical enviormental phase shift.. working on
  
  init_number Rec_type
  //==1 stock-recruit relationship assumes an average value based on R_ave
  //==2 Beverton-Holt population-recruit functions based on population-specific input steepness, R0 (R_ave), M, and weight
  //==3 environmental recruitment - sine fucntion based on amplitude and frequency

  init_number use_stock_comp_info_survey
  //Determines whether it is assumed that info (stock composition data) is available determine natal origin for age composition data
  //==0 calc OBS survey age comps by area (summed across natal population)
  //==1 calc OBS survey age comps by natal population within each area

  init_number use_stock_comp_info_catch
  //Determines whether it is assumed that info (stock composition data) is available determine natal origin for age composition data
  //==0 calc OBS catch age comps by area (summed across natal population)
  //==1 calc OBS catch age comps by natal population within each area

  init_number number_tags_switch
   //determines how total tags are determined and how they are distributed across regions
   //==0 use input tags by year, pop, reg
   //==1 total tags are based on the input fract of total  abund tagged and is distributed evely across regions
   //==2 regional tags are based on the input fract of regional abund tagged
   //==3 total tags are based on the input of total abund tagged and is distributed based on survey prop of abundance in each region
   //==4 regional tags are randomly (uniformly) distributed based on max and min number of tags per region
   //==5 regional tags are randomly (uniformly) distributed based on max and min number of tags per region AND
   // prob that tagging occurs in a given year/region is based on opport_tag_prob and a uniform distributed random number
   
  init_number age_dist_tags
   //determines how tags are distributed across ages
   //==0 regional tags distributed evely across all ages
   //==1 regional tags distributed according to survey proportions at age
   //==2 regional tags distributed according to regional catch proportions at age
   
  init_number tag_fit_ages_switch
   //determines whether tags are fit by age cohorts or by region-only cohorts
   //the latter is for situations where age is unknown and so assume no tag age dynamics and that all fish fully selected
   //in this case OM assumes normal tagging age-based dynamics, but EM ignores age structure in tags and assumes all tagged fish are fully selected
   //therefore creates inherent process error in tag dynamics because EM uses T, F, M dynamics of fully selected age (error magnified if OM or EM has age-based M or T)
   //==0, fit by age-based cohorts
   //==1, fit by region-based cohorts

   init_number tag_fit_ages_switch_OM
   //  #DOES NOT WORK WITH TAG MIXING!!!
   
   //determines whether OM should use same dynamics as EM when fitting tags by cohort
   //if ==1 then OM tag dynamics use only fully selected F instead of maintaining age-based dynamics then summing across ages
   //==0, OM differs from EM (OM maintains age-based dynamics)
   //==1, OM matches EM and uses only fully selected F to calculate tag recaps (use tag_age_sel to define age of full selection)

  init_number sim_tag_mixing_switch
   //determines whether tags have a different F or T in first year of release compared to rest of population (i.e., if there is incomplete mixing)
   //==0 F and T same as rest of pop (complete mixing)
   //==1 F and/or T have are different from rest of population (incomplete mixing)
  init_number sim_tag_mixing_T_switch
   //determines whether tags have a differentT in first year of release compared to rest of population (i.e., if there is incomplete mixing)
   //==0 T same as rest of pop (complete mixing)
   //==1 T have are different from rest of population (incomplete mixing)
  
  init_number est_tag_mixing_switch
   //determines whether different F or T compared to rest of population in first year of release are estimated  (i.e., estimate F and/or T to account for  incomplete mixing)
   //==0 F and T same as rest of pop (assume complete mixing)
   //==1 F is estimated different from rest of population (incomplete mixing F only) in first year of release
   //==2 T is estimated different from rest of population (incomplete mixing T only) in first year of release
   //==3 F AND T are estimated different from rest of population (incomplete mixing F AND T) in first year of release

///////////////////////////////////////////////////////////////////////////////
//////// ADDITIONAL PARAMETERS FROM DAT FILE //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/////tagging data parameters
  init_int nyrs_release //number of years with tag release events
  !! int ny_rel=nyrs_release;
  init_vector yrs_releases(1,ny_rel) //vector containing the model years with releases
  init_3darray input_ntags(1,np,1,nreg,1,ny_rel)
  init_vector input_total_tags(1,ny_rel)
  init_vector frac_abund_tagged(1,ny_rel) //proportion of abundance that is tagged in each release year
  init_matrix max_tags(1,np,1,nreg) //maximum tags for uniform distribution if tags_switch==4 or 5
  init_matrix min_tags(1,np,1,nreg) //maximum tags for uniform distribution if tags_switch==4 or 5
  init_number opport_tag_prob //cutoff value defining whether uniform RNG represents no tagging(<tag_prob) or tagging (>tag_prob) in a given area
  init_number opport_tag_prob_year //cutoff value defining whether uniform RNG represents no tagging(<tag_prob) or tagging (>tag_prob) in a given year
  init_int max_life_tags //number of years that tag recaptures will be tallied for after release (assume proportional to longevity of the species)...use this to avoid calculating tag recaptures for all remaining model years after release since # recaptures are often extremely limited after a few years after release
  init_int tag_age_sel //age used to compute tag dynamics (F and T) when fitting by cohort and OM tag dynamics match EM
  
  init_number SIM_ntag //the ESS used to simulate multinomial tagging dat
  init_number myseed_ntags //just the seed for RNG
  init_number myseed_prob_tag //just the seed for RNG
  init_number myseed_prob_tag_year //just the seed for RNG
  init_3darray report_rate_TRUE(1,np,1,ny_rel,1,nreg) //tag reporting rate (assume constant for all recaptures within a given release cohort, but can be variable across populations or regions)...could switch to allow variation across fleets instead
//reporting rate is assumed to be function of release event and recap location (not a function of recap year...could expand to this, but not priority at momement)

  init_vector F_tag_scalar(1,ny_rel)
  init_vector T_tag_res(1,ny_rel)

///////////
  init_vector sigma_T(1,np) //variance term for random variation in T
///////Density-dependent movement parameters
  init_3darray input_Bstar(1,np,1,nreg,1,na) //used with move_switch==8
  init_matrix SSB_zero_appor(1,np,1,nreg) //used with move_switch==8
  init_3darray A(1,np,1,nreg,1,na) //used with move_switch==8
  init_3darray DD_residency(1,np,1,nreg,1,na) //used with move_switch==8
/////////////
  init_number return_age // used if move_swith ==6
  init_vector return_probability(1,np) // used if move_swith==6
  init_vector spawn_return_prob(1,np) // used if natal_homing_swith==2
  init_number phase_F //must be turned on (==1) if F_type==3
  init_number phase_dummy //must be turned on (==1) if F_type!=3
  init_vector tspawn(1,np) //time of spawning in proportion of year (0-1)
  init_vector steep(1,np) //B-H steepness
  init_vector ln_R_ave(1,np) //Average Recruitment or R0 for B-H S-R curve
  init_vector amplitude(1,np) //amplitude of periodic recruitment in % of R_ave 
  init_vector freq(1,np) //frequency of recruitment in years (ie 10 for peak every 10 years)
  init_5darray input_T(1,np,1,nreg,1,na,1,np,1,nreg)  //movement matrix
  init_5darray input_T_year(1,np,1,nreg,1,ny,1,np,1,nreg)  //movement matrix
  init_matrix input_residency_larval(1,np,1,nreg)  //larval residency probability
  init_3darray input_residency(1,np,1,nreg,1,na) //
  init_3darray sel_beta1(1,np,1,nreg,1,nf)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta2(1,np,1,nreg,1,nf)   //selectivity inflection parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta3(1,np,1,nreg,1,nf)  //selectivity slope parameter 2 for double selectivity
  init_3darray sel_beta4(1,np,1,nreg,1,nf)  //selectivity inflection parameter 2 for double logistic selectivity
  init_3darray sel_beta1_survey(1,np,1,nreg,1,nfs)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta2_survey(1,np,1,nreg,1,nfs)   //selectivity inflection parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta3_survey(1,np,1,nreg,1,nfs)  //selectivity slope parameter 2 for double selectivity
  init_3darray sel_beta4_survey(1,np,1,nreg,1,nfs)  //selectivity inflection parameter 2 for double logistic selectivity
  init_4darray input_selectivity(1,np,1,nreg,1,na,1,nf) //fishery selectivity by area/region/age/fleet
  init_4darray input_survey_selectivity(1,np,1,nreg,1,na,1,nfs)//survey selectivity
  init_3darray q_survey(1,np,1,nreg,1,nfs) // catchability for different surveys(fleets)operating in different areas
  init_3darray input_F(1,np,1,nreg,1,nf)
  init_3darray dunce_F(1,np,1,nreg,1,3) //min and max F for dunce cap F alternative
  init_3darray F_rho(1,np,1,nreg,1,nf) //degree of autocorrelation (0-1) if F switch = 8; random walk in F
  init_number input_F_MSY
  init_matrix input_M_TRUE(1,np,1,na)
  init_vector sigma_recruit(1,np)
  init_vector sigma_rec_prop(1,np) //error around recruit apportionment
  init_3darray sigma_F(1,np,1,nreg,1,nf)

//##########################################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################
//#########################################################################################################################################
//######### FOR NATAL HOMING THERE IS NO ACCOUNTING OF REGIONAL DIFFERENCES IN VITAL RATES ACROSS REGIONS WITHIN A POPULATION
//######### IE BECAUSE GENETICS DEFINE VITAL RATES, THEY MUST ALL BE THE SAME
//######### **********DO NOT INPUT REGIONALLY VARYING VITAL RATES, NATAL REGION WILL NOT BE PROPERLY TRACKED IN SSB CALCS #############
//#########################################################################################################################################
  init_3darray input_weight(1,np,1,nreg,1,na)
  init_3darray input_catch_weight(1,np,1,nreg,1,na)
  init_3darray fecundity(1,np,1,nreg,1,na)
  init_3darray maturity(1,np,1,nreg,1,na)
  init_matrix prop_fem(1,np,1,nreg) //proportion of population assumed to be female for SSB calcs (typically use 0.5)

//##########################################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################
//#########################################################################################################################################

  init_matrix input_Rec_prop(1,np,1,nreg)
  init_matrix equil_ssb_apport(1,np,1,nreg)

  init_3darray frac_natal(1,np,1,np,1,nreg)
  init_4darray input_init_abund(1,np,1,np,1,nreg,1,na)
  init_matrix rec_index_sigma(1,np,1,nreg)
  init_3darray sigma_survey(1,np,1,nreg,1,nfs)
  init_4darray sigma_survey_overlap(1,np,1,np,1,nreg,1,nfs)
  init_3darray sigma_catch(1,np,1,nreg,1,nf)
  init_4darray sigma_catch_overlap(1,np,1,np,1,nreg,1,nf)
  init_3darray SIM_ncatch(1,np,1,nreg,1,nf) //cannot exceed 2000, otherwise change dimension of temp vector below
  init_4darray SIM_ncatch_overlap(1,np,1,np,1,nreg,1,nf) //cannot exceed 2000, otherwise change dimension of temp vector below
  init_3darray SIM_nsurvey(1,np,1,nreg,1,nfs) //cannot exceed 2000, otherwise change dimension of temp vector below
  init_4darray SIM_nsurvey_overlap(1,np,1,np,1,nreg,1,nfs) //cannot exceed 2000, otherwise change dimension of temp vector below

//################################################################################################################
//################################################################################################################
//################################################################################################################
//### IF PARSE TAC OR U USING OBS DATA THEN MAKE SURE THAT FULL TAC OR U FOR THAT AREA IS INPUT FOR EACH FLEET IN THAT AREA ###
//########################################################################################################
  init_3darray input_TAC(1,np,1,nreg,1,nf)
  init_3darray input_u(1,np,1,nreg,1,nf)

//################################################################################################################
//################################################################################################################
//################################################################################################################

//NR parameters
  init_number max_Fnew //
  init_number Fnew_start
  init_number NR_iterations
  init_number NR_dev
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //#### EM inputs
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
  //population structure for formatting inputs for EM
  init_number EM_structure //
  // ==0 the RM is panmictic
  // ==1 the EM is metamictic
  // ==2 the EM is metapop  
  init_number npops_EM
  !! int np_em=npops_EM;

//////////////////////////////////////////////////////
  init_ivector nregions_EM(1,np_em) //number of regions within a population - for metamictic regions = areas, populations = 1
  !! ivector nreg_em=nregions_EM;
  init_ivector nfleets_EM(1,np_em) //number of fleets in each region by each population
  !! ivector nf_em=nfleets_EM;
  init_ivector nfleets_survey_EM(1,np_em) //number of fleets in each region by each population
  !! ivector nfs_em=nfleets_survey_EM;

  init_matrix tsurvey_EM(1,np_em,1,nreg_em)

  init_number diagnostics_switch
  //#==0 allow OBS data to be used for estimation
  //#==1 allow TRUE data from without error to be used for estimation

  init_number fleets_as_areas_switch // bypasses some challlenging portions of the code  - temporary
  //==0 not fleets as areas
  //==1 fleets as areas approach for EM

  init_number move_switch_EM
  ///// Sets the type of adult movement pattern (sets age class>1 movements)
  //==0 no movement, set T phases=-1
  //==1 input movement, set T phases=-1
  //>1 estimate movement among all regions/pops...set either of the movement phases to >0
  
  init_number report_rate_switch_EM
  //==(-1), use input_report_rate_EM, reporting rate may differ from TRUE report rate
  //==0, use TRUE reporting rate from OM
  //==1, estimate time-invariant (pop, reg varying) reporting rate, ph_rep_rate_CNST>0
  //==2, estimate a time-, pop-, reg-varying reporting rate ph_rep_rate_YR>0
  
////// Population Structure switches
  init_number natal_homing_switch_EM
  //==0 no natal homing (SSB is sum of SSB in population regardless of natal origin; weight/mat/fecund/ are based on current population not natal population) - Metapopulation/metamictic
  //==1 do natal homing (a fish only adds to SSB if it is in its natal population at spawning time; weight/mat/fecund/ are based on natal population) - Natal homing
  //natal homing  assumes genetic based life history and contribution to SSB (i.e., natal homing and no demographic mixing), natal homing==0 assumes demographic mixing (e.g. metapopulations where life history is more location based)

  init_number spawn_return_switch_EM
   //==0 if natal_homing_switch==1 then only fish that are in natal population add to SSB
   //==1 natal_homing_switch==1 a fraction of fish return to natal population to spawn (inpopsantaneous migration to natal population and back at time of spawning) based spawn_return_prob; weight/mat/fecund/ are based on natal population)

  init_number select_switch_EM
  //==0 input selectivity
  //==1 logistic selectivity based on input sel_beta1 and sel_beta2
  //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4

  init_number select_switch_survey_EM
  //==0 input selectivity
  //==1 logistic selectivity based on input sel_beta1 and sel_beta2
  //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4

 //determine how to estimate R0 when there are multiple regions within a population that have different vital rates
  init_number maturity_switch_equil_EM
  //==0 for equal by area or average
  //==1 weighted average using equil_ssb_apportion to determine proportional contribution to equil vital rates by region
  //SSB0 must be calculated to determine stock-recruit function (if only know steepness and R0 for the population)
  //Use equilibrium SPR calcs to get SSB0, but to do so requires vital rates (maturity, weight), which are typically constant across a population
  //With multiple regions within a pop each with different vitals, must make assumption regarding the proportional contribution of each region's demograhics to equil SSB
  //When ==1 just assume equal (average) contributions, when ==1 input a proportional contribution (ie assume one region has higher carrying capacity and contributes more to equil SSB)
  
  init_number SSB_type_EM
  //==1 fecundity based SSB
  //==2 weight based SSB
  
  init_number Rec_type_EM
  //==1 stock-recruit relationship assumes an average value based on R_ave
  //==2 Beverton-Holt population-recruit functions based on population-specific input steepness, R0 (R_ave), M, and weight

  init_number apportionment_type_EM
  //==-1 no recruitment apportionment to regions within a population (each region within a population gets full amount of recruits from SR curve)
  //==0 apportionment to each region is based on relative SSB in region compared to population SSB
  //==1 input apportionment
  //==2 recruits are apportioned equally to each region within a population
  //==3 estimate pop/reg apportionment
  //==4 estimate pop/reg/year apportionment

  init_number use_stock_comp_info_survey_EM //for likelihood calcs
  //Determines whether it is assumed that info (stock composition data) is available determine natal origin for age composition data
  //==0 calc survey age comps by area (summed across natal population)
  //==1 calc survey age comps by natal population within each area

  init_number use_stock_comp_info_catch_EM //for likelihood calcss
  //Determines whether it is assumed that info (stock composition data) is available determine natal origin for age composition data
  //==0 calc  catch age comps by area (summed across natal population)
  //==1 calc  catch age comps by natal population within each area
   init_number F_switch_EM
  //==1 estimate yearly F
  //==2 random walk in F
   init_number M_switch_EM
   //all M types are constant across regions
  //==(-1) use input_M (can differ from TRUE M)
  //==0 use input_M_TRUE
  //==1 estimate constant M (const across pop and age)
  //==2 estimate population-based M (const across ages)
  //==3 estimate age-based M (const across pop)
  //==4 estimate age- and population-varying M

  init_number recruit_devs_switch_EM
  //==0 use stock-recruit relationphip directly (make sure to set ph_rec=0), also assumes initial abund for all ages=R0
  //==1 allow lognormal error around SR curve (i.e., include randomness based on input sigma_recruit)

   init_number recruit_randwalk_switch_EM
  //==0 no random walk recruitment deviations
  //==1 have random walk lognormal recruitment deviations (requirs recruit_devs_switch==1)....NEEDS WORK!!!!!

  init_number init_abund_switch_EM
  //==0 est init_abund at age, where age-1 abund in year 1 is R_ave (to avoid overparametrization of recruitment)
  //==1 decay from R_ave. If phase is negative, input_init_abund is used.
  init_number est_dist_init_abund_EM
  //==(-2) use input_dist_init_abund specified for EM (can differ from true)
  //==(-1) assume all fish in a pop are equally distributed across regions in that pop (no fish start outside natal pop)
  //==0 use true distribution of init_abundance
  //==1 estimate the spatial distribution of init abundance, need to make ph_non_natal_init OR ph_reg_init non-negative depending on spatial structure used

  init_vector tspawn_EM(1,np_em) //timing of spawn for EM needs to match npops_EM
  init_vector return_probability_EM(1,np_em)
  init_vector spawn_return_prob_EM(1,np_em)

  init_int do_tag_EM
  init_int do_tag_mult //if==0 assume neg binomial, if==1 assume multinomial (same as OM)

  init_5darray input_T_EM(1,np_em,1,nreg_em,1,na,1,np_em,1,nreg_em)// input T matrix for EM
  init_vector input_R_ave_EM(1,np_em)
  init_3darray input_rec_prop_EM(1,np_em,1,nreg_em,1,nyrs)//input recruit apportionment for EM
  init_vector sigma_recruit_EM(1,np_em)
  init_vector input_steep_EM(1,np_em)
  init_4darray init_abund_EM(1,np_em,1,np_em,1,nreg_em,1,na)
  init_3darray input_dist_init_abund_EM(1,np_em,1,np_em,1,nreg_em)
  init_matrix input_M_EM(1,np_em,1,na)

  init_4darray input_selectivity_EM(1,np,1,nreg,1,na,1,nf) //fishery selectivity by area/region/age/fleet
  init_4darray input_survey_selectivity_EM(1,np,1,nreg,1,na,1,nfs)//survey selectivity
  init_matrix input_report_rate_EM(1,np_em,1,nreg_em)

  init_int ph_lmr
  init_number lb_R_ave
  init_number ub_R_ave
  init_int ph_rec
  init_number lb_rec_devs
  init_number ub_rec_devs
  init_int ph_rec_app_CNST
  init_int ph_rec_app_YR
  init_number lb_rec_app
  init_number ub_rec_app
  init_int ph_init_abund
  init_int ph_reg_init
  init_int ph_non_natal_init
  init_number lb_init_dist
  init_number ub_init_dist
  init_number lb_init_abund
  init_number ub_init_abund
  init_int ph_F
  init_number lb_F
  init_number ub_F 
  init_int ph_steep
  init_number lb_steep
  init_number ub_steep
  init_int ph_M_CNST
  init_int ph_M_pop_CNST
  init_int ph_M_age_CNST
  init_int ph_M_pop_age
  init_number lb_M
  init_number ub_M
  init_int ph_sel_log
  init_number lb_sel_beta1
  init_number ub_sel_beta1
  init_number lb_sel_beta2
  init_number ub_sel_beta2
  init_number lb_sel_beta3
  init_number ub_sel_beta3
  init_number lb_sel_beta4
  init_number ub_sel_beta4
  init_number lb_sel_beta1_surv
  init_number ub_sel_beta1_surv
  init_number lb_sel_beta2_surv
  init_number ub_sel_beta2_surv
  init_number lb_sel_beta3_surv
  init_number ub_sel_beta3_surv
  init_number lb_sel_beta4_surv
  init_number ub_sel_beta4_surv
  init_int ph_sel_log_surv
  init_int ph_sel_dubl
  init_int ph_sel_dubl_surv
  init_int ph_q
  init_number lb_q
  init_number ub_q
  init_int ph_F_rho // if we want random walk F, not implemented
  init_number lb_F_rho
  init_number ub_F_rho
 // MOVEMENT PARAMETERS ///////////
  init_int phase_T_YR
  init_int phase_T_YR_ALT_FREQ
  init_int T_est_freq
  init_int phase_T_YR_AGE_ALT_FREQ
  init_int T_est_age_freq
  init_int juv_age
  init_int phase_T_CNST
  init_int phase_T_CNST_AGE
  init_int phase_T_YR_AGE
  init_int phase_T_CNST_AGE_no_AG1
  init_int phase_T_YR_AGE_no_AG1
  init_int phase_T_YR_AGE_ALT_FREQ_no_AG1
  init_number lb_T
  init_number ub_T
  init_number T_start  //starting value for T parameters in log space, for 2 pops ==-2 gives ~90% residency, for 3 pops ==-3 gives ~60% residency
//////////////////////////////////////

 // Tagging Parameters //////////////
  init_int phase_rep_rate_YR
  init_int phase_rep_rate_CNST
  init_number lb_B
  init_number ub_B
  init_number ph_T_tag
  init_number ph_F_tag
  init_number lb_scalar
  init_number ub_scalar
 ///////////////////////////

  init_int ph_dummy
 // Likelihood weights //////////
  init_number wt_surv
  init_number wt_catch
  init_number wt_fish_age
  init_number wt_srv_age 
  init_number wt_rec
  init_number wt_tag
  init_number wt_F_pen
  init_number wt_M_pen
  init_number wt_B_pen
 /////////////////////////////////////

 // Penalty Function Inputs /////////////////////////////
  init_number report_rate_sigma
  init_number report_rate_ave
  init_number abund_pen_switch // include penalty (norm2) on init_abund_devs?  0==no, 1==yes
  init_number wt_abund_pen
  init_number move_pen_switch //inlcude movement penalty in log space?  0==no, 1==yes
  init_number wt_move_pen
  init_number Tpen
  init_number Tpen2
  init_number sigma_Tpen_EM 
 /////////////////////////////////////////////////////////

 //error for the EM
  init_4darray OBS_survey_fleet_bio_se_EM(1,np_em,1,nreg_em,1,ny,1,nfs_em)
  init_4darray OBS_yield_fleet_se_EM(1,np_em,1,nreg_em,1,ny,1,nf_em)
  init_4darray OBS_survey_prop_N_EM(1,np_em,1,nreg_em,1,ny,1,nfs_em) //cannot exceed 2000, otherwise change dimension of temp vector below
  init_4darray OBS_catch_prop_N_EM(1,np_em,1,nreg_em,1,ny,1,nf_em) //cannot exceed 2000, otherwise change dimension of temp vector below
  init_4darray tag_N_EM(1,np_em,1,nreg_em,1,ny_rel,1,na)
  //summing tag prop array

  !! int xn=na*ny;
  init_5darray T_Full_Input(1,np,1,nreg,1,xn,1,np,1,nreg)
//###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################
 //###################################################################################################################################


//////////////////////////////////////////////////////////////
  init_int debug
  init_number myseed_yield
  init_number myseed_survey
  init_number myseed_F
  init_number myseed_rec_devs
  init_number myseed_rec_apport
  init_number myseed_rec_index
  init_number myseed_survey_age
  init_number myseed_catch_age
  init_number myseed_tag
  init_number myseed_T

  //fill in a vector of years
  vector years(1,nyrs)
  !!years.fill_seqadd(double(1),1.0);
  
    init_imatrix nregions_temp(1,np,1,np) //used to fill tag_recap matrices

  !! for(int j=1;j<=npops;j++) //recap stock
  !! {
  !!  for (int r=1;r<=npops;r++) //recap region
  !!  {
  !!    if(j<=r)
  !!     {
  !!     nregions_temp(j,r)=0;
  !!     }
  !!    if(j>r)
  !!     {
  !!     nregions_temp(j,r)=nreg(r); //create temp matrix that holds the number of regions that exist in all previous populations (so can sum for use in calcs below)
  !!     }
  !!   }
  !!  }
  ivector nreg_temp(1,np)

  int a
  int y
  int z
  int k
  int j
  int i
  int s
  int r
  int n
  int w
  int p
  int v
  int x
  int u
  int d
  int xx

//Add counters to enumerate the regions for the report out section for 6d arrays
  int region_counter
 !! cout << "debug = " << debug << endl;
 !! cout << "If debug != 1541 then .dat file not setup correctly" << endl;
 !! cout << "input read" << endl;

 
PARAMETER_SECTION

  !! ivector nr=nregions;
  !! int nps=npops;
  !! int nyr=nyrs;
  !! int nag=nages;
  !! ivector nfl=nfleets;
  !! ivector nfls=nfleets_survey;  


 init_matrix F_est(1,nps,1,nr,phase_F)

 //For dunce cap F
 matrix Fstartyr(1,nps,1,nr)
 matrix minF(1,nps,1,nr)
 matrix maxF(1,nps,1,nr)
 matrix stepF(1,nps,1,nr)
 vector R_ave(1,nps)
 
 // vitals
 6darray T(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr)
 5darray T_year(1,nps,1,nr,1,nyr,1,nps,1,nr)
 6darray rel_bio(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr)
 3darray Bstar(1,nps,1,nr,1,nag)
 3darray c(1,nps,1,nr,1,nag)
 4darray Fract_Move_DD(1,nps,1,nr,1,nyr,1,nag)
 5darray selectivity(1,nps,1,nr,1,nyr,1,nag,1,nfl)
 4darray selectivity_age(1,nps,1,nr,1,nag,1,nfl)
 4darray F_year(1,nps,1,nr,1,nyr,1,nfl)
 5darray F_fleet(1,nps,1,nr,1,nyr,1,nag,1,nfl)
 4darray F(1,nps,1,nr,1,nyr,1,nag)
 4darray M(1,nps,1,nr,1,nyr,1,nag)
 matrix rec_devs(1,nps,1,nyr)
 matrix rec_devs_randwalk(1,nps,1,nyr)
 4darray init_abund(1,nps,1,nps,1,nr,1,nag)
 4darray weight_population(1,nps,1,nr,1,nyr,1,nag)
 4darray weight_catch(1,nps,1,nr,1,nyr,1,nag)
 3darray wt_mat_mult(1,nps,1,nyr,1,nag)
 4darray wt_mat_mult_reg(1,nps,1,nr,1,nyr,1,nag)
 3darray ave_mat_temp(1,nps,1,nag,1,nr) //to calc average maturity
 matrix ave_mat(1,nps,1,nag) //to calc average maturity
 matrix SPR_N(1,nps,1,nag)
 matrix SPR_SSB(1,nps,1,nag)
 vector SPR(1,nps)
 vector SSB_zero(1,nps)
 vector alpha(1,nps)
 vector beta(1,nps)

//recruitment 
 3darray recruits_BM(1,nps,1,nr,1,nyr)
 3darray recruits_AM(1,nps,1,nr,1,nyr)
 3darray Rec_Prop(1,nps,1,nr,1,nyr-1)
 3darray Rec_prop_temp1(1,nps,1,nyr-1,1,nr)
 3darray Rec_prop_temp2(1,nps,1,nyr-1,1,nr)

 3darray rec_index_BM(1,nps,1,nr,1,nyr)
 3darray rec_index_prop_BM(1,nps,1,nr,1,nyr)
 3darray rec_index_BM_temp(1,nps,1,nyr,1,nr)
 3darray rec_index_AM(1,nps,1,nr,1,nyr)
 3darray rec_index_prop_AM(1,nps,1,nr,1,nyr)
 3darray rec_index_AM_temp(1,nps,1,nyr,1,nr)

 vector env_rec(1,nyr)

//abundance 
 4darray abundance_at_age_BM(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_at_age_AM(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_in(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_res(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_leave(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_spawn(1,nps,1,nr,1,nyr,1,nag)

//biomass
 4darray biomass_BM_age(1,nps,1,nr,1,nyr,1,nag)
 4darray biomass_AM_age(1,nps,1,nr,1,nyr,1,nag)
 3darray biomass_BM(1,nps,1,nr,1,nyr)
 3darray biomass_AM(1,nps,1,nr,1,nyr)
 4darray bio_in(1,nps,1,nr,1,nyr,1,nag)
 4darray bio_res(1,nps,1,nr,1,nyr,1,nag)
 4darray bio_leave(1,nps,1,nr,1,nyr,1,nag)

 //tagging data
  !! int nyr_rel=nyrs_release;
  !! ivector xy(1,nyr_rel);
  !! ivector nt(1,nyr_rel);
  !! ivector nt2(1,nyr_rel);
  !! int nt3=max_life_tags*sum(nregions)+1;
  !! ivector tag_age(1,nyr_rel);

  !!  for(int x=1; x<=nyrs_release; x++)
  !!   {
  !!    xx=yrs_releases(x);
  !!    xy(x)=min(max_life_tags,nyrs-xx+1);
  !!    nt(x)=xy(x)*sum(nregions)+1;
  !!    nt2(x)=nt(x)-1;
  !!    tag_age(x)=xy(x);
  !!   }

 7darray tags_avail(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr)
 7darray recaps(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr)
 7darray tag_prop(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr)
 5darray tag_prop_final(1,nps,1,nr,1,nyr_rel,1,nag,1,nt)
 5darray SIM_tag_prop(1,nps,1,nr,1,nyr_rel,1,nag,1,nt)
 5darray OBS_tag_prop_final(1,nps,1,nr,1,nyr_rel,1,nag,1,nt)
 3darray total_recap_temp(1,nps,1,nr,1,tag_age)
 vector rand_tag_prop_temp2(1,nt3)  
 5darray tag_prop_temp2(1,nps,1,nr,1,nyr_rel,1,nag,1,nt2)
 4darray tag_prop_temp2_no_age(1,nps,1,nr,1,nyr_rel,1,nt2)
 5darray rand_tag_prop_temp(1,nps,1,nr,1,nyr_rel,1,nag,1,2000) //should make function of max(ncatch) but had issues making an index
 4darray rand_tag_prop_temp_no_age(1,nps,1,nr,1,nyr_rel,1,2000) //should make function of max(ncatch) but had issues making an index

 matrix tags_avail_temp(1,nps,1,nr)
 3darray tag_prop_temp(1,nps,1,max_life_tags,1,nr)

 vector ntags_total(1,nyr_rel)
 matrix ntags_total_temp(1,nps,1,nr)
 4darray ntags(1,nps,1,nr,1,nyr_rel,1,nag)
 3darray ntags_region(1,nps,1,nr,1,nyr_rel)
 3darray prob_tag(1,nps,1,nr,1,nyr_rel)
 vector prob_tag_year(1,nyr_rel)

 3darray prob_tag_RN(1,nps,1,nr,1,nyr_rel)
 3darray ntags_RN(1,nps,1,nr,1,nyr_rel)
 vector prob_tag_year_RN(1,nyr_rel)
 4darray total_rec(1,nps,1,nr,1,nyr_rel,1,nag)
 4darray not_rec(1,nps,1,nr,1,nyr_rel,1,nag)
 4darray tag_prop_not_rec(1,nps,1,nr,1,nyr_rel,1,nag)

  3darray total_rec_no_age(1,nps,1,nr,1,nyr_rel)
  3darray not_rec_no_age(1,nps,1,nr,1,nyr_rel)
  3darray ntags_no_age(1,nps,1,nr,1,nyr_rel)
  4darray tag_prop_final_no_age(1,nps,1,nr,1,nyr_rel,1,nt)
  6darray tag_prop_no_age(1,nps,1,nr,1,nyr_rel,1,tag_age,1,nps,1,nr)
  3darray tag_prop_not_rec_no_age(1,nps,1,nr,1,nyr_rel)
  7darray tag_recap_no_age_temp(1,nps,1,nr,1,nyr_rel,1,tag_age,1,nps,1,nr,1,nag) //recaps
  4darray SIM_tag_prop_no_age(1,nps,1,nr,1,nyr_rel,1,nt)
  4darray OBS_tag_prop_final_no_age(1,nps,1,nr,1,nyr_rel,1,nt)
  3darray age_full_selection(1,nps,1,nr,1,nyr)
  vector age_full_selection_temp(1,nag)
  4darray F_tag(1,nps,1,nr,1,nyr_rel,1,nag)
  6darray T_tag(1,nps,1,nr,1,nyr_rel,1,nag,1,nps,1,nr)
  
 //survey index
 5darray survey_selectivity(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray survey_selectivity_age(1,nps,1,nr,1,nag,1,nfls)
 5darray survey_selectivity_temp(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 6darray true_survey_fleet_overlap_age(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag)
 6darray survey_at_age_region_fleet_overlap_prop(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag)
 6darray SIM_survey_prop_overlap(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag)
 // 6darray OBS_survey_prop_overlap(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag)
 6darray OBS_survey_prop_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag) //**fixed KB now matches dims for EM

 6darray true_survey_fleet_overlap_age_bio(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag)
 5darray true_survey_fleet_bio_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls)
 4darray true_survey_region_bio_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray true_survey_population_bio_overlap(1,nps,1,nyr,1,nps)
 matrix true_survey_natal_bio_overlap(1,nyr,1,nps)
 vector true_survey_total_bio_overlap(1,nyr)
 5darray true_survey_fleet_age(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 
 // some new data stuff to do true survey abundances for tag releases
 5darray true_survey_fleet_age_temp(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray true_survey_region_abundance(1,nps,1,nyr,1,nr,1,nag)
 4darray true_survey_region_abundance_temp(1,nps,1,nyr,1,nag,1,nr) 
 3darray true_survey_population_abundance_temp(1,nyr,1,nag,1,nps)
 3darray true_survey_population_abundance(1,nyr,1,nps,1,nag)
 3darray true_survey_total_abundance_temp(1,nyr,1,nag,1,nps)
 matrix true_survey_total_abundance(1,nyr,1,nag)
 // end new tag stuff


 5darray survey_at_age_fleet_prop(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 5darray SIM_survey_prop(1,nps,1,nr,1,nfls,1,nyr,1,nag)
  //5darray OBS_survey_prop(1,nps,1,nr,1,nfls,1,nyr,1,nag)
 5darray OBS_survey_prop(1,nps,1,nr,1,nyr,1,nfls,1,nag) //fixed **KB

 5darray true_survey_fleet_age_bio(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 4darray true_survey_fleet_bio(1,nps,1,nr,1,nyr,1,nfls)
 3darray true_survey_region_bio(1,nps,1,nyr,1,nr)
 matrix true_survey_population_bio(1,nyr,1,nps)
 vector true_survey_total_bio(1,nyr)
 5darray OBS_survey_fleet_bio_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls)
 4darray OBS_survey_region_bio_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray OBS_survey_population_bio_overlap(1,nps,1,nyr,1,nps)
 matrix OBS_survey_natal_bio_overlap(1,nyr,1,nps)
 vector OBS_survey_total_bio_overlap(1,nyr)
 4darray OBS_survey_fleet_bio(1,nps,1,nr,1,nyr,1,nfls)
 3darray OBS_survey_region_bio(1,nps,1,nyr,1,nr)
 matrix OBS_survey_population_bio(1,nyr,1,nps)
 vector OBS_survey_total_bio(1,nyr)
 3darray apport_region_survey_biomass(1,nps,1,nr,1,nyr)

 //yield & BRP calcs 
 5darray catch_at_age_fleet(1,nps,1,nr,1,nyr,1,nag,1,nfl)
 5darray catch_at_age_fleet_prop(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 5darray SIM_catch_prop(1,nps,1,nr,1,nfl,1,nyr,1,nag)
  //5darray OBS_catch_prop(1,nps,1,nr,1,nfl,1,nyr,1,nag)
 5darray OBS_catch_prop(1,nps,1,nr,1,nyr,1,nfl,1,nag)//**fixed KB

 4darray yield_fleet(1,nps,1,nr,1,nyr,1,nfl)
 4darray catch_at_age_region(1,nps,1,nr,1,nyr,1,nag)
 4darray catch_at_age_region_prop(1,nps,1,nr,1,nyr,1,nag)
 3darray yield_region(1,nps,1,nr,1,nyr)
 3darray catch_at_age_population(1,nps,1,nyr,1,nag)
 3darray catch_at_age_population_prop(1,nps,1,nyr,1,nag)
 matrix yield_population(1,nps,1,nyr)
 3darray SSB_region(1,nps,1,nr,1,nyr)
 matrix SSB_population(1,nps,1,nyr)
 vector SSB_total(1,nyr)
 3darray abundance_population(1,nps,1,nyr,1,nag)
 matrix abundance_total(1,nyr,1,nag)
 matrix biomass_population(1,nps,1,nyr)
 vector biomass_total(1,nyr)
 matrix catch_at_age_total(1,nyr,1,nag)
 matrix catch_at_age_total_prop(1,nyr,1,nag)
 vector yield_total(1,nyr)
 4darray harvest_rate_region_num(1,nps,1,nr,1,nyr,1,nag)
 3darray harvest_rate_population_num(1,nps,1,nyr,1,nag)
 matrix harvest_rate_total_num(1,nyr,1,nag)
 3darray harvest_rate_region_bio(1,nps,1,nr,1,nyr)
 matrix harvest_rate_population_bio(1,nps,1,nyr)
 vector harvest_rate_total_bio(1,nyr)
 3darray depletion_region(1,nps,1,nr,1,nyr)
 matrix depletion_population(1,nps,1,nyr)
 vector depletion_total(1,nyr)
 5darray abundance_at_age_BM_overlap_region(1,nps,1,nps,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_BM_overlap_population(1,nps,1,nps,1,nyr,1,nag)
 5darray abundance_at_age_AM_overlap_region(1,nps,1,nps,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_AM_overlap_population(1,nps,1,nps,1,nyr,1,nag)
 4darray abundance_AM_overlap_region_all_natal(1,nps,1,nr,1,nyr,1,nag)
 5darray abundance_spawn_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray SSB_region_overlap(1,nps,1,nps,1,nr,1,nyr)
 3darray SSB_population_overlap(1,nps,1,nps,1,nyr)
 matrix SSB_natal_overlap(1,nps,1,nyr)

 6darray catch_at_age_region_fleet_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray catch_at_age_region_fleet_overlap_prop(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray SIM_catch_prop_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray OBS_catch_prop_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 5darray catch_at_age_region_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 5darray catch_at_age_region_overlap_prop(1,nps,1,nps,1,nr,1,nyr,1,nag)
 5darray yield_region_fleet_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr)
 4darray yield_region_overlap(1,nps,1,nps,1,nr,1,nyr)
 4darray catch_at_age_population_overlap(1,nps,1,nps,1,nyr,1,nag)
 4darray catch_at_age_population_overlap_prop(1,nps,1,nps,1,nyr,1,nag)
 3darray yield_population_overlap(1,nps,1,nps,1,nyr)
 3darray abundance_natal_overlap(1,nps,1,nyr,1,nag)
 5darray biomass_BM_age_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 5darray biomass_AM_age_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray biomass_BM_overlap_region(1,nps,1,nps,1,nr,1,nyr)
 4darray biomass_AM_overlap_region(1,nps,1,nps,1,nr,1,nyr)
 5darray biomass_AM_overlap_region_all_natal_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 4darray biomass_AM_overlap_age_region_all_natal(1,nps,1,nr,1,nyr,1,nag)
 3darray biomass_AM_overlap_region_all_natal(1,nps,1,nr,1,nyr)
 3darray biomass_population_overlap(1,nps,1,nps,1,nyr)
 matrix biomass_natal_overlap(1,nps,1,nyr)
 3darray catch_at_age_natal_overlap(1,nps,1,nyr,1,nag)
 3darray catch_at_age_natal_overlap_prop(1,nps,1,nyr,1,nag)
 matrix yield_natal_overlap(1,nps,1,nyr)
 5darray harvest_rate_region_fleet_bio_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr)
 4darray harvest_rate_region_bio_overlap(1,nps,1,nps,1,nr,1,nyr)
 3darray harvest_rate_population_bio_overlap(1,nps,1,nps,1,nyr)
 matrix harvest_rate_natal_bio_overlap(1,nps,1,nyr)
 4darray depletion_region_overlap(1,nps,1,nps,1,nr,1,nyr)
 3darray depletion_population_overlap(1,nps,1,nps,1,nyr)
 matrix depletion_natal_overlap(1,nps,1,nyr)
 3darray Bratio_population_overlap(1,nps,1,nps,1,nyr)
 matrix Bratio_natal_overlap(1,nps,1,nyr)
 matrix Bratio_population(1,nps,1,nyr)
 vector Bratio_total(1,nyr)

 //Observed Yield
 5darray OBS_yield_region_fleet_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfl)
 4darray OBS_yield_region_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray OBS_yield_population_overlap(1,nps,1,nyr,1,nps)
 matrix OBS_yield_natal_overlap(1,nyr,1,nps)
 vector OBS_yield_total_overlap(1,nyr)
 4darray OBS_yield_fleet(1,nps,1,nr,1,nyr,1,nfl)
 3darray OBS_yield_region(1,nps,1,nyr,1,nr)
 matrix OBS_yield_population(1,nyr,1,nps)
 vector OBS_yield_total(1,nyr)
 3darray apport_yield_region(1,nps,1,nr,1,nyr)

 matrix biomass_BM_temp(1,nps,1,nr)
 number biomass_BM_temp2
 5darray biomass_BM_overlap_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 4darray init_abund_temp(1,nps,1,nr,1,nag,1,nps)
 5darray rand_SIM_survey_prop_temp(1,nps,1,nr,1,nyr,1,nfls,1,2000) //should make function of max(ncatch) but had issues making an index, used 2000 as placeholder since nsurvey unlikely to exceed 2000
 6darray rand_SIM_survey_prop_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,2000) //should make function of max(ncatch) but had issues making an index, used 2000 as placeholder since nsurvey unlikely to exceed 2000
 5darray rand_SIM_catch_prop_temp(1,nps,1,nr,1,nyr,1,nfl,1,2000) //should make function of max(ncatch) but had issues making an index
 6darray rand_SIM_catch_prop_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfl,1,2000) //should make function of max(ncatch) but had issues making an index
 vector rand_SIM_survey_prop_temp2(1,nages)
 vector rand_SIM_catch_prop_temp2(1,nages)
 5darray OBS_survey_fleet_bio_temp(1,nps,1,nr,1,nyr,1,nfls,1,nps)
 5darray true_survey_fleet_bio_overlap_temp(1,nps,1,nr,1,nyr,1,nfls,1,nps)
 5darray catch_at_age_fleet_prop_temp(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 matrix abundance_move_temp(1,nps,1,nr)
 matrix bio_move_temp(1,nps,1,nr)
 5darray yield_fleet_temp(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 4darray yield_region_temp(1,nps,1,nr,1,nyr,1,nag)
 3darray yield_population_temp(1,nps,1,nyr,1,nag)
 4darray SSB_region_temp(1,nps,1,nr,1,nyr,1,nag)
 matrix SSB_total_temp(1,nyr,1,nps)
 4darray abundance_population_temp(1,nps,1,nyr,1,nag,1,nr)
 3darray abundance_total_temp(1,nyr,1,nag,1,nps)
 3darray biomass_population_temp(1,nps,1,nyr,1,nr)
 matrix biomass_total_temp(1,nyr,1,nps)
 3darray catch_at_age_total_temp(1,nyr,1,nag,1,nps)
 matrix yield_total_temp(1,nyr,1,nps)
 4darray catch_at_age_population_temp(1,nps,1,nyr,1,nag,1,nr)
 5darray SSB_region_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 matrix abundance_move_overlap_temp(1,nps,1,nr)
 5darray OBS_yield_fleet_temp(1,nps,1,nr,1,nyr,1,nfl,1,nps)
 6darray yield_region_fleet_temp_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 5darray yield_region_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray yield_population_temp_overlap(1,nps,1,nps,1,nyr,1,nag)
 4darray abundance_natal_temp_overlap(1,nps,1,nyr,1,nag,1,nps)
 4darray biomass_population_temp_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray biomass_natal_temp_overlap(1,nps,1,nyr,1,nps)
 4darray catch_at_age_natal_temp_overlap(1,nps,1,nyr,1,nag,1,nps)
 3darray yield_natal_temp_overlap(1,nps,1,nyr,1,nps)
 5darray catch_at_age_population_temp_overlap(1,nps,1,nps,1,nyr,1,nag,1,nr)
 3darray SSB_natal_overlap_temp(1,nps,1,nyr,1,nps)
 matrix SSB_overlap_natal(1,nps,1,nr)
 5darray abundance_AM_overlap_region_all_natal_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 3darray SSB_population_temp(1,nps,1,nyr,1,nr)
 4darray SSB_population_temp_overlap(1,nps,1,nps,1,nyr,1,nr)
 
 4darray res_TAC(1,nps,1,nr,1,nfl,1,nyr)
 3darray res_u(1,nps,1,nr,1,nyr)
 number Fnew
 number delt
 number fofF
 number fprimeF
 vector fofFvect(1,nag)
 vector fprimeFhigh(1,nag)
 vector fprimeFlow(1,nag)
 4darray TAC(1,nps,1,nr,1,nfl,1,nyr)
 3darray u(1,nps,1,nr,1,nfl)
 4darray yield_RN(1,nps,1,nr,1,nyr,1,nfl)
 5darray yield_RN_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfl)
 4darray survey_RN(1,nps,1,nr,1,nyr,1,nfls)
 5darray survey_RN_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls)
 4darray F_RN(1,nps,1,nr,1,nyr,1,nfl)
 4darray T_RN(1,nps,1,nr,1,nyr,1,nag)
 matrix rec_devs_RN(1,nps,1,nyr)
 3darray Rec_apport_RN(1,nps,1,nyr-1,1,nr)
 3darray rec_index_RN(1,nps,1,nr,1,nyr)

 4darray yield_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_flts)
 5darray yield_RN_temp_overlap(1,max_pops,1,max_pops,1,max_regs,1,max_yrs,1,max_flts)
 4darray survey_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_surv_flts)
 5darray survey_RN_temp_overlap(1,max_pops,1,max_pops,1,max_regs,1,max_yrs,1,max_surv_flts)
 4darray F_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_flts)
 4darray T_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_ages)
 matrix rec_devs_RN_temp(1,max_pops,1,max_yrs)
 3darray Rec_apport_RN_temp(1,max_pops,1,max_yrs,1,max_regs)
 3darray rec_index_RN_temp(1,max_pops,1,max_regs,1,max_yrs)
 3darray prob_tag_RN_temp(1,max_pops,1,max_regs,1,max_tag_yrs)
 3darray ntags_RN_temp(1,max_pops,1,max_regs,1,max_tag_yrs)
 vector prob_tag_year_RN_temp(1,max_tag_yrs)
 
 init_number dummy(phase_dummy)

   //Setting up parameters for export to mismatch EM model

 4darray abund_frac_age_region(1,nps,1,nr,1,nyr,1,nag)
 3darray abund_frac_region_year(1,nps,1,nr,1,nyr)
 matrix abund_frac_region(1,nps,1,nr)

 3darray input_weight_region_temp(1,nps,1,nag,1,nr)
 matrix input_weight_region(1,nps,1,nag)
 matrix input_weight_population_temp(1,nag,1,nps)
 vector input_weight_population(1,nag)
 3darray input_catch_weight_region_temp(1,nps,1,nag,1,nr)
 matrix input_catch_weight_region(1,nps,1,nag)
 matrix input_catch_weight_population_temp(1,nag,1,nps)
 vector input_catch_weight_population(1,nag)
 3darray fecundity_region_temp(1,nps,1,nag,1,nr)
 matrix fecundity_region(1,nps,1,nag)
 matrix fecundity_population_temp(1,nag,1,nps)
 vector fecundity_population(1,nag)
 3darray maturity_region_temp(1,nps,1,nag,1,nr)
 matrix maturity_region(1,nps,1,nag)
 matrix maturity_population_temp(1,nag,1,nps)
 vector maturity_population(1,nag)
 
 number prop_fem_pan
 matrix prop_fem_temp(1,nps,1,nr)

 matrix rec_index_BM_population(1,nps,1,nyr)
 vector rec_index_pan(1,nyr)
 3darray rec_index_temp(1,nps,1,nyr,1,nr)
 matrix rec_index_temp2(1,nyr,1,nps)

 4darray selectivity_age_temp(1,nps,1,nag,1,nfl,1,nr)
 3darray selectivity_age_pop(1,nps,1,nag,1,nfl)
 4darray survey_selectivity_age_temp(1,nps,1,nag,1,nfl,1,nr)
 3darray survey_selectivity_age_pop(1,nps,1,nag,1,nfl)

//weighted average of age comps
 5darray OBS_survey_prop_temp(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray OBS_survey_prop_temp2(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_survey_prop_temp3(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_survey_prop_temp4(1,nps,1,nyr,1,nag,1,nr)
 3darray OBS_survey_prop_population(1,nps,1,nyr,1,nag)
 3darray OBS_survey_prop_pan_temp(1,nyr,1,nag,1,nps)
 matrix  OBS_survey_prop_pan(1,nyr,1,nag)

 5darray OBS_catch_prop_temp(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray OBS_catch_prop_temp2(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_catch_prop_temp3(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_catch_prop_temp4(1,nps,1,nyr,1,nag,1,nr)
 3darray OBS_catch_prop_population(1,nps,1,nyr,1,nag)
 3darray OBS_catch_prop_pan_temp(1,nyr,1,nag,1,nps)
 matrix  OBS_catch_prop_pan(1,nyr,1,nag)


//summing ntags
 4darray ntags_temp1(1,nps,1,nyr_rel,1,nag,1,nr)
 3darray ntags_population(1,nps,1,nyr_rel,1,nag)
 3darray ntags_pan_temp(1,nyr_rel,1,nag,1,nps)
 matrix  ntags_pan(1,nyr_rel,1,nag)


//summing OBS tag prop 
5darray OBS_tag_prop_population_temp(1,nps,1,nyr_rel,1,nag,1,nt,1,nr)
4darray OBS_tag_prop_population_temp2(1,nps,1,nyr_rel,1,nag,1,nt)
4darray OBS_tag_prop_population_final(1,nps,1,nyr_rel,1,nag,1,nt)
4darray OBS_tag_prop_pan_temp(1,nyr_rel,1,nag,1,nt,1,nps)
3darray OBS_tag_prop_pan_temp2(1,nyr_rel,1,nag,1,nt)
3darray OBS_tag_prop_pan_final_temp(1,nyr_rel,1,nag,1,max_life_tags+1)//set up a new array for summed proportion
3darray OBS_tag_prop_pan_final(1,nyr_rel,1,nag,1,max_life_tags+1)//set up a new array for summed proportion

//4darray OBS_tag_prop_population_temp_no_age(1,nps,1,nyr_rel,1,nt,1,nr)
//3darray OBS_tag_prop_population_final_no_age(1,nps,1,nyr_rel,1,nt)
//      !!cout<<ntags_region<<endl;

//3darray OBS_tag_prop_population_temp2_no_age(1,nps,1,nyr_rel,1,nt)

//3darray OBS_tag_prop_pan_temp_no_age(1,nyr_rel,1,nt,1,nps)
//matrix OBS_tag_prop_pan_temp2_no_age(1,nyr_rel,1,nt)
//matrix OBS_tag_prop_pan_final_temp_no_age(1,nyr_rel,1,max_life_tags+1)//set up a new array for summed proportion
//matrix OBS_tag_prop_pan_final_no_age(1,nyr_rel,1,max_life_tags+1)//set up a new array for summed proportion

 //calculating true fraction natal
 3darray init_abund_reg_temp(1,nps,1,nps,1,nr)
 matrix init_abund_pop_temp(1,nps,1,nps)
 3darray frac_natal_true(1,nps,1,nps,1,nr)

///////////////////////////////////////////////////////////////////////////////////////////////
//Temporary (reorganized) 6d arrary parameters 
//Reorganize so that region is the second dimension. 
//////////////////////////////////////////////////////////////////////////////////////////////
//survey index  
 //6darray ro_true_survey_fleet_overlap_age(1,nps,1,nr,1,nps,1,nyr,1,nfls,1,nag) 
 //6darray ro_survey_at_age_region_fleet_overlap_prop(1,nps,1,nr,1,nps,1,nfls,1,nyr,1,nag) 
// 6darray ro_SIM_survey_prop_overlap(1,nps,1,nr,1,nps,1,nfls,1,nyr,1,nag)
// 6darray ro_OBS_survey_prop_overlap(1,nps,1,nr,1,nps,1,nfls,1,nyr,1,nag)
// 6darray ro_true_survey_fleet_overlap_age_bio(1,nps,1,nr,1,nps,1,nyr,1,nfls,1,nag) 


//yield & BRP calcs 
// 6darray ro_catch_at_age_region_fleet_overlap(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag) 
// 6darray ro_catch_at_age_region_fleet_overlap_prop(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag)
// 6darray ro_SIM_catch_prop_overlap(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag) 
// 6darray ro_OBS_catch_prop_overlap(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag) 


  objective_function_value f

 !! cout << "parameters set" << endl;
 
   
PROCEDURE_SECTION

  get_random_numbers();

  get_movement();

  get_selectivity();

  get_F_age();

  get_vitals();

  get_SPR();

  get_env_Rec();

  get_DD_move_parameters();

  get_abundance();

  get_rand_survey_CAA_prop();

  get_rand_CAA_prop();

  get_tag_recaptures();

  get_observed_tag_recaptures();

  evaluate_the_objective_function();

FUNCTION get_random_numbers

 // NOTE THAT IT APPEARS THE RNG JUST CONTINUES A SINGLE SET OF RANDOM NUMBERS WITHIN A DATA OBJECT
 // THUS FOR EACH DATA SET OR SIMULATED PARAMETER THE VALUES WILL NOT REPEAT SO DON'T NEED MULTIPLE RANDOM
 // NUMBERS FOR EACH DATA SET, ETC...

   random_number_generator myrand_yield(myseed_yield);
   random_number_generator myrand_survey(myseed_survey);
   random_number_generator myrand_F(myseed_F);
   random_number_generator myrand_T(myseed_T);
   random_number_generator myrand_rec_devs(myseed_rec_devs);
   random_number_generator myrand_rec_apport(myseed_rec_apport);
   random_number_generator myrand_rec_index(myseed_rec_index);
   random_number_generator myrand_ntags(myseed_ntags);
   random_number_generator myrand_prob_tag(myseed_prob_tag);
   random_number_generator myrand_prob_tag_year(myseed_prob_tag_year);

 //to avoid issues of varying RNGs across scenarios/runs (due to different input index lengths) we generate the RNG
 //and put them into temp arrays that using fixed index legnths meant to be greater than the maximum possible
 //values that would be used for a given index (note that there is no error protection built in so that if someone
 //uses a higher index than input as the max value, then RNGs will repeat)
 //the values from the temp arrays are then read into the actual RN arrays
 //this way if one run has a different index it will have the same RNGs up until the new dimension (e.g., if one run
 //has 3 fleets and another has 4, the one with 4 fleets with have the same RNGs for the first 3 fleets, then will
 //keep reading the RNG string for the 4th fleet)
  for (int p=1;p<=max_pops;p++)
   {
    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {       
        for (int y=1;y<=max_yrs;y++)
         {
          for (int z=1;z<=max_flts;z++)
           {
            for (int x=1;x<=max_surv_flts;x++)
             {
            for(int s=1; s<=max_tag_yrs; s++)
               {
             for (int a=1;a<=max_ages;a++)
                {
                 yield_RN_temp_overlap(p,j,r,y,z)=randn(myrand_yield);
                 survey_RN_temp_overlap(p,j,r,y,x)=randn(myrand_survey);
                 yield_RN_temp(j,r,y,z)=randn(myrand_yield);
                 survey_RN_temp(j,r,y,x)=randn(myrand_survey);
                 F_RN_temp(j,r,y,z)=randn(myrand_F);
                 rec_devs_RN_temp(j,y)=randn(myrand_rec_devs);
                 rec_index_RN_temp(j,r,y)=randn(myrand_rec_index);
                 T_RN_temp(j,r,y,a)=randn(myrand_T);
                 ntags_RN_temp(j,r,s)=randu(myrand_ntags);
                 prob_tag_RN_temp(j,r,s)=randu(myrand_prob_tag);
                 prob_tag_year_RN_temp(s)=randu(myrand_prob_tag_year);
               if(apportionment_type==3)//completely random apportionment
                {
                 Rec_apport_RN_temp(j,y,r)=randu(myrand_rec_apport);//generate a positive random number bw 0-1
                }
               if(apportionment_type==4)//completely random apportionment
                {
                 Rec_apport_RN_temp(j,y,r)=randn(myrand_rec_apport);//generate a positive random number bw 0-1
               }
         }
        }
       }
      }
     }
    }
   }
  }
   

 //once RNGs are created they are now input into the actual RN arrays used in the calculations based on the true
 //index lengths (
 //*******NOTE*********
 //if index length > input associated max value used in previous loop then will have issues with repeating RNGs
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
          for (int z=1;z<=nfleets(j);z++)
           {
            for (int x=1;x<=nfleets_survey(j);x++)
             {
              for(int s=1; s<=nyrs_release; s++)
               {
                for (int a=1;a<=nages;a++)
                {
                 yield_RN(j,r,y,z)=yield_RN_temp(j,r,y,z);
                 yield_RN_overlap(p,j,r,y,z)=yield_RN_temp_overlap(p,j,r,y,z);
                 survey_RN(j,r,y,x)=survey_RN_temp(j,r,y,x);
                 survey_RN_overlap(p,j,r,y,x)=survey_RN_temp_overlap(p,j,r,y,x);
                 F_RN(j,r,y,z)=F_RN_temp(j,r,y,z);
                 rec_devs_RN(j,y)=rec_devs_RN_temp(j,y);
                 rec_index_RN(j,r,y)=rec_index_RN_temp(j,r,y);
                 T_RN(j,r,y,a)=T_RN_temp(j,r,y,a);
                 ntags_RN(j,r,s)=ntags_RN_temp(j,r,s);
                 prob_tag_RN(j,r,s)=prob_tag_RN_temp(j,r,s);
                 prob_tag_year_RN(s)=prob_tag_year_RN_temp(s);
             if(y>1) //rec apport is of length nyrs-1, so need to make sure y>1
              {
               if(apportionment_type==3)//completely random apportionment
                {
                 Rec_apport_RN(j,y-1,r)=Rec_apport_RN_temp(j,y-1,r);//generate a positive random number bw 0-1
                }
               if(apportionment_type==4)//completely random apportionment
                {
                 Rec_apport_RN(j,y-1,r)=Rec_apport_RN_temp(j,y-1,r);//generate a positive random number bw 0-1
               }
              }
         }
        }
       }
      }
     }
    }
   }
  }
///////BUILD MOVEMENT MATRIX////////
FUNCTION get_movement

//POSSIBLE ADDITIONS:
  //new functional forms
  //simulate yearly differences and age-based movement
  //random movement same as a random F value
  //random walk
  //yearly deviationp around an input average value (e.g., unidirectional movement with yearly deviationp based on a variance term and lognormal distribution)
  //need to enpure random values still sum to 1 for an area (relativization of values)
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
                 if(move_switch==0)  //fix at no movement
                  {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1;
                    }
                   else
                    {
                     T(j,r,y,a,k,n)=0;
                    }
                  }
                 if(move_switch==1) // use input movement
                  {
                   T(j,r,y,a,k,n)=input_T(j,r,a,k,n);
                  }
                 if(move_switch==2) // only allow movement within a population (ie regions within a population) not across populations based on input residency term
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency(j,r,a);
                   }
                   if(j==k && r!=n)
                   {
                    T(j,r,y,a,k,n)=(1-input_residency(j,r,a))/(nregions(j)-1);
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(move_switch==3) //symmetric movement but only allow movement within a population (ie regionp within a population) not across populations based on input residency term
                  {
                   if(j==k)
                   {
                    T(j,r,y,a,k,n)=1/(nregions(j));
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(move_switch==4) //symmetric movement across all populations and regions
                  {
                   T(j,r,y,a,k,n)=1/sum(nregions);
                  }
                 if(move_switch==5) // allow movement across all regions and populations, based on population/region specific residency
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency(j,r,a);
                   }
                   else
                   {
                    T(j,r,y,a,k,n)=(1-input_residency(j,r,a))/(sum(nregions)-1);
                   }
                  }
                 if(move_switch==9) // use input yearly  movement
                  {
                   T(j,r,y,a,k,n)=input_T_year(j,r,y,k,n);
                  }

                 if(move_switch==10) // use input yearly  movement
                  {
                   T(j,r,y,a,k,n)=T_Full_Input(j,r,a+(y-1)*nages,k,n);
                  }
              if(larval_move_switch>(-1)) //fill in age-1 movement if age-1 movement is different than other ages
              {
               if(a==1 && larval_move_switch==0) // allow different movement from adults
                {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1;
                    }
                   else
                    {
                     T(j,r,y,a,k,n)=0;
                    }
                  }
                 if(a==1 && larval_move_switch==1) // use input movement
                  {
                   T(j,r,y,a,k,n)=input_T(j,r,a,k,n);
                  }
                 if(a==1 && larval_move_switch==2) // only allow movement within a population (ie regionp within a population) not across populations based on input residency term
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency_larval(j,r);
                   }
                   if(j==k && r!=n)
                   {
                    T(j,r,y,a,k,n)=(1-input_residency_larval(j,r))/(nregions(j)-1);
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(a==1 && larval_move_switch==3) //symmetric movement but only allow movement within a population (ie regionp within a population) not across populations
                  {
                   if(j==k)
                   {
                    T(j,r,y,a,k,n)=1/(nregions(j));
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(a==1 && larval_move_switch==4) //symmetric movement across all populations and regionp
                  {
                   T(j,r,y,a,k,n)=1/sum(nregions);
                  }
                 if(larval_move_switch==5) // allow movement across all regions and populations, based on population/region specific residency
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency_larval(j,r);
                   }
                   else
                   {
                    T(j,r,y,a,k,n)=(1-input_residency_larval(j,r))/(sum(nregions)-1);
                   }
                  }
                 }
              if(T(j,r,y,a,k,n)>1) //ensure T doesn't exceed bounds
                {
                 T(j,r,y,a,k,n)=1;
                }
               if(T(j,r,y,a,k,n)<0)
                {
                 T(j,r,y,a,k,n)=0;
                }
       }
      } 
     }
    }
   }
  }

 if(rand_move==1)
  {
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              if((a>1 && move_switch!=0) || (a==1 && larval_move_switch>0) || (a==1 && larval_move_switch==(-1) &&  move_switch!=0)) 
              {
               T(j,r,y,a,k,n)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j)));
              }
               if(move_switch==0)
                {
                 T(j,r,y,a,k,n)=T(j,r,y,a,k,n);
                }
               if(a==1 && larval_move_switch==0)
                {
                 T(j,r,y,a,k,n)=T(j,r,y,a,k,n);
                }
               if(T(j,r,y,a,k,n)>1) //ensure T doesn't exceed bounds
                {
                 T(j,r,y,a,k,n)=1;
                }
               if(T(j,r,y,a,k,n)<0)
                {
                 T(j,r,y,a,k,n)=0;
                }
              }
             }
            }
           }
          }
     }
    }

///////SELECTIVITY CALCULATIONS///////
FUNCTION get_selectivity
//POSSIBLE ADDITIONS:
  //yearly selectivity
  
//fishery selectivity
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
        {
        age_full_selection_temp=0;
         for (int a=1;a<=nages;a++)
           {
            for (int z=1;z<=nfleets(j);z++)
              {
               if(select_switch==2) //4 parameter double logistic selectivity
                {
                 selectivity(j,r,y,a,z)=1/((1+mfexp(-sel_beta1(j,r,z)*(a-sel_beta2(j,r,z))))*(1+mfexp(-sel_beta3(j,r,z)*(a-sel_beta4(j,r,z)))));
                }
                if(select_switch==1) //two parameter logistic selectivity
                {
                 selectivity(j,r,y,a,z)=1/(1+mfexp(-sel_beta1(j,r,z)*(a-sel_beta2(j,r,z))));
                //selectivity(j,r,y,a,z)=1/(1+mfexp(-log(19)*(a-(sel_beta1(j,r,z)))/(sel_beta2(j,r,z)))); 
                }
                if(select_switch==0) //input selectivity at age constant by year
                {
                 selectivity(j,r,y,a,z)=input_selectivity(j,r,a,z);
                }
               // age_full_selection_temp(a)=selectivity(j,r,y,a,z);
               }
              }
   //    if(tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
     //    {
       //  for (int a=1;a<=nages;a++)
         //  {
           // for (int z=1;z<=nfleets(j);z++)
             // {
               //if(selectivity(j,r,y,a,z)==max(age_full_selection_temp)); //really only getting max age of last fleet but good enough for now
                //{
                 age_full_selection(j,r,y)=tag_age_sel;
                }
               }
              }


//survey selectivity 
 for (int j=1;j<=npops;j++)
    {
     for (int r=1;r<=nregions(j);r++)
      {
       for (int y=1;y<=nyrs;y++)
         {
          for (int a=1;a<=nages;a++)
            {
             for (int z=1;z<=nfleets_survey(j);z++)
               {

               if(select_switch_survey==2) //4 parameter double logistic selectivity
                {
                 survey_selectivity(j,r,y,a,z)=1/((1+mfexp(-sel_beta1_survey(j,r,z)*(a-sel_beta2_survey(j,r,z))))*(1+mfexp(-sel_beta3_survey(j,r,z)*(a-sel_beta4_survey(j,r,z)))));
                }
                if(select_switch_survey==1) //two parameter logistic selectivity
                {
                 survey_selectivity(j,r,y,a,z)=1/(1+mfexp(-sel_beta1_survey(j,r,z)*(a-sel_beta2_survey(j,r,z))));
                //selectivity(j,r,y,a,z)=1/(1+mfexp(-log(19)*(a-(sel_beta1(j,r,z)))/(sel_beta2(j,r,z)))); 
                }
                if(select_switch_survey==0) //input selectivity
                {
                 survey_selectivity(j,r,y,a,z)=input_survey_selectivity(j,r,a,z);
                }
              }
             }
            }
          }
        }

  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
          for (int a=1;a<=nages;a++)
            {
            for (int z=1;z<=nfleets(j);z++)
              {
                selectivity_age(j,r,a,z)=selectivity(j,r,nyrs,a,z);
              }
             for (int z=1;z<=nfleets_survey(j);z++)
               {
                survey_selectivity_age(j,r,a,z)=survey_selectivity(j,r,nyrs,a,z);
               }
              }
             }
            }


 
///////FISHING MORTALITY CALCULATIONS///////
FUNCTION get_F_age

  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int z=1;z<=nfleets(j);z++)
           { 
             if(F_switch==1) //input F directly
              {
               F_year(j,r,y,z)=input_F(j,r,z);
              }
             if(F_switch==2) //input single yearly F (i.e., using single population FMSY for all populations)
              {
               F_year(j,r,y,z)=input_F_MSY;
              }
             if(F_switch==3) //estimate F that achieves desired Total F
              {
               F_year(j,r,y,z)=F_est(j,r);
              }
             if(F_switch==4) //split F_MSY by npops
              {
               F_year(j,r,y,z)=input_F_MSY/npops;
              }
             if(F_switch==5) //split F_MSY by total number of regionp
              {
               F_year(j,r,y,z)=input_F_MSY/sum(nregions);
              }
             if(F_switch==6) //split F_MSY by total number of fleets
              {
               F_year(j,r,y,z)=input_F_MSY/sum(nfleets);
              }
             if(F_switch==7) //F devs about input F based on sigma_F
              {
               F_year(j,r,y,z)=input_F(j,r,z)*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
              }
             if(F_switch==8) //random walk in F
              {
               F_year(j,r,y,z)=F_rho(j,r,z)*F_year(j,r,y-1,z)*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
              }
             if(F_switch==9)  //Dunce cap F
              {
               Fstartyr(j,r)=dunce_F(j,r,1);
               minF(j,r)=dunce_F(j,r,2);
               maxF(j,r)=dunce_F(j,r,3);
               stepF(j,r)=(maxF(j,r)-minF(j,r))/((nyrs-Fstartyr(j,r))/2);
              if(y<Fstartyr(j,r))
               {
                F_year(j,r,y,z)=0;
               }
               if(y>=Fstartyr(j,r))
               {
               if(y<((nyrs-Fstartyr(j,r))/2+Fstartyr(j,r)))
                {
                 F_year(j,r,y,z)=minF(j,r)+(y-Fstartyr(j,r))*stepF(j,r)*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
                }
               }
               if(y>=((nyrs-Fstartyr(j,r))/2+Fstartyr(j,r)))
                {
                 F_year(j,r,y,z)=maxF(j,r)-((y-Fstartyr(j,r))-((nyrs-Fstartyr(j,r))/2))*stepF(j,r)*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
                  if(F_year(j,r,y,z)<minF(j,r)) //needed because the stepF decrease can be randomly be greater than preceding F, and so F goes negative
                  {
                  F_year(j,r,y,z)=0.5*minF(j,r);
                  }
                }
              }
             F_fleet(j,r,y,a,z)=F_year(j,r,y,z)*selectivity(j,r,y,a,z);
             F(j,r,y,a)=sum(F_fleet(j,r,y,a)); 
             M(j,r,y,a)=input_M_TRUE(j,a);
           }
          }
         }
        }
       }
FUNCTION get_vitals
//POSSIBLE ADDITIONS:
  //random walk in apportionment or random to give time-varying
  //switch for input recruitment devs by year to recreate a given population trajectory
 R_ave=ln_R_ave; ///this is annoying...//a quick fix
 
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
          for (int a=1;a<=nages;a++)
           {
            for (int z=1;z<=nfleets(j);z++)
             {
              weight_population(j,r,y,a)=input_weight(j,r,a);
              weight_catch(j,r,y,a)=input_catch_weight(j,r,a);

              if(maturity_switch_equil==0) // for SPR calculations when maturity across areas is equal or if want a straight average of maturity across areas
               {
                if(SSB_type==1) //fecundity based SSB
                 {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*fecundity(j,r,a)*maturity(j,r,a);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                 }
               if(SSB_type==2) //weight based SSB
                {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*weight_population(j,r,y,a)*maturity(j,r,a);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                }
               }
              if(maturity_switch_equil==1)
               {// calculates the weighted average matruity based on equilibrium apportionment of SSB - allows for unequal influence of maturity/weight
                if(SSB_type==1) //fecundity based SSB
                 {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*fecundity(j,r,a)*maturity(j,r,a)*equil_ssb_apport(j,r);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                 }
               if(SSB_type==2) //weight based SSB
                {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*weight_population(j,r,y,a)*maturity(j,r,a);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                }
               }                        
               if(recruit_devs_switch==0)  //use population recruit relationship directly
                {
                 rec_devs(j,y)=1;
                }
               if(recruit_devs_switch==1)  // allow lognormal error around SR curve
                {
                 rec_devs(j,y)=mfexp(rec_devs_RN(j,y)*sigma_recruit(j)-.5*square(sigma_recruit(j)));

                 if(recruit_randwalk_switch==1)
                 {
                  rec_devs_randwalk(j,y)=rec_devs(j,y);
                  if(y!=1)
                   {
                    rec_devs(j,y)=rec_devs(j,y-1)*rec_devs_randwalk(j,y); //is this correct?
                   }
                 }
                }
               if(SSB_type==1) //fecundity based SSB
                {
                 wt_mat_mult_reg(j,r,y,a)=prop_fem(j,r)*fecundity(j,r,a)*maturity(j,r,a);// for yearly SSB calcs
                }
               if(SSB_type==2) //weight based SSB
                {
                 wt_mat_mult_reg(j,r,y,a)=prop_fem(j,r)*weight_population(j,r,y,a)*maturity(j,r,a);
                }
             }
            }
           }

        for (int y=1;y<=nyrs-1;y++)
         {
               if(apportionment_type==1) //input recruitment apportionment directly by population and region
                {
                 Rec_Prop(j,r,y)=input_Rec_prop(j,r);
                }
               if(apportionment_type==2) //equal apportionment by nregions
                {
               Rec_Prop(j,r,y)=1.0/nregions(j);
                }
               if(apportionment_type==(-1)) // no apportionment
                {
                 Rec_Prop(j,r,y)=1;
                }
               if(apportionment_type==3)//completely random apportionment
                {
                Rec_prop_temp1(j,y,r)=Rec_apport_RN(j,y,r);//generate a positive random number bw 0-1
                }                 
               if(apportionment_type==4) //add input obersvation error to input recruit proportions following Schnute and Richards 1995 (as implemented in Cox and Kronlund 2008)
                {
                 Rec_prop_temp1(j,y,r)=log(input_Rec_prop(j,r))+sigma_rec_prop(j)*Rec_apport_RN(j,y,r);//applying the additive error in log space; this equals "log(u) + error" in Cox and Kronlund Table 1
                }
               }   
             }
        
  for (int r=1;r<=nregions(j);r++)
       {
        for (int y=1;y<=nyrs-1;y++)
         {
          if(apportionment_type==3)
           { //need to standardize year by region matrix to ensure new proportions sum to one
            Rec_Prop(j,r,y)=Rec_prop_temp1(j,y,r)/sum(Rec_prop_temp1(j,y));
           }
          if(apportionment_type==4)
           { //need to run through region by year matrix to calculate second half of Schnute and Richards 1995 random mult equations and to standardize randomized apportioments to total one
            Rec_prop_temp2(j,y,r)=Rec_prop_temp1(j,y,r)-(sum(Rec_prop_temp1(j,y))/nregions(j)); //finish equation T1.21 in Cox and Kronlund Table 1 (note that here nregions is the same as A (number of ages) in Table 1 paper
            Rec_Prop(j,r,y)=mfexp(Rec_prop_temp2(j,y,r))/sum(mfexp(Rec_prop_temp2(j,y))); // back-transform and standardize
           }
         }
        }
       }
     }
//SPR calcs are done with either  average maturity/weight across all the regions within a population or assuming an input population fraction at equilibrium
// while the full SSB calcs use the region specific maturity/weight
FUNCTION get_SPR
     for (int k=1;k<=npops;k++)
     {
      for (int n=1;n<=nages;n++)
       {
        if(n==1)
         {
          SPR_N(k,n)=1000;
          SPR_SSB(k,n)=wt_mat_mult(k,1,n)*SPR_N(k,n);
         }
        if(n>1 && n<nages)
         {
          SPR_N(k,n)=SPR_N(k,n-1)*mfexp(-M(k,1,1,n-1));
          SPR_SSB(k,n)=wt_mat_mult(k,1,n)*SPR_N(k,n);
         }
        if(n==nages)
         {
          SPR_N(k,n)=SPR_N(k,n-1)*mfexp(-M(k,1,1,n))*(1/(1-mfexp(-M(k,1,1,n))));
          SPR_SSB(k,n)=wt_mat_mult(k,1,n)*SPR_N(k,n);
         }
       }
     SPR(k)=sum(SPR_SSB(k))/1000;
     SSB_zero(k)=SPR(k)*R_ave(k);
      if(Rec_type==2) //BH recruitment
      {
      alpha(k)=(SSB_zero(k)/R_ave(k))*((1-steep(k))/(4*steep(k)));//alternate parameterization
      beta(k)=(5*steep(k)-1)/(4*steep(k)*R_ave(k));
      }
    }


FUNCTION get_env_Rec // calculate autocorrelated recruitment - input period and amplitude
       for (int p=1;p<=npops;p++)
        {
        for (int j=1;j<=npops;j++)
         {
         for (int r=1;r<=nregions(j);r++)
          {
          for (int y=1;y<=nyrs;y++)
           {
           if(y==1)
             {
              env_rec(y)=input_init_abund(p,j,r,1);
             }           
           if(y>1)
            {
             env_rec(y) = R_ave(j)+(R_ave(j)*amplitude)*sin(((2*M_PI)/freq)*years(y)+freq);
            }
           }
         }
        }
       }
FUNCTION get_DD_move_parameters
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////DD Movement///////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there           
      if(move_switch==8)
       {
        for (int s=1;s<=npops;s++)
         {
         for (int a=1;a<=nages;a++)
          {
          for (int n=1;n<=nregions(s);n++)
           {       
            if(use_input_Bstar==0 && Rec_type==2) //set Bstar==SSB0
             {
              if(nregions(s)==1)
               {
                Bstar(s,n,a)=SSB_zero(s);
               }
              if(nregions(s)>1)
               {
                Bstar(s,n,a)=SSB_zero(s)*SSB_zero_appor(s,n);
               }
             }
            if(use_input_Bstar==1 || Rec_type!=2)
             {
              Bstar(s,n,a)=input_Bstar(s,n,a);
             }
            if(DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age)
             {
              c(s,n,a)=-log(A(s,n,a))/sum(Bstar(s,n));  ///Bstar is the point of inflection in T
             }
            if(DD_move_age_switch==1) //DD movement is not age based (based on total bio not bio at age)
             {
              c(s,n,a)=-log(A(s,n,a))/Bstar(s,n,a);  ///Bstar is the point of inflection in T
             }
          }
         }
        }
       }
FUNCTION get_abundance


       for (int y=1;y<=nyrs;y++)
        {

//need to add all_natal calcs for all values that might want to compare across an area, because regular non-natal homing calcs do not account for different weights across natal populationp so if if any part of calc involves weight it needs to use the biomass_all_natal value not just biomass_AM

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// Y==1 Overlap Calcs ///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       if(y==1)
         {
         for (int a=1;a<=nages;a++)
          {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {

                  if(init_abund_switch==1)
                  {
                  init_abund(p,j,r,a)=R_ave(p)*pow(mfexp(-(M(p,r,y,a))),a-1)*frac_natal(p,j,r);
                  }
                  
                  if(init_abund_switch==0) {
                  init_abund(p,j,r,a)=input_init_abund(p,j,r,a);
                  }
                
                    abundance_at_age_BM_overlap_region(p,j,y,a,r)=init_abund(p,j,r,a);
                    abundance_at_age_BM_overlap_population(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                    biomass_BM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                    biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////////////NON-NATAL Homing calcs /////////////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////    
                   init_abund_temp(j,r,a,p)=init_abund(p,j,r,a);
                   abundance_at_age_BM(j,r,y,a)=sum(init_abund_temp(j,r,a));
                   recruits_BM(j,r,y)=abundance_at_age_BM(j,r,y,1);
////////////////////////////////////////////////////////////////////////////////

                  if(natal_homing_switch==0)
                   {
                    biomass_BM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
                  if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {
                    biomass_BM_overlap_temp(j,r,y,a,p)=biomass_BM_age_overlap(p,j,r,y,a);
                    biomass_BM_age(j,r,y,a)=sum(biomass_BM_overlap_temp(j,r,y,a));
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
                   
 ///get year one recruitment index
               rec_index_BM(j,r,y) = recruits_BM(j,r,y)*mfexp(rec_index_RN(j,r,y)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));
               rec_index_BM_temp(j,y,r)=rec_index_BM(j,r,y);
               rec_index_prop_BM(j,r,y)=rec_index_BM(j,r,y)/sum(rec_index_BM_temp(j,y));
              }
             }
            }
           }
          } //close loops so have full biomass vectors filled in at start of DD movement calcs
         for (int a=1;a<=nages;a++)
          {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Density-Dependent Movement Calcs///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
            
                    if(move_switch==8 && DD_move_age_switch==1) //DD movement is age based (based on age based biomass not total bio)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM_age(j,r,y,a)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n,a);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r,a)); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                          if(a>1 || (a==1 && larval_move_switch!=0))
                           {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                           }
                            if(a==1 && larval_move_switch==0)
                             {
                              T(j,r,y,a,s,n)=T(j,r,y,a,s,n);
                             }
                          }
                        }
                    }

                    if(move_switch==8 && DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM(j,r,y)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,1)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM(s,n,y)/sum(Bstar(s,n));
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM(j,r,y)/sum(Bstar(j,r))); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                          if(a>1 || (a==1 && larval_move_switch!=0))
                           {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                           }
                            if(a==1 && larval_move_switch==0)
                             {
                              T(j,r,y,a,s,n)=T(j,r,y,a,s,n);
                             }
                          }
                        }
                    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {              
                    if(move_switch!=6 || move_switch!=7 || a==1)  //if movement is not type=6 or a==1 (and movement type 6)
                     {
                       abundance_move_overlap_temp(k,n)=init_abund(p,k,n,a)*T(p,n,y,a,j,r); //with overlap always use natal population movement rates  (i.e., use p inpsead of k)
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal population movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=init_abund(p,k,n,a)*return_probability(p); //with overlap always use natal population movement rates
                      }
                     }
                   }
                  }
                    if(move_switch!=6 || move_switch!=7  || a==1)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
                     }
                    if(move_switch==7)  //all fish stay where they were (i.e., no return migration)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=init_abund(p,j,r,a);
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a<return_age || a>return_age)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=init_abund(p,j,r,a); //with overlap always use natal population movement rates                     
                       }
                      if(a==return_age && p==j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=init_abund(p,j,r,a)+(sum(abundance_move_overlap_temp)/nregions(p));
                       }
                      if(a==return_age && p!=j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=(1-return_probability(p))*init_abund(p,j,r,a);
                       }
                      }
                      
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_population(p,r,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_population(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));
                biomass_AM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_population_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_population_overlap(p,j,y)=sum(biomass_population_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_population_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

///////////MIGHT WANT TO DO FOLLOWING CALCS FOR NATAL HOMING AS WELL (DO METAPOP TYPE CALCS BELOW)
              //  abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
              //  abundance_res(j,r,y,a)=abundance_move_temp(j,r);
              //  abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
              //  bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
              //  bio_res(j,r,y,a)=bio_move_temp(j,r);
              //  bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,y,a)-bio_res(j,r,y,a);
                
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance/////////////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////
    
                abundance_move_temp=0;
                bio_move_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=abundance_at_age_BM(k,n,y,a)*T(k,n,y,a,j,r);
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_population(k,n,y,a);                   
                   }
                  }
                  if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);
                   }
                  if(natal_homing_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_AM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }
                   
                abundance_population_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_population(j,y,a)=sum(abundance_population_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_population(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                biomass_population_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_population(j,y)=sum(biomass_population_temp(j,y));
                biomass_total_temp(y,j)=biomass_population(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                recruits_AM(j,r,y)=abundance_at_age_AM(j,r,y,a);                
                rec_index_AM(j,r,y)=recruits_AM(j,r,y)*mfexp(rec_index_RN(j,r,y)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));
                rec_index_AM_temp(j,y,r)=rec_index_AM(j,r,y);
                rec_index_prop_AM(j,r,y)=rec_index_AM(j,r,y)/sum(rec_index_AM_temp(j,y));

                abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                bio_res(j,r,y,a)=bio_move_temp(j,r);
                bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,r,y,a)-bio_res(j,r,y,a);
    
          } //end fleets loop

             for (int z=1;z<=nfleets_survey(j);z++)    /// survey index  1. Currently set up for more than 1 survey fleet
              {
               if(tsurvey(j,r)==0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                {
                  true_survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*q_survey(j,r,z);
                  true_survey_fleet_overlap_age_bio(p,j,r,y,z,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  true_survey_fleet_bio_overlap(p,j,r,y,z)=sum(true_survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  true_survey_fleet_bio_overlap_temp(j,r,y,z,p)=true_survey_fleet_bio_overlap(p,j,r,y,z);
                  OBS_survey_fleet_bio_overlap(p,j,r,y,z)=true_survey_fleet_bio_overlap(p,j,r,y,z)*mfexp(survey_RN_overlap(p,j,r,y,z)*sigma_survey_overlap(p,j,r,z)-.5*square(sigma_survey_overlap(p,j,r,z)));
                  OBS_survey_fleet_bio_temp(j,r,y,z,p)=OBS_survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  true_survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*q_survey(j,r,z);
                  true_survey_fleet_age_bio(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_age_bio(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=true_survey_fleet_bio(j,r,y,z)*mfexp(survey_RN(j,r,y,z)*sigma_survey(j,r,z)-.5*square(sigma_survey(j,r,z)));
                 }
                if(natal_homing_switch==1)
                 {
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_bio_overlap_temp(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=sum(OBS_survey_fleet_bio_temp(j,r,y,z));  
                 }
             
                  true_survey_region_bio_overlap(p,j,y,r)=sum(true_survey_fleet_bio_overlap(p,j,r,y));               
                  true_survey_population_bio_overlap(p,y,j)=sum(true_survey_region_bio_overlap(p,j,y));               
                  true_survey_natal_bio_overlap(y,p)=sum(true_survey_population_bio_overlap(p,y));               
                  true_survey_total_bio_overlap(y)=sum(true_survey_natal_bio_overlap(y));
                  OBS_survey_region_bio_overlap(p,j,y,r)=sum(OBS_survey_fleet_bio_overlap(p,j,r,y));
                  OBS_survey_population_bio_overlap(p,y,j)=sum(OBS_survey_region_bio_overlap(p,j,y));
                  OBS_survey_natal_bio_overlap(y,p)=sum(OBS_survey_population_bio_overlap(p,y));
                  OBS_survey_total_bio_overlap(y)=sum(OBS_survey_natal_bio_overlap(y));
                  
                         // calculate true survey abundance for tagging simulation **dh**
                  true_survey_fleet_age_temp(j,r,y,a,z)=true_survey_fleet_age(j,r,y,z,a);
                  true_survey_region_abundance(j,y,r,a)=sum(true_survey_fleet_age_temp(j,r,y,a));
                  true_survey_region_abundance_temp(j,y,a,r)=true_survey_region_abundance(j,y,r,a);
                  true_survey_population_abundance(y,j,a)=sum(true_survey_region_abundance_temp(j,y,a));
                  true_survey_population_abundance_temp(y,a,j)=true_survey_population_abundance(y,j,a);
                  true_survey_total_abundance(y,a)=sum(true_survey_population_abundance_temp(y,a));
                   // end abundance for tagging *dh
                  true_survey_region_bio(j,y,r)=sum(true_survey_fleet_bio(j,r,y));
                  true_survey_population_bio(y,j)=sum(true_survey_region_bio(j,y));
                  true_survey_total_bio(y)=sum(true_survey_population_bio(y));
                  OBS_survey_region_bio(j,y,r)=sum(OBS_survey_fleet_bio(j,r,y));
                  OBS_survey_population_bio(y,j)=sum(OBS_survey_region_bio(j,y));
                  OBS_survey_total_bio(y)=sum(OBS_survey_population_bio(y));

                  apport_region_survey_biomass(j,r,y)=OBS_survey_region_bio(j,y,r)/OBS_survey_population_bio(y,j);
                  
                }  //tsurvey==0
               } //end survey_fleets
            }
           }
          }
         } //end age loop

 ///////////////////////////////////////////////////////////////////////////////////
 ////////////////NEWTON-RAPHSON for YEAR 1////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                   if(model_type_switch>0)  //Find the fully selected F that would result in the TAC
                     {
                      for (int a=1;a<=nages;a++)
                       {
                         for(int x=1;x<=nfleets(j);x++)
                           {
                             if(model_type_switch==1)
                                {
                                 if(natal_homing_switch>0)
                                   {
                                               if(parse_TAC==0) //use input TAC
                                                {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x);
                                                   }
                                                   if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y);
                                                   }
                                                   if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y);
                                                   }
                                                }
                                             if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==2)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //apportion TAC equally among fleets if y<timelag
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_TAC(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                   }
                                                 if(parse_TAC_source==3)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  }
                                                 }
                                     if(TAC(j,r,x,y)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(TAC(j,r,x,y)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                             {
                                              if(s==1)
                                               {
                                                 fofFvect(s)=((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                                 fprimeFhigh(s)=(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                                 fprimeFlow(s)=(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));                                              
                                               }
                                              if(s>1)
                                               {
                                                fofFvect(s)=((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                                fprimeFhigh(s)=(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                                fprimeFlow(s)=(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                               }
                                             }
                                              fofF=sum(fofFvect)-TAC(j,r,x,y);
                                              fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                              Fnew=Fnew-(fofF/fprimeF);                                           
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
                                 if(natal_homing_switch==0)
                                   {
                                               if(parse_TAC==0) //use input TAC
                                                {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x);
                                                   }
                                                   if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y);
                                                   }
                                                   if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y);
                                                   }
                                                }
                                               if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==2)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //apportion TAC equally among fleets if y<timelag
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_TAC(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                   }
                                                 if(parse_TAC_source==3)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  }
                                                 }
                                     if(TAC(j,r,x,y)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(TAC(j,r,x,y)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFhigh(s)=weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFlow(s)=weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFhigh(s)=weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFlow(s)=weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              }
                                             }
                                            fofF=sum(fofFvect)-TAC(j,r,x,y);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                       } 
                                      }
                                     }

                             if(model_type_switch==2)
                                {
                                 if(natal_homing_switch>0)
                                   {
                                               if(parse_TAC==0) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 u(j,r,x)=input_u(j,r,x);
                                                }
                                               if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                   }
                                                 if(parse_TAC_source==2)
                                                  {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       u(j,r,x)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                 if(parse_TAC_source==3)
                                                  {
                                                   u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                  }
                                                 }
                                     if(u(j,r,x)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(u(j,r,x)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                             {
                                              if(s==1)
                                               {  
                                                 fofFvect(s)=(((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                 fprimeFhigh(s)=((((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                 fprimeFlow(s)=((((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM_overlap_region_all_natal(j,r,y);                                              
                                               }
                                              if(s>1)
                                               {
                                                fofFvect(s)=(((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                fprimeFhigh(s)=((((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                fprimeFlow(s)=((((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                               }
                                             } 
                                            fofF=sum(fofFvect)-u(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
                                 if(natal_homing_switch==0)
                                   {
                                               if(parse_TAC==0) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 u(j,r,x)=input_u(j,r,x);
                                                }
                                               if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                   }
                                                 if(parse_TAC_source==2)
                                                  {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       u(j,r,x)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                 if(parse_TAC_source==3)
                                                  {
                                                   u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                  }
                                                 }
                                     if(u(j,r,x)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(u(j,r,x)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=(weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFhigh(s)=(weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFlow(s)=(weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=(weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFhigh(s)=(weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFlow(s)=(weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              }
                                             }
                                            fofF=sum(fofFvect)-u(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                       } 
                                      }
                                     }
                             F_fleet(j,r,y,a,x)=Fnew*selectivity(j,r,y,a,x);                                         
                         }
                       F(j,r,y,a)=sum(F_fleet(j,r,y,a)); 
                      }                     
                     }
                    }
                   }
                  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         for (int a=1;a<=nages;a++)
           {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
                abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));
                SSB_region_temp_overlap(p,j,r,y,a)=abundance_spawn_overlap(p,j,r,y,a)*wt_mat_mult_reg(p,r,y,a); //changed mat by region
                SSB_region_overlap(p,j,r,y)=sum(SSB_region_temp_overlap(p,j,r,y));

                  SSB_overlap_natal=0;
                  if(natal_homing_switch==1 && spawn_return_switch==1)
                   {
                    for(int k=1;k<=npops;k++)
                     {
                      for (int n=1;n<=nregions(k);n++)
                       {
                        if(p==k && j==k)
                         {
                          SSB_overlap_natal(k,n)=0;  // already account for all fish already in natal population in calc below
                         }   
                        if(p==j && j!=k)
                         {
                          SSB_overlap_natal(k,n)=spawn_return_prob(p)*SSB_region_overlap(p,k,n,y);
                         } 
                       }
                      } 
                      if(p==j)
                      {
                       SSB_region_overlap(p,j,r,y)=SSB_region_overlap(p,j,r,y)+(sum(SSB_overlap_natal)/nregions(p));  //reutrning SSB is split evenly across natal regionp
                       }
                   }
                   
                SSB_population_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                SSB_population_overlap(p,j,y)=sum(SSB_population_temp_overlap(p,j,y));
                SSB_natal_overlap_temp(p,y,j)=SSB_population_overlap(p,j,y);
                SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));
                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_population(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                if(a==1)
                 {
                  catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F(j,r,y,a)+M(j,r,y,a))*(1-tspawn(p))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                  catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))*(1-tspawn(p))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a)));              
                 }
                if(a>1)
                 {
                  catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a))); //
                  catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a)));              
                 }
                yield_region_fleet_temp_overlap(p,j,r,z,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_fleet_overlap(p,j,r,z,y,a);
                yield_region_fleet_overlap(p,j,r,z,y)=sum(yield_region_fleet_temp_overlap(p,j,r,z,y));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_population_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_population_overlap(p,j,y,a)=sum(catch_at_age_population_temp_overlap(p,j,y,a));
                yield_population_temp_overlap(p,j,y,a)=weight_catch(p,r,y,a)*catch_at_age_population_overlap(p,j,y,a);
                yield_population_overlap(p,j,y)=sum(yield_population_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_population_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_population_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_fleet_bio_overlap(p,j,r,z,y)=yield_region_fleet_overlap(p,j,r,z,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_population_bio_overlap(p,j,y)=yield_population_overlap(p,j,y)/biomass_population_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_population_overlap(p,j,y)=biomass_population_overlap(p,j,y)/biomass_population_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
                Bratio_population_overlap(p,j,y)=SSB_population_overlap(p,j,y)/SSB_zero(p);
                Bratio_natal_overlap(p,y)=SSB_natal_overlap(p,y)/SSB_zero(p);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////Y==1 Abundance Calcs//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  CONTINUATION of NON-NATAL HOMING CALCS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                if(a==1)
                 {
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))*(1-tspawn(j))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));
                 }
                if(a>1)
                 {
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));
                 }
                yield_fleet_temp(j,r,y,z,a)=weight_catch(j,r,y,a)*catch_at_age_fleet(j,r,y,a,z);
                yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                yield_region_temp(j,r,y,a)=weight_catch(j,r,y,a)*catch_at_age_region(j,r,y,a);
                yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                catch_at_age_population_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                catch_at_age_population(j,y,a)=sum(catch_at_age_population_temp(j,y,a));
                yield_population_temp(j,y,a)=weight_catch(j,r,y,a)*catch_at_age_population(j,y,a);
                yield_population(j,y)=sum(yield_population_temp(j,y));
                catch_at_age_total_temp(y,a,j)=catch_at_age_population(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_population(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_population_num(j,y,a)=catch_at_age_population(j,y,a)/abundance_population(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_population_bio(j,y)=yield_population(j,y)/biomass_population(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_population(j,y)=biomass_population(j,y)/biomass_population(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);

              if(natal_homing_switch>0)
               {
               if(p==j)
               {
                SSB_region(j,r,y)=SSB_region_overlap(p,j,r,y);  //if natal homing only account for SSB that is in its natal populations area, don't sum across natal populations
               }
              }
              if(natal_homing_switch==0)
              {
                SSB_region_temp(j,r,y,a)=abundance_spawn(j,r,y,a)*wt_mat_mult_reg(j,r,y,a); // changed mat by region
                SSB_region(j,r,y)=sum(SSB_region_temp(j,r,y));
              }
                SSB_population_temp(j,y,r)=SSB_region(j,r,y); 
                SSB_population(j,y)=sum(SSB_population_temp(j,y)); 
                SSB_total_temp(y,j)=SSB_population(j,y);
                SSB_total(y)=sum(SSB_total_temp(y));
                Bratio_population(j,y)=SSB_population(j,y)/SSB_zero(j);
                Bratio_total(y)=SSB_total(y)/sum(SSB_zero);

                //YIELD observation error
                OBS_yield_region_fleet_overlap(p,j,r,y,z)=yield_region_fleet_overlap(p,j,r,z,y)*mfexp(yield_RN_overlap(p,j,r,y,z)*sigma_catch_overlap(p,j,r,z)-.5*square(sigma_catch_overlap(p,j,r,z)));
                OBS_yield_fleet_temp(j,r,y,z,p)=OBS_yield_region_fleet_overlap(p,j,r,y,z);
               if(natal_homing_switch==0)
                {
                 OBS_yield_fleet(j,r,y,z)=yield_fleet(j,r,y,z)*mfexp(yield_RN(j,r,y,z)*sigma_catch(j,r,z)-.5*square(sigma_catch(j,r,z)));
                }
               if(natal_homing_switch==1)
                {
                 OBS_yield_fleet(j,r,y,z)=sum(OBS_yield_fleet_temp(j,r,y,z));  
                }
                OBS_yield_region_overlap(p,j,y,r)=sum(OBS_yield_region_fleet_overlap(p,j,r,y));
                OBS_yield_population_overlap(p,y,j)=sum(OBS_yield_region_overlap(p,j,y));
                OBS_yield_natal_overlap(y,p)=sum(OBS_yield_population_overlap(p,y));
                OBS_yield_total_overlap(y)=sum(OBS_yield_natal_overlap(y));
                OBS_yield_region(j,y,r)=sum(OBS_yield_fleet(j,r,y));
                OBS_yield_population(y,j)=sum(OBS_yield_region(j,y));
                OBS_yield_total(y)=sum(OBS_yield_population(y));
             //apportion variables
                apport_yield_region(j,r,y)=OBS_yield_region(j,y,r)/OBS_yield_population(y,j);

          } //end fleets loop

             for (int z=1;z<=nfleets_survey(j);z++)    /// survey index  1. Currently set up for more than 1 survey fleet
              {
               if(tsurvey(j,r)>0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                {
                  true_survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  true_survey_fleet_overlap_age_bio(p,j,r,y,z,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  true_survey_fleet_bio_overlap(p,j,r,y,z)=sum(true_survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  true_survey_fleet_bio_overlap_temp(j,r,y,z,p)=true_survey_fleet_bio_overlap(p,j,r,y,z);
                  OBS_survey_fleet_bio_overlap(p,j,r,y,z)=true_survey_fleet_bio_overlap(p,j,r,y,z)*mfexp(survey_RN_overlap(p,j,r,y,z)*sigma_survey_overlap(p,j,r,z)-.5*square(sigma_survey_overlap(p,j,r,z)));
                  OBS_survey_fleet_bio_temp(j,r,y,z,p)=OBS_survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  true_survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  true_survey_fleet_age_bio(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_age_bio(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=true_survey_fleet_bio(j,r,y,z)*mfexp(survey_RN(j,r,y,z)*sigma_survey(j,r,z)-.5*square(sigma_survey(j,r,z)));
                 }


                if(natal_homing_switch==1)
                 {
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_bio_overlap_temp(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=sum(OBS_survey_fleet_bio_temp(j,r,y,z));  
                 }
             
                  true_survey_region_bio_overlap(p,j,y,r)=sum(true_survey_fleet_bio_overlap(p,j,r,y));               
                  true_survey_population_bio_overlap(p,y,j)=sum(true_survey_region_bio_overlap(p,j,y));               
                  true_survey_natal_bio_overlap(y,p)=sum(true_survey_population_bio_overlap(p,y));               
                  true_survey_total_bio_overlap(y)=sum(true_survey_natal_bio_overlap(y));
                  OBS_survey_region_bio_overlap(p,j,y,r)=sum(OBS_survey_fleet_bio_overlap(p,j,r,y));
                  OBS_survey_population_bio_overlap(p,y,j)=sum(OBS_survey_region_bio_overlap(p,j,y));
                  OBS_survey_natal_bio_overlap(y,p)=sum(OBS_survey_population_bio_overlap(p,y));
                  OBS_survey_total_bio_overlap(y)=sum(OBS_survey_natal_bio_overlap(y));
                  
                  
                   // calculate true survey abundance for tagging simulation **dh**
                 true_survey_fleet_age_temp(j,r,y,a,z)=true_survey_fleet_age(j,r,y,z,a);
                 true_survey_region_abundance(j,y,r,a)=sum(true_survey_fleet_age_temp(j,r,y,a));
                 true_survey_region_abundance_temp(j,y,a,r)=true_survey_region_abundance(j,y,r,a);
                 true_survey_population_abundance(y,j,a)=sum(true_survey_region_abundance_temp(j,y,a));
                 true_survey_population_abundance_temp(y,a,j)=true_survey_population_abundance(y,j,a);
                 true_survey_total_abundance(y,a)=sum(true_survey_population_abundance_temp(y,a));

             // end abundance for tagging *dh

                  true_survey_region_bio(j,y,r)=sum(true_survey_fleet_bio(j,r,y));
                  true_survey_population_bio(y,j)=sum(true_survey_region_bio(j,y));
                  true_survey_total_bio(y)=sum(true_survey_population_bio(y));
                  OBS_survey_region_bio(j,y,r)=sum(OBS_survey_fleet_bio(j,r,y));
                  OBS_survey_population_bio(y,j)=sum(OBS_survey_region_bio(j,y));
                  OBS_survey_total_bio(y)=sum(OBS_survey_population_bio(y));
                

                  apport_region_survey_biomass(j,r,y)=OBS_survey_region_bio(j,y,r)/OBS_survey_population_bio(y,j);


         
                } //tsurvey>0
               }  //end survey_fleets
             }
            }
           }
          }//end age loop
        } //end yr 1

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////Recruitment Calcs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////need to fix and use SSB and R_ave for BH calcs.

    if(y>1)
     {
        for (int a=1;a<=nages;a++)
          {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
                   if(a==1)
                   {
                 if(natal_homing_switch>0)
                 {
                 if(p==j)
                 {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a population
                    {
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }

                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a population
                    {
                    recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a population
                    {
                     recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1)));
                    }
                  }
                  
                if(Rec_type==3) //environmental recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a population
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }
                 }
                 }

             if(natal_homing_switch==0)
              {
                if(Rec_type==1) //average recruitment
                  {
                  if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                    }
                    if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                    }
                   }

                if(Rec_type==3) //average recruitment
                  {
                  if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=env_rec(y)* rec_devs(j,y)*Rec_Prop(j,r,y-1);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                   {
                    recruits_BM(j,r,y)=env_rec(y)* R_ave(j)*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                   }
                 }
                 }
               rec_index_BM(j,r,y)=recruits_BM(j,r,y)*mfexp(rec_index_RN(j,r,y)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));
               rec_index_BM_temp(j,y,r)=rec_index_BM(j,r,y);
               rec_index_prop_BM(j,r,y)=rec_index_BM(j,r,y)/sum(rec_index_BM_temp(j,y));

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////a==1 overlap calcs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                 if(p==j)
                 {
                   abundance_at_age_BM_overlap_region(p,j,y,a,r)=recruits_BM(j,r,y);
                 }
                 if(p!=j)
                 {
                   abundance_at_age_BM_overlap_region(p,j,y,a,r)=0;
                 }
                abundance_at_age_BM_overlap_population(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));

                biomass_BM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////// metapop calcs /////////////////////////////////////////////////////////////////////////

                 abundance_at_age_BM(j,r,y,a)=recruits_BM(j,r,y);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

                  if(natal_homing_switch==0)
                   {
                    biomass_BM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
                  if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {
                  if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {
                    biomass_BM_overlap_temp(j,r,y,a,p)=biomass_BM_age_overlap(p,j,r,y,a);
                    biomass_BM_age(j,r,y,a)=sum(biomass_BM_overlap_temp(j,r,y,a));
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
                   }
              }
             }
            }
           }
          } //close loops so have full biomass vectors filled in at start of DD movement calcs
   
        if(a==1)
         {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Density-Dependent Movement Calcs///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8 && DD_move_age_switch==1) //DD movement is age based (based on age based biomass not total bio)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM_age(j,r,y,a)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }       
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n,a);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r,a)); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                          if(a>1 || (a==1 && larval_move_switch!=0))
                           {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                           }
                            if(a==1 && larval_move_switch==0)
                             {
                              T(j,r,y,a,s,n)=T(j,r,y,a,s,n);
                             }
                          }
                        }
                    }

                    if(move_switch==8 && DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_AM(j,r,y-1)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,1)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }       
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_AM(s,n,y-1)/sum(Bstar(s,n));
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_AM(j,r,y-1)/sum(Bstar(j,r))); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                          if(a>1 || (a==1 && larval_move_switch!=0))
                           {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                           }
                            if(a==1 && larval_move_switch==0)
                             {
                              T(j,r,y,a,s,n)=T(j,r,y,a,s,n);
                             }
                          }
                        }
                    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_move_overlap_temp=0;
               
                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {

  ///////////NEED TO CHECK WHY USING FULL RECRUITMENT FUNCTIOnp TO CALCULATE MOVEMENT OF AGE-1 FISH HERE, MIGHT BE ABLE TO CONDEnpE CODE AND USE RECRUITS_BM OR ABUNDANCE_BM


                 if(natal_homing_switch>0)
                 {
                 if(p==k)
                 {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }

                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_population_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_population_overlap(p,k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_population_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_population_overlap(p,k,y-1)))*rec_devs(k,y)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r);
                    }
                  }

                if(Rec_type==3) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }
                 }
                 }
                 
             if(natal_homing_switch==0)
              {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                    }
                 if(Rec_type==2) //BH recruitment
                  {
                 if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    abundance_move_overlap_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }

                if(Rec_type==3) //environmental recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }

                }
                   }
                  }

                abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_population(p,r,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_population(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));
                biomass_AM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_population_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_population_overlap(p,j,y)=sum(biomass_population_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_population_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////y>1, A==1 METAPOP TYPE CALCS (MOVEMENT)//////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

                    abundance_move_temp=0;
                    bio_move_temp=0;

        for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
             if(natal_homing_switch==0)
              {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(k,n,y,a,j,r);
                     }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    abundance_move_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                if(Rec_type==3) //env recruitment
                  {
                  if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y)*Rec_Prop(k,n,y-1)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                }
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_population(k,n,y,a);
              }
             }
                  
                  if(natal_homing_switch>0)
                   {                  
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)=biomass_AM_overlap_region_all_natal(j,r,y);                   
                   }
                  if(natal_homing_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_AM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }

               recruits_AM(j,r,y)=abundance_at_age_AM(j,r,y,a);
               rec_index_AM(j,r,y)=recruits_AM(j,r,y)*mfexp(rec_index_RN(j,r,y)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));
               rec_index_AM_temp(j,y,r)=rec_index_AM(j,r,y);
               rec_index_prop_AM(j,r,y)=rec_index_AM(j,r,y)/sum(rec_index_AM_temp(j,y));

                biomass_population_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_population(j,y)=sum(biomass_population_temp(j,y));
                biomass_total_temp(y,j)=biomass_population(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

               abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
               abundance_res(j,r,y,a)=abundance_move_temp(j,r);
               abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
               bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
               bio_res(j,r,y,a)=bio_move_temp(j,r);
               bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,r,y,a)-bio_res(j,r,y,a);
           }
          }
         }
        }
       } //end a==1 if statement

 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       if(a==2) //account for partial year mortality during spawning year
        {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
                  abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p)));
                  abundance_at_age_BM_overlap_population(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                  biomass_BM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                  biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));

/////////////////////////metapop/////////////////////////////////////////////////////////
            abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(j))); //account for time of spawning (i.e., if born midway only experience a half year of mortality from age-1 to age-2)
/////////////////////////////////////////////////////////////////////////////////////////

                  if(natal_homing_switch==0)
                   {
                    biomass_BM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
                  if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {
                    biomass_BM_overlap_temp(j,r,y,a,p)=biomass_BM_age_overlap(p,j,r,y,a);
                    biomass_BM_age(j,r,y,a)=sum(biomass_BM_overlap_temp(j,r,y,a));
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
              }
             }
            }
           }
          } //close loops so have full biomass vectors filled in at start of DD movement calcs
          
         if(a==2)
          { 
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Density-Dependent Movement Calcs///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8 && DD_move_age_switch==1) //DD movement is age based (based on age based biomass not total bio)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM_age(j,r,y,a)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }       
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n,a);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r,a)); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                          }
                        }
                    }

                    if(move_switch==8 && DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_AM(j,r,y-1)));
                       if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,1)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }      
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_AM(s,n,y-1)/sum(Bstar(s,n));
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_AM(j,r,y-1)/sum(Bstar(j,r))); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                          }
                        }
                    }
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7  ||a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(p)))*T(p,n,y,a,j,r); //with overlap always use natal population movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal population movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(p)))*return_probability(p); //with overlap always use natal population movement rates
                      }
                     }
                   }
                  }
                    if(move_switch!=6 || move_switch!=7  || a==1)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
                     }
                    if(move_switch==7)  //all fish stay where they were (i.e., no return migration)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p)));
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a<return_age || a>return_age)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p))); //with overlap always use natal population movement rates                     
                       }
                      if(a==return_age && p==j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p)))+(sum(abundance_move_overlap_temp)/nregions(p));
                       }
                      if(a==return_age && p!=j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=(1-return_probability(p))*abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p)));
                       }
                      }
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_population(p,r,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_population(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));
                biomass_AM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_population_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_population_overlap(p,j,y)=sum(biomass_population_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_population_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

 //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
   ///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    

                    abundance_move_temp=0;
                     bio_move_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=abundance_at_age_AM(k,n,y-1,a-1)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(k)))*T(k,n,y,a,j,r);
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_population(k,n,y,a);
                   }
                  }

                  if(natal_homing_switch>0)
                   {                  
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);                  
                   }
                  if(natal_homing_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_AM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }

                biomass_population_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_population(j,y)=sum(biomass_population_temp(j,y));
                biomass_total_temp(y,j)=biomass_population(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                   abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                   abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                   abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                   bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                   bio_res(j,r,y,a)=bio_move_temp(j,r);
                   bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,r,y,a)-bio_res(j,r,y,a);
             }
            }
           }
          }
         }//end a==2 if statement

 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       if(a>2 && a<nages)
        {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
                  abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
                  abundance_at_age_BM_overlap_population(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                  biomass_BM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                  biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));

/////////////////////////metapop/////////////////////////////////////////////////////////
                 abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
/////////////////////////////////////////////////////////////////////////////////////////

                  if(natal_homing_switch==0)
                   {
                    biomass_BM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
                  if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {
                    biomass_BM_overlap_temp(j,r,y,a,p)=biomass_BM_age_overlap(p,j,r,y,a);
                    biomass_BM_age(j,r,y,a)=sum(biomass_BM_overlap_temp(j,r,y,a));
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
              }
             }
            }
           }
          } //close loops so have full biomass vectors filled in at start of DD movement calcs


       if(a>2 && a<nages)
        {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Density-Dependent Movement Calcs///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
            
                    if(move_switch==8 && DD_move_age_switch==1) //DD movement is age based (based on age based biomass not total bio)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM_age(j,r,y,a)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }       
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n,a);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r,a)); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                          }
                        }
                    }

                    if(move_switch==8 && DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_AM(j,r,y-1)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,1)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }       
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_AM(s,n,y-1)/sum(Bstar(s,n));
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_AM(j,r,y-1)/sum(Bstar(j,r))); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                          }
                        }
                    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7 || a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(p,n,y,a,j,r); //with overlap always use natal population movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal population movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*return_probability(p); //with overlap always use natal population movement rates
                      }
                     }
                   }
                  }
                    if(move_switch!=6 || move_switch!=7  || a==1)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
                     }
                    if(move_switch==7)  //all fish stay where they were (i.e., no return migration)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a<return_age || a>return_age)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))); //with overlap always use natal population movement rates                     
                       }
                      if(a==return_age && p==j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+(sum(abundance_move_overlap_temp)/nregions(p));
                       }
                      if(a==return_age && p!=j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=(1-return_probability(p))*abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
                       }
                      }
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_population(p,r,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_population(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));
                biomass_AM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_population_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_population_overlap(p,j,y)=sum(biomass_population_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_population_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));
                
 //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
   ///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    

                   abundance_move_temp=0;
                   bio_move_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=abundance_at_age_AM(k,n,y-1,a-1)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(k,n,y,a,j,r);
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_population(k,n,y,a);
                   }
                  }
                  
                  if(natal_homing_switch>0)
                   {
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);
                   }
                  if(natal_homing_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_AM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }

                biomass_population_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_population(j,y)=sum(biomass_population_temp(j,y));
                biomass_total_temp(y,j)=biomass_population(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                bio_res(j,r,y,a)=bio_move_temp(j,r);
                bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,r,y,a)-bio_res(j,r,y,a);
           }
          }
         }
        }
       }//end a>2 <nages if statement

 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       if(a==nages) //account for fish already in plus group
        {
          for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
                  abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
                  abundance_at_age_BM_overlap_population(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                  biomass_BM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                  biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                
/////////////////////////metapop/////////////////////////////////////////////////////////
               abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM(j,r,y-1,a)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
/////////////////////////////////////////////////////////////////////////////////////////

                  if(natal_homing_switch==0)
                   {
                    biomass_BM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
                  if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {
                    biomass_BM_overlap_temp(j,r,y,a,p)=biomass_BM_age_overlap(p,j,r,y,a);
                    biomass_BM_age(j,r,y,a)=sum(biomass_BM_overlap_temp(j,r,y,a));
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                   }
              }
             }
            }
           }
          } //close loops so have full biomass vectors filled in at start of DD movement calcs
          
       if(a==nages) //account for fish already in plus group
        {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Density-Dependent Movement Calcs///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8 && DD_move_age_switch==1) //DD movement is age based (based on age based biomass not total bio)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM_age(j,r,y,a)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }       
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n,a);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r,a)); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                          }
                        }
                    }

                    if(move_switch==8 && DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_AM(j,r,y-1)));
                      if(rand_move==1)
                       {
                        Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,1)*sigma_T(j)-0.5*square(sigma_T(j)));
                       }
                       if(Fract_Move_DD(j,r,y,a)>1)
                        {
                         Fract_Move_DD(j,r,y,a)=1;
                        }
                       if(Fract_Move_DD(j,r,y,a)<0)
                        {
                         Fract_Move_DD(j,r,y,a)=0;
                        }       
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_AM(s,n,y-1)/sum(Bstar(s,n));
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_AM(j,r,y-1)/sum(Bstar(j,r))); //do not include population/region moving from in summation
                        for (int s=1;s<=npops;s++)
                         {
                          for (int n=1;n<=nregions(s);n++)
                           {
                            if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                             {
                              rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                             }
                            if(j!=s || r!=n) 
                             {
                              rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                             }
                            if(sum(nregions)==2)
                             {
                              rel_bio(j,r,y,a,s,n)=1;
                             }
                            }
                           }                          
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {
                           if(j==s && r==n) 
                            {
                             T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                            }
                           if(j!=s || r!=n) 
                            {
                             T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                            }
                          }
                        }
                    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7  || a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(p,n,y,a,j,r)+abundance_at_age_AM_overlap_region(p,k,y-1,a,n)*mfexp(-(M(k,n,y-1,a)+F(k,n,y-1,a)))*T(p,n,y,a,j,r); //with overlap always use natal population movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal population movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*return_probability(p)+abundance_at_age_AM_overlap_region(p,k,y-1,a,n)*mfexp(-(M(k,n,y-1,a)+F(k,n,y-1,a)))*return_probability(p); //with overlap always use natal population movement rates
                      }
                     }
                   }
                  }
                    if(move_switch!=6 || move_switch!=7  || a==1)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
                     }
                    if(move_switch==7)  //all fish stay where they were (i.e., no return migration)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a<return_age || a>return_age)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a))); //with overlap always use natal population movement rates                     
                       }
                      if(a==return_age && p==j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)))+(sum(abundance_move_overlap_temp)/nregions(p));
                       }
                      if(a==return_age && p!=j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=(1-return_probability(p))*abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+(1-return_probability(p))*abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
                       }
                      }
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_population(p,r,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_population(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));
                biomass_AM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_population_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_population_overlap(p,j,y)=sum(biomass_population_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_population_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

 //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
   ///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
                  
                   abundance_move_temp=0;
                   bio_move_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=abundance_at_age_AM(k,n,y-1,a-1)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(k,n,y,a,j,r)+abundance_at_age_AM(k,n,y-1,a)*mfexp(-(M(k,n,y-1,a)+F(k,n,y-1,a)))*T(k,n,y,a,j,r);;
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_population(k,n,y,a);
                   }
                  }
                  if(natal_homing_switch>0)
                   {
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);
                   }
                  if(natal_homing_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_AM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }

                biomass_population_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_population(j,y)=sum(biomass_population_temp(j,y));
                biomass_total_temp(y,j)=biomass_population(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                   abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                   abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                   abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                   bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                   bio_res(j,r,y,a)=bio_move_temp(j,r);
                   bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,r,y,a)-bio_res(j,r,y,a);
            }
           }
          }
         }
        } //end nages if statement

           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets_survey(j);z++)    /// survey index  1. Currently set up for more than 1 survey fleet
                  {
                   if(tsurvey(j,r)==0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                   {
                  true_survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*q_survey(j,r,z);
                  true_survey_fleet_overlap_age_bio(p,j,r,y,z,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  true_survey_fleet_bio_overlap(p,j,r,y,z)=sum(true_survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  true_survey_fleet_bio_overlap_temp(j,r,y,z,p)=true_survey_fleet_bio_overlap(p,j,r,y,z);
                  OBS_survey_fleet_bio_overlap(p,j,r,y,z)=true_survey_fleet_bio_overlap(p,j,r,y,z)*mfexp(survey_RN_overlap(p,j,r,y,z)*sigma_survey_overlap(p,j,r,z)-.5*square(sigma_survey_overlap(p,j,r,z)));
                  OBS_survey_fleet_bio_temp(j,r,y,z,p)=OBS_survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  true_survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*q_survey(j,r,z);
                  true_survey_fleet_age_bio(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_age_bio(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=true_survey_fleet_bio(j,r,y,z)*mfexp(survey_RN(j,r,y,z)*sigma_survey(j,r,z)-.5*square(sigma_survey(j,r,z)));
                 }
                if(natal_homing_switch==1)
                 {
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_bio_overlap_temp(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=sum(OBS_survey_fleet_bio_temp(j,r,y,z));  
                 }
             
                  true_survey_region_bio_overlap(p,j,y,r)=sum(true_survey_fleet_bio_overlap(p,j,r,y));               
                  true_survey_population_bio_overlap(p,y,j)=sum(true_survey_region_bio_overlap(p,j,y));               
                  true_survey_natal_bio_overlap(y,p)=sum(true_survey_population_bio_overlap(p,y));               
                  true_survey_total_bio_overlap(y)=sum(true_survey_natal_bio_overlap(y));
                  OBS_survey_region_bio_overlap(p,j,y,r)=sum(OBS_survey_fleet_bio_overlap(p,j,r,y));
                  OBS_survey_population_bio_overlap(p,y,j)=sum(OBS_survey_region_bio_overlap(p,j,y));
                  OBS_survey_natal_bio_overlap(y,p)=sum(OBS_survey_population_bio_overlap(p,y));
                  OBS_survey_total_bio_overlap(y)=sum(OBS_survey_natal_bio_overlap(y));
                  
                   // calculate true survey abundance for tagging simulation **dh**
                 true_survey_fleet_age_temp(j,r,y,a,z)=true_survey_fleet_age(j,r,y,z,a);
                 true_survey_region_abundance(j,y,r,a)=sum(true_survey_fleet_age_temp(j,r,y,a));
                 true_survey_region_abundance_temp(j,y,a,r)=true_survey_region_abundance(j,y,r,a);
                 true_survey_population_abundance(y,j,a)=sum(true_survey_region_abundance_temp(j,y,a));
                 true_survey_population_abundance_temp(y,a,j)=true_survey_population_abundance(y,j,a);
                 true_survey_total_abundance(y,a)=sum(true_survey_population_abundance_temp(y,a));

             // end abundance for tagging *dh



                  true_survey_region_bio(j,y,r)=sum(true_survey_fleet_bio(j,r,y));
                  true_survey_population_bio(y,j)=sum(true_survey_region_bio(j,y));
                  true_survey_total_bio(y)=sum(true_survey_population_bio(y));
                  OBS_survey_region_bio(j,y,r)=sum(OBS_survey_fleet_bio(j,r,y));
                  OBS_survey_population_bio(y,j)=sum(OBS_survey_region_bio(j,y));
                  OBS_survey_total_bio(y)=sum(OBS_survey_population_bio(y));

                  apport_region_survey_biomass(j,r,y)=OBS_survey_region_bio(j,y,r)/OBS_survey_population_bio(y,j);
                  
                } //tsurvey==0
               } //end survey_fleets 
      }
     }
    }
   } //end age loop
 ///////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////
 /////////////////NEWTON RAPHSON CALCS//////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                   if(model_type_switch>0)  //Find the fully selected F that would result in the TAC
                     {
                      for (int a=1;a<=nages;a++)
                       {
                         for(int x=1;x<=nfleets(j);x++)
                           {
                             if(model_type_switch==1)
                                {
                                 if(natal_homing_switch>0)
                                   {
                                               if(parse_TAC==0) //use input TAC
                                                {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x);
                                                   }
                                                   if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y);
                                                   }
                                                   if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y);
                                                   }
                                                }
                                             if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==2)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //apportion TAC equally among fleets if y<timelag
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_TAC(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                   }
                                                 if(parse_TAC_source==3)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  }
                                                 }
                                     if(TAC(j,r,x,y)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(TAC(j,r,x,y)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                             {
                                              if(s==1)
                                               {
                                                 fofFvect(s)=((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                                 fprimeFhigh(s)=(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                                 fprimeFlow(s)=(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));                                              
                                               }
                                              if(s>1)
                                               {
                                                fofFvect(s)=((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                                fprimeFhigh(s)=(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                                fprimeFlow(s)=(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                               }
                                             }
                                              fofF=sum(fofFvect)-TAC(j,r,x,y);
                                              fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                              Fnew=Fnew-(fofF/fprimeF);                                           
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
                                 if(natal_homing_switch==0)
                                   {
                                               if(parse_TAC==0) //use input TAC
                                                {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x);
                                                   }
                                                   if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y);
                                                   }
                                                   if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y);
                                                   }
                                                }
                                             if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(y==1)
                                                    {
                                                     TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                                                    }
                                                   if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                                                    {
                                                     TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                                                    }
                                                   }
                                                  }
                                                 if(parse_TAC_source==2)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //apportion TAC equally among fleets if y<timelag
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_TAC(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_TAC(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                   }
                                                 if(parse_TAC_source==3)
                                                  {
                                                  if(calc_TAC_from_uMSY==0)
                                                   {
                                                    TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==1)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  if(calc_TAC_from_uMSY==2)
                                                   {
                                                    TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                                                   }
                                                  }
                                                 }
                                     if(TAC(j,r,x,y)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(TAC(j,r,x,y)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFhigh(s)=weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFlow(s)=weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFhigh(s)=weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFlow(s)=weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              }
                                             }
                                            fofF=sum(fofFvect)-TAC(j,r,x,y);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                       } 
                                      }
                                     }

                             if(model_type_switch==2)
                                {
                                 if(natal_homing_switch>0)
                                   {
                                               if(parse_TAC==0) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 u(j,r,x)=input_u(j,r,x);
                                                }
                                               if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                   }
                                                 if(parse_TAC_source==2)
                                                  {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       u(j,r,x)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                 if(parse_TAC_source==3)
                                                  {
                                                   u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                  }
                                                 }
                                     if(u(j,r,x)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(u(j,r,x)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                             {
                                              if(s==1)
                                               {  
                                                 fofFvect(s)=(((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                 fprimeFhigh(s)=((((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                 fprimeFlow(s)=((((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM_overlap_region_all_natal(j,r,y);                                              
                                               }
                                              if(s>1)
                                               {
                                                fofFvect(s)=(((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                fprimeFhigh(s)=((((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                                fprimeFlow(s)=((((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                                               }
                                             } 
                                            fofF=sum(fofFvect)-u(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
                                 if(natal_homing_switch==0)
                                   {
                                               if(parse_TAC==0) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 u(j,r,x)=input_u(j,r,x);
                                                }
                                               if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                                                {
                                                 if(parse_TAC_source==0)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                  }
                                                 if(parse_TAC_source==1)
                                                  {
                                                   if(y==1)
                                                    {
                                                     u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                    }
                                                   if(y>1)
                                                    {
                                                     u(j,r,x)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                    }
                                                   }
                                                 if(parse_TAC_source==2)
                                                  {
                                                   if(TAC_survey_parse_timelag_switch==1) //use timelag
                                                    {
                                                       if(y<=TAC_survey_parse_timelag) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>TAC_survey_parse_timelag)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)/nfleets(j);
                                                        }
                                                     }
                                                  if(TAC_survey_parse_timelag_switch==0) //no timelag
                                                    {
                                                     if(tsurvey(j,r)==0)
                                                      {
                                                       u(j,r,x)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)/nfleets(j);
                                                      }
                                                     if(tsurvey(j,r)>0)
                                                      {
                                                       if(y==1) //first year apportion TAC equally among fleets
                                                        {                                                     
                                                         u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                        }
                                                       if(y>1)
                                                        {                                                     
                                                         u(j,r,x)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                                                        }
                                                       }
                                                     }
                                                    }
                                                 if(parse_TAC_source==3)
                                                  {
                                                   u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                                                  }
                                                 }
                                     if(u(j,r,x)==0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                       {
                                        Fnew=0;
                                       }
                                     if(u(j,r,x)>0) //iterationp have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterationp
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=(weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFhigh(s)=(weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFlow(s)=(weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=(weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFhigh(s)=(weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFlow(s)=(weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              }
                                             }
                                            fofF=sum(fofFvect)-u(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                       } 
                                      }
                                     }
                             F_fleet(j,r,y,a,x)=Fnew*selectivity(j,r,y,a,x);                                         
                         }
                       F(j,r,y,a)=sum(F_fleet(j,r,y,a)); 
                      }                     
                     }
                    }
                   }
                  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //a==1 natal homing
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       for (int a=1;a<=nages;a++)
         {
           for (int p=1;p<=npops;p++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
          if(a==1)
            {

                abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));
                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_population(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))*(1-tspawn(p))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a)));              
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))*(1-tspawn(j))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                yield_region_fleet_temp_overlap(p,j,r,z,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_fleet_overlap(p,j,r,z,y,a);
                yield_region_fleet_overlap(p,j,r,z,y)=sum(yield_region_fleet_temp_overlap(p,j,r,z,y));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_population_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_population_overlap(p,j,y,a)=sum(catch_at_age_population_temp_overlap(p,j,y,a));
                yield_population_temp_overlap(p,j,y,a)=weight_catch(p,r,y,a)*catch_at_age_population_overlap(p,j,y,a);
                yield_population_overlap(p,j,y)=sum(yield_population_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_population_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_population_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_fleet_bio_overlap(p,j,r,z,y)=yield_region_fleet_overlap(p,j,r,z,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_population_bio_overlap(p,j,y)=yield_population_overlap(p,j,y)/biomass_population_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_population_overlap(p,j,y)=biomass_population_overlap(p,j,y)/biomass_population_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
               
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////a==1 metapop type abundance calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
                   abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                   catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))*(1-tspawn(j))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a)); //account for time of spawning in catch (tspawn divides out in F/(F+M

                   //catchsum(j,r,y,z)=sum(catch_at_age_fleet(j,r,y,z));
                   
                   yield_fleet_temp(j,r,y,z,a)=weight_catch(j,r,y,a)*catch_at_age_fleet(j,r,y,a,z);
                   yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                   catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                   yield_region_temp(j,r,y,a)=weight_catch(j,r,y,a)*catch_at_age_region(j,r,y,a);
                   yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                   catch_at_age_population_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                   catch_at_age_population(j,y,a)=sum(catch_at_age_population_temp(j,y,a));
                   yield_population_temp(j,y,a)=weight_catch(j,r,y,a)*catch_at_age_population(j,y,a);
                   yield_population(j,y)=sum(yield_population_temp(j,y));

                abundance_population_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_population(j,y,a)=sum(abundance_population_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_population(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_population(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_population(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_population_num(j,y,a)=catch_at_age_population(j,y,a)/abundance_population(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_population_bio(j,y)=yield_population(j,y)/biomass_population(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_population(j,y)=biomass_population(j,y)/biomass_population(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);
            } //end a==1 loop


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////a==2 overlap calcs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if(a==2)
    {
                abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));
                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_population(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a)));              
                yield_region_fleet_temp_overlap(p,j,r,z,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_fleet_overlap(p,j,r,z,y,a);
                yield_region_fleet_overlap(p,j,r,z,y)=sum(yield_region_fleet_temp_overlap(p,j,r,z,y));
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));              
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_population_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_population_overlap(p,j,y,a)=sum(catch_at_age_population_temp_overlap(p,j,y,a));
                yield_population_temp_overlap(p,j,y,a)=weight_catch(p,r,y,a)*catch_at_age_population_overlap(p,j,y,a);
                yield_population_overlap(p,j,y)=sum(yield_population_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_population_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_population_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_fleet_bio_overlap(p,j,r,z,y)=yield_region_fleet_overlap(p,j,r,z,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_population_bio_overlap(p,j,y)=yield_population_overlap(p,j,y)/biomass_population_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_population_overlap(p,j,y)=biomass_population_overlap(p,j,y)/biomass_population_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
                            
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////a==2 metapop type abundance calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                   
                  abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));                  
                  yield_fleet_temp(j,r,y,z,a)=weight_catch(j,r,y,a)*catch_at_age_fleet(j,r,y,a,z);
                  yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                  catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                  yield_region_temp(j,r,y,a)=weight_catch(j,r,y,a)*catch_at_age_region(j,r,y,a);
                  yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                  catch_at_age_population_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                  catch_at_age_population(j,y,a)=sum(catch_at_age_population_temp(j,y,a));
                  yield_population_temp(j,y,a)=weight_catch(j,r,y,a)*catch_at_age_population(j,y,a);
                  yield_population(j,y)=sum(yield_population_temp(j,y));

                abundance_population_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_population(j,y,a)=sum(abundance_population_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_population(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_population(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_population(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_population_num(j,y,a)=catch_at_age_population(j,y,a)/abundance_population(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_population_bio(j,y)=yield_population(j,y)/biomass_population(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_population(j,y)=biomass_population(j,y)/biomass_population(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);
          } //end a==2 loop
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////a>2 a<nages overlap calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
             if(a>2 && a<nages)
               {
                abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));
                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_population(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a)));              
                yield_region_fleet_temp_overlap(p,j,r,z,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_fleet_overlap(p,j,r,z,y,a);
                yield_region_fleet_overlap(p,j,r,z,y)=sum(yield_region_fleet_temp_overlap(p,j,r,z,y));
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_population_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_population_overlap(p,j,y,a)=sum(catch_at_age_population_temp_overlap(p,j,y,a));
                yield_population_temp_overlap(p,j,y,a)=weight_catch(p,r,y,a)*catch_at_age_population_overlap(p,j,y,a);
                yield_population_overlap(p,j,y)=sum(yield_population_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_population_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_population_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_fleet_bio_overlap(p,j,r,z,y)=yield_region_fleet_overlap(p,j,r,z,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_population_bio_overlap(p,j,y)=yield_population_overlap(p,j,y)/biomass_population_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_population_overlap(p,j,y)=biomass_population_overlap(p,j,y)/biomass_population_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
                                
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////a>2 a<nages metapop type abundance calcs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
                  abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));
                  yield_fleet_temp(j,r,y,z,a)=weight_catch(j,r,y,a)*catch_at_age_fleet(j,r,y,a,z);
                  yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                  catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                  yield_region_temp(j,r,y,a)=weight_catch(j,r,y,a)*catch_at_age_region(j,r,y,a);
                  yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                  catch_at_age_population_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                  catch_at_age_population(j,y,a)=sum(catch_at_age_population_temp(j,y,a));
                  yield_population_temp(j,y,a)=weight_catch(j,r,y,a)*catch_at_age_population(j,y,a);
                  yield_population(j,y)=sum(yield_population_temp(j,y));

                abundance_population_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_population(j,y,a)=sum(abundance_population_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_population(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_population(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_population(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_population_num(j,y,a)=catch_at_age_population(j,y,a)/abundance_population(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_population_bio(j,y)=yield_population(j,y)/biomass_population(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_population(j,y)=biomass_population(j,y)/biomass_population(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);
       } //end a>2 <nages loop

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////a==nages overlap calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           if(a==nages) //account for fish already in plus group
            {
                abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));
                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_population(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a)));              
                yield_region_fleet_temp_overlap(p,j,r,z,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_fleet_overlap(p,j,r,z,y,a);
                yield_region_fleet_overlap(p,j,r,z,y)=sum(yield_region_fleet_temp_overlap(p,j,r,z,y));
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_population_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_population_overlap(p,j,y,a)=sum(catch_at_age_population_temp_overlap(p,j,y,a));
                yield_population_temp_overlap(p,j,y,a)=weight_catch(p,r,y,a)*catch_at_age_population_overlap(p,j,y,a);
                yield_population_overlap(p,j,y)=sum(yield_population_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_population_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_population_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_fleet_bio_overlap(p,j,r,z,y)=yield_region_fleet_overlap(p,j,r,z,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_population_bio_overlap(p,j,y)=yield_population_overlap(p,j,y)/biomass_population_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_population_overlap(p,j,y)=biomass_population_overlap(p,j,y)/biomass_population_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);               
                 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////a==nages metapop type abundance calcs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
                  abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));
                  yield_fleet_temp(j,r,y,z,a)=weight_catch(j,r,y,a)*catch_at_age_fleet(j,r,y,a,z);
                  yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                  catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                  yield_region_temp(j,r,y,a)=weight_catch(j,r,y,a)*catch_at_age_region(j,r,y,a);
                  yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                  catch_at_age_population_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                  catch_at_age_population(j,y,a)=sum(catch_at_age_population_temp(j,y,a));
                  yield_population_temp(j,y,a)=weight_catch(j,r,y,a)*catch_at_age_population(j,y,a);
                  yield_population(j,y)=sum(yield_population_temp(j,y));

                abundance_population_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_population(j,y,a)=sum(abundance_population_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_population(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_population(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_population(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_population_num(j,y,a)=catch_at_age_population(j,y,a)/abundance_population(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_population_bio(j,y)=yield_population(j,y)/biomass_population(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_population(j,y)=biomass_population(j,y)/biomass_population(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);
        } //end nages if statement

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////SSB calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   //don't adjust SSB for fish not in natal population area because SR calculationp do this automatically by only using
   //SSB that is in natal population area (i.e., p==j)
   //for spawning migrationp scenarios the SSB is adjusted outside the loop to remove SSB of fish that actually returned
   // to natal population (i.e., remove this SSB from the non-natal areas...doesn't impact SR calcs so can do this outside
   //loops without conpequence to model
   
                SSB_region_temp_overlap(p,j,r,y,a)=abundance_spawn_overlap(p,j,r,y,a)*wt_mat_mult_reg(p,r,y,a); //added region
                SSB_region_overlap(p,j,r,y)=sum(SSB_region_temp_overlap(p,j,r,y));
                 SSB_overlap_natal=0;
                  if(natal_homing_switch==1 && spawn_return_switch==1) //spawning return calculationp
                   {
                    for(int k=1;k<=npops;k++)
                     {
                      for (int n=1;n<=nregions(k);n++)
                       {
                        if(p==k && j==k)
                         {
                          SSB_overlap_natal(k,n)=0;
                         }   
                        if(p==j && j!=k)
                         {
                          SSB_overlap_natal(k,n)=spawn_return_prob(p)*SSB_region_overlap(p,k,n,y);
                         } 
                       }
                      } 
                      if(p==j)
                      {
                       SSB_region_overlap(p,j,r,y)=SSB_region_overlap(p,j,r,y)+(sum(SSB_overlap_natal)/nregions(p));
                       }
                   } 
                SSB_population_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                SSB_population_overlap(p,j,y)=sum(SSB_population_temp_overlap(p,j,y));
                SSB_natal_overlap_temp(p,y,j)=SSB_population_overlap(p,j,y);
                SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));  /// this is adjusted below outside y loop to account for fish not spawning

              if(natal_homing_switch>0)
               {
               if(p==j)  //accounts for not being in natal area
               {
                SSB_region(j,r,y)=SSB_region_overlap(p,j,r,y);
               }
              }
              if(natal_homing_switch==0)
              {
                SSB_region_temp(j,r,y,a)=abundance_spawn(j,r,y,a)*wt_mat_mult_reg(j,r,y,a); 
                SSB_region(j,r,y)=sum(SSB_region_temp(j,r,y));
              }
                SSB_population_temp(j,y,r)=SSB_region(j,r,y);
                SSB_population(j,y)=sum(SSB_population_temp(j,y));
                SSB_total_temp(y,j)=SSB_population(j,y);
                SSB_total(y)=sum(SSB_total_temp(y));

                Bratio_population_overlap(p,j,y)=SSB_population_overlap(p,j,y)/SSB_zero(p);
                Bratio_natal_overlap(p,y)=SSB_natal_overlap(p,y)/SSB_zero(p);
                Bratio_population(j,y)=SSB_population(j,y)/SSB_zero(j);
                Bratio_total(y)=SSB_total(y)/sum(SSB_zero);

                //YIELD observation error
                OBS_yield_region_fleet_overlap(p,j,r,y,z)=yield_region_fleet_overlap(p,j,r,z,y)*mfexp(yield_RN_overlap(p,j,r,y,z)*sigma_catch_overlap(p,j,r,z)-.5*square(sigma_catch_overlap(p,j,r,z)));
                OBS_yield_fleet_temp(j,r,y,z,p)=OBS_yield_region_fleet_overlap(p,j,r,y,z);
               if(natal_homing_switch==0)
                {
                 OBS_yield_fleet(j,r,y,z)=yield_fleet(j,r,y,z)*mfexp(yield_RN(j,r,y,z)*sigma_catch(j,r,z)-.5*square(sigma_catch(j,r,z)));
                }
               if(natal_homing_switch==1)
                {
                 OBS_yield_fleet(j,r,y,z)=sum(OBS_yield_fleet_temp(j,r,y,z));  
                }
                OBS_yield_region_overlap(p,j,y,r)=sum(OBS_yield_region_fleet_overlap(p,j,r,y));
                OBS_yield_population_overlap(p,y,j)=sum(OBS_yield_region_overlap(p,j,y));
                OBS_yield_natal_overlap(y,p)=sum(OBS_yield_population_overlap(p,y));
                OBS_yield_total_overlap(y)=sum(OBS_yield_natal_overlap(y));
                OBS_yield_region(j,y,r)=sum(OBS_yield_fleet(j,r,y));
                OBS_yield_population(y,j)=sum(OBS_yield_region(j,y));
                OBS_yield_total(y)=sum(OBS_yield_population(y));
             //apportion variables
                apport_yield_region(j,r,y)=OBS_yield_region(j,y,r)/OBS_yield_population(y,j);

          } //end fleets loop
             for (int z=1;z<=nfleets_survey(j);z++)    /// survey index  1. Currently set up for more than 1 survey fleet
              {
               if(tsurvey(j,r)>0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                {
                  true_survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  true_survey_fleet_overlap_age_bio(p,j,r,y,z,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  true_survey_fleet_bio_overlap(p,j,r,y,z)=sum(true_survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  true_survey_fleet_bio_overlap_temp(j,r,y,z,p)=true_survey_fleet_bio_overlap(p,j,r,y,z);
                  OBS_survey_fleet_bio_overlap(p,j,r,y,z)=true_survey_fleet_bio_overlap(p,j,r,y,z)*mfexp(survey_RN_overlap(p,j,r,y,z)*sigma_survey_overlap(p,j,r,z)-.5*square(sigma_survey_overlap(p,j,r,z)));
                  OBS_survey_fleet_bio_temp(j,r,y,z,p)=OBS_survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  true_survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  true_survey_fleet_age_bio(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_age_bio(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=true_survey_fleet_bio(j,r,y,z)*mfexp(survey_RN(j,r,y,z)*sigma_survey(j,r,z)-.5*square(sigma_survey(j,r,z)));
                 }
                if(natal_homing_switch==1)
                 {
                  true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_bio_overlap_temp(j,r,y,z));
                  OBS_survey_fleet_bio(j,r,y,z)=sum(OBS_survey_fleet_bio_temp(j,r,y,z));  
                 }
             
                  true_survey_region_bio_overlap(p,j,y,r)=sum(true_survey_fleet_bio_overlap(p,j,r,y));               
                  true_survey_population_bio_overlap(p,y,j)=sum(true_survey_region_bio_overlap(p,j,y));               
                  true_survey_natal_bio_overlap(y,p)=sum(true_survey_population_bio_overlap(p,y));               
                  true_survey_total_bio_overlap(y)=sum(true_survey_natal_bio_overlap(y));
                  OBS_survey_region_bio_overlap(p,j,y,r)=sum(OBS_survey_fleet_bio_overlap(p,j,r,y));
                  OBS_survey_population_bio_overlap(p,y,j)=sum(OBS_survey_region_bio_overlap(p,j,y));
                  OBS_survey_natal_bio_overlap(y,p)=sum(OBS_survey_population_bio_overlap(p,y));
                  OBS_survey_total_bio_overlap(y)=sum(OBS_survey_natal_bio_overlap(y));
                  
                   // calculate true survey abundance for tagging simulation **dh**
                 true_survey_fleet_age_temp(j,r,y,a,z)=true_survey_fleet_age(j,r,y,z,a);
                 true_survey_region_abundance(j,y,r,a)=sum(true_survey_fleet_age_temp(j,r,y,a));
                 true_survey_region_abundance_temp(j,y,a,r)=true_survey_region_abundance(j,y,r,a);
                 true_survey_population_abundance(y,j,a)=sum(true_survey_region_abundance_temp(j,y,a));
                 true_survey_population_abundance_temp(y,a,j)=true_survey_population_abundance(y,j,a);
                 true_survey_total_abundance(y,a)=sum(true_survey_population_abundance_temp(y,a));

             // end abundance for tagging *dh

                  
                  true_survey_region_bio(j,y,r)=sum(true_survey_fleet_bio(j,r,y));
                  true_survey_population_bio(y,j)=sum(true_survey_region_bio(j,y));
                  true_survey_total_bio(y)=sum(true_survey_population_bio(y));
                  OBS_survey_region_bio(j,y,r)=sum(OBS_survey_fleet_bio(j,r,y));
                  OBS_survey_population_bio(y,j)=sum(OBS_survey_region_bio(j,y));
                  OBS_survey_total_bio(y)=sum(OBS_survey_population_bio(y));

                  apport_region_survey_biomass(j,r,y)=OBS_survey_region_bio(j,y,r)/OBS_survey_population_bio(y,j);
                  
                }  //tsurvey>0
               } //end survey_fleets
      }
     }
    }
   } // end age loop
  } //end yr>1 loop
 }  //end y loop

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////now adjusting natal homing ssb for spawning migration scenarios where fish left non-natal area
       for (int p=1;p<=npops;p++)
        {
         for (int j=1;j<=npops;j++)
          {
           for (int r=1;r<=nregions(p);r++)
            {
             for (int t=1;t<=nregions(j);t++)
              {
               for (int y=1;y<=nyrs;y++)
                {
                 SSB_overlap_natal=0;
                  if(natal_homing_switch==1 && spawn_return_switch==1)
                    {
                    if(p!=j)  //update SSB that doesn't spawn (ie doesn't return to natal population)
                    {
                     SSB_region_overlap(p,j,r,y)=(1-spawn_return_prob(p))*SSB_region_overlap(p,j,r,y);
                    }
                   SSB_population_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                   SSB_population_overlap(p,j,y)=sum(SSB_population_temp_overlap(p,j,y));
                   SSB_natal_overlap_temp(p,y,j)=SSB_population_overlap(p,j,y);
                   SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));
                   Bratio_population_overlap(p,j,y)=SSB_population_overlap(p,j,y)/SSB_zero(p);
                   Bratio_natal_overlap(p,y)=SSB_natal_overlap(p,y)/SSB_zero(p);
                   Bratio_population(j,y)=SSB_population(j,y)/SSB_zero(j);
                   Bratio_total(y)=SSB_total(y)/sum(SSB_zero);
                    }
                   }
                  }
                 }
                }
               }

  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for(int y=1;y<=nyrs;y++)
       {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T_year(j,r,y,k,n)=T(j,r,y,4,k,n);            
              } 
             }
            }
           }
          }

 ///manually calculating the true distribution of initial abundance for use in EM (so don't have to estimate)
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
    {
     for (int r=1;r<=nregions(j);r++)
      {
        init_abund_reg_temp(p,j,r)=sum(init_abund(p,j,r));
      }
     }
    }
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
    {
        init_abund_pop_temp(p,j)=sum(init_abund_reg_temp(p,j));
      }
     }
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
    {
     for (int r=1;r<=nregions(j);r++)
      {
       frac_natal_true(p,j,r)=sum(init_abund(p,j,r))/sum(init_abund_pop_temp(p));
      }
     }
    }

FUNCTION get_rand_survey_CAA_prop
 random_number_generator myrand_survey_age(myseed_survey_age);
 
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets_survey(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            for (int a=1;a<=nages;a++)
             {
              survey_at_age_region_fleet_overlap_prop(p,j,r,z,y,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)/sum(true_survey_fleet_overlap_age(p,j,r,y,z));              
              survey_at_age_fleet_prop(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)/sum(true_survey_fleet_age(j,r,y,z));
             }
            }
           }
          }
         }
        }

  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets_survey(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            if(use_stock_comp_info_survey==0)
             {
              rand_SIM_survey_prop_temp(j,r,y,z).fill_multinomial(myrand_survey_age,value(survey_at_age_fleet_prop(j,r,y,z)));
             }
            if(use_stock_comp_info_survey==1)
             {
              rand_SIM_survey_prop_temp_overlap(p,j,r,y,z).fill_multinomial(myrand_survey_age,value(survey_at_age_region_fleet_overlap_prop(p,j,r,z,y)));
             }
           }
         }
       }
     }
   }

  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets_survey(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
             rand_SIM_survey_prop_temp2=0;
            if(use_stock_comp_info_survey==0)
             {
               for(int n=1;n<=SIM_nsurvey(j,r,z);n++) 
                {
                 rand_SIM_survey_prop_temp2(value(rand_SIM_survey_prop_temp(j,r,y,z,n)))+= 1.0;
                }
               SIM_survey_prop(j,r,z,y)=rand_SIM_survey_prop_temp2;
             }
            if(use_stock_comp_info_survey==1)
             {              
               for(int n=1;n<=SIM_nsurvey_overlap(p,j,r,z);n++) /// look into changing this so can have ncatch change by year (ie different sample sizes for beginning and end of timeseries)
                {
                 rand_SIM_survey_prop_temp2(value(rand_SIM_survey_prop_temp_overlap(p,j,r,y,z,n)))+= 1.0;
                }
               SIM_survey_prop_overlap(p,j,r,z,y)=rand_SIM_survey_prop_temp2;
             }
             
        for(int a=1;a<=nages;a++)
         {
          if(use_stock_comp_info_survey==0)
           {
            OBS_survey_prop(j,r,y,z,a)=SIM_survey_prop(j,r,z,y,a)/SIM_nsurvey(j,r,z);
           }
          if(use_stock_comp_info_survey==1)
           {
            OBS_survey_prop_overlap(p,j,r,y,z,a)=SIM_survey_prop_overlap(p,j,r,z,y,a)/SIM_nsurvey_overlap(p,j,r,z);
           }
         }
       }
      }
     }
    }
   }

FUNCTION get_rand_CAA_prop
 random_number_generator myrand_catch_age(myseed_catch_age);
 
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            for (int a=1;a<=nages;a++)
             {
                 catch_at_age_region_fleet_overlap_prop(p,j,r,z,y,a)=catch_at_age_region_fleet_overlap(p,j,r,z,y,a)/sum(catch_at_age_region_fleet_overlap(p,j,r,z,y));              
                 catch_at_age_region_overlap_prop(p,j,r,y,a)=catch_at_age_region_overlap(p,j,r,y,a)/sum(catch_at_age_region_overlap(p,j,r,y));
                 catch_at_age_population_overlap_prop(p,j,y,a)=catch_at_age_population_overlap(p,j,y,a)/sum(catch_at_age_population_overlap(p,j,y));
                 catch_at_age_natal_overlap_prop(p,y,a)=catch_at_age_natal_overlap(p,y,a)/sum(catch_at_age_natal_overlap(p,y));

                 catch_at_age_fleet_prop_temp(j,r,y,z,a)=catch_at_age_fleet(j,r,y,a,z);
                 catch_at_age_fleet_prop(j,r,y,z,a)=catch_at_age_fleet_prop_temp(j,r,y,z,a)/sum(catch_at_age_fleet_prop_temp(j,r,y,z));
                 catch_at_age_region_prop(j,r,y,a)=catch_at_age_region(j,r,y,a)/sum(catch_at_age_region(j,r,y));
                 catch_at_age_population_prop(j,y,a)=catch_at_age_population(j,y,a)/sum(catch_at_age_population(j,y));
                 catch_at_age_total_prop(y,a)=catch_at_age_total(y,a)/sum(catch_at_age_total(y));
              }
            }
           }
          }
         }
        }
        
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            if(use_stock_comp_info_catch==0)
             {
              rand_SIM_catch_prop_temp(j,r,y,z).fill_multinomial(myrand_catch_age,value(catch_at_age_fleet_prop(j,r,y,z)));
             }
            if(use_stock_comp_info_catch==1)
             {
              rand_SIM_catch_prop_temp_overlap(p,j,r,y,z).fill_multinomial(myrand_catch_age,value(catch_at_age_region_fleet_overlap_prop(p,j,r,z,y)));
             }
           }
         }
       }
     }
   }
            
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
             rand_SIM_catch_prop_temp2=0;

            if(use_stock_comp_info_catch==0)
             {
               for(int n=1;n<=SIM_ncatch(j,r,z);n++) /// look into changing this so can have ncatch change by year (ie different sample sizes for beginning and end of timeseries)
                {
                 rand_SIM_catch_prop_temp2(value(rand_SIM_catch_prop_temp(j,r,y,z,n)))+= 1.0;
                }
               SIM_catch_prop(j,r,z,y)=rand_SIM_catch_prop_temp2;
             }
            if(use_stock_comp_info_catch==1)
             {              
               for(int n=1;n<=SIM_ncatch_overlap(p,j,r,z);n++) /// look into changing this so can have ncatch change by year (ie different sample sizes for beginning and end of timeseries)
                {
                 rand_SIM_catch_prop_temp2(value(rand_SIM_catch_prop_temp_overlap(p,j,r,y,z,n)))+= 1.0;
                }
               SIM_catch_prop_overlap(p,j,r,z,y)=rand_SIM_catch_prop_temp2;
             }
             
        for(int a=1;a<=nages;a++)
         {
          if(use_stock_comp_info_catch==0)
           {
            OBS_catch_prop(j,r,y,z,a)=SIM_catch_prop(j,r,z,y,a)/SIM_ncatch(j,r,z);
           }
          if(use_stock_comp_info_catch==1)
           {
            OBS_catch_prop_overlap(p,j,r,y,z,a)=SIM_catch_prop_overlap(p,j,r,z,y,a)/SIM_ncatch_overlap(p,j,r,z);
           }
         }
       }
      }
     }
    }
   }

FUNCTION get_tag_recaptures
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  // In some calcs (mortality and movement) need to account for true age (age+time at large),
  // because multiple cohorts this means need to factor in release year
  // so true age becomes age+year-release year
  // using subscript notation this equals= (a+y-xx)
  // similarly because account for release age, do not need to worry about plus group calcs as carry recaptures out to max_life_tags and never assume plus group (just use plus group mortality and movement values in calcs where a>=max_age)
  /////////////////////////////////////////////////////////////////////////////
  //reporting rate is assumed to be function of release event and recap location (not a function of recap year...could expand to this, but not priority at momement)
  ///////////////////////////////////////////////////////////////////////////////////

 if(do_tag==1)
  {
   for(int x=1; x<=nyrs_release; x++)
    {
     xx=yrs_releases(x);
     ntags_total_temp=0;
      for (int i=1;i<=npops;i++)
       {
         for (int n=1;n<=nregions(i);n++)
          {
           if(prob_tag_RN(i,n,x)>=opport_tag_prob) //set prob tag (for opportunistic tagging) to 1 or 0 (occurs/doesn't occur) based on uniform random number, if opport_Tag_prop==0.5 then 50/50 chance of tagging
           {
            prob_tag(i,n,x)=1.0;
           }
            if(prob_tag_RN(i,n,x)<opport_tag_prob)
           {
            prob_tag(i,n,x)=0.0;
           }
           if(prob_tag_year_RN(x)>=opport_tag_prob_year) //set prob tag (for opportunistic tagging) to 1 or 0 (occurs/doesn't occur) based on uniform random number, if opport_Tag_prop==0.5 then 50/50 chance of tagging
           {
            prob_tag_year(x)=1.0;
           }
            if(prob_tag_year_RN(x)<opport_tag_prob_year)
           {
            prob_tag_year(x)=0.0;
           }
         if(number_tags_switch==-2) //input total tags is distributed evenly across regions
          {
            ntags_region(i,n,x)=input_total_tags(x)/(sum(nregions)); //sum abundance across ages and all areas/pops
          }
         if(number_tags_switch==-1) //input total tags is distributed across regions according to survey abund
          {
            ntags_region(i,n,x)=input_total_tags(x)*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i))); //sum abundance across ages and all areas/pops; //sum abundance across ages by region
          }
         if(number_tags_switch==0) // tags by area is a fixed number based on input # tags
          {
            ntags_region(i,n,x)=input_ntags(i,n,x); 
          }
         if(number_tags_switch==1) //total tags is a fraction of total abundance and tags assigned evely to each region
          {
            ntags_region(i,n,x)=frac_abund_tagged(x)*sum(abundance_total(xx))/(sum(nregions)); //sum abundance across ages and all areas/pops
          }
         if(number_tags_switch==2) //total tags is a fraction of pop/reg abundance
          {
            ntags_region(i,n,x)=frac_abund_tagged(x)*sum(abundance_at_age_AM(i,n,xx)); //sum abundance across ages by region
          }
         if(number_tags_switch==3) //total tags is a fraction of total abundance then assigned to area based on survey abundance
          {
            ntags_region(i,n,x)=frac_abund_tagged(x)*sum(abundance_total(xx))*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i))); //sum abundance across ages and all areas/pops
          }
         if(number_tags_switch==4) // tags by area is a random number based on uniform dist with mean==((max+min)/2)
          {
           //number of tags is random number, but all areas receive tag releases (unless max and min set to 0 then no tags)
           ntags_region(i,n,x)=ntags_RN(i,n,x)*(max_tags(i,n)-min_tags(i,n))+min_tags(i,n); //scale uniform number ([0,1)) to interval that want
          }
         if(number_tags_switch==5) //number of tags is input and distributed proportional to survey abundance AND probability of tagging in a given AREA is random number (i.e., meant to simulate completely opportunistic tagging)
          {
           ntags_region(i,n,x)=input_total_tags(x)*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i)))*prob_tag(i,n,x); //if prob_tag==0 then no tagging in that release year/area
          }
         if(number_tags_switch==6) //number of tags is input and distributed proportional to survey abundance AND probability of tagging in a given YEAR is random number
          {
           ntags_region(i,n,x)=input_total_tags(x)*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i)))*prob_tag_year(x); //if prob_tag==0 then no tagging in that release year/area
          }
         if(number_tags_switch==7) //number of tags is random number AND probability of tagging in a given YEAR is random number AND probability of tagging in a given AREA is a random number
          {
           ntags_region(i,n,x)=(ntags_RN(i,n,x)*(max_tags(i,n)-min_tags(i,n))+min_tags(i,n))*prob_tag_year(x)*prob_tag(i,n,x); //if prob_tag==0 then no tagging in that release year/area
          }
         ntags_total_temp(i,n)=ntags_region(i,n,x);
          }
        }
         ntags_total(x)=sum(ntags_total_temp);
       }

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
       for (int a=1;a<=nages;a++) //release age 
        {
         if(age_dist_tags==0)  //distribute tags evenly across ages
          {
           ntags(i,n,x,a)=ntags_region(i,n,x)/(nages);
          }
         if(age_dist_tags==1)  //distribute tags across ages based on regional survey proportions at age (i.e., based on the proportions you would expect to catch fish in a survey...use primary survey in a region)
          {
           ntags(i,n,x,a)=ntags_region(i,n,x)*survey_at_age_fleet_prop(i,n,xx,1,a);
          }
         if(age_dist_tags==2)  //distribute tags across ages based on regional catch proportions at age (i.e., based on the proportions you would expect to catch fish in the fishery...use average across fleets in a region)
          {
           ntags(i,n,x,a)=ntags_region(i,n,x)*catch_at_age_region_prop(i,n,xx,a);
          }
        }
       }
      }
     }


////////////////////////////////////////////////////////////////////////////////////////////
 // for the EM .dat when have mismatches
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age 
        {
        ntags_temp1(i,x,a,n)=ntags(i,n,x,a);
        ntags_population(i,x,a)=sum(ntags_temp1(i,x,a));
        ntags_pan_temp(x,a,i)=ntags_population(i,x,a);
        ntags_pan(x,a)=sum(ntags_pan_temp(x,a));
        }
       }
      }
     }
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////     
//tag dynamics
//////////////////////////////////////////////////

 //assume tags released in natal population
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              if(y==1) //year of release for a cohort
               {
                if(sim_tag_mixing_switch==0) //assume complete mixing of tagged and untagged fish
                 {
                  tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T(i,n,xx,a,j,r); 
                  recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,xx,a)*(1.-mfexp(-(F(j,r,xx,a)+M(j,r,xx,a))))/(F(j,r,xx,a)+(M(j,r,xx,a)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)

               if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                {
                 tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T(i,n,xx,tag_age_sel,j,r); 
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,xx,tag_age_sel)*(1.-mfexp(-(F(j,r,xx,tag_age_sel)+M(j,r,xx,tag_age_sel))))/(F(j,r,xx,tag_age_sel)+(M(j,r,xx,tag_age_sel)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                                
                }
                
                 }
                if(sim_tag_mixing_switch==1) //assume incomplete mixing of tagged and untagged fish
                 {
                 if(sim_tag_mixing_T_switch==0) //assume incomplete mixing of tagged and untagged fish only for F
                  {
                   tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T(i,n,xx,a,j,r); 
                  }
                 if(sim_tag_mixing_T_switch==1) //assume incomplete mixing of tagged and untagged fish 
                  {                 
                   if(i==j && n==r)
                    {
                      T_tag(i,n,x,a,j,r)=T_tag_res(x);
                    }
                   if(i!=j || n!=r)
                    {
                     if(move_switch==8)
                     {
                      T_tag(i,n,x,a,j,r)=(1-T_tag_res(x))*rel_bio(i,n,xx,a,j,r);
                     }
                     if(move_switch!=8)
                     {
                      T_tag(i,n,x,a,j,r)=(1-T_tag_res(x))/(sum(nregions)-1); //split movement evenly across remaining regions
                     }
                    }
                      if(T_tag(i,n,x,a,j,r)>1)
                       {
                        T_tag(i,n,x,a,j,r)=1;
                       }
                      if(T_tag(i,n,x,a,j,r)<0)
                       {
                        T_tag(i,n,x,a,j,r)=0;
                       }
                    tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T_tag(i,n,x,a,j,r);                        
                   }
                  F_tag(j,r,x,a)=F_tag_scalar(x)*F(j,r,xx,a);
                  recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*(F_tag_scalar(x)*F(j,r,xx,a))*(1.-mfexp(-((F_tag_scalar(x)*F(j,r,xx,a))+M(j,r,xx,a))))/((F_tag_scalar(x)*F(j,r,xx,a))+(M(j,r,xx,a)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)
                 }
               }
            if(y==2) // must account for the maximum age so use min function to ensure that not exceeding the max age
             {
             if(sim_tag_mixing_switch==0) //assume complete mixing of tagged and untagged fish
              {
               tags_avail_temp=0;
               for(int p=1;p<=npops;p++)
               {
                for (int s=1;s<=nregions(p);s++)
                {
                 if(natal_homing_switch==0) //if no natal homing
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }                        
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current population only origin population and destination population
                //########################################################################################################              
                 if(natal_homing_switch==1) //if natal homing 
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(i,n,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }               
                }
               }
                 tags_avail(i,n,x,a,y,j,r)=sum(tags_avail_temp); //sum across all pops/regs of tags that moved into pop j reg r
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),min((a+y),nages))*(1.-mfexp(-(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))))))/(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)   

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),tag_age_sel)*(1.-mfexp(-(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)))))/(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
                    }    
               }
             if(sim_tag_mixing_switch==1) //assume incomplete mixing of tagged and untagged fish
              {
               tags_avail_temp=0;
               for(int p=1;p<=npops;p++)
               {
                for (int s=1;s<=nregions(p);s++)
                {
                 if(natal_homing_switch==0) //if no natal homing
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),min((a+y),nages),j,r)*mfexp(-((F_tag_scalar(x)*F(p,s,(xx+y-2),min(((a+y)-1),nages)))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                 }                        
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current population only origin population and destination population
                //########################################################################################################              
                 if(natal_homing_switch==1) //if natal homing 
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(i,n,(xx+y-1),min((a+y),nages),j,r)*mfexp(-((F_tag_scalar(x)*F(p,s,(xx+y-2),min(((a+y)-1),nages)))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)             
                 }               
                }
               }
                 tags_avail(i,n,x,a,y,j,r)=sum(tags_avail_temp); //sum across all pops/regs of tags that moved into pop j reg r
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),min((a+y),nages))*(1.-mfexp(-(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))))))/(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
               }
              }
              if(y>2) // must account for the maximum age so use min function to ensure that not exceeding the max age
              {
               tags_avail_temp=0;
               for(int p=1;p<=npops;p++)
               {
                for (int s=1;s<=nregions(p);s++)
                {
                 if(natal_homing_switch==0) //if no natal homing
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }                        
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current population only origin population and destination population
                //########################################################################################################              
                 if(natal_homing_switch==1) //if natal homing 
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(i,n,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }               
                }
               }
                 tags_avail(i,n,x,a,y,j,r)=sum(tags_avail_temp); //sum across all pops/regs of tags that moved into pop j reg r
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),min((a+y),nages))*(1.-mfexp(-(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))))))/(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),tag_age_sel)*(1.-mfexp(-(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)))))/(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
                    }  
               }
             }
            }
           }
          }
         }
        }
       }

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
        total_recap_temp.initialize();
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              total_recap_temp(j,r,y)=recaps(i,n,x,a,y,j,r);
             }
            }
           }
             total_rec(i,n,x,a)=sum(total_recap_temp);
             not_rec(i,n,x,a)=ntags(i,n,x,a)-total_rec(i,n,x,a);  //for ntags  at a given age all entries represent all tags released so can just use any of the entries (hence the i,x,a,1 subscripts)
           }
          }
         }
        }
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              if(ntags(i,n,x,a)>0)
               {
                tag_prop(i,n,x,a,y,j,r)=recaps(i,n,x,a,y,j,r)/ntags(i,n,x,a);
                tag_prop_not_rec(i,n,x,a)=not_rec(i,n,x,a)/ntags(i,n,x,a);
               }
              if(ntags(i,n,x,a)==0)
               {                   
                tag_prop(i,n,x,a,y,j,r)=0;
                tag_prop_not_rec(i,n,x,a)=0;
               }
             } 
            }
           }
          }
         }
        }
       }

  
 //in order to use the multinomial RNG, need to have a vector of probabilities
 //this essentially requires stacking the tag_prop array into vectors for each release cohort covering all recap states (recap year, population, region)
 //need to extract each recap prob vector and store, then combine into array where can extract the last index (essentially the columns of the array)
 //can then use the fill.multinomial with that vector of probabilities
 if(tag_fit_ages_switch==0) // fit cohorts by region and age
 {
 nreg_temp=rowsum(nregions_temp);
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
       for (int a=1;a<=nages;a++) //release age 
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              tag_prop_temp(j,y,r)=0; //for whatever reason ADMB won't let a 3darray=double...this is workaround for total_recap_temp=0;, ie setting temp 3darray to 0 after loops through all years 
             }
            }
           }
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++) //recap region
              {
               tag_prop_temp(j,y,r)=tag_prop(i,n,x,a,y,j,r); //fill temp array with a single release cohort's recapture probabilities (excluding not recaptured)
              }
            }
           }
         //tag_prop_temp2=0;
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
              {
               tag_prop_temp2(i,n,x,a,((y-1)*sum(nregions)+nreg_temp(j)+r))=tag_prop_temp(j,y,r); //stack array into single vector by year, population, region
              }
             }
            }
             for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
              {
              if(s<(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1))
              {
               tag_prop_final(i,n,x,a,s)=tag_prop_temp2(i,n,x,a,s);
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags<=(nyrs-xx+1)) //add not recap probability to final entry of temp array
              {
               tag_prop_final(i,n,x,a,s)=tag_prop_not_rec(i,n,x,a);  //for estimation model will use this version of tag_prop in likelihood
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags>(nyrs-xx+1)) //add not recap probability to final entry of temp array with adjustment for release events where model ends before max_life_tags (so NR remains in last state)
              {
               tag_prop_final(i,n,x,a,(max_life_tags*sum(nregions)+1))=tag_prop_not_rec(i,n,x,a);  //for estimation model will use this version of tag_prop in likelihood
              }
            }
           }
          }
         }
        }
       }
 if(tag_fit_ages_switch==1) //only fit cohorts by region not age
 {
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
               tag_recap_no_age_temp(i,n,x,y,j,r,a)=recaps(i,n,x,a,y,j,r);
               total_rec_no_age(i,n,x)=sum(total_rec(i,n,x));
               not_rec_no_age(i,n,x)=sum(not_rec(i,n,x));
             }
            }
           }
          }
        ntags_no_age(i,n,x)=sum(ntags(i,n,x));
       }
      }
     }

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              if(ntags_no_age(i,n,x)>0)
               {
                tag_prop_no_age(i,n,x,y,j,r)=sum(tag_recap_no_age_temp(i,n,x,y,j,r))/ntags_no_age(i,n,x);
                tag_prop_not_rec_no_age(i,n,x)=not_rec_no_age(i,n,x)/ntags_no_age(i,n,x);
               }
              if(ntags_no_age(i,n,x)==0)
               {                   
                tag_prop_no_age(i,n,x,y,j,r)=0;
                tag_prop_not_rec_no_age(i,n,x)=0;
               }
             } 
            }
           }
          }
         }
        }
       }
 nreg_temp=rowsum(nregions_temp);
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              tag_prop_temp(j,y,r)=0; //for whatever reason ADMB won't let a 3darray=double...this is workaround for total_recap_temp=0;, ie setting temp 3darray to 0 after loops through all years 
             }
            }
           }
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++) //recap region
              {
               tag_prop_temp(j,y,r)=tag_prop_no_age(i,n,x,y,j,r); //fill temp array with a single release cohort's recapture probabilities (excluding not recaptured)
              }
            }
           }
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
              {
               tag_prop_temp2_no_age(i,n,x,((y-1)*sum(nregions)+nreg_temp(j)+r))=tag_prop_temp(j,y,r); //stack array into single vector by year, population, region
              }
             }
            }
             for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
              {
              if(s<(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1))
              {
               tag_prop_final_no_age(i,n,x,s)=tag_prop_temp2_no_age(i,n,x,s);
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags<=(nyrs-xx+1)) //add not recap probability to final entry of temp array
              {
               tag_prop_final_no_age(i,n,x,s)=tag_prop_not_rec_no_age(i,n,x);  //for estimation model will use this version of tag_prop in likelihood
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags>(nyrs-xx+1)) //add not recap probability to final entry of temp array with adjustment for release events where model ends before max_life_tags (so NR remains in last state)
              {
               tag_prop_final_no_age(i,n,x,(max_life_tags*sum(nregions)+1))=tag_prop_not_rec_no_age(i,n,x);  //for estimation model will use this version of tag_prop in likelihood
              }
            }
           }
          }
         }
        } //end tag_age_switch_loop
       } //end_do_tag loop
FUNCTION get_observed_tag_recaptures
 if(do_tag==1)
  {
  
 random_number_generator myrand_tag(myseed_tag);
 
 if(tag_fit_ages_switch==0) // fit cohorts by region and age
 {
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
       for (int a=1;a<=nages;a++) //release age 
        {
         rand_tag_prop_temp(i,n,x,a).fill_multinomial(myrand_tag,value(tag_prop_final(i,n,x,a))); //fill multinomial using recap prop for each cohort
        }
       }
      }
     }
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
       for (int a=1;a<=nages;a++) //release age 
        {
          rand_tag_prop_temp2=0;
            for(int s=1;s<=SIM_ntag;s++) /// look into changing this so can have ntag change by year (ie different sample sizes for beginning and end of timeseries)
             {
               rand_tag_prop_temp2(value(rand_tag_prop_temp(i,n,x,a,s)))+= 1.0; //count up number in each recap state
              }
           SIM_tag_prop(i,n,x,a)=rand_tag_prop_temp2; //put number in each recap state into appropriate release cohort
        }
       }
      }
     }

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
       for (int a=1;a<=nages;a++) //release age 
        {
         for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
          {           
           OBS_tag_prop_final(i,n,x,a,s)=SIM_tag_prop(i,n,x,a,s)/SIM_ntag;
              if(ntags(i,n,x,a)==0)
               {                   
                OBS_tag_prop_final(i,n,x,a,s)=0; 
               }
           //summarize by population for panmictic EM
             OBS_tag_prop_population_temp(i,x,a,s,n)=SIM_tag_prop(i,n,x,a,s);
             OBS_tag_prop_population_temp2(i,x,a,s)=sum(OBS_tag_prop_population_temp(i,x,a,s));
             OBS_tag_prop_population_final(i,x,a,s)=OBS_tag_prop_population_temp2(i,x,a,s)/(nregions(i)*SIM_ntag);
           //  OBS_tag_prop_pan_temp(x,a,s,i)=OBS_tag_prop_population_temp2(i,x,a,s);
         //    OBS_tag_prop_pan_temp2(x,a,s)=sum(OBS_tag_prop_pan_temp(x,a,s));

     //        if(OM_structure>0 && EM_structure==0){
      //            for (int k=1;k<=max_life_tags+1;k++)
       //               {
       //                OBS_tag_prop_pan_final_temp(x,a,k)=OBS_tag_prop_pan_temp2(x,a,k)+OBS_tag_prop_pan_temp2(x,a,k+max_life_tags);
       //                OBS_tag_prop_pan_final_temp(x,a,max_life_tags+1)=OBS_tag_prop_pan_temp2(x,a,max_life_tags*sum(nregions)+1);// this is a nightmare
        //               OBS_tag_prop_pan_final(x,a,k)=OBS_tag_prop_pan_final_temp(x,a,k)/(SIM_ntag*sum(nregions));//this probably wont work with 3 areas
                       //OBS_tag_prop_pan_final(x,a,k)=1; //for place holding only
       //                }
        //              }
                    }
                    }
                   }
                  }
                 }
                }

 if(tag_fit_ages_switch==1) // fit cohorts by region and age
 {
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
         rand_tag_prop_temp_no_age(i,n,x).fill_multinomial(myrand_tag,value(tag_prop_final_no_age(i,n,x))); //fill multinomial using recap prop for each cohort
        }
       }
      }
     
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
          rand_tag_prop_temp2=0;
            for(int s=1;s<=SIM_ntag;s++) /// look into changing this so can have ntag change by year (ie different sample sizes for beginning and end of timeseries)
             {
               rand_tag_prop_temp2(value(rand_tag_prop_temp_no_age(i,n,x,s)))+= 1.0; //count up number in each recap state
              }
           SIM_tag_prop_no_age(i,n,x)=rand_tag_prop_temp2; //put number in each recap state into appropriate release cohort
        }
       }
      }
     

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
   for(int x=1; x<=nyrs_release; x++)
    {
     xx=yrs_releases(x);
         for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
          {           
           OBS_tag_prop_final_no_age(i,n,x,s)=SIM_tag_prop_no_age(i,n,x,s)/SIM_ntag;
              if(ntags_no_age(i,n,x)==0)
               {                   
                OBS_tag_prop_final_no_age(i,n,x,s)=0; 
              }
           //summarize by population for panmictic EM
//             OBS_tag_prop_population_temp_no_age(i,x,s,n)=SIM_tag_prop_no_age(i,n,x,s);
//             OBS_tag_prop_population_temp2_no_age(i,x,s)=sum(OBS_tag_prop_population_temp_no_age(i,x,s));
 //            OBS_tag_prop_population_final_no_age(i,x,s)=OBS_tag_prop_population_temp2_no_age(i,x,s)/(nregions(i)*SIM_ntag);
           //  OBS_tag_prop_pan_temp(x,a,s,i)=OBS_tag_prop_population_temp2(i,x,a,s);
         //    OBS_tag_prop_pan_temp2(x,a,s)=sum(OBS_tag_prop_pan_temp(x,a,s));

     //        if(OM_structure>0 && EM_structure==0){
      //            for (int k=1;k<=max_life_tags+1;k++)
       //               {
       //                OBS_tag_prop_pan_final_temp(x,a,k)=OBS_tag_prop_pan_temp2(x,a,k)+OBS_tag_prop_pan_temp2(x,a,k+max_life_tags);
       //                OBS_tag_prop_pan_final_temp(x,a,max_life_tags+1)=OBS_tag_prop_pan_temp2(x,a,max_life_tags*sum(nregions)+1);// this is a nightmare
        //               OBS_tag_prop_pan_final(x,a,k)=OBS_tag_prop_pan_final_temp(x,a,k)/(SIM_ntag*sum(nregions));//this probably wont work with 3 areas
                       //OBS_tag_prop_pan_final(x,a,k)=1; //for place holding only
       //                }
                     // }
                    }
                    }
                   }
                  }
                 
                } //end do_age_tag_switch loop
  } //end do_tag loop

FUNCTION evaluate_the_objective_function
   f=0.0;

           for (int y=1;y<=nyrs;y++)
            {
             for (int j=1;j<=npops;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j);z++)
                  {
                    res_TAC(j,r,z,y)=input_TAC(j,r,z)-yield_fleet(j,r,y,z);
                  }
                 res_u(j,r,y)=sum(input_u(j,r))-harvest_rate_region_bio(j,r,y);
                }
               }
              }

REPORT_SECTION

 //for the mismatch calculate the abundance fraction by region/population
  for (int y=1;y<=nyrs;y++)
   {
   for (int p=1;p<=npops;p++)
    {
    for (int r=1;r<=nregions(p);r++)
      {
      for (int z=1;z<=nfleets(p);z++)
        {
         for(int a=1;a<=nages;a++)
          {
        abund_frac_age_region(p,r,y,a)=abundance_at_age_AM(p,r,y,a)/abundance_total(y,a);
        abund_frac_region_year(p,r,y)=sum(abundance_at_age_AM(p,r,y))/sum(abundance_total(y));
        abund_frac_region(p,r)=sum(abund_frac_region_year(p,r))/nyrs;//average of all years
        }}}}}


 //Aggregating OBS values and vitals for panmictic EM
 if(EM_structure==0 && OM_structure>0){ 
  for (int y=1;y<=nyrs;y++)
   {
   for (int p=1;p<=npops;p++)
    {
    for (int r=1;r<=nregions(p);r++)
      {
      for (int z=1;z<=nfleets(p);z++)
        {
         for(int a=1;a<=nages;a++)
          {
       //aggregating weight at age
        input_weight_region_temp(p,a,r)=input_weight(p,r,a)*abund_frac_region(p,r);//rearrange to summarize and weight for output
        input_weight_region(p,a)=sum(input_weight_region_temp(p,a));
        input_weight_population_temp(a,p)=input_weight_region(p,a);
        input_weight_population(a)=sum(input_weight_population_temp(a));

        //aggregating catch weight at age
        input_catch_weight_region_temp(p,a,r)=input_catch_weight(p,r,a)*abund_frac_region(p,r);//sum by region
        input_catch_weight_region(p,a)=sum(input_catch_weight_region_temp(p,a));
        input_catch_weight_population_temp(a,p)=input_catch_weight_region(p,a);
        input_catch_weight_population(a)=sum(input_catch_weight_population_temp(a));

        //aggregating fecundity
        fecundity_region_temp(p,a,r)=fecundity(p,r,a)*abund_frac_region(p,r);//sum by region
        fecundity_region(p,a)=sum(fecundity_region_temp(p,a));
        fecundity_population_temp(a,p)=fecundity_region(p,a);
        fecundity_population(a)=sum(fecundity_population_temp(a));

        //aggregating maturity
        maturity_region_temp(p,a,r)=maturity(p,r,a)*abund_frac_region(p,r);//sum by region
        maturity_region(p,a)=sum(maturity_region_temp(p,a));
        maturity_population_temp(a,p)=maturity_region(p,a);
        maturity_population(a)=sum(maturity_population_temp(a));

        //aggregating selectivity
        selectivity_age_temp(p,a,z,r)=selectivity_age(p,r,a,z)*abund_frac_region(p,r);
        selectivity_age_pop(p,a,z)=sum(selectivity_age_temp(p,a,z));
        survey_selectivity_age_temp(p,a,z,r)=survey_selectivity_age(p,r,a,z)*abund_frac_region(p,r);
        survey_selectivity_age_pop(p,a,z)=sum(survey_selectivity_age_temp(p,a,z));

        //aggregating the age comps

        //survey
        OBS_survey_prop_temp(p,r,y,a,z)=OBS_survey_prop(p,r,z,y,a);
        OBS_survey_prop_temp2(p,r,y,a)=sum(OBS_survey_prop_temp(p,r,y,a));
        OBS_survey_prop_temp3(p,r,y,a)=OBS_survey_prop_temp2(p,r,y,a)*abund_frac_age_region(p,r,y,a);
        OBS_survey_prop_temp4(p,y,a,r)=OBS_survey_prop_temp3(p,r,y,a);
        OBS_survey_prop_population(p,y,a)=sum(OBS_survey_prop_temp4(p,y,a));
        OBS_survey_prop_pan_temp(y,a,p)= OBS_survey_prop_population(p,y,a);
        OBS_survey_prop_pan(y,a)=sum(OBS_survey_prop_pan_temp(y,a));

        //catch
        OBS_catch_prop_temp(p,r,y,a,z)= OBS_catch_prop(p,r,z,y,a);
        OBS_catch_prop_temp2(p,r,y,a)=sum(OBS_catch_prop_temp(p,r,y,a));
        OBS_catch_prop_temp3(p,r,y,a)= OBS_catch_prop_temp2(p,r,y,a)*abund_frac_age_region(p,r,y,a);
        OBS_catch_prop_temp4(p,y,a,r)= OBS_catch_prop_temp3(p,r,y,a);
        OBS_catch_prop_population(p,y,a)=sum(OBS_catch_prop_temp4(p,y,a));
        OBS_catch_prop_pan_temp(y,a,p)= OBS_catch_prop_population(p,y,a);
        OBS_catch_prop_pan(y,a)=sum(OBS_catch_prop_pan_temp(y,a));


        } //end age loop
        } //end fleets loop

        //proportion female
        prop_fem_temp(p,r)= prop_fem(p,r)*abund_frac_region(p,r); 
        prop_fem_pan=sum(prop_fem_temp);
        rec_index_temp(p,y,r)=rec_index_BM(p,r,y)*abund_frac_region_year(p,r,y); //rearrange and weight for summing

       } //end reg loop

       //rec index
        rec_index_BM_population(p,y)=sum(rec_index_temp(p,y));// combined by region
        rec_index_temp2(y,p)=rec_index_BM_population(p,y);
        rec_index_pan(y)=sum(rec_index_temp2(y));//npops; combined by populations
        
      } //end pop loop          
     } //end year loop
  }


//Additional model structure parameters
  report<<"#nages"<<endl;
  report<<nages<<endl;
  report<<"#nyrs"<<endl;
  report<<nyrs<<endl;

//EM structure 
  report<<"#npops_EM"<<endl;
  report<<npops_EM<<endl;
  report<<"#nregions_EM"<<endl;
  report<<nregions_EM<<endl;
  report<<"#nfleets_EM"<<endl;
  report<<nfleets_EM<<endl;
  report<<"#nfleets_survey_EM"<<endl;
  report<<nfleets_survey_EM<<endl;
  
//OM structure 
  report<<"#npops_OM"<<endl;
  report<<npops<<endl;
  report<<"#nregions_OM"<<endl;
  report<<nregions<<endl;
  report<<"#nfleets_OM"<<endl;
  report<<nfleets<<endl;
  report<<"#nfleets_survey_OM"<<endl;
  report<<nfleets_survey<<endl;

//EM parameters input from OM .dat
  report<<"#tsurvey_EM"<<endl;
  report<<tsurvey_EM<<endl;
  report<<"#diagnostics_switch"<<endl;
  report<<diagnostics_switch<<endl;
  report<<"#fleets_as_areas_switch"<<endl;
  report<<fleets_as_areas_switch<<endl;
  report<<"#move_switch"<<endl;
  report<<move_switch_EM<<endl;
  report<<"#report_rate_switch"<<endl;
  report<<report_rate_switch_EM<<endl;
  report<<"#natal_homing_switch"<<endl;
  report<<natal_homing_switch_EM<<endl;
  report<<"#spawn_return_switch"<<endl;
  report<<spawn_return_switch_EM<<endl;
  report<<"#select_switch"<<endl;
  report<<select_switch_EM<<endl;
  report<<"#select_switch_survey"<<endl;
  report<<select_switch_survey_EM<<endl;
  report<<"#maturity_switch_equil"<<endl;
  report<<maturity_switch_equil_EM<<endl;
  report<<"#SSB_type"<<endl;
  report<<SSB_type_EM<<endl;
  report<<"#Rec_type"<<endl;
  report<<Rec_type_EM<<endl;
  report<<"#apportionment_type"<<endl;
  report<<apportionment_type_EM<<endl;
  report<<"#use_stock_comp_info_survey"<<endl;
  report<<use_stock_comp_info_survey_EM<<endl;
  report<<"#use_stock_comp_info_catch"<<endl;
  report<<use_stock_comp_info_catch_EM<<endl;
  report<<"#F_switch"<<endl;
  report<<F_switch_EM<<endl;
  report<<"#M_switch"<<endl;
  report<<M_switch_EM<<endl;
  report<<"#recruit_devs_switch"<<endl;
  report<<recruit_devs_switch_EM<<endl;
  report<<"#recruit_randwalk_switch"<<endl;
  report<<recruit_randwalk_switch_EM<<endl;
  report<<"#init_abund_switch"<<endl;
  report<<init_abund_switch_EM<<endl;
  report<<"#est_dist_init_abund"<<endl;
  report<<est_dist_init_abund_EM<<endl;
  report<<"#tspawn_EM"<<endl;
  report<<tspawn_EM<<endl;
  report<<"#return_age"<<endl;
  report<<return_age<<endl;
  report<<"#return_probability"<<endl;
  report<<return_probability_EM<<endl;
  report<<"#spawn_return_prob"<<endl;
  report<<spawn_return_prob_EM<<endl;
  report<<"#do_tag"<<endl;
  report<<do_tag_EM<<endl;
  report<<"#fit_tag_age_switch"<<endl;
  report<<tag_fit_ages_switch<<endl;
  report<<"#do_tag_mult"<<endl;
  report<<do_tag_mult<<endl;
  report<<"#est_tag_mixing_switch"<<endl;
  report<<est_tag_mixing_switch<<endl;

  report<<"#sigma_recruit_EM"<<endl;
  report<<sigma_recruit_EM<<endl;

  report<<"#ph_lmr"<<endl;
  report<<ph_lmr<<endl;
  report<<"#lb_R_ave"<<endl;
  report<<lb_R_ave<<endl;
  report<<"#ub_R_ave"<<endl;
  report<<ub_R_ave<<endl;
  report<<"#ph_rec"<<endl;
  report<<ph_rec<<endl;
  report<<"#lb_rec_devs"<<endl;
  report<<lb_rec_devs<<endl;
  report<<"#ub_rec_devs"<<endl;
  report<<ub_rec_devs<<endl;
  report<<"#ph_rec_app_CNST"<<endl;
  report<<ph_rec_app_CNST<<endl;
  report<<"#ph_rec_app_YR"<<endl;
  report<<ph_rec_app_YR<<endl;
  report<<"#lb_rec_app"<<endl;
  report<<lb_rec_app<<endl;
  report<<"#ub_rec_app"<<endl;
  report<<ub_rec_app<<endl;
  report<<"#ph_init_abund"<<endl;
  report<<ph_init_abund<<endl;
  report<<"#ph_reg_init"<<endl;
  report<<ph_reg_init<<endl;
  report<<"#ph_non_natal_init"<<endl;
  report<<ph_non_natal_init<<endl;
  report<<"#lb_init_dist"<<endl;
  report<<lb_init_dist<<endl;
  report<<"#ub_init_dist"<<endl;
  report<<ub_init_dist<<endl;
  report<<"#lb_init_abund"<<endl;
  report<<lb_init_abund<<endl;
  report<<"#ub_init_abund"<<endl;
  report<<ub_init_abund<<endl;
  report<<"#ph_F"<<endl;
  report<<ph_F<<endl;
  report<<"#lb_F"<<endl;
  report<<lb_F<<endl;
  report<<"#ub_F"<<endl;
  report<<ub_F<<endl;
  report<<"#ph_steep"<<endl;
  report<<ph_steep<<endl;
  report<<"#lb_steep"<<endl;
  report<<lb_steep<<endl;
  report<<"#ub_steep"<<endl;
  report<<ub_steep<<endl;
  report<<"#ph_M_CNST"<<endl;
  report<<ph_M_CNST<<endl;
  report<<"#ph_M_pop_CNST"<<endl;
  report<<ph_M_pop_CNST<<endl;
  report<<"#ph_M_age_CNST"<<endl;
  report<<ph_M_age_CNST<<endl;
  report<<"#ph_M_pop_age"<<endl;
  report<<ph_M_pop_age<<endl;
  report<<"#lb_M"<<endl;
  report<<lb_M<<endl;
  report<<"#ub_M"<<endl;
  report<<ub_M<<endl;
  report<<"#ph_sel_log"<<endl;
  report<<ph_sel_log<<endl;
  report<<"#lb_sel_beta1"<<endl;
  report<<lb_sel_beta1<<endl;
  report<<"#ub_sel_beta1"<<endl;
  report<<ub_sel_beta1<<endl;
  report<<"#lb_sel_beta2"<<endl;
  report<<lb_sel_beta2<<endl;
  report<<"#ub_sel_beta2"<<endl;
  report<<ub_sel_beta2<<endl;
  report<<"#lb_sel_beta3"<<endl;
  report<<lb_sel_beta3<<endl;
  report<<"#ub_sel_beta3"<<endl;
  report<<ub_sel_beta3<<endl;
  report<<"#lb_sel_beta4"<<endl;
  report<<lb_sel_beta4<<endl;
  report<<"#ub_sel_beta4"<<endl;
  report<<ub_sel_beta4<<endl;
  report<<"#lb_sel_beta1_surv"<<endl;
  report<<lb_sel_beta1_surv<<endl;
  report<<"#ub_sel_beta1_surv"<<endl;
  report<<ub_sel_beta1_surv<<endl;
  report<<"#lb_sel_beta2_surv"<<endl;
  report<<lb_sel_beta2_surv<<endl;
  report<<"#ub_sel_beta2_surv"<<endl;
  report<<ub_sel_beta2_surv<<endl;
  report<<"#lb_sel_beta3_surv"<<endl;
  report<<lb_sel_beta3_surv<<endl;
  report<<"#ub_sel_beta3_surv"<<endl;
  report<<ub_sel_beta3_surv<<endl;
  report<<"#lb_sel_beta4_surv"<<endl;
  report<<lb_sel_beta4_surv<<endl;
  report<<"#ub_sel_beta4_surv"<<endl;
  report<<ub_sel_beta4_surv<<endl;
  report<<"#ph_sel_log_surv"<<endl;
  report<<ph_sel_log_surv<<endl;
  report<<"#ph_sel_dubl"<<endl;
  report<<ph_sel_dubl<<endl;
  report<<"#ph_sel_dubl_surv"<<endl;
  report<<ph_sel_dubl_surv<<endl;
  report<<"#ph_q"<<endl;
  report<<ph_q<<endl;
  report<<"#lb_q"<<endl;
  report<<lb_q<<endl;
  report<<"#ub_q"<<endl;
  report<<ub_q<<endl;
  report<<"#ph_F_rho"<<endl;
  report<<ph_F_rho<<endl;
  report<<"#lb_F_rho"<<endl;
  report<<lb_F_rho<<endl;
  report<<"#ub_F_rho"<<endl;
  report<<ub_F_rho<<endl;
  report<<"#phase_T_YR"<<endl;
  report<<phase_T_YR<<endl;
  report<<"#phase_T_YR_ALT_FREQ"<<endl;
  report<<phase_T_YR_ALT_FREQ<<endl;
  report<<"#T_est_freq"<<endl;
  report<<T_est_freq<<endl;
  report<<"#phase_T_YR_AGE_ALT_FREQ"<<endl;
  report<<phase_T_YR_AGE_ALT_FREQ<<endl;
  report<<"#T_est_age_freq"<<endl;
  report<<T_est_age_freq<<endl;
  report<<"#juv_age"<<endl;
  report<<juv_age<<endl;  
  report<<"#phase_T_CNST"<<endl;
  report<<phase_T_CNST<<endl;
  report<<"#phase_T_CNST_AGE"<<endl;
  report<<phase_T_CNST_AGE<<endl;
  report<<"#phase_T_YR_AGE"<<endl;
  report<<phase_T_YR_AGE<<endl;
  report<<"#phase_T_CNST_AGE_no_AG1"<<endl;
  report<<phase_T_CNST_AGE_no_AG1<<endl;
  report<<"#phase_T_YR_AGE_no_AG1"<<endl;
  report<<phase_T_YR_AGE_no_AG1<<endl;
  report<<"#phase_T_YR_AGE_ALT_FREQ_no_AG1"<<endl;
  report<<phase_T_YR_AGE_ALT_FREQ_no_AG1<<endl;
  report<<"#lb_T"<<endl;
  report<<lb_T<<endl;
  report<<"#ub_T"<<endl;
  report<<ub_T<<endl;
  report<<"#T_start"<<endl;
  report<<T_start<<endl;
  report<<"#phase_rep_rate_YR"<<endl;
  report<<phase_rep_rate_YR<<endl;
  report<<"#phase_rep_rate_CNST"<<endl;
  report<<phase_rep_rate_CNST<<endl;
  report<<"#lb_B"<<endl;
  report<<lb_B<<endl;
  report<<"#ub_B"<<endl;
  report<<ub_B<<endl;
  
  report<<"#ph_T_tag"<<endl;
  report<<ph_T_tag<<endl;
  report<<"#ph_F_tag"<<endl;
  report<<ph_F_tag<<endl;
  report<<"#lb_scalar"<<endl;
  report<<lb_scalar<<endl;
  report<<"#ub_scalar"<<endl;
  report<<ub_scalar<<endl;
  
  report<<"#ph_dummy"<<endl;
  report<<ph_dummy<<endl;
  report<<"#wt_surv"<<endl;
  report<<wt_surv<<endl;
  report<<"#wt_catch"<<endl;
  report<<wt_catch<<endl;
  report<<"#wt_fish_age"<<endl;
  report<<wt_fish_age<<endl;
  report<<"#wt_srv_age"<<endl;
  report<<wt_srv_age<<endl;
  report<<"#wt_rec"<<endl;
  report<<wt_rec<<endl;
  report<<"#wt_tag"<<endl;
  report<<wt_tag<<endl;
  report<<"#wt_F_pen"<<endl;
  report<<wt_F_pen<<endl;
  report<<"#wt_M_pen"<<endl;
  report<<wt_M_pen<<endl;
  report<<"#wt_B_pen"<<endl;
  report<<wt_B_pen<<endl;
  report<<"#report_rate_sigma"<<endl;
  report<<report_rate_sigma<<endl;
  report<<"#report_rate_ave"<<endl;
  report<<report_rate_ave<<endl;

  report<<"#abund_pen_switch"<<endl;
  report<<abund_pen_switch<<endl;
  report<<"#wt_abund_pen"<<endl;
  report<<wt_abund_pen<<endl;
  report<<"#move_pen_switch"<<endl;
  report<<move_pen_switch<<endl;
  report<<"#wt_move_pen"<<endl;
  report<<wt_move_pen<<endl;
  report<<"#Tpen"<<endl;
  report<<Tpen<<endl;
  report<<"#Tpen2"<<endl;
  report<<Tpen2<<endl;
  report<<"#sigma_Tpen_EM"<<endl;
  report<<sigma_Tpen_EM<<endl;
 
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
///REPORTING THE CORRECT EM PARAMETERS///////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

//Spatial to panmictic EM inputs
   if(EM_structure==0 && OM_structure>=0){ 
      report<<"#input_weight"<<endl;
      report<<input_weight_population<<endl;
      report<<"#input_catch_weight"<<endl;
      report<<input_catch_weight_population<<endl;
      report<<"#fecundity"<<endl;
      report<<fecundity_population<<endl;
      report<<"#maturity"<<endl;
      report<<maturity_population<<endl;
      report<<"#prop_fem"<<endl; 
      report<<prop_fem_pan<<endl;
      report<<"#OBS_rec_index_BM"<<endl;
      report<<rec_index_pan<<endl;

 //for straight panmictic EM 
    if(sum(nfleets_EM)==1){
      report<<"#OBS_survey_fleet"<<endl;
      report<<OBS_survey_total_bio<<endl;
      report<<"#OBS_survey_fleet_bio_se_EM"<<endl;
      report<<OBS_survey_fleet_bio_se_EM<<endl;
      report<<"#OBS_survey_prop"<<endl;
      report<<OBS_survey_prop_pan<<endl;
      report<<"#OBS_survey_prop_N_EM"<<endl;
      report<<OBS_survey_prop_N_EM<<endl;
      report<<"#OBS_yield_fleet"<<endl;
      report<<OBS_yield_total<<endl;
      report<<"#OBS_yield_fleet_se_EM"<<endl;
      report<<OBS_yield_fleet_se_EM<<endl;    
      report<<"#OBS_catch_prop"<<endl;
      report<<OBS_catch_prop_pan<<endl; 
      report<<"#OBS_catch_prop_N_EM"<<endl;
      report<<OBS_catch_prop_N_EM<<endl;
      
//tagging information
      report<<"#nyrs_release"<<endl;
      report<<nyrs_release<<endl;
      report<<"#years_of_tag_releases "<<endl;
      report<<yrs_releases<<endl;
      report<<"#max_life_tags"<<endl;
      report<<max_life_tags<<endl;
      report<<"#age_full_selection"<<endl;
      report<<age_full_selection<<endl;
      report<<"#input_report_rate_EM"<<endl;
      report<<input_report_rate_EM<<endl;
      report<<"#ntags"<<endl;
      report<<ntags_pan<<endl;
      report<<"#ntags_total"<<endl;
      report<<ntags_total<<endl;
      report<<"#tag_N_EM"<<endl;
      report<<tag_N_EM<<endl;
      report<<"#input_T_EM"<<endl;
      report<<input_T_EM<<endl;
      report<<"#OBS_tag_prop_pan_final"<<endl;
      report<<OBS_tag_prop_pan_final<<endl;
     // report<<"#OBS_tag_prop_pan_final_no_age"<<endl;
     // report<<OBS_tag_prop_pan_final_no_age<<endl;
      }


 //for fleets-as-areas approach
      if(sum(nfleets_EM)>1){
      //fleet specific outputs for fishery, panmictic for survey
      report<<"#OBS_survey_fleet"<<endl;
      report<<OBS_survey_total_bio<<endl;
      report<<"#OBS_survey_fleet_bio_se_EM"<<endl;
      report<<OBS_survey_fleet_bio_se_EM<<endl;
      report<<"#OBS_survey_prop"<<endl;
      report<<OBS_survey_prop_pan<<endl;
      report<<"#OBS_survey_prop_N_EM"<<endl;
      report<<OBS_survey_prop_N_EM<<endl;
      report<<"#OBS_yield_fleet"<<endl;
      report<<OBS_yield_fleet<<endl;
      report<<"#OBS_yield_fleet_se_EM"<<endl;
      report<<OBS_yield_fleet_se_EM<<endl;
      report<<"#OBS_catch_prop"<<endl;
      report<<OBS_catch_prop<<endl;
      report<<"#OBS_catch_prop_N_EM"<<endl;
      report<<OBS_catch_prop_N_EM<<endl;
      
//tagging information
      report<<"#nyrs_release"<<endl;
      report<<nyrs_release<<endl;
      report<<"#years_of_tag_releases "<<endl;
      report<<yrs_releases<<endl;
      report<<"#max_life_tags"<<endl;
      report<<max_life_tags<<endl;
      report<<"#age_full_selection"<<endl;
      report<<age_full_selection<<endl;
      report<<"#input_report_rate_EM"<<endl;
      report<<input_report_rate_EM<<endl;
      report<<"#ntags"<<endl;
      report<<ntags_pan<<endl;
      report<<"#ntags_total"<<endl;
      report<<ntags_total<<endl;
      report<<"#tag_N_EM"<<endl;
      report<<tag_N_EM<<endl;
      report<<"#input_T_EM"<<endl;
      report<<input_T_EM<<endl;
      report<<"#OBS_tag_prop_pan_final"<<endl;
      report<<OBS_tag_prop_pan_final<<endl;
     // report<<"#OBS_tag_prop_pan_final_no_age"<<endl;
     // report<<OBS_tag_prop_pan_final_no_age<<endl;
      }
      }

//spatial to spatial EM inputs or panmictic matching - no aggregation needed
     if((EM_structure>0 && OM_structure>0) || (EM_structure==0 && OM_structure==0)){
      report<<"#input_weight"<<endl;
      report<<input_weight<<endl;
      report<<"#input_catch_weight"<<endl;
      report<<input_catch_weight<<endl;
      report<<"#fecundity"<<endl;
      report<<fecundity<<endl;
      report<<"#maturity"<<endl;
      report<<maturity<<endl;
      report<<"#prop_fem"<<endl; 
      report<<prop_fem<<endl;
      report<<"#OBS_rec_index_BM"<<endl;
      report<<rec_index_BM<<endl;
      report<<"#OBS_survey_fleet_bio"<<endl;
      report<<OBS_survey_fleet_bio<<endl;
      report<<"#OBS_survey_fleet_bio_se_EM"<<endl;
      report<<OBS_survey_fleet_bio_se_EM<<endl;
      report<<"#OBS_survey_prop"<<endl;
      report<<OBS_survey_prop<<endl;
      report<<"#OBS_survey_prop_N_EM"<<endl;
      report<<OBS_survey_prop_N_EM<<endl;
      report<<"#OBS_yield_fleet"<<endl;
      report<<OBS_yield_fleet<<endl;
      report<<"#OBS_yield_fleet_se_EM"<<endl;
      report<<OBS_yield_fleet_se_EM<<endl;
      report<<"#OBS_catch_prop"<<endl;
      report<<OBS_catch_prop<<endl;
      report<<"#OBS_catch_prop_N_EM"<<endl;
      report<<OBS_catch_prop_N_EM<<endl;

//tagging information
      report<<"#nyrs_release"<<endl;
      report<<nyrs_release<<endl;
      report<<"#years_of_tag_releases "<<endl;
      report<<yrs_releases<<endl;
      report<<"#max_life_tags"<<endl;
      report<<max_life_tags<<endl;
      report<<"#age_full_selection"<<endl;
      report<<age_full_selection<<endl;
      report<<"#input_report_rate_EM"<<endl;
      report<<input_report_rate_EM<<endl;
      report<<"#ntags"<<endl;
      report<<ntags<<endl;
      report<<"#ntags_total"<<endl;
      report<<ntags_total<<endl;
      report<<"#tag_N_EM"<<endl;
      report<<tag_N_EM<<endl;
      report<<"#input_T_EM"<<endl;
      report<<input_T_EM<<endl;
      report<<"#OBS_tag_prop_final"<<endl;
      report<<OBS_tag_prop_final<<endl;
      report<<"#OBS_tag_prop_final_no_age"<<endl;
      report<<OBS_tag_prop_final_no_age<<endl;
      }


///////////
// panmictic to spatial will go here eventually
/////////

//Additional inputs for EM specified in OM .dat
  report<<"#input_M_EM"<<endl;
  report<<input_M_EM<<endl;
  report<<"#input_R_ave_EM"<<endl;
  report<<input_R_ave_EM<<endl;
  report<<"#input_Rec_Prop_EM"<<endl;
  report<<input_rec_prop_EM<<endl;
  report<<"#input_selectivity_EM"<<endl;
  report<<input_selectivity_EM<<endl;
  report<<"#input_survey_selectivity_EM"<<endl;
  report<<input_survey_selectivity_EM<<endl;
  report<<"#input_steep_EM"<<endl;
  report<<input_steep_EM<<endl;
  report<<"#input_dist_init_abund"<<endl;
  report<<input_dist_init_abund_EM<<endl;

/// TRUE VALUES FROM OM
  report<<"#frac_natal_true"<<endl;
  report<<frac_natal_true<<endl;  
  report<<"#input_M_TRUE"<<endl;
  report<<input_M_TRUE<<endl;
  report<<"#init_abund_TRUE"<<endl;
  report<<init_abund<<endl;
  report<<"#q_survey"<<endl;
  report<<q_survey<<endl;
  report<<"#sel_beta1"<<endl;
  report<<sel_beta1<<endl;
  report<<"#sel_beta2"<<endl;
  report<<sel_beta2<<endl;
  report<<"#sel_beta3"<<endl;
  report<<sel_beta3<<endl;
  report<<"#sel_beta4"<<endl;
  report<<sel_beta4<<endl;
  report<<"#sel_beta1_survey"<<endl;
  report<<sel_beta1_survey<<endl;
  report<<"#sel_beta2_survey"<<endl;
  report<<sel_beta2_survey<<endl;
  report<<"#sel_beta3_survey"<<endl;
  report<<sel_beta3_survey<<endl;
  report<<"#sel_beta4_survey"<<endl;
  report<<sel_beta4_survey<<endl;
  report<<"#steep"<<endl;
  report<<steep<<endl;
  report<<"#R_ave"<<endl;
  report<<R_ave<<endl;
  report<<"#SSB_zero"<<endl;
  report<<SSB_zero<<endl;
  report<<"#rec_devs"<<endl;
  report<<rec_devs<<endl;
  report<<"#Rec_Prop"<<endl;
  report<<Rec_Prop<<endl;
  report<<"#recruits_BM"<<endl;
  report<<recruits_BM<<endl;
  report<<"#F"<<endl;
  report<<F<<endl;
  report<<"#Fyear"<<endl;
  report<<F_year<<endl;
  report<<"#biomass_AM"<<endl;
  report<<biomass_AM<<endl;
  report<<"#biomass_population"<<endl;
  report<<biomass_population<<endl;

//add to EM input
  report<<"#catch_at_age_fleet_prop"<<endl;
  report<<catch_at_age_fleet_prop<<endl;
  report<<"#yield_fleet"<<endl;
  report<<yield_fleet<<endl;
  report<<"#survey_at_age_fleet_prop"<<endl;
  report<<survey_at_age_fleet_prop<<endl;
  report<<"#true_survey_fleet_bio"<<endl;
  report<<true_survey_fleet_bio<<endl;
//

  report<<"#harvest_rate_region_bio"<<endl;
  report<<harvest_rate_region_bio<<endl;
  report<<"#depletion_region"<<endl;
  report<<depletion_region<<endl;
  report<<"#SSB_region"<<endl;
  report<<SSB_region<<endl;
  report<<"#Bratio_population"<<endl;
  report<<Bratio_population<<endl;
  report<<"#T_year"<<endl;
  report<<T_year<<endl;

//reporting true aggregated selectivity if fleets ==1
  if(EM_structure==0 && sum(nfleets_EM)==1)
  {
  report<<"#selectivity_age"<<endl;
  report<<selectivity_age_pop<<endl;
  report<<"#survey_selectivity_age"<<endl;
  report<<survey_selectivity_age_pop<<endl;
  }
  
  else{
  report<<"#selectivity_age"<<endl;
  report<<selectivity_age<<endl;
  report<<"#survey_selectivity_age"<<endl;
  report<<survey_selectivity_age<<endl;
  }
  
 //true tag information
  report<<"#TRUE_tag_prop"<<endl;
  report<<tag_prop_final<<endl;
  report<<"#TRUE_tag_prop_no_age"<<endl;
  report<<tag_prop_final_no_age<<endl;
  
  report<<"#T"<<endl;
  report<<T<<endl;
  report<<"#report_rate_TRUE"<<endl;
  report<<report_rate_TRUE<<endl;
  
  report<<"#move_switch_OM"<<endl;
  report<<move_switch<<endl;
  report<<"#DD_move_age_switch_OM"<<endl;
  report<<DD_move_age_switch<<endl;
  
//report the abundance_frac for later calcs if needed
  report<<"#abund_frac_age_region"<<endl;
  report<<abund_frac_age_region<<endl;
  report<<"#abund_frac_year"<<endl;
  report<<abund_frac_region_year<<endl;
  report<<"#abund_frac_region"<<endl;
  report<<abund_frac_region<<endl;

  report<<"#F_tag_scalar"<<endl;
  report<<F_tag_scalar<<endl;
  report<<"#T_tag_res"<<endl;
  report<<T_tag_res<<endl;
  report<<"#sim_tag_mixing_switch"<<endl;
  report<<sim_tag_mixing_switch<<endl;
  
  report<<"#debug"<<endl;
  report<<debug<<endl;


//to get abundance after movement for init abund
// report<<"#Abund_BM"<<endl;
// report<<abundance_at_age_BM<<endl;
// report<<"#Abund_AM"<<endl;
// report<<abundance_at_age_AM<<endl;

 // Some stuff for looking at tag calculations 
  /*
  report<<"#OBS_tag_prop_final"<<endl;
  report<<OBS_tag_prop_final(1,1,1)<<endl;
  report<<"nt"<<endl;
  report<<"true_survey_population_bio"<<true_survey_population_bio<<endl;
  report<<"true_survey_region_bio"<<true_survey_region_bio<<endl;
  report<<"abundance_total"<<abundance_total<<endl<<"abundance_population"<<abundance_population<<endl;
  report<<"true fleet_age"<<true_survey_fleet_age<<endl;
  report<<"survey_abundance"<<endl<<true_survey_region_abundance<<endl;
  report<<"temp_fleet_age"<<endl<<true_survey_fleet_age_temp<<endl;


  report<<"ntags(1,1,1,1)"<<endl<<ntags(1,1,1,9)<<endl<<"ntags_total(x)"<<endl<<ntags_total(1)<<endl;
  report<<"true_survey_population_abundance(xx,i,a)"<<endl<<true_survey_population_abundance(1,1,9)<<endl;
  report<<"true_survey_total_abundance(xx,a)"<<endl<<sum(true_survey_total_abundance(1))<<endl;
  report<<"true_survey_region_abundance(i,xx,n,a)"<<endl<<true_survey_region_abundance(1,1,1,9)<<endl;
  report<<"true_survey_population_abundance(xx,i,a)"<<endl<<true_survey_population_abundance(1,1,9)<<endl;
  report<<"survey_selectivity(i,n,xx,a,1)"<<endl<<survey_selectivity(1,1,1,9,1)<<endl;
  report<<"sum(survey_selectivity_temp(i,n,xx,1)))"<<endl<<sum(survey_selectivity_temp(1,1,1,1))<<endl;
 */

RUNTIME_SECTION
  convergence_criteria .001,.0001, 1.0e-4, 1.0e-7
  maximum_function_evaluations 100000
  


