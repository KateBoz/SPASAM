
/////////////////////////////////////////////////////////
// Spatial Estimation Model based on Operating Model by Daniel Goethel (NMFS SEFSC)  
// started by Dana Hanselman (AFSC) //Further screwed up by Jon Deroba (NEFSC)//Terrorized by Dan Goethel (SEFSC)
// Agonized over again by Dana Hanselman between bloodlettings...
//////////////////////////////////////////////////////////

GLOBALS_SECTION
  #include "admodel.h"
  #include <contrib.h>
  #define EOUT(var) cout <<#var<<" "<<var<<endl;

TOP_OF_MAIN_SECTION
  arrmblsize=500000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);

DATA_SECTION
////////////////////////////////////////////////////////////////////////////////////
/////MODEL STRUCTURE INPUTS/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
  init_int nages //number of ages
  init_int nyrs //number of years of catch
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
  ivector nfs(1,np)
  !! nfs=nfleets_survey;

////////////////////////////////////////////////////////////////////////////////////
//////////////SWITCHES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  init_matrix tsurvey(1,np,1,nreg) //time of survey in proportion of year (0,1)
  init_number larval_move_switch
  ///// Changes the type of larval movement pattern (sets age class 1 movements)
  //==0 no movement
  //==1 input movement
  //==2 movement within population only based on residency (symmetric)
  //==3 symmetric movement but only allow movement within a population (ie regions within a population) not across populations
  //==4 symmetric movement across all populations and regions
  //==5 allow movement across all regions and populations, based on population/region specific residency (symmetric off-diag)

  init_number move_switch
  ///// Sets the type of adult movement pattern (sets age class>1 movements)
  //==0 no movement, set T phases=-1
  //==1 input movement, set T phases=-1
  //==2 movement among populations, set phase_T_pop>1...est yearly movement rates
  //==3 movement among regions (only 1 pop), set phase_T_reg>1...est yearly movement rates
  //==4 time-invariant movement multiple pops, set phase_T_CNST_pop>1
  //==5 time-invariant movement multiple regs, set phase_T_CNST_reg>1
  
////// Population Structure switches
  init_number natal_homing_switch
  //==0 no natal homing (SSB is sum of SSB in population regardless of natal origin; weight/mat/fecund/ are based on current population not natal population) - Metapopulation/metamictic
  //==1 do natal homing (a fish only adds to SSB if it is in its natal population at spawning time; weight/mat/fecund/ are based on natal population) - Natal homing
  //natal homing  assumes genetic based life history and contribution to SSB (i.e., natal homing and no demographic mixing), natal homing==0 assumes demographic mixing (e.g. metapopulations where life history is more location based)

  init_number spawn_return_switch
   //==0 if natal_homing_switch==1 then only fish that are in natal population add to SSB
   //==1 natal_homing_switch==1 a fraction of fish return to natal population to spawn (inpopsantaneous migration to natal population and back at time of spawning) based spawn_return_prob; weight/mat/fecund/ are based on natal population)
//////////////////////////////////////////////////////

  init_number select_switch
  //==0 input selectivity
  //==1 logistic selectivity based on input sel_beta1 and sel_beta2
  //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4
/////////////////////////////////////////////////////
init_number select_switch_survey

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
  
  init_number Rec_type
  //==1 stock-recruit relationship assumes an average value based on R_ave
  //==2 Beverton-Holt population-recruit functions based on population-specific input steepness, R0 (R_ave), M, and weight

  init_number apportionment_type
  //==-1 no recruitment apportionment to regions within a population (each region within a population gets full amount of recruits from SR curve)
  //==0 apportionment to each region is based on relative SSB in region compared to population SSB
  //==1 input apportionment
  //==2 recruits are apportioned equally to each region within a population
  //==3 estimate pop/reg apportionment
  //==4 estimate pop/reg/year apportionment

  init_number use_stock_comp_info_survey //for likelihood calcs
  //Determines whether it is assumed that info (stock composition data) is available determine natal origin for age composition data
  //==0 calc survey age comps by area (summed across natal population)
  //==1 calc survey age comps by natal population within each area

  init_number use_stock_comp_info_catch //for likelihood calcss
  //Determines whether it is assumed that info (stock composition data) is available determine natal origin for age composition data
  //==0 calc  catch age comps by area (summed across natal population)
  //==1 calc  catch age comps by natal population within each area
   init_number F_switch
  //negative phase == use input_F_TRUE
  //==1 estimate yearly F
  //==2 random walk in F
 init_number recruit_devs_switch
  //==0 use stock-recruit relationphip directly (make sure to set ph_rec=0), also assumes initial abund for all ages=R0
  //==1 allow lognormal error around SR curve (i.e., include randomness based on input sigma_recruit)


   init_number recruit_randwalk_switch
  //==0 no random walk recruitment deviations
  //==1 have random walk lognormal recruitment deviations (requirs recruit_devs_switch==1)....NEEDS WORK!!!!!
  
  init_vector tspawn(1,np) //time of spawning in proportion of year (0-1)
  init_number return_age // used if move_swith ==6
  init_vector return_probability(1,np) // used if move_swith==6
  init_vector spawn_return_prob(1,np) // used if natal_homing_swith==2
  init_int do_tag
  init_int do_tag_mult //if==0 assume neg binomial, if==1 assume multinomial (same as OM)
  init_vector sigma_recruit(1,np)

///////////////////////////////////////////////////////////////////////////////
////////READ IN THE SPECS  //////////////////////////////////
//////////////////////////////////
/////////////////////////////////////////////
//phases and bounds for est parameters

  init_int ph_lmr
  init_int ph_rec
  init_int ph_abund_devs
  init_int ph_F
  init_int ph_steep
  init_int ph_M
  init_int ph_sel_log
  init_number lb_sel_beta1
  init_number ub_sel_beta1
  init_number lb_sel_beta2
  init_number ub_sel_beta2
  init_int ph_sel_log_surv
  init_int ph_sel_dubl
  init_int ph_sel_dubl_surv
  init_int ph_q
  init_int ph_F_rho // if we want random walk F
  init_int phase_T_pop //use if mult pops
  init_int phase_T_reg //use if mult regs
  init_int phase_T_CNST_pop //use if mult pops
  init_int phase_T_CNST_reg //use if mult regs
  init_int ph_dummy
  number ph_theta
 // !!if(do_tag_mult==0)
 // !! {
 // !!  ph_theta==2
 // !! }
 // !!if(do_tag_mult==2)
//  !! {
//  !!  ph_theta==-2
//  !! }
 // likleihood weights
   init_number wt_srv
   init_number wt_catch
   init_number wt_fish_age
   init_number wt_srv_age 
   init_number wt_rec
   init_number wt_tag
   init_int abund_pen_switch
   init_int move_pen_switch
   init_number Tpen
   init_number Tpen2
//###########READ BIO DATA###############################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################
//#########################################################################################################################################
//######### FOR NATAL HOMING THERE IS NO ACCOUNTING OF REGIONAL DIFFERENCES IN VITAL RATES ACROSS REGIONS WITHIN A POPULATION
//######### IE BECAUSE GENETICS DEFINE VITAL RATES, THEY MUST ALL BE THE SAME
//######### **********DO NOT INPUT REGIONALLY VARYING VITAL RATES, NATAL REGION WILL NOT BE PROPERLY TRACKED IN SSB CALCS #############
//#########################################################################################################################################
  ////////////////BIOLOGICAL PARAMETERS////////////////
  /////////////////////////////////////////////////////
  init_matrix input_Rec_prop(1,np,1,nreg)
  init_3darray input_weight(1,np,1,nreg,1,na)  
  init_3darray input_catch_weight(1,np,1,nreg,1,na)
  init_3darray fecundity(1,np,1,nreg,1,na)
  init_3darray maturity(1,np,1,nreg,1,na)
  init_matrix prop_fem(1,np,1,nreg) //proportion of population assumed to be female for SSB calcs (typically use 0.5)
  //##########################################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################
//#########################################################################################################################################

///////////////////////////////////////////////////////////////////////////////////////READ IN THE OBS DATA  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//Recruit indices, I made one up, but didn't fit, would be the same as survey -lnL
  init_3darray OBS_rec_index_BM(1,np,1,nreg,1,ny)

// 3darray OBS_rec_index_AM(1,npops,1,nr,1,nyr)

//Survey data and age comps
  init_4darray OBS_survey_fleet_bio(1,np,1,nreg,1,ny,1,nfs)  
  init_4darray OBS_survey_fleet_bio_se(1,np,1,nreg,1,ny,1,nfs) //survey standard erros
  init_5darray OBS_survey_prop(1,np,1,nreg,1,ny,1,nfs,1,na)  
  init_4darray OBS_survey_prop_N(1,np,1,nreg,1,ny,1,nfs)   //sample size

//Catch and age compos
 init_4darray OBS_yield_fleet(1,np,1,nreg,1,ny,1,nf)
 init_4darray OBS_yield_fleet_se(1,np,1,nreg,1,ny,1,nf)  // standard error on catch
 init_5darray OBS_catch_at_age_fleet_prop(1,np,1,nreg,1,ny,1,nf,1,na) 
 init_4darray OBS_catch_at_age_fleet_prop_N(1,np,1,nreg,1,ny,1,nf) //sample size

//Catch Prop
//tagging data parameters
  init_int nyrs_release //number of years with tag release events  
  !! int ny_rel=nyrs_release;
  init_vector yrs_releases(1,ny_rel) //vector containing the model years with releases 
  init_int max_life_tags //number of years that tag recaptures will be tallied for after release (assume proportional to longevity of the species)...use this to avoid calculating tag recaptures for all remaining model years after release since # recaptures are often extremely limited after a few years after release
    !! int tag_age=max_life_tags;
  init_3darray report_rate(1,np,1,ny_rel,1,nreg) //tag reporting rate (assume constant for all recaptures within a given release cohort, but can be variable across populations or regions)...could switch to allow variation across fleets instead
  init_4darray ntags(1,np,1,nreg,1,ny_rel,1,na) //releases
  init_vector ntags_total(1,ny_rel)  
  init_4darray OBS_tag_prop_N(1,np,1,nreg,1,ny_rel,1,na) //eff_N for tag mult tag_prop
     !! int recap_index=np*nreg(1)*ny_rel; // using this to dimension down arrays
 //  init_5darray OBS_recaps_temp(1,recap_index,1,na,1,tag_age,1,np,1,3) // not sure how to have 1,nreg here
   init_5darray input_T(1,np,1,nreg,1,na,1,np,1,nreg)

  !! int nyr_rel=nyrs_release;
  !! ivector xy(1,nyr_rel);
  !! ivector nt(1,nyr_rel);

  !!  for(int x=1; x<=nyrs_release; x++)
  !!   {
  !!    xx=yrs_releases(x);
  !!    xy(x)=min(max_life_tags,nyrs-xx+1);
  !!    nt(x)=xy(x)*sum(nregions)+1;
  !!   }
   init_5darray OBS_tag_prop_final(1,np,1,nreg,1,nyr_rel,1,na,1,nt)
   matrix input_residency_larval(1,np,1,nreg)  //larval residency probability
   3darray input_residency(1,np,1,nreg,1,na) //
// probably need to calculate quantity below  or omit *dh
   vector frac_total_abund_tagged(1,ny_rel) //proportion of total abundance that is tagged in each 
 //  7darray OBS_recaps(1,np,1,nreg,1,ny_rel,1,na,1,tag_age,1,np,1,nreg) // for filling for calcs later

   init_4darray init_abund_TRUE(1,np,1,np,1,nreg,1,na);  //input true initial abundance; just used for reporting
   init_matrix input_M(1,np,1,na); // input for now, if we estimate we will want to limit how many Ms
  //##########################################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################
 // TRUE VALUES
  //##########################################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################
  init_3darray q_survey_TRUE(1,np,1,nreg,1,nfs) // catchability for different surveys(fleets)operating in different areas
  init_3darray sel_beta1_TRUE(1,np,1,nreg,1,nf)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta2_TRUE(1,np,1,nreg,1,nf)   //selectivity inflection parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta3_TRUE(1,np,1,nreg,1,nf)  //selectivity slope parameter 2 for double selectivity
  init_3darray sel_beta4_TRUE(1,np,1,nreg,1,nf)  //selectivity inflection parameter 2 for double logistic selectivity
  init_3darray sel_beta1_survey_TRUE(1,np,1,nreg,1,nfs)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta2_survey_TRUE(1,np,1,nreg,1,nfs)   //selectivity inflection parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta3_survey_TRUE(1,np,1,nreg,1,nfs)  //selectivity slope parameter 2 for double selectivity
  init_3darray sel_beta4_survey_TRUE(1,np,1,nreg,1,nfs)  //selectivity inflection parameter 2 for double logistic selectivity
  init_vector steep_TRUE(1,np) //B-H steepness
  init_vector R_ave_TRUE(1,np) //Average Recruitment or R0 for B-H S-R curve
  init_vector SSB_zero_TRUE(1,np)
  init_matrix rec_devs_TRUE(1,np,1,ny)
  init_3darray Rec_Prop_TRUE(1,np,1,nreg,1,ny)
  init_3darray recruits_BM_TRUE(1,np,1,nreg,1,ny)
  init_4darray F_TRUE(1,np,1,nreg,1,ny,1,na)
  init_4darray F_year_TRUE(1,np,1,nreg,1,ny,1,nfs)
  init_3darray biomass_AM_TRUE(1,np,1,nreg,1,ny)
  init_matrix biomass_population_TRUE(1,np,1,ny)
  init_3darray harvest_rate_region_bio_TRUE(1,np,1,nreg,1,ny)
  init_3darray depletion_region_TRUE(1,np,1,nreg,1,ny)
  init_3darray SSB_region_TRUE(1,np,1,nreg,1,ny)
  init_matrix Bratio_population_TRUE(1,np,1,ny)
  init_5darray T_year_TRUE(1,np,1,nreg,1,ny,1,np,1,nreg)
  init_4darray selectivity_age_TRUE(1,np,1,nreg,1,na,1,nf)
  init_4darray survey_selectivity_age_TRUE(1,np,1,nreg,1,na,1,nfs)
  //##########################################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################

// end of file marker
   init_int debug;
  //##########################################################################################################################################
//#########################################################################################################################################
//##########################################################################################################################################
//////////////////////////////////////
// END DATA SECTION /////////////////////
///////////////////////////////////////
  //fill in a vector of years
  // Not sure what this shit does, but doesn't hurt so I left it
  vector years(1,nyrs)
  !!years.fill_seqadd(double(1),1.0);
  
    imatrix nregions_temp(1,np,1,np) //used to fill tag_recap matrices

  !! for(int j=1;j<=np;j++) //recap stock
  !! {
  !!  for (int r=1;r<=np;r++) //recap region
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

  int region_counter
 
 !! cout << "debug = " << debug << endl;
 !! cout << "If debug != 1541 then .dat file is meh" << endl;
 !! cout << "input read" << endl;

 !! cout << "end setup of containers" << endl;


PARAMETER_SECTION
 !! cout << "begin parameter section" << endl;

  !! ivector nr=nregions;
  !! int parpops=npops;
  !! int nps=npops;
  !! int nyr=nyrs;
  !! int nag=nages;
  !! ivector nfl=nfleets;
  !! ivector nfls=nfleets_survey;  

  !! int k;
// need these for some if the init parameter arrays below
  !! int fishfleet=nfleets(1); 
  !! int survfleet=nfleets_survey(1);
  !! int parreg=nregions(1);

 //movement paramters

 //##################################################################################################
 // NO IDEA HOW TO MAKE THESE ESTIMATED OBJECTS WHEN THERE IS ONLY 1 POPULATION
 // THE SAME IS TRUE FOR T_REG IF NUMBER OF REGIONS IS 1
 // DECIDED TO JUST MAKE 2 EM, ONE FOR ESTIMATING MOVEMENT BETWEEN POPS, 1 FOR EST MOVEMENT BETWEEN REGIONS
 // CAN'T HAVE MOVEMENT BETWEEN MULTIPLE REGIONS AND MULTIPLE POPULATIONS BECAUSE MY HEAD WOULD EXPLODE AND THE UNIVERSE WOULD IMPLODE
 // DG
 //#########################################################################################################
 //!!if(nps>1)
 // !!{
     init_3darray ln_T_est_pop(1,nyrs,1,nps,1,nps,phase_T_pop);  //DO NOT USE
     init_matrix ln_T_CNST_pop(1,nps,1,nps,phase_T_CNST_pop);  //DO NOT USE
 // !!}
 ///#####################################################################################################
   matrix G_pop(1,nps,1,nps);
   vector G_temp_pop(1,nps);

   init_3darray ln_T_est_reg(1,nyrs,1,parreg,1,parreg-1,phase_T_reg);  //movement parameters time-varying
   matrix G_reg(1,parreg,1,parreg);
   vector G_temp_reg(1,parreg);

   init_matrix ln_T_CNST_reg(1,parreg,1,parreg-1,phase_T_CNST_reg);  //movement parameters, time-invariant
  
  // selectivity parameters

   !! int sel_lgth=parpops*sum(nr);

   init_bounded_matrix log_sel_beta1(1,sel_lgth,1,fishfleet,lb_sel_beta1,ub_sel_beta1,ph_sel_log);   //selectivity slope parameter 1 for logistic selectivity/double logistic
   init_bounded_matrix log_sel_beta2(1,sel_lgth,1,fishfleet,lb_sel_beta2,ub_sel_beta2,ph_sel_log);   //selectivity inflection parameter 1 for logistic selectivity/double logistic
   init_bounded_matrix log_sel_beta3(1,sel_lgth,1,fishfleet,-10,5,ph_sel_dubl);  //selectivity slope parameter 2 for double selectivity
   init_bounded_matrix log_sel_beta4(1,sel_lgth,1,fishfleet,-10,5,ph_sel_dubl);//selectivity inflection parameter 2 for double logistic selectivity
   init_bounded_matrix log_sel_beta1surv(1,parpops,1,survfleet,-10,5,ph_sel_log_surv);   //selectivity slope parameter 1 for logistic selectivity/double logistic
   init_bounded_matrix log_sel_beta2surv(1,parpops,1,survfleet,-10,5,ph_sel_log_surv) ;  //selectivity inflection parameter 1 for logistic selectivity/double logistic
   init_bounded_matrix log_sel_beta3surv(1,parpops,1,survfleet,-10,5,ph_sel_dubl_surv);   //selectivity slope parameter 1 for logistic selectivity/double logistic
   init_bounded_matrix log_sel_beta4surv(1,parpops,1,survfleet,-10,5,ph_sel_dubl_surv) ;  //selectivity inflection parameter 1 for logistic selectivity/double logistic
  

  3darray sel_beta1(1,nps,1,nr,1,nfl)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  3darray sel_beta2(1,nps,1,nr,1,nfl)   //selectivity inflection parameter 1 for logistic selectivity/double logistic
  3darray sel_beta3(1,nps,1,nr,1,nfl)  //selectivity slope parameter 2 for double selectivity
  3darray sel_beta4(1,nps,1,nr,1,nfl)  //selectivity inflection parameter 2 for double logistic selectivity
  3darray sel_beta1surv(1,nps,1,nr,1,nfls)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  3darray sel_beta2surv(1,nps,1,nr,1,nfls)  //selectivity inflection parameter 1 for logistic selectivity/double logistic
  3darray sel_beta3surv(1,nps,1,nr,1,nfls)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  3darray sel_beta4surv(1,nps,1,nr,1,nfls)  //selectivity inflection parameter 1 for logistic selectivity/double logistic


  5darray survey_selectivity(1,nps,1,nr,1,nyr,1,nag,1,nfls)  //param
  5darray selectivity(1,nps,1,nr,1,nyr,1,nag,1,nfl)
  4darray survey_selectivity_age(1,nps,1,nr,1,nag,1,nfls)  //param
  4darray selectivity_age(1,nps,1,nr,1,nag,1,nfl)
  
  init_bounded_matrix ln_q(1,parpops,1,survfleet,-30,15,ph_q)
 
  3darray q_survey(1,parpops,1,nr,1,nfls)  //
 //###########WHY HAVE Q estimated by pop by applied by region?  can't q_survey just be a matrix?
//F parameters
   !! int F_lgth=parpops*sum(nr)*nyr;

  init_bounded_matrix ln_F(1,F_lgth,1,fishfleet,-7,2,ph_F) //the actual parameters  
  init_bounded_matrix F_rho(1,parpops,1,fishfleet,0,2,ph_F_rho) //random walk params*
  5darray F_fleet(1,nps,1,nr,1,nyr,1,nag,1,nfl) //derived quantity in which we likely have interest in precision
  4darray F_year(1,nps,1,nr,1,nyr,1,nfl) //derived quantity in which we likely have interest in precision
  4darray F(1,nps,1,nr,1,nyr,1,nag) //derived quantity in which we likely have interest in precision
  4darray M(1,nps,1,nr,1,nyr,1,nag)
 //tagging data and parameters

  !! int nyr_rel=nyrs_release;
  !! ivector xy(1,nyr_rel);
  !! ivector nt(1,nyr_rel);
  !! ivector nt2(1,nyr_rel);
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
  4darray total_rec(1,nps,1,nr,1,nyr_rel,1,nag)
  4darray not_rec(1,nps,1,nr,1,nyr_rel,1,nag)
  7darray tag_prop(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr)
  4darray tag_prop_not_rec(1,nps,1,nr,1,nyr_rel,1,nag)
  5darray tag_prop_final(1,nps,1,nr,1,nyr_rel,1,nag,1,nt)
  6darray T(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr) 
  5darray T_terminal(1,nps,1,nr,1,nag,1,nps,1,nr) 
  5darray T_year(1,nps,1,nr,1,nyr,1,nps,1,nr) 
  7darray recaps(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr) //recaps
 
//recruitment parameters
  init_bounded_matrix ln_rec_prop(1,parpops,1,parreg,-4,0,-1)  // do we need this?
  //is phase suppose to be negative below?? --DG
  // for now... not sure what these are used for in estimation model? Are they necessary?
  //YES..if we want to estimate recruit apportionment (typical in SS)..need to implement logit transform in later code
  init_3darray ln_rec_prop_year(1,parpops,1,parreg,1,nyr,-1) 
 
 //#########################################################################################################################
 //## rec devs and initial devs for abundance...this approach will be problematic if we don't want to estimate recruit devs at any point (or fix them)
 //## because if rec_devs_switch==0 we set them to 1 so will need to turn off both rec_devs and initial_abund_devs
 //############################################################################################################################
  //init_bounded_matrix ln_rec_devs_RN(1,parpops,1,nyr+nages-1,-40,40,ph_rec) //actual parameters (log scale devs)
   !! int dev_lgth=nps*(nyr-1);
  init_bounded_dev_vector ln_rec_devs_RN(1,dev_lgth,-40,40,ph_rec)
  init_bounded_matrix ln_abund_devs(1,nps,1,nag,-10,10,ph_abund_devs)
 //###############################################################################################################################
 //###################################################################################################################################
  matrix abund_devs(1,nps,1,nages)
  matrix rec_devs(1,nps,1,nyr-1) //derived quantity as exp(rec_devs_RN)
  init_bounded_vector ln_R_ave(1,parpops,2,20,ph_lmr) //estimated parameter Average Recruitment
  vector R_ave(1,parpops) // switch to log scale
  vector SSB_zero(1,nps) //derived quantity
  init_bounded_vector steep(1,parpops,0.2,1,ph_steep) //B-H steepness //could be estimated parameter or input value
  vector alpha(1,nps) //derived quantity
  vector beta(1,nps) //derived quantity
  matrix SR(1,nps,1,nyr-1)
  matrix total_recruits(1,nps,1,nyr-1)
//end recruitment parameters


 4darray init_abund(1,nps,1,nps,1,nr,1,nag)  //estimated as devs from exponential decline from R_ave
 
 // 6-d arrays
 6darray survey_fleet_overlap_age(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag) 
 6darray survey_at_age_region_fleet_overlap_prop(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag)
 6darray survey_fleet_overlap_age_bio(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag)
 6darray catch_at_age_region_fleet_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray catch_at_age_region_fleet_overlap_prop(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray yield_region_fleet_temp_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)

 4darray weight_population(1,nps,1,nr,1,nyr,1,nag)
 4darray weight_catch(1,nps,1,nr,1,nyr,1,nag)
 3darray wt_mat_mult(1,nps,1,nyr,1,nag)
 4darray wt_mat_mult_reg(1,nps,1,nr,1,nyr,1,nag)
 3darray ave_mat_temp(1,nps,1,nag,1,nr) //to calc average maturity 
 matrix ave_mat(1,nps,1,nag) //to calc average maturity
 matrix SPR_N(1,nps,1,nag) 
 matrix SPR_SSB(1,nps,1,nag) 
 vector SPR(1,nps) 
 
//recruitment 
 3darray recruits_BM(1,nps,1,nr,1,nyr) //param
 3darray recruits_AM(1,nps,1,nr,1,nyr) //param
 3darray rec_index_BM(1,nps,1,nr,1,nyr)
 3darray rec_index_AM(1,nps,1,nr,1,nyr)
 3darray rec_index_prop_BM(1,nps,1,nr,1,nyr) //
 3darray rec_index_BM_temp(1,nps,1,nyr,1,nr) //
 3darray rec_index_prop_AM(1,nps,1,nr,1,nyr) //
 3darray rec_index_AM_temp(1,nps,1,nyr,1,nr) //
 matrix rec_devs_randwalk(1,nps,1,nyr-1)

 3darray Rec_Prop(1,nps,1,nr,1,nyr)
 3darray Rec_prop_temp1(1,nps,1,nr,1,nyr)
 3darray Rec_prop_temp2(1,nps,1,nr,1,nyr)
 vector env_rec(1,nyr)

//abundance 
 4darray abundance_at_age_BM(1,nps,1,nr,1,nyr,1,nag) //param
 4darray abundance_at_age_AM(1,nps,1,nr,1,nyr,1,nag) //param
 4darray abundance_in(1,nps,1,nr,1,nyr,1,nag) //param
 4darray abundance_res(1,nps,1,nr,1,nyr,1,nag) //param
 4darray abundance_leave(1,nps,1,nr,1,nyr,1,nag) //param
 4darray abundance_spawn(1,nps,1,nr,1,nyr,1,nag) //param

//biomass
 4darray biomass_BM_age(1,nps,1,nr,1,nyr,1,nag) //param
 4darray biomass_AM_age(1,nps,1,nr,1,nyr,1,nag) //param
 3darray biomass_BM(1,nps,1,nr,1,nyr) //param
 3darray biomass_AM(1,nps,1,nr,1,nyr) //param
 4darray bio_in(1,nps,1,nr,1,nyr,1,nag) //param
 4darray bio_res(1,nps,1,nr,1,nyr,1,nag) //param
 4darray bio_leave(1,nps,1,nr,1,nyr,1,nag) //param

 //yield & BRP calcs 
 5darray catch_at_age_fleet(1,nps,1,nr,1,nyr,1,nag,1,nfl)
 5darray catch_at_age_fleet_prop(1,nps,1,nr,1,nyr,1,nfl,1,nag) 
 4darray yield_fleet(1,nps,1,nr,1,nyr,1,nfl)
 4darray catch_at_age_region(1,nps,1,nr,1,nyr,1,nag)
 4darray catch_at_age_region_prop(1,nps,1,nr,1,nyr,1,nag)
 3darray yield_region(1,nps,1,nr,1,nyr)
 3darray catch_at_age_population(1,nps,1,nyr,1,nag)
 3darray catch_at_age_population_prop(1,nps,1,nyr,1,nag)
 matrix yield_population(1,nps,1,nyr)
 3darray SSB_region(1,nps,1,nr,1,nyr) //param
 matrix SSB_population(1,nps,1,nyr) //param
 vector SSB_total(1,nyr) //param
 3darray abundance_population(1,nps,1,nyr,1,nag) //param
 matrix abundance_total(1,nyr,1,nag) //param
 matrix biomass_population(1,nps,1,nyr) //param
 vector biomass_total(1,nyr) //param
 matrix catch_at_age_total(1,nyr,1,nag) //
 matrix catch_at_age_total_prop(1,nyr,1,nag) //
 vector yield_total(1,nyr) //
 4darray harvest_rate_region_num(1,nps,1,nr,1,nyr,1,nag) //
 3darray harvest_rate_population_num(1,nps,1,nyr,1,nag) //
 matrix harvest_rate_total_num(1,nyr,1,nag) //
 3darray harvest_rate_region_bio(1,nps,1,nr,1,nyr) //
 matrix harvest_rate_population_bio(1,nps,1,nyr) //
 vector harvest_rate_total_bio(1,nyr) //
 3darray depletion_region(1,nps,1,nr,1,nyr) //
 matrix depletion_population(1,nps,1,nyr) //
 vector depletion_total(1,nyr) //

 5darray abundance_at_age_BM_overlap_region(1,nps,1,nps,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_BM_overlap_population(1,nps,1,nps,1,nyr,1,nag)
 5darray abundance_at_age_AM_overlap_region(1,nps,1,nps,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_AM_overlap_population(1,nps,1,nps,1,nyr,1,nag)
 4darray abundance_AM_overlap_region_all_natal(1,nps,1,nr,1,nyr,1,nag)
 5darray abundance_spawn_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray SSB_region_overlap(1,nps,1,nps,1,nr,1,nyr)
 3darray SSB_population_overlap(1,nps,1,nps,1,nyr)
 matrix SSB_natal_overlap(1,nps,1,nyr)

 5darray biomass_BM_overlap_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 4darray init_abund_temp(1,nps,1,nr,1,nag,1,nps)
 5darray survey_fleet_bio_overlap_temp(1,nps,1,nr,1,nyr,1,nfls,1,nps)
 5darray catch_at_age_fleet_prop_temp(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 matrix abundance_move_temp(1,nps,1,nr)
 matrix bio_move_temp(1,nps,1,nr)
 matrix abundance_move_overlap_temp(1,nps,1,nr)
 5darray abundance_AM_overlap_region_all_natal_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 5darray biomass_AM_overlap_region_all_natal_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 3darray SSB_natal_overlap_temp(1,nps,1,nyr,1,nps)
 4darray abundance_natal_temp_overlap(1,nps,1,nyr,1,nag,1,nps)
 4darray SSB_population_temp_overlap(1,nps,1,nps,1,nyr,1,nr)
 5darray SSB_region_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)

//survey index
 5darray survey_fleet_bio_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls)
 4darray survey_region_bio_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray survey_population_bio_overlap(1,nps,1,nyr,1,nps)
 matrix survey_natal_bio_overlap(1,nyr,1,nps)
 vector survey_total_bio_overlap(1,nyr)
 5darray survey_fleet_age(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 5darray survey_at_age_fleet_prop(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 5darray survey_fleet_age_bio(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 4darray survey_fleet_bio(1,nps,1,nr,1,nyr,1,nfls)
 3darray survey_region_bio(1,nps,1,nyr,1,nr)
 matrix survey_population_bio(1,nyr,1,nps)
 vector survey_total_bio(1,nyr) 
 //caa
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
 matrix SSB_overlap_natal(1,nps,1,nr)

 5darray yield_region_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray yield_population_temp_overlap(1,nps,1,nps,1,nyr,1,nag)
 3darray yield_natal_temp_overlap(1,nps,1,nyr,1,nps)
 5darray catch_at_age_population_temp_overlap(1,nps,1,nps,1,nyr,1,nag,1,nr)
 4darray catch_at_age_natal_temp_overlap(1,nps,1,nyr,1,nag,1,nps)
 5darray yield_fleet_temp(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 4darray yield_region_temp(1,nps,1,nr,1,nyr,1,nag)
 3darray yield_population_temp(1,nps,1,nyr,1,nag)
 3darray catch_at_age_total_temp(1,nyr,1,nag,1,nps)
 4darray catch_at_age_population_temp(1,nps,1,nyr,1,nag,1,nr)
 matrix yield_total_temp(1,nyr,1,nps)
 4darray SSB_region_temp(1,nps,1,nr,1,nyr,1,nag)
 matrix SSB_total_temp(1,nyr,1,nps)
 3darray SSB_population_temp(1,nps,1,nyr,1,nr)
 3darray biomass_population_temp(1,nps,1,nyr,1,nr)
 matrix biomass_total_temp(1,nyr,1,nps)
 4darray biomass_population_temp_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray biomass_natal_temp_overlap(1,nps,1,nyr,1,nps)
 4darray abundance_population_temp(1,nps,1,nyr,1,nag,1,nr)
 3darray abundance_total_temp(1,nyr,1,nag,1,nps)
 3darray total_recap_temp(1,nps,1,nr,1,tag_age)
 matrix tags_avail_temp(1,nps,1,nr)
 3darray tag_prop_temp(1,nps,1,nyr_rel,1,nr)
 5darray tag_prop_temp2(1,nps,1,nr,1,nyr_rel,1,nag,1,nt2)

 // likelihood components
  number survey_age_like
  number fish_age_like
  number rec_like
  number tag_like
  number tag_like_temp
  number catch_like 
  number survey_like
  number Tpen_like
 // init_bounded_number  theta(1,100,ph_theta);   // for negbinomial -lnL

 init_number dummy(ph_dummy)
  objective_function_value f;
  
  !! cout << "parameters set" << endl;

 //INITIALIZATION_SECTION  //set initial values
 //  steep .814;
 //  ln_q 0;
 //  ln_R_ave 7;
 //  log_sel_beta1 0;
 //  log_sel_beta2 2;
 //  log_sel_beta1surv 0;
 //  log_sel_beta2surv 2;
 //  ln_F -.7
 //  ln_rec_devs_RN 0;

PROCEDURE_SECTION
 
   get_movement();
   get_selectivity();
   get_F_age();
   get_vitals();
   get_SPR();
   get_abundance();
   get_survey_CAA_prop();
   get_CAA_prop();
   get_tag_recaptures();
   evaluate_the_objective_function();
///////BUILD MOVEMENT MATRIX////////
FUNCTION get_movement
 if(move_switch==0)
  {
     for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
     for(int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              if(j==k & r==n)
              {
               T(j,r,y,a,k,n)=1;
              }
              if(j!=k || r!=n)
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

 if(move_switch==1 || (phase_T_pop<0 && phase_T_reg<0 && phase_T_CNST_pop<0 && phase_T_CNST_reg<0 && move_switch!=0)) // if T fixed set it to input T
  {
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
     for(int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T(j,r,y,a,k,n)=T_year_TRUE(j,r,y,k,n);            
       }
      } 
     }
    }
   }
  }
  }
  
  if(move_switch>1)
  {
 if(npops>1 && sum(nregions)>npops && phase_T_pop>0) //not coded to do multiple regions AND multiple populations 
  {
      cout << "model not setup to estimate T for this pop structure" << endl;
  }

 if(npops==1 && sum(nregions)==1) //if panmictic then movement is 100%
  {

  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
     for(int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T(j,r,y,a,k,n)=1;            
       }
      } 
     }
    }
   }
  }
  }


  
 if(move_switch==2 && npops>1 && sum(nregions)==npops && phase_T_pop>0) //for mult pops but not mult regions
  {
   for(int y=1;y<=nyrs;y++)
    {
    G_pop=0;
    G_temp_pop=0;
     for (int j=1;j<=npops;j++)
      {
      for (int i=1;i<=npops;i++) 
       {
            if(j==i)
            {
            G_pop(j,i)=1;
            }
            if(i>j)
            {
            G_pop(j,i)=mfexp(ln_T_est_pop(y,j,i-1));
            }
            if(j!=i && i<j)
            {
            G_pop(j,i)=mfexp(ln_T_est_pop(y,j,i));
            }
        }
       }    
        G_temp_pop=rowsum(G_pop);     
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T(j,r,y,a,k,n)=G_pop(j,k)/G_temp_pop(j);
            
       }
      } 
     }
    }
   }
  }
  }

 if(move_switch==3 && npops==1 && sum(nregions)>npops && phase_T_reg>0) //for mult reg but not mult pops
  {
   for(int y=1;y<=nyrs;y++)
    {
    G_reg=0;
    G_temp_reg=0;
     for (int j=1;j<=nregions(1);j++)
      {
      for (int i=1;i<=nregions(1);i++) 
       {
            if(j==i)
            {
            G_reg(j,i)=1;
            }
            if(i>j)
            {
            G_reg(j,i)=mfexp(ln_T_est_reg(y,j,i-1));
            }
            if(j!=i && i<j)
            {
            G_reg(j,i)=mfexp(ln_T_est_reg(y,j,i));
            }
        }
       }    
        G_temp_reg=rowsum(G_reg);     
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T(j,r,y,a,k,n)=G_reg(r,n)/G_temp_reg(r);
       }
      } 
     }
    }
   }
  }
  }


 if(move_switch==4 && npops>1 && sum(nregions)==npops && phase_T_CNST_pop>0) //for mult pops but not mult regions
  {
    G_pop=0;
    G_temp_pop=0;
     for (int j=1;j<=npops;j++)
      {
      for (int i=1;i<=npops;i++) 
       {
            if(j==i)
            {
            G_pop(j,i)=1;
            }
            if(i>j)
            {
            G_pop(j,i)=mfexp(ln_T_CNST_pop(j,i-1));
            }
            if(j!=i && i<j)
            {
            G_pop(j,i)=mfexp(ln_T_CNST_pop(j,i));
            }
        }
       }    
        G_temp_pop=rowsum(G_pop);     
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for(int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T(j,r,y,a,k,n)=G_pop(j,k)/G_temp_pop(j);
            
       }
      } 
     }
    }
   }
  }
  }

 if(move_switch==5 && npops==1 && sum(nregions)>npops && phase_T_CNST_reg>0) //for mult reg but not mult pops
  {
    G_reg=0;
    G_temp_reg=0;
     for (int j=1;j<=nregions(1);j++)
      {
      for (int i=1;i<=nregions(1);i++) 
       {
            if(j==i)
            {
            G_reg(j,i)=1;
            }
            if(i>j)
            {
            G_reg(j,i)=mfexp(ln_T_CNST_reg(j,i-1));
            }
            if(j!=i && i<j)
            {
            G_reg(j,i)=mfexp(ln_T_CNST_reg(j,i));
            }
        }
       }    
        G_temp_reg=rowsum(G_reg);     
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for(int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T(j,r,y,a,k,n)=G_reg(r,n)/G_temp_reg(r);
       }
      } 
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
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T_terminal(j,r,a,k,n)=T(j,r,nyrs,a,k,n);            
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
///////SELECTIVITY CALCULATIONS///////
FUNCTION get_selectivity
 nreg_temp=rowsum(nregions_temp);

 if (ph_sel_log<0 && ph_sel_dubl<0)
 {
 sel_beta1=sel_beta1_TRUE;   //selectivity slope parameter 1 for logistic selectivity/double logistic
 sel_beta2=sel_beta2_TRUE;   //selectivity inflection parameter 1 for logistic selectivity/double logistic
 sel_beta3=sel_beta3_TRUE;  //selectivity slope parameter 2 for double selectivity
 sel_beta4=sel_beta4_TRUE;
 }
 else
 {
    for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
      {
        for (int z=1;z<=nfleets(j);z++)                    
      {
 // get betas on their arithmetic scale
 sel_beta1(j,r,z)=mfexp(log_sel_beta1(r+nreg_temp(j),z));
 sel_beta2(j,r,z)=mfexp(log_sel_beta2(r+nreg_temp(j),z));
 sel_beta3(j,r,z)=mfexp(log_sel_beta3(r+nreg_temp(j),z));
 sel_beta4(j,r,z)=mfexp(log_sel_beta4(r+nreg_temp(j),z));
  }}}
 } //close else for positive phase

 if (ph_sel_log_surv<0 && ph_sel_dubl_surv<0)
 {
 sel_beta1surv=sel_beta1_survey_TRUE;
 sel_beta2surv=sel_beta2_survey_TRUE;
 sel_beta3surv=sel_beta3_survey_TRUE;
 sel_beta4surv=sel_beta4_survey_TRUE;
 }
 else
 {
     for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
      {
        for (int z=1;z<=nfleets_survey(j);z++)                    
      {
 // get betas on their arithmetic scale
 sel_beta1surv(j,r,z)=mfexp(log_sel_beta1surv(j,z));
 sel_beta2surv(j,r,z)=mfexp(log_sel_beta2surv(j,z));
 sel_beta3surv(j,r,z)=mfexp(log_sel_beta3surv(j,z));
 sel_beta4surv(j,r,z)=mfexp(log_sel_beta4surv(j,z));
  }}}
 } //close else for survey select positive phase
 
 //fishery selectivity
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
               if(select_switch==2) //4 parameter double logistic selectivity
                {
                selectivity(j,r,y,a,z)=1/((1+mfexp(-sel_beta1(j,r,z)*(a-sel_beta2(j,r,z))))*(1+mfexp(-sel_beta3(j,r,z)*(a-sel_beta4(j,r,z)))));
                }
                if(select_switch==1) //two parameter logistic selectivity
                {
                selectivity(j,r,y,a,z)=1/(1+mfexp(-sel_beta1(j,r,z)*(a-sel_beta2(j,r,z)))); 
                }
        //       if(select_switch==0) //input selectivity at age constant by year
           //     {
               //  selectivity(j,r,y,a,z)=input_selectivity(j,r,a,z);
             //   }
              }
            }
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
                 survey_selectivity(j,r,y,a,z)=1/((1+mfexp(-sel_beta1surv(j,r,z)*(a-sel_beta2surv(j,r,z))))*(1+mfexp(-sel_beta3surv(j,r,z)*(a-sel_beta4surv(j,r,z)))));
                }
                if(select_switch_survey==1) //two parameter logistic selectivity
                {
                 survey_selectivity(j,r,y,a,z)=1/(1+mfexp(-sel_beta1surv(j,r,z)*(a-sel_beta2surv(j,r,z)))); //
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
             if(ph_F<0) //inefficient place for this because it'll reassign in all these loops but F_fleet needs to be in loops below...
             {
              F_year=F_year_TRUE;
             }
             else {
             if(F_switch==1) //estimate annual deviations
              {
               F_year(j,r,y,z)=mfexp(ln_F(y+nreg_temp(j)*nyrs+(r-1)*nyrs,z)); 
              }
             if(F_switch==2) //random walk or AR1  in F if we ever want to try it.
              {
               F_year(j,r,y,z)=mfexp(ln_F(y+nreg_temp(j)*nyrs+(r-1)*nyrs,z));  
               if(y>1)
               {
               F_year(j,r,y,z)=F_rho(j,z)*mfexp(((y-1)+nreg_temp(j)*nyrs+(r-1)*nyrs,z));            
               }
              }
             } //end else for negative F phase
             F_fleet(j,r,y,a,z)=F_year(j,r,y,z)*selectivity(j,r,y,a,z);    
           }
               F(j,r,y,a)=sum(F_fleet(j,r,y,a)); // moved down one loop I think was correct *dh*
           }
         }
        }
       }

FUNCTION get_vitals
 //set q to true if neg phase
   if(ph_q<0){
      q_survey=q_survey_TRUE;
      }
      
   else{
       for(int i=1;i<=npops;i++) 
    {
        for(int k=1;k<=nregions(i);k++) 
        {
          for(int j=1;j<=nfleets_survey(i);j++) 
          {
           q_survey(i,k,j) = mfexp(ln_q(i,j));}
    }
    }
    }
    
//set R_ave to true if neg phase
  if(ph_lmr<0){
  R_ave=R_ave_TRUE;
  }
  else{
  R_ave=mfexp(ln_R_ave+square(sigma_recruit)*0.5);
   }
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
              M(j,r,y,a)=input_M(j,a);
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

               if(SSB_type==1) //fecundity based SSB
                {
                 wt_mat_mult_reg(j,r,y,a)=prop_fem(j,r)*fecundity(j,r,a)*maturity(j,r,a);// for yearly SSB calcs
                }
               if(SSB_type==2) //weight based SSB
                {
                 wt_mat_mult_reg(j,r,y,a)=prop_fem(j,r)*weight_population(j,r,y,a)*maturity(j,r,a);
                }
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
               if(apportionment_type==(3)) //no yearly  NEED LOGIT TRANSFORM
                {
                 Rec_Prop(j,r,y)=ln_rec_prop(j,r);
                }
               if(apportionment_type==(4)) //Yearly  NEED LOGIT TRANSFORM
                {
                 Rec_Prop(j,r,y)=ln_rec_prop_year(j,r,y);
                }  
               }   
             }
           }         
         }
       }
     }

    for (int j=1;j<=npops;j++)
     { 
        for (int a=1;a<=nages;a++)
         {
          abund_devs(j,a)=mfexp(ln_abund_devs(j,a));
          if(ph_abund_devs<0)
          {
           abund_devs(j,a)=1;
          }
         }
      }
    for (int j=1;j<=npops;j++)
     { 
        for (int y=1;y<=nyrs-1;y++)
         {
               if(recruit_devs_switch==0)  //use population recruit relationship directly
                {
                 rec_devs(j,y)=1;
                }
               if(recruit_devs_switch==1)  // allow lognormal error around SR curve
                {
                 rec_devs(j,y)=mfexp(ln_rec_devs_RN(y+(j-1)*(nyrs-1))*sigma_recruit(j)-.5*square(sigma_recruit(j)));

                  if(ph_rec<0){
                   
                         rec_devs(j,y)=rec_devs_TRUE(j,y+1); //rec_devs vector in OM has length=nyrs, but only begins being used in year 2
                           }     
             if(recruit_randwalk_switch==1)
                 {
                  rec_devs_randwalk(j,y)=rec_devs(j,y);
                  if(y>(nages+1)) //start random walk in year3 based on year 2 devs as starting point (don't use equilibrium devs from year 1)
                   {
                    rec_devs(j,y)=rec_devs(j,y-1)*rec_devs_randwalk(j,y);  //is this correct?
                   }
                 }
                }
               }
              }
              
             
//SPR calcs are done with eitehr  average maturity/weight across all the regions within a population or assuming an input population fraction at equilibrium
// while the full SSB calcs use the region specific maturity/weight
FUNCTION get_SPR

  if(ph_steep<0){
       steep=steep_TRUE;
       }
       
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

FUNCTION get_abundance

       for (int y=1;y<=nyrs;y++)
        {

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
                    init_abund(p,j,r,a)=R_ave(p)*abund_devs(j,a)*pow(mfexp(-(M(p,j,r,a))),a);  //not sure what this is doing; JJD: I think this is estimating init_abundance as deviations from an exponential decline from Rave
                    if(ph_abund_devs<0)
                     {
                      init_abund(p,j,r,a)=init_abund_TRUE(p,j,r,a);
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
               rec_index_BM(j,r,y) = recruits_BM(j,r,y);
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

                 abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {              
                    if(move_switch!=6 || move_switch!=7 || a==1)  //if movement is not type=6 or a==1 (and movement type 6)
                     {
                       abundance_move_overlap_temp(k,n)=init_abund(p,k,n,a)*T(p,n,y,a,j,r); //with overlap always use natal population movement rates  (i.e., use p inpopsead of k)
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
                rec_index_AM(j,r,y)=recruits_AM(j,r,y);
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
                  survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*q_survey(j,r,z);
                  survey_fleet_overlap_age_bio(p,j,r,y,z,a)=survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  survey_fleet_bio_overlap(p,j,r,y,z)=sum(survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  survey_fleet_bio_overlap_temp(j,r,y,z,p)=survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*q_survey(j,r,z);
                  survey_fleet_age_bio(j,r,y,z,a)=survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_age_bio(j,r,y,z));
                 }
                if(natal_homing_switch==1)
                 {
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_bio_overlap_temp(j,r,y,z));
                 }
             
                  survey_region_bio_overlap(p,j,y,r)=sum(survey_fleet_bio_overlap(p,j,r,y));               
                  survey_population_bio_overlap(p,y,j)=sum(survey_region_bio_overlap(p,j,y));               
                  survey_natal_bio_overlap(y,p)=sum(survey_population_bio_overlap(p,y));               
                  survey_total_bio_overlap(y)=sum(survey_natal_bio_overlap(y));

                  survey_region_bio(j,y,r)=sum(survey_fleet_bio(j,r,y));
                  survey_population_bio(y,j)=sum(survey_region_bio(j,y));
                  survey_total_bio(y)=sum(survey_population_bio(y));
                  
                }  //tsurvey==0
               } //end survey_fleets
            }
           }
          }
         } //end age loop

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
                  catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F(j,r,y,a)+M(j,r,y,a))*(1-tspawn(p))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a))); // took out regional M *dh*
                  catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))*(1-tspawn(p))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a))); // took out regional M *dh*
                 }
                if(a>1)
                 {
                  catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a))); //// took out regional M *dh*
                  catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a))); // took out regional M *dh*            
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
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))*(1-tspawn(j))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a)); // took out regional M *dh*
                 }
                if(a>1)
                 {
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a)); // took out regional M *dh*
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

          } //end fleets loop
             for (int z=1;z<=nfleets_survey(j);z++)    /// survey index  1. Currently set up for more than 1 survey fleet
              {
               if(tsurvey(j,r)>0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                {
                  survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  survey_fleet_overlap_age_bio(p,j,r,y,z,a)=survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  survey_fleet_bio_overlap(p,j,r,y,z)=sum(survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  survey_fleet_bio_overlap_temp(j,r,y,z,p)=survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  survey_fleet_age_bio(j,r,y,z,a)=survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_age_bio(j,r,y,z));
                 }
                if(natal_homing_switch==1)
                 {
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_bio_overlap_temp(j,r,y,z));
                 }
             
                  survey_region_bio_overlap(p,j,y,r)=sum(survey_fleet_bio_overlap(p,j,r,y));               
                  survey_population_bio_overlap(p,y,j)=sum(survey_region_bio_overlap(p,j,y));               
                  survey_natal_bio_overlap(y,p)=sum(survey_population_bio_overlap(p,y));               
                  survey_total_bio_overlap(y)=sum(survey_natal_bio_overlap(y));

                  survey_region_bio(j,y,r)=sum(survey_fleet_bio(j,r,y));
                  survey_population_bio(y,j)=sum(survey_region_bio(j,y));
                  survey_total_bio(y)=sum(survey_population_bio(y));
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
                     SR(j,y-1)=R_ave(j);
                     total_recruits(j,y-1)=R_ave(j)*rec_devs(j,y-1);
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     SR(j,y-1)=R_ave(j);
                     total_recruits(j,y-1)=R_ave(j)*rec_devs(j,y-1);
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y-1)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }

                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a population
                    {
                    SR(j,y-1)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)));
                    total_recruits(j,y-1)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y-1);
                    recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a population
                    {
                     SR(j,y-1)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)));
                     total_recruits(j,y-1)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y-1);
                     recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y-1)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1)));
                    }
                  }
                  
                if(Rec_type==3) //environmental recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     SR(j,y-1)=env_rec(y);
                     total_recruits(j,y-1)=env_rec(y)*rec_devs(j,y-1);
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a population
                    {
                     SR(j,y-1)=env_rec(y);
                     total_recruits(j,y-1)=env_rec(y)*rec_devs(j,y-1);
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y-1)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
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
                    SR(j,y-1)=R_ave(j);
                    total_recruits(j,y-1)=R_ave(j)*rec_devs(j,y-1);
                    recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y-1)*Rec_Prop(j,r,y);
                    }
                    if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    SR(j,y-1)=R_ave(j);
                    total_recruits(j,y-1)=R_ave(j)*rec_devs(j,y-1);
                    recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y-1)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     SR(j,y-1)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)));
                     total_recruits(j,y-1)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y-1);
                     recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     SR(j,y-1)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)));
                     total_recruits(j,y-1)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y-1);
                     recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y-1)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                    }
                   }

                if(Rec_type==3) //average recruitment
                  {
                  if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     SR(j,y-1)=env_rec(y);
                     total_recruits(j,y-1)=env_rec(y)* rec_devs(j,y-1);
                     recruits_BM(j,r,y)=env_rec(y)* rec_devs(j,y-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                   {
                    SR(j,y-1)=env_rec(y);
                    total_recruits(j,y-1)=env_rec(y)* rec_devs(j,y-1);
                    recruits_BM(j,r,y)=env_rec(y)* R_ave(j)*rec_devs(j,y-1)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                   }
                 }
                 }
               rec_index_BM(j,r,y)=recruits_BM(j,r,y);
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_move_overlap_temp=0;
               
                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {

                 if(natal_homing_switch>0)
                 {
                 if(p==k)
                 {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y-1)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }

                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_population_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_population_overlap(p,k,y-1)))*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_population_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_population_overlap(p,k,y-1)))*rec_devs(k,y-1)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r);
                    }
                  }

                if(Rec_type==3) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y-1)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
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
                    abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                    }
                 if(Rec_type==2) //BH recruitment
                  {
                 if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    abundance_move_overlap_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }

                if(Rec_type==3) //environmental recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
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
                     abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                     }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    abundance_move_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                if(Rec_type==3) //env recruitment
                  {
                  if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
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
               rec_index_AM(j,r,y)=recruits_AM(j,r,y);
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
                  abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p))); // took out regional M
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

                abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7  ||a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,r,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(p)))*T(p,n,y,a,j,r); //with overlap always use natal population movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal population movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,r,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(p)))*return_probability(p); //with overlap always use natal population movement rates
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

                abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7 || a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,r,y-1,a-1)+F(k,n,y-1,a-1)))*T(p,n,y,a,j,r); //with overlap always use natal population movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal population movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,r,y-1,a-1)+F(k,n,y-1,a-1)))*return_probability(p); //with overlap always use natal population movement rates
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

                abundance_move_overlap_temp=0;

                for (int k=1;k<=npops;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7  || a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,r,y-1,a-1)+F(k,n,y-1,a-1)))*T(p,n,y,a,j,r)+abundance_at_age_AM_overlap_region(p,k,y-1,a,n)*mfexp(-(M(k,r,y-1,a)+F(k,n,y-1,a)))*T(p,n,y,a,j,r); //with overlap always use natal population movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal population movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,r,y-1,a-1)+F(k,n,y-1,a-1)))*return_probability(p)+abundance_at_age_AM_overlap_region(p,k,y-1,a,n)*mfexp(-(M(k,r,y-1,a)+F(k,n,y-1,a)))*return_probability(p); //with overlap always use natal population movement rates
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
                  survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*q_survey(j,r,z);
                  survey_fleet_overlap_age_bio(p,j,r,y,z,a)=survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  survey_fleet_bio_overlap(p,j,r,y,z)=sum(survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  survey_fleet_bio_overlap_temp(j,r,y,z,p)=survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*q_survey(j,r,z);
                  survey_fleet_age_bio(j,r,y,z,a)=survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_age_bio(j,r,y,z));
                 }
                if(natal_homing_switch==1)
                 {
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_bio_overlap_temp(j,r,y,z));
                 }
             
                  survey_region_bio_overlap(p,j,y,r)=sum(survey_fleet_bio_overlap(p,j,r,y));               
                  survey_population_bio_overlap(p,y,j)=sum(survey_region_bio_overlap(p,j,y));               
                  survey_natal_bio_overlap(y,p)=sum(survey_population_bio_overlap(p,y));               
                  survey_total_bio_overlap(y)=sum(survey_natal_bio_overlap(y));
                  
                  survey_region_bio(j,y,r)=sum(survey_fleet_bio(j,r,y));
                  survey_population_bio(y,j)=sum(survey_region_bio(j,y));
                  survey_total_bio(y)=sum(survey_population_bio(y));

                } //tsurvey==0
               } //end survey_fleets 
      }
     }
    }
   } //end age loop

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

          } //end fleets loop
             for (int z=1;z<=nfleets_survey(j);z++)    /// survey index  1. Currently set up for more than 1 survey fleet
              {
               if(tsurvey(j,r)>0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                {
                  survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  survey_fleet_overlap_age_bio(p,j,r,y,z,a)=survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                  survey_fleet_bio_overlap(p,j,r,y,z)=sum(survey_fleet_overlap_age_bio(p,j,r,y,z));  
                  survey_fleet_bio_overlap_temp(j,r,y,z,p)=survey_fleet_bio_overlap(p,j,r,y,z);

                if(natal_homing_switch==0)
                 {
                  survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                  survey_fleet_age_bio(j,r,y,z,a)=survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_age_bio(j,r,y,z));
                 }
                if(natal_homing_switch==1)
                 {
                  survey_fleet_bio(j,r,y,z)=sum(survey_fleet_bio_overlap_temp(j,r,y,z));
                 }
             
                  survey_region_bio_overlap(p,j,y,r)=sum(survey_fleet_bio_overlap(p,j,r,y));               
                  survey_population_bio_overlap(p,y,j)=sum(survey_region_bio_overlap(p,j,y));               
                  survey_natal_bio_overlap(y,p)=sum(survey_population_bio_overlap(p,y));               
                  survey_total_bio_overlap(y)=sum(survey_natal_bio_overlap(y));

                  survey_region_bio(j,y,r)=sum(survey_fleet_bio(j,r,y));
                  survey_population_bio(y,j)=sum(survey_region_bio(j,y));
                  survey_total_bio(y)=sum(survey_population_bio(y));

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
FUNCTION get_survey_CAA_prop      
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
              survey_at_age_region_fleet_overlap_prop(p,j,r,z,y,a)=survey_fleet_overlap_age(p,j,r,y,z,a)/sum(survey_fleet_overlap_age(p,j,r,y,z));              
              survey_at_age_fleet_prop(j,r,y,z,a)=survey_fleet_age(j,r,y,z,a)/sum(survey_fleet_age(j,r,y,z));
             }
            }
           }
          }
         }
        }


FUNCTION get_CAA_prop

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
                 catch_at_age_fleet_prop_temp(j,r,y,z,a)=catch_at_age_fleet(j,r,y,a,z);
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
            for (int a=1;a<=nages;a++)
             {
                 catch_at_age_region_fleet_overlap_prop(p,j,r,z,y,a)=catch_at_age_region_fleet_overlap(p,j,r,z,y,a)/sum(catch_at_age_region_fleet_overlap(p,j,r,z,y));              
                 catch_at_age_region_overlap_prop(p,j,r,y,a)=catch_at_age_region_overlap(p,j,r,y,a)/sum(catch_at_age_region_overlap(p,j,r,y));
                 catch_at_age_population_overlap_prop(p,j,y,a)=catch_at_age_population_overlap(p,j,y,a)/sum(catch_at_age_population_overlap(p,j,y));
                 catch_at_age_natal_overlap_prop(p,y,a)=catch_at_age_natal_overlap(p,y,a)/sum(catch_at_age_natal_overlap(p,y));

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
 if(do_tag==1)
  {
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
                 tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T(i,n,xx,a,j,r); 
                 recaps(i,n,x,a,y,j,r)=report_rate(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,xx,a)*(1.-mfexp(-(F(j,r,xx,a)+M(j,r,xx,a))))/(F(j,r,xx,a)+(M(j,r,xx,a)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                                
               }
              if(y>1) // must account for the maximum age so use min function to ensure that not exceeding the max age
              {
               tags_avail_temp=0;
               for(int p=1;p<=npops;p++)
               {
                for (int s=1;s<=nregions(p);s++)
                {
                 if(natal_homing_switch==0) //if no natal homing
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                 }                        
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current population only origin population and destination population
                //########################################################################################################              
                 if(natal_homing_switch==1) //if natal homing 
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(i,n,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)             
                 }               
                }
               }
                 tags_avail(i,n,x,a,y,j,r)=sum(tags_avail_temp); //sum across all pops/regs of tags that moved into pop j reg r
                 recaps(i,n,x,a,y,j,r)=report_rate(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),min((a+y),nages))*(1.-mfexp(-(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))))))/(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
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

FUNCTION evaluate_the_objective_function
  // f=dummy; //in case all the estimated parameters are turned off
   f=0.0;
  Tpen_like.initialize();
 if(move_pen_switch==1)
 {
 if(active(ln_T_est_reg))
  {
  for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
   {
     for (int j=1;j<=nregions(1);j++)
      {
      for (int i=1;i<=nregions(1)-1;i++) 
       {
          if(ln_T_est_reg(y,j,i)<Tpen)
           {
           Tpen_like+=(Tpen-ln_T_est_reg(y,j,i))*(Tpen-ln_T_est_reg(y,j,i));
          }
          if(ln_T_est_reg(y,j,i)>Tpen2)
           {
           Tpen_like+=(Tpen2-ln_T_est_reg(y,j,i))*(Tpen2-ln_T_est_reg(y,j,i));
          }
         }
        }
       }
      }

 if(active(ln_T_est_pop))
  {
  for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
   {
     for (int j=1;j<=npops;j++)
      {
      for (int i=1;i<=npops-1;i++) 
       {
          if(ln_T_est_pop(y,j,i)<Tpen)
           {
           Tpen_like+=(Tpen-ln_T_est_pop(y,j,i))*(Tpen-ln_T_est_pop(y,j,i));
          }
          if(ln_T_est_pop(y,j,i)>Tpen2)
           {
           Tpen_like+=(Tpen2-ln_T_est_pop(y,j,i))*(Tpen2-ln_T_est_pop(y,j,i));
          }
         }
        }
       }
      }

 if(active(ln_T_CNST_reg))
  {
     for (int j=1;j<=nregions(1);j++)
      {
      for (int i=1;i<=nregions(1)-1;i++) 
       {
          if(ln_T_CNST_reg(j,i)<Tpen)
           {
           Tpen_like+=(Tpen-ln_T_CNST_reg(j,i))*(Tpen-ln_T_CNST_reg(j,i));
          }
          if(ln_T_CNST_reg(j,i)>Tpen2)
           {
           Tpen_like+=(Tpen2-ln_T_CNST_reg(j,i))*(Tpen2-ln_T_CNST_reg(j,i));
          }
         }
        }
      }

 if(active(ln_T_CNST_pop))
  {
     for (int j=1;j<=npops;j++)
      {
      for (int i=1;i<=npops-1;i++) 
       {
          if(ln_T_CNST_pop(j,i)<Tpen)
           {
           Tpen_like+=(Tpen-ln_T_CNST_pop(j,i))*(Tpen-ln_T_CNST_pop(j,i));
          }
          if(ln_T_est_pop(y,j,i)>Tpen2)
           {
           Tpen_like+=(Tpen2-ln_T_CNST_pop(j,i))*(Tpen2-ln_T_CNST_pop(j,i));
          }
         }
        }
       }
   f+=Tpen_like;   
 }
 catch_like.initialize();
 fish_age_like.initialize();
 survey_age_like.initialize();
 survey_like.initialize();
 tag_like.initialize();
 tag_like_temp.initialize();
 rec_like.initialize();
 // Calculate multinomial likelihoods for compositions (Fournier style)
 // and survey biomass lognormal likelihood
      for (int j=1;j<=npops;j++)
     {
       for (int r=1;r<=nregions(j);r++)
       {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
             for (int z=1;z<=nfleets_survey(j);z++)
              {
            survey_like +=   wt_srv*square((log(OBS_survey_fleet_bio(j,r,y,z)+0.0001)-log(survey_fleet_bio(j,r,y,z)+0.0001) )/ (2.*square(OBS_survey_fleet_bio_se(j,r,y,z)))); //OBS_survey_fleet_bio(j,r,y,z))));
              survey_age_like -= OBS_survey_prop_N(j,r,y,z) * ((OBS_survey_prop(j,r,y,z)+0.001)*log(survey_at_age_fleet_prop(j,r,y,z)+0.001));
             }
            }
           }
          }
        
     
     // catch likelihood and multinomial fishery ages
           for (int j=1;j<=npops;j++)
     {
       for (int r=1;r<=nregions(j);r++)
       {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
             for (int z=1;z<=nfleets(j);z++)
              {
           catch_like+= square((log(OBS_yield_fleet(j,r,y,z)+0.0001)-log(yield_fleet(j,r,y,z)+0.0001) )/ (2.*square(OBS_yield_fleet_se(j,r,y,z)))); //OBS_yield_fleet(j,r,y,z))));
           fish_age_like -= OBS_catch_at_age_fleet_prop_N(j,r,y,z)*((OBS_catch_at_age_fleet_prop(j,r,y,z)+0.001)*log(catch_at_age_fleet_prop(j,r,y,z)+0.001));
             }
             }
          }
         }

 if(active(ln_rec_devs_RN))
  {
   for(int j=1;j<=npops;j++)
    {  //is this correct?  we aren't penalizing against devs from Rave or SR function?  also we do want to do calcs across all years right (since not including year index in rec_devs here)?
     for(int y=1;y<=nyrs-1;y++)
      {
       rec_like +=log(sigma_recruit(j))+(square(log(total_recruits(j,y)+.00001)-log(SR(j,y)+.00001)))/(2*square(sigma_recruit(j)));
      }
    }
  }

 //if(active(ln_rec_devs_RN))
// {
// for(int j=1;j<=npops;j++) {  //is this correct?  we aren't penalizing against devs from Rave or SR function?  also we do want to do calcs across all years right (since not including year index in rec_devs here)?
// rec_like += (norm2(ln_rec_devs_RN(j)+sigma_recruit(j)*sigma_recruit(j)/2.)/(2.*square(sigma_recruit(j))) + (size_count(ln_rec_devs_RN(j)))*log(sigma_recruit(j))); 
// }}

// if(do_tag_mult==0)
//  {
//   for (int i=1;i<=npops;i++)
//  {
//   for (int n=1;n<=nregions(i);n++)
//    {
//    for(int x=1; x<=nyrs_release; x++)
//     {
//           xx=yrs_releases(x);
///      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
 //       {
 ///        for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
  //        {  
   //       for(int j=1;j<=npops;j++) //recap stock
    //        {
    //         for (int r=1;r<=nregions(j);r++)
     //        {
      
    //         tag_like -= log_negbinomial_density(OBS_recaps(i,n,x,a,y,j,r),recaps(i,n,x,a,y,j,r)+0.00001,theta);  //negative binomial tag likelihood
     //       }}}}}}} 
  // I have done Poisson and negative-binomial, if you want multinomial than that's up to you
  // Right now this has a 7d array for recaps and a 6d for OBS
 // }

 if(do_tag_mult==1)
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
        if(ntags(i,n,x,a)==0)
         {
          OBS_tag_prop_N(i,n,x,a)==0; //make tag likelihood==0 if there are no releases in a given cohort (i.e., mainly if no releases at young ages due to selectivity==0)
         }
         for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
           {
            if(max_life_tags<=(nyrs-xx+1)) //complete cohorts
             {
              tag_like -= OBS_tag_prop_N(i,n,x,a) * ((tag_prop_final(i,n,x,a,s)+0.001)*log(OBS_tag_prop_final(i,n,x,a,s)+0.001));
             }
            if(max_life_tags>(nyrs-xx+1)) //need special calcs for incomplete cohorts (ie model ends before end max_life_tags reached)
             {
               if(s<((nyrs-xx+1)*sum(nregions)+1))
                {
                 tag_like_temp +=((tag_prop_final(i,n,x,a,s)+0.001)*log(OBS_tag_prop_final(i,n,x,a,s)+0.001));
                }
               if(s==((nyrs-xx+1)*sum(nregions)+1))
                {
                 tag_like_temp += ((tag_prop_final(i,n,x,a,(max_life_tags*sum(nregions)+1))+0.001)*log(OBS_tag_prop_final(i,n,x,a,(max_life_tags*sum(nregions)+1))+0.001));
                }
               tag_like -= OBS_tag_prop_N(i,n,x,a)*tag_like_temp;
             }


            }}}}}
  }


  /// Early penalty to keep F under wraps
 if (current_phase()<last_phase()) f+= 1000*(norm2(mfexp(ln_F)));

 //may be better way to penalize....
 if(abund_pen_switch==1)
 {
  if (active(ln_abund_devs)) f+= norm2(ln_abund_devs); 
 } 
// Sum objective function
   f           += survey_like*wt_srv;
   f           += catch_like*wt_catch;
   f           += fish_age_like*wt_fish_age;
   f           += survey_age_like*wt_srv_age;
   f           += rec_like*wt_rec;
   f           += tag_like*wt_tag;
  
REPORT_SECTION
  report<<"$nages"<<endl;
  report<<nages<<endl;
  report<<"$nyrs"<<endl;
  report<<nyrs<<endl;
  report<<"$npops"<<endl;
  report<<npops<<endl;
  report<<"$nregions"<<endl;
  report<<nregions<<endl;
  report<<"$nfleets"<<endl;
  report<<nfleets<<endl;
  report<<"$nfleets_survey"<<endl;
  report<<nfleets_survey<<endl;
  
  report<<"$sigma_recruit"<<endl;
  report<<sigma_recruit<<endl;

  report<<"$M"<<endl;
  report<<M<<endl;

 //EST values
  report<<"$T_terminal"<<endl;
  report<<T_terminal<<endl;
  report<<"$T_year"<<endl;
  report<<T_year<<endl;
  report<<"$init_abund"<<endl;
  report<<init_abund<<endl;
  report<<"$alpha"<<endl;
  report<<alpha<<endl;
  report<<"$beta"<<endl;
  report<<beta<<endl;
  report<<"$input_T"<<endl;
  report<<input_T<<endl;
  report<<"$q_survey"<<endl;
  report<<q_survey<<endl;
  report<<"$sel_beta1"<<endl;
  report<<sel_beta1<<endl;
  report<<"$sel_beta2"<<endl;
  report<<sel_beta2<<endl;
  report<<"$sel_beta3"<<endl;
  report<<sel_beta3<<endl;
  report<<"$sel_beta4"<<endl;
  report<<sel_beta4<<endl;
  report<<"$sel_beta1_survey"<<endl;
  report<<sel_beta1surv<<endl;
  report<<"$sel_beta2_survey"<<endl;
  report<<sel_beta2surv<<endl;
  report<<"$sel_beta3_survey"<<endl;
  report<<sel_beta3surv<<endl;
  report<<"$sel_beta4_survey"<<endl;
  report<<sel_beta4surv<<endl;
  report<<"$steep"<<endl;
  report<<steep<<endl;
  report<<"$R_ave"<<endl;
  report<<R_ave<<endl;
  report<<"$SSB_zero"<<endl;
  report<<SSB_zero<<endl; 
  report<<"$rec_devs"<<endl;
  report<<rec_devs<<endl;
  report<<"$Rec_Prop"<<endl;
  report<<Rec_Prop<<endl;
  report<<"$recruits_BM"<<endl;
  report<<recruits_BM<<endl;
  report<<"$F"<<endl;
  report<<F<<endl;
  report<<"$F_year"<<endl;
  report<<F_year<<endl;
  
  report<<"$biomass_AM"<<endl;
  report<<biomass_AM<<endl;
  report<<"$biomass_population"<<endl;
  report<<biomass_population<<endl;
  report<<"$harvest_rate_region_bio"<<endl;
  report<<harvest_rate_region_bio<<endl;
  report<<"$depletion_region"<<endl;
  report<<depletion_region<<endl;
  report<<"$SSB_region"<<endl;
  report<<SSB_region<<endl;
  report<<"$Bratio_population"<<endl;
  report<<Bratio_population<<endl;
  
 /// TRUE VALUES
  report<<"$input_T"<<endl;
  report<<input_T<<endl;
  report<<"$T_year_TRUE"<<endl;
  report<<T_year_TRUE<<endl;
  report<<"$q_survey_TRUE"<<endl;
  report<<q_survey_TRUE<<endl;
  report<<"$sel_beta1_TRUE"<<endl;
  report<<sel_beta1_TRUE<<endl;
  report<<"$sel_beta2_TRUE"<<endl;
  report<<sel_beta2_TRUE<<endl;
  report<<"$sel_beta3_TRUE"<<endl;
  report<<sel_beta3_TRUE<<endl;
  report<<"$sel_beta4_TRUE"<<endl;
  report<<sel_beta4_TRUE<<endl;
  report<<"$sel_beta1_survey_TRUE"<<endl;
  report<<sel_beta1_survey_TRUE<<endl;
  report<<"$sel_beta2_survey_TRUE"<<endl;
  report<<sel_beta2_survey_TRUE<<endl;
  report<<"$sel_beta3_survey_TRUE"<<endl;
  report<<sel_beta3_survey_TRUE<<endl;
  report<<"$sel_beta4_survey_TRUE"<<endl;
  report<<sel_beta4_survey_TRUE<<endl; 
  report<<"$steep_TRUE"<<endl;
  report<<steep_TRUE<<endl;
  report<<"$R_ave_TRUE"<<endl;
  report<<R_ave_TRUE<<endl;
  report<<"$SSB_zero_TRUE"<<endl;
  report<<SSB_zero_TRUE<<endl; 
  report<<"$rec_devs_TRUE"<<endl;
  report<<rec_devs_TRUE<<endl;
  report<<"$Rec_Prop_TRUE"<<endl;
  report<<Rec_Prop_TRUE<<endl;
  report<<"$recruits_BM_TRUE"<<endl;
  report<<recruits_BM_TRUE<<endl;
  report<<"$F_TRUE"<<endl;
  report<<F_TRUE<<endl;
  report<<"$F_year_TRUE"<<endl;
  report<<F_year_TRUE<<endl;
  report<<"$Init_Abund_TRUE"<<endl;
  report<<init_abund_TRUE<<endl;
  report<<"$Init_Abund"<<endl;
  report<<init_abund<<endl;
  report<<"$Init_Abund_Devs"<<endl;
  report<<abund_devs<<endl;
  report<<"$biomass_AM_TRUE"<<endl;
  report<<biomass_AM_TRUE<<endl;
  report<<"$biomass_population_TRUE"<<endl;
  report<<biomass_population_TRUE<<endl;
  report<<"$harvest_rate_region_bio_TRUE"<<endl;
  report<<harvest_rate_region_bio_TRUE<<endl;
  report<<"$depletion_region_TRUE"<<endl;
  report<<depletion_region_TRUE<<endl;
  report<<"$SSB_region_TRUE"<<endl;
  report<<SSB_region_TRUE<<endl;
  report<<"$Bratio_population_TRUE"<<endl;
  report<<Bratio_population_TRUE<<endl;

 //OBS and Pred Values
  report<<"$OBS_survey_fleet_bio"<<endl;
  report<<OBS_survey_fleet_bio<<endl;
  report<<"$survey_fleet_bio"<<endl;
  report<<survey_fleet_bio<<endl;
  report<<"$OBS_survey_fleet_bio_se"<<endl;
  report<<OBS_survey_fleet_bio_se<<endl;
  report<<"$OBS_survey_prop"<<endl;
  report<<OBS_survey_prop<<endl;
  report<<"$survey_prop"<<endl;
  report<<survey_at_age_fleet_prop<<endl;

  report<<"$OBS_yield_fleet"<<endl;
  report<<OBS_yield_fleet<<endl;
  report<<"$yield_fleet"<<endl;
  report<<yield_fleet<<endl;
  report<<"$OBS_yield_fleet_se"<<endl;
  report<<OBS_yield_fleet_se<<endl;
  report<<"$OBS_catch_prop"<<endl;
  report<<OBS_catch_at_age_fleet_prop<<endl;
  report<<"$catch_at_age_fleet_prop"<<endl;
  report<<catch_at_age_fleet_prop<<endl;

/// TAG INFORMATION
  report<<"$nyrs_release"<<endl;
  report<<nyrs_release<<endl;
  report<<"$years_of_tag_releases "<<endl;
  report<<yrs_releases<<endl;
  report<<"$max_life_tags"<<endl;
  report<<max_life_tags<<endl;
  report<<"$report_rate"<<endl;
  report<<report_rate<<endl;
  report<<"$ntags_total"<<endl;
  report<<ntags_total<<endl;
  report<<"$ntags"<<endl;
  report<<ntags<<endl;
  report<<"$OBS_tag_prop_final"<<endl;
  report<<OBS_tag_prop_final<<endl;
  report<<"$tag_prop_final"<<endl;
  report<<tag_prop_final<<endl;

  report<<"$likelihood components"<<endl;
  report<<"$tag_like"<<endl;
  report<<tag_like<<endl;
  report<<"$fish_age_like"<<endl;
  report<<fish_age_like<<endl;
  report<<"$survey_age_like"<<endl;
  report<<survey_age_like<<endl;
  report<<"$survey_like"<<endl;
  report<<survey_like<<endl;
  report<<"$catch_like"<<endl;
  report<<catch_like<<endl;
  report<<"$rec_like"<<endl;
  report<<rec_like<<endl;

  report<<"$selectivity_age"<<endl;
  report<<selectivity_age<<endl;
  report<<"$selectivity_age_TRUE"<<endl;
  report<<selectivity_age_TRUE<<endl;
  report<<"$survey_selectivity_age"<<endl;
  report<<survey_selectivity_age<<endl;
  report<<"$survey_selectivity_age_TRUE"<<endl;
  report<<survey_selectivity_age_TRUE<<endl;

  report<<"$total_recruit"<<endl;
  report<<total_recruits<<endl;
  report<<"$SR"<<endl;
  report<<SR<<endl;
  
  save_gradients(gradients);
RUNTIME_SECTION
  convergence_criteria .001,.0001, 1.0e-4, 1.0e-7
  maximum_function_evaluations 100000
  


