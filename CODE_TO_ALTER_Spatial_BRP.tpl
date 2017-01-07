/////////////////////////////////////////////////////////
// Spatial Operating Model built by Daniel Goethel (NMFS SEFSC)  
// edited by Katelyn Bosley
//////////////////////////////////////////////////////////

//need to add all_natal calcs for all values that might want to compare across an area, because regular non-natal
//homing calcs do not account for different weights across natal populations so if if any part of calc involves weight
//it needs to use the biomass_all_natal value not just biomass_AM


GLOBALS_SECTION
  #include "admodel.h"
  #define EOUT(var) cout <<#var<<" "<<var<<endl;
 // enum MovementType{
  //     NO_MOVEMENT = 0,
  //     }

TOP_OF_MAIN_SECTION
  arrmblsize=500000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);


//input data file contains values for switches in the model

DATA_SECTION
////////////////////////////////////////////////////////////////////////////////////
/////MODEL STRUCTURE INPUTS/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  init_int nages //number of ages
  init_int nyrs //number of years for simulation
  init_int nstocks //number of stocks...likely to swtich to populations

///////////// DESCRIBE ALL OF THESE RAGGED ARRAYS and population setups (regions per area fleets per region)
  !! int ns=nstocks;
  !! int ny=nyrs;
  !! int na=nages;
  //////////////////////////////////////////////////////
  
  init_ivector nregions(1,ns) //number of regions within a population - for metamictic regions = areas, stocks = 1

  !! ivector nreg=nregions;

  init_imatrix nfleets(1,ns,1,nreg) //number of fleets in each region by each stock (population)

  !! imatrix nf=nfleets;


////////////////////////////////////////////////////////////////////////////////////
//////////////SWITCHES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

///// Changes the type of harvest model 
  init_number model_type_switch
  //==0 use to simulate catch based on input F
  //==1 use input TAC to set F
  //==2 use input harvest rate


///// Changes the type of larval movement pattern (sets age class 1 movements)
  init_number larval_move_switch
  //==0 no movement
  //==1 input movement
  //==2 movement within stock only based on residency (symmetric)
  //==3 symmetric movement but only allow movement within a stock (ie regions within a stock) not across stocks
  //==4 symmetric movement across all stocks and regions
  //==5 allow movement across all regions and stocks, based on stock/region specific residency (symmetric off-diag)


////// ADD MOVEMENT SWITCHES DENSITY-DEPENDENT, possibly for metamictic with ontogenetic movement
  init_number move_switch
  //==0 no movement
  //==1 input movement
  //==2 movement within stock only based on input residency (symmetric off diagnol)
  //==3 symmetric movement but only allow movement within a stock (ie regions within a stock) not across stocks
  //==4 symmetric movement across all stocks and regions
  //==5 allow movement across all regions and stocks, based on stock/region specific residency (symmetric off diagnol)
  //==6 natal return, based on age of return and return probability (certain fraction of fish make return migration to natal stock eg ontogenetic migration)
  //==7 larvae stay in the stock that they move to (i.e., for overlap, do not return to natal stock if adult movement==0...otherwise with natal
  //    homing would return to natal stock because natal residency is 100% and use natal movement rates (not current stock movement rates like with metapopulation/random movement))
//////////////////////////////////////////////////////

////// SEPARATE INTO POP STRUCTURE AND SPAWN RETURN SWITCHES
  init_number overlap_switch
  //==0 no overlap (SSB is sum of SSB in stock regardless of natal origin; weight/mat/fecund/ are based on current stock not natal stock) - Metapopulation/metamictic
  //==1 do overlap (a fish only adds to SSB if it is in its natal stock at spawning time; weight/mat/fecund/ are based on natal stock) - Natal homing
  //==2 do overlap with natal migration (a fraction of fish return to natal stock to spawn (instantaneous migration to natal population and back at time of spawning) based spawn_return_prob; weight/mat/fecund/ are based on natal stock)
  //overlap>0 assumes genetic based life history and contribution to SSB (i.e., natal homing and no demographic mixing), overlap==0 assumes demographic mixing (e.g. metapopulations where life history is more location based)
//////////////////////////////////////////////////////

////// ADD DOUBLE SELECTIVITY ////////////////////////
  init_number select_switch
  //==0 input selectivity
  //==1 logistic selectivity based on input sel_beta1 and sel_beta2
  //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4
/////////////////////////////////////////////////////

  init_number F_switch
  //==1 input F
  //==2 input single overall F (FMSY)
  //==3 Estimate F that minimizes difference in input and estimated total harvest rate
  //==4 overall F (FMSY) is split evenly among stocks (each fleet uses stock F)
  //==5 overall F (FMSY) is is split evenly among all regions (each fleet uses region F)
  //==6 overall F (FMSY) is split evenly among fleets
  //==7 random walk in F (NOT YET IMPLEMENTED)
  
  init_number recruit_devs_switch
  //==0 use stock-recruit relationship directly
  //==1 allow lognormal error around SR curve (i.e., include randomness based on input sigma_recruit)
  
  init_number SSB_type
  //==1 fecundity based SSB
  //==2 weight based SSB
  
  init_number apportionment_type
  //==-1 no recruitment apportionment to regions within a stock (each region within a stock gets full amount of recruits from SR curve)
  //==0 apportionment to each region is based on relative SSB in region compared to stock SSB
  //==1 input apportionment
  //==2 recruits are apportioned equally to each region within a stock
  
  init_number Rec_type
  //==1 stock-recruit relationship assumes an average value based on R_ave
  //==2 Beverton-Holt stock-recruit functions based on stock-specific input steepness, R0 (R_ave), M, and weight
  //==3 input recruitment - use for 


///////////////////////////////////////////////////////////////////////////////
//////// ADDITIONAL PARAMETERS FROM DAT FILE //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  init_number return_age // used if move_swith ==6
  init_vector return_probability(1,ns) // used if overlap_swith==6
  init_vector spawn_return_prob(1,ns) // used if overlap_swith==2
  init_number phase_F //must be turned on (==1) if F_type==3
  init_number phase_dummy //must be turned on (==1) if F_type!=3
  init_vector tspawn(1,ns) //time of spawning
  init_vector steep(1,ns) //B-H steepness
  init_vector R_ave(1,ns) //Average Recruitment or R0 for B-H S-R curve
  init_vector amplitude(1,ns) //amplitude of periodic recruitment in % of R_ave 
  init_vector freq(1,ns) //frequency of recruitment in years (ie 10 for peak every 10 years)
  //check that this ragged array is set up properly
  init_5darray input_T(1,ns,1,nreg,1,na,1,ns,1,nreg)  //movement matrix 
  init_matrix input_residency_larval(1,ns,1,nreg)  //larval residency probability
  init_3darray input_residency(1,ns,1,nreg,1,na) //
  init_3darray sel_beta1(1,ns,1,nreg,1,nf)   //selectivity slope parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta2(1,ns,1,nreg,1,nf)   //selectivity inflection parameter 1 for logistic selectivity/double logistic
  init_3darray sel_beta3(1,ns,1,nreg,1,nf)  //selectivity slope parameter 2 for double selectivity
  init_3darray sel_beta4(1,ns,1,nreg,1,nf)  //selectivity inflection parameter 2 for double logistic selectivity
  init_4darray input_selectivity(1,ns,1,nreg,1,na,1,nf) //
  init_3darray input_F(1,ns,1,nreg,1,nf)
  init_number input_F_MSY
  init_matrix input_M(1,ns,1,na)
  init_vector sigma_recruit(1,ns)
  init_matrix input_weight(1,ns,1,na)
  init_matrix input_catch_weight(1,ns,1,na)
  init_matrix fecundity(1,ns,1,na)
  init_matrix maturity(1,ns,1,na)
  init_matrix input_Rec_prop(1,ns,1,nreg)
  init_3darray init_abund(1,ns,1,nreg,1,na)
  init_matrix rec_index_sigma(1,ns,1,nreg)
  init_vector caa_sigma(1,na)
 
  init_3darray input_TAC(1,ns,1,nreg,1,nf)
  init_3darray input_u(1,ns,1,nreg,1,nf)
  init_number max_Fnew //
  init_number Fnew_start
  init_number NR_iterations
  init_number NR_dev
  
  init_int debug
  init_number myseed

  vector years(1,nyrs)
  !!years.fill_seqadd(double(1),1.0);
  
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

 !! cout << "debug = " << debug << endl;
 !! cout << "If debug != 1541 then .dat file not setup correctly" << endl;
 !! cout << "input read" << endl;

 
 
 //!!cout << rec_index_sigma << endl;
 //!!cout << freq << endl;
 //!!exit(43);

PARAMETER_SECTION

  !! ivector nr=nregions;
  !! int nst=nstocks;
  !! int nyr=nyrs;
  !! int nag=nages;
  !! imatrix nfl=nfleets;

 init_matrix F_est(1,nst,1,nr,phase_F)


 6darray T(1,nst,1,nr,1,nyr,1,nag,1,nst,1,nr)
 5darray selectivity(1,nst,1,nr,1,nyr,1,nag,1,nfl)
 4darray F_year(1,nst,1,nr,1,nyr,1,nfl)
 5darray F_fleet(1,nst,1,nr,1,nyr,1,nag,1,nfl)
 4darray F(1,nst,1,nr,1,nyr,1,nag)
 4darray M(1,nst,1,nr,1,nyr,1,nag)
 matrix rec_devs(1,nst,1,nyr)
 3darray weight_stock(1,nst,1,nyr,1,nag)
 3darray weight_catch(1,nst,1,nyr,1,nag)
 3darray wt_mat_mult(1,nst,1,nyr,1,nag)
 3darray weight_stock_overlap(1,nst,1,nyr,1,nag)
 3darray weight_catch_overlap(1,nst,1,nyr,1,nag)
 3darray wt_mat_mult_overlap(1,nst,1,nyr,1,nag)
 3darray Rec_Prop(1,nst,1,nr,1,nyr)
 matrix SPR_N(1,nst,1,nag)
 matrix SPR_SSB(1,nst,1,nag)
 vector SPR(1,nst)
 vector SSB_zero(1,nst)
 vector alpha(1,nst)
 vector beta(1,nst)

 3darray recruits_BM(1,nst,1,nr,1,nyr)
 3darray recruits_AM(1,nst,1,nr,1,nyr)
 
 3darray rec_index_BM(1,nst,1,nr,1,nyr)

 vector env_rec(1,nyr)
 
 4darray abundance_at_age_BM(1,nst,1,nr,1,nyr,1,nag)
 4darray abundance_at_age_AM(1,nst,1,nr,1,nyr,1,nag)
 4darray abundance_in(1,nst,1,nr,1,nyr,1,nag)
 4darray abundance_res(1,nst,1,nr,1,nyr,1,nag)
 4darray abundance_leave(1,nst,1,nr,1,nyr,1,nag)
 4darray abundance_spawn(1,nst,1,nr,1,nyr,1,nag)
 
 4darray biomass_BM_age(1,nst,1,nr,1,nyr,1,nag)
 4darray biomass_AM_age(1,nst,1,nr,1,nyr,1,nag)
 3darray biomass_BM(1,nst,1,nr,1,nyr)
 3darray biomass_AM(1,nst,1,nr,1,nyr)
 4darray bio_in(1,nst,1,nr,1,nyr,1,nag)
 4darray bio_res(1,nst,1,nr,1,nyr,1,nag)
 4darray bio_leave(1,nst,1,nr,1,nyr,1,nag)

 5darray catch_at_age_fleet(1,nst,1,nr,1,nyr,1,nag,1,nfl)
 4darray yield_fleet(1,nst,1,nr,1,nyr,1,nfl)
 4darray catch_at_age_region(1,nst,1,nr,1,nyr,1,nag)
 3darray yield_region(1,nst,1,nr,1,nyr)
 3darray catch_at_age_stock(1,nst,1,nyr,1,nag)
 matrix yield_stock(1,nst,1,nyr)
 3darray SSB_region(1,nst,1,nr,1,nyr)
 matrix SSB_stock(1,nst,1,nyr)
 vector SSB_total(1,nyr)
 3darray abundance_stock(1,nst,1,nyr,1,nag)
 matrix abundance_total(1,nyr,1,nag)
 matrix biomass_stock(1,nst,1,nyr)
 vector biomass_total(1,nyr)
 matrix catch_at_age_total(1,nyr,1,nag)
 vector yield_total(1,nyr)
 4darray harvest_rate_region_num(1,nst,1,nr,1,nyr,1,nag)
 3darray harvest_rate_stock_num(1,nst,1,nyr,1,nag)
 matrix harvest_rate_total_num(1,nyr,1,nag)
 3darray harvest_rate_region_bio(1,nst,1,nr,1,nyr)
 matrix harvest_rate_stock_bio(1,nst,1,nyr)
 vector harvest_rate_total_bio(1,nyr)
 3darray depletion_region(1,nst,1,nr,1,nyr)
 matrix depletion_stock(1,nst,1,nyr)
 vector depletion_total(1,nyr)
 5darray abundance_at_age_BM_overlap_region(1,nst,1,nst,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_BM_overlap_stock(1,nst,1,nst,1,nyr,1,nag)
 5darray abundance_at_age_AM_overlap_region(1,nst,1,nst,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_AM_overlap_stock(1,nst,1,nst,1,nyr,1,nag)
 4darray abundance_AM_overlap_region_all_natal(1,nst,1,nr,1,nyr,1,nag)
 5darray abundance_spawn_overlap(1,nst,1,nst,1,nr,1,nyr,1,nag)
 4darray SSB_region_overlap(1,nst,1,nst,1,nr,1,nyr)
 3darray SSB_stock_overlap(1,nst,1,nst,1,nyr)
 matrix SSB_natal_overlap(1,nst,1,nyr)

 4darray true_caa_fleet_total(1,nst,1,nr,1,nyr,1,nfl)
 4darray true_caa_overlap_reg_total(1,nst,1,nst,1,nr,1,nyr)
 5darray true_caa_fleet_temp(1,nst,1,nr,1,nyr,1,nfl,1,nag)
 5darray true_caa_overlap_reg_temp(1,nst,1,nst,1,nr,1,nyr,1,nag)
 
 5darray CAA_prop_fleet(1,nst,1,nr,1,nyr,1,nfl,1,nag)
 5darray CAA_prop_overlap(1,nst,1,nr,1,nyr,1,nfl,1,nag)
 
 5darray obs_caa_fleet(1,nst,1,nr,1,nyr,1,nfl,1,nag)
 4darray obs_caa_fleet_total(1,nst,1,nr,1,nyr,1,nfl)
 5darray obs_prop_fleet(1,nst,1,nr,1,nyr,1,nfl,1,nag)

 5darray obs_caa_overlap_reg(1,nst,1,nst,1,nr,1,nyr,1,nag)
 4darray obs_caa_overlap_reg_total(1,nst,1,nst,1,nr,1,nyr)
 5darray obs_prop_overlap_reg(1,nst,1,nst,1,nr,1,nyr,1,nag)

 5darray catch_at_age_region_overlap(1,nst,1,nst,1,nr,1,nyr,1,nag)
 4darray yield_region_overlap(1,nst,1,nst,1,nr,1,nyr)
 4darray catch_at_age_stock_overlap(1,nst,1,nst,1,nyr,1,nag)
 3darray yield_stock_overlap(1,nst,1,nst,1,nyr)
 3darray abundance_natal_overlap(1,nst,1,nyr,1,nag)
 5darray biomass_BM_age_overlap(1,nst,1,nst,1,nr,1,nyr,1,nag)
 5darray biomass_AM_age_overlap(1,nst,1,nst,1,nr,1,nyr,1,nag)
 4darray biomass_BM_overlap_region(1,nst,1,nst,1,nr,1,nyr)
 4darray biomass_AM_overlap_region(1,nst,1,nst,1,nr,1,nyr)
 5darray biomass_AM_overlap_region_all_natal_temp(1,nst,1,nr,1,nyr,1,nag,1,nst)
 4darray biomass_AM_overlap_age_region_all_natal(1,nst,1,nr,1,nyr,1,nag)
 3darray biomass_AM_overlap_region_all_natal(1,nst,1,nr,1,nyr)
 3darray biomass_stock_overlap(1,nst,1,nst,1,nyr)
 matrix biomass_natal_overlap(1,nst,1,nyr)
 3darray catch_at_age_natal_overlap(1,nst,1,nyr,1,nag)
 matrix yield_natal_overlap(1,nst,1,nyr)
 4darray harvest_rate_region_bio_overlap(1,nst,1,nst,1,nr,1,nyr)
 3darray harvest_rate_stock_bio_overlap(1,nst,1,nst,1,nyr)
 matrix harvest_rate_natal_bio_overlap(1,nst,1,nyr)
 4darray depletion_region_overlap(1,nst,1,nst,1,nr,1,nyr)
 3darray depletion_stock_overlap(1,nst,1,nst,1,nyr)
 matrix depletion_natal_overlap(1,nst,1,nyr)
 3darray Bratio_stock_overlap(1,nst,1,nst,1,nyr)
 matrix Bratio_natal_overlap(1,nst,1,nyr)
 matrix Bratio_stock(1,nst,1,nyr)
 vector Bratio_total(1,nyr)
 
 matrix abundance_move_temp(1,nst,1,nr)
 matrix bio_move_temp(1,nst,1,nr)
 5darray yield_fleet_temp(1,nst,1,nr,1,nyr,1,nfl,1,nag)
 4darray yield_region_temp(1,nst,1,nr,1,nyr,1,nag)
 3darray yield_stock_temp(1,nst,1,nyr,1,nag)
 4darray SSB_region_temp(1,nst,1,nr,1,nyr,1,nag)
 matrix SSB_total_temp(1,nyr,1,nst)
 4darray abundance_stock_temp(1,nst,1,nyr,1,nag,1,nr)
 3darray abundance_total_temp(1,nyr,1,nag,1,nst)
 3darray biomass_stock_temp(1,nst,1,nyr,1,nr)
 matrix biomass_total_temp(1,nyr,1,nst)
 3darray catch_at_age_total_temp(1,nyr,1,nag,1,nst)
 matrix yield_total_temp(1,nyr,1,nst)
 4darray catch_at_age_stock_temp(1,nst,1,nyr,1,nag,1,nr)
 5darray SSB_region_temp_overlap(1,nst,1,nst,1,nr,1,nyr,1,nag)
 matrix abundance_move_overlap_temp(1,nst,1,nr)

 5darray yield_region_temp_overlap(1,nst,1,nst,1,nr,1,nyr,1,nag)
 4darray yield_stock_temp_overlap(1,nst,1,nst,1,nyr,1,nag)
 4darray abundance_natal_temp_overlap(1,nst,1,nyr,1,nag,1,nst)
 4darray biomass_stock_temp_overlap(1,nst,1,nst,1,nyr,1,nr)
 3darray biomass_natal_temp_overlap(1,nst,1,nyr,1,nst)
 4darray catch_at_age_natal_temp_overlap(1,nst,1,nyr,1,nag,1,nst)
 3darray yield_natal_temp_overlap(1,nst,1,nyr,1,nst)
 5darray catch_at_age_stock_temp_overlap(1,nst,1,nst,1,nyr,1,nag,1,nr)
 3darray SSB_natal_overlap_temp(1,nst,1,nyr,1,nst)
 matrix SSB_overlap_natal(1,nst,1,nr)
 5darray abundance_AM_overlap_region_all_natal_temp(1,nst,1,nr,1,nyr,1,nag,1,nst)
 3darray SSB_stock_temp(1,nst,1,nyr,1,nr)
 4darray SSB_stock_temp_overlap(1,nst,1,nst,1,nyr,1,nr)
 4darray res_TAC(1,nst,1,nr,1,nfl,1,nyr)
 3darray res_u(1,nst,1,nr,1,nyr)
 number Fnew
 number delt
 number fofF
 number fprimeF
 vector fofFvect(1,nag)
 vector fprimeFhigh(1,nag)
 vector fprimeFlow(1,nag)

 init_number dummy(phase_dummy)

  objective_function_value f

 !! cout << "parameters set" << endl;
 

INITIALIZATION_SECTION  //set initial values
     
  F_est .7;
  dummy 1;

PROCEDURE_SECTION

  get_movement();

  get_selectivity();

  get_F_age();

  get_vitals();

  get_SPR();
  
  get_env_Rec();

  //cout<<env_rec<<endl;

  get_abundance();

  //cout<<recruits_BM<<endl;
  //exit(43);
  
  get_rec_index();

  get_rand_CAA_prop();

  evaluate_the_objective_function();


///////BUILD MOVEMENT MATRIX////////
FUNCTION get_movement

//POSSIBLE ADDITIONS:
  //new functional forms
  //simulate yearly differences and age-based movement
  //random movement same as a random F value
  //random walk
  //density dependent
  //yearly deviations around an input average value (e.g., unidirectional movement with yearly deviations based on a variance term and lognormal distribution)
  //need to ensure random values still sum to 1 for an area (relativization of values)

  for (int j=1;j<=nstocks;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=nstocks;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
               if(a==1) // allow different movement from adults
                {
                 if(larval_move_switch==0)  //fix at no movement
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
                 if(larval_move_switch==1) // use input movement
                  {
                   T(j,r,y,a,k,n)=input_T(j,r,a,k,n);
                  }
                 if(larval_move_switch==2) // only allow movement within a stock (ie regions within a stock) not across stocks based on input residency term
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
                 if(larval_move_switch==3) //symmetric movement but only allow movement within a stock (ie regions within a stock) not across stocks
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
                 if(larval_move_switch==4) //symmetric movement across all stocks and regions
                  {
                   T(j,r,y,a,k,n)=1/sum(nregions);
                  }
                 if(larval_move_switch==5) // allow movement across all regions and stocks, based on stock/region specific residency
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
               else
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
                 if(move_switch==2) // only allow kovement within a stock (ie regions within a stock) not across stocks based on input residency term
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
                 if(move_switch==3) //symmetric movement but only allow movement within a stock (ie regions within a stock) not across stocks based on input residency term
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
                 if(move_switch==4) //symmetric movement across all stocks and regions
                  {
                   T(j,r,y,a,k,n)=1/sum(nregions);
                  }
                 if(move_switch==5) // allow movement across all regions and stocks, based on stock/region specific residency
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

  for (int j=1;j<=nstocks;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
        {
         for (int a=1;a<=nages;a++)
           {
            for (int z=1;z<=nfleets(j,r);z++)
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
              }
            }
          }
        }
      }


///////FISHING MORTALITY CALCULATIONS///////
FUNCTION get_F_age
   random_number_generator myrand(myseed);

 /// GOING TO NEED TIME-VARYING F
 // Year component to input F
 // Random variation/noise (e.g., TAC model) around an average value (or around input F values)
 //  rand_F.fill_randn(myrand1);
 //  for (int n=1;n<=nyrs_F_pre_TAC;n++) //prior to TAC fishery is assumed to fluctuate around scalar*FMSY if scalar=1.0 then fluctuates around FMSY otherwise can fluctuate around some hi or low F value to increase/decrease stock size relative to BMSY
 //  {
 //   SIM_F_devs(n)=mfexp(rand_F(n)*sigma_F-0.5*square(sigma_F));
 //   if(SIM_Fmsy_switch==1) //set F=input F for all years
 //   {
 //    SIM_F_pre_TAC(n)=Fmsy;
 //   }
 //   if(SIM_Fmsy_switch==0) //set F=scalar*Fmsy
 //   {    
 //    SIM_F_pre_TAC(n)=Fmsy_scalar*Fmsy*SIM_F_devs(n);
 //   }
 //  }
 // random walk


  for (int j=1;j<=nstocks;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int z=1;z<=nfleets(j,r);z++)
           { 
             if(F_switch==1) //input F directly
              {
               F_year(j,r,y,z)=input_F(j,r,z);
              }
             if(F_switch==2) //input single yearly F (i.e., using single stock FMSY for all stocks)
              {
               F_year(j,r,y,z)=input_F_MSY;
              }
             if(F_switch==3) //estimate F that achieves desired Total F
              {
               F_year(j,r,y,z)=F_est(j,r);
              }
             if(F_switch==4) //split F_MSY by nstocks
              {
               F_year(j,r,y,z)=input_F_MSY/nstocks;
              }
             if(F_switch==5) //split F_MSY by total number of regions
              {
               F_year(j,r,y,z)=input_F_MSY/sum(nregions);
              }
             if(F_switch==6) //split F_MSY by total number of fleets
              {
               F_year(j,r,y,z)=input_F_MSY/sum(nfleets);
              }
         //  if(F_switch==7) //random walk in F
         //   {
         //    F_year(j,r,y,z)=F_est(j,r,z);
         //   }
             F_fleet(j,r,y,a,z)=F_year(j,r,y,z)*selectivity(j,r,y,a,z);
             F(j,r,y,a)=sum(F_fleet(j,r,y,a)); 
             M(j,r,y,a)=input_M(j,a);
           }
          }
         }
        }
       }



FUNCTION get_vitals
//POSSIBLE ADDITIONS:
  //random walk in apportionment or random to give time-varying
  //switch for input recruitment devs by year to recreate a given stock trajectory
  //what if recruit index came before or after movement
 
  random_number_generator myrand(myseed);
  
  for (int p=1;p<=nstocks;p++)
   {
    for (int j=1;j<=nstocks;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int y=1;y<=nyrs;y++)
         {
          for (int a=1;a<=nages;a++)
           {
            for (int z=1;z<=nfleets(j,r);z++)
             {
              weight_stock(j,y,a)=input_weight(j,a);
              weight_catch(j,y,a)=input_catch_weight(j,a);
              
               if(recruit_devs_switch==0)  //use stock recruit relationship directly
                {
                 rec_devs(j,y)=1;
                }
               if(recruit_devs_switch==1)  // allow lognormal error around SR curve
                {
////// CHECK THIS CALCULATION MAKE SURE USING NEW RANDN every time through loop
////// rand_rec.fill_randn(myrand);
                 rec_devs(j,y)=mfexp(randn(myrand)*sigma_recruit(j)-.5*square(sigma_recruit(j)));
                }
               if(SSB_type==1) //fecundity based SSB
                {
                 wt_mat_mult(j,y,a)=0.5*fecundity(j,a)*maturity(j,a);
                }
               if(SSB_type==2) //weight based SSB
                {
                 wt_mat_mult(j,y,a)=0.5*weight_stock(j,y,a)*maturity(j,a);
                }
               if(apportionment_type==1) //input recruitment apportionment directly by stock and region
                {
                 Rec_Prop(j,r,y)=input_Rec_prop(j,r);
                }
               if(apportionment_type==2) //equal apportionment by nregions
                {
                 Rec_Prop(j,r,y)=1/nregions(j);
                }
               if(apportionment_type==(-1)) // no apportionment
                {
                 Rec_Prop(j,r,y)=1;
                }
               //if(apportionment_type==3)
                //{
                // Rec_Prop(j,r,y) = rec_index_BM(j,r,y)/sum(rec_index_BM(j,r,y));
                //}
             }
           }
         }
       }
     }
   }


///////SPR CALCS///////
FUNCTION get_SPR
  if(Rec_type==2) //BH recruitment
   {
    for (int k=1;k<=nstocks;k++)
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
     alpha(k)=SPR(k)*(1-steep(k))/(4*steep(k));
     beta(k)=(5*steep(k)-1)/(4*steep(k)*R_ave(k));
     SSB_zero(k)=SPR(k)*R_ave(k);
     }
   }


FUNCTION get_env_Rec // calculate temporally autocorrelated recruitment - input period and amplitude

   for (int y=1;y<=nyrs;y++)
         {
         for (int j=1;j<=nstocks;j++)
             {
            for (int r=1;r<=nregions(j);r++)
               {
          if(y==1)
             {
             env_rec(y)=init_abund(j,r,1);
             }
             
          if(y>1)
            {
            env_rec(y) = R_ave(j)+(R_ave(j)*amplitude)*sin(((2*PI)/freq)*years(y)+freq);
             }
           }
         }

       }


 //cout<<env_rec<<endl;
 //exit(43);


FUNCTION get_abundance
  random_number_generator myrand(myseed);

       for (int y=1;y<=nyrs;y++)
        {

//need to add all_natal calcs for all values that might want to compare across an area, because regular non-natal homing calcs do not account for different weights across natal populations so if if any part of calc involves weight it needs to use the biomass_all_natal value not just biomass_AM

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// Y==1 Overlap Calcs ///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

             if(y==1)
              {
               for (int a=1;a<=nages;a++)
                 {
           for (int p=1;p<=nstocks;p++)
            {
             for (int j=1;j<=nstocks;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j,r);z++)
                  {
                  if(p==j)
                   {
         /// ADD NATAL HOMING SPECIFIC ABUNDANCE CAN BE DISTRIBUTED ACROSS AREA (ie start in non-natal area)
         
                    abundance_at_age_BM_overlap_region(p,j,y,a,r)=init_abund(j,r,a);
                    abundance_at_age_BM_overlap_stock(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                   }
                  if(p!=j)
                   {
                    abundance_at_age_BM_overlap_region(p,j,y,a,r)=0;
                    abundance_at_age_BM_overlap_stock(p,j,y,a)=0;
                   }

                abundance_move_overlap_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                  
                    if(move_switch!=6  || a==1)  //if movement is not type=6 or a==1 (and movement type 6)
                     {
                      if(p==k)
                      {
                       abundance_move_overlap_temp(k,n)=init_abund(k,n,a)*T(p,n,y,a,j,r); //with overlap always use natal stock movement rates  (i.e., use p instead of k)
                      }
                     }
                   }
                  }

  
  
                    if(move_switch!=6 || move_switch!=7  || a==1) //if movement is not type=6,7 or a==1 (and movement type 6)
                     {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(p==j)
                       {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=init_abund(j,r,a);                    
                       }
                      }
                    if(move_switch==7  && a!=1)
                     {
                      if(p==j)
                      {
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=init_abund(j,r,a);
                      }
                     }
 ////////DO FOLLOWING CALCS FOR BM AS WELL                
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_stock(p,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
/////////////////////////////////////////////////////////////

                abundance_at_age_AM_overlap_stock(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));

                biomass_BM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                biomass_AM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_stock_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_stock_overlap(p,j,y)=sum(biomass_stock_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_stock_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

///////////MIGHT WANT TO DO FOLLOWING CALCS FOR NATAL HOMING AS WELL (DO METAPOP TYPE CALCS BELOW)
              //  abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
              //  abundance_res(j,r,y,a)=abundance_move_temp(j,r);
              //  abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
              //  bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
              //  bio_res(j,r,y,a)=bio_move_temp(j,r);
              //  bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_stock(j,y,a)-bio_res(j,r,y,a);
                


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance/////////////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////
    
                abundance_at_age_BM(j,r,y,a)=init_abund(j,r,a);
                recruits_BM(j,r,y)=abundance_at_age_BM(j,r,y,1);

            ////get year one recruitment index
               // rec_index_BM(j,r,y) = recruits_BM(j,r,y)*mfexp(randn(myrand)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));


                abundance_move_temp=0;
                bio_move_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=init_abund(k,n,a)*T(k,n,y,a,j,r);
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_stock(k,y,a);                   
                   }
                  }
                  if(overlap_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
                   {

 //// ADD IN BM CALCS TO ENSURE USING RIGHT VALUES FOR NATAL HOMING SCENARIOS
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);
/////////////////////////////
                   }
                  if(overlap_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_BM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_AM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }
                   
                abundance_stock_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_stock(j,y,a)=sum(abundance_stock_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_stock(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                biomass_stock_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_stock(j,y)=sum(biomass_stock_temp(j,y));
                biomass_total_temp(y,j)=biomass_stock(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                recruits_AM(j,r,y)=abundance_at_age_AM(j,r,y,a);
                abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                bio_res(j,r,y,a)=bio_move_temp(j,r);
                bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_stock(j,y,a)-bio_res(j,r,y,a);
                
              }
              }
              }
              }
              } //end age loop



 ///////////////////////////////////////////////////////////////////////////////////
 ////////////////NEWTON-RAPHSON for YEAR 1////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////
           for (int p=1;p<=nstocks;p++)
            {
             for (int j=1;j<=nstocks;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                   if(model_type_switch>0)  //Find the fully selected F that would result in the TAC
                     {
                      for (int a=1;a<=nages;a++)
                       {
                         for(int x=1;x<=nfleets(j,r);x++)
                           {
                             if(model_type_switch==1)
                                {
                                 if(overlap_switch>0)
                                   {
                                   if(input_TAC(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_TAC(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
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
 
                                            fofF=sum(fofFvect)-input_TAC(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
                                 if(overlap_switch==0)
                                   {
                                   if(input_TAC(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_TAC(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFhigh(s)=weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFlow(s)=weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFhigh(s)=weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFlow(s)=weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              }
                                             }

                                            fofF=sum(fofFvect)-input_TAC(j,r,x);
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
                                 if(overlap_switch>0)
                                   {
                                   if(input_u(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_u(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
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
 
                                            fofF=sum(fofFvect)-input_u(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
                                 if(overlap_switch==0)
                                   {
                                   if(input_u(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_u(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=(weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFhigh(s)=(weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFlow(s)=(weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=(weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFhigh(s)=(weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFlow(s)=(weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              }
                                             }

                                            fofF=sum(fofFvect)-input_u(j,r,x);
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
           for (int p=1;p<=nstocks;p++)
            {
             for (int j=1;j<=nstocks;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j,r);z++)
                  {
                abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));
                SSB_region_temp_overlap(p,j,r,y,a)=abundance_spawn_overlap(p,j,r,y,a)*wt_mat_mult(p,y,a); 
                SSB_region_overlap(p,j,r,y)=sum(SSB_region_temp_overlap(p,j,r,y));

                  SSB_overlap_natal=0;
                  if(overlap_switch==2)
                   {
                    for(int k=1;k<=nstocks;k++)
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
                       SSB_region_overlap(p,j,r,y)=SSB_region_overlap(p,j,r,y)+(sum(SSB_overlap_natal)/nregions(p));  //reutrning SSB is split evenly across natal regions
                       }
                   }
                   
                SSB_stock_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                SSB_stock_overlap(p,j,y)=sum(SSB_stock_temp_overlap(p,j,y));
                SSB_natal_overlap_temp(p,y,j)=SSB_stock_overlap(p,j,y);
                SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));

                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_stock(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                if(a==1)
                 {
                  catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F(j,r,y,a)+M(j,r,y,a))*(1-tspawn(p))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                 }
                if(a>1)
                 {
                  catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a))); //
                 }
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_stock_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_stock_overlap(p,j,y,a)=sum(catch_at_age_stock_temp_overlap(p,j,y,a));
                yield_stock_temp_overlap(p,j,y,a)=weight_catch(p,y,a)*catch_at_age_stock_overlap(p,j,y,a);
                yield_stock_overlap(p,j,y)=sum(yield_stock_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_stock_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_stock_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_stock_bio_overlap(p,j,y)=yield_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_stock_overlap(p,j,y)=biomass_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
                Bratio_stock_overlap(p,j,y)=SSB_stock_overlap(p,j,y)/SSB_zero(p);
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
                yield_fleet_temp(j,r,y,z,a)=weight_catch(j,y,a)*catch_at_age_fleet(j,r,y,a,z);
                yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                yield_region_temp(j,r,y,a)=weight_catch(j,y,a)*catch_at_age_region(j,r,y,a);
                yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                catch_at_age_stock_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                catch_at_age_stock(j,y,a)=sum(catch_at_age_stock_temp(j,y,a));
                yield_stock_temp(j,y,a)=weight_catch(j,y,a)*catch_at_age_stock(j,y,a);
                yield_stock(j,y)=sum(yield_stock_temp(j,y));
                catch_at_age_total_temp(y,a,j)=catch_at_age_stock(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_stock(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_stock_num(j,y,a)=catch_at_age_stock(j,y,a)/abundance_stock(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_stock_bio(j,y)=yield_stock(j,y)/biomass_stock(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_stock(j,y)=biomass_stock(j,y)/biomass_stock(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);

              if(overlap_switch>0)
               {
               if(p==j)
               {
                SSB_region(j,r,y)=SSB_region_overlap(p,j,r,y);  //if natal homing only account for SSB that is in its natal populations area, don't sum across natal populations
               }
              }
              if(overlap_switch==0)
              {
                SSB_region_temp(j,r,y,a)=abundance_spawn(j,r,y,a)*wt_mat_mult(j,y,a); 
                SSB_region(j,r,y)=sum(SSB_region_temp(j,r,y));
              }
                SSB_stock_temp(j,y,r)=SSB_region(j,r,y); 
                SSB_stock(j,y)=sum(SSB_stock_temp(j,y)); 
                SSB_total_temp(y,j)=SSB_stock(j,y);
                SSB_total(y)=sum(SSB_total_temp(y));
                Bratio_stock(j,y)=SSB_stock(j,y)/SSB_zero(j);
                Bratio_total(y)=SSB_total(y)/sum(SSB_zero);
               }
              }
             }
            }
              } //end age loop
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
           for (int p=1;p<=nstocks;p++)
            {
             for (int j=1;j<=nstocks;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j,r);z++)
                  {
               if(a==1)
                {
                 if(overlap_switch==1 || overlap_switch==2)
                 {
                 if(p==j)
                 {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to stock SSB but not region SSB)
                    }
                  }

 ////////////add recruits_BM with error (recruitment index) - and apportion switch 3
 
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=((SSB_stock_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_stock_overlap(p,j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=((SSB_stock_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_stock_overlap(p,j,y-1)))*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1)));
                    }
                  }
                  
                 if(Rec_type==3) //cyclical recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to stock SSB but not region SSB)
                    }
                  }
                 }
                 }

             if(overlap_switch==0)
              {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_stock(j,y-1));
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=((SSB_stock(j,y-1))/(alpha(j)+beta(j)*SSB_stock(j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=((SSB_stock(j,y-1))/(alpha(j)+beta(j)*SSB_stock(j,y-1)))*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_stock(j,y-1));
                    }

                 if(Rec_type==3) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_stock(j,y-1));
                    }
                  }
                  }
                 }


//////////////////finish filling out the recruitment index
    
                 //rec_index_BM(j,r,y) = recruits_BM(j,r,y)*mfexp(randn(myrand)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));


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


                abundance_at_age_BM_overlap_stock(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));

                abundance_move_overlap_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {

  ///////////NEED TO CHECK WHY USING FULL RECRUITMENT FUNCTIONS TO CALCULATE MOVEMENT OF AGE-1 FISH HERE, MIGHT BE ABLE TO CONDENSE CODE AND USE RECRUITS_BM OR ABUNDANCE_BM


                 if(overlap_switch==1 || overlap_switch==2)
                 {
                 if(p==k)
                 {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to stock SSB but not region SSB)
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_stock_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_stock_overlap(p,k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a stock
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_stock_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_stock_overlap(p,k,y-1)))*rec_devs(k,y)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r);
                    }
                  }
                 }
                 }
             if(overlap_switch==0)
              {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                     abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_stock(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                    abundance_move_overlap_temp(k,n)=((SSB_stock(k,y-1))/(alpha(k)+beta(k)*SSB_stock(k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a stock
                    {
                     abundance_move_overlap_temp(k,n)=((SSB_stock(k,y-1))/(alpha(k)+beta(k)*SSB_stock(k,y-1)))*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_stock(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                }
                   }
                  }

                abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_stock(p,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_stock(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));

                biomass_BM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                biomass_AM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_stock_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_stock_overlap(p,j,y)=sum(biomass_stock_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_stock_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////y>1, A==1 METAPOP TYPE CALCS (MOVEMENT)//////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
                 abundance_at_age_BM(j,r,y,a)=recruits_BM(j,r,y);

                    abundance_move_temp=0;
                    bio_move_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
             if(overlap_switch==0)
              {
                 if(Rec_type==1) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                     abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_stock(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                    abundance_move_temp(k,n)=((SSB_stock(k,y-1))/(alpha(k)+beta(k)*SSB_stock(k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a stock
                    {
                     abundance_move_temp(k,n)=((SSB_stock(k,y-1))/(alpha(k)+beta(k)*SSB_stock(k,y-1)))*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_stock(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }

                if(Rec_type==3) //cyclical recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 || apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regions within a stock
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a stock
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y)*(SSB_region(k,n,y-1)/SSB_stock(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                }
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_stock(k,y,a);
                   }
                  }
                  if(overlap_switch>0)
                   {
  //// ADD IN BM CALCS TO ENSURE USING RIGHT VALUES FOR NATAL HOMING SCENARIOS
                  
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)=biomass_AM_overlap_region_all_natal(j,r,y);
                    
                   }
                  if(overlap_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_BM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_AM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));

                   }


                biomass_stock_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_stock(j,y)=sum(biomass_stock_temp(j,y));
                biomass_total_temp(y,j)=biomass_stock(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                   abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                   abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                   abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                   bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                   bio_res(j,r,y,a)=bio_move_temp(j,r);
                   bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_stock(j,y,a)-bio_res(j,r,y,a);


            } //end a==1 if statement

 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(a==2) //account for partial year mortality during spawning year
                 {
                  abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p)));
                  abundance_at_age_BM_overlap_stock(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));

                abundance_move_overlap_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7  ||a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(p)))*T(p,n,y,a,j,r); //with overlap always use natal stock movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal stock movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(p)))*return_probability(p); //with overlap always use natal stock movement rates
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
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(p))); //with overlap always use natal stock movement rates                     
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
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_stock(p,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_stock(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));

                biomass_BM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                biomass_AM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_stock_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_stock_overlap(p,j,y)=sum(biomass_stock_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_stock_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

 //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
   ///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
            abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(j))); //account for time of spawning (i.e., if born midway only experience a half year of mortality from age-1 to age-2)

                    abundance_move_temp=0;
                     bio_move_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=abundance_at_age_AM(k,n,y-1,a-1)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1))*(1-tspawn(k)))*T(k,n,y,a,j,r);
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_stock(k,y,a);
                   }
                  }

                  if(overlap_switch>0)
                   {
  //// ADD IN BM CALCS TO ENSURE USING RIGHT VALUES FOR NATAL HOMING SCENARIOS
                  
                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);
                  
                   }
                  if(overlap_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_BM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_AM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));

                   }

                biomass_stock_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_stock(j,y)=sum(biomass_stock_temp(j,y));
                biomass_total_temp(y,j)=biomass_stock(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));


                   abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                   abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                   abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                   bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                   bio_res(j,r,y,a)=bio_move_temp(j,r);
                   bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_stock(j,y,a)-bio_res(j,r,y,a);

             } //end a==2 if statement

 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(a>2 && a<nages)
                 {
                  abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
                  abundance_at_age_BM_overlap_stock(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                  
                abundance_move_overlap_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7 || a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(p,n,y,a,j,r); //with overlap always use natal stock movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal stock movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*return_probability(p); //with overlap always use natal stock movement rates
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
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))); //with overlap always use natal stock movement rates                     
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
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_stock(p,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_stock(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));

                biomass_BM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                biomass_AM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_stock_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_stock_overlap(p,j,y)=sum(biomass_stock_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_stock_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));
                
 //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
   ///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
                  abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));

                   abundance_move_temp=0;
                   bio_move_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=abundance_at_age_AM(k,n,y-1,a-1)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(k,n,y,a,j,r);
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_stock(k,y,a);
                   }
                  }
                  
                  if(overlap_switch>0)
                   {

  //// ADD IN BM CALCS TO ENSURE USING RIGHT VALUES FOR NATAL HOMING SCENARIOS

                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);

                   }
                  if(overlap_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                    biomass_BM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_BM(j,r,y,a);
                    biomass_AM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_AM(j,r,y,a);
                    biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                    biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }


                biomass_stock_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_stock(j,y)=sum(biomass_stock_temp(j,y));
                biomass_total_temp(y,j)=biomass_stock(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                   abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                   abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                   abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                   bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                   bio_res(j,r,y,a)=bio_move_temp(j,r);
                   bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_stock(j,y,a)-bio_res(j,r,y,a);

          } //end a>2 <nages if statement

 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(a==nages) //account for fish already in plus group
                 {
                  abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
                  abundance_at_age_BM_overlap_stock(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                  
                abundance_move_overlap_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                    if(move_switch!=6  || move_switch!=7  || a==1)
                     {
                      abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(p,n,y,a,j,r)+abundance_at_age_AM_overlap_region(p,k,y-1,a,n)*mfexp(-(M(k,n,y-1,a)+F(k,n,y-1,a)))*T(p,n,y,a,j,r); //with overlap always use natal stock movement rates
                     }
                    if(move_switch==6 && a>1)
                     {
                      if(a==return_age && p==j && p==k && j==k)
                      {
                       abundance_move_overlap_temp(k,n)=0; //with overlap always use natal stock movement rates
                      }
                      if(a==return_age && p==j && j!=k)
                      {
                       abundance_move_overlap_temp(k,n)=abundance_at_age_AM_overlap_region(p,k,y-1,a-1,n)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*return_probability(p)+abundance_at_age_AM_overlap_region(p,k,y-1,a,n)*mfexp(-(M(k,n,y-1,a)+F(k,n,y-1,a)))*return_probability(p); //with overlap always use natal stock movement rates
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
                        abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a))); //with overlap always use natal stock movement rates                     
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
                biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_stock(p,y,a);
                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a));
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                abundance_at_age_AM_overlap_stock(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));

                biomass_BM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                biomass_AM_age_overlap(p,j,r,y,a)=weight_stock(p,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r);
                biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_stock_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_stock_overlap(p,j,y)=sum(biomass_stock_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_stock_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

 //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
   ///////////////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
               abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM(j,r,y-1,a)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
                  
                   abundance_move_temp=0;
                   bio_move_temp=0;

                for (int k=1;k<=nstocks;k++)
                 {
                  for (int n=1;n<=nregions(k);n++)
                   {
                     abundance_move_temp(k,n)=abundance_at_age_AM(k,n,y-1,a-1)*mfexp(-(M(k,n,y-1,a-1)+F(k,n,y-1,a-1)))*T(k,n,y,a,j,r)+abundance_at_age_AM(k,n,y-1,a)*mfexp(-(M(k,n,y-1,a)+F(k,n,y-1,a)))*T(k,n,y,a,j,r);;
                     bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_stock(k,y,a);
                   }
                  }
                  if(overlap_switch>0)
                   {
  //// ADD IN BM CALCS TO ENSURE USING RIGHT VALUES FOR NATAL HOMING SCENARIOS

                    abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
                    biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);

                   }
                  if(overlap_switch==0)
                   {
                    abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
                biomass_BM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_BM(j,r,y,a);
                biomass_AM_age(j,r,y,a)=weight_stock(j,y,a)*abundance_at_age_AM(j,r,y,a);
                biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
                biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
                   }

                biomass_stock_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_stock(j,y)=sum(biomass_stock_temp(j,y));
                biomass_total_temp(y,j)=biomass_stock(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

                   abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                   abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                   abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                   bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                   bio_res(j,r,y,a)=bio_move_temp(j,r);
                   bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_stock(j,y,a)-bio_res(j,r,y,a);

        } //end nages if statement
       }
       }
       }
       }
     } //end age loop


 ///////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////
 /////////////////NEWTON RAPHSON CALCS//////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////

          for (int p=1;p<=nstocks;p++)
            {
             for (int j=1;j<=nstocks;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                   if(model_type_switch>0)  //Find the fully selected F that would result in the TAC
                     {
                      for (int a=1;a<=nages;a++)
                       {
                         for(int x=1;x<=nfleets(j,r);x++)
                           {
                             if(model_type_switch==1)
                                {
                                 if(overlap_switch>0)
                                   {
                                   if(input_TAC(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_TAC(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
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
                                            fofF=sum(fofFvect)-input_TAC(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
    ///////////////////////////////////////////////
                                 if(overlap_switch==0)
                                   {
                                   if(input_TAC(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_TAC(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFhigh(s)=weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));
                                               fprimeFlow(s)=weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j))));                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFhigh(s)=weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              fprimeFlow(s)=weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                                              }
                                             }

                                            fofF=sum(fofFvect)-input_TAC(j,r,x);
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
                                 if(overlap_switch>0)
                                   {
                                   if(input_u(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_u(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
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
                                            
                                            fofF=sum(fofFvect)-input_u(j,r,x);
                                            fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                                            Fnew=Fnew-(fofF/fprimeF);
                                          if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                                          if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                                         }
                                         } 
                                       }
    ////////////////////////////////                               
                                 if(overlap_switch==0)
                                   {
                                   if(input_u(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                                      {
                                       Fnew=0;
                                      }
                                   if(input_u(j,r,x)>0)
                                      {
                                        Fnew=Fnew_start;
                                        for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
                                         {
                                          delt=Fnew*NR_dev;  // NR_dev~0.001
                                           for(int s=1;s<=nages;s++)
                                            {
                                              if(s==1)
                                              {
                                               fofFvect(s)=(weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFhigh(s)=(weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);
                                               fprimeFlow(s)=(weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))*(1-tspawn(j)))))/biomass_AM(j,r,y);                                              
                                              }
                                              if(s>1)
                                              {
                                              fofFvect(s)=(weight_catch(j,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFhigh(s)=(weight_catch(j,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              fprimeFlow(s)=(weight_catch(j,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                                              }
                                             }

                                            fofF=sum(fofFvect)-input_u(j,r,x);
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
           for (int p=1;p<=nstocks;p++)
            {
             for (int j=1;j<=nstocks;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j,r);z++)
                  {
          if(a==1)
            {

                abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));

                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_stock(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))*(1-tspawn(j))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_stock_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_stock_overlap(p,j,y,a)=sum(catch_at_age_stock_temp_overlap(p,j,y,a));
                yield_stock_temp_overlap(p,j,y,a)=weight_catch(p,y,a)*catch_at_age_stock_overlap(p,j,y,a);
                yield_stock_overlap(p,j,y)=sum(yield_stock_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_stock_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_stock_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_stock_bio_overlap(p,j,y)=yield_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_stock_overlap(p,j,y)=biomass_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////a==1 metapop type abundance calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



                   recruits_AM(j,r,y)=abundance_at_age_AM(j,r,y,a);
 
                   abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));

                   catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))*(1-tspawn(j))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a)); //account for time of spawning in catch (tspawn divides out in F/(F+M

                   //catchsum(j,r,y,z)=sum(catch_at_age_fleet(j,r,y,z));
                   
                   yield_fleet_temp(j,r,y,z,a)=weight_catch(j,y,a)*catch_at_age_fleet(j,r,y,a,z);
                   yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                   catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                   yield_region_temp(j,r,y,a)=weight_catch(j,y,a)*catch_at_age_region(j,r,y,a);
                   yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                   catch_at_age_stock_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                   catch_at_age_stock(j,y,a)=sum(catch_at_age_stock_temp(j,y,a));
                   yield_stock_temp(j,y,a)=weight_catch(j,y,a)*catch_at_age_stock(j,y,a);
                   yield_stock(j,y)=sum(yield_stock_temp(j,y));

                abundance_stock_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_stock(j,y,a)=sum(abundance_stock_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_stock(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_stock(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_stock(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_stock_num(j,y,a)=catch_at_age_stock(j,y,a)/abundance_stock(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_stock_bio(j,y)=yield_stock(j,y)/biomass_stock(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_stock(j,y)=biomass_stock(j,y)/biomass_stock(j,1);
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

                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_stock(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_stock_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_stock_overlap(p,j,y,a)=sum(catch_at_age_stock_temp_overlap(p,j,y,a));
                yield_stock_temp_overlap(p,j,y,a)=weight_catch(p,y,a)*catch_at_age_stock_overlap(p,j,y,a);
                yield_stock_overlap(p,j,y)=sum(yield_stock_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_stock_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_stock_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_stock_bio_overlap(p,j,y)=yield_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_stock_overlap(p,j,y)=biomass_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
              
                 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////a==2 metapop type abundance calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
                   
                  abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));
                  yield_fleet_temp(j,r,y,z,a)=weight_catch(j,y,a)*catch_at_age_fleet(j,r,y,a,z);
                  yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                  catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                  yield_region_temp(j,r,y,a)=weight_catch(j,y,a)*catch_at_age_region(j,r,y,a);
                  yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                  catch_at_age_stock_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                  catch_at_age_stock(j,y,a)=sum(catch_at_age_stock_temp(j,y,a));
                  yield_stock_temp(j,y,a)=weight_catch(j,y,a)*catch_at_age_stock(j,y,a);
                  yield_stock(j,y)=sum(yield_stock_temp(j,y));

                abundance_stock_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_stock(j,y,a)=sum(abundance_stock_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_stock(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_stock(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_stock(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_stock_num(j,y,a)=catch_at_age_stock(j,y,a)/abundance_stock(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_stock_bio(j,y)=yield_stock(j,y)/biomass_stock(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_stock(j,y)=biomass_stock(j,y)/biomass_stock(j,1);
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

                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_stock(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_stock_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_stock_overlap(p,j,y,a)=sum(catch_at_age_stock_temp_overlap(p,j,y,a));
                yield_stock_temp_overlap(p,j,y,a)=weight_catch(p,y,a)*catch_at_age_stock_overlap(p,j,y,a);
                yield_stock_overlap(p,j,y)=sum(yield_stock_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_stock_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_stock_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_stock_bio_overlap(p,j,y)=yield_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_stock_overlap(p,j,y)=biomass_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
                
                 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////a>2 a<nages metapop type abundance calcs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   
                  abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));

                  //catchsum(j,r,y,z)=sum(catch_at_age_fleet(j,r,y,z));

                  yield_fleet_temp(j,r,y,z,a)=weight_catch(j,y,a)*catch_at_age_fleet(j,r,y,a,z);
                  yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                  catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                  yield_region_temp(j,r,y,a)=weight_catch(j,y,a)*catch_at_age_region(j,r,y,a);
                  yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                  catch_at_age_stock_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                  catch_at_age_stock(j,y,a)=sum(catch_at_age_stock_temp(j,y,a));
                  yield_stock_temp(j,y,a)=weight_catch(j,y,a)*catch_at_age_stock(j,y,a);
                  yield_stock(j,y)=sum(yield_stock_temp(j,y));

                abundance_stock_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_stock(j,y,a)=sum(abundance_stock_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_stock(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_stock(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_stock(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_stock_num(j,y,a)=catch_at_age_stock(j,y,a)/abundance_stock(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_stock_bio(j,y)=yield_stock(j,y)/biomass_stock(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_stock(j,y)=biomass_stock(j,y)/biomass_stock(j,1);
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

                abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_stock(p,j,y,a);
                abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));
                catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-exp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a)));
                yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,y,a)*catch_at_age_region_overlap(p,j,r,y,a);
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                catch_at_age_stock_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
                catch_at_age_stock_overlap(p,j,y,a)=sum(catch_at_age_stock_temp_overlap(p,j,y,a));
                yield_stock_temp_overlap(p,j,y,a)=weight_catch(p,y,a)*catch_at_age_stock_overlap(p,j,y,a);
                yield_stock_overlap(p,j,y)=sum(yield_stock_temp_overlap(p,j,y));
                catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_stock_overlap(p,j,y,a);
                catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
                yield_natal_temp_overlap(p,y,j)=yield_stock_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_stock_bio_overlap(p,j,y)=yield_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1);
                depletion_stock_overlap(p,j,y)=biomass_stock_overlap(p,j,y)/biomass_stock_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);               
                 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////a==nages metapop type abundance calcs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

           
                   
                  abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));
                  catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));

                  //catchsum(j,r,y,z)=sum(catch_at_age_fleet(j,r,y,z));
                   
                  yield_fleet_temp(j,r,y,z,a)=weight_catch(j,y,a)*catch_at_age_fleet(j,r,y,a,z);
                  yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));
                  catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                  yield_region_temp(j,r,y,a)=weight_catch(j,y,a)*catch_at_age_region(j,r,y,a);
                  yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                  catch_at_age_stock_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                  catch_at_age_stock(j,y,a)=sum(catch_at_age_stock_temp(j,y,a));
                  yield_stock_temp(j,y,a)=weight_catch(j,y,a)*catch_at_age_stock(j,y,a);
                  yield_stock(j,y)=sum(yield_stock_temp(j,y));

                abundance_stock_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_stock(j,y,a)=sum(abundance_stock_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_stock(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));
                catch_at_age_total_temp(y,a,j)=catch_at_age_stock(j,y,a);
                catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                yield_total_temp(y,j)=yield_stock(j,y);
                yield_total(y)=sum(yield_total_temp(y));
                harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                harvest_rate_stock_num(j,y,a)=catch_at_age_stock(j,y,a)/abundance_stock(j,y,a);
                harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_stock_bio(j,y)=yield_stock(j,y)/biomass_stock(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1);
                depletion_stock(j,y)=biomass_stock(j,y)/biomass_stock(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);
        } //end nages if statement


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////SSB calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   //don't adjust SSB for fish not in natal stock area because SR calculations do this automatically by only using
   //SSB that is in natal population area (i.e., p==j)
   //for spawning migrations scenarios the SSB is adjusted outside the loop to remove SSB of fish that actually returned
   // to natal population (i.e., remove this SSB from the non-natal areas...doesn't impact SR calcs so can do this outside
   //loops without consequence to model
   
                SSB_region_temp_overlap(p,j,r,y,a)=abundance_spawn_overlap(p,j,r,y,a)*wt_mat_mult(p,y,a); 
                SSB_region_overlap(p,j,r,y)=sum(SSB_region_temp_overlap(p,j,r,y));
                 SSB_overlap_natal=0;
                  if(overlap_switch==2) //spawning return calculations
                   {
                    for(int k=1;k<=nstocks;k++)
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
                SSB_stock_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                SSB_stock_overlap(p,j,y)=sum(SSB_stock_temp_overlap(p,j,y));
                SSB_natal_overlap_temp(p,y,j)=SSB_stock_overlap(p,j,y);
                SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));  /// this is adjusted below outside y loop to account for fish not spawning

              if(overlap_switch>0)
               {
               if(p==j)  //accounts for not being in natal area
               {
                SSB_region(j,r,y)=SSB_region_overlap(p,j,r,y);
               }
              }
              if(overlap_switch==0)
              {
                SSB_region_temp(j,r,y,a)=abundance_spawn(j,r,y,a)*wt_mat_mult(j,y,a); 
                SSB_region(j,r,y)=sum(SSB_region_temp(j,r,y));
              }
                SSB_stock_temp(j,y,r)=SSB_region(j,r,y);
                SSB_stock(j,y)=sum(SSB_stock_temp(j,y));
                SSB_total_temp(y,j)=SSB_stock(j,y);
                SSB_total(y)=sum(SSB_total_temp(y));

                Bratio_stock_overlap(p,j,y)=SSB_stock_overlap(p,j,y)/SSB_zero(p);
                Bratio_natal_overlap(p,y)=SSB_natal_overlap(p,y)/SSB_zero(p);
                Bratio_stock(j,y)=SSB_stock(j,y)/SSB_zero(j);
                Bratio_total(y)=SSB_total(y)/sum(SSB_zero);
   }
   }
   }
   }
  } // end age loop
 } //end yr>1 loop
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     } //end y loop


  /////now adjusting natal homing ssb for spawning migration scenarios where fish left non-natal area
       for (int p=1;p<=nstocks;p++)
        {
         for (int j=1;j<=nstocks;j++)
          {
           for (int r=1;r<=nregions(p);r++)
            {
             for (int t=1;t<=nregions(j);t++)
              {
               for (int y=1;y<=nyrs;y++)
                {
                 SSB_overlap_natal=0;
                  if(overlap_switch==2)
                    {
                    if(p!=j)  //update SSB that doesn't spawn (ie doesn't return to natal stock)
                    {
                     SSB_region_overlap(p,j,r,y)=(1-spawn_return_prob(p))*SSB_region_overlap(p,j,r,y);
                    }
                   SSB_stock_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                   SSB_stock_overlap(p,j,y)=sum(SSB_stock_temp_overlap(p,j,y));
                   SSB_natal_overlap_temp(p,y,j)=SSB_stock_overlap(p,j,y);
                   SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));
                Bratio_stock_overlap(p,j,y)=SSB_stock_overlap(p,j,y)/SSB_zero(p);
                Bratio_natal_overlap(p,y)=SSB_natal_overlap(p,y)/SSB_zero(p);
                Bratio_stock(j,y)=SSB_stock(j,y)/SSB_zero(j);
                Bratio_total(y)=SSB_total(y)/sum(SSB_zero);
                    }
                   }
                  }
                 }
                }
               }



 

FUNCTION get_rec_index //function to calculate recruitment index

// can add after movement too for larval drift scenario
  random_number_generator myrand(myseed);

  for (int j=1;j<=nstocks;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        rec_index_BM(j,r,y) = recruits_BM(j,r,y)*mfexp(randn(myrand)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));
       }
     }
   }




FUNCTION get_rand_CAA_prop
 random_number_generator myrand(myseed);
 //calculate the total CAA for each stock/area/region/fleet for non-natal and natal calcs
 
  for (int p=1;p<=nstocks;p++)
   {
    for (int j=1;j<=nstocks;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int y=1;y<=nyrs;y++)
         {
          for (int z=1;z<=nfleets(j,r);z++)
             {
           for (int a=1;a<=nages;a++)
               {
               
            if(overlap_switch==0)//for non-natal homing -
              {
               true_caa_fleet_temp(j,r,y,z,a)=catch_at_age_fleet(j,r,y,a,z); //summing catch across age
               true_caa_fleet_total(j,r,y,z) = sum(true_caa_fleet_temp(j,r,y,z));
               //calc random catches assuming true are means with log-normal age-specific error
               obs_caa_fleet(j,r,y,z,a)= true_caa_fleet_temp(j,r,y,z,a)*mfexp(randn(myrand)*caa_sigma(a)-0.5*square(caa_sigma(a)));
               obs_caa_fleet_total(j,r,y,z) = sum(obs_caa_fleet(j,r,y,z));
               
             for (int a=1;a<=nages;a++)
               {
               CAA_prop_fleet(j,r,y,z,a) =  true_caa_fleet_temp(j,r,y,z,a)/true_caa_fleet_total(j,r,y,z);
               //observed age comp
               obs_prop_fleet(j,r,y,z,a) = obs_caa_fleet(j,r,y,z,a)/obs_caa_fleet_total(j,r,y,z);
               }
              }
               
            if(overlap_switch>0)//for natal-homing...havent really tested this
              {
               true_caa_overlap_reg_temp(p,j,r,y,a) = catch_at_age_region_overlap(p,j,r,y,a);
               true_caa_overlap_reg_total(p,j,r,y) = sum(true_caa_overlap_reg_temp(p,j,r,y));
               //calc random catches assuming true are means with log-normal age-specific error
               obs_caa_overlap_reg(p,j,r,y,a)= true_caa_overlap_reg_temp(p,j,r,y,a)*mfexp(randn(myrand)*caa_sigma(a)-0.5*square(caa_sigma(a)));
               obs_caa_overlap_reg_total(p,j,r,y) = sum(obs_caa_overlap_reg(p,j,r,y));


              for (int a=1;a<=nages;a++)
               {
               CAA_prop_overlap(p,j,r,y,a) = true_caa_overlap_reg_temp(p,j,r,y,a)/true_caa_overlap_reg_total(p,j,r,y);
               //observed age comp
               obs_prop_overlap_reg(p,j,r,y,a)= obs_caa_overlap_reg(p,j,r,y,a)/obs_caa_overlap_reg_total(p,j,r,y);
               }
              }
             }
            }
          }  
        }
      }
    }


  //cout<<obs_caa_fleet<<endl;
  //cout<<obs_caa_fleet_total<<endl;
  //exit(43);



FUNCTION evaluate_the_objective_function
   f=0.0;

           for (int y=1;y<=nyrs;y++)
            {
             for (int j=1;j<=nstocks;j++)
              {
               for (int r=1;r<=nregions(j);r++)
                {
                 for (int z=1;z<=nfleets(j,r);z++)
                  {
                    res_TAC(j,r,z,y)=input_TAC(j,r,z)-yield_fleet(j,r,y,z);
                  }
                 res_u(j,r,y)=sum(input_u(j,r))-harvest_rate_region_bio(j,r,y);
                }
               }
              }


REPORT_SECTION
  report<<"$res_TAC"<<endl;
  report<<res_TAC<<endl;
  report<<"$res_u"<<endl;
  report<<res_u<<endl;
  report<<"$myseed"<<endl;
  report<<myseed<<endl;
  report<<"$nages"<<endl;
  report<<nages<<endl;
  report<<"$nyrs"<<endl;
  report<<nyrs<<endl;
  report<<"$nstocks"<<endl;
  report<<nstocks<<endl;
  report<<"$nregions"<<endl;
  report<<nregions<<endl;
  report<<"$nages"<<endl;
  report<<nages<<endl;
  report<<"$nfleets"<<endl;
  report<<nfleets<<endl;

  report<<"$larval_move_switch"<<endl;
  report<<larval_move_switch<<endl;
  report<<"$move_switch"<<endl;
  report<<move_switch<<endl;
  report<<"$overlap_switch"<<endl;
  report<<overlap_switch<<endl;
  report<<"$F_switch"<<endl;
  report<<F_switch<<endl;
  report<<"$recruit_devs_switch"<<endl;
  report<<recruit_devs_switch<<endl;
  report<<"$SSB_type"<<endl;
  report<<SSB_type<<endl;
  report<<"$apportionment_type"<<endl;
  report<<apportionment_type<<endl;
  report<<"$Rec_type"<<endl;
  report<<Rec_type<<endl;
  report<<"$return_age"<<endl;
  report<<return_age<<endl;
  report<<"$return_probability"<<endl;
  report<<return_probability<<endl;
  report<<"$spawn_return_prob"<<endl;
  report<<spawn_return_prob<<endl;

  report<<"$tspawn"<<endl;
  report<<tspawn<<endl;
  report<<"$steep"<<endl;
  report<<steep<<endl;
  report<<"$R_ave"<<endl;
  report<<R_ave<<endl;
  report<<"$SSB_zero"<<endl;
  report<<SSB_zero<<endl;




//  report<<"$input_T"<<endl;
//  report<<input_T<<endl;
//  report<<"$input_residency"<<endl;
//  report<<input_residency<<endl;
//  report<<"$input_residency_larval"<<endl;
//  report<<input_residency_larval<<endl;
//  report<<"$sel_beta1"<<endl;
//  report<<sel_beta1<<endl;
//  report<<"$sel_beta2"<<endl;
//  report<<sel_beta2<<endl;
//  report<<"$sel_beta3"<<endl;
//  report<<sel_beta3<<endl;
//  report<<"$sel_beta4"<<endl;
//  report<<sel_beta4<<endl;
//  report<<"$input_selectivity"<<endl;
//  report<<input_selectivity<<endl;
//  report<<"$input_F"<<endl;
//  report<<input_F<<endl;
//  report<<"$input_F_MSY"<<endl;
//  report<<input_F_MSY<<endl;
//  report<<"$input_M"<<endl;
//  report<<input_M<<endl;
//  report<<"$sigma_recruit"<<endl;
//  report<<sigma_recruit<<endl;
//  report<<"$input_weight"<<endl;
//  report<<input_weight<<endl;
//  report<<"$input_catch_weight"<<endl;
//  report<<input_catch_weight<<endl;
//  report<<"$fecundity"<<endl;
//  report<<fecundity<<endl;
//  report<<"$maturity"<<endl;
//  report<<maturity<<endl;
//  report<<"$input_Rec_prop"<<endl;
//  report<<input_Rec_prop<<endl;
//  report<<"$init_abund"<<endl;
//  report<<init_abund<<endl;
//  report<<$"rec_index_sigma"<<end;
//  report<<rec_index_sigma<<end;

//  report<<"$F_est"<<endl;
//  report<<F_est<<endl;


//  report<<"$rec_devs"<<endl;
 // report<<rec_devs<<endl;
//  report<<"$weight_stock"<<endl;
//  report<<weight_stock<<endl;
//  report<<"$weight_catch"<<endl;
//  report<<weight_catch<<endl;
//  report<<"$wt_mat_mult"<<endl;
//  report<<wt_mat_mult<<endl;
//  report<<"$weight_stock_overlap"<<endl;
//  report<<weight_stock_overlap<<endl;
//  report<<"$weight_catch_overlap"<<endl;
//  report<<weight_catch_overlap<<endl;
//  report<<"$wt_mat_mult_overlap"<<endl;
//  report<<wt_mat_mult_overlap<<endl;
//  report<<"$Rec_Prop"<<endl;
//  report<<Rec_Prop<<endl;
//  report<<"$SPR_N"<<endl;
//  report<<SPR_N<<endl;
//  report<<"$SPR_SSB"<<endl;
// report<<SPR_SSB<<endl;

  report<<"$alpha"<<endl;
  report<<alpha<<endl;
  report<<"$beta"<<endl;
  report<<beta<<endl;


//  report<<"$abundance_at_age_BM"<<endl;
//  report<<abundance_at_age_BM<<endl;
//  report<<"$abundance_at_age_AM"<<endl;
//  report<<abundance_at_age_AM<<endl;
//  report<<"T"<<endl;
//  report<<T<<endl;
//  report<<"$abundance_AM_overlap_region_all_natal"<<endl;
//  report<<abundance_AM_overlap_region_all_natal<<endl;
//  report<<"$abundance_stock"<<endl;
//  report<<abundance_stock<<endl;
//  report<<"$abundance_total"<<endl;
//  report<<abundance_total<<endl;
//  report<<"$abundance_in"<<endl;
//  report<<abundance_in<<endl;
//  report<<"$abundance_res"<<endl;
//  report<<abundance_res<<endl;
//  report<<"$abundance_leave"<<endl;
//  report<<abundance_leave<<endl;
//  report<<"$abundance_spawn"<<endl;
//  report<<abundance_spawn<<endl;
//  report<<"$biomass_BM"<<endl;
//  report<<biomass_BM<<endl;
//  report<<"$biomass_BM_age"<<endl;
//  report<<biomass_BM_age<<endl;

//  report<<"$biomass_AM_age"<<endl;
//  report<<biomass_AM_age<<endl;

//  report<<"$bio_in"<<endl;
//  report<<bio_in<<endl;
//  report<<"$bio_res"<<endl;
//  report<<bio_res<<endl;
//  report<<"$bio_leave"<<endl;
//  report<<bio_leave<<endl;

//  report<<"$catch_at_age_fleet"<<endl;
//  report<<catch_at_age_fleet<<endl;
//  report<<"$catch_at_age_region"<<endl;
//  report<<catch_at_age_region<<endl;
//  report<<"$catch_at_age_stock"<<endl;
//  report<<catch_at_age_stock<<endl;
//  report<<"$catch_at_age_total"<<endl;
//  report<<catch_at_age_total<<endl;

//  report<<"$obs_caa_fleet"<<endl;
//  report<<obs_caa_fleet<<endl;
//  report<<"$obs_prop_fleet"<<endl;
//  report<<obs_prop_fleet<<endl;

//  report<<"$obs_caa_overlap_reg"<<endl;
//  report<<obs_caa_overlap_reg<<endl;
//  report<<"$obs_prop_overlap_reg"<<endl;
//  report<<obs_prop_overlap_reg<<endl;

//  report<<"$recruits_BM"<<endl;
//  report<<recruits_BM<<endl;
//  report<<"$rec_index_BM"<<endl;
//  report<<rec_index_BM<<end;


//  report<<"$harvest_rate_region_num"<<endl;
//  report<<harvest_rate_region_num<<endl;
//  report<<"$harvest_rate_stock_num"<<endl;
//  report<<harvest_rate_stock_num<<endl;
//  report<<"$harvest_rate_total_num"<<endl;
//  report<<harvest_rate_total_num<<endl;



//  report<<"$residuals"<<endl;
//  report<<res<<endl;
//  report<<"$harvest_rate_pen"<<endl;
//  report<<harvest_rate_pen<<endl;

//  report<<"$abundance_at_age_BM_overlap_region"<<endl;
//  report<<abundance_at_age_BM_overlap_region<<endl;
//  report<<"$abundance_at_age_BM_overlap_stock"<<endl;
//  report<<abundance_at_age_BM_overlap_stock<<endl;
//  report<<"$abundance_at_age_AM_overlap_region"<<endl;
//  report<<abundance_at_age_AM_overlap_region<<endl;
//  report<<"$abundance_at_age_AM_overlap_stock"<<endl;
//  report<<abundance_at_age_AM_overlap_stock<<endl;
//  report<<"$abundance_natal_overlap"<<endl;
//  report<<abundance_natal_overlap<<endl;
//  report<<"$abundance_spawn_overlap"<<endl;
//  report<<abundance_spawn_overlap<<endl;



//  report<<"$catch_at_age_region_overlap"<<endl;
//  report<<catch_at_age_region_overlap<<endl;
//  report<<"$catch_at_age_stock_overlap"<<endl;
//  report<<catch_at_age_stock_overlap<<endl;
//  report<<"$catch_at_age_natal_overlap"<<endl;
//  report<<catch_at_age_natal_overlap<<endl;



//  report<<"$biomass_BM_overlap_region"<<endl;
//  report<<biomass_BM_overlap_region<<endl;
//  report<<"$biomass_BM_age_overlap"<<endl;
//  report<<biomass_BM_age_overlap<<endl;

 // report<<"$biomass_AM_age_overlap"<<endl;
//  report<<biomass_AM_age_overlap<<endl;


  report<<"$F"<<endl;
  report<<F<<endl;
  report<<"$M"<<endl;
  report<<M<<endl;

  report<<"$biomass_AM"<<endl;
  report<<biomass_AM<<endl;
  report<<"$biomass_stock"<<endl;
  report<<biomass_stock<<endl;
  report<<"$biomass_total"<<endl;
  report<<biomass_total<<endl;

  report<<"$yield_fleet"<<endl;
  report<<yield_fleet<<endl;
  report<<"$yield_region"<<endl;
  report<<yield_region<<endl;
  report<<"$yield_stock"<<endl;
  report<<yield_stock<<endl;
  report<<"$yield_total"<<endl;
  report<<yield_total<<endl;

  report<<"$harvest_rate_region_bio"<<endl;
  report<<harvest_rate_region_bio<<endl;
  report<<"$harvest_rate_stock_bio"<<endl;
  report<<harvest_rate_stock_bio<<endl;
  report<<"$harvest_rate_total_bio"<<endl;
  report<<harvest_rate_total_bio<<endl;

  report<<"$depletion_region"<<endl;
  report<<depletion_region<<endl;
  report<<"$depletion_stock"<<endl;
  report<<depletion_stock<<endl;
  report<<"$depletion_total"<<endl;
  report<<depletion_total<<endl;

  report<<"$SSB_region"<<endl;
  report<<SSB_region<<endl;
  report<<"$SSB_stock"<<endl;
  report<<SSB_stock<<endl;
  report<<"$SSB_total"<<endl;
  report<<SSB_total<<endl;

  report<<"$SSB_region_overlap"<<endl;
  report<<SSB_region_overlap<<endl;
  report<<"$SSB_stock_overlap"<<endl;
  report<<SSB_stock_overlap<<endl;
  report<<"$SSB_natal_overlap"<<endl;
  report<<SSB_natal_overlap<<endl;

  report<<"$yield_region_overlap"<<endl;
  report<<yield_region_overlap<<endl;
  report<<"$yield_stock_overlap"<<endl;
  report<<yield_stock_overlap<<endl;
  report<<"$yield_natal_overlap"<<endl;
  report<<yield_natal_overlap<<endl;

  report<<"$biomass_AM_overlap_region"<<endl;
  report<<biomass_AM_overlap_region<<endl;
  report<<"$biomass_stock_overlap"<<endl;
  report<<biomass_stock_overlap<<endl;
  report<<"$biomass_natal_overlap"<<endl;
  report<<biomass_natal_overlap<<endl;

  report<<"$harvest_rate_region_bio_overlap"<<endl;
  report<<harvest_rate_region_bio_overlap<<endl;
  report<<"$harvest_rate_stock_bio_overlap"<<endl;
  report<<harvest_rate_stock_bio_overlap<<endl;
  report<<"$harvest_rate_natal_bio_overlap"<<endl;
  report<<harvest_rate_natal_bio_overlap<<endl;

  report<<"$depletion_region_overlap"<<endl;
  report<<depletion_region_overlap<<endl;
  report<<"$depletion_stock_overlap"<<endl;
  report<<depletion_stock_overlap<<endl;
  report<<"$depletion_natal_overlap"<<endl;
  report<<depletion_natal_overlap<<endl;

  report<<"$Bratio_stock"<<endl;
  report<<Bratio_stock<<endl;
  report<<"$Bratio_total"<<endl;
  report<<Bratio_total<<endl;

  report<<"$Bratio_stock_overlap"<<endl;
  report<<Bratio_stock_overlap<<endl;
  report<<"$Bratio_natal_overlap"<<endl;
  report<<Bratio_natal_overlap<<endl;


  
RUNTIME_SECTION
  convergence_criteria .001,.0001, 1.0e-4, 1.0e-7
  maximum_function_evaluations 100000
  


