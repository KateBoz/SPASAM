#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"
  #define EOUT(var) cout <<#var<<" "<<var<<endl;
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <tim_om.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nages.allocate("nages");
  nyrs.allocate("nyrs");
  npops.allocate("npops");
 int np=npops;
 int ny=nyrs;
 int na=nages;
  nregions.allocate(1,np,"nregions");
 ivector nreg=nregions;
  nfleets.allocate(1,np,"nfleets");
 ivector nf=nfleets;
  nfleets_survey.allocate(1,np,"nfleets_survey");
 ivector nfs=nfleets_survey;
  model_type_switch.allocate("model_type_switch");
  do_tag.allocate("do_tag");
  parse_TAC.allocate("parse_TAC");
  calc_TAC_from_uMSY.allocate("calc_TAC_from_uMSY");
  parse_TAC_source.allocate("parse_TAC_source");
  TAC_survey_parse_timelag_switch.allocate("TAC_survey_parse_timelag_switch");
  TAC_survey_parse_timelag.allocate("TAC_survey_parse_timelag");
  tsurvey.allocate(1,np,1,nreg,"tsurvey");
  larval_move_switch.allocate("larval_move_switch");
  move_switch.allocate("move_switch");
  use_input_Bstar.allocate("use_input_Bstar");
  natal_homing_switch.allocate("natal_homing_switch");
  spawn_return_switch.allocate("spawn_return_switch");
  select_switch.allocate("select_switch");
  select_switch_survey.allocate("select_switch_survey");
  F_switch.allocate("F_switch");
  recruit_devs_switch.allocate("recruit_devs_switch");
  recruit_randwalk_switch.allocate("recruit_randwalk_switch");
  maturity_switch_equil.allocate("maturity_switch_equil");
  SSB_type.allocate("SSB_type");
  apportionment_type.allocate("apportionment_type");
  Rec_type.allocate("Rec_type");
  use_stock_comp_info_survey.allocate("use_stock_comp_info_survey");
  use_stock_comp_info_catch.allocate("use_stock_comp_info_catch");
  nyrs_release.allocate("nyrs_release");
 int ny_rel=nyrs_release;
  yrs_releases.allocate(1,ny_rel,"yrs_releases");
  frac_total_abund_tagged.allocate(1,ny_rel,"frac_total_abund_tagged");
  max_life_tags.allocate("max_life_tags");
  SIM_ntag.allocate("SIM_ntag");
  report_rate.allocate(1,np,1,ny_rel,1,nreg,"report_rate");
  input_Bstar.allocate(1,np,1,nreg,"input_Bstar");
  SSB_zero_appor.allocate(1,np,1,nreg,"SSB_zero_appor");
  A.allocate(1,np,1,nreg,"A");
  return_age.allocate("return_age");
  return_probability.allocate(1,np,"return_probability");
  spawn_return_prob.allocate(1,np,"spawn_return_prob");
  phase_F.allocate("phase_F");
  phase_dummy.allocate("phase_dummy");
  tspawn.allocate(1,np,"tspawn");
  steep.allocate(1,np,"steep");
  ln_R_ave.allocate(1,np,"ln_R_ave");
  amplitude.allocate(1,np,"amplitude");
  freq.allocate(1,np,"freq");
  input_T.allocate(1,np,1,nreg,1,na,1,np,1,nreg,"input_T");
  input_residency_larval.allocate(1,np,1,nreg,"input_residency_larval");
  input_residency.allocate(1,np,1,nreg,1,na,"input_residency");
  sel_beta1.allocate(1,np,1,nreg,1,nf,"sel_beta1");
  sel_beta2.allocate(1,np,1,nreg,1,nf,"sel_beta2");
  sel_beta3.allocate(1,np,1,nreg,1,nf,"sel_beta3");
  sel_beta4.allocate(1,np,1,nreg,1,nf,"sel_beta4");
  sel_beta1_survey.allocate(1,np,1,nreg,1,nfs,"sel_beta1_survey");
  sel_beta2_survey.allocate(1,np,1,nreg,1,nfs,"sel_beta2_survey");
  sel_beta3_survey.allocate(1,np,1,nreg,1,nfs,"sel_beta3_survey");
  sel_beta4_survey.allocate(1,np,1,nreg,1,nfs,"sel_beta4_survey");
  input_selectivity.allocate(1,np,1,nreg,1,na,1,nf,"input_selectivity");
  input_survey_selectivity.allocate(1,np,1,nreg,1,na,1,nfs,"input_survey_selectivity");
  q_survey.allocate(1,np,1,nreg,1,nfs,"q_survey");
  input_F.allocate(1,np,1,nreg,1,nf,"input_F");
  dunce_F.allocate(1,np,1,nreg,1,3,"dunce_F");
  F_rho.allocate(1,np,1,nreg,1,nf,"F_rho");
  input_F_MSY.allocate("input_F_MSY");
  input_M.allocate(1,np,1,na,"input_M");
  sigma_recruit.allocate(1,np,"sigma_recruit");
  sigma_rec_prop.allocate(1,np,"sigma_rec_prop");
  sigma_F.allocate(1,np,1,nreg,1,nf,"sigma_F");
  input_weight.allocate(1,np,1,nreg,1,na,"input_weight");
  input_catch_weight.allocate(1,np,1,nreg,1,na,"input_catch_weight");
  fecundity.allocate(1,np,1,nreg,1,na,"fecundity");
  maturity.allocate(1,np,1,nreg,1,na,"maturity");
  prop_fem.allocate(1,np,1,nreg,"prop_fem");
  input_Rec_prop.allocate(1,np,1,nreg,"input_Rec_prop");
  equil_ssb_apport.allocate(1,np,1,nreg,"equil_ssb_apport");
  init_abund.allocate(1,np,1,np,1,nreg,1,na,"init_abund");
  rec_index_sigma.allocate(1,np,1,nreg,"rec_index_sigma");
  sigma_survey.allocate(1,np,1,nreg,1,nfs,"sigma_survey");
  sigma_survey_overlap.allocate(1,np,1,np,1,nreg,1,nfs,"sigma_survey_overlap");
  sigma_catch.allocate(1,np,1,nreg,1,nf,"sigma_catch");
  sigma_catch_overlap.allocate(1,np,1,np,1,nreg,1,nf,"sigma_catch_overlap");
  SIM_ncatch.allocate(1,np,1,nreg,1,nf,"SIM_ncatch");
  SIM_ncatch_overlap.allocate(1,np,1,np,1,nreg,1,nf,"SIM_ncatch_overlap");
  SIM_nsurvey.allocate(1,np,1,nreg,1,nf,"SIM_nsurvey");
  SIM_nsurvey_overlap.allocate(1,np,1,np,1,nreg,1,nf,"SIM_nsurvey_overlap");
  input_TAC.allocate(1,np,1,nreg,1,nf,"input_TAC");
  input_u.allocate(1,np,1,nreg,1,nf,"input_u");
  max_Fnew.allocate("max_Fnew");
  Fnew_start.allocate("Fnew_start");
  NR_iterationp.allocate("NR_iterationp");
  NR_dev.allocate("NR_dev");
  larval_move_switch_EM.allocate("larval_move_switch_EM");
  move_switch_EM.allocate("move_switch_EM");
  natal_homing_switch_EM.allocate("natal_homing_switch_EM");
  spawn_return_switch_EM.allocate("spawn_return_switch_EM");
  select_switch_EM.allocate("select_switch_EM");
  maturity_switch_equil_EM.allocate("maturity_switch_equil_EM");
  SSB_type_EM.allocate("SSB_type_EM");
  Rec_type_EM.allocate("Rec_type_EM");
  apportionment_type_EM.allocate("apportionment_type_EM");
  use_stock_comp_info_survey_EM.allocate("use_stock_comp_info_survey_EM");
  use_stock_comp_info_catch_EM.allocate("use_stock_comp_info_catch_EM");
  F_switch_EM.allocate("F_switch_EM");
  recruit_devs_switch_EM.allocate("recruit_devs_switch_EM");
  recruit_randwalk_switch_EM.allocate("recruit_randwalk_switch_EM");
  do_tag_EM.allocate("do_tag_EM");
  do_tag_mult.allocate("do_tag_mult");
  ph_lmr.allocate("ph_lmr");
  ph_rec.allocate("ph_rec");
  ph_abund_devs.allocate("ph_abund_devs");
  ph_F.allocate("ph_F");
  ph_steep.allocate("ph_steep");
  ph_M.allocate("ph_M");
  ph_sel_log.allocate("ph_sel_log");
  lb_sel_beta1.allocate("lb_sel_beta1");
  ub_sel_beta1.allocate("ub_sel_beta1");
  lb_sel_beta2.allocate("lb_sel_beta2");
  ub_sel_beta2.allocate("ub_sel_beta2");
  ph_sel_log_surv.allocate("ph_sel_log_surv");
  ph_sel_dubl.allocate("ph_sel_dubl");
  ph_sel_dubl_surv.allocate("ph_sel_dubl_surv");
  ph_q.allocate("ph_q");
  ph_F_rho.allocate("ph_F_rho");
  ph_T_YR.allocate("ph_T_YR");
  ph_T_CNST.allocate("ph_T_CNST");
  ph_dummy.allocate("ph_dummy");
  wt_surv.allocate("wt_surv");
  wt_catch.allocate("wt_catch");
  wt_fish_age.allocate("wt_fish_age");
  wt_srv_age.allocate("wt_srv_age");
  wt_rec.allocate("wt_rec");
  wt_tag.allocate("wt_tag");
  abund_pen_switch.allocate("abund_pen_switch");
  move_pen_switch.allocate("move_pen_switch");
  Tpen.allocate("Tpen");
  Tpen2.allocate("Tpen2");
  OBS_survey_fleet_bio_se.allocate(1,np,1,nreg,1,ny,1,nfs,"OBS_survey_fleet_bio_se");
  OBS_yield_fleet_se.allocate(1,np,1,nreg,1,ny,1,nf,"OBS_yield_fleet_se");
  OBS_survey_prop_N.allocate(1,np,1,nreg,1,ny,1,nfs,"OBS_survey_prop_N");
  OBS_catch_prop_N.allocate(1,np,1,nreg,1,ny,1,nf,"OBS_catch_prop_N");
  tag_N.allocate(1,np,1,nreg,1,ny_rel,1,na,"tag_N");
  debug.allocate("debug");
  myseed_yield.allocate("myseed_yield");
  myseed_survey.allocate("myseed_survey");
  myseed_F.allocate("myseed_F");
  myseed_rec_devs.allocate("myseed_rec_devs");
  myseed_rec_apport.allocate("myseed_rec_apport");
  myseed_rec_index.allocate("myseed_rec_index");
  myseed_survey_age.allocate("myseed_survey_age");
  myseed_catch_age.allocate("myseed_catch_age");
  myseed_tag.allocate("myseed_tag");
  years.allocate(1,nyrs);
years.fill_seqadd(double(1),1.0);
  nregions_temp.allocate(1,np,1,np,"nregions_temp");
 for(int j=1;j<=npops;j++) //recap stock
 {
  for (int r=1;r<=npops;r++) //recap region
  {
    if(j<=r)
     {
     nregions_temp(j,r)=0;
     }
    if(j>r)
     {
     nregions_temp(j,r)=nreg(r); //create temp matrix that holds the number of regions that exist in all previous populations (so can sum for use in calcs below)
     }
   }
  }
  nreg_temp.allocate(1,np);
 cout << "debug = " << debug << endl;
 cout << "If debug != 1541 then .dat file not setup correctly" << endl;
 cout << "input read" << endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
 ivector nr=nregions;
 int nps=npops;
 int nyr=nyrs;
 int nag=nages;
 ivector nfl=nfleets;
 ivector nfls=nfleets_survey;  
  F_est.allocate(1,nps,1,nr,phase_F,"F_est");
  Fstartyr.allocate(1,nps,1,nr,"Fstartyr");
  #ifndef NO_AD_INITIALIZE
    Fstartyr.initialize();
  #endif
  minF.allocate(1,nps,1,nr,"minF");
  #ifndef NO_AD_INITIALIZE
    minF.initialize();
  #endif
  maxF.allocate(1,nps,1,nr,"maxF");
  #ifndef NO_AD_INITIALIZE
    maxF.initialize();
  #endif
  stepF.allocate(1,nps,1,nr,"stepF");
  #ifndef NO_AD_INITIALIZE
    stepF.initialize();
  #endif
  R_ave.allocate(1,nps,"R_ave");
  #ifndef NO_AD_INITIALIZE
    R_ave.initialize();
  #endif
  T.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr,"T");
  #ifndef NO_AD_INITIALIZE
    T.initialize();
  #endif
  T_year.allocate(1,nps,1,nr,1,nyr,1,nps,1,nr,"T_year");
  #ifndef NO_AD_INITIALIZE
    T_year.initialize();
  #endif
  rel_bio.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr,"rel_bio");
  #ifndef NO_AD_INITIALIZE
    rel_bio.initialize();
  #endif
  Bstar.allocate(1,nps,1,nr,"Bstar");
  #ifndef NO_AD_INITIALIZE
    Bstar.initialize();
  #endif
  c.allocate(1,nps,1,nr,"c");
  #ifndef NO_AD_INITIALIZE
    c.initialize();
  #endif
  Fract_Move_DD.allocate(1,nps,1,nr,1,nyr,1,nag,"Fract_Move_DD");
  #ifndef NO_AD_INITIALIZE
    Fract_Move_DD.initialize();
  #endif
  selectivity.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfl,"selectivity");
  #ifndef NO_AD_INITIALIZE
    selectivity.initialize();
  #endif
  selectivity_age.allocate(1,nps,1,nr,1,nag,1,nfl,"selectivity_age");
  #ifndef NO_AD_INITIALIZE
    selectivity_age.initialize();
  #endif
  F_year.allocate(1,nps,1,nr,1,nyr,1,nfl,"F_year");
  #ifndef NO_AD_INITIALIZE
    F_year.initialize();
  #endif
  F_fleet.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfl,"F_fleet");
  #ifndef NO_AD_INITIALIZE
    F_fleet.initialize();
  #endif
  F.allocate(1,nps,1,nr,1,nyr,1,nag,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  M.allocate(1,nps,1,nr,1,nyr,1,nag,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  rec_devs.allocate(1,nps,1,nyr,"rec_devs");
  #ifndef NO_AD_INITIALIZE
    rec_devs.initialize();
  #endif
  rec_devs_randwalk.allocate(1,nps,1,nyr,"rec_devs_randwalk");
  #ifndef NO_AD_INITIALIZE
    rec_devs_randwalk.initialize();
  #endif
  weight_population.allocate(1,nps,1,nr,1,nyr,1,nag,"weight_population");
  #ifndef NO_AD_INITIALIZE
    weight_population.initialize();
  #endif
  weight_catch.allocate(1,nps,1,nr,1,nyr,1,nag,"weight_catch");
  #ifndef NO_AD_INITIALIZE
    weight_catch.initialize();
  #endif
  wt_mat_mult.allocate(1,nps,1,nyr,1,nag,"wt_mat_mult");
  #ifndef NO_AD_INITIALIZE
    wt_mat_mult.initialize();
  #endif
  wt_mat_mult_reg.allocate(1,nps,1,nr,1,nyr,1,nag,"wt_mat_mult_reg");
  #ifndef NO_AD_INITIALIZE
    wt_mat_mult_reg.initialize();
  #endif
  ave_mat_temp.allocate(1,nps,1,nag,1,nr,"ave_mat_temp");
  #ifndef NO_AD_INITIALIZE
    ave_mat_temp.initialize();
  #endif
  ave_mat.allocate(1,nps,1,nag,"ave_mat");
  #ifndef NO_AD_INITIALIZE
    ave_mat.initialize();
  #endif
  SPR_N.allocate(1,nps,1,nag,"SPR_N");
  #ifndef NO_AD_INITIALIZE
    SPR_N.initialize();
  #endif
  SPR_SSB.allocate(1,nps,1,nag,"SPR_SSB");
  #ifndef NO_AD_INITIALIZE
    SPR_SSB.initialize();
  #endif
  SPR.allocate(1,nps,"SPR");
  #ifndef NO_AD_INITIALIZE
    SPR.initialize();
  #endif
  SSB_zero.allocate(1,nps,"SSB_zero");
  #ifndef NO_AD_INITIALIZE
    SSB_zero.initialize();
  #endif
  alpha.allocate(1,nps,"alpha");
  #ifndef NO_AD_INITIALIZE
    alpha.initialize();
  #endif
  beta.allocate(1,nps,"beta");
  #ifndef NO_AD_INITIALIZE
    beta.initialize();
  #endif
  recruits_BM.allocate(1,nps,1,nr,1,nyr,"recruits_BM");
  #ifndef NO_AD_INITIALIZE
    recruits_BM.initialize();
  #endif
  recruits_AM.allocate(1,nps,1,nr,1,nyr,"recruits_AM");
  #ifndef NO_AD_INITIALIZE
    recruits_AM.initialize();
  #endif
  Rec_Prop.allocate(1,nps,1,nr,1,nyr,"Rec_Prop");
  #ifndef NO_AD_INITIALIZE
    Rec_Prop.initialize();
  #endif
  Rec_prop_temp1.allocate(1,nps,1,nyr,1,nr,"Rec_prop_temp1");
  #ifndef NO_AD_INITIALIZE
    Rec_prop_temp1.initialize();
  #endif
  Rec_prop_temp2.allocate(1,nps,1,nyr,1,nr,"Rec_prop_temp2");
  #ifndef NO_AD_INITIALIZE
    Rec_prop_temp2.initialize();
  #endif
  rec_index_BM.allocate(1,nps,1,nr,1,nyr,"rec_index_BM");
  #ifndef NO_AD_INITIALIZE
    rec_index_BM.initialize();
  #endif
  rec_index_prop_BM.allocate(1,nps,1,nr,1,nyr,"rec_index_prop_BM");
  #ifndef NO_AD_INITIALIZE
    rec_index_prop_BM.initialize();
  #endif
  rec_index_BM_temp.allocate(1,nps,1,nyr,1,nr,"rec_index_BM_temp");
  #ifndef NO_AD_INITIALIZE
    rec_index_BM_temp.initialize();
  #endif
  rec_index_AM.allocate(1,nps,1,nr,1,nyr,"rec_index_AM");
  #ifndef NO_AD_INITIALIZE
    rec_index_AM.initialize();
  #endif
  rec_index_prop_AM.allocate(1,nps,1,nr,1,nyr,"rec_index_prop_AM");
  #ifndef NO_AD_INITIALIZE
    rec_index_prop_AM.initialize();
  #endif
  rec_index_AM_temp.allocate(1,nps,1,nyr,1,nr,"rec_index_AM_temp");
  #ifndef NO_AD_INITIALIZE
    rec_index_AM_temp.initialize();
  #endif
  env_rec.allocate(1,nyr,"env_rec");
  #ifndef NO_AD_INITIALIZE
    env_rec.initialize();
  #endif
  abundance_at_age_BM.allocate(1,nps,1,nr,1,nyr,1,nag,"abundance_at_age_BM");
  #ifndef NO_AD_INITIALIZE
    abundance_at_age_BM.initialize();
  #endif
  abundance_at_age_AM.allocate(1,nps,1,nr,1,nyr,1,nag,"abundance_at_age_AM");
  #ifndef NO_AD_INITIALIZE
    abundance_at_age_AM.initialize();
  #endif
  abundance_in.allocate(1,nps,1,nr,1,nyr,1,nag,"abundance_in");
  #ifndef NO_AD_INITIALIZE
    abundance_in.initialize();
  #endif
  abundance_res.allocate(1,nps,1,nr,1,nyr,1,nag,"abundance_res");
  #ifndef NO_AD_INITIALIZE
    abundance_res.initialize();
  #endif
  abundance_leave.allocate(1,nps,1,nr,1,nyr,1,nag,"abundance_leave");
  #ifndef NO_AD_INITIALIZE
    abundance_leave.initialize();
  #endif
  abundance_spawn.allocate(1,nps,1,nr,1,nyr,1,nag,"abundance_spawn");
  #ifndef NO_AD_INITIALIZE
    abundance_spawn.initialize();
  #endif
  biomass_BM_age.allocate(1,nps,1,nr,1,nyr,1,nag,"biomass_BM_age");
  #ifndef NO_AD_INITIALIZE
    biomass_BM_age.initialize();
  #endif
  biomass_AM_age.allocate(1,nps,1,nr,1,nyr,1,nag,"biomass_AM_age");
  #ifndef NO_AD_INITIALIZE
    biomass_AM_age.initialize();
  #endif
  biomass_BM.allocate(1,nps,1,nr,1,nyr,"biomass_BM");
  #ifndef NO_AD_INITIALIZE
    biomass_BM.initialize();
  #endif
  biomass_AM.allocate(1,nps,1,nr,1,nyr,"biomass_AM");
  #ifndef NO_AD_INITIALIZE
    biomass_AM.initialize();
  #endif
  bio_in.allocate(1,nps,1,nr,1,nyr,1,nag,"bio_in");
  #ifndef NO_AD_INITIALIZE
    bio_in.initialize();
  #endif
  bio_res.allocate(1,nps,1,nr,1,nyr,1,nag,"bio_res");
  #ifndef NO_AD_INITIALIZE
    bio_res.initialize();
  #endif
  bio_leave.allocate(1,nps,1,nr,1,nyr,1,nag,"bio_leave");
  #ifndef NO_AD_INITIALIZE
    bio_leave.initialize();
  #endif
 int nyr_rel=nyrs_release;
 ivector xy(1,nyr_rel);
 ivector nt(1,nyr_rel);
 ivector nt2(1,nyr_rel);
 int nt3=max_life_tags*sum(nregions)+1;
 ivector tag_age(1,nyr_rel);
  for(int x=1; x<=nyrs_release; x++)
   {
    xx=yrs_releases(x);
    xy(x)=min(max_life_tags,nyrs-xx+1);
    nt(x)=xy(x)*sum(nregions)+1;
    nt2(x)=nt(x)-1;
    tag_age(x)=xy(x);
   }
  tags_avail.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr,"tags_avail");
  #ifndef NO_AD_INITIALIZE
    tags_avail.initialize();
  #endif
  recaps.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr,"recaps");
  #ifndef NO_AD_INITIALIZE
    recaps.initialize();
  #endif
  tag_prop.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr,"tag_prop");
  #ifndef NO_AD_INITIALIZE
    tag_prop.initialize();
  #endif
  tag_prop_final.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,nt,"tag_prop_final");
  #ifndef NO_AD_INITIALIZE
    tag_prop_final.initialize();
  #endif
  SIM_tag_prop.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,nt,"SIM_tag_prop");
  #ifndef NO_AD_INITIALIZE
    SIM_tag_prop.initialize();
  #endif
  OBS_tag_prop_final.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,nt,"OBS_tag_prop_final");
  #ifndef NO_AD_INITIALIZE
    OBS_tag_prop_final.initialize();
  #endif
  total_recap_temp.allocate(1,nps,1,nr,1,tag_age,"total_recap_temp");
  #ifndef NO_AD_INITIALIZE
    total_recap_temp.initialize();
  #endif
  rand_tag_prop_temp2.allocate(1,nt3,"rand_tag_prop_temp2");
  #ifndef NO_AD_INITIALIZE
    rand_tag_prop_temp2.initialize();
  #endif
  tag_prop_temp2.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,nt2,"tag_prop_temp2");
  #ifndef NO_AD_INITIALIZE
    tag_prop_temp2.initialize();
  #endif
  rand_tag_prop_temp.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,2000,"rand_tag_prop_temp");
  #ifndef NO_AD_INITIALIZE
    rand_tag_prop_temp.initialize();
  #endif
  tags_avail_temp.allocate(1,nps,1,nr,"tags_avail_temp");
  #ifndef NO_AD_INITIALIZE
    tags_avail_temp.initialize();
  #endif
  tag_prop_temp.allocate(1,nps,1,nyr_rel,1,nr,"tag_prop_temp");
  #ifndef NO_AD_INITIALIZE
    tag_prop_temp.initialize();
  #endif
  ntags_total.allocate(1,nyr_rel,"ntags_total");
  #ifndef NO_AD_INITIALIZE
    ntags_total.initialize();
  #endif
  ntags.allocate(1,nps,1,nr,1,nyr_rel,1,nag,"ntags");
  #ifndef NO_AD_INITIALIZE
    ntags.initialize();
  #endif
  total_rec.allocate(1,nps,1,nr,1,nyr_rel,1,nag,"total_rec");
  #ifndef NO_AD_INITIALIZE
    total_rec.initialize();
  #endif
  not_rec.allocate(1,nps,1,nr,1,nyr_rel,1,nag,"not_rec");
  #ifndef NO_AD_INITIALIZE
    not_rec.initialize();
  #endif
  tag_prop_not_rec.allocate(1,nps,1,nr,1,nyr_rel,1,nag,"tag_prop_not_rec");
  #ifndef NO_AD_INITIALIZE
    tag_prop_not_rec.initialize();
  #endif
  survey_selectivity.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfls,"survey_selectivity");
  #ifndef NO_AD_INITIALIZE
    survey_selectivity.initialize();
  #endif
  survey_selectivity_age.allocate(1,nps,1,nr,1,nag,1,nfls,"survey_selectivity_age");
  #ifndef NO_AD_INITIALIZE
    survey_selectivity_age.initialize();
  #endif
  survey_selectivity_temp.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nag,"survey_selectivity_temp");
  #ifndef NO_AD_INITIALIZE
    survey_selectivity_temp.initialize();
  #endif
  true_survey_fleet_overlap_age.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag,"true_survey_fleet_overlap_age");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_overlap_age.initialize();
  #endif
  survey_at_age_region_fleet_overlap_prop.allocate(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag,"survey_at_age_region_fleet_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    survey_at_age_region_fleet_overlap_prop.initialize();
  #endif
  SIM_survey_prop_overlap.allocate(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag,"SIM_survey_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    SIM_survey_prop_overlap.initialize();
  #endif
  OBS_survey_prop_overlap.allocate(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag,"OBS_survey_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_prop_overlap.initialize();
  #endif
  true_survey_fleet_overlap_age_bio.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag,"true_survey_fleet_overlap_age_bio");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_overlap_age_bio.initialize();
  #endif
  true_survey_fleet_bio_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,"true_survey_fleet_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_bio_overlap.initialize();
  #endif
  true_survey_region_bio_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"true_survey_region_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    true_survey_region_bio_overlap.initialize();
  #endif
  true_survey_population_bio_overlap.allocate(1,nps,1,nyr,1,nps,"true_survey_population_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    true_survey_population_bio_overlap.initialize();
  #endif
  true_survey_natal_bio_overlap.allocate(1,nyr,1,nps,"true_survey_natal_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    true_survey_natal_bio_overlap.initialize();
  #endif
  true_survey_total_bio_overlap.allocate(1,nyr,"true_survey_total_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    true_survey_total_bio_overlap.initialize();
  #endif
  true_survey_fleet_age.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nag,"true_survey_fleet_age");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_age.initialize();
  #endif
  true_survey_fleet_age_temp.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfls,"true_survey_fleet_age_temp");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_age_temp.initialize();
  #endif
  true_survey_region_abundance.allocate(1,nps,1,nyr,1,nr,1,nag,"true_survey_region_abundance");
  #ifndef NO_AD_INITIALIZE
    true_survey_region_abundance.initialize();
  #endif
  true_survey_region_abundance_temp.allocate(1,nps,1,nyr,1,nag,1,nr,"true_survey_region_abundance_temp");
  #ifndef NO_AD_INITIALIZE
    true_survey_region_abundance_temp.initialize();
  #endif
  true_survey_population_abundance.allocate(1,nyr,1,nps,1,nag,"true_survey_population_abundance");
  #ifndef NO_AD_INITIALIZE
    true_survey_population_abundance.initialize();
  #endif
  true_survey_population_abundance_temp.allocate(1,nyr,1,nag,1,nps,"true_survey_population_abundance_temp");
  #ifndef NO_AD_INITIALIZE
    true_survey_population_abundance_temp.initialize();
  #endif
  true_survey_total_abundance.allocate(1,nyr,1,nag,"true_survey_total_abundance");
  #ifndef NO_AD_INITIALIZE
    true_survey_total_abundance.initialize();
  #endif
  survey_at_age_fleet_prop.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nag,"survey_at_age_fleet_prop");
  #ifndef NO_AD_INITIALIZE
    survey_at_age_fleet_prop.initialize();
  #endif
  SIM_survey_prop.allocate(1,nps,1,nr,1,nfls,1,nyr,1,nag,"SIM_survey_prop");
  #ifndef NO_AD_INITIALIZE
    SIM_survey_prop.initialize();
  #endif
  OBS_survey_prop.allocate(1,nps,1,nr,1,nfls,1,nyr,1,nag,"OBS_survey_prop");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_prop.initialize();
  #endif
  true_survey_fleet_age_bio.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nag,"true_survey_fleet_age_bio");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_age_bio.initialize();
  #endif
  true_survey_fleet_bio.allocate(1,nps,1,nr,1,nyr,1,nfls,"true_survey_fleet_bio");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_bio.initialize();
  #endif
  true_survey_region_bio.allocate(1,nps,1,nyr,1,nr,"true_survey_region_bio");
  #ifndef NO_AD_INITIALIZE
    true_survey_region_bio.initialize();
  #endif
  true_survey_population_bio.allocate(1,nyr,1,nps,"true_survey_population_bio");
  #ifndef NO_AD_INITIALIZE
    true_survey_population_bio.initialize();
  #endif
  true_survey_total_bio.allocate(1,nyr,"true_survey_total_bio");
  #ifndef NO_AD_INITIALIZE
    true_survey_total_bio.initialize();
  #endif
  OBS_survey_fleet_bio_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,"OBS_survey_fleet_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_fleet_bio_overlap.initialize();
  #endif
  OBS_survey_region_bio_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"OBS_survey_region_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_region_bio_overlap.initialize();
  #endif
  OBS_survey_population_bio_overlap.allocate(1,nps,1,nyr,1,nps,"OBS_survey_population_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_population_bio_overlap.initialize();
  #endif
  OBS_survey_natal_bio_overlap.allocate(1,nyr,1,nps,"OBS_survey_natal_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_natal_bio_overlap.initialize();
  #endif
  OBS_survey_total_bio_overlap.allocate(1,nyr,"OBS_survey_total_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_total_bio_overlap.initialize();
  #endif
  OBS_survey_fleet_bio.allocate(1,nps,1,nr,1,nyr,1,nfls,"OBS_survey_fleet_bio");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_fleet_bio.initialize();
  #endif
  OBS_survey_region_bio.allocate(1,nps,1,nyr,1,nr,"OBS_survey_region_bio");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_region_bio.initialize();
  #endif
  OBS_survey_population_bio.allocate(1,nyr,1,nps,"OBS_survey_population_bio");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_population_bio.initialize();
  #endif
  OBS_survey_total_bio.allocate(1,nyr,"OBS_survey_total_bio");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_total_bio.initialize();
  #endif
  apport_region_survey_biomass.allocate(1,nps,1,nr,1,nyr,"apport_region_survey_biomass");
  #ifndef NO_AD_INITIALIZE
    apport_region_survey_biomass.initialize();
  #endif
  catch_at_age_fleet.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfl,"catch_at_age_fleet");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_fleet.initialize();
  #endif
  catch_at_age_fleet_prop.allocate(1,nps,1,nr,1,nyr,1,nfl,1,nag,"catch_at_age_fleet_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_fleet_prop.initialize();
  #endif
  SIM_catch_prop.allocate(1,nps,1,nr,1,nfl,1,nyr,1,nag,"SIM_catch_prop");
  #ifndef NO_AD_INITIALIZE
    SIM_catch_prop.initialize();
  #endif
  OBS_catch_prop.allocate(1,nps,1,nr,1,nfl,1,nyr,1,nag,"OBS_catch_prop");
  #ifndef NO_AD_INITIALIZE
    OBS_catch_prop.initialize();
  #endif
  yield_fleet.allocate(1,nps,1,nr,1,nyr,1,nfl,"yield_fleet");
  #ifndef NO_AD_INITIALIZE
    yield_fleet.initialize();
  #endif
  catch_at_age_region.allocate(1,nps,1,nr,1,nyr,1,nag,"catch_at_age_region");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region.initialize();
  #endif
  catch_at_age_region_prop.allocate(1,nps,1,nr,1,nyr,1,nag,"catch_at_age_region_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region_prop.initialize();
  #endif
  yield_region.allocate(1,nps,1,nr,1,nyr,"yield_region");
  #ifndef NO_AD_INITIALIZE
    yield_region.initialize();
  #endif
  catch_at_age_population.allocate(1,nps,1,nyr,1,nag,"catch_at_age_population");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population.initialize();
  #endif
  catch_at_age_population_prop.allocate(1,nps,1,nyr,1,nag,"catch_at_age_population_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population_prop.initialize();
  #endif
  yield_population.allocate(1,nps,1,nyr,"yield_population");
  #ifndef NO_AD_INITIALIZE
    yield_population.initialize();
  #endif
  SSB_region.allocate(1,nps,1,nr,1,nyr,"SSB_region");
  #ifndef NO_AD_INITIALIZE
    SSB_region.initialize();
  #endif
  SSB_population.allocate(1,nps,1,nyr,"SSB_population");
  #ifndef NO_AD_INITIALIZE
    SSB_population.initialize();
  #endif
  SSB_total.allocate(1,nyr,"SSB_total");
  #ifndef NO_AD_INITIALIZE
    SSB_total.initialize();
  #endif
  abundance_population.allocate(1,nps,1,nyr,1,nag,"abundance_population");
  #ifndef NO_AD_INITIALIZE
    abundance_population.initialize();
  #endif
  abundance_total.allocate(1,nyr,1,nag,"abundance_total");
  #ifndef NO_AD_INITIALIZE
    abundance_total.initialize();
  #endif
  biomass_population.allocate(1,nps,1,nyr,"biomass_population");
  #ifndef NO_AD_INITIALIZE
    biomass_population.initialize();
  #endif
  biomass_total.allocate(1,nyr,"biomass_total");
  #ifndef NO_AD_INITIALIZE
    biomass_total.initialize();
  #endif
  catch_at_age_total.allocate(1,nyr,1,nag,"catch_at_age_total");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_total.initialize();
  #endif
  catch_at_age_total_prop.allocate(1,nyr,1,nag,"catch_at_age_total_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_total_prop.initialize();
  #endif
  yield_total.allocate(1,nyr,"yield_total");
  #ifndef NO_AD_INITIALIZE
    yield_total.initialize();
  #endif
  harvest_rate_region_num.allocate(1,nps,1,nr,1,nyr,1,nag,"harvest_rate_region_num");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_region_num.initialize();
  #endif
  harvest_rate_population_num.allocate(1,nps,1,nyr,1,nag,"harvest_rate_population_num");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_population_num.initialize();
  #endif
  harvest_rate_total_num.allocate(1,nyr,1,nag,"harvest_rate_total_num");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_total_num.initialize();
  #endif
  harvest_rate_region_bio.allocate(1,nps,1,nr,1,nyr,"harvest_rate_region_bio");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_region_bio.initialize();
  #endif
  harvest_rate_population_bio.allocate(1,nps,1,nyr,"harvest_rate_population_bio");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_population_bio.initialize();
  #endif
  harvest_rate_total_bio.allocate(1,nyr,"harvest_rate_total_bio");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_total_bio.initialize();
  #endif
  depletion_region.allocate(1,nps,1,nr,1,nyr,"depletion_region");
  #ifndef NO_AD_INITIALIZE
    depletion_region.initialize();
  #endif
  depletion_population.allocate(1,nps,1,nyr,"depletion_population");
  #ifndef NO_AD_INITIALIZE
    depletion_population.initialize();
  #endif
  depletion_total.allocate(1,nyr,"depletion_total");
  #ifndef NO_AD_INITIALIZE
    depletion_total.initialize();
  #endif
  abundance_at_age_BM_overlap_region.allocate(1,nps,1,nps,1,nyr,1,nag,1,nr,"abundance_at_age_BM_overlap_region");
  #ifndef NO_AD_INITIALIZE
    abundance_at_age_BM_overlap_region.initialize();
  #endif
  abundance_at_age_BM_overlap_population.allocate(1,nps,1,nps,1,nyr,1,nag,"abundance_at_age_BM_overlap_population");
  #ifndef NO_AD_INITIALIZE
    abundance_at_age_BM_overlap_population.initialize();
  #endif
  abundance_at_age_AM_overlap_region.allocate(1,nps,1,nps,1,nyr,1,nag,1,nr,"abundance_at_age_AM_overlap_region");
  #ifndef NO_AD_INITIALIZE
    abundance_at_age_AM_overlap_region.initialize();
  #endif
  abundance_at_age_AM_overlap_population.allocate(1,nps,1,nps,1,nyr,1,nag,"abundance_at_age_AM_overlap_population");
  #ifndef NO_AD_INITIALIZE
    abundance_at_age_AM_overlap_population.initialize();
  #endif
  abundance_AM_overlap_region_all_natal.allocate(1,nps,1,nr,1,nyr,1,nag,"abundance_AM_overlap_region_all_natal");
  #ifndef NO_AD_INITIALIZE
    abundance_AM_overlap_region_all_natal.initialize();
  #endif
  abundance_spawn_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"abundance_spawn_overlap");
  #ifndef NO_AD_INITIALIZE
    abundance_spawn_overlap.initialize();
  #endif
  SSB_region_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,"SSB_region_overlap");
  #ifndef NO_AD_INITIALIZE
    SSB_region_overlap.initialize();
  #endif
  SSB_population_overlap.allocate(1,nps,1,nps,1,nyr,"SSB_population_overlap");
  #ifndef NO_AD_INITIALIZE
    SSB_population_overlap.initialize();
  #endif
  SSB_natal_overlap.allocate(1,nps,1,nyr,"SSB_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    SSB_natal_overlap.initialize();
  #endif
  catch_at_age_region_fleet_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"catch_at_age_region_fleet_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region_fleet_overlap.initialize();
  #endif
  catch_at_age_region_fleet_overlap_prop.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"catch_at_age_region_fleet_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region_fleet_overlap_prop.initialize();
  #endif
  SIM_catch_prop_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"SIM_catch_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    SIM_catch_prop_overlap.initialize();
  #endif
  OBS_catch_prop_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"OBS_catch_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_catch_prop_overlap.initialize();
  #endif
  catch_at_age_region_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"catch_at_age_region_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region_overlap.initialize();
  #endif
  catch_at_age_region_overlap_prop.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"catch_at_age_region_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region_overlap_prop.initialize();
  #endif
  yield_region_fleet_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,"yield_region_fleet_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_region_fleet_overlap.initialize();
  #endif
  yield_region_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,"yield_region_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_region_overlap.initialize();
  #endif
  catch_at_age_population_overlap.allocate(1,nps,1,nps,1,nyr,1,nag,"catch_at_age_population_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population_overlap.initialize();
  #endif
  catch_at_age_population_overlap_prop.allocate(1,nps,1,nps,1,nyr,1,nag,"catch_at_age_population_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population_overlap_prop.initialize();
  #endif
  yield_population_overlap.allocate(1,nps,1,nps,1,nyr,"yield_population_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_population_overlap.initialize();
  #endif
  abundance_natal_overlap.allocate(1,nps,1,nyr,1,nag,"abundance_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    abundance_natal_overlap.initialize();
  #endif
  biomass_BM_age_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"biomass_BM_age_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_BM_age_overlap.initialize();
  #endif
  biomass_AM_age_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"biomass_AM_age_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_AM_age_overlap.initialize();
  #endif
  biomass_BM_overlap_region.allocate(1,nps,1,nps,1,nr,1,nyr,"biomass_BM_overlap_region");
  #ifndef NO_AD_INITIALIZE
    biomass_BM_overlap_region.initialize();
  #endif
  biomass_AM_overlap_region.allocate(1,nps,1,nps,1,nr,1,nyr,"biomass_AM_overlap_region");
  #ifndef NO_AD_INITIALIZE
    biomass_AM_overlap_region.initialize();
  #endif
  biomass_AM_overlap_region_all_natal_temp.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,"biomass_AM_overlap_region_all_natal_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_AM_overlap_region_all_natal_temp.initialize();
  #endif
  biomass_AM_overlap_age_region_all_natal.allocate(1,nps,1,nr,1,nyr,1,nag,"biomass_AM_overlap_age_region_all_natal");
  #ifndef NO_AD_INITIALIZE
    biomass_AM_overlap_age_region_all_natal.initialize();
  #endif
  biomass_AM_overlap_region_all_natal.allocate(1,nps,1,nr,1,nyr,"biomass_AM_overlap_region_all_natal");
  #ifndef NO_AD_INITIALIZE
    biomass_AM_overlap_region_all_natal.initialize();
  #endif
  biomass_population_overlap.allocate(1,nps,1,nps,1,nyr,"biomass_population_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_population_overlap.initialize();
  #endif
  biomass_natal_overlap.allocate(1,nps,1,nyr,"biomass_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_natal_overlap.initialize();
  #endif
  catch_at_age_natal_overlap.allocate(1,nps,1,nyr,1,nag,"catch_at_age_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_natal_overlap.initialize();
  #endif
  catch_at_age_natal_overlap_prop.allocate(1,nps,1,nyr,1,nag,"catch_at_age_natal_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_natal_overlap_prop.initialize();
  #endif
  yield_natal_overlap.allocate(1,nps,1,nyr,"yield_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_natal_overlap.initialize();
  #endif
  harvest_rate_region_fleet_bio_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,"harvest_rate_region_fleet_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_region_fleet_bio_overlap.initialize();
  #endif
  harvest_rate_region_bio_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,"harvest_rate_region_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_region_bio_overlap.initialize();
  #endif
  harvest_rate_population_bio_overlap.allocate(1,nps,1,nps,1,nyr,"harvest_rate_population_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_population_bio_overlap.initialize();
  #endif
  harvest_rate_natal_bio_overlap.allocate(1,nps,1,nyr,"harvest_rate_natal_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    harvest_rate_natal_bio_overlap.initialize();
  #endif
  depletion_region_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,"depletion_region_overlap");
  #ifndef NO_AD_INITIALIZE
    depletion_region_overlap.initialize();
  #endif
  depletion_population_overlap.allocate(1,nps,1,nps,1,nyr,"depletion_population_overlap");
  #ifndef NO_AD_INITIALIZE
    depletion_population_overlap.initialize();
  #endif
  depletion_natal_overlap.allocate(1,nps,1,nyr,"depletion_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    depletion_natal_overlap.initialize();
  #endif
  Bratio_population_overlap.allocate(1,nps,1,nps,1,nyr,"Bratio_population_overlap");
  #ifndef NO_AD_INITIALIZE
    Bratio_population_overlap.initialize();
  #endif
  Bratio_natal_overlap.allocate(1,nps,1,nyr,"Bratio_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    Bratio_natal_overlap.initialize();
  #endif
  Bratio_population.allocate(1,nps,1,nyr,"Bratio_population");
  #ifndef NO_AD_INITIALIZE
    Bratio_population.initialize();
  #endif
  Bratio_total.allocate(1,nyr,"Bratio_total");
  #ifndef NO_AD_INITIALIZE
    Bratio_total.initialize();
  #endif
  OBS_yield_region_fleet_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfl,"OBS_yield_region_fleet_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_region_fleet_overlap.initialize();
  #endif
  OBS_yield_region_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"OBS_yield_region_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_region_overlap.initialize();
  #endif
  OBS_yield_population_overlap.allocate(1,nps,1,nyr,1,nps,"OBS_yield_population_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_population_overlap.initialize();
  #endif
  OBS_yield_natal_overlap.allocate(1,nyr,1,nps,"OBS_yield_natal_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_natal_overlap.initialize();
  #endif
  OBS_yield_total_overlap.allocate(1,nyr,"OBS_yield_total_overlap");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_total_overlap.initialize();
  #endif
  OBS_yield_fleet.allocate(1,nps,1,nr,1,nyr,1,nfl,"OBS_yield_fleet");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_fleet.initialize();
  #endif
  OBS_yield_region.allocate(1,nps,1,nyr,1,nr,"OBS_yield_region");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_region.initialize();
  #endif
  OBS_yield_population.allocate(1,nyr,1,nps,"OBS_yield_population");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_population.initialize();
  #endif
  OBS_yield_total.allocate(1,nyr,"OBS_yield_total");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_total.initialize();
  #endif
  apport_yield_region.allocate(1,nps,1,nr,1,nyr,"apport_yield_region");
  #ifndef NO_AD_INITIALIZE
    apport_yield_region.initialize();
  #endif
  biomass_BM_temp.allocate(1,nps,1,nr,"biomass_BM_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_BM_temp.initialize();
  #endif
  biomass_BM_temp2.allocate("biomass_BM_temp2");
  #ifndef NO_AD_INITIALIZE
  biomass_BM_temp2.initialize();
  #endif
  biomass_BM_overlap_temp.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,"biomass_BM_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_BM_overlap_temp.initialize();
  #endif
  init_abund_temp.allocate(1,nps,1,nr,1,nag,1,nps,"init_abund_temp");
  #ifndef NO_AD_INITIALIZE
    init_abund_temp.initialize();
  #endif
  rand_SIM_survey_prop_temp.allocate(1,nps,1,nr,1,nyr,1,nfls,1,2000,"rand_SIM_survey_prop_temp");
  #ifndef NO_AD_INITIALIZE
    rand_SIM_survey_prop_temp.initialize();
  #endif
  rand_SIM_survey_prop_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,2000,"rand_SIM_survey_prop_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    rand_SIM_survey_prop_temp_overlap.initialize();
  #endif
  rand_SIM_catch_prop_temp.allocate(1,nps,1,nr,1,nyr,1,nfl,1,2000,"rand_SIM_catch_prop_temp");
  #ifndef NO_AD_INITIALIZE
    rand_SIM_catch_prop_temp.initialize();
  #endif
  rand_SIM_catch_prop_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfl,1,2000,"rand_SIM_catch_prop_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    rand_SIM_catch_prop_temp_overlap.initialize();
  #endif
  rand_SIM_survey_prop_temp2.allocate(1,nages,"rand_SIM_survey_prop_temp2");
  #ifndef NO_AD_INITIALIZE
    rand_SIM_survey_prop_temp2.initialize();
  #endif
  rand_SIM_catch_prop_temp2.allocate(1,nages,"rand_SIM_catch_prop_temp2");
  #ifndef NO_AD_INITIALIZE
    rand_SIM_catch_prop_temp2.initialize();
  #endif
  OBS_survey_fleet_bio_temp.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nps,"OBS_survey_fleet_bio_temp");
  #ifndef NO_AD_INITIALIZE
    OBS_survey_fleet_bio_temp.initialize();
  #endif
  true_survey_fleet_bio_overlap_temp.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nps,"true_survey_fleet_bio_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    true_survey_fleet_bio_overlap_temp.initialize();
  #endif
  catch_at_age_fleet_prop_temp.allocate(1,nps,1,nr,1,nyr,1,nfl,1,nag,"catch_at_age_fleet_prop_temp");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_fleet_prop_temp.initialize();
  #endif
  abundance_move_temp.allocate(1,nps,1,nr,"abundance_move_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_move_temp.initialize();
  #endif
  bio_move_temp.allocate(1,nps,1,nr,"bio_move_temp");
  #ifndef NO_AD_INITIALIZE
    bio_move_temp.initialize();
  #endif
  yield_fleet_temp.allocate(1,nps,1,nr,1,nyr,1,nfl,1,nag,"yield_fleet_temp");
  #ifndef NO_AD_INITIALIZE
    yield_fleet_temp.initialize();
  #endif
  yield_region_temp.allocate(1,nps,1,nr,1,nyr,1,nag,"yield_region_temp");
  #ifndef NO_AD_INITIALIZE
    yield_region_temp.initialize();
  #endif
  yield_population_temp.allocate(1,nps,1,nyr,1,nag,"yield_population_temp");
  #ifndef NO_AD_INITIALIZE
    yield_population_temp.initialize();
  #endif
  SSB_region_temp.allocate(1,nps,1,nr,1,nyr,1,nag,"SSB_region_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_region_temp.initialize();
  #endif
  SSB_total_temp.allocate(1,nyr,1,nps,"SSB_total_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_total_temp.initialize();
  #endif
  abundance_population_temp.allocate(1,nps,1,nyr,1,nag,1,nr,"abundance_population_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_population_temp.initialize();
  #endif
  abundance_total_temp.allocate(1,nyr,1,nag,1,nps,"abundance_total_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_total_temp.initialize();
  #endif
  biomass_population_temp.allocate(1,nps,1,nyr,1,nr,"biomass_population_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_population_temp.initialize();
  #endif
  biomass_total_temp.allocate(1,nyr,1,nps,"biomass_total_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_total_temp.initialize();
  #endif
  catch_at_age_total_temp.allocate(1,nyr,1,nag,1,nps,"catch_at_age_total_temp");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_total_temp.initialize();
  #endif
  yield_total_temp.allocate(1,nyr,1,nps,"yield_total_temp");
  #ifndef NO_AD_INITIALIZE
    yield_total_temp.initialize();
  #endif
  catch_at_age_population_temp.allocate(1,nps,1,nyr,1,nag,1,nr,"catch_at_age_population_temp");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population_temp.initialize();
  #endif
  SSB_region_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"SSB_region_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    SSB_region_temp_overlap.initialize();
  #endif
  abundance_move_overlap_temp.allocate(1,nps,1,nr,"abundance_move_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_move_overlap_temp.initialize();
  #endif
  OBS_yield_fleet_temp.allocate(1,nps,1,nr,1,nyr,1,nfl,1,nps,"OBS_yield_fleet_temp");
  #ifndef NO_AD_INITIALIZE
    OBS_yield_fleet_temp.initialize();
  #endif
  yield_region_fleet_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"yield_region_fleet_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_region_fleet_temp_overlap.initialize();
  #endif
  yield_region_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"yield_region_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_region_temp_overlap.initialize();
  #endif
  yield_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nag,"yield_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_population_temp_overlap.initialize();
  #endif
  abundance_natal_temp_overlap.allocate(1,nps,1,nyr,1,nag,1,nps,"abundance_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    abundance_natal_temp_overlap.initialize();
  #endif
  biomass_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"biomass_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_population_temp_overlap.initialize();
  #endif
  biomass_natal_temp_overlap.allocate(1,nps,1,nyr,1,nps,"biomass_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_natal_temp_overlap.initialize();
  #endif
  catch_at_age_natal_temp_overlap.allocate(1,nps,1,nyr,1,nag,1,nps,"catch_at_age_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_natal_temp_overlap.initialize();
  #endif
  yield_natal_temp_overlap.allocate(1,nps,1,nyr,1,nps,"yield_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_natal_temp_overlap.initialize();
  #endif
  catch_at_age_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nag,1,nr,"catch_at_age_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population_temp_overlap.initialize();
  #endif
  SSB_natal_overlap_temp.allocate(1,nps,1,nyr,1,nps,"SSB_natal_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_natal_overlap_temp.initialize();
  #endif
  SSB_overlap_natal.allocate(1,nps,1,nr,"SSB_overlap_natal");
  #ifndef NO_AD_INITIALIZE
    SSB_overlap_natal.initialize();
  #endif
  abundance_AM_overlap_region_all_natal_temp.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,"abundance_AM_overlap_region_all_natal_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_AM_overlap_region_all_natal_temp.initialize();
  #endif
  SSB_population_temp.allocate(1,nps,1,nyr,1,nr,"SSB_population_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_population_temp.initialize();
  #endif
  SSB_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"SSB_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    SSB_population_temp_overlap.initialize();
  #endif
  res_TAC.allocate(1,nps,1,nr,1,nfl,1,nyr,"res_TAC");
  #ifndef NO_AD_INITIALIZE
    res_TAC.initialize();
  #endif
  res_u.allocate(1,nps,1,nr,1,nyr,"res_u");
  #ifndef NO_AD_INITIALIZE
    res_u.initialize();
  #endif
  Fnew.allocate("Fnew");
  #ifndef NO_AD_INITIALIZE
  Fnew.initialize();
  #endif
  delt.allocate("delt");
  #ifndef NO_AD_INITIALIZE
  delt.initialize();
  #endif
  fofF.allocate("fofF");
  #ifndef NO_AD_INITIALIZE
  fofF.initialize();
  #endif
  fprimeF.allocate("fprimeF");
  #ifndef NO_AD_INITIALIZE
  fprimeF.initialize();
  #endif
  fofFvect.allocate(1,nag,"fofFvect");
  #ifndef NO_AD_INITIALIZE
    fofFvect.initialize();
  #endif
  fprimeFhigh.allocate(1,nag,"fprimeFhigh");
  #ifndef NO_AD_INITIALIZE
    fprimeFhigh.initialize();
  #endif
  fprimeFlow.allocate(1,nag,"fprimeFlow");
  #ifndef NO_AD_INITIALIZE
    fprimeFlow.initialize();
  #endif
  TAC.allocate(1,nps,1,nr,1,nfl,1,nyr,"TAC");
  #ifndef NO_AD_INITIALIZE
    TAC.initialize();
  #endif
  u.allocate(1,nps,1,nr,1,nfl,"u");
  #ifndef NO_AD_INITIALIZE
    u.initialize();
  #endif
  yield_RN.allocate(1,nps,1,nr,1,nyr,1,nfl,"yield_RN");
  #ifndef NO_AD_INITIALIZE
    yield_RN.initialize();
  #endif
  yield_RN_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfl,"yield_RN_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_RN_overlap.initialize();
  #endif
  survey_RN.allocate(1,nps,1,nr,1,nyr,1,nfls,"survey_RN");
  #ifndef NO_AD_INITIALIZE
    survey_RN.initialize();
  #endif
  survey_RN_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,"survey_RN_overlap");
  #ifndef NO_AD_INITIALIZE
    survey_RN_overlap.initialize();
  #endif
  F_RN.allocate(1,nps,1,nr,1,nyr,1,nfl,"F_RN");
  #ifndef NO_AD_INITIALIZE
    F_RN.initialize();
  #endif
  rec_devs_RN.allocate(1,nps,1,nyr,"rec_devs_RN");
  #ifndef NO_AD_INITIALIZE
    rec_devs_RN.initialize();
  #endif
  Rec_apport_RN.allocate(1,nps,1,nyr,1,nr,"Rec_apport_RN");
  #ifndef NO_AD_INITIALIZE
    Rec_apport_RN.initialize();
  #endif
  rec_index_RN.allocate(1,nps,1,nr,1,nyr,"rec_index_RN");
  #ifndef NO_AD_INITIALIZE
    rec_index_RN.initialize();
  #endif
  dummy.allocate(phase_dummy,"dummy");
  ro_true_survey_fleet_overlap_age.allocate(1,nps,1,nr,1,nps,1,nyr,1,nfls,1,nag,"ro_true_survey_fleet_overlap_age");
  #ifndef NO_AD_INITIALIZE
    ro_true_survey_fleet_overlap_age.initialize();
  #endif
  ro_survey_at_age_region_fleet_overlap_prop.allocate(1,nps,1,nr,1,nps,1,nfls,1,nyr,1,nag,"ro_survey_at_age_region_fleet_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    ro_survey_at_age_region_fleet_overlap_prop.initialize();
  #endif
  ro_SIM_survey_prop_overlap.allocate(1,nps,1,nr,1,nps,1,nfls,1,nyr,1,nag,"ro_SIM_survey_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    ro_SIM_survey_prop_overlap.initialize();
  #endif
  ro_OBS_survey_prop_overlap.allocate(1,nps,1,nr,1,nps,1,nfls,1,nyr,1,nag,"ro_OBS_survey_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    ro_OBS_survey_prop_overlap.initialize();
  #endif
  ro_true_survey_fleet_overlap_age_bio.allocate(1,nps,1,nr,1,nps,1,nyr,1,nfls,1,nag,"ro_true_survey_fleet_overlap_age_bio");
  #ifndef NO_AD_INITIALIZE
    ro_true_survey_fleet_overlap_age_bio.initialize();
  #endif
  ro_catch_at_age_region_fleet_overlap.allocate(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag,"ro_catch_at_age_region_fleet_overlap");
  #ifndef NO_AD_INITIALIZE
    ro_catch_at_age_region_fleet_overlap.initialize();
  #endif
  ro_catch_at_age_region_fleet_overlap_prop.allocate(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag,"ro_catch_at_age_region_fleet_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    ro_catch_at_age_region_fleet_overlap_prop.initialize();
  #endif
  ro_SIM_catch_prop_overlap.allocate(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag,"ro_SIM_catch_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    ro_SIM_catch_prop_overlap.initialize();
  #endif
  ro_OBS_catch_prop_overlap.allocate(1,nps,1,nr,1,nps,1,nfl,1,nyr,1,nag,"ro_OBS_catch_prop_overlap");
  #ifndef NO_AD_INITIALIZE
    ro_OBS_catch_prop_overlap.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
 cout << "parameters set" << endl;
}

void model_parameters::userfunction(void)
{
  f =0.0;
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
}

void model_parameters::get_random_numbers(void)
{
   random_number_generator myrand_yield(myseed_yield);
   random_number_generator myrand_survey(myseed_survey);
   random_number_generator myrand_F(myseed_F);
   random_number_generator myrand_rec_devs(myseed_rec_devs);
   random_number_generator myrand_rec_apport(myseed_rec_apport);
   random_number_generator myrand_rec_index(myseed_rec_index);
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
              yield_RN(j,r,y,z)=randn(myrand_yield);
              yield_RN_overlap(p,j,r,y,z)=randn(myrand_yield);
              survey_RN(j,r,y,x)=randn(myrand_survey);
              survey_RN_overlap(p,j,r,y,x)=randn(myrand_survey);
              F_RN(j,r,y,z)=randn(myrand_F);
              rec_devs_RN(j,y)=randn(myrand_rec_devs);
               if(apportionment_type==3)//completely random apportionment
                {
                 Rec_apport_RN(j,y,r)=randu(myrand_rec_apport);//generate a positive random number bw 0-1
                }
               if(apportionment_type==4)//completely random apportionment
                {
                 Rec_apport_RN(j,y,r)=randn(myrand_rec_apport);//generate a positive random number bw 0-1
               }
              rec_index_RN(j,r,y)=randn(myrand_rec_index);
             }
           }
         }
       }
     }
    }
}

void model_parameters::get_movement(void)
{
  //new functional forms
  //simulate yearly differences and age-based movement
  //random movement same as a random F value
  //random walk
  //denpity dependent
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
                 if(larval_move_switch==2) // only allow movement within a population (ie regionp within a population) not across populations based on input residency term
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
                 if(larval_move_switch==3) //symmetric movement but only allow movement within a population (ie regionp within a population) not across populations
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
                 if(larval_move_switch==4) //symmetric movement across all populations and regionp
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
                }
       }
      } 
     }
    }
   }
  }
}

void model_parameters::get_selectivity(void)
{
  //yearly selectivity
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
}

void model_parameters::get_F_age(void)
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
                  if(F_year(j,r,y,z)<0) //needed because the stepF decrease can be randomly be greater than preceding F, and so F goes negative
                  {
                  F_year(j,r,y,z)=minF(j,r);
                  }
                }
              }
             F_fleet(j,r,y,a,z)=F_year(j,r,y,z)*selectivity(j,r,y,a,z);
             F(j,r,y,a)=sum(F_fleet(j,r,y,a)); 
             M(j,r,y,a)=input_M(j,a);
           }
          }
         }
        }
       }
}

void model_parameters::get_vitals(void)
{
  //random walk in apportionment or random to give time-varying
  //switch for input recruitment devs by year to recreate a given population trajectory
 R_ave=mfexp(ln_R_ave); ///this is annoying...
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
               //if(apportionment_type==5)    /// add in the two switches for shifting approtionment based on the enviromment. and random with normal dis
               }   
             }
           }         
         }
  for (int r=1;r<=nregions(j);r++)
       {
        for (int y=1;y<=nyrs;y++)
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
}

void model_parameters::get_SPR(void)
{
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
}

void model_parameters::get_env_Rec(void)
{
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
              env_rec(y)=init_abund(p,j,r,1);
             }           
           if(y>1)
            {
             env_rec(y) = R_ave(j)+(R_ave(j)*amplitude)*sin(((2*M_PI)/freq)*years(y)+freq);
            }
           }
         }
        }
       }
}

void model_parameters::get_DD_move_parameters(void)
{
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////DD Movement///////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there           
      if(move_switch==8)
       {
        for (int s=1;s<=npops;s++)
         {
          for (int n=1;n<=nregions(s);n++)
           {       
            if(use_input_Bstar==0 && Rec_type==2) //set Bstar==SSB0
             {
              if(nregions(s)==1)
               {
                Bstar(s,n)=SSB_zero(s);
               }
              if(nregions(s)>1)
               {
                Bstar(s,n)=SSB_zero(s)*SSB_zero_appor(s,n);
               }
             }
            if(use_input_Bstar==1 || Rec_type!=2)
             {
              Bstar(s,n)=input_Bstar(s,n);
             }
           c(s,n)=-log(A(s,n))/Bstar(s,n);  ///Bstar is the point of inflection in T
          }
         }
        }
}

void model_parameters::get_abundance(void)
{
       for (int y=1;y<=nyrs;y++)
        {
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
                    abundance_at_age_BM_overlap_region(p,j,y,a,r)=init_abund(p,j,r,a);
                    abundance_at_age_BM_overlap_population(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                    biomass_BM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                    biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y));
                   init_abund_temp(j,r,a,p)=init_abund(p,j,r,a);
                   abundance_at_age_BM(j,r,y,a)=sum(init_abund_temp(j,r,a));
                   recruits_BM(j,r,y)=abundance_at_age_BM(j,r,y,1);
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
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8)
                     {
                      Fract_Move_DD(j,r,y,a)=(1-input_residency(j,r,a))/(1+A(j,r)*mfexp(c(j,r)*biomass_BM_age(j,r,y,a)));
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r)); //do not include population/region moving from in summation
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
              //  abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
              //  abundance_res(j,r,y,a)=abundance_move_temp(j,r);
              //  abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
              //  bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
              //  bio_res(j,r,y,a)=bio_move_temp(j,r);
              //  bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,y,a)-bio_res(j,r,y,a);
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y);
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
                    recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y);
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
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*Rec_Prop(j,r,y);
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
                    recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y);
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
                     recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y);
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
                     recruits_BM(j,r,y)=env_rec(y)* rec_devs(j,y)*Rec_Prop(j,r,y);
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
                 abundance_at_age_BM(j,r,y,a)=recruits_BM(j,r,y);
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
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8)
                     {
                      Fract_Move_DD(j,r,y,a)=((1-input_residency(j,r,a))/(1+A(j,r)*mfexp(c(j,r)*biomass_BM_age(j,r,y,a))));
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r)); //do not include population/region moving from in summation
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
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
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
                      abundance_move_overlap_temp(k,n)=((SSB_population_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_population_overlap(p,k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
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
                      abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
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
                    abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
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
                    abundance_move_overlap_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
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
                     abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
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
                     abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
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
                    abundance_move_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
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
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
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
            abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1))*(1-tspawn(j))); //account for time of spawning (i.e., if born midway only experience a half year of mortality from age-1 to age-2)
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
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8)
                     {
                      Fract_Move_DD(j,r,y,a)=((1-input_residency(j,r,a))/(1+A(j,r)*mfexp(c(j,r)*biomass_BM_age(j,r,y,a))));
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r)); //do not include population/region moving from in summation
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
                 abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
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
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8)
                     {
                      Fract_Move_DD(j,r,y,a)=((1-input_residency(j,r,a))/(1+A(j,r)*mfexp(c(j,r)*biomass_BM_age(j,r,y,a))));
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r)); //do not include population/region moving from in summation
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
               abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM(j,r,y-1,a)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
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
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                    if(move_switch==8)
                     {
                      Fract_Move_DD(j,r,y,a)=((1-input_residency(j,r,a))/(1+A(j,r)*mfexp(c(j,r)*biomass_BM_age(j,r,y,a))));
                      biomass_BM_temp=0;
                      biomass_BM_temp2=0;
                       for (int s=1;s<=npops;s++)
                        {
                         for (int n=1;n<=nregions(s);n++)
                          {       
                           biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n);
                          }
                        }
                       biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r)); //do not include population/region moving from in summation
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                                        for(int i=1;i<=NR_iterationp;i++)  //newton-raphson iterationp
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
                  true_survey_region_abundance(j,y,r)=sum(true_survey_fleet_age_temp(j,r,y,a));
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
}

void model_parameters::get_rand_survey_CAA_prop(void)
{
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
            OBS_survey_prop(j,r,z,y,a)=SIM_survey_prop(j,r,z,y,a)/SIM_nsurvey(j,r,z);
           }
          if(use_stock_comp_info_survey==1)
           {
            OBS_survey_prop_overlap(p,j,r,z,y,a)=SIM_survey_prop_overlap(p,j,r,z,y,a)/SIM_nsurvey_overlap(p,j,r,z);
           }
         }
       }
      }
     }
    }
   }
}

void model_parameters::get_rand_CAA_prop(void)
{
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
            OBS_catch_prop(j,r,z,y,a)=SIM_catch_prop(j,r,z,y,a)/SIM_ncatch(j,r,z);
           }
          if(use_stock_comp_info_catch==1)
           {
            OBS_catch_prop_overlap(p,j,r,z,y,a)=SIM_catch_prop_overlap(p,j,r,z,y,a)/SIM_ncatch_overlap(p,j,r,z);
           }
         }
       }
      }
     }
    }
   }
}

void model_parameters::get_tag_recaptures(void)
{
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
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age 
        {
          survey_selectivity_temp(i,n,xx,1,a)=survey_selectivity(i,n,xx,a,1);
          ntags_total(x)=frac_total_abund_tagged(x)*sum(abundance_total(xx));
         // ntags(i,n,x,a)=ntags_total(x)*(true_survey_population_bio(xx,i)/true_survey_total_bio(xx))*(true_survey_region_bio(i,xx,n)/true_survey_population_bio(xx,i))*(survey_selectivity(i,n,xx,a,1)/sum(survey_selectivity_temp(i,n,xx,1)));
          ntags(i,n,x,a)=ntags_total(x)*(true_survey_population_abundance(xx,i,a)/sum(true_survey_total_abundance(xx)))*(true_survey_region_abundance(i,xx,n,a)/true_survey_population_abundance(xx,i,a)); //*(survey_selectivity(i,n,xx,a,1)/sum(survey_selectivity_temp(i,n,xx,1)));
          }
       }
      }
     }
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
}

void model_parameters::get_observed_tag_recaptures(void)
{
 if(do_tag==1)
  {
 random_number_generator myrand_tag(myseed_tag);
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
          }
        }
       }
      }
     }
 }
}

void model_parameters::evaluate_the_objective_function(void)
{
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report<<"#nages"<<endl;
  report<<nages<<endl;
  report<<"#nyrs"<<endl;
  report<<nyrs<<endl;
  report<<"#npops"<<endl;
  report<<npops<<endl;
  report<<"#nregions"<<endl;
  report<<nregions<<endl;
  report<<"#nfleets"<<endl;
  report<<nfleets<<endl;
  report<<"#nfleets_survey"<<endl;
  report<<nfleets_survey<<endl;
  report<<"#tsurvey"<<endl;
  report<<tsurvey<<endl;
  report<<"#larval_move_switch"<<endl;
  report<<larval_move_switch_EM<<endl;
  report<<"#move_switch"<<endl;
  report<<move_switch_EM<<endl;
  report<<"#natal_homing_switch"<<endl;
  report<<natal_homing_switch_EM<<endl;
  report<<"#spawn_return_switch"<<endl;
  report<<spawn_return_switch_EM<<endl;
  report<<"#select_switch"<<endl;
  report<<select_switch_EM<<endl;
  report<<"#select_switch_survey"<<endl;
  report<<select_switch_survey<<endl;
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
  report<<"#recruit_devs_switch"<<endl;
  report<<recruit_devs_switch_EM<<endl;
  report<<"#recruit_randwalk_switch"<<endl;
  report<<recruit_randwalk_switch_EM<<endl;
  report<<"#tspawn"<<endl;
  report<<tspawn<<endl;
  report<<"#return_age"<<endl;
  report<<return_age<<endl;
  report<<"#return_probability"<<endl;
  report<<return_probability<<endl;
  report<<"#spawn_return_prob"<<endl;
  report<<spawn_return_prob<<endl;
  report<<"#do_tag"<<endl;
  report<<do_tag_EM<<endl;
  report<<"#do_tag_mult"<<endl;
  report<<do_tag_mult<<endl;
  report<<"#sigma_recruit"<<endl;
  report<<sigma_recruit<<endl;
  report<<"#ph_lmr"<<endl;
  report<<ph_lmr<<endl;
  report<<"#ph_rec"<<endl;
  report<<ph_rec<<endl;
  report<<"#ph_abund_devs"<<endl;
  report<<ph_abund_devs<<endl;
  report<<"#ph_F"<<endl;
  report<<ph_F<<endl;
  report<<"#ph_steep"<<endl;
  report<<ph_steep<<endl;
  report<<"#ph_M"<<endl;
  report<<ph_M<<endl;
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
  report<<"#ph_sel_log_surv"<<endl;
  report<<ph_sel_log_surv<<endl;
  report<<"#ph_sel_dubl"<<endl;
  report<<ph_sel_dubl<<endl;
  report<<"#ph_sel_dubl_surv"<<endl;
  report<<ph_sel_dubl_surv<<endl;
  report<<"#ph_q"<<endl;
  report<<ph_q<<endl;
  report<<"#ph_F_rho"<<endl;
  report<<ph_F_rho<<endl;
  report<<"#ph_T_YR"<<endl;
  report<<ph_T_YR<<endl;
  report<<"#ph_T_CNST"<<endl;
  report<<ph_T_CNST<<endl;
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
  report<<"#abund_pen_switch"<<endl;
  report<<abund_pen_switch<<endl;
  report<<"#move_pen_switch"<<endl;
  report<<move_pen_switch<<endl;
  report<<"#Tpen"<<endl;
  report<<Tpen<<endl;
  report<<"#Tpen2"<<endl;
  report<<Tpen2<<endl;
  report<<"#input_Rec_prop"<<endl;
  report<<input_Rec_prop<<endl;
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
  report<<"#OBS_survey_fleet_bio_se"<<endl;
  report<<OBS_survey_fleet_bio_se<<endl;
  report<<"#OBS_survey_prop"<<endl;
  report<<OBS_survey_prop<<endl;
  report<<"#OBS_survey_prop_N"<<endl;
  report<<OBS_survey_prop_N<<endl;
  report<<"#OBS_yield_fleet"<<endl;
  report<<OBS_yield_fleet<<endl;
  report<<"#OBS_yield_fleet_se"<<endl;
  report<<OBS_yield_fleet_se<<endl;
  report<<"#OBS_catch_prop"<<endl;
  report<<OBS_catch_prop<<endl;
  report<<"#OBS_catch_prop_N"<<endl;
  report<<OBS_catch_prop_N<<endl;
  report<<"#nyrs_release"<<endl;
  report<<nyrs_release<<endl;
  report<<"#years_of_tag_releases "<<endl;
  report<<yrs_releases<<endl;
  report<<"#max_life_tags"<<endl;
  report<<max_life_tags<<endl;
  report<<"#report_rate"<<endl;
  report<<report_rate<<endl;
  report<<"#ntags"<<endl;
  report<<ntags<<endl;
  report<<"#ntags_total"<<endl;
  report<<ntags_total<<endl;
  report<<"#tag_N"<<endl;
  report<<tag_N<<endl;
  report<<"#input_T"<<endl;
  report<<input_T<<endl;  
  report<<"#OBS_tag_prop_final"<<endl;
  report<<OBS_tag_prop_final<<endl;
  report<<"#init_abund"<<endl;
  report<<init_abund<<endl;
  report<<"#input_M"<<endl;
  report<<input_M<<endl;
 /// TRUE VALUES
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
  report<<"#selectivity_age"<<endl;
  report<<selectivity_age<<endl;
  report<<"#survey_selectivity_age"<<endl;
  report<<survey_selectivity_age<<endl;
  report<<"#debug"<<endl;
  report<<debug<<endl;
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
}

void model_parameters::set_runtime(void)
{
  dvector temp("{.001,.0001, 1.0e-4, 1.0e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{100000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize=500000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
