  #include "admodel.h"
  #include <contrib.h>
  #define EOUT(var) cout <<#var<<" "<<var<<endl;
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <TIM_EM.htp>

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
  nfs.allocate(1,np);
 nfs=nfleets_survey;
  tsurvey.allocate(1,np,1,nreg,"tsurvey");
  larval_move_switch.allocate("larval_move_switch");
  move_switch.allocate("move_switch");
  natal_homing_switch.allocate("natal_homing_switch");
  spawn_return_switch.allocate("spawn_return_switch");
  select_switch.allocate("select_switch");
  select_switch_survey.allocate("select_switch_survey");
  maturity_switch_equil.allocate("maturity_switch_equil");
  SSB_type.allocate("SSB_type");
  Rec_type.allocate("Rec_type");
  apportionment_type.allocate("apportionment_type");
  use_stock_comp_info_survey.allocate("use_stock_comp_info_survey");
  use_stock_comp_info_catch.allocate("use_stock_comp_info_catch");
  F_switch.allocate("F_switch");
  recruit_devs_switch.allocate("recruit_devs_switch");
  recruit_randwalk_switch.allocate("recruit_randwalk_switch");
  tspawn.allocate(1,np,"tspawn");
  return_age.allocate("return_age");
  return_probability.allocate(1,np,"return_probability");
  spawn_return_prob.allocate(1,np,"spawn_return_prob");
  do_tag.allocate("do_tag");
  do_tag_mult.allocate("do_tag_mult");
  sigma_recruit.allocate(1,np,"sigma_recruit");
  ph_lmr.allocate("ph_lmr");
  ph_rec.allocate("ph_rec");
  ph_F.allocate("ph_F");
  ph_steep.allocate("ph_steep");
  ph_M.allocate("ph_M");
  ph_sel_log.allocate("ph_sel_log");
  ph_sel_log_surv.allocate("ph_sel_log_surv");
  ph_sel_dubl.allocate("ph_sel_dubl");
  ph_sel_dubl_surv.allocate("ph_sel_dubl_surv");
  ph_q.allocate("ph_q");
  ph_F_rho.allocate("ph_F_rho");
  phase_T_pop.allocate("phase_T_pop");
  phase_T_reg.allocate("phase_T_reg");
  phase_T_CNST_pop.allocate("phase_T_CNST_pop");
  phase_T_CNST_reg.allocate("phase_T_CNST_reg");
  wt_srv.allocate("wt_srv");
  wt_catch.allocate("wt_catch");
  wt_fish_age.allocate("wt_fish_age");
  wt_srv_age.allocate("wt_srv_age");
  wt_rec.allocate("wt_rec");
  wt_tag.allocate("wt_tag");
  move_pen_switch.allocate("move_pen_switch");
  Tpen.allocate("Tpen");
  Tpen2.allocate("Tpen2");
  input_Rec_prop.allocate(1,np,1,nreg,"input_Rec_prop");
  input_weight.allocate(1,np,1,nreg,1,na,"input_weight");
  input_catch_weight.allocate(1,np,1,nreg,1,na,"input_catch_weight");
  fecundity.allocate(1,np,1,nreg,1,na,"fecundity");
  maturity.allocate(1,np,1,nreg,1,na,"maturity");
  prop_fem.allocate(1,np,1,nreg,"prop_fem");
  OBS_rec_index_BM.allocate(1,np,1,nreg,1,ny,"OBS_rec_index_BM");
  OBS_survey_fleet_bio.allocate(1,np,1,nreg,1,ny,1,nfs,"OBS_survey_fleet_bio");
  OBS_survey_fleet_bio_se.allocate(1,np,1,nreg,1,ny,1,nfs,"OBS_survey_fleet_bio_se");
  OBS_survey_prop.allocate(1,np,1,nreg,1,ny,1,nfs,1,na,"OBS_survey_prop");
  OBS_survey_prop_N.allocate(1,np,1,nreg,1,ny,1,nfs,"OBS_survey_prop_N");
  OBS_yield_fleet.allocate(1,np,1,nreg,1,ny,1,nf,"OBS_yield_fleet");
  OBS_yield_fleet_se.allocate(1,np,1,nreg,1,ny,1,nf,"OBS_yield_fleet_se");
  OBS_catch_at_age_fleet_prop.allocate(1,np,1,nreg,1,ny,1,nf,1,na,"OBS_catch_at_age_fleet_prop");
  OBS_catch_at_age_fleet_prop_N.allocate(1,np,1,nreg,1,ny,1,nf,"OBS_catch_at_age_fleet_prop_N");
  nyrs_release.allocate("nyrs_release");
 int ny_rel=nyrs_release;
  yrs_releases.allocate(1,ny_rel,"yrs_releases");
  max_life_tags.allocate("max_life_tags");
 int tag_age=max_life_tags;
  report_rate.allocate(1,np,1,ny_rel,1,nreg,"report_rate");
  ntags.allocate(1,np,1,nreg,1,ny_rel,1,na,"ntags");
  ntags_total.allocate(1,ny_rel,"ntags_total");
  OBS_tag_prop_N.allocate(1,np,1,nreg,1,ny_rel,1,na,"OBS_tag_prop_N");
 int recap_index=np*nreg(1)*ny_rel; // using this to dimension down arrays
  input_T.allocate(1,np,1,nreg,1,na,1,np,1,nreg,"input_T");
 int nyr_rel=nyrs_release;
 ivector xy(1,nyr_rel);
 ivector nt(1,nyr_rel);
  for(int x=1; x<=nyrs_release; x++)
   {
    xx=yrs_releases(x);
    xy(x)=min(max_life_tags,nyrs-xx+1);
    nt(x)=xy(x)*sum(nregions)+1;
   }
  OBS_tag_prop_final.allocate(1,np,1,nreg,1,nyr_rel,1,na,1,nt,"OBS_tag_prop_final");
  input_residency_larval.allocate(1,np,1,nreg);
  input_residency.allocate(1,np,1,nreg,1,na);
  frac_total_abund_tagged.allocate(1,ny_rel);
  init_abund2.allocate(1,np,1,np,1,nreg,1,na,"init_abund2");
  input_M.allocate(1,np,1,na,"input_M");
  q_survey_TRUE.allocate(1,np,1,nreg,1,nfs,"q_survey_TRUE");
  sel_beta1_TRUE.allocate(1,np,1,nreg,1,nf,"sel_beta1_TRUE");
  sel_beta2_TRUE.allocate(1,np,1,nreg,1,nf,"sel_beta2_TRUE");
  sel_beta3_TRUE.allocate(1,np,1,nreg,1,nf,"sel_beta3_TRUE");
  sel_beta4_TRUE.allocate(1,np,1,nreg,1,nf,"sel_beta4_TRUE");
  sel_beta1_survey_TRUE.allocate(1,np,1,nreg,1,nfs,"sel_beta1_survey_TRUE");
  sel_beta2_survey_TRUE.allocate(1,np,1,nreg,1,nfs,"sel_beta2_survey_TRUE");
  sel_beta3_survey_TRUE.allocate(1,np,1,nreg,1,nfs,"sel_beta3_survey_TRUE");
  sel_beta4_survey_TRUE.allocate(1,np,1,nreg,1,nfs,"sel_beta4_survey_TRUE");
  steep_TRUE.allocate(1,np,"steep_TRUE");
  R_ave_TRUE.allocate(1,np,"R_ave_TRUE");
  SSB_zero_TRUE.allocate(1,np,"SSB_zero_TRUE");
  rec_devs_TRUE.allocate(1,np,1,ny,"rec_devs_TRUE");
  Rec_Prop_TRUE.allocate(1,np,1,nreg,1,ny,"Rec_Prop_TRUE");
  recruits_BM_TRUE.allocate(1,np,1,nreg,1,ny,"recruits_BM_TRUE");
  F_TRUE.allocate(1,np,1,nreg,1,ny,1,na,"F_TRUE");
  F_year_TRUE.allocate(1,np,1,nreg,1,ny,1,nfs,"F_year_TRUE");
  biomass_AM_TRUE.allocate(1,np,1,nreg,1,ny,"biomass_AM_TRUE");
  biomass_population_TRUE.allocate(1,np,1,ny,"biomass_population_TRUE");
  harvest_rate_region_bio_TRUE.allocate(1,np,1,nreg,1,ny,"harvest_rate_region_bio_TRUE");
  depletion_region_TRUE.allocate(1,np,1,nreg,1,ny,"depletion_region_TRUE");
  SSB_region_TRUE.allocate(1,np,1,nreg,1,ny,"SSB_region_TRUE");
  Bratio_population_TRUE.allocate(1,np,1,ny,"Bratio_population_TRUE");
  T_year_TRUE.allocate(1,np,1,nreg,1,ny,1,np,1,nreg,"T_year_TRUE");
  selectivity_age_TRUE.allocate(1,np,1,nreg,1,na,1,nf,"selectivity_age_TRUE");
  survey_selectivity_age_TRUE.allocate(1,np,1,nreg,1,na,1,nfs,"survey_selectivity_age_TRUE");
  debug.allocate("debug");
  years.allocate(1,nyrs);
years.fill_seqadd(double(1),1.0);
  nregions_temp.allocate(1,np,1,np);
 for(int j=1;j<=np;j++) //recap stock
 {
  for (int r=1;r<=np;r++) //recap region
  {
    if(j<=r)
     {
     nregions_temp(j,r)=0;
     }
    if(j>r)
     {
     nregions_temp(j,r)=nreg(j-1); //create temp matrix that holds the number of regions that exist in all previous populations (so can sum for use in calcs below)
     }
   }
  }
  nreg_temp.allocate(1,np);
 cout << "debug = " << debug << endl;
 cout << "If debug != 1541 then .dat file is meh" << endl;
 cout << "input read" << endl;
 cout << "end setup of containers" << endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
 cout << "begin parameter section" << endl;
 ivector nr=nregions;
 int parpops=npops;
 int nps=npops;
 int nyr=nyrs;
 int nag=nages;
 ivector nfl=nfleets;
 ivector nfls=nfleets_survey;  
 int k;
 int fishfleet=nfleets(1); 
 int survfleet=nfleets_survey(1);
 int parreg=nregions(1);
  ln_T_est_pop.allocate(1,nyrs,1,nps,1,nps,phase_T_pop,"ln_T_est_pop");
  ln_T_CNST_pop.allocate(1,nps,1,nps,phase_T_CNST_pop,"ln_T_CNST_pop");
  G_pop.allocate(1,nps,1,nps,"G_pop");
  #ifndef NO_AD_INITIALIZE
    G_pop.initialize();
  #endif
  G_temp_pop.allocate(1,nps,"G_temp_pop");
  #ifndef NO_AD_INITIALIZE
    G_temp_pop.initialize();
  #endif
  ln_T_est_reg.allocate(1,nyrs,1,parreg,1,parreg-1,phase_T_reg,"ln_T_est_reg");
  G_reg.allocate(1,parreg,1,parreg,"G_reg");
  #ifndef NO_AD_INITIALIZE
    G_reg.initialize();
  #endif
  G_temp_reg.allocate(1,parreg,"G_temp_reg");
  #ifndef NO_AD_INITIALIZE
    G_temp_reg.initialize();
  #endif
  ln_T_CNST_reg.allocate(1,parreg,1,parreg-1,phase_T_CNST_reg,"ln_T_CNST_reg");
  log_sel_beta1.allocate(1,parpops,1,nr,1,fishfleet,ph_sel_log,"log_sel_beta1");
  log_sel_beta2.allocate(1,parpops,1,nr,1,fishfleet,ph_sel_log,"log_sel_beta2");
  log_sel_beta3.allocate(1,parpops,1,nr,1,fishfleet,ph_sel_dubl,"log_sel_beta3");
  log_sel_beta4.allocate(1,parpops,1,nr,1,fishfleet,ph_sel_dubl,"log_sel_beta4");
  log_sel_beta1surv.allocate(1,parpops,1,survfleet,-10,5,ph_sel_log_surv,"log_sel_beta1surv");
  log_sel_beta2surv.allocate(1,parpops,1,survfleet,-10,5,ph_sel_log_surv,"log_sel_beta2surv");
  log_sel_beta3surv.allocate(1,parpops,1,survfleet,-10,5,ph_sel_dubl_surv,"log_sel_beta3surv");
  log_sel_beta4surv.allocate(1,parpops,1,survfleet,-10,5,ph_sel_dubl_surv,"log_sel_beta4surv");
  sel_beta1.allocate(1,nps,1,nr,1,nfl,"sel_beta1");
  #ifndef NO_AD_INITIALIZE
    sel_beta1.initialize();
  #endif
  sel_beta2.allocate(1,nps,1,nr,1,nfl,"sel_beta2");
  #ifndef NO_AD_INITIALIZE
    sel_beta2.initialize();
  #endif
  sel_beta3.allocate(1,nps,1,nr,1,nfl,"sel_beta3");
  #ifndef NO_AD_INITIALIZE
    sel_beta3.initialize();
  #endif
  sel_beta4.allocate(1,nps,1,nr,1,nfl,"sel_beta4");
  #ifndef NO_AD_INITIALIZE
    sel_beta4.initialize();
  #endif
  sel_beta1surv.allocate(1,nps,1,nr,1,nfls,"sel_beta1surv");
  #ifndef NO_AD_INITIALIZE
    sel_beta1surv.initialize();
  #endif
  sel_beta2surv.allocate(1,nps,1,nr,1,nfls,"sel_beta2surv");
  #ifndef NO_AD_INITIALIZE
    sel_beta2surv.initialize();
  #endif
  sel_beta3surv.allocate(1,nps,1,nr,1,nfls,"sel_beta3surv");
  #ifndef NO_AD_INITIALIZE
    sel_beta3surv.initialize();
  #endif
  sel_beta4surv.allocate(1,nps,1,nr,1,nfls,"sel_beta4surv");
  #ifndef NO_AD_INITIALIZE
    sel_beta4surv.initialize();
  #endif
  survey_selectivity.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfls,"survey_selectivity");
  #ifndef NO_AD_INITIALIZE
    survey_selectivity.initialize();
  #endif
  selectivity.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfl,"selectivity");
  #ifndef NO_AD_INITIALIZE
    selectivity.initialize();
  #endif
  survey_selectivity_age.allocate(1,nps,1,nr,1,nag,1,nfls,"survey_selectivity_age");
  #ifndef NO_AD_INITIALIZE
    survey_selectivity_age.initialize();
  #endif
  selectivity_age.allocate(1,nps,1,nr,1,nag,1,nfl,"selectivity_age");
  #ifndef NO_AD_INITIALIZE
    selectivity_age.initialize();
  #endif
  ln_q.allocate(1,parpops,1,survfleet,-30,15,ph_q,"ln_q");
  q_survey.allocate(1,parpops,1,nr,1,nfls,"q_survey");
  #ifndef NO_AD_INITIALIZE
    q_survey.initialize();
  #endif
  ln_F.allocate(1,parpops,1,nyr,1,fishfleet,ph_F,"ln_F");
  F_rho.allocate(1,parpops,1,fishfleet,0,2,ph_F_rho,"F_rho");
  F_fleet.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfl,"F_fleet");
  #ifndef NO_AD_INITIALIZE
    F_fleet.initialize();
  #endif
  F_year.allocate(1,nps,1,nr,1,nyr,1,nfl,"F_year");
  #ifndef NO_AD_INITIALIZE
    F_year.initialize();
  #endif
  F.allocate(1,nps,1,nr,1,nyr,1,nag,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  M.allocate(1,nps,1,nr,1,nyr,1,nag,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
 int nyr_rel=nyrs_release;
 ivector xy(1,nyr_rel);
 ivector nt(1,nyr_rel);
 ivector nt2(1,nyr_rel);
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
  total_rec.allocate(1,nps,1,nr,1,nyr_rel,1,nag,"total_rec");
  #ifndef NO_AD_INITIALIZE
    total_rec.initialize();
  #endif
  not_rec.allocate(1,nps,1,nr,1,nyr_rel,1,nag,"not_rec");
  #ifndef NO_AD_INITIALIZE
    not_rec.initialize();
  #endif
  tag_prop.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr,"tag_prop");
  #ifndef NO_AD_INITIALIZE
    tag_prop.initialize();
  #endif
  tag_prop_not_rec.allocate(1,nps,1,nr,1,nyr_rel,1,nag,"tag_prop_not_rec");
  #ifndef NO_AD_INITIALIZE
    tag_prop_not_rec.initialize();
  #endif
  tag_prop_final.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,nt,"tag_prop_final");
  #ifndef NO_AD_INITIALIZE
    tag_prop_final.initialize();
  #endif
  T.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr,"T");
  #ifndef NO_AD_INITIALIZE
    T.initialize();
  #endif
  T_terminal.allocate(1,nps,1,nr,1,nag,1,nps,1,nr,"T_terminal");
  #ifndef NO_AD_INITIALIZE
    T_terminal.initialize();
  #endif
  T_year.allocate(1,nps,1,nr,1,nyr,1,nps,1,nr,"T_year");
  #ifndef NO_AD_INITIALIZE
    T_year.initialize();
  #endif
  recaps.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr,"recaps");
  #ifndef NO_AD_INITIALIZE
    recaps.initialize();
  #endif
  ln_rec_prop.allocate(1,parpops,1,parreg,-4,0,-1,"ln_rec_prop");
  ln_rec_prop_year.allocate(1,parpops,1,parreg,1,nyr,-1,"ln_rec_prop_year");
  ln_rec_devs_RN.allocate(1,parpops,1,nyr+nages-1,-40,40,ph_rec,"ln_rec_devs_RN");
  rec_devs.allocate(1,nps,1,nyr+nages-1,"rec_devs");
  #ifndef NO_AD_INITIALIZE
    rec_devs.initialize();
  #endif
  ln_R_ave.allocate(1,parpops,2,20,ph_lmr,"ln_R_ave");
  R_ave.allocate(1,parpops,"R_ave");
  #ifndef NO_AD_INITIALIZE
    R_ave.initialize();
  #endif
  SSB_zero.allocate(1,nps,"SSB_zero");
  #ifndef NO_AD_INITIALIZE
    SSB_zero.initialize();
  #endif
  steep.allocate(1,parpops,0.2,1,ph_steep,"steep");
  alpha.allocate(1,nps,"alpha");
  #ifndef NO_AD_INITIALIZE
    alpha.initialize();
  #endif
  beta.allocate(1,nps,"beta");
  #ifndef NO_AD_INITIALIZE
    beta.initialize();
  #endif
  init_abund.allocate(1,nps,1,nps,1,nr,1,nag,"init_abund");
  #ifndef NO_AD_INITIALIZE
    init_abund.initialize();
  #endif
  survey_fleet_overlap_age.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag,"survey_fleet_overlap_age");
  #ifndef NO_AD_INITIALIZE
    survey_fleet_overlap_age.initialize();
  #endif
  survey_at_age_region_fleet_overlap_prop.allocate(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag,"survey_at_age_region_fleet_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    survey_at_age_region_fleet_overlap_prop.initialize();
  #endif
  survey_fleet_overlap_age_bio.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag,"survey_fleet_overlap_age_bio");
  #ifndef NO_AD_INITIALIZE
    survey_fleet_overlap_age_bio.initialize();
  #endif
  catch_at_age_region_fleet_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"catch_at_age_region_fleet_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region_fleet_overlap.initialize();
  #endif
  catch_at_age_region_fleet_overlap_prop.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"catch_at_age_region_fleet_overlap_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_region_fleet_overlap_prop.initialize();
  #endif
  yield_region_fleet_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag,"yield_region_fleet_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_region_fleet_temp_overlap.initialize();
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
  recruits_BM.allocate(1,nps,1,nr,1,nyr,"recruits_BM");
  #ifndef NO_AD_INITIALIZE
    recruits_BM.initialize();
  #endif
  recruits_AM.allocate(1,nps,1,nr,1,nyr,"recruits_AM");
  #ifndef NO_AD_INITIALIZE
    recruits_AM.initialize();
  #endif
  rec_index_BM.allocate(1,nps,1,nr,1,nyr,"rec_index_BM");
  #ifndef NO_AD_INITIALIZE
    rec_index_BM.initialize();
  #endif
  rec_index_AM.allocate(1,nps,1,nr,1,nyr,"rec_index_AM");
  #ifndef NO_AD_INITIALIZE
    rec_index_AM.initialize();
  #endif
  rec_index_prop_BM.allocate(1,nps,1,nr,1,nyr,"rec_index_prop_BM");
  #ifndef NO_AD_INITIALIZE
    rec_index_prop_BM.initialize();
  #endif
  rec_index_BM_temp.allocate(1,nps,1,nyr,1,nr,"rec_index_BM_temp");
  #ifndef NO_AD_INITIALIZE
    rec_index_BM_temp.initialize();
  #endif
  rec_index_prop_AM.allocate(1,nps,1,nr,1,nyr,"rec_index_prop_AM");
  #ifndef NO_AD_INITIALIZE
    rec_index_prop_AM.initialize();
  #endif
  rec_index_AM_temp.allocate(1,nps,1,nyr,1,nr,"rec_index_AM_temp");
  #ifndef NO_AD_INITIALIZE
    rec_index_AM_temp.initialize();
  #endif
  rec_devs_randwalk.allocate(1,nps,1,nyr-1,"rec_devs_randwalk");
  #ifndef NO_AD_INITIALIZE
    rec_devs_randwalk.initialize();
  #endif
  Rec_Prop.allocate(1,nps,1,nr,1,nyr,"Rec_Prop");
  #ifndef NO_AD_INITIALIZE
    Rec_Prop.initialize();
  #endif
  Rec_prop_temp1.allocate(1,nps,1,nr,1,nyr,"Rec_prop_temp1");
  #ifndef NO_AD_INITIALIZE
    Rec_prop_temp1.initialize();
  #endif
  Rec_prop_temp2.allocate(1,nps,1,nr,1,nyr,"Rec_prop_temp2");
  #ifndef NO_AD_INITIALIZE
    Rec_prop_temp2.initialize();
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
  catch_at_age_fleet.allocate(1,nps,1,nr,1,nyr,1,nag,1,nfl,"catch_at_age_fleet");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_fleet.initialize();
  #endif
  catch_at_age_fleet_prop.allocate(1,nps,1,nr,1,nyr,1,nfl,1,nag,"catch_at_age_fleet_prop");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_fleet_prop.initialize();
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
  biomass_BM_overlap_temp.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,"biomass_BM_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_BM_overlap_temp.initialize();
  #endif
  init_abund_temp.allocate(1,nps,1,nr,1,nag,1,nps,"init_abund_temp");
  #ifndef NO_AD_INITIALIZE
    init_abund_temp.initialize();
  #endif
  survey_fleet_bio_overlap_temp.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nps,"survey_fleet_bio_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    survey_fleet_bio_overlap_temp.initialize();
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
  abundance_move_overlap_temp.allocate(1,nps,1,nr,"abundance_move_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_move_overlap_temp.initialize();
  #endif
  abundance_AM_overlap_region_all_natal_temp.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,"abundance_AM_overlap_region_all_natal_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_AM_overlap_region_all_natal_temp.initialize();
  #endif
  biomass_AM_overlap_region_all_natal_temp.allocate(1,nps,1,nr,1,nyr,1,nag,1,nps,"biomass_AM_overlap_region_all_natal_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_AM_overlap_region_all_natal_temp.initialize();
  #endif
  SSB_natal_overlap_temp.allocate(1,nps,1,nyr,1,nps,"SSB_natal_overlap_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_natal_overlap_temp.initialize();
  #endif
  abundance_natal_temp_overlap.allocate(1,nps,1,nyr,1,nag,1,nps,"abundance_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    abundance_natal_temp_overlap.initialize();
  #endif
  SSB_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"SSB_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    SSB_population_temp_overlap.initialize();
  #endif
  SSB_region_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"SSB_region_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    SSB_region_temp_overlap.initialize();
  #endif
  survey_fleet_bio_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nfls,"survey_fleet_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    survey_fleet_bio_overlap.initialize();
  #endif
  survey_region_bio_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"survey_region_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    survey_region_bio_overlap.initialize();
  #endif
  survey_population_bio_overlap.allocate(1,nps,1,nyr,1,nps,"survey_population_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    survey_population_bio_overlap.initialize();
  #endif
  survey_natal_bio_overlap.allocate(1,nyr,1,nps,"survey_natal_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    survey_natal_bio_overlap.initialize();
  #endif
  survey_total_bio_overlap.allocate(1,nyr,"survey_total_bio_overlap");
  #ifndef NO_AD_INITIALIZE
    survey_total_bio_overlap.initialize();
  #endif
  survey_fleet_age.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nag,"survey_fleet_age");
  #ifndef NO_AD_INITIALIZE
    survey_fleet_age.initialize();
  #endif
  survey_at_age_fleet_prop.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nag,"survey_at_age_fleet_prop");
  #ifndef NO_AD_INITIALIZE
    survey_at_age_fleet_prop.initialize();
  #endif
  survey_fleet_age_bio.allocate(1,nps,1,nr,1,nyr,1,nfls,1,nag,"survey_fleet_age_bio");
  #ifndef NO_AD_INITIALIZE
    survey_fleet_age_bio.initialize();
  #endif
  survey_fleet_bio.allocate(1,nps,1,nr,1,nyr,1,nfls,"survey_fleet_bio");
  #ifndef NO_AD_INITIALIZE
    survey_fleet_bio.initialize();
  #endif
  survey_region_bio.allocate(1,nps,1,nyr,1,nr,"survey_region_bio");
  #ifndef NO_AD_INITIALIZE
    survey_region_bio.initialize();
  #endif
  survey_population_bio.allocate(1,nyr,1,nps,"survey_population_bio");
  #ifndef NO_AD_INITIALIZE
    survey_population_bio.initialize();
  #endif
  survey_total_bio.allocate(1,nyr,"survey_total_bio");
  #ifndef NO_AD_INITIALIZE
    survey_total_bio.initialize();
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
  SSB_overlap_natal.allocate(1,nps,1,nr,"SSB_overlap_natal");
  #ifndef NO_AD_INITIALIZE
    SSB_overlap_natal.initialize();
  #endif
  yield_region_temp_overlap.allocate(1,nps,1,nps,1,nr,1,nyr,1,nag,"yield_region_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_region_temp_overlap.initialize();
  #endif
  yield_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nag,"yield_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_population_temp_overlap.initialize();
  #endif
  yield_natal_temp_overlap.allocate(1,nps,1,nyr,1,nps,"yield_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    yield_natal_temp_overlap.initialize();
  #endif
  catch_at_age_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nag,1,nr,"catch_at_age_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population_temp_overlap.initialize();
  #endif
  catch_at_age_natal_temp_overlap.allocate(1,nps,1,nyr,1,nag,1,nps,"catch_at_age_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_natal_temp_overlap.initialize();
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
  catch_at_age_total_temp.allocate(1,nyr,1,nag,1,nps,"catch_at_age_total_temp");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_total_temp.initialize();
  #endif
  catch_at_age_population_temp.allocate(1,nps,1,nyr,1,nag,1,nr,"catch_at_age_population_temp");
  #ifndef NO_AD_INITIALIZE
    catch_at_age_population_temp.initialize();
  #endif
  yield_total_temp.allocate(1,nyr,1,nps,"yield_total_temp");
  #ifndef NO_AD_INITIALIZE
    yield_total_temp.initialize();
  #endif
  SSB_region_temp.allocate(1,nps,1,nr,1,nyr,1,nag,"SSB_region_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_region_temp.initialize();
  #endif
  SSB_total_temp.allocate(1,nyr,1,nps,"SSB_total_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_total_temp.initialize();
  #endif
  SSB_population_temp.allocate(1,nps,1,nyr,1,nr,"SSB_population_temp");
  #ifndef NO_AD_INITIALIZE
    SSB_population_temp.initialize();
  #endif
  biomass_population_temp.allocate(1,nps,1,nyr,1,nr,"biomass_population_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_population_temp.initialize();
  #endif
  biomass_total_temp.allocate(1,nyr,1,nps,"biomass_total_temp");
  #ifndef NO_AD_INITIALIZE
    biomass_total_temp.initialize();
  #endif
  biomass_population_temp_overlap.allocate(1,nps,1,nps,1,nyr,1,nr,"biomass_population_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_population_temp_overlap.initialize();
  #endif
  biomass_natal_temp_overlap.allocate(1,nps,1,nyr,1,nps,"biomass_natal_temp_overlap");
  #ifndef NO_AD_INITIALIZE
    biomass_natal_temp_overlap.initialize();
  #endif
  abundance_population_temp.allocate(1,nps,1,nyr,1,nag,1,nr,"abundance_population_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_population_temp.initialize();
  #endif
  abundance_total_temp.allocate(1,nyr,1,nag,1,nps,"abundance_total_temp");
  #ifndef NO_AD_INITIALIZE
    abundance_total_temp.initialize();
  #endif
  total_recap_temp.allocate(1,nps,1,nr,1,tag_age,"total_recap_temp");
  #ifndef NO_AD_INITIALIZE
    total_recap_temp.initialize();
  #endif
  tags_avail_temp.allocate(1,nps,1,nr,"tags_avail_temp");
  #ifndef NO_AD_INITIALIZE
    tags_avail_temp.initialize();
  #endif
  tag_prop_temp.allocate(1,nps,1,nyr_rel,1,nr,"tag_prop_temp");
  #ifndef NO_AD_INITIALIZE
    tag_prop_temp.initialize();
  #endif
  tag_prop_temp2.allocate(1,nps,1,nr,1,nyr_rel,1,nag,1,nt2,"tag_prop_temp2");
  #ifndef NO_AD_INITIALIZE
    tag_prop_temp2.initialize();
  #endif
  survey_age_like.allocate("survey_age_like");
  #ifndef NO_AD_INITIALIZE
  survey_age_like.initialize();
  #endif
  fish_age_like.allocate("fish_age_like");
  #ifndef NO_AD_INITIALIZE
  fish_age_like.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  tag_like.allocate("tag_like");
  #ifndef NO_AD_INITIALIZE
  tag_like.initialize();
  #endif
  tag_like_temp.allocate("tag_like_temp");
  #ifndef NO_AD_INITIALIZE
  tag_like_temp.initialize();
  #endif
  catch_like.allocate("catch_like");
  #ifndef NO_AD_INITIALIZE
  catch_like.initialize();
  #endif
  survey_like.allocate("survey_like");
  #ifndef NO_AD_INITIALIZE
  survey_like.initialize();
  #endif
  Tpen_like.allocate("Tpen_like");
  #ifndef NO_AD_INITIALIZE
  Tpen_like.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
 cout << "parameters set" << endl;
}

void model_parameters::initializationfunction(void)
{
  steep.set_initial_value(.814);
}

void model_parameters::userfunction(void)
{
  f =0.0;
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
}

void model_parameters::get_movement(void)
{
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
 if(move_switch==1) // if T fixed set it to input T
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
              T(j,r,y,a,k,n)=input_T(j,r,a,k,n);            
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
}

void model_parameters::get_selectivity(void)
{
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
 sel_beta1(j,r,z)=mfexp(log_sel_beta1(j,r,z));
 sel_beta2(j,r,z)=mfexp(log_sel_beta2(j,r,z));
 sel_beta3(j,r,z)=mfexp(log_sel_beta3(j,r,z));
 sel_beta4(j,r,z)=mfexp(log_sel_beta4(j,r,z));
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
             if(ph_F<0) //inefficient place for this because it'll reassign in all these loops but F_fleet needs to be in loops below...
             {
              F_year=F_year_TRUE;
             }
             else {
             if(F_switch==1) //estimate annual deviations
              {
               F_year(j,r,y,z)=mfexp(ln_F(j,y,z)); 
              }
             if(F_switch==2) //random walk or AR1  in F if we ever want to try it.
              {
               if(y==1) //year one separate because no y-1  
               {
               F_year(j,r,y,z)=mfexp(ln_F(j,y,z));  
               }
               if(y>1)
               {
               F_year(j,r,y,z)=F_rho(j,z)*mfexp(ln_F(j,y-1,z));            
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
}

void model_parameters::get_vitals(void)
{
       for(int i=1;i<=npops;i++) 
    {
        for(int k=1;k<=nregions(i);k++) 
        {
          for(int j=1;j<=nfleets_survey(i);j++) 
          {
           q_survey(i,k,j) = mfexp(ln_q(i,j));
    }
    }
    }
  R_ave=mfexp(ln_R_ave+square(sigma_recruit)*0.5);
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
        for (int y=1;y<=nyrs+nages-1;y++)
         {
               if(recruit_devs_switch==0)  //use population recruit relationship directly
                {
                 rec_devs(j,y)=1;
                }
               if(recruit_devs_switch==1)  // allow lognormal error around SR curve
                {
                 rec_devs(j,y)=mfexp(ln_rec_devs_RN(j,y)*sigma_recruit(j)-.5*square(sigma_recruit(j)));
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
                    init_abund(p,j,r,a)=R_ave(p)*rec_devs(j,a)*pow(mfexp(-(M(p,j,r,a))),a);  //not sure what this is doing; JJD: I think this is estimating init_abundance as deviations from an exponential decline from Rave
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
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y+nages-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y+nages-1)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regions within a population
                    {
                    recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y+nages-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a population
                    {
                     recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y+nages-1)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1)));
                    }
                  }
                if(Rec_type==3) //environmental recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y+nages-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a population
                    {
                     recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y+nages-1)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
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
                    recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y+nages-1)*Rec_Prop(j,r,y);
                    }
                    if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y+nages-1)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y+nages-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y+nages-1)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                    }
                   }
                if(Rec_type==3) //average recruitment
                  {
                  if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     recruits_BM(j,r,y)=env_rec(y)* rec_devs(j,y+nages-1)*Rec_Prop(j,r,y);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                   {
                    recruits_BM(j,r,y)=env_rec(y)* R_ave(j)*rec_devs(j,y+nages-1)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                   }
                 }
                 }
               rec_index_BM(j,r,y)=recruits_BM(j,r,y);
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
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y+nages-1)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_population_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_population_overlap(p,k,y-1)))*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=((SSB_population_overlap(p,k,y-1))/(alpha(k)+beta(k)*SSB_population_overlap(p,k,y-1)))*rec_devs(k,y+nages-1)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r);
                    }
                  }
                if(Rec_type==3) //average recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(p,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                      abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y+nages-1)*(SSB_region_overlap(p,k,n,y-1)/sum(SSB_region_overlap(p,k,y-1)))*T(p,n,y,a,j,r); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
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
                    abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    abundance_move_overlap_temp(k,n)=R_ave(k)*rec_devs(k,y+nages-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                    }
                 if(Rec_type==2) //BH recruitment
                  {
                 if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    abundance_move_overlap_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y+nages-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                if(Rec_type==3) //environmental recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_overlap_temp(k,n)=env_rec(y)*rec_devs(k,y+nages-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
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
                     abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                     }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                    abundance_move_temp(k,n)=R_ave(k)*rec_devs(k,y+nages-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                 if(Rec_type==2) //BH recruitment
                  {
                   if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1))  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                    abundance_move_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=((SSB_population(k,y-1))/(alpha(k)+beta(k)*SSB_population(k,y-1)))*rec_devs(k,y+nages-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
                    }
                  }
                if(Rec_type==3) //env recruitment
                  {
                  if(apportionment_type==1 || apportionment_type==2 ||apportionment_type==3||apportionment_type==4||apportionment_type==(-1)) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y+nages-1)*Rec_Prop(k,n,y)*T(k,n,y,a,j,r);
                    }
                   if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                    {
                     abundance_move_temp(k,n)=env_rec(y)*rec_devs(k,y+nages-1)*(SSB_region(k,n,y-1)/SSB_population(k,y-1))*T(k,n,y,a,j,r);
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
}

void model_parameters::get_survey_CAA_prop(void)
{
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
}

void model_parameters::get_CAA_prop(void)
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

void model_parameters::evaluate_the_objective_function(void)
{
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
 for(int j=1;j<=npops;j++) {  //is this correct?  we aren't penalizing against devs from Rave or SR function?  also we do want to do calcs across all years right (since not including year index in rec_devs here)?
 rec_like += (norm2(ln_rec_devs_RN(j)+sigma_recruit(j)*sigma_recruit(j)/2.)/(2.*square(sigma_recruit(j))) + (size_count(ln_rec_devs_RN(j)))*log(sigma_recruit(j))); 
 }
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
 if (current_phase()<4) f+= 1000*(norm2(mfexp(ln_F)));
   f           += survey_like*wt_srv;
   f           += catch_like*wt_catch;
   f           += fish_age_like*wt_fish_age;
   f           += survey_age_like*wt_srv_age;
   f           += rec_like*wt_rec;
   f           += tag_like*wt_tag;
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
  report<<"$Init_Abund_TRUE"<<endl;
  report<<init_abund2<<endl;
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
  save_gradients(gradients);
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
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
