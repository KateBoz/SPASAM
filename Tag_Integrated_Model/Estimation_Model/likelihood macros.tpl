
//NEED TO:  generalize using the correct loop values (i.e., j, r, z, y, a, etc.; sum over all options to add to likelihood function;
//          make sure correct when putting into the EM (once coded)
  f_survey(j,r)=lk_lognormal(pred_survey(j,r), obs_survey(j,r), survey_cv(j,r), w_survey(j,r));
  f_catch(j,r,z)=lk_lognormal(pred_catch(j,r,z), obs_catch(j,r,z), catch_cv(j,r,z), w_catch(j,r,z));
  
  f_catch_comps(j,r,y,z,a)=lk_multinomial(nsamp_catch_comps(j,r,y,z,a), pred_catch_comps(j,r,y,z,a), obs_catch_comps(j,r,y,z,a), nyrs_catch_comps(j,r,y,z,a), minSS_catch_comps(j,r,y,z,a), w_catch_comps(j,r,z));
  f_survey_comps(j,r,y,z,a)=lk_multinomial(nsamp_survey_comps(j,r,y,z,a), pred_survey_comps(j,r,y,z,a), obs_survey_comps(j,r,y,z,a), nyrs_survey_comps(j,r,q,z,a), minSS_survey_comps(j,r,y,z,a), w_survey_comps(j,r,z));

   
//-----------------------------------------------------------------------------------
//Likelihood contribution: lognormal
FUNCTION dvariable lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
  //small_number is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+small_number),(obs+small_number)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  LkvalTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+small_number),
               log(elem_div((pred_comp(ii)+small_number), (obs_comp(ii)+small_number)))));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

