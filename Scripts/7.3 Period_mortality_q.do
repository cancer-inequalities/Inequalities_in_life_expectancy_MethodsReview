/****************************************************************************************
* Mortality bounds by education quintile and sex
* Rank-based bounds (Novosad et al.)
*
* OUTPUT:
*   age sex quintile mort_lb mort_ub}
*
********************************************************************************
* NOTE ON PATHS
********************************************************************************
*
* This script assumes that the working directory is set to the root of the
* repository. All file paths are defined relative to this location.
*
****************************************************************************************/

clear
set more off
version 17

********************************************************************************
* 0. LOAD DATA WITH MORTALITY + RANKS
********************************************************************************
use "Out/mortality_rates_2017_ranked_period.dta", clear

confirm variable age
confirm variable sex
confirm variable rate
confirm variable rank

********************************************************************************
* 1. CONVERT MORTALITY → SURVIVAL
********************************************************************************
gen surv = 100000 - rate
label variable surv "Survival per 100,000"

********************************************************************************
* 2. DEFINE QUINTILE CUTS
********************************************************************************
local qcuts 0 20 40 60 80 100

********************************************************************************
* 3. CREATE BOUND VARIABLES (WIDE FORMAT)
********************************************************************************
forvalues q = 1/5 {
    gen surv_lb_q`q' = .
    gen surv_ub_q`q' = .
    gen mort_lb_q`q' = .
    gen mort_ub_q`q' = .
}

********************************************************************************
* 4. COMPUTE SURVIVAL BOUNDS (AGE × SEX × QUINTILE)
********************************************************************************
levelsof age, local(ages)
levelsof sex, local(sexes)

foreach s of local sexes {
    forvalues q = 1/5 {

        local s_cut : word `q' of `qcuts'
        local t_cut : word `=`q'+1' of `qcuts'

        foreach a of local ages {

            quietly bound_param if age == `a' & sex == `s', ///
                 xvar(rank) ///
				yvar(surv) ///
				s(`s_cut') t(`t_cut') ///
				minmom(0) maxmom(100000) ///
				forcemono
            replace surv_lb_q`q' = real(r(mu_lb)) if age == `a' & sex == `s'
            replace surv_ub_q`q' = real(r(mu_ub)) if age == `a' & sex == `s'
        }
    }
}

********************************************************************************
* 5. SURVIVAL → MORTALITY BOUNDS
********************************************************************************
forvalues q = 1/5 {
    replace mort_lb_q`q' = 100000 - surv_ub_q`q'
    replace mort_ub_q`q' = 100000 - surv_lb_q`q'
}

********************************************************************************
* 6. COLLAPSE TO ONE ROW PER AGE × SEX
*    (required before reshape)
********************************************************************************
collapse (mean) mort_lb_q* mort_ub_q*, by(age sex)

isid age sex   // must pass

********************************************************************************
* 7. RESHAPE TO LONG FORMAT (AGE × SEX × QUINTILE)
********************************************************************************
reshape long mort_lb_q mort_ub_q, i(age sex) j(quintile)

rename mort_lb_q mort_lb
rename mort_ub_q mort_ub

order age sex quintile mort_lb mort_ub
sort sex age quintile

********************************************************************************
* 8. LABELS
********************************************************************************
label define sexlbl 1 "Men" 2 "Women"
label values sex sexlbl
label variable mort_lb "Lower bound mortality (per 100,000)"
label variable mort_ub "Upper bound mortality (per 100,000)"
label variable quintile "Education quintile"

********************************************************************************
* 9. FINAL SANITY CHECKS
********************************************************************************
* Bounds ordered correctly
assert mort_lb <= mort_ub if !missing(mort_lb)

* One row per age × sex × quintile
isid age sex quintile

********************************************************************************
* 10. SAVE OUTPUTS
********************************************************************************
save ///
"Out/mortality_bounds_by_age_sex_quintile_period.dta", ///
replace

********************************************************************************
* END OF DO-FILE
********************************************************************************

