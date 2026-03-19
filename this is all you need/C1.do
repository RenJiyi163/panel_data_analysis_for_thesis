preserve
keep if sample_sdm == 1
xtset province_id year

xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    , wmat(W_econ_raw) model(sdm) fe vce(robust)
est store sdm_econ

xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    , wmat(W_econ_raw) model(sdm) fe impacts vce(robust)

xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    , wmat(W_econ_raw) model(sar) fe vce(robust)
est store sar_econ

xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    , wmat(W_adj_raw) model(sdm) fe vce(robust)
est store sdm_adj

xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    , wmat(W_adj_raw) model(sdm) fe impacts vce(robust)

esttab sdm_econ sar_econ using "spatial_main.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间计量模型（经济距离权重矩阵，2011—2017年）") ///
    mtitles("SDM（杜宾）" "SAR（自回归）") ///
    note("经济距离矩阵：w_ij=1/|人均GDP差距|，行标准化。N=217（31省×7年）。")

esttab sdm_econ sdm_adj using "spatial_robust.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间权重矩阵稳健性：经济距离 vs 地理邻接") ///
    mtitles("SDM经济距离" "SDM地理邻接")

restore