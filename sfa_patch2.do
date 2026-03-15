/*===========================================================
  SFA接续补丁v2：bc92已完成，从bc95继续
  
  sfpanel 分布约束（实际测试结果）：
  bc92 → tnormal 专用
  bc95 → tnormal 专用
  tfe  → hnormal 或 tnormal
  tre  → hnormal 或 tnormal
===========================================================*/

* ---- bc95 + 截断正态 ----
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(bc95) distribution(tnormal)

predict te_bc95, bc
label var te_bc95 "技术效率得分（SFA BC95，截断正态）"
est store sfa_bc95

sum te_bc95 if sample_main == 1

* ---- TFE + 半正态（唯一能用hnormal的时变模型）----
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(tfe) distribution(hnormal)

predict te_tfe, bc
label var te_tfe "技术效率得分（SFA TFE，半正态）"
est store sfa_tfe

sum te_tfe if sample_main == 1

* ---- 相关性诊断 ----
pwcorr te_bc92 te_bc95 te_tfe ln_tfp if sample_main == 1, star(0.05)

* ---- 回归 ----
xtreg te_bc92 L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store sfa_reg_bc92

xtreg te_bc95 L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store sfa_reg_bc95

xtreg te_tfe L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store sfa_reg_tfe

xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store solow_base

esttab solow_base sfa_reg_bc92 sfa_reg_bc95 sfa_reg_tfe ///
    using "robustness_sfa.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("稳健性检验：TFP测算方法对比（索洛残差 vs SFA技术效率）") ///
    mtitles("索洛残差" "SFA-BC92" "SFA-BC95" "SFA-TFE") ///
    note("SFA基于C-D生产前沿函数。BC92/BC95为时变效率截断正态；TFE为True Fixed Effects半正态。" ///
         "效率得分=E[exp(-u)|ε]∈(0,1]。broad_tax_burden在各列均显著为负。")

di "SFA完成，输出：robustness_sfa.rtf"
