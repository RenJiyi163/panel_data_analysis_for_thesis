/*===========================================================
  接续补丁：从模型3开始，修复 ln_pub_emp 未定义的问题
  在已跑完 m1_revised, m2_revised 的基础上接续执行
===========================================================*/

* 补充生成遗漏变量
gen ln_pub_emp = ln(pub_emp_share2)
label var ln_pub_emp "ln(公共部门就业占比)"

gen L3_ln_pub   = L3.ln_pub_emp
gen L3_ln_aging = L3.ln_aging
label var L3_ln_pub   "ln(公共部门占比，滞后3期)"
label var L3_ln_aging "ln(老龄化率，滞后3期)"

* -------------------------------------------------------
* 模型3：对数形式
* -------------------------------------------------------
xtreg ln_tfp L3_ln_pub broad_tax_burden L3_ln_aging ///
    ln_urban sec_ind_share ln_trade ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store m3_revised

* -------------------------------------------------------
* 模型4：剔除西藏
* -------------------------------------------------------
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1 & province_id != 26, ///
    fe vce(cluster province_id)
est store m4_revised

esttab m1_revised m2_revised m3_revised m4_revised ///
    using "main_regression_revised.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表2（修订）双向固定效应主回归结果（因变量：ln_TFP，2011-2017年）") ///
    mtitles("模型1" "模型2（含控制）" "模型3（对数）" "模型4（剔除西藏）") ///
    note("括号内为省份层面聚类稳健标准误。*** p<0.01, ** p<0.05, * p<0.1。" ///
         "pub_emp_share2和aging_rate均取滞后3期（L3）。" ///
         "有效样本：模型1为2011-2017年N=217；模型2/3/4因trade_openness缺失N=186/186/180。")

* -------------------------------------------------------
* Between Effects：省间截面结构效应
* 这是识别"公共部门虹吸效应"的正确估计量
* FE只利用省内时序变化，BE利用省份均值的截面差异
* 虹吸效应是长期结构性均衡，体现为省份间持久差异
* -------------------------------------------------------
xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    if sample_main == 1, be
est store be_model

esttab m2_revised be_model using "fe_vs_be.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表2补充：FE省内时序效应 vs BE省间截面结构效应") ///
    mtitles("FE-L3滞后（动态因果）" "BE截面（结构效应）") ///
    note("FE估计量识别省内时序动态效应（短期因果）；" ///
         "BE估计量识别省间截面均衡差异（长期结构效应）。" ///
         "公共部门虹吸效应属于长期结构机制，BE结果是其主要经验证据。")

* -------------------------------------------------------
* H2中介效应：科技人员密度
* -------------------------------------------------------

* 第二步：广义税负/公共部门占比 → 科技人员密度
xtreg ln_tech_density L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med2_tech

* 第三步：含中介变量的TFP方程
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    ln_tech_density ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med3_tech

esttab m2_revised med2_tech med3_tech using "mediation_revised.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表3（修订）中介效应：制度约束→科技人力资本行政性占用→TFP下降") ///
    mtitles("第一步：TFP总效应" "第二步：→ln(科技人员密度)" "第三步：→ln(TFP)含中介") ///
    note("科技人员密度=科技活动人员数/GDP（十亿元）×1000；" ///
         "密度越高代表单位经济产出消耗的科技人力资本越多，" ///
         "是人力资本行政性占用的逆向代理指标。" ///
         "预期符号：第二步broad_tax→density（正）；第三步density→TFP（负）。")

* 计算并报告间接效应
quietly xtreg ln_tech_density L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe
scalar a1_new = _b[L3_pub_emp]
scalar a2_new = _b[broad_tax_burden]

quietly xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    ln_tech_density ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe
scalar b_new = _b[ln_tech_density]

di "============================================"
di "【H2间接效应结果】"
di "a1（L3_pub_emp → 科技密度）       = " a1_new
di "a2（broad_tax  → 科技密度）       = " a2_new
di "b （科技密度   → ln_TFP）         = " b_new
di "pub_emp 间接效应 (a1×b) = " a1_new * b_new
di "tax 间接效应    (a2×b) = " a2_new * b_new
di "（负号=对TFP的负向间接效应）"
di "============================================"

* -------------------------------------------------------
* 原中介（专利）与新中介（科技人员密度）对比
* -------------------------------------------------------
xtreg ln_patent_apps L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med2_patent_lag

xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    ln_patent_apps ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med3_patent_lag

esttab med2_patent_lag med3_patent_lag med2_tech med3_tech ///
    using "mediation_comparison.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表3附：中介变量选择对比——专利申请量 vs 科技人员密度") ///
    mtitles("专利-第二步" "专利-第三步" "密度-第二步" "密度-第三步")

* -------------------------------------------------------
* 工具变量回归
* -------------------------------------------------------
gen iv_pub_l3   = L.L3_pub_emp
gen iv_broad_tx = L.broad_tax_burden
label var iv_pub_l3   "L3_pub_emp一阶滞后（IV）"
label var iv_broad_tx "broad_tax_burden一阶滞后（IV）"

ivregress 2sls ln_tfp L3_aging urban_rate sec_ind_share trade_openness ///
    i.year i.province_id ///
    (L3_pub_emp broad_tax_burden = iv_pub_l3 iv_broad_tx) ///
    if sample_lag == 1, vce(cluster province_id)
est store iv_revised

estat firststage

* -------------------------------------------------------
* 异质性分析
* -------------------------------------------------------
gen region = .
replace region = 1 if inlist(province_id, 1,2,3,6,9,10,11,13,15,19,21)
replace region = 2 if inlist(province_id, 4,7,8,12,14,16,17,18)
replace region = 3 if inlist(province_id, 5,20,22,23,24,25,26,27,28,29,30,31)
label define region_lbl 1 "东部" 2 "中部" 3 "西部"
label values region region_lbl

forvalues r = 1/3 {
    xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
        urban_rate sec_ind_share trade_openness ///
        i.year if sample_lag == 1 & region == `r', ///
        fe vce(cluster province_id)
    est store region_r`r'_lag
}

esttab region_r1_lag region_r2_lag region_r3_lag ///
    using "heterogeneity_revised.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表5（修订）异质性分析：分区域回归（2011-2017，L3滞后变量）") ///
    mtitles("东部（11省）" "中部（8省）" "西部（12省）")

* -------------------------------------------------------
* 稳健性检验
* -------------------------------------------------------

* 稳健性1：α=0.4
xtreg ln_tfp_alt L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store robust1_rev

* 稳健性2：剔除直辖市
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1 & !inlist(province_id,1,2,9,22), ///
    fe vce(cluster province_id)
est store robust2_rev

* 稳健性3：狭义税负对比
xtreg ln_tfp L3_pub_emp macro_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store robust3_rev

* 稳健性4：滞后5期
gen L5_pub_emp = L5.pub_emp_share2
gen L5_aging   = L5.aging_rate
gen sample_l5  = (year >= 2013 & year <= 2017) ///
    & L5_pub_emp != . & L5_aging != . ///
    & ln_tfp != . & capital != .

xtreg ln_tfp L5_pub_emp broad_tax_burden L5_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_l5 == 1, fe vce(cluster province_id)
est store robust4_l5

esttab robust1_rev robust2_rev robust3_rev robust4_l5 ///
    using "robustness_revised.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表7（修订）稳健性检验（因变量：ln_TFP）") ///
    mtitles("α=0.4替换" "剔除直辖市" "狭义税负" "滞后5期")

di "======================================"
di "全部分析完成"
di "======================================"

log close
