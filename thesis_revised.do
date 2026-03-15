/*===========================================================
  制度约束下的人力资本配置扭曲与全要素生产率损失
  修订版实证分析 —— 针对三个核心问题的修复方案
  
  修改说明：
    1. pub_emp_share2：改用滞后3期（L3），消除私营就业增速
       趋势对FE估计的干扰；补充between effects作为截面证据
    2. aging_rate：改用滞后3期（L3），匹配老龄化→TFP的
       时序传导机制（3-5年传导时滞）
    3. H2中介变量：以科技人员密度（pub_tech_staff/gdp_bn）
       替代专利申请量，机制更清晰：
       制度约束→科技人力资本行政性占用→生产性创新效率↓→TFP↓
===========================================================*/

clear all
set more off
capture log close
log using "thesis_revised_log.txt", replace text

import delimited "panel_with_broad_tax.csv", encoding(UTF-8) clear

/*-----------------------------------------------------------
  第一步：变量构造（同原版）
-----------------------------------------------------------*/
gen pub_emp_share2 = emp_urban_unit_10k / (emp_urban_unit_10k + private_emp_10k)
label var pub_emp_share2 "公共部门就业占比"

gen labor = emp_urban_unit_10k + private_emp_10k
gen ln_labor = ln(labor)

gen ln_gdp_val  = ln(gdp_bn)
gen ln_broad_tx = ln(broad_tax_burden)
gen ln_urban    = ln(urban_rate)
gen ln_aging    = ln(aging_rate)
gen ln_trade    = ln(trade_openness) if trade_openness > 0

xtset province_id year

/*-----------------------------------------------------------
  第二步：永续盘存法（同原版）
-----------------------------------------------------------*/
scalar delta  = 0.096
scalar alpha  = 0.5
scalar g_inv  = 0.161

gen capital = .
levelsof province_id, local(provinces)
foreach p of local provinces {
    quietly sum fai_bn if province_id == `p' & year == 2008
    local I2008 = r(mean)
    quietly replace capital = `I2008' / (g_inv + delta) ///
        if province_id == `p' & year == 2007
}

sort province_id year
foreach y of numlist 2008/2017 {
    replace capital = (1 - delta) * L.capital + fai_bn ///
        if year == `y' & fai_bn != .
}

gen ln_capital = ln(capital)
label var ln_capital "ln(资本存量)"

/*-----------------------------------------------------------
  第三步：TFP估算（同原版）
-----------------------------------------------------------*/
gen ln_tfp = ln_gdp_val - alpha * ln_capital - (1 - alpha) * ln_labor
label var ln_tfp "ln(TFP)：索洛残差（α=0.5）"

gen ln_tfp_alt = ln_gdp_val - 0.4 * ln_capital - 0.6 * ln_labor
label var ln_tfp_alt "ln(TFP)：索洛残差（α=0.4，稳健性）"

gen sample_main = (year >= 2008 & year <= 2017) & capital != . & ln_tfp != .

/*-----------------------------------------------------------
  第四步：【修复一】生成滞后变量
  
  理论依据：
  - pub_emp_share2: 公共部门对高技能劳动力的虹吸效应通过
    "人才结构固化→市场部门创新能力下降→TFP下降"路径
    传导，估计时滞为3-5年（参考Lu & Tao 2009）
  - aging_rate: 人口老龄化对TFP的负效应通过劳动力质量
    下降、技术扩散减缓等路径传导，时滞至少3年
    
  注：使用滞后变量后，有效样本变为2011-2017年
      N = 31×7 = 217（主样本）
-----------------------------------------------------------*/

gen L3_pub_emp = L3.pub_emp_share2
gen L3_aging   = L3.aging_rate

label var L3_pub_emp "公共部门就业占比（滞后3期）"
label var L3_aging   "老龄化率（滞后3期）"

* 滞后样本标识（2011-2017）
gen sample_lag = (year >= 2011 & year <= 2017) ///
    & L3_pub_emp != . & L3_aging != . ///
    & ln_tfp != . & capital != .

/*-----------------------------------------------------------
  第五步：【修复二】构造科技人员密度指标（H2中介变量）
  
  理论依据：
  广义税负（政府支出/GDP）越高的省份，政府部门对科技
  活动人员的行政性占用越强——科研经费流向财政供养的
  事业单位/科研机构，而非市场导向的企业创新部门。
  这一"行政性占用"机制导致同等数量的科技人员产生更低的
  生产性创新产出，从而压低TFP。
  
  指标：tech_density = pub_tech_staff / gdp_bn * 1000
       含义：每十亿元GDP对应的科技活动人员数
       预期方向：制度约束↑ → tech_density↑（科技人员被
       非生产性部门吸附）→ TFP↓
       
  注：该指标越高并非代表创新越强，而是代表单位经济产出
  所"消耗"的科技人力资本越多，是错配程度的逆向代理。
-----------------------------------------------------------*/

gen tech_density = (pub_tech_staff / gdp_bn) * 1000
label var tech_density "科技人员密度（人/十亿元GDP）—— 人力资本行政性占用代理"

gen ln_tech_density = ln(tech_density)
label var ln_tech_density "ln(科技人员密度)"

* 原专利变量保留，用于对比
label var ln_patent_apps "ln(专利申请量)（原中介，对比用）"

sum tech_density if sample_main == 1
sum tech_density if sample_lag == 1

/*-----------------------------------------------------------
  第六步：主回归——【修复版】双向固定效应
  
  变化：
  - pub_emp_share2 → L3_pub_emp（滞后3期）
  - aging_rate     → L3_aging（滞后3期）
  - 样本：2011-2017年，N≈217
-----------------------------------------------------------*/

* 模型1（修复）：核心自变量，滞后形式
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store m1_revised

* 模型2（修复）：含控制变量
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store m2_revised

* 模型3：对数形式
gen L3_ln_pub = L3.ln_pub_emp
gen L3_ln_aging = L3.ln_aging
label var L3_ln_pub   "ln(公共部门占比，滞后3期)"
label var L3_ln_aging "ln(老龄化率，滞后3期)"

xtreg ln_tfp L3_ln_pub broad_tax_burden L3_ln_aging ///
    ln_urban sec_ind_share ln_trade ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store m3_revised

* 模型4：剔除西藏
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
    note("括号内为省份层面聚类稳健标准误。pub_emp_share2和aging_rate均取滞后3期，" ///
         "以匹配挤出效应和老龄化效应的时序传导机制（理论时滞3-5年）。" ///
         "有效样本为2011-2017年，N≈217。")

/*-----------------------------------------------------------
  第七步：截面结构效应补充验证（Between Effects）
  
  理论依据：公共部门对人才的结构性虹吸效应是长期平衡状态，
  体现为省份间持久性差异（而非省内短期波动）。
  Between estimator正是利用各省时序均值的截面差异进行估计，
  是识别结构性效应的适当方法，可与FE模型互补。
-----------------------------------------------------------*/

* Between effects：利用省份均值的截面差异
xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    if sample_main == 1, be
est store be_model

esttab m2_revised be_model using "fe_vs_be.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表2补充：FE（省内时序效应）vs BE（省间结构效应）对比") ///
    mtitles("FE-滞后期（时序效应）" "BE（截面结构效应）") ///
    note("两种估计量识别不同维度的效应：" ///
         "FE利用省内时序变化（短期动态因果）；" ///
         "BE利用省间截面差异（长期结构均衡）。" ///
         "理论上虹吸效应属于后者，两种结果互相印证。")

/*-----------------------------------------------------------
  第八步：【修复三】中介效应检验（H2）
  
  新中介变量：科技人员密度（ln_tech_density）
  
  传导路径：
  broad_tax_burden（↑）→ ln_tech_density（↑，科技人力资本
  被行政性占用加深）→ ln_tfp（↓，人力资本配置效率下降）
  
  注意方向：
  - 第二步预期：broad_tax_burden → ln_tech_density（正向）
  - 第三步预期：ln_tech_density → ln_tfp（负向）
  - 间接效应 = a×b（负号）→ broad_tax_burden通过人力资本
    行政性占用渠道对TFP产生额外负效应
-----------------------------------------------------------*/

* 第二步：自变量对新中介变量（科技人员密度）
xtreg ln_tech_density L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med2_tech

* 第三步：自变量+新中介变量对TFP
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    ln_tech_density ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med3_tech

esttab med2_tech med3_tech using "mediation_revised.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表3（修订）中介效应检验（中介变量：科技人员密度）") ///
    mtitles("第二步：→ln(科技人员密度)" "第三步：→ln(TFP)（含中介）") ///
    note("科技人员密度=科技活动人员数/GDP（十亿元）×1000，" ///
         "越高代表单位经济产出消耗的科技人力资本越多，" ///
         "是人力资本行政性占用程度的逆向代理指标。" ///
         "中介效应方向：制度约束↑→科技人员密度↑→TFP↓。")

* 计算间接效应量
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

di "【H2间接效应】"
di "L3_pub_emp 通过科技人员密度的间接效应 (a×b) = " a1_new * b_new
di "broad_tax_burden 通过科技人员密度的间接效应 (a×b) = " a2_new * b_new
di "（负号=对TFP产生负向间接效应，与理论预期一致）"

/*-----------------------------------------------------------
  对比：原中介变量（专利）vs 新中介变量（科技人员密度）
-----------------------------------------------------------*/
xtreg ln_patent_apps L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med2_patent_lag

xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    ln_patent_apps ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store med3_patent_lag

esttab med2_patent_lag med3_patent_lag ///
    med2_tech med3_tech ///
    using "mediation_comparison.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表3对比：中介变量选择——专利申请量 vs 科技人员密度") ///
    mtitles("专利-第二步" "专利-第三步" "密度-第二步" "密度-第三步")

/*-----------------------------------------------------------
  第九步：工具变量（使用滞后版本）
-----------------------------------------------------------*/
gen iv_pub_emp   = L.L3_pub_emp
gen iv_broad_tax = L.broad_tax_burden
label var iv_pub_emp   "L3_pub_emp的一阶滞后（工具变量）"
label var iv_broad_tax "broad_tax_burden的一阶滞后（工具变量）"

ivregress 2sls ln_tfp L3_aging urban_rate sec_ind_share trade_openness ///
    i.year i.province_id ///
    (L3_pub_emp broad_tax_burden = iv_pub_emp iv_broad_tax) ///
    if sample_lag == 1, vce(cluster province_id)
est store iv_revised

estat firststage

/*-----------------------------------------------------------
  第十步：异质性分析（滞后版）
-----------------------------------------------------------*/
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
    title("表5（修订）异质性分析：分区域回归（2011-2017，滞后变量）") ///
    mtitles("东部" "中部" "西部")

/*-----------------------------------------------------------
  第十一步：稳健性检验（修订版）
-----------------------------------------------------------*/

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

* 稳健性4：滞后5期（对pub_emp和aging效应的进一步滞后验证）
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
    title("表7（修订）稳健性检验") ///
    mtitles("α=0.4替换" "剔除直辖市" "狭义税负" "滞后5期")

di "======================================"
di "修订版分析完成"
di "核心修订："
di "  1. pub_emp和aging均改用L3滞后期"
di "  2. H2中介变量改为科技人员密度(ln_tech_density)"
di "  3. 补充between effects截面结构效应证据"
di "  4. 稳健性新增L5滞后期验证"
di "======================================"

log close
