/*===========================================================
  随机前沿分析（SFA）—— TFP分解与效率损耗机制检验
  
  理论逻辑：
  索洛残差法估算的TFP是技术效率与技术进步的混合体，
  无法区分"配置扭曲导致的效率损耗"与"技术本身停滞"。
  SFA将生产前沿以下的残差分解为：
    - 技术无效率项 u_it（管理与配置效率，单侧分布）
    - 随机误差项 v_it（测量误差与随机冲击，双侧正态）
  
  如果广义税负和公共部门扩张主要压制的是技术效率分量
  （而非技术进步分量），则直接验证了"配置扭曲→效率损耗"
  的机制，比索洛残差法更精确地回应论文的核心理论命题。
  
  生产函数设定：柯布-道格拉斯（C-D）
  ln(GDP_it) = α + β₁·ln(K_it) + β₂·ln(L_it) + v_it - u_it
  
  u_it = 技术无效率项，服从截断正态分布
  效率得分 TE_it = exp(-u_it)，范围[0,1]，1为完全有效
  
  样本：31省，2008-2017年，N=310
===========================================================*/

clear all
set more off
capture log close
log using "sfa_log.txt", replace text

* ---- 安装SFA包 ----
* ssc install sfpanel, replace   /* 面板SFA */
* ssc install frontier, replace  /* 截面SFA（备用）*/

import delimited "panel_with_broad_tax.csv", encoding(UTF-8) clear

/*-----------------------------------------------------------
  第一步：重建基础变量
-----------------------------------------------------------*/
gen labor = emp_urban_unit_10k + private_emp_10k
gen pub_emp_share2 = emp_urban_unit_10k / (emp_urban_unit_10k + private_emp_10k)

scalar delta = 0.096
scalar g_inv = 0.161

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

gen ln_gdp_val = ln(gdp_bn)
gen ln_capital = ln(capital)
gen ln_labor   = ln(labor)
gen ln_gdp_pc  = ln(gdp_per_cap) if gdp_per_cap > 0

xtset province_id year

gen sample_main = (year >= 2008 & year <= 2017) ///
    & capital != . & labor != . & gdp_bn != .

/*-----------------------------------------------------------
  第二步：面板SFA估计
  使用 sfpanel 命令（Belotti et al., 2013）
  
  模型选择：时变效率模型（TFE / TRE）
  - bc92: Battese-Coelli 1992 模型（u_it = u_i * exp(η(t-T))）
    效率随时间变化，但单调，适合本文10年短面板
  - tfe:  True Fixed Effects，允许任意时变效率
  - kumb90: Kumbhakar 1990 时变模型
  
  推荐：先用 bc92，再用 tfe 作稳健性对比
-----------------------------------------------------------*/

* 模型1：Battese-Coelli 1992（时变效率，效率指数分布）
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(bc92) ///
    distribution(hn)   /* hn=半正态，tn=截断正态，e=指数 */
    
* 保存效率得分
predict te_bc92, te
label var te_bc92 "技术效率得分（BC92，半正态）"

est store sfa_bc92

di "BC92模型技术效率得分描述统计："
sum te_bc92 if sample_main == 1

* 模型2：True Fixed Effects（允许任意时变效率）
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(tfe) ///
    distribution(hn)

predict te_tfe, te
label var te_tfe "技术效率得分（TFE，半正态）"

est store sfa_tfe

di "TFE模型技术效率得分描述统计："
sum te_tfe if sample_main == 1

* 模型3：截断正态分布（稳健性）
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(bc92) ///
    distribution(tn)

predict te_tn, te
label var te_tn "技术效率得分（BC92，截断正态）"

* 相关性诊断：SFA效率得分 vs 索洛残差TFP
gen ln_tfp_solow = ln_gdp_val - 0.5 * ln_capital - 0.5 * ln_labor

pwcorr te_bc92 te_tfe ln_tfp_solow ///
    if sample_main == 1, star(0.05)

di "诊断：若SFA效率得分与索洛残差TFP高度正相关（r>0.6），"
di "说明两种测算方法捕捉的是同一底层变量，互相印证。"

/*-----------------------------------------------------------
  第三步：以技术效率得分作为因变量，重跑主回归
  
  核心问题：广义税负和公共部门扩张
  是否主要损害了技术效率（配置扭曲渠道）
  而非技术进步（技术停滞渠道）？
  
  如果税负对TE的回归系数显著为负，
  说明损耗机制主要是配置效率渠道，
  与论文"人力资本错配→效率损耗"的理论高度吻合。
-----------------------------------------------------------*/

gen L3_pub_emp = L3.pub_emp_share2
gen L3_aging   = L3.aging_rate
label var L3_pub_emp "公共部门就业占比（滞后3期）"
label var L3_aging   "老龄化率（滞后3期）"

gen sample_lag = (year >= 2011 & year <= 2017) ///
    & L3_pub_emp != . & L3_aging != . & te_bc92 != .

* 技术效率得分作为因变量
xtreg te_bc92 L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store sfa_reg_bc92

xtreg te_tfe L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store sfa_reg_tfe

/*-----------------------------------------------------------
  第四步：并行对比——索洛残差 vs SFA效率得分
  这是论文中"TFP测算方法稳健性检验"的核心表格
-----------------------------------------------------------*/

* 重建索洛残差模型（若在新会话，需先跑thesis_patch.do）
xtreg ln_tfp_solow L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store solow_check

esttab solow_check sfa_reg_bc92 sfa_reg_tfe ///
    using "robustness_sfa.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("稳健性检验：TFP测算方法对比（索洛残差 vs SFA技术效率）") ///
    mtitles("索洛残差TFP" "SFA效率得分（BC92）" "SFA效率得分（TFE）") ///
    note("括号内为省份层面聚类稳健标准误。" ///
         "SFA基于柯布-道格拉斯前沿函数（ln(GDP)=f(K,L)）估算，" ///
         "技术效率得分TE∈(0,1]，1=完全有效。" ///
         "若broad_tax_burden在三列中均显著为负，" ///
         "说明税负抑制效应不依赖TFP测算方法的参数假设；" ///
         "若SFA效率得分（配置效率分量）的系数绝对值更大，" ///
         "则进一步验证损耗主要来自配置扭曲渠道而非技术停滞。")

/*-----------------------------------------------------------
  第五步：效率得分的空间分布可视化（描述性）
  按省份均值排名，识别高效率vs低效率省份格局
-----------------------------------------------------------*/

bysort province_id: egen mean_te = mean(te_bc92) ///
    if sample_main == 1
bysort province_id: egen mean_tax = mean(broad_tax_burden) ///
    if sample_main == 1

* 低效率省份是否系统性地具有更高税负？
pwcorr mean_te mean_tax if year == 2008, star(0.05)

di "截面相关性诊断：省份平均效率得分 vs 平均广义税负"
di "预期符号：负相关（高税负省份效率更低）"

/*-----------------------------------------------------------
  注意事项：
  1. sfpanel需要Stata 13+，若报错请更新
  2. 若sfpanel未安装，先运行 ssc install sfpanel
  3. distribution()选项：hn（半正态）最保守，
     tn（截断正态）允许非零众数，e（指数）最灵活
     三种分布假设的结果应方向一致，否则需讨论
  4. 在论文中：
     - SFA部分放入第六章6.2.4稳健性检验
     - 核心展示：broad_tax_burden在三种TFP下均显著
     - 进一步论述：SFA可分离配置效率分量，
       税负对TE的负效应证实损耗来自配置扭曲渠道
===========================================================*/

log close
