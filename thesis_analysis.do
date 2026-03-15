/*===========================================================
  制度约束下的人力资本配置扭曲与全要素生产率损失
  实证分析主do文件 —— 基于索洛残差法估算TFP
  样本：31省，2008-2017年，N=310（平衡面板）
  核心更新：
    1. 永续盘存法估算省级资本存量
    2. 索洛残差法构建真实TFP
    3. 广义税负口径（政府支出/GDP）
===========================================================*/

clear all
set more off
capture log close
log using "thesis_log.txt", replace text

/*-----------------------------------------------------------
  第一步：数据导入与基础整理
-----------------------------------------------------------*/

import delimited "panel_with_broad_tax.csv", encoding(UTF-8) clear

* 生成核心变量
gen pub_emp_share2 = emp_urban_unit_10k / (emp_urban_unit_10k + private_emp_10k)
label var pub_emp_share2   "公共部门就业占比（代理）"
label var broad_tax_burden "广义税负（政府支出/GDP，%）"

* 生成对数变量
gen ln_pub_emp  = ln(pub_emp_share2)
gen ln_aging    = ln(aging_rate)
gen ln_urban    = ln(urban_rate)
gen ln_trade    = ln(trade_openness) if trade_openness > 0
gen ln_broad_tx = ln(broad_tax_burden)

* 劳动投入：城镇单位+私营个体就业（万人）
gen labor = emp_urban_unit_10k + private_emp_10k
label var labor "劳动投入（城镇单位+私营个体，万人）"
gen ln_labor = ln(labor)

* 设定面板结构
xtset province_id year

/*-----------------------------------------------------------
  第二步：永续盘存法估算省级资本存量
  方法：K_t = (1-δ)*K_{t-1} + I_t
  参数：δ=0.096（折旧率，参考张军等2004）
        α=0.5（资本产出弹性，参考中国文献惯例）
  初始资本存量：K_2007 = I_2008 / (g + δ)
        g=0.161（2008-2017年固定资产投资平均增速）
  注：fai_bn单位为亿元，gdp_bn单位为亿元，保持一致
-----------------------------------------------------------*/

* 折旧率和资本弹性
scalar delta = 0.096
scalar alpha = 0.5
scalar g_inv = 0.161

* 初始资本存量（2007年）= 2008年投资 / (平均增速 + 折旧率)
* 按省份生成初始值
gen capital = .
label var capital "省级资本存量（永续盘存法，亿元）"

* 设定2007年初始资本存量
bysort province_id (year): replace capital = fai_bn[2] / (g_inv + delta) ///
    if year == 2007
* 注：fai_bn[2]是2008年的值（数据从2005年开始，2007年行是第3行，但fai从2008起有数）
* 更稳健的写法：
drop capital
gen capital = .

* 用2008年投资初始化2007年资本存量
levelsof province_id, local(provinces)
foreach p of local provinces {
    quietly sum fai_bn if province_id == `p' & year == 2008
    local I2008 = r(mean)
    quietly replace capital = `I2008' / (g_inv + delta) ///
        if province_id == `p' & year == 2007
}

* 永续盘存递推：2008年起
sort province_id year
foreach y of numlist 2008/2017 {
    replace capital = (1 - delta) * L.capital + fai_bn ///
        if year == `y' & fai_bn != .
}

* 检查资本存量
sum capital if year >= 2008 & year <= 2017
gen ln_capital = ln(capital)
label var ln_capital "ln(资本存量)"

/*-----------------------------------------------------------
  第三步：索洛残差法估算TFP
  ln(TFP) = ln(GDP) - α*ln(K) - (1-α)*ln(L)
  α = 0.5（资本产出弹性）
-----------------------------------------------------------*/

gen ln_gdp_val = ln(gdp_bn)

gen ln_tfp = ln_gdp_val - alpha * ln_capital - (1 - alpha) * ln_labor
label var ln_tfp "ln(TFP)：索洛残差（α=0.5）"

* 稳健性备用：α=0.4
gen ln_tfp_alt = ln_gdp_val - 0.4 * ln_capital - 0.6 * ln_labor
label var ln_tfp_alt "ln(TFP)：索洛残差（α=0.4，稳健性）"

* 限定主回归样本
gen sample_main = (year >= 2008 & year <= 2017) & capital != . & ln_tfp != .

/*-----------------------------------------------------------
  第四步：描述性统计
-----------------------------------------------------------*/

estpost summarize ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    ln_patent_apps urban_rate sec_ind_share trade_openness ///
    ln_capital ln_labor ///
    if sample_main == 1, detail

esttab using "desc_stats.rtf", replace ///
    cells("mean(fmt(3)) sd(fmt(3)) min(fmt(3)) max(fmt(3)) count(fmt(0))") ///
    title("表1 描述性统计（2008-2017年，N=310）") ///
    note("注：TFP采用索洛残差法估算，资本存量采用永续盘存法（δ=0.096，α=0.5）。广义税负为一般公共预算支出/GDP，遵循奥地利学派政府支出即税收口径。") ///
    nonumber nomtitle

* 相关系数矩阵
pwcorr ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    ln_patent_apps if sample_main == 1, star(0.05)

/*-----------------------------------------------------------
  第五步：主回归——双向固定效应
  因变量：ln_tfp（索洛残差TFP）
  H1：pub_emp_share2、broad_tax_burden、aging_rate
      与ln_tfp显著负相关
  注：不纳入ln_gdp作控制变量（TFP由GDP构造，会产生
      机械相关）；改用urban_rate、sec_ind_share、
      trade_openness作为控制变量
-----------------------------------------------------------*/

* 模型1：仅核心自变量
xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store m1

* 模型2：加入控制变量
xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store m2

* 模型3：对数形式
xtreg ln_tfp ln_pub_emp broad_tax_burden ln_aging ///
    ln_urban sec_ind_share ln_trade ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store m3

* 模型4：剔除西藏（广义税负异常省份）
xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1 & province_id != 26, fe vce(cluster province_id)
est store m4

* 输出主回归结果
esttab m1 m2 m3 m4 using "main_regression.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) ///
    b(3) se(3) ///
    title("表2 双向固定效应主回归结果（因变量：ln_TFP）") ///
    mtitles("模型1" "模型2" "模型3（对数）" "模型4（剔除西藏）") ///
    note("括号内为聚类标准误（省份层面）。所有模型包含省份固定效应和年份固定效应。TFP为索洛残差法估算（α=0.5，δ=0.096）。广义税负为一般公共预算支出/GDP。")

/*-----------------------------------------------------------
  第六步：中介效应检验
  H2：制度约束通过抑制创新（ln_patent_apps）损害TFP
  路径：pub_emp_share2 → ln_patent_apps → ln_tfp
        broad_tax_burden → ln_patent_apps → ln_tfp
-----------------------------------------------------------*/

* 第二步：自变量对中介变量
xtreg ln_patent_apps pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store med_step2

* 第三步：自变量+中介变量同时对因变量
xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    ln_patent_apps ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store med_step3

esttab med_step2 med_step3 using "mediation.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表3 中介效应检验（广义税负口径，因变量：ln_TFP）") ///
    mtitles("第二步：→ln(专利申请)" "第三步：→ln(TFP)（含中介）")

* 间接效应
quietly xtreg ln_patent_apps pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe
scalar a1 = _b[pub_emp_share2]
scalar a2 = _b[broad_tax_burden]

quietly xtreg ln_tfp ln_patent_apps pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe
scalar b_med = _b[ln_patent_apps]

di "pub_emp_share2 间接效应（a×b）= " a1 * b_med
di "broad_tax_burden 间接效应（a×b）= " a2 * b_med

/*-----------------------------------------------------------
  第七步：工具变量回归（处理内生性）
  工具变量：滞后一期的pub_emp_share2和broad_tax_burden
-----------------------------------------------------------*/

xtset province_id year
gen iv_pub_emp   = L.pub_emp_share2
gen iv_broad_tax = L.broad_tax_burden
label var iv_pub_emp   "pub_emp_share2滞后一期"
label var iv_broad_tax "broad_tax_burden滞后一期"

ivregress 2sls ln_tfp aging_rate urban_rate sec_ind_share trade_openness ///
    i.year i.province_id ///
    (pub_emp_share2 broad_tax_burden = iv_pub_emp iv_broad_tax) ///
    if sample_main == 1, vce(cluster province_id)
est store iv_model

estat firststage

/*-----------------------------------------------------------
  第八步：异质性分析
  H3：挤出效应在市场化程度低的地区更强
-----------------------------------------------------------*/

gen region = .
replace region = 1 if inlist(province_id, 1,2,3,6,9,10,11,13,15,19,21)
replace region = 2 if inlist(province_id, 4,7,8,12,14,16,17,18)
replace region = 3 if inlist(province_id, 5,20,22,23,24,25,26,27,28,29,30,31)
label define region_lbl 1 "东部" 2 "中部" 3 "西部"
label values region region_lbl

forvalues r = 1/3 {
    xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
        urban_rate sec_ind_share trade_openness ///
        i.year if sample_main == 1 & region == `r', fe vce(cluster province_id)
    est store region_`r'
}

esttab region_1 region_2 region_3 using "heterogeneity.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表5 异质性分析：分区域回归（因变量：ln_TFP）") ///
    mtitles("东部" "中部" "西部")

* 交互项检验
gen pub_emp_west  = pub_emp_share2   * (region == 3)
gen broad_tx_west = broad_tax_burden * (region == 3)
label var pub_emp_west  "公共部门占比×西部"
label var broad_tx_west "广义税负×西部"

xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    pub_emp_west broad_tx_west ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store interaction

esttab interaction using "interaction.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表6 交互项检验（西部虚拟变量）")

/*-----------------------------------------------------------
  第九步：稳健性检验
-----------------------------------------------------------*/

* 稳健性1：替换α=0.4（资本弹性敏感性）
xtreg ln_tfp_alt pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store robust1

* 稳健性2：剔除直辖市
xtreg ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1 & !inlist(province_id,1,2,9,22), ///
    fe vce(cluster province_id)
est store robust2

* 稳健性3：狭义税负口径对比
xtreg ln_tfp pub_emp_share2 macro_tax_burden aging_rate ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_main == 1, fe vce(cluster province_id)
est store robust3

esttab robust1 robust2 robust3 using "robustness.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("表7 稳健性检验（因变量：ln_TFP）") ///
    mtitles("α=0.4替换" "剔除直辖市" "狭义税负对比")

/*-----------------------------------------------------------
  第十步：输出汇总
-----------------------------------------------------------*/

di "======================================"
di "所有分析完成"
di "样本：31省，2008-2017年，N=310"
di "因变量：ln_TFP（索洛残差法，α=0.5，δ=0.096）"
di "核心解释变量：pub_emp_share2，broad_tax_burden"
di "生成文件："
di "  desc_stats.rtf       描述性统计"
di "  main_regression.rtf  主回归"
di "  mediation.rtf        中介效应"
di "  iv_regression.rtf    工具变量"
di "  heterogeneity.rtf    异质性分析"
di "  interaction.rtf      交互项检验"
di "  robustness.rtf       稳健性检验"
di "======================================"

log close
