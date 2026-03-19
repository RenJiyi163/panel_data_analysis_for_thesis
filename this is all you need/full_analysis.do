/*===========================================================
  完整实证分析主文件
  从CSV从头开始，不依赖任何先前状态
  
  包含：
  1. 主回归（TWFE）
  2. 中介效应
  3. 工具变量
  4. 异质性分析
  5. 稳健性检验（含α替换、剔除直辖市、狭义税负、L5）
  6. DEA稳健性（需Python：shell python dea_solve.py）
  7. SFA随机前沿分析
  8. 空间杜宾模型（SDM）
===========================================================*/

clear all
set more off
capture log close
log using "full_analysis_log.txt", replace text

/*-----------------------------------------------------------
  第一步：数据构建（完整版，一次性生成所有变量）
-----------------------------------------------------------*/
import delimited "panel_with_broad_tax.csv", encoding(UTF-8) clear

* 面板设定（必须最先做）
xtset province_id year

* 基础变量
gen pub_emp_share2 = emp_urban_unit_10k / (emp_urban_unit_10k + private_emp_10k)
gen labor          = emp_urban_unit_10k + private_emp_10k

* 永续盘存法（xtset在前，L.可用）
scalar delta = 0.096
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

* 对数变量
gen ln_gdp_val      = ln(gdp_bn)
gen ln_capital      = ln(capital)
gen ln_labor        = ln(labor)
gen ln_pub_emp      = ln(pub_emp_share2)
gen ln_aging        = ln(aging_rate)
gen ln_urban        = ln(urban_rate)
gen ln_trade        = ln(trade_openness) if trade_openness > 0
gen ln_tech_density = ln((pub_tech_staff / gdp_bn) * 1000)

* TFP
gen ln_tfp     = ln_gdp_val - 0.5*ln_capital - 0.5*ln_labor
gen ln_tfp_alt = ln_gdp_val - 0.4*ln_capital - 0.6*ln_labor

* 滞后变量
gen L3_pub_emp = L3.pub_emp_share2
gen L5_pub_emp = L5.pub_emp_share2
gen L3_aging   = L3.aging_rate
gen L5_aging   = L5.aging_rate

* 区域分组
gen region = .
replace region = 1 if inlist(province_id,1,2,3,6,9,10,11,13,15,19,21)
replace region = 2 if inlist(province_id,4,7,8,12,14,16,17,18)
replace region = 3 if inlist(province_id,5,20,22,23,24,25,26,27,28,29,30,31)
label define region_lbl 1 "东部" 2 "中部" 3 "西部"
label values region region_lbl

* 样本标识
* main：2008-2017，有资本和TFP
gen sample_main = (year>=2008 & year<=2017) & capital!=. & ln_tfp!=.
* lag：2011-2017，L3变量有值，不要求trade非缺失（让各回归自行处理）
gen sample_lag  = (year>=2011 & year<=2017) & L3_pub_emp!=. & L3_aging!=. & ln_tfp!=.
* l5：2013-2017，L5变量有值
gen sample_l5   = (year>=2013 & year<=2017) & L5_pub_emp!=. & L5_aging!=. & ln_tfp!=.
* SDM专用：严格平衡面板，31省×7年=217，不含trade_openness
* 检查sample_lag下ln_tfp等核心变量是否完整
gen sample_sdm  = sample_lag & ln_tfp!=. & L3_pub_emp!=. & broad_tax_burden!=. & L3_aging!=. & urban_rate!=. & sec_ind_share!=.

* 确认SDM样本是否严格平衡
quietly tab province_id if sample_sdm==1
di "SDM样本省份数：" r(r)
quietly tab year if sample_sdm==1
di "SDM样本年份数：" r(r)
count if sample_sdm==1
di "SDM样本总观测：" r(N) "（应为217=31×7）"

/*-----------------------------------------------------------
  第二步：描述性统计
-----------------------------------------------------------*/
estpost summarize ln_tfp pub_emp_share2 broad_tax_burden aging_rate ///
    ln_tech_density urban_rate sec_ind_share trade_openness ///
    ln_capital ln_labor if sample_main==1, detail

esttab using "desc_stats.rtf", replace ///
    cells("mean(fmt(3)) sd(fmt(3)) min(fmt(3)) max(fmt(3)) count(fmt(0))") ///
    title("描述性统计（2008-2017年，N=310）") nonumber nomtitle

/*-----------------------------------------------------------
  第三步：主回归——双向固定效应
-----------------------------------------------------------*/
* 模型1：核心变量
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store m1

* 模型2：含控制变量
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store m2

* 模型3：对数形式
xtreg ln_tfp ln_pub_emp broad_tax_burden ln_aging ///
    ln_urban sec_ind_share ln_trade ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store m3

* 模型4：剔除西藏
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1 & province_id!=26, fe vce(cluster province_id)
est store m4

esttab m1 m2 m3 m4 using "main_regression.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("主回归：双向固定效应（因变量：ln_TFP，2011-2017年）") ///
    mtitles("模型1" "模型2（含控制）" "模型3（对数）" "模型4（剔除西藏）") ///
    note("括号内为省份层面聚类稳健标准误。pub_emp和aging取滞后3期（L3）。")

/*-----------------------------------------------------------
  第四步：中介效应（Baron-Kenny三步法）
-----------------------------------------------------------*/
* 第二步：自变量→中介变量
xtreg ln_tech_density L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store med_step2

* 第三步：含中介变量→TFP
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    ln_tech_density ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store med_step3

esttab m2 med_step2 med_step3 using "mediation.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("中介效应：制度约束→科技人力资本行政性占用→TFP") ///
    mtitles("第一步：TFP总效应" "第二步：→ln(科技人员密度)" "第三步：→ln(TFP)含中介")

* 间接效应量
quietly xtreg ln_tech_density L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness i.year if sample_lag==1, fe
scalar a_tax = _b[broad_tax_burden]
scalar a_pub = _b[L3_pub_emp]
quietly xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    ln_tech_density urban_rate sec_ind_share trade_openness i.year if sample_lag==1, fe
scalar b_med = _b[ln_tech_density]
di "间接效应：broad_tax = " a_tax*b_med "  pub_emp = " a_pub*b_med

/*-----------------------------------------------------------
  第五步：工具变量（2SLS）
-----------------------------------------------------------*/
gen iv_pub_emp   = L.L3_pub_emp
gen iv_broad_tax = L.broad_tax_burden
label var iv_pub_emp   "L3_pub_emp一阶滞后（IV）"
label var iv_broad_tax "broad_tax_burden一阶滞后（IV）"

ivregress 2sls ln_tfp L3_aging urban_rate sec_ind_share trade_openness ///
    i.year i.province_id ///
    (L3_pub_emp broad_tax_burden = iv_pub_emp iv_broad_tax) ///
    if sample_lag==1, vce(cluster province_id)
est store iv_model
estat firststage

/*-----------------------------------------------------------
  第六步：异质性分析
-----------------------------------------------------------*/
forvalues r = 1/3 {
    xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
        urban_rate sec_ind_share trade_openness ///
        i.year if sample_lag==1 & region==`r', fe vce(cluster province_id)
    est store region_`r'
}

esttab region_1 region_2 region_3 using "heterogeneity.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("异质性分析：分区域回归（因变量：ln_TFP）") ///
    mtitles("东部" "中部" "西部")

/*-----------------------------------------------------------
  第七步：稳健性检验（原有四组）
-----------------------------------------------------------*/
* 稳1：α=0.4
xtreg ln_tfp_alt L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store robust1

* 稳2：剔除直辖市
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1 & !inlist(province_id,1,2,9,22), ///
    fe vce(cluster province_id)
est store robust2

* 稳3：狭义税负
xtreg ln_tfp L3_pub_emp macro_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store robust3

* 稳4：滞后5期
xtreg ln_tfp L5_pub_emp broad_tax_burden L5_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_l5==1, fe vce(cluster province_id)
est store robust4

esttab robust1 robust2 robust3 robust4 using "robustness.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("稳健性检验（因变量：ln_TFP）") ///
    mtitles("α=0.4替换" "剔除直辖市" "狭义税负" "滞后5期")

/*-----------------------------------------------------------
  第八步：DEA稳健性
  需要：dea_efficiency.csv（由dea_solve.py生成）
  如未生成：先运行 shell python dea_solve.py
-----------------------------------------------------------*/
* 导出DEA投入产出数据
preserve
    keep if sample_main==1
    keep province_id year capital labor gdp_bn
    sort year province_id
    export delimited "dea_input.csv", replace
restore

* 调用Python计算DEA
shell python dea_solve.py

* 导入效率得分
capture confirm file "dea_efficiency.csv"
if _rc == 0 {
    preserve
        import delimited "dea_efficiency.csv", clear varnames(1)
        save "dea_efficiency.dta", replace
    restore
    merge m:1 province_id year using "dea_efficiency.dta", ///
        keepusing(dea_eff) nogenerate
    gen ln_dea_eff = ln(dea_eff)
    label var dea_eff    "DEA技术效率得分（VRS）"
    label var ln_dea_eff "ln(DEA效率得分)"

    pwcorr ln_tfp ln_dea_eff if sample_main==1, star(0.05)

    xtreg ln_dea_eff L3_pub_emp broad_tax_burden L3_aging ///
        urban_rate sec_ind_share trade_openness ///
        i.year if sample_lag==1, fe vce(cluster province_id)
    est store dea_reg

    esttab m2 dea_reg using "robustness_dea.rtf", replace ///
        star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
        title("稳健性：索洛残差 vs DEA效率得分") ///
        mtitles("索洛残差TFP" "DEA效率得分（VRS）")
}
else {
    di as err "dea_efficiency.csv未找到，跳过DEA模块"
    di as err "请确认 dea_solve.py 在同一文件夹，且Python已安装scipy"
}

/*-----------------------------------------------------------
  第九步：SFA随机前沿分析
-----------------------------------------------------------*/
capture which sfpanel
if _rc != 0 {
    di "安装 sfpanel..."
    ssc install sfpanel
}

* BC92：时变效率，截断正态
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main==1, model(bc92) distribution(tnormal)
predict te_bc92, bc
label var te_bc92 "SFA技术效率（BC92，截断正态）"
est store sfa_bc92

* TFE：True Fixed Effects，半正态
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main==1, model(tfe) distribution(hnormal)
predict te_tfe, bc
label var te_tfe "SFA技术效率（TFE，半正态）"
est store sfa_tfe

sum te_bc92 te_tfe if sample_main==1
pwcorr te_bc92 te_tfe ln_tfp if sample_main==1, star(0.05)

xtreg te_bc92 L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store sfa_reg_bc92

xtreg te_tfe L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag==1, fe vce(cluster province_id)
est store sfa_reg_tfe

esttab m2 sfa_reg_bc92 sfa_reg_tfe using "robustness_sfa.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("稳健性：TFP测算方法对比（索洛残差 vs SFA技术效率）") ///
    mtitles("索洛残差" "SFA-BC92" "SFA-TFE") ///
    note("SFA基于C-D前沿函数。BC92=时变效率截断正态；TFE=True FE半正态。")

/*-----------------------------------------------------------
  第十步：空间杜宾模型（SDM）
  
  关键：xsmle要求严格平衡面板
  解法：不含trade_openness（它有缺失），用sample_sdm
        sample_sdm = sample_lag且核心变量无缺失 = 217观测值
-----------------------------------------------------------*/

* ---- 构建经济距离权重矩阵 ----
bysort province_id: egen mean_gdppc = mean(gdp_per_cap) ///
    if sample_main==1 & gdp_per_cap > 0

preserve
    keep if sample_main==1 & year==2010
    keep province_id mean_gdppc
    sort province_id
    quietly sum mean_gdppc
    replace mean_gdppc = r(mean) if mean_gdppc==. | mean_gdppc==0
    mkmat mean_gdppc, matrix(GDPPC)
restore

mata:
function spatial_lag(string scalar y, string scalar w, string scalar out) {
    st_matrix(out, st_matrix(w) * st_matrix(y))
}

gdppc = st_matrix("GDPPC")
n = rows(gdppc)
W_e = J(n, n, 0)
for (i=1; i<=n; i++) {
    for (j=1; j<=n; j++) {
        if (i != j) {
            d = abs(gdppc[i,1] - gdppc[j,1])
            if (d > 0) W_e[i,j] = 1/d
        }
    }
}
for (i=1; i<=n; i++) {
    rs = sum(W_e[i,.])
    if (rs > 0) W_e[i,.] = W_e[i,.] / rs
}
printf("经济距离矩阵行和均值: %f\n", mean(rowsum(W_e)))
st_matrix("W_econ_raw", W_e)
end

* ---- 构建地理邻接权重矩阵 ----
mata:
n = 31
A = J(n, n, 0)
A[1,2]=1;  A[1,3]=1
A[2,1]=1;  A[2,3]=1
A[3,1]=1;  A[3,2]=1;  A[3,4]=1;  A[3,5]=1;  A[3,6]=1;  A[3,15]=1; A[3,16]=1
A[4,3]=1;  A[4,5]=1;  A[4,16]=1; A[4,27]=1
A[5,3]=1;  A[5,4]=1;  A[5,6]=1;  A[5,7]=1;  A[5,8]=1;  A[5,27]=1; A[5,28]=1; A[5,30]=1
A[6,3]=1;  A[6,5]=1;  A[6,7]=1
A[7,5]=1;  A[7,6]=1;  A[7,8]=1
A[8,5]=1;  A[8,7]=1
A[9,10]=1; A[9,11]=1
A[10,9]=1; A[10,11]=1; A[10,12]=1; A[10,15]=1
A[11,9]=1; A[11,10]=1; A[11,12]=1; A[11,13]=1; A[11,14]=1
A[12,10]=1; A[12,11]=1; A[12,14]=1; A[12,15]=1; A[12,16]=1; A[12,17]=1
A[13,11]=1; A[13,14]=1; A[13,19]=1
A[14,11]=1; A[14,12]=1; A[14,13]=1; A[14,17]=1; A[14,18]=1; A[14,19]=1
A[15,3]=1; A[15,10]=1; A[15,12]=1; A[15,16]=1
A[16,3]=1; A[16,4]=1; A[16,12]=1; A[16,15]=1; A[16,17]=1; A[16,27]=1
A[17,12]=1; A[17,14]=1; A[17,16]=1; A[17,18]=1; A[17,22]=1; A[17,27]=1
A[18,14]=1; A[18,17]=1; A[18,19]=1; A[18,20]=1; A[18,22]=1; A[18,24]=1
A[19,13]=1; A[19,14]=1; A[19,18]=1; A[19,20]=1
A[20,18]=1; A[20,19]=1; A[20,24]=1; A[20,25]=1
A[22,17]=1; A[22,18]=1; A[22,23]=1; A[22,24]=1; A[22,27]=1
A[23,22]=1; A[23,24]=1; A[23,25]=1; A[23,26]=1; A[23,27]=1; A[23,28]=1; A[23,29]=1
A[24,18]=1; A[24,20]=1; A[24,22]=1; A[24,23]=1; A[24,25]=1
A[25,20]=1; A[25,23]=1; A[25,24]=1; A[25,26]=1
A[26,23]=1; A[26,25]=1; A[26,29]=1; A[26,31]=1
A[27,4]=1; A[27,5]=1; A[27,16]=1; A[27,17]=1; A[27,22]=1; A[27,23]=1; A[27,28]=1; A[27,30]=1
A[28,5]=1; A[28,23]=1; A[28,27]=1; A[28,29]=1; A[28,30]=1; A[28,31]=1
A[29,23]=1; A[29,26]=1; A[29,28]=1; A[29,30]=1; A[29,31]=1
A[30,5]=1; A[30,27]=1; A[30,28]=1; A[30,29]=1
A[31,26]=1; A[31,28]=1; A[31,29]=1
printf("对称性检验（应为0）: %f\n", max(abs(A - A')))
for (i=1; i<=n; i++) {
    rs = sum(A[i,.])
    if (rs > 0) A[i,.] = A[i,.] / rs
}
st_matrix("W_adj_raw", A)
end

* ---- 空间自相关诊断 ----
gen W_ln_tfp = .
forvalues y = 2011/2017 {
    preserve
        keep if year==`y' & sample_sdm==1
        sort province_id
        mkmat ln_tfp, matrix(Y_`y')
    restore
    mata: spatial_lag("Y_`y'", "W_econ_raw", "WY_`y'")
    local row = 1
    forvalues p = 1/31 {
        quietly replace W_ln_tfp = el("WY_`y'", `row', 1) ///
            if year==`y' & province_id==`p' & sample_sdm==1
        local ++row
    }
}
pwcorr ln_tfp W_ln_tfp if sample_sdm==1, star(0.05)
di "空间自相关诊断完成"

* ---- 验证SDM样本严格平衡 ----
* xsmle对缺失值零容忍，在此强制检查
quietly count if sample_sdm==1 & (ln_tfp==. | L3_pub_emp==. | broad_tax_burden==. | L3_aging==. | urban_rate==. | sec_ind_share==.)
if r(N) > 0 {
    di as err "SDM样本中存在 " r(N) " 个缺失值，请检查"
}
else {
    di "SDM样本无缺失值，可以运行"
}

* ---- SDM估计：经济距离矩阵 ----
* 注意：xsmle用sample_sdm（不含trade_openness）确保平衡
xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_sdm==1, ///
    wmat(W_econ_raw) model(sdm) fe vce(robust)
est store sdm_econ

* 效应分解
xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_sdm==1, ///
    wmat(W_econ_raw) model(sdm) fe impacts vce(robust)

* SAR基准
xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_sdm==1, ///
    wmat(W_econ_raw) model(sar) fe vce(robust)
est store sar_econ

* SEM基准
xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_sdm==1, ///
    wmat(W_econ_raw) model(sem) fe vce(robust)
est store sem_econ

* 地理邻接矩阵稳健性
xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_sdm==1, ///
    wmat(W_adj_raw) model(sdm) fe vce(robust)
est store sdm_adj

xsmle ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_sdm==1, ///
    wmat(W_adj_raw) model(sdm) fe impacts vce(robust)

* 输出
esttab sdm_econ sar_econ sem_econ ///
    using "spatial_main.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间计量模型（经济距离权重矩阵，2011—2017年）") ///
    mtitles("SDM" "SAR" "SEM") ///
    note("经济距离矩阵：w_ij=1/|人均GDP差距|，行标准化。" ///
         "SDM自动包含W·Y和W·X项。省份固定效应，稳健标准误。" ///
         "控制变量不含trade_openness以确保面板严格平衡（N=217）。")

esttab sdm_econ sdm_adj ///
    using "spatial_robust.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间权重矩阵稳健性：经济距离 vs 地理邻接（SDM）") ///
    mtitles("SDM经济距离" "SDM地理邻接")

/*-----------------------------------------------------------
  完成
-----------------------------------------------------------*/
di as txt _newline "========================================"
di as txt "全部分析完成"
di as txt "输出文件："
di as txt "  desc_stats.rtf       描述性统计"
di as txt "  main_regression.rtf  主回归"
di as txt "  mediation.rtf        中介效应"
di as txt "  heterogeneity.rtf    异质性分析"
di as txt "  robustness.rtf       稳健性（原有四组）"
di as txt "  robustness_dea.rtf   DEA稳健性"
di as txt "  robustness_sfa.rtf   SFA稳健性"
di as txt "  spatial_main.rtf     SDM空间计量主表"
di as txt "  spatial_robust.rtf   SDM权重矩阵稳健性"
di as txt "========================================"

log close
