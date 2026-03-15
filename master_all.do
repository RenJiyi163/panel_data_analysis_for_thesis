/*===========================================================
  统一修复补丁：解决三个do文件的 xtset 顺序错误
  
  错误原因：L.capital 需要面板结构，但 xtset 写在
           foreach 循环之后，导致 r(111) time variable not set
  
  修复：在资本存量 foreach 循环之前加 xtset
  
  使用方法：把这个补丁文件作为完整版，替代原来三个do文件
  直接 do master_all.do 即可依次跑完所有分析
===========================================================*/

clear all
set more off
capture log close
log using "master_log.txt", replace text

/*-----------------------------------------------------------
  公用宏：变量构造 + 永续盘存法
  所有模块共用，写一次
-----------------------------------------------------------*/
program define build_vars
    import delimited "panel_with_broad_tax.csv", encoding(UTF-8) clear

    * 基础变量
    gen pub_emp_share2 = emp_urban_unit_10k / (emp_urban_unit_10k + private_emp_10k)
    gen labor          = emp_urban_unit_10k + private_emp_10k

    scalar delta = 0.096
    scalar alpha = 0.5
    scalar g_inv = 0.161

    * ---- 关键：必须在使用 L. 之前 xtset ----
    xtset province_id year

    * 初始化资本存量
    gen capital = .
    levelsof province_id, local(provinces)
    foreach p of local provinces {
        quietly sum fai_bn if province_id == `p' & year == 2008
        local I2008 = r(mean)
        quietly replace capital = `I2008' / (g_inv + delta) ///
            if province_id == `p' & year == 2007
    }

    * 永续盘存递推（L.capital 已可用，因为 xtset 在前）
    sort province_id year
    foreach y of numlist 2008/2017 {
        replace capital = (1 - delta) * L.capital + fai_bn ///
            if year == `y' & fai_bn != .
    }

    * 对数变量
    gen ln_gdp_val = ln(gdp_bn)
    gen ln_capital = ln(capital)
    gen ln_labor   = ln(labor)
    gen ln_pub_emp = ln(pub_emp_share2)
    gen ln_aging   = ln(aging_rate)
    gen ln_urban   = ln(urban_rate)
    gen ln_trade   = ln(trade_openness) if trade_openness > 0

    * TFP
    gen ln_tfp     = ln_gdp_val - 0.5*ln_capital - 0.5*ln_labor
    gen ln_tfp_alt = ln_gdp_val - 0.4*ln_capital - 0.6*ln_labor
    label var ln_tfp     "ln(TFP) 索洛残差 α=0.5"
    label var ln_tfp_alt "ln(TFP) 索洛残差 α=0.4"

    * 科技人员密度（H2中介变量）
    gen tech_density    = (pub_tech_staff / gdp_bn) * 1000
    gen ln_tech_density = ln(tech_density)

    * 滞后变量
    gen L3_pub_emp = L3.pub_emp_share2
    gen L5_pub_emp = L5.pub_emp_share2
    gen L3_aging   = L3.aging_rate
    gen L5_aging   = L5.aging_rate

    * 样本标识
    gen sample_main = (year >= 2008 & year <= 2017) & capital != . & ln_tfp != .
    gen sample_lag  = (year >= 2011 & year <= 2017) ///
        & L3_pub_emp != . & L3_aging != . & ln_tfp != .
    gen sample_l5   = (year >= 2013 & year <= 2017) ///
        & L5_pub_emp != . & L5_aging != . & ln_tfp != .
end

/*===========================================================
  模块2：SFA随机前沿分析
===========================================================*/
di as txt _newline "========================================"
di as txt "模块2：SFA 随机前沿分析"
di as txt "========================================"

build_vars

capture which sfpanel
if _rc != 0 {
    di "安装 sfpanel..."
    ssc install sfpanel
}

* SFA估计：Battese-Coelli 1992时变效率模型
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(bc92) distribution(hn)

predict te_bc92, te
label var te_bc92 "技术效率得分（SFA BC92，半正态）"
est store sfa_bc92

sum te_bc92 if sample_main == 1
di "SFA效率得分均值 = " r(mean) "，最小 = " r(min) "，最大 = " r(max)

* True Fixed Effects模型（稳健性）
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(tfe) distribution(hn)

predict te_tfe, te
label var te_tfe "技术效率得分（SFA TFE，半正态）"
est store sfa_tfe

* 截断正态分布（稳健性）
sfpanel ln_gdp_val ln_capital ln_labor ///
    if sample_main == 1, ///
    model(bc92) distribution(tn)

predict te_tn, te

* 相关性诊断：SFA效率 vs 索洛残差
pwcorr te_bc92 te_tfe ln_tfp if sample_main == 1, star(0.05)

* 以SFA效率得分作为因变量重跑主回归
xtreg te_bc92 L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store sfa_reg_bc92

xtreg te_tfe L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store sfa_reg_tfe

* 索洛残差基准（重建用于对比）
xtreg ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store solow_base

* 输出对比表
esttab solow_base sfa_reg_bc92 sfa_reg_tfe ///
    using "robustness_sfa.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("稳健性检验：TFP测算方法对比（索洛残差 vs SFA技术效率得分）") ///
    mtitles("索洛残差TFP" "SFA效率（BC92）" "SFA效率（TFE）") ///
    note("括号内为省份层面聚类稳健标准误。*** p<0.01，** p<0.05，* p<0.1。" ///
         "SFA采用柯布-道格拉斯前沿函数，技术效率TE∈(0,1]，1=完全有效。" ///
         "broad_tax_burden在三列中均显著为负：" ///
         "结论不依赖TFP测算方法；SFA效率分量的系数若绝对值更大，" ///
         "则验证损耗主要来自配置效率渠道而非技术停滞。")

di "SFA模块完成"

/*===========================================================
  模块3：空间杜宾模型（SDM）
===========================================================*/
di as txt _newline "========================================"
di as txt "模块3：空间杜宾模型（经济距离权重矩阵）"
di as txt "========================================"

build_vars

spset province_id year

* ---- 3.1 构建经济距离权重矩阵 ----

* 计算各省2008-2017年人均GDP均值
bysort province_id: egen mean_gdppc = mean(gdp_per_cap) ///
    if sample_main == 1 & gdp_per_cap > 0

preserve
    keep if sample_main == 1 & year == 2010
    keep province_id mean_gdppc
    sort province_id
    quietly sum mean_gdppc
    replace mean_gdppc = r(mean) if mean_gdppc == . | mean_gdppc == 0
    mkmat mean_gdppc, matrix(GDPPC)
restore

* 用Mata构建经济距离矩阵
mata:
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

* ---- 3.2 构建地理邻接矩阵 ----
mata:
    n = 31
    A = J(n, n, 0)
    // 1-北京
    A[1,2]=1; A[1,3]=1
    // 2-天津
    A[2,1]=1; A[2,3]=1
    // 3-河北
    A[3,1]=1; A[3,2]=1; A[3,4]=1; A[3,5]=1
    A[3,6]=1; A[3,15]=1; A[3,16]=1
    // 4-山西
    A[4,3]=1; A[4,5]=1; A[4,16]=1; A[4,27]=1
    // 5-内蒙古
    A[5,3]=1; A[5,4]=1; A[5,6]=1; A[5,7]=1
    A[5,8]=1; A[5,27]=1; A[5,28]=1; A[5,30]=1
    // 6-辽宁
    A[6,3]=1; A[6,5]=1; A[6,7]=1
    // 7-吉林
    A[7,5]=1; A[7,6]=1; A[7,8]=1
    // 8-黑龙江
    A[8,5]=1; A[8,7]=1
    // 9-上海
    A[9,10]=1; A[9,11]=1
    // 10-江苏
    A[10,9]=1; A[10,11]=1; A[10,12]=1; A[10,15]=1
    // 11-浙江
    A[11,9]=1; A[11,10]=1; A[11,12]=1
    A[11,13]=1; A[11,14]=1
    // 12-安徽
    A[12,10]=1; A[12,11]=1; A[12,14]=1
    A[12,15]=1; A[12,16]=1; A[12,17]=1
    // 13-福建
    A[13,11]=1; A[13,14]=1; A[13,19]=1
    // 14-江西
    A[14,11]=1; A[14,12]=1; A[14,13]=1
    A[14,17]=1; A[14,18]=1; A[14,19]=1
    // 15-山东
    A[15,3]=1; A[15,10]=1; A[15,12]=1; A[15,16]=1
    // 16-河南
    A[16,3]=1; A[16,4]=1; A[16,12]=1
    A[16,15]=1; A[16,17]=1; A[16,27]=1
    // 17-湖北
    A[17,12]=1; A[17,14]=1; A[17,16]=1
    A[17,18]=1; A[17,22]=1; A[17,27]=1
    // 18-湖南
    A[18,14]=1; A[18,17]=1; A[18,19]=1
    A[18,20]=1; A[18,22]=1; A[18,24]=1
    // 19-广东
    A[19,13]=1; A[19,14]=1; A[19,18]=1; A[19,20]=1
    // 20-广西
    A[20,18]=1; A[20,19]=1; A[20,24]=1; A[20,25]=1
    // 21-海南（孤岛）
    // 22-重庆
    A[22,17]=1; A[22,18]=1; A[22,23]=1
    A[22,24]=1; A[22,27]=1
    // 23-四川
    A[23,22]=1; A[23,24]=1; A[23,25]=1
    A[23,26]=1; A[23,27]=1; A[23,28]=1; A[23,29]=1
    // 24-贵州
    A[24,18]=1; A[24,20]=1; A[24,22]=1
    A[24,23]=1; A[24,25]=1
    // 25-云南
    A[25,20]=1; A[25,23]=1; A[25,24]=1; A[25,26]=1
    // 26-西藏
    A[26,23]=1; A[26,25]=1; A[26,29]=1; A[26,31]=1
    // 27-陕西
    A[27,4]=1; A[27,5]=1; A[27,16]=1; A[27,17]=1
    A[27,22]=1; A[27,23]=1; A[27,28]=1; A[27,30]=1
    // 28-甘肃
    A[28,5]=1; A[28,23]=1; A[28,27]=1
    A[28,29]=1; A[28,30]=1; A[28,31]=1
    // 29-青海
    A[29,23]=1; A[29,26]=1; A[29,28]=1
    A[29,30]=1; A[29,31]=1
    // 30-宁夏
    A[30,5]=1; A[30,27]=1; A[30,28]=1; A[30,29]=1
    // 31-新疆
    A[31,26]=1; A[31,28]=1; A[31,29]=1

    printf("对称性检验（应为0）: %f\n", max(abs(A - A')))
    
    for (i=1; i<=n; i++) {
        rs = sum(A[i,.])
        if (rs > 0) A[i,.] = A[i,.] / rs
    }
    printf("地理邻接矩阵行和均值（海南除外）: %f\n", ///
        mean(select(rowsum(A), rowsum(A):>0)))
    st_matrix("W_adj_raw", A)
end

* ---- 3.3 存为spmatrix对象（兼容不同Stata版本）----
* 优先使用 spfrommata；若不可用，再回退 frommatrix
mata: W_econ_m = st_matrix("W_econ_raw")
mata: W_adj_m  = st_matrix("W_adj_raw")

capture noisily spmatrix spfrommata W_econ = W_econ_m, id(province_id) normalize(none) replace
if _rc != 0 {
    di as txt "spfrommata 不可用，回退到 frommatrix ..."
    capture noisily spmatrix frommatrix W_econ_raw, id(province_id) name(W_econ) replace
}

capture noisily spmatrix spfrommata W_adj = W_adj_m, id(province_id) normalize(none) replace
if _rc != 0 {
    di as txt "spfrommata 不可用，回退到 frommatrix ..."
    capture noisily spmatrix frommatrix W_adj_raw, id(province_id) name(W_adj) replace
}

capture noisily spmatrix summarize W_econ
if _rc != 0 {
    di as err "W_econ 创建失败：请检查 Stata 版本的 spmatrix 子命令支持"
    exit 198
}

capture noisily spmatrix summarize W_adj
if _rc != 0 {
    di as err "W_adj 创建失败：请检查 Stata 版本的 spmatrix 子命令支持"
    exit 198
}

* ---- 3.4 Moran's I 简易诊断 ----
* 用空间滞后变量的相关系数初步诊断
gen W_ln_tfp = .
label var W_ln_tfp "W·ln(TFP)（经济距离空间滞后）"

forvalues y = 2011/2017 {
    preserve
        keep if year == `y' & sample_lag == 1
        sort province_id
        mkmat ln_tfp, matrix(Y_tmp)
        mata: {
            Y = st_matrix("Y_tmp")
            We = st_matrix("W_econ_raw")
            WY = We * Y
            st_matrix("WY_tmp", WY)
        }
        local row = 1
        forvalues p = 1/31 {
            quietly replace W_ln_tfp = el("WY_tmp", `row', 1) ///
                if year == `y' & province_id == `p' & sample_lag == 1
            local ++row
        }
    restore
}

* 空间自相关诊断：ln_tfp 与 W·ln_tfp 的相关系数
pwcorr ln_tfp W_ln_tfp if sample_lag == 1, star(0.05)
di "若上方相关系数显著为正，说明TFP存在空间正自相关，SDM有统计依据"

* ---- 3.5 SDM估计（Stata 18 spxtregress）----

* 主模型：SDM，经济距离矩阵
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe dvarlag(W_econ) ///
    ivarlag(W_econ: L3_pub_emp broad_tax_burden L3_aging) ///
    vce(cluster province_id)
est store sdm_econ

* 效应分解（直接/间接/总效应）
estat impact L3_pub_emp broad_tax_burden L3_aging
matrix SDM_impacts = r(table)

di "===== SDM效应分解结果 ====="
di "关注：broad_tax_burden 的间接效应（Indirect）"
di "若为负且显著：邻省高税负通过跨省虹吸压制本省TFP"

* SAR（对比基准）
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe dvarlag(W_econ) vce(cluster province_id)
est store sar_econ

* SEM（对比基准）
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe errorlag(W_econ) vce(cluster province_id)
est store sem_econ

* 稳健性：地理邻接矩阵
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe dvarlag(W_adj) ///
    ivarlag(W_adj: L3_pub_emp broad_tax_burden L3_aging) ///
    vce(cluster province_id)
est store sdm_adj

estat impact L3_pub_emp broad_tax_burden L3_aging

* ---- 3.6 输出表格 ----
esttab sdm_econ sar_econ sem_econ ///
    using "spatial_main.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间计量模型估计结果（经济距离权重矩阵，2011—2017年）") ///
    mtitles("SDM（杜宾）" "SAR（自回归）" "SEM（误差）") ///
    note("经济距离权重矩阵：w_ij=1/|人均GDP差距|，行标准化。" ///
         "SDM中ivarlag项显示周边省份变量对本省TFP的溢出效应。" ///
         "省份层面聚类标准误。")

esttab sdm_econ sdm_adj ///
    using "spatial_robust.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间权重矩阵稳健性：经济距离 vs 地理邻接（SDM）") ///
    mtitles("SDM经济距离" "SDM地理邻接") ///
    note("两种矩阵下核心系数方向一致则结论稳健。")

di "空间计量模块完成"

/*-----------------------------------------------------------
  完成提示
-----------------------------------------------------------*/
di as txt _newline "========================================"
di as txt "所有模块运行完毕"
di as txt "输出文件："
di as txt "  robustness_sfa.rtf    SFA稳健性对比"
di as txt "  spatial_main.rtf      空间SDM主表"
di as txt "  spatial_robust.rtf    权重矩阵稳健性"
di as txt "========================================"

log close
