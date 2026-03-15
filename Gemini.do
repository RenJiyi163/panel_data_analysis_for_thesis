/*===========================================================
   修复版：master_all_v2.do
   修复点：
   1. DEA：改为自动安装
   2. SFA：bc92 只接受 tnormal；半正态改用 bc95 模型
   3. SDM：修复空间权重矩阵导入逻辑，移除冗余命令
===========================================================*/

clear all
set more off
capture log close
log using "master_log.txt", replace text

program define build_vars
    import delimited "panel_with_broad_tax.csv", encoding(UTF-8) clear
    gen pub_emp_share2 = emp_urban_unit_10k / (emp_urban_unit_10k + private_emp_10k)
    gen labor          = emp_urban_unit_10k + private_emp_10k
    scalar delta = 0.096
    scalar alpha = 0.5
    scalar g_inv = 0.161
    xtset province_id year
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
    gen ln_pub_emp = ln(pub_emp_share2)
    gen ln_aging   = ln(aging_rate)
    gen ln_urban   = ln(urban_rate)
    gen ln_trade   = ln(trade_openness) if trade_openness > 0
    gen ln_tfp     = ln_gdp_val - 0.5*ln_capital - 0.5*ln_labor
    gen ln_tfp_alt = ln_gdp_val - 0.4*ln_capital - 0.6*ln_labor
    gen tech_density    = (pub_tech_staff / gdp_bn) * 1000
    gen ln_tech_density = ln(tech_density)
    gen L3_pub_emp = L3.pub_emp_share2
    gen L5_pub_emp = L5.pub_emp_share2
    gen L3_aging   = L3.aging_rate
    gen L5_aging   = L5.aging_rate
    gen sample_main = (year >= 2008 & year <= 2017) & capital != . & ln_tfp != .
    gen sample_lag  = (year >= 2011 & year <= 2017) ///
        & L3_pub_emp != . & L3_aging != . & ln_tfp != .
    gen sample_l5   = (year >= 2013 & year <= 2017) ///
        & L5_pub_emp != . & L5_aging != . & ln_tfp != .
end

/*===========================================================
  模块3：空间杜宾模型（SDM）—— 最终修正版
  修正：spmatrix import 正确语法（去掉 matrix() 包裹）
===========================================================*/
di _newline "======== 模块3：SDM ========"

build_vars  // 确保变量已生成（若前面已调用可省略）

* ---- 经济距离矩阵（基于人均GDP均值） ----
bysort province_id: egen mean_gdppc = mean(gdp_per_cap) ///
    if sample_main == 1 & gdp_per_cap > 0

preserve
    keep if sample_main == 1 & year == 2010  // 任选一年，仅用于构建截面权重
    keep province_id mean_gdppc
    sort province_id
    quietly sum mean_gdppc
    replace mean_gdppc = r(mean) if mean_gdppc == . | mean_gdppc == 0
    mkmat mean_gdppc, matrix(GDPPC)
restore

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
    printf("经济距离矩阵行和均值（应为1）: %f\n", mean(rowsum(W_e)))
    st_matrix("W_econ_raw", W_e)
end

* ---- 地理邻接矩阵 ----
mata:
    n = 31
    A = J(n, n, 0)
    A[1,2]=1; A[1,3]=1
    A[2,1]=1; A[2,3]=1
    A[3,1]=1; A[3,2]=1; A[3,4]=1; A[3,5]=1; A[3,6]=1; A[3,15]=1; A[3,16]=1
    A[4,3]=1; A[4,5]=1; A[4,16]=1; A[4,27]=1
    A[5,3]=1; A[5,4]=1; A[5,6]=1; A[5,7]=1; A[5,8]=1; A[5,27]=1; A[5,28]=1; A[5,30]=1
    A[6,3]=1; A[6,5]=1; A[6,7]=1
    A[7,5]=1; A[7,6]=1; A[7,8]=1
    A[8,5]=1; A[8,7]=1
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

* 声明空间数据（需先设定面板）
xtset province_id year
spset province_id

* 导入权重矩阵（关键修正：直接使用矩阵名，不加 matrix()）
spmatrix import W_econ = W_econ_raw, id(province_id) replace
spmatrix import W_adj  = W_adj_raw,  id(province_id) replace

* 可选：使用 estat moran 进行空间自相关检验（需先运行一个不含空间项的模型）
reg ln_tfp L3_pub_emp broad_tax_burden L3_aging urban_rate sec_ind_share i.year if sample_lag == 1
spmatrix summarize W_econ  // 查看权重矩阵信息
estat moran, errorlag(W_econ)  // 对残差进行莫兰检验

* ---- SDM主模型：经济距离 ----
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe dvarlag(W_econ) ///
    ivarlag(W_econ: L3_pub_emp broad_tax_burden L3_aging) ///
    vce(cluster province_id)
est store sdm_econ

* 直接效应、间接效应、总效应
estat impact L3_pub_emp broad_tax_burden L3_aging

* ---- SAR基准 ----
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe dvarlag(W_econ) vce(cluster province_id)
est store sar_econ

* ---- SEM基准 ----
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe errorlag(W_econ) vce(cluster province_id)
est store sem_econ

* ---- 稳健性：地理邻接矩阵 ----
spxtregress ln_tfp L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share ///
    if sample_lag == 1, ///
    fe dvarlag(W_adj) ///
    ivarlag(W_adj: L3_pub_emp broad_tax_burden L3_aging) ///
    vce(cluster province_id)
est store sdm_adj
estat impact L3_pub_emp broad_tax_burden L3_aging

* 输出结果
esttab sdm_econ sar_econ sem_econ ///
    using "spatial_main.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间计量模型（经济距离权重矩阵，2011—2017年）") ///
    mtitles("SDM" "SAR" "SEM") ///
    note("经济距离权重矩阵：w_ij=1/|人均GDP差距|，行标准化。")

esttab sdm_econ sdm_adj ///
    using "spatial_robust.rtf", replace ///
    star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
    title("空间权重矩阵稳健性：经济距离 vs 地理邻接") ///
    mtitles("SDM经济距离" "SDM地理邻接")

di "SDM模块完成"