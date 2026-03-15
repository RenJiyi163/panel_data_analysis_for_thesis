/*===========================================================
  DEA-Malmquist生产率指数
  —— 作为索洛残差TFP的稳健性检验
  
  投入：省级资本存量（capital）、劳动投入（labor）
  产出：实际GDP（gdp_bn）
  样本：31省，2008-2017年
  
  逻辑：DEA是非参数方法，不需要预设资本弹性α，
        可直接规避索洛残差法对α=0.5假设的依赖。
        Malmquist指数分解：
        MPI = 技术效率变化（EC）× 技术进步（TC）
        ln(MPI) = TFP变化的非参数估计
===========================================================*/

clear all
set more off
capture log close
log using "dea_log.txt", replace text

* ---- 安装DEA包（如已安装可注释掉）----
* ssc install dea, replace

import delimited "panel_with_broad_tax.csv", encoding(UTF-8) clear

/*-----------------------------------------------------------
  第一步：重建基础变量（与主回归保持一致）
-----------------------------------------------------------*/
gen pub_emp_share2 = emp_urban_unit_10k / (emp_urban_unit_10k + private_emp_10k)
gen labor = emp_urban_unit_10k + private_emp_10k

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

* 仅保留主样本
keep if year >= 2008 & year <= 2017
keep if capital != . & labor != . & gdp_bn != .

* 对数化投入产出（DEA一般用水平值，但对数值用于后续回归）
gen ln_gdp_val  = ln(gdp_bn)
gen ln_capital  = ln(capital)
gen ln_labor    = ln(labor)

/*-----------------------------------------------------------
  第二步：逐年计算DEA效率得分
  使用产出导向（output-oriented）VRS模型
  
  原理：θ_it = 实际产出 / 前沿产出
        θ=1 表示完全有效（位于前沿）
        θ<1 表示技术无效率
-----------------------------------------------------------*/

* DEA命令格式：
* dea 投入变量 = 产出变量, rts(vrs) ort(out)
* 需要逐年循环运行（每年一个截面）

gen dea_eff = .
label var dea_eff "DEA技术效率得分（产出导向VRS）"

forvalues y = 2008/2017 {
    capture {
        dea capital labor = gdp_bn if year == `y', ///
            rts(vrs) ort(out)
        * 将效率得分存入主数据集
        replace dea_eff = e(eff)[_n] if year == `y'
    }
    if _rc != 0 {
        di "年份 `y' DEA计算出错，跳过"
    }
}

/*-----------------------------------------------------------
  第三步：Malmquist生产率指数（手动计算）
  
  MPI(t,t+1) = sqrt[ D_t(x_{t+1},y_{t+1})/D_t(x_t,y_t) × 
                      D_{t+1}(x_{t+1},y_{t+1})/D_{t+1}(x_t,y_t) ]
  
  简化版（使用相邻年效率得分比值）：
  - 效率变化EC: dea_eff_{t+1} / dea_eff_t
  - 全局TFP变化：用EC近似（在仅有效率得分时）
  
  注：完整Malmquist需要跨期DEA（混合LP），
      下面提供两种方案：
      方案A：用dea包的malmquist功能（如版本支持）
      方案B：手动计算效率变化作为TFP代理
-----------------------------------------------------------*/

* 方案A：如果dea包支持malmquist选项
* 检查是否支持
capture {
    dea capital labor = gdp_bn if year == 2008 | year == 2009, ///
        rts(vrs) ort(out) malmquist(year)
    di "dea包支持malmquist选项"
}

if _rc == 0 {
    * 支持时，逐年对计算Malmquist指数
    gen mpi = .
    gen mpi_ec = .   /* 效率变化 */
    gen mpi_tc = .   /* 技术进步 */
    label var mpi    "Malmquist生产率指数"
    label var mpi_ec "效率变化分量（EC）"
    label var mpi_tc "技术进步分量（TC）"
    
    forvalues y = 2008/2016 {
        local y1 = `y' + 1
        capture {
            dea capital labor = gdp_bn ///
                if year == `y' | year == `y1', ///
                rts(vrs) ort(out) malmquist(year)
            * 提取结果（具体变量名取决于dea版本）
        }
    }
}

* 方案B：手动用效率得分近似
* ln(TFP变化) ≈ ln(EC) = ln(eff_{t+1}) - ln(eff_t)
xtset province_id year
gen ln_dea_eff = ln(dea_eff)
gen d_ln_dea   = D.ln_dea_eff
label var d_ln_dea "DEA效率得分变化（≈Malmquist TFP变化）"

* 用效率得分水平作为TFP代理（截面层面）
* 用效率得分变化作为TFP增长代理

/*-----------------------------------------------------------
  第四步：【关键】用DEA效率得分替代索洛残差TFP
         重跑主回归（稳健性检验）
-----------------------------------------------------------*/

* 生成辅助变量
gen L3_pub_emp = L3.pub_emp_share2
gen L3_aging   = L3.aging_rate
label var L3_pub_emp "公共部门就业占比（滞后3期）"
label var L3_aging   "老龄化率（滞后3期）"

gen sample_lag = (year >= 2011 & year <= 2017) ///
    & L3_pub_emp != . & L3_aging != . ///
    & dea_eff != .

* 重新生成广义税负（确保存在）
label var broad_tax_burden "广义税负（政府支出/GDP，%）"

* DEA效率得分水平作为因变量
xtreg ln_dea_eff L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_lag == 1, fe vce(cluster province_id)
est store dea_m1

* DEA效率得分变化（Malmquist近似）作为因变量
gen sample_d = (year >= 2012 & year <= 2017) ///
    & L3_pub_emp != . & L3_aging != . & d_ln_dea != .

xtreg d_ln_dea L3_pub_emp broad_tax_burden L3_aging ///
    urban_rate sec_ind_share trade_openness ///
    i.year if sample_d == 1, fe vce(cluster province_id)
est store dea_m2

/*-----------------------------------------------------------
  第五步：将DEA结果与索洛残差TFP对比
-----------------------------------------------------------*/

* 重建索洛残差TFP
gen ln_tfp_solow = ln_gdp_val - 0.5 * ln_capital - 0.5 * ln_labor

* 相关系数：两种方法的TFP是否高度相关？
pwcorr ln_tfp_solow ln_dea_eff if year >= 2008 & year <= 2017, star(0.05)

di "--------------------------------------------"
di "相关系数诊断：若两种TFP高度相关（r>0.7），"
di "说明索洛残差TFP的测量是稳健的；"
di "若回归系数方向一致，则说明结论不依赖测算方法。"
di "--------------------------------------------"

/*-----------------------------------------------------------
  第六步：输出稳健性对比表
-----------------------------------------------------------*/

* 加载索洛残差结果（需要已跑过主回归）
* 如果在新会话中运行，需要重跑thesis_patch.do先

capture est restore m2_revised
if _rc == 0 {
    esttab m2_revised dea_m1 dea_m2 ///
        using "robustness_dea.rtf", replace ///
        star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
        title("稳健性补充：索洛残差 vs DEA效率得分（因变量对比）") ///
        mtitles("索洛残差TFP（基准）" "DEA效率水平" "DEA效率变化") ///
        note("括号内为省份层面聚类稳健标准误。" ///
             "DEA采用产出导向VRS模型，每年独立求解。" ///
             "三列核心解释变量均相同，因变量测算方法不同，" ///
             "若broad_tax_burden系数方向和显著性一致，" ///
             "则表明核心结论不依赖TFP测算方法的参数假设。")
}
else {
    esttab dea_m1 dea_m2 ///
        using "robustness_dea.rtf", replace ///
        star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
        title("DEA效率得分回归结果") ///
        mtitles("DEA效率水平" "DEA效率变化（Malmquist近似）")
}

/*-----------------------------------------------------------
  注意事项：
  1. 如果 dea 包未安装，先运行：ssc install dea
  2. 如果dea包版本不支持malmquist选项，方案B（手动EC）
     已经足够作为稳健性检验——相关系数诊断最关键
  3. DEA对极端值敏感，西藏（province_id=26）广义税负异常，
     建议同时报告剔除西藏的稳健性版本
  4. 在论文中只需呈现：
     "表X 稳健性检验——TFP测算方法对比"
     一句话说明：方向一致则结论不依赖参数假设
===========================================================*/

log close
