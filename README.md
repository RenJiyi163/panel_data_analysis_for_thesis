# 制度约束下的人力资本配置扭曲与全要素生产率损失
本存储库包含用于毕业论文的面板数据和分析脚本。数据文件包含多期面板数据，脚本文件（.do 文件）包括用于数据清洗、估算和计量经济学模型分析的 Stata 代码。分析方法涵盖了数据包络分析（DEA）和随机前沿分析（SFA）以及空间杜宾矩阵（SDM）。本项目旨在探索与制度约束相关的经济学或统计学问题，所有代码和数据已公开以便进一步研究和复现。
# 数据来源
本研究使用的宏观经济数据主要来源于国家统计局官方网站（stats.gov.cn）公开发布的年度统计数据。部分指标的定义和整理参考了学术研究中通用的数据库标准。所有原始数据均可通过公开渠道获取，数据处理代码完全开源，以保证研究的可复现性。
# 具体使用方法
下载"panel_with_broad_tax.csv""full_analysis.do""C1.do",将三个文件置于同一文件夹下。打开STATA18（只在该环境下测试过，其他版本请自行测试是否可行），cd切换到前述文件夹目录，使用do指令先后执行"full_analysis.do""C1.do"即可得出论文实证部分相同的数据结果。
# 注意
full_analysis.do 的SDM模块需要使用 preserve/keep 方式运行，该部分代码已单独放在 C1.do 中。请先执行 full_analysis.do，再执行 C1.do 以完成空间计量部分。

Institutional Constraints, Human Capital Misallocation, and Total Factor Productivity Loss: Evidence from China's Provincial Panel
This repository contains the replication package for a provincial panel study (31 provinces, 2008–2017) examining how institutional constraints — measured by fiscal burden, public sector employment share, and population aging — suppress total factor productivity (TFP) in China.
Methods: Two-way fixed effects, IV (2SLS), mediation analysis, DEA, SFA, and Spatial Durbin Model (SDM)
Key finding: A 1 percentage point increase in the broad fiscal burden ratio (government expenditure/GDP) reduces provincial TFP by approximately 1.0–1.8%, robust across all specifications.
Data: Primary sources are official annual statistics published by China's National Bureau of Statistics (stats.gov.cn). All data values are publicly available.
Replication: Download panel_with_broad_tax.csv, full_analysis.do, and C1.do to the same directory. In Stata 18, cd to that directory, then run full_analysis.do followed by C1.do. (The SDM module requires a separate script due to panel balancing constraints — this is by design, not an error.)
