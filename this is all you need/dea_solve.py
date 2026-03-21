"""
dea_solve.py  ——  在Stata外独立运行，或由Stata用shell调用
用法：python dea_solve.py
输入：当前目录下的 dea_input.csv
输出：当前目录下的 dea_efficiency.csv
"""
import pandas as pd
import numpy as np
from scipy.optimize import linprog
import sys, os

input_path  = "dea_input.csv"
output_path = "dea_efficiency.csv"

if not os.path.exists(input_path):
    print(f"错误：找不到 {input_path}，请先在Stata中运行export步骤")
    sys.exit(1)

df = pd.read_csv(input_path)
df = df.dropna(subset=["capital","labor","gdp_bn"])
print(f"读入数据：{len(df)} 行，年份 {sorted(df['year'].unique())}")

results = []
for year in sorted(df["year"].unique()):
    yr = df[df["year"]==year].sort_values("province_id").reset_index(drop=True)
    X  = yr[["capital","labor"]].values   # n×2 投入
    Y  = yr[["gdp_bn"]].values            # n×1 产出
    n  = len(yr)

    for k in range(n):
        xk, yk = X[k], Y[k]
        # 变量：[θ, λ_1,...,λ_n]
        c     = np.zeros(1+n); c[0] = -1.0  # min -θ
        # 投入约束：X·λ ≤ xk
        A1    = np.zeros((2, 1+n))
        for i in range(2): A1[i, 1:] = X[:, i]
        b1    = xk
        # 产出约束：θ·yk - Y·λ ≤ 0
        A2    = np.zeros((1, 1+n))
        A2[0,0] = yk[0]; A2[0,1:] = -Y[:,0]
        b2    = np.array([0.0])
        # VRS：Σλ=1
        A_eq  = np.zeros((1, 1+n)); A_eq[0,1:] = 1.0
        b_eq  = np.array([1.0])
        bounds = [(1.0, None)] + [(0.0, None)] * n
        res   = linprog(c,
                        A_ub=np.vstack([A1,A2]),
                        b_ub=np.concatenate([b1,b2]),
                        A_eq=A_eq, b_eq=b_eq,
                        bounds=bounds, method="highs")
        eff = 1.0/res.x[0] if res.success else float("nan")
        results.append({
            "province_id": int(yr.loc[k, "province_id"]),
            "year":        int(year),
            "dea_eff":     eff
        })

out = pd.DataFrame(results)
out.to_csv(output_path, index=False)
print(f"完成！保存至 {output_path}")
print(out["dea_eff"].describe().round(4).to_string())
