# 細胞計數與種盤稀釋計算器

血球計數盤加 trypan blue 計數，算出**原液濃度、存活率**，以及種盤要**取多少原液、加多少 DMEM**。

**線上版（手機也能用）：https://cell-seeding-calc.vercel.app**

輸入活細胞、死細胞數，選種盤設定（96 孔 / 24 孔，或自填每孔細胞數），即時得到配方、警告與速算提示。

## 三種用法

| 版本 | 檔案 | 啟動 |
|---|---|---|
| 網頁版 | `web/index.html` | 直接用瀏覽器開，或用上面的線上版 |
| 桌面視窗（tkinter） | `cell_seeding_calc.py` | `python cell_seeding_calc.py --gui`，Windows 可雙擊 `Start-CellSeeding.cmd` |
| 命令列 | `cell_seeding_calc.py` | `python cell_seeding_calc.py`（互動，每題有預設值）或一行跑完 |

一行跑完的範例：

```bash
python cell_seeding_calc.py --a-live 200 --a-dead 8 --b-live 208 --b-dead 10 --preset 96d1
```

Python 版只用標準函式庫（3.7+），網頁版沒有任何外部相依。

## 公式

- 總稀釋倍數 `DF` = 預稀釋倍數 × trypan 稀釋倍數，其中倍數 = (細胞液 + 稀釋液) / 細胞液
- 每槽濃度 (cells/mL) = 該槽活細胞總數 / 格數 × `DF` × 10⁴（一個大格 0.1 µL）
- 原液濃度 `C` = 兩槽平均；兩槽差異 % = |上 − 下| / 平均 × 100；存活率 % = 活 / (活 + 死) × 100
- 目標濃度 `Ct` = 每孔細胞數 / 每孔體積(µL) × 1000
- 配製總體積 `Vt` = 孔數 × 每孔體積 × (1 + 餘量%)，進位到 500 µL 的倍數
- 取原液 `V1` = `Ct` × `Vt` / `C`；加 DMEM = `Vt` − `V1`；`V1` < 50 µL 時改用兩段稀釋

判斷規則：每格 > 100 顆太濃、< 20 顆太稀；兩槽差異 > 20% 或存活率 ≤ 90% 視為錯誤；細胞數不夠或 `V1` > `Vt` 視為錯誤；`V1` > 0.6 × `Vt` 提醒原液偏稀。

## 測試

```bash
python -m unittest test_cell_seeding_calc -v   # Python 版（含 tkinter 視窗的煙霧測試）
node --test test_web_core.mjs                  # 網頁版的計算核心
```

兩邊用同一組案例驗證（96d1、24tw、原液 1.2×10⁷ 的兩段稀釋）。

## 紀錄

- Python 版每次執行會把輸入與結果附加到程式旁的 `seeding_log.csv`（已列入 `.gitignore`，不會上傳）。
- 網頁版的紀錄只存在使用者自己的瀏覽器，不會上傳；可用「下載 CSV」匯出，欄位與 Python 版相同。

## 注意

結果僅供實驗參考，請自行核對公式與你的實驗條件。
