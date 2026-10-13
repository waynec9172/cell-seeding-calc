# 細胞計數與種盤稀釋計算器

血球計數盤加 trypan blue 計數，支援**種盤稀釋配方**與**切孔／收孔細胞量推算**（支援 PWA 離線使用）。

**線上版（手機也能用）：https://cell-seeding-calc.vercel.app**

## 兩大核心功能

1. **🌱 種盤配方（Seeding）**：輸入計數顆數與種盤設定（96 孔 / 24 孔或自訂），自動計算**原液取量 $V_1$、加多少 DMEM**，並支援 $V_1 < 50\ \mu\text{L}$ 時的自動兩段序列稀釋。
2. **🔬 切孔／收孔細胞量（Harvest）**：從 96 孔盤、24 孔盤、12 孔、6 孔或 Dish 切下 1 孔（或多孔 pooled）下來打散計數，即時換算：
   - **單孔平均活細胞數**（cells/well）與總細胞數
   - **表面生長密度**（cells/cm²）
   - **整盤推估總產量**（全盤孔數同密度推算）
   - **增殖倍數與倍增世代**（相較當初種植細胞數）
   - 可一鍵將收孔算得之濃度帶入種盤配方

## PWA 離線支援（手機可安裝）

網頁版已完整配置 PWA（Progressive Web App）與 Service Worker 快取：
- **離線使用**：在無 Wi-Fi 或無訊號之細胞房（P2 實驗室 / 無菌操作室）中仍可 100% 正常離線計算。
- **安裝至桌面**：
  - Android / Chrome：直接點擊頁面上方「📱 安裝 App」按鈕。
  - iOS / Safari：點擊分享按鈕 ➔ 選擇「加入主畫面」。

## 三種用法

| 版本 | 檔案 | 啟動 |
|---|---|---|
| 網頁版（支援 PWA） | `web/index.html` | 直接用瀏覽器開，或用線上版 |
| 桌面視窗（tkinter） | `cell_seeding_calc.py` | `python cell_seeding_calc.py --gui`，Windows 可雙擊 `Start-CellSeeding.cmd` |
| 命令列 | `cell_seeding_calc.py` | `python cell_seeding_calc.py`（互動模式）或一行跑完 |

一行跑完的範例：

```bash
# 種盤配方
python cell_seeding_calc.py --a-live 200 --a-dead 8 --b-live 208 --b-dead 10 --preset 96d1

# 切孔／收孔產量分析（例如 96 孔切 1 孔，100 µL DMEM 打散）
python cell_seeding_calc.py --mode harvest --harvest-plate 96 --harvest-vol-ul 100 --a-live 80 --a-dead 4 --seeded-cells 10000
```


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
