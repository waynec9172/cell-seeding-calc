#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""細胞計數與種盤稀釋計算器（cell_seeding_calc.py）

血球計數盤 + trypan blue 計數 → 原液濃度、存活率 → 種盤要取多少原液、加多少 DMEM。
只用 Python 標準函式庫（3.7+）。

使用說明
--------
圖形介面（輸入即時重算；Windows 可直接雙擊同資料夾的 Start-CellSeeding.cmd）::

    python cell_seeding_calc.py --gui

互動模式（每題顯示預設值，直接按 Enter 採用）::

    python cell_seeding_calc.py

一行跑完（有給 --a-live 就不會提問，沒給的欄位用預設值）::

    python cell_seeding_calc.py --a-live 200 --a-dead 8 --b-live 208 --b-dead 10 --preset 96d1

    --dmem-ml 2            打散細胞用 DMEM 體積 (mL)
    --pre-cell 10 --pre-dmem 0   計數前預稀釋：細胞液 µL、DMEM µL（預設不稀釋）
    --tb-cell 10 --tb-dye 10     Trypan blue 混合：細胞液 µL、trypan blue µL
    --grids 4              每槽數幾個大格（4 或 9）
    --a-live/--a-dead      上槽活、死細胞總數
    --b-live/--b-dead      下槽活、死細胞總數（省略就只用上槽）
    --preset               96d1 / 96d3 / 96d5 / 24tw / 24d1 / custom
    --cells-per-well --well-ul --wells   自訂設定（也可用來覆寫 preset）
    --margin 25            多抓的餘量 %
    --no-log / --log-file  不寫紀錄 / 指定紀錄檔（預設為本程式旁的 seeding_log.csv）

每次執行結束，會把日期時間、所有輸入與結果附加到 seeding_log.csv（不存在就自動建立）。
測試：python -m unittest test_cell_seeding_calc -v
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import sys
from dataclasses import dataclass
from typing import Callable, Dict, List, NamedTuple, Optional, Tuple

# ---------------------------------------------------------------- 常數 ----

COUNT_FACTOR = 1e4            # 一個大格 0.1 µL = 10^-4 mL → 濃度 = 每格顆數 × 10^4 /mL
VALID_GRIDS = (4, 9)
DENSE_PER_GRID = 100          # 每槽平均每格 > 100：太濃
SPARSE_PER_GRID = 20          # 每槽平均每格 < 20：太稀
MAX_DIFF_PCT = 20.0           # 兩槽差異 > 20%：錯誤
MIN_VIABILITY_PCT = 90.0      # 存活率 ≤ 90%：錯誤
THIN_STOCK_FRACTION = 0.6     # V1 > 0.6 × Vt：原液偏稀
MIN_PIPETTE_UL = 50.0         # V1 < 50 µL：改兩段稀釋
VT_STEP_UL = 500.0            # 配製總體積進位單位
TWO_STEP_FACTORS = (2, 5, 10, 20, 50)
MIN_TRANSFER_UL = 100.0       # 兩段稀釋中，希望每次至少吸取的體積
MIN_INTERMEDIATE_UL = 1000.0  # 中間稀釋液最小總體積
MIN_STOCK_TAKE_UL = 20.0      # 第一段至少取的原液體積
WELLS_PER_MIX = 3             # 每加幾孔回管子混勻一次

OK, WARN, ERR = "OK", "!", "X"


class Preset(NamedTuple):
    cells_per_well: float
    well_ul: float
    wells: int


PRESETS: Dict[str, Preset] = {
    "96d1": Preset(10000, 100, 39),
    "96d3": Preset(4000, 100, 30),
    "96d5": Preset(1000, 100, 30),
    "24tw": Preset(5000, 500, 12),
    "24d1": Preset(20000, 500, 12),
}
CUSTOM = "custom"
DEFAULT_PRESET = "96d1"


class HarvestPreset(NamedTuple):
    name: str
    area: Optional[float]  # cm²
    wells: int
    default_vol_ul: float


HARVEST_PRESETS: Dict[str, HarvestPreset] = {
    "96": HarvestPreset("96 孔盤", 0.32, 96, 100.0),
    "48": HarvestPreset("48 孔盤", 0.95, 48, 250.0),
    "24": HarvestPreset("24 孔盤", 1.9, 24, 500.0),
    "12": HarvestPreset("12 孔盤", 3.8, 12, 1000.0),
    "6": HarvestPreset("6 孔盤", 9.5, 6, 1500.0),
    "10cm": HarvestPreset("10 cm Dish", 55.0, 1, 5000.0),
    "custom": HarvestPreset("自訂規格", None, 1, 500.0),
}


class Finding(NamedTuple):

    level: str      # OK / ! / X
    message: str


# ------------------------------------------------------------ 數字與驗證 ----

def round_half_up(value: float) -> int:
    """四捨五入（0.5 一律進位，不用 Python 的銀行家捨入）。"""
    return int(math.floor(value + 0.5))


def ceil_to_multiple(value: float, step: float) -> float:
    """無條件進位到 step 的倍數（先 round 掉浮點誤差，避免 5000.0000001 變 5500）。"""
    return math.ceil(round(value / step, 9)) * step


def fmt_num(value: float) -> str:
    """一般數字：整數不帶小數，其餘最多 2 位。"""
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.2f}".rstrip("0").rstrip(".")


def sci(value: float, digits: int = 2) -> str:
    """科學記號，例如 1020000 → 1.02×10^6。"""
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    mag = abs(value)
    exp = math.floor(math.log10(mag))
    mant = round(mag / 10 ** exp, digits)
    if mant >= 10:
        mant /= 10
        exp += 1
    text = f"{mant:.{digits}f}".rstrip("0").rstrip(".")
    return f"{sign}{text}×10^{exp}"


def _finite(value: float, name: str) -> float:
    if value is None or not math.isfinite(value):
        raise ValueError(f"{name}必須是有效的數字")
    return value


def check_non_negative(value: float, name: str) -> float:
    _finite(value, name)
    if value < 0:
        raise ValueError(f"{name}不能是負數（收到 {fmt_num(value)}）")
    return value


def check_positive(value: float, name: str) -> float:
    check_non_negative(value, name)
    if value == 0:
        raise ValueError(f"{name}必須大於 0")
    return value


def check_dilution_factor(df: float) -> float:
    _finite(df, "稀釋倍數")
    if df < 1:
        raise ValueError(f"稀釋倍數不能 < 1（收到 {fmt_num(df)}）")
    return df


def check_grids(grids: float) -> int:
    if grids not in VALID_GRIDS:
        raise ValueError(f"格數只能是 4 或 9（收到 {fmt_num(grids)}）")
    return int(grids)


def number_parser(name: str, positive: bool = False) -> Callable[[str], float]:
    """文字 → float 的解析器（同時給互動模式與 argparse 用）。"""
    def parse(text: str) -> float:
        try:
            value = float(str(text).strip())
        except ValueError:
            raise ValueError(f"{name}要輸入數字（收到 {text!r}）") from None
        return check_positive(value, name) if positive else check_non_negative(value, name)
    return parse


def parse_grids(text: str) -> int:
    try:
        value = float(str(text).strip())
    except ValueError:
        raise ValueError(f"格數只能是 4 或 9（收到 {text!r}）") from None
    return check_grids(value)


def parse_wells(text: str) -> int:
    value = number_parser("孔數", positive=True)(text)
    if not value.is_integer():
        raise ValueError(f"孔數必須是整數（收到 {fmt_num(value)}）")
    return int(value)


def parse_preset(text: str) -> str:
    name = str(text).strip().lower()
    if name in PRESETS or name == CUSTOM:
        return name
    choices = " / ".join(list(PRESETS) + [CUSTOM])
    raise ValueError(f"實驗設定只能是 {choices}（收到 {text!r}）")


# ------------------------------------------------------------ 核心計算 ----
# 以下全為純函式：不讀寫檔案、不 print、不 input。

def dilution_factor(sample_ul: float, diluent_ul: float) -> float:
    """(樣品 + 稀釋液) / 樣品。"""
    check_positive(sample_ul, "細胞液體積")
    check_non_negative(diluent_ul, "稀釋液體積")
    return (sample_ul + diluent_ul) / sample_ul


def total_dilution(pre_cell_ul: float, pre_diluent_ul: float,
                   tb_cell_ul: float, tb_dye_ul: float) -> float:
    """總稀釋倍數 DF = 預稀釋倍數 × trypan 稀釋倍數。"""
    return (dilution_factor(pre_cell_ul, pre_diluent_ul)
            * dilution_factor(tb_cell_ul, tb_dye_ul))


def chamber_concentration(live: float, grids: int, df: float) -> float:
    """單槽濃度 (cells/mL) = 活細胞總數 / 格數 × DF × 10^4。"""
    check_non_negative(live, "活細胞總數")
    check_grids(grids)
    check_dilution_factor(df)
    return live / grids * df * COUNT_FACTOR


def difference_pct(conc_a: float, conc_b: float) -> float:
    """兩槽差異 % = |a − b| / 平均 × 100。"""
    mean = (conc_a + conc_b) / 2
    if mean == 0:
        return 0.0
    return abs(conc_a - conc_b) / mean * 100


def viability_pct(live: float, dead: float) -> float:
    total = live + dead
    if total <= 0:
        raise ValueError("活細胞與死細胞都是 0，無法計算存活率")
    return live / total * 100


def target_concentration(cells_per_well: float, well_ul: float) -> float:
    """目標濃度 Ct (cells/mL) = 每孔細胞數 / 每孔體積 (µL) × 1000。"""
    check_positive(well_ul, "每孔體積")
    return cells_per_well / well_ul * 1000


def prepare_volume(wells: int, well_ul: float, margin_pct: float,
                   step_ul: float = VT_STEP_UL) -> float:
    """配製總體積 Vt (µL) = 孔數 × 每孔體積 × (1 + 餘量/100)，進位到 step 的倍數。"""
    check_positive(wells, "孔數")
    check_positive(well_ul, "每孔體積")
    check_non_negative(margin_pct, "餘量 %")
    return ceil_to_multiple(wells * well_ul * (1 + margin_pct / 100), step_ul)


def quick_k(df: float, grids: int, chambers: int = 2) -> float:
    """速算常數 K：濃度 = (各槽活細胞總數相加) × K。"""
    return df * COUNT_FACTOR / (chambers * grids)


@dataclass(frozen=True)
class TwoStepPlan:
    factor: int              # 中間稀釋倍數 f
    intermediate_ul: float   # 中間稀釋液總體積 Vi
    stock_ul: float          # 第一段：原液 Vi/f
    diluent1_ul: float       # 第一段：DMEM Vi − Vi/f
    transfer_ul: float       # 第二段：中間稀釋液 V1×f
    diluent2_ul: float       # 第二段：DMEM Vt − V1×f
    transfer_ok: bool        # V1×f 是否達到 100 µL（f 用到 50 仍不足時為 False）


def two_step_dilution(v1_ul: float, total_ul: float) -> TwoStepPlan:
    """V1 太小時的兩段稀釋。"""
    check_positive(v1_ul, "V1")
    factor = TWO_STEP_FACTORS[-1]
    for candidate in TWO_STEP_FACTORS:
        if v1_ul * candidate >= MIN_TRANSFER_UL:
            factor = candidate
            break
    transfer = v1_ul * factor
    intermediate = max(MIN_INTERMEDIATE_UL,
                       ceil_to_multiple(transfer + MIN_TRANSFER_UL, 100))
    if intermediate / factor < MIN_STOCK_TAKE_UL:
        intermediate = MIN_STOCK_TAKE_UL * factor
    stock = intermediate / factor
    return TwoStepPlan(
        factor=factor,
        intermediate_ul=intermediate,
        stock_ul=stock,
        diluent1_ul=intermediate - stock,
        transfer_ul=transfer,
        diluent2_ul=total_ul - transfer,
        transfer_ok=transfer >= MIN_TRANSFER_UL,
    )


@dataclass(frozen=True)
class SeedingPlan:
    target_conc: float             # Ct, cells/mL
    total_ul: float                # Vt
    cells_needed: float            # Ct × Vt / 1000
    tube_cells: float              # 整管活細胞 = C × 打散體積
    v1_ul: float                   # 要取的原液
    dmem_ul: float                 # Vt − V1
    two_step: Optional[TwoStepPlan]

    @property
    def feasible(self) -> bool:
        return self.v1_ul <= self.total_ul and self.cells_needed <= self.tube_cells


def plan_seeding(stock_conc: float, cells_per_well: float, well_ul: float,
                 wells: int, margin_pct: float = 25.0,
                 tube_ml: float = 2.0) -> SeedingPlan:
    """由原液濃度與實驗設定算出配方。"""
    check_positive(stock_conc, "原液濃度")
    check_positive(cells_per_well, "每孔細胞數")
    check_positive(tube_ml, "打散用 DMEM 體積")
    ct = target_concentration(cells_per_well, well_ul)
    vt = prepare_volume(wells, well_ul, margin_pct)
    needed = ct * vt / 1000
    v1 = ct * vt / stock_conc
    two_step = two_step_dilution(v1, vt) if v1 < MIN_PIPETTE_UL else None
    return SeedingPlan(
        target_conc=ct, total_ul=vt, cells_needed=needed,
        tube_cells=stock_conc * tube_ml, v1_ul=v1, dmem_ul=vt - v1,
        two_step=two_step,
    )


@dataclass(frozen=True)
class HarvestPlan:
    plate_type: str
    plate_name: str
    harvested_wells: int
    harvest_vol_ul: float
    total_harvest_live: float
    cells_per_well: float
    total_per_well: float
    area: Optional[float]
    density: Optional[float]
    wells_total: Optional[int]
    plate_total: Optional[float]
    seeded_per_well: Optional[float]
    fold: Optional[float]
    doublings: Optional[float]


def plan_harvest(stock_conc: float, viability_pct: float,
                 harvest_vol_ul: float, harvested_wells: int,
                 plate_type: str = "96",
                 custom_area: Optional[float] = None,
                 custom_wells: Optional[int] = None,
                 seeded_per_well: Optional[float] = None) -> HarvestPlan:
    """由懸浮液濃度與收孔條件計算單孔細胞數、生長密度與整盤推估。"""
    check_positive(stock_conc, "懸浮液濃度")
    check_positive(harvest_vol_ul, "收孔打散體積")
    check_positive(harvested_wells, "收取的孔數")

    preset = HARVEST_PRESETS.get(plate_type, HARVEST_PRESETS["custom"])
    area = custom_area if plate_type == "custom" else preset.area
    wells_total = custom_wells if plate_type == "custom" else preset.wells

    total_harvest_live = stock_conc * (harvest_vol_ul / 1000.0)
    cells_per_well = total_harvest_live / harvested_wells
    total_per_well = (cells_per_well / (viability_pct / 100.0)) if viability_pct > 0 else cells_per_well
    density = (cells_per_well / area) if area and area > 0 else None
    plate_total = (cells_per_well * wells_total) if wells_total and wells_total > 0 else None

    fold = None
    doublings = None
    if seeded_per_well and seeded_per_well > 0:
        fold = cells_per_well / seeded_per_well
        doublings = math.log2(fold) if fold > 0 else 0.0

    return HarvestPlan(
        plate_type=plate_type,
        plate_name=preset.name,
        harvested_wells=harvested_wells,
        harvest_vol_ul=harvest_vol_ul,
        total_harvest_live=total_harvest_live,
        cells_per_well=cells_per_well,
        total_per_well=total_per_well,
        area=area,
        density=density,
        wells_total=wells_total,
        plate_total=plate_total,
        seeded_per_well=seeded_per_well,
        fold=fold,
        doublings=doublings,
    )


# ------------------------------------------------------------- 整體分析 ----


@dataclass(frozen=True)
class Inputs:
    tube_ml: float
    pre_cell_ul: float
    pre_diluent_ul: float
    tb_cell_ul: float
    tb_dye_ul: float
    grids: int
    a_live: float
    a_dead: float
    b_live: Optional[float]
    b_dead: Optional[float]
    preset: str
    cells_per_well: float
    well_ul: float
    wells: int
    margin_pct: float

    def __post_init__(self) -> None:
        check_positive(self.tube_ml, "打散用 DMEM 體積")
        check_positive(self.pre_cell_ul, "預稀釋的細胞液體積")
        check_non_negative(self.pre_diluent_ul, "預稀釋的 DMEM 體積")
        check_positive(self.tb_cell_ul, "Trypan 混合的細胞液體積")
        check_non_negative(self.tb_dye_ul, "Trypan blue 體積")
        check_grids(self.grids)
        check_non_negative(self.a_live, "上槽活細胞數")
        check_non_negative(self.a_dead, "上槽死細胞數")
        if (self.b_live is None) != (self.b_dead is None):
            raise ValueError("下槽的活細胞數與死細胞數要一起輸入，或兩者都留空")
        if self.b_live is not None:
            check_non_negative(self.b_live, "下槽活細胞數")
            check_non_negative(self.b_dead, "下槽死細胞數")
        check_positive(self.cells_per_well, "每孔細胞數")
        check_positive(self.well_ul, "每孔體積")
        check_positive(self.wells, "孔數")
        check_non_negative(self.margin_pct, "餘量 %")
        if self.a_live + self.a_dead + (self.b_live or 0) + (self.b_dead or 0) <= 0:
            raise ValueError("沒有輸入任何細胞數，無法計算")


@dataclass(frozen=True)
class Analysis:
    inputs: Inputs
    df: float
    conc_a: float
    conc_b: Optional[float]
    per_grid_a: float
    per_grid_b: Optional[float]
    stock_conc: float
    diff_pct: Optional[float]
    viability_pct: float
    tube_cells: float
    k: float
    plan: Optional[SeedingPlan]
    findings: Tuple[Finding, ...]
    harvest_plan: Optional[HarvestPlan] = None

    @property
    def has_error(self) -> bool:
        return any(f.level == ERR for f in self.findings)


def _density_finding(label: str, per_grid: float, grids: int) -> Finding:
    if per_grid > DENSE_PER_GRID:
        return Finding(WARN, f"{label}平均每格 {fmt_num(per_grid)} 顆 > {DENSE_PER_GRID}：太濃，"
                             "建議先稀釋 2 倍再重新計數")
    if per_grid < SPARSE_PER_GRID:
        hint = "打散時用少一點 DMEM" if grids == 9 else "改數 9 格，或打散時用少一點 DMEM"
        return Finding(WARN, f"{label}平均每格 {fmt_num(per_grid)} 顆 < {SPARSE_PER_GRID}：太稀，建議{hint}")
    return Finding(OK, f"{label}平均每格 {fmt_num(per_grid)} 顆"
                       f"（{SPARSE_PER_GRID}–{DENSE_PER_GRID} 之間，適合計數）")


def analyze(inp: Inputs) -> Analysis:
    """輸入 → 全部計算結果與警告。"""
    df = total_dilution(inp.pre_cell_ul, inp.pre_diluent_ul, inp.tb_cell_ul, inp.tb_dye_ul)
    conc_a = chamber_concentration(inp.a_live, inp.grids, df)
    per_grid_a = inp.a_live / inp.grids
    if inp.b_live is not None:
        conc_b: Optional[float] = chamber_concentration(inp.b_live, inp.grids, df)
        per_grid_b: Optional[float] = inp.b_live / inp.grids
        stock = (conc_a + conc_b) / 2
        diff: Optional[float] = difference_pct(conc_a, conc_b)
        live_total = inp.a_live + inp.b_live
        dead_total = inp.a_dead + inp.b_dead
    else:
        conc_b = per_grid_b = diff = None
        stock = conc_a
        live_total, dead_total = inp.a_live, inp.a_dead
    viability = viability_pct(live_total, dead_total)
    tube_cells = stock * inp.tube_ml

    findings: List[Finding] = [_density_finding("上槽", per_grid_a, inp.grids)]
    if per_grid_b is not None:
        findings.append(_density_finding("下槽", per_grid_b, inp.grids))

    if diff is None:
        findings.append(Finding(WARN, "只用上槽計數，無法評估兩槽差異"))
    elif diff > MAX_DIFF_PCT:
        findings.append(Finding(ERR, f"兩槽差異 {diff:.2f}% > {MAX_DIFF_PCT:g}%：細胞可能沒打散，"
                                     "請重新吹打、重新計數"))
    else:
        findings.append(Finding(OK, f"兩槽差異 {diff:.2f}%（≤ {MAX_DIFF_PCT:g}%）"))

    if viability <= MIN_VIABILITY_PCT:
        findings.append(Finding(ERR, f"存活率 {viability:.2f}% ≤ {MIN_VIABILITY_PCT:g}%：這批細胞不建議使用"))
    else:
        findings.append(Finding(OK, f"存活率 {viability:.2f}%（> {MIN_VIABILITY_PCT:g}%）"))
    if dead_total == 0:
        findings.append(Finding(WARN, "死細胞總數為 0，存活率 100% 僅供參考，請確認是否漏填死細胞數"))

    plan: Optional[SeedingPlan] = None
    if stock <= 0:
        findings.append(Finding(ERR, "原液濃度為 0（沒有計到活細胞），無法計算配方"))
    else:
        plan = plan_seeding(stock, inp.cells_per_well, inp.well_ul, inp.wells,
                            inp.margin_pct, inp.tube_ml)
        if plan.cells_needed > plan.tube_cells:
            findings.append(Finding(ERR, f"需要 {sci(plan.cells_needed)} 顆 > 整管活細胞 "
                                         f"{sci(plan.tube_cells)} 顆：細胞不夠"))
        else:
            findings.append(Finding(OK, f"需要 {sci(plan.cells_needed)} 顆 ≤ 整管活細胞 "
                                        f"{sci(plan.tube_cells)} 顆"))
        if plan.v1_ul > plan.total_ul:
            findings.append(Finding(ERR, f"V1 {fmt_num(plan.v1_ul)} µL > 總體積 {fmt_num(plan.total_ul)} µL：原液比目標濃度"
                                         "還稀，請重新離心、用較少 DMEM 打散"))
        elif plan.v1_ul > THIN_STOCK_FRACTION * plan.total_ul:
            findings.append(Finding(WARN, f"V1 {fmt_num(plan.v1_ul)} µL 超過總體積的 "
                                          f"{THIN_STOCK_FRACTION:.0%}：原液偏稀，仍可配製"))
        else:
            findings.append(Finding(OK, f"V1 {fmt_num(plan.v1_ul)} µL，約為總體積的 "
                                        f"{plan.v1_ul / plan.total_ul:.0%}"))
        if plan.two_step is not None:
            findings.append(Finding(WARN, f"V1 {fmt_num(plan.v1_ul)} µL < {MIN_PIPETTE_UL:g} µL：吸量誤差太大，"
                                          f"改用兩段稀釋（先稀釋 {plan.two_step.factor} 倍）"))
            if not plan.two_step.transfer_ok:
                findings.append(Finding(WARN, f"即使稀釋 {plan.two_step.factor} 倍，第二段仍只需吸 "
                                              f"{fmt_num(plan.two_step.transfer_ul)} µL "
                                              f"（< {MIN_TRANSFER_UL:g} µL），建議再多稀釋一次或增加配製量"))

    return Analysis(
        inputs=inp, df=df, conc_a=conc_a, conc_b=conc_b,
        per_grid_a=per_grid_a, per_grid_b=per_grid_b, stock_conc=stock,
        diff_pct=diff, viability_pct=viability, tube_cells=tube_cells,
        k=quick_k(df, inp.grids, 2 if conc_b is not None else 1),
        plan=plan, findings=tuple(findings),
    )


def analyze_harvest(inp: Inputs, plate_type: str = "96",
                    harvest_vol_ul: float = 100.0,
                    harvested_wells: int = 1,
                    custom_area: Optional[float] = None,
                    custom_wells: Optional[int] = None,
                    seeded_per_well: Optional[float] = None) -> Analysis:
    """收孔細胞分析：結合計數結果與收孔條件。"""
    base = analyze(inp)
    if base.stock_conc <= 0:
        return base
    h = plan_harvest(
        base.stock_conc, base.viability_pct, harvest_vol_ul,
        harvested_wells, plate_type, custom_area, custom_wells, seeded_per_well
    )
    extra_findings: List[Finding] = []
    if h.density is not None:
        if h.density > 2.5e5:
            extra_findings.append(Finding(WARN, f"生長密度 {sci(h.density)} cells/cm² 偏高（> 2.5×10⁵）：細胞可能已 100% 滿盤（過密），建議縮短培養時間"))
        elif h.density < 2e4:
            extra_findings.append(Finding(WARN, f"生長密度 {sci(h.density)} cells/cm² 偏低（< 2×10⁴）：約 < 20% confluence，若要做收樣細胞量可能不足"))
        else:
            extra_findings.append(Finding(OK, f"生長密度 {sci(h.density)} cells/cm²：約在對數生長期／適合收樣區間"))
    if h.fold is not None:
        if h.fold < 1:
            extra_findings.append(Finding(WARN, f"收孔細胞數小於當初種植數（增殖倍數 {fmt_num(h.fold)}× < 1）：細胞可能未貼附或生長停滯"))
        else:
            doublings_str = fmt_num(round(h.doublings, 1)) if h.doublings is not None else "0"
            extra_findings.append(Finding(OK, f"較當初種植增殖 {fmt_num(h.fold)} 倍（約經歷 {doublings_str} 次分裂倍增）"))

    return Analysis(
        inputs=base.inputs, df=base.df, conc_a=base.conc_a, conc_b=base.conc_b,
        per_grid_a=base.per_grid_a, per_grid_b=base.per_grid_b, stock_conc=base.stock_conc,
        diff_pct=base.diff_pct, viability_pct=base.viability_pct, tube_cells=base.tube_cells,
        k=base.k, plan=base.plan, findings=tuple(list(base.findings) + extra_findings),
        harvest_plan=h
    )


# ---------------------------------------------------------------- 輸出 ----


def _ul(value: float) -> str:
    return f"{round_half_up(value)} µL"


def recipe_steps(inp: Inputs, plan: SeedingPlan) -> List[str]:
    """配方步驟（µL 四捨五入到整數；DMEM 用總體積相減，確保加總吻合）。"""
    vt = round_half_up(plan.total_ul)
    steps: List[str] = []
    if plan.two_step is None:
        v1 = round_half_up(plan.v1_ul)
        steps.append(f"取原液 {v1} µL 到新管")
        steps.append(f"加 DMEM {vt - v1} µL（總體積 {vt} µL），輕輕吹打混勻")
    else:
        t = plan.two_step
        vi = round_half_up(t.intermediate_ul)
        stock = round_half_up(t.stock_ul)
        transfer = round_half_up(t.transfer_ul)
        steps.append(f"【第一段，稀釋 {t.factor} 倍】取原液 {stock} µL 到新管")
        steps.append(f"加 DMEM {vi - stock} µL（中間稀釋液共 {vi} µL），輕輕吹打混勻")
        steps.append(f"【第二段】取中間稀釋液 {transfer} µL 到另一支新管")
        steps.append(f"加 DMEM {vt - transfer} µL（總體積 {vt} µL），輕輕吹打混勻")
    steps.append(f"每孔加 {fmt_num(inp.well_ul)} µL，每加 {WELLS_PER_MIX} 孔回管子混勻一次"
                 f"（共 {inp.wells} 孔）")
    return steps


def render_report(a: Analysis) -> str:
    inp = a.inputs
    pre = (inp.pre_cell_ul + inp.pre_diluent_ul) / inp.pre_cell_ul
    tb = (inp.tb_cell_ul + inp.tb_dye_ul) / inp.tb_cell_ul
    lines: List[str] = []
    add = lines.append

    add("=== 計數結果 ===")
    add(f"總稀釋倍數 DF：{fmt_num(round(a.df, 4))}"
        f"（預稀釋 {fmt_num(round(pre, 4))} × trypan {fmt_num(round(tb, 4))}）")
    add(f"上槽濃度：{sci(a.conc_a)} cells/mL（平均每格 {fmt_num(a.per_grid_a)} 顆）")
    if a.conc_b is not None:
        add(f"下槽濃度：{sci(a.conc_b)} cells/mL（平均每格 {fmt_num(a.per_grid_b)} 顆）")
        add(f"兩槽差異：{a.diff_pct:.2f} %")
    else:
        add("下槽：未輸入，只用上槽")
    add(f"存活率：{a.viability_pct:.2f} %")
    add(f"原液濃度 C：{sci(a.stock_conc)} cells/mL")
    add(f"整管活細胞（{fmt_num(inp.tube_ml)} mL）：{sci(a.tube_cells)} 顆")

    add("")
    add("=== 判斷 ===")
    for f in a.findings:
        add(f"[{f.level}] {f.message}")

    add("")
    add("=== 配方 ===")
    p = a.plan
    if p is None:
        add("無法配製（原液濃度為 0）。")
    else:
        add(f"目標：每孔 {fmt_num(inp.cells_per_well)} cells / {fmt_num(inp.well_ul)} µL"
            f"（{inp.preset}）→ Ct = {sci(p.target_conc)} cells/mL")
        add(f"配製總體積 Vt = {inp.wells} 孔 × {fmt_num(inp.well_ul)} µL × (1 + {fmt_num(inp.margin_pct)}%)"
            f"，進位到 {VT_STEP_UL:g} 的倍數 = {_ul(p.total_ul)}；需要細胞 {sci(p.cells_needed)} 顆")
        if not p.feasible:
            add("此批細胞無法依目標配製（見上方 [X] 項目）。以下數字僅供參考：")
        for i, step in enumerate(recipe_steps(inp, p), 1):
            add(f"  {i}. {step}")

    add("")
    add("=== 速算提示 ===")
    if a.conc_b is not None:
        add(f"濃度 = (上槽 + 下槽) × K，K = DF × 10^4 / (2 × 格數) = "
            f"{fmt_num(round(a.k, 4))}")
        live_sum = inp.a_live + inp.b_live
        add(f"本次：({fmt_num(inp.a_live)} + {fmt_num(inp.b_live)}) × {fmt_num(round(a.k, 4))}"
            f" = {sci(live_sum * a.k)} cells/mL")
    else:
        add(f"濃度 = 上槽 × K，K = DF × 10^4 / 格數 = {fmt_num(round(a.k, 4))}")
        add(f"本次：{fmt_num(inp.a_live)} × {fmt_num(round(a.k, 4))} = {sci(inp.a_live * a.k)} cells/mL")
    return "\n".join(lines)


def render_harvest_report(a: Analysis) -> str:
    inp = a.inputs
    h = a.harvest_plan
    pre = (inp.pre_cell_ul + inp.pre_diluent_ul) / inp.pre_cell_ul
    tb = (inp.tb_cell_ul + inp.tb_dye_ul) / inp.tb_cell_ul
    lines: List[str] = []
    add = lines.append

    add("=== 計數結果 ===")
    add(f"總稀釋倍數 DF：{fmt_num(round(a.df, 4))}"
        f"（預稀釋 {fmt_num(round(pre, 4))} × trypan {fmt_num(round(tb, 4))}）")
    add(f"上槽濃度：{sci(a.conc_a)} cells/mL（平均每格 {fmt_num(a.per_grid_a)} 顆）")
    if a.conc_b is not None:
        add(f"下槽濃度：{sci(a.conc_b)} cells/mL（平均每格 {fmt_num(a.per_grid_b)} 顆）")
        add(f"兩槽差異：{a.diff_pct:.2f} %")
    else:
        add("下槽：未輸入，只用上槽")
    add(f"存活率：{a.viability_pct:.2f} %")
    add(f"懸浮液濃度 C：{sci(a.stock_conc)} cells/mL")
    if h is not None:
        add(f"收下總活細胞（{fmt_num(h.harvest_vol_ul)} µL）：{sci(h.total_harvest_live)} 顆")

    add("")
    add("=== 判斷 ===")
    for f in a.findings:
        add(f"[{f.level}] {f.message}")

    if h is not None:
        add("")
        add("=== 收孔細胞量分析 ===")
        area_info = f"（單孔面積 {fmt_num(h.area)} cm²，整盤 {h.wells_total} 孔）" if h.area else ""
        add(f"來源孔盤：{h.plate_name}{area_info}")
        add(f"收取孔數：{h.harvested_wells} 孔，打散體積 {fmt_num(h.harvest_vol_ul)} µL")
        add(f"單孔活細胞數：{sci(h.cells_per_well)} cells/well（約 {fmt_num(round_half_up(h.cells_per_well))} 顆）")
        add(f"單孔總細胞數（含死細胞）：{sci(h.total_per_well)} cells/well")
        if h.density is not None:
            add(f"表面生長密度：{sci(h.density)} cells/cm²")
        if h.plate_total is not None:
            add(f"整盤推估總量：{sci(h.plate_total)} 顆（全盤 {h.wells_total} 孔同密度推估）")
        if h.fold is not None:
            doublings_str = fmt_num(round(h.doublings, 1)) if h.doublings is not None else "0"
            add(f"增殖倍數：{fmt_num(h.fold)} 倍（由 {fmt_num(h.seeded_per_well)} 顆種植，約經歷 {doublings_str} 次分裂倍增）")

    add("")
    add("=== 速算提示 ===")
    if a.conc_b is not None:
        add(f"濃度 = (上槽 + 下槽) × K，K = DF × 10^4 / (2 × 格數) = {fmt_num(round(a.k, 4))}")
        live_sum = inp.a_live + inp.b_live
        add(f"本次：({fmt_num(inp.a_live)} + {fmt_num(inp.b_live)}) × {fmt_num(round(a.k, 4))} = {sci(live_sum * a.k)} cells/mL")
    else:
        add(f"濃度 = 上槽 × K，K = DF × 10^4 / 格數 = {fmt_num(round(a.k, 4))}")
        add(f"本次：{fmt_num(inp.a_live)} × {fmt_num(round(a.k, 4))} = {sci(inp.a_live * a.k)} cells/mL")

    return "\n".join(lines)


# ---------------------------------------------------------------- 紀錄 ----

LOG_COLUMNS = [
    "日期時間", "打散DMEM_mL", "預稀釋細胞液_uL", "預稀釋DMEM_uL", "TB細胞液_uL", "TB_uL",
    "總稀釋倍數DF", "格數", "上槽活", "上槽死", "下槽活", "下槽死",
    "實驗設定", "每孔細胞數", "每孔體積_uL", "孔數", "餘量_%",
    "上槽濃度_per_mL", "下槽濃度_per_mL", "原液濃度C_per_mL", "兩槽差異_%", "存活率_%",
    "整管活細胞", "目標濃度Ct_per_mL", "配製總體積Vt_uL", "需要細胞數",
    "原液V1_uL", "DMEM_uL", "稀釋方式", "中間稀釋倍數", "中間液總體積Vi_uL", "第二段中間液_uL",
    "判斷",
]


def build_log_row(a: Analysis, when: dt.datetime) -> Dict[str, object]:
    inp, p = a.inputs, a.plan
    blank = ""
    t = p.two_step if p else None
    if a.harvest_plan is not None:
        h = a.harvest_plan
        return {
            "日期時間": when.strftime("%Y-%m-%d %H:%M:%S"),
            "打散DMEM_mL": round(h.harvest_vol_ul / 1000.0, 4),
            "預稀釋細胞液_uL": inp.pre_cell_ul,
            "預稀釋DMEM_uL": inp.pre_diluent_ul,
            "TB細胞液_uL": inp.tb_cell_ul,
            "TB_uL": inp.tb_dye_ul,
            "總稀釋倍數DF": round(a.df, 6),
            "格數": inp.grids,
            "上槽活": inp.a_live,
            "上槽死": inp.a_dead,
            "下槽活": blank if inp.b_live is None else inp.b_live,
            "下槽死": blank if inp.b_dead is None else inp.b_dead,
            "實驗設定": f"收孔:{h.plate_name}",
            "每孔細胞數": round_half_up(h.cells_per_well),
            "每孔體積_uL": h.harvest_vol_ul,
            "孔數": h.harvested_wells,
            "餘量_%": round(h.fold, 2) if h.fold is not None else blank,
            "上槽濃度_per_mL": f"{a.conc_a:.6g}",
            "下槽濃度_per_mL": blank if a.conc_b is None else f"{a.conc_b:.6g}",
            "原液濃度C_per_mL": f"{a.stock_conc:.6g}",
            "兩槽差異_%": blank if a.diff_pct is None else round(a.diff_pct, 2),
            "存活率_%": round(a.viability_pct, 2),
            "整管活細胞": f"{h.total_harvest_live:.6g}",
            "目標濃度Ct_per_mL": blank,
            "配製總體積Vt_uL": blank,
            "需要細胞數": f"{h.plate_total:.6g}" if h.plate_total else blank,
            "原液V1_uL": blank,
            "DMEM_uL": blank,
            "稀釋方式": "收孔分析",
            "中間稀釋倍數": blank,
            "中間液總體積Vi_uL": blank,
            "第二段中間液_uL": blank,
            "判斷": " | ".join(f"[{f.level}] {f.message}" for f in a.findings if f.level != OK),
        }
    return {
        "日期時間": when.strftime("%Y-%m-%d %H:%M:%S"),
        "打散DMEM_mL": inp.tube_ml,
        "預稀釋細胞液_uL": inp.pre_cell_ul,
        "預稀釋DMEM_uL": inp.pre_diluent_ul,
        "TB細胞液_uL": inp.tb_cell_ul,
        "TB_uL": inp.tb_dye_ul,
        "總稀釋倍數DF": round(a.df, 6),
        "格數": inp.grids,
        "上槽活": inp.a_live,
        "上槽死": inp.a_dead,
        "下槽活": blank if inp.b_live is None else inp.b_live,
        "下槽死": blank if inp.b_dead is None else inp.b_dead,
        "實驗設定": inp.preset,
        "每孔細胞數": inp.cells_per_well,
        "每孔體積_uL": inp.well_ul,
        "孔數": inp.wells,
        "餘量_%": inp.margin_pct,
        "上槽濃度_per_mL": f"{a.conc_a:.6g}",
        "下槽濃度_per_mL": blank if a.conc_b is None else f"{a.conc_b:.6g}",
        "原液濃度C_per_mL": f"{a.stock_conc:.6g}",
        "兩槽差異_%": blank if a.diff_pct is None else round(a.diff_pct, 2),
        "存活率_%": round(a.viability_pct, 2),
        "整管活細胞": f"{a.tube_cells:.6g}",
        "目標濃度Ct_per_mL": blank if p is None else f"{p.target_conc:.6g}",
        "配製總體積Vt_uL": blank if p is None else p.total_ul,
        "需要細胞數": blank if p is None else f"{p.cells_needed:.6g}",
        "原液V1_uL": blank if p is None else round(p.v1_ul, 2),
        "DMEM_uL": blank if p is None else round(p.dmem_ul, 2),
        "稀釋方式": blank if p is None else ("兩段" if t else "一段"),
        "中間稀釋倍數": blank if t is None else t.factor,
        "中間液總體積Vi_uL": blank if t is None else t.intermediate_ul,
        "第二段中間液_uL": blank if t is None else round(t.transfer_ul, 2),
        "判斷": " | ".join(f"[{f.level}] {f.message}" for f in a.findings if f.level != OK),
    }



def append_log(path: str, row: Dict[str, object]) -> None:
    """附加一列；檔案不存在或是空的就先寫欄位名稱（UTF-8 BOM，Excel 直接開得了中文）。"""
    new_file = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", newline="", encoding="utf-8-sig" if new_file else "utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=LOG_COLUMNS)
        if new_file:
            writer.writeheader()
        writer.writerow(row)


# ------------------------------------------------------------ 輸入介面 ----

InputFn = Callable[[str], str]


def ask(label: str, default, parser: Callable[[str], object],
        input_fn: InputFn = input, optional: bool = False, shown: Optional[str] = None):
    """顯示預設值並重複詢問直到輸入合法；Enter = 預設值。"""
    if shown is None:
        if default is not None:
            shown = f"預設 {fmt_num(default) if isinstance(default, (int, float)) else default}"
        else:
            shown = "留空略過" if optional else "必填"
    while True:
        raw = input_fn(f"{label} [{shown}]: ").strip()
        if raw == "":
            if default is not None:
                return default
            if optional:
                return None
            print("  這一題必須輸入。")
            continue
        try:
            return parser(raw)
        except ValueError as exc:
            print(f"  輸入錯誤：{exc}")


def build_parser() -> argparse.ArgumentParser:
    def conv(parser: Callable[[str], object]):
        def wrapped(text: str):
            try:
                return parser(text)
            except ValueError as exc:
                raise argparse.ArgumentTypeError(str(exc)) from None
        return wrapped

    p = argparse.ArgumentParser(
        description="細胞計數與種盤稀釋計算器。不給 --a-live 就進入互動模式。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--dmem-ml", type=conv(number_parser("打散用 DMEM 體積", True)), default=None,
                   help="打散細胞用的 DMEM 體積 (mL)（預設 2）")
    p.add_argument("--pre-cell", type=conv(number_parser("預稀釋的細胞液體積", True)), default=None,
                   help="計數前預稀釋：細胞液 µL（預設 10）")
    p.add_argument("--pre-dmem", type=conv(number_parser("預稀釋的 DMEM 體積")), default=None,
                   help="計數前預稀釋：DMEM µL（預設 0 = 不稀釋）")
    p.add_argument("--tb-cell", type=conv(number_parser("Trypan 混合的細胞液體積", True)), default=None,
                   help="Trypan 混合：細胞液 µL（預設 10）")
    p.add_argument("--tb-dye", type=conv(number_parser("Trypan blue 體積")), default=None,
                   help="Trypan 混合：trypan blue µL（預設 10）")
    p.add_argument("--grids", type=conv(parse_grids), default=None, help="每槽數幾個大格：4 或 9（預設 4）")
    p.add_argument("--a-live", type=conv(number_parser("上槽活細胞數")), default=None, help="上槽活細胞總數")
    p.add_argument("--a-dead", type=conv(number_parser("上槽死細胞數")), default=None, help="上槽死細胞總數（預設 0）")
    p.add_argument("--b-live", type=conv(number_parser("下槽活細胞數")), default=None, help="下槽活細胞總數（可省略）")
    p.add_argument("--b-dead", type=conv(number_parser("下槽死細胞數")), default=None, help="下槽死細胞總數（預設 0）")
    p.add_argument("--preset", type=conv(parse_preset), default=None,
                   help="實驗設定：" + " / ".join(list(PRESETS) + [CUSTOM]) + f"（預設 {DEFAULT_PRESET}）")
    p.add_argument("--cells-per-well", type=conv(number_parser("每孔細胞數", True)), default=None,
                   help="每孔細胞數（custom 必填；也可覆寫 preset）")
    p.add_argument("--well-ul", type=conv(number_parser("每孔體積", True)), default=None,
                   help="每孔體積 µL（custom 必填；也可覆寫 preset）")
    p.add_argument("--wells", type=conv(parse_wells), default=None,
                   help="孔數（custom 必填；也可覆寫 preset）")
    p.add_argument("--margin", type=conv(number_parser("餘量 %")), default=None, help="多抓的餘量 %%（預設 25）")
    p.add_argument("--mode", choices=["seeding", "harvest"], default="seeding",
                   help="模式：seeding (種盤配方) 或 harvest (切孔／收孔細胞量分析)")
    p.add_argument("--harvest-plate", choices=list(HARVEST_PRESETS), default="96",
                   help="收孔來源盤別（預設 96）")
    p.add_argument("--harvest-wells", type=conv(parse_wells), default=1,
                   help="收取孔數（預設 1）")
    p.add_argument("--harvest-vol-ul", type=conv(number_parser("收孔打散體積", True)), default=None,
                   help="收孔打散用 DMEM 體積 µL（96 孔預設 100，24 孔預設 500）")
    p.add_argument("--seeded-cells", type=conv(number_parser("當初每孔種植數", True)), default=None,
                   help="當初每孔種植細胞數（選填，用來算增殖倍數）")
    p.add_argument("--custom-area", type=conv(number_parser("自訂單孔面積", True)), default=None,
                   help="自訂單孔面積 cm²（harvest-plate 為 custom 時使用）")
    p.add_argument("--custom-wells", type=conv(parse_wells), default=None,
                   help="自訂整盤孔數（harvest-plate 為 custom 時使用）")
    p.add_argument("--gui", action="store_true", help="開啟圖形介面（命令列參數的值會當作初始值）")
    p.add_argument("-i", "--interactive", action="store_true", help="即使有給 --a-live 也逐題詢問（參數值當預設）")
    p.add_argument("--log-file", default=None, help="紀錄檔路徑（預設為本程式旁的 seeding_log.csv）")
    p.add_argument("--no-log", action="store_true", help="這次不寫紀錄")
    return p



def default_log_path() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "seeding_log.csv")


def _setup_from_args(args: argparse.Namespace) -> Tuple[str, Optional[float], Optional[float], Optional[int]]:
    """preset 名稱 + 可能由參數覆寫的三個值（未給為 None）。"""
    name = args.preset
    if name is None:
        given = any(v is not None for v in (args.cells_per_well, args.well_ul, args.wells))
        name = CUSTOM if given else DEFAULT_PRESET
    return name, args.cells_per_well, args.well_ul, args.wells


def collect_batch(args: argparse.Namespace) -> Inputs:
    """非互動：沒給的欄位用預設值；custom 必須給齊三個值。"""
    if args.a_live is None:
        raise ValueError("請用 --a-live 提供上槽活細胞數（或不帶參數進入互動模式）")
    if args.b_live is None and args.b_dead is not None:
        raise ValueError("有 --b-dead 就必須有 --b-live")
    name, cells, well_ul, wells = _setup_from_args(args)
    if name != CUSTOM:
        base = PRESETS[name]
        cells = base.cells_per_well if cells is None else cells
        well_ul = base.well_ul if well_ul is None else well_ul
        wells = base.wells if wells is None else wells
    elif None in (cells, well_ul, wells):
        raise ValueError("custom 需要同時給 --cells-per-well、--well-ul、--wells")
    b_live = args.b_live
    b_dead = None if b_live is None else (args.b_dead or 0.0)
    return Inputs(
        tube_ml=2.0 if args.dmem_ml is None else args.dmem_ml,
        pre_cell_ul=10.0 if args.pre_cell is None else args.pre_cell,
        pre_diluent_ul=0.0 if args.pre_dmem is None else args.pre_dmem,
        tb_cell_ul=10.0 if args.tb_cell is None else args.tb_cell,
        tb_dye_ul=10.0 if args.tb_dye is None else args.tb_dye,
        grids=4 if args.grids is None else args.grids,
        a_live=args.a_live, a_dead=args.a_dead or 0.0,
        b_live=b_live, b_dead=b_dead,
        preset=name, cells_per_well=cells, well_ul=well_ul, wells=wells,
        margin_pct=25.0 if args.margin is None else args.margin,
    )


def collect_interactive(args: argparse.Namespace, input_fn: InputFn = input) -> Inputs:
    """互動：逐題詢問；命令列參數的值會變成該題的預設值。"""
    def d(value, fallback):
        return fallback if value is None else value

    print("=== 細胞計數與種盤稀釋計算器 ===（每題直接按 Enter 採用 [ ] 內的預設值）")
    tube_ml = ask("打散細胞用的 DMEM 體積 (mL)", d(args.dmem_ml, 2.0),
                  number_parser("打散用 DMEM 體積", True), input_fn)
    print("計數前預先稀釋（預設不稀釋：10 + 0 → 稀釋倍數 1）")
    pre_cell = ask("  預稀釋：細胞液 (µL)", d(args.pre_cell, 10.0),
                   number_parser("預稀釋的細胞液體積", True), input_fn)
    pre_dil = ask("  預稀釋：DMEM (µL)", d(args.pre_dmem, 0.0),
                  number_parser("預稀釋的 DMEM 體積"), input_fn)
    print("Trypan blue 混合（預設 10 + 10 → 稀釋倍數 2）")
    tb_cell = ask("  細胞液 (µL)", d(args.tb_cell, 10.0),
                  number_parser("Trypan 混合的細胞液體積", True), input_fn)
    tb_dye = ask("  Trypan blue (µL)", d(args.tb_dye, 10.0),
                 number_parser("Trypan blue 體積"), input_fn)
    grids = ask("每槽數幾個大格 (4 或 9)", d(args.grids, 4), parse_grids, input_fn)

    a_live = ask("上槽：活細胞總數", args.a_live, number_parser("上槽活細胞數"), input_fn)
    a_dead = ask("上槽：死細胞總數", d(args.a_dead, 0.0), number_parser("上槽死細胞數"), input_fn)
    b_live = ask("下槽：活細胞總數（留空 = 只用上槽）", args.b_live,
                 number_parser("下槽活細胞數"), input_fn, optional=True)
    b_dead = None
    if b_live is not None:
        b_dead = ask("下槽：死細胞總數", d(args.b_dead, 0.0), number_parser("下槽死細胞數"), input_fn)

    name0, cells, well_ul, wells = _setup_from_args(args)
    print("實驗設定：" + "；".join(
        f"{k} = {fmt_num(v.cells_per_well)} cells/孔, {fmt_num(v.well_ul)} µL, {v.wells} 孔"
        for k, v in PRESETS.items()) + "；custom = 自訂")
    name = ask("實驗設定", name0, parse_preset, input_fn)
    if name == CUSTOM:
        base = PRESETS[DEFAULT_PRESET]
        cells = ask("  每孔細胞數", d(cells, base.cells_per_well), number_parser("每孔細胞數", True), input_fn)
        well_ul = ask("  每孔體積 (µL)", d(well_ul, base.well_ul), number_parser("每孔體積", True), input_fn)
        wells = ask("  孔數", d(wells, base.wells), parse_wells, input_fn)
    else:
        base = PRESETS[name]
        cells, well_ul, wells = d(cells, base.cells_per_well), d(well_ul, base.well_ul), d(wells, base.wells)
    margin = ask("多抓的餘量 (%)", d(args.margin, 25.0), number_parser("餘量 %"), input_fn)

    return Inputs(tube_ml=tube_ml, pre_cell_ul=pre_cell, pre_diluent_ul=pre_dil,
                  tb_cell_ul=tb_cell, tb_dye_ul=tb_dye, grids=grids,
                  a_live=a_live, a_dead=a_dead, b_live=b_live, b_dead=b_dead,
                  preset=name, cells_per_well=cells, well_ul=well_ul, wells=wells,
                  margin_pct=margin)


# ------------------------------------------------------------ 圖形介面 ----
# 表單解析與摘要是純函式（可測試）；視窗本身只負責呈現。

FORM_FIELDS = ("dmem_ml", "pre_cell", "pre_dmem", "tb_cell", "tb_dye", "grids",
               "a_live", "a_dead", "b_live", "b_dead",
               "preset", "cells_per_well", "well_ul", "wells", "margin")


def form_defaults(args: Optional[argparse.Namespace] = None) -> Dict[str, str]:
    """表單的初始字串；命令列參數有給的值會當作初始值。"""
    def pick(name: str, fallback) -> str:
        value = getattr(args, name, None) if args is not None else None
        return fmt_num(fallback if value is None else value)

    values = {
        "dmem_ml": pick("dmem_ml", 2), "pre_cell": pick("pre_cell", 10), "pre_dmem": pick("pre_dmem", 0),
        "tb_cell": pick("tb_cell", 10), "tb_dye": pick("tb_dye", 10), "grids": pick("grids", 4),
        "a_live": "" if getattr(args, "a_live", None) is None else pick("a_live", 0),
        "a_dead": pick("a_dead", 0),
        "b_live": "" if getattr(args, "b_live", None) is None else pick("b_live", 0),
        "b_dead": "" if getattr(args, "b_dead", None) is None else pick("b_dead", 0),
        "margin": pick("margin", 25),
    }
    name = getattr(args, "preset", None) or DEFAULT_PRESET
    base = PRESETS.get(name, PRESETS[DEFAULT_PRESET])
    cells = pick("cells_per_well", base.cells_per_well)
    well_ul = pick("well_ul", base.well_ul)
    wells = pick("wells", base.wells)
    if name != CUSTOM and (cells, well_ul, wells) != tuple(
            fmt_num(v) for v in (base.cells_per_well, base.well_ul, base.wells)):
        name = CUSTOM      # 參數覆寫了 preset 的數值 → 視為自訂
    values.update(preset=name, cells_per_well=cells, well_ul=well_ul, wells=wells)
    return values


def inputs_from_form(form: Dict[str, str]) -> Inputs:
    """把表單的字串轉成 Inputs；有問題就丟 ValueError（訊息直接顯示給使用者）。"""
    def text(key: str) -> str:
        return form.get(key, "").strip()

    def num(key: str, name: str, positive: bool = False) -> float:
        if text(key) == "":
            raise ValueError(f"請輸入{name}")
        return number_parser(name, positive)(text(key))

    if text("grids") == "":
        raise ValueError("請選擇格數（4 或 9）")
    grids = parse_grids(text("grids"))
    cells = num("cells_per_well", "每孔細胞數", True)
    well_ul = num("well_ul", "每孔體積", True)
    if text("wells") == "":
        raise ValueError("請輸入孔數")
    wells = parse_wells(text("wells"))
    # 種盤設定只是「快速帶入」；實際用表單上的數字，與 preset 不同就記為 custom
    preset = parse_preset(text("preset") or CUSTOM)
    if preset != CUSTOM and (cells, well_ul, wells) != tuple(PRESETS[preset]):
        preset = CUSTOM

    b_live: Optional[float] = None
    b_dead: Optional[float] = None
    if text("b_live") != "":
        b_live = num("b_live", "下槽活細胞數")
        b_dead = num("b_dead", "下槽死細胞數") if text("b_dead") != "" else 0.0
    elif text("b_dead") != "":
        raise ValueError("下槽有死細胞數時，也要輸入下槽活細胞數（或清空下槽）")

    return Inputs(
        tube_ml=num("dmem_ml", "打散用 DMEM 體積", True),
        pre_cell_ul=num("pre_cell", "預稀釋的細胞液體積", True),
        pre_diluent_ul=num("pre_dmem", "預稀釋的 DMEM 體積"),
        tb_cell_ul=num("tb_cell", "Trypan 混合的細胞液體積", True),
        tb_dye_ul=num("tb_dye", "Trypan blue 體積"),
        grids=grids,
        a_live=num("a_live", "上槽活細胞數"),
        a_dead=num("a_dead", "上槽死細胞數") if text("a_dead") != "" else 0.0,
        b_live=b_live, b_dead=b_dead,
        preset=preset, cells_per_well=cells, well_ul=well_ul, wells=wells,
        margin_pct=num("margin", "餘量 %"),
    )


def summary_line(a: Analysis) -> Tuple[str, str]:
    """GUI 最上方一行重點：(嚴重程度, 文字)。"""
    if a.has_error:
        level = ERR
    elif any(f.level == WARN for f in a.findings):
        level = WARN
    else:
        level = OK
    text = f"原液 {sci(a.stock_conc)} cells/mL　存活率 {a.viability_pct:.1f}%"
    p = a.plan
    if p is None:
        text += "　無法配製"
    elif not p.feasible:
        text += "　無法依目標配製"
    elif p.two_step is not None:
        text += f"　V1 僅 {fmt_num(round(p.v1_ul, 1))} µL → 兩段稀釋"
    else:
        v1 = round_half_up(p.v1_ul)
        text += f"　取原液 {v1} µL ＋ DMEM {round_half_up(p.total_ul) - v1} µL"
    return level, text


try:  # 沒有 Tcl/Tk 的 Python 仍可使用命令列模式
    import tkinter as tk
    from tkinter import ttk
    from tkinter import font as tkfont
except ImportError:  # pragma: no cover
    tk = ttk = tkfont = None  # type: ignore[assignment]


class SeedingApp:
    """tkinter 視窗：左邊輸入、右邊即時結果；輸入每改一次就重算。"""

    COLORS = {OK: "#1a7f37", WARN: "#b45309", ERR: "#c62828", "muted": "#6b7280"}
    PAD = 6

    def __init__(self, root, initial: Dict[str, str], log_path: str) -> None:
        self.root = root
        self.log_path = log_path
        self.analysis: Optional[Analysis] = None
        self.vars = {k: tk.StringVar(master=root, value=initial.get(k, "")) for k in FORM_FIELDS}
        self.setup_entries: List = []
        self._filling = False   # 正在由 preset 帶入數值（此時不要把 preset 切成 custom）

        root.title("細胞計數與種盤稀釋計算器")
        self._setup_fonts()
        self._build_inputs()
        self._build_output()
        self._on_preset()
        self.vars["preset"].trace_add("write", self._on_preset)   # 先換 preset 數值，再重算
        for key in ("cells_per_well", "well_ul", "wells"):
            self.vars[key].trace_add("write", self._on_setup_edit)
        for var in self.vars.values():
            var.trace_add("write", self.refresh)
        self.refresh()
        self.a_live_entry.focus_set()

    # --- 建構視窗 ---
    def _setup_fonts(self) -> None:
        families = set(tkfont.families(self.root))
        for name in ("TkDefaultFont", "TkTextFont", "TkFixedFont"):
            font = tkfont.nametofont(name)
            if "Microsoft JhengHei UI" in families:
                font.configure(family="Microsoft JhengHei UI")
            font.configure(size=10)
        base = tkfont.nametofont("TkDefaultFont")
        family = base.actual("family")
        self.headline_font = tkfont.Font(self.root, family=family, size=13, weight="bold")
        self.section_font = tkfont.Font(self.root, family=family, size=11, weight="bold")

    def _build_inputs(self) -> None:
        left = ttk.Frame(self.root, padding=self.PAD)
        left.grid(row=0, column=0, sticky="ns")
        pad = self.PAD

        def box(title: str):
            frame = ttk.LabelFrame(left, text=title, padding=pad)
            frame.pack(fill="x", pady=(0, pad))
            return frame

        def entry(parent, row: int, label: str, key: str, unit: str = "", col: int = 0):
            ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", pady=2)
            widget = ttk.Entry(parent, textvariable=self.vars[key], width=9, justify="right")
            widget.grid(row=row, column=col + 1, padx=pad, pady=2)
            if unit:
                ttk.Label(parent, text=unit).grid(row=row, column=col + 2, sticky="w")
            return widget

        dil = box("打散與稀釋")
        muted = self.COLORS["muted"]

        def caption(row: int, text: str):
            ttk.Label(dil, text=text, foreground=muted, wraplength=420, justify="left"
                      ).grid(row=row, column=0, columnspan=4, sticky="w", pady=(6, 0))

        entry(dil, 0, "打散用 DMEM", "dmem_ml", "mL（離心後打散整管用）")

        caption(2, "① 計數前預稀釋（細胞太濃才用；不稀釋就維持 10＋0）")
        entry(dil, 3, "取細胞液", "pre_cell", "µL")
        entry(dil, 4, "加 DMEM", "pre_dmem", "µL")
        self.pre_factor = ttk.Label(dil, text="", font=self.headline_font)
        self.pre_factor.grid(row=3, column=3, rowspan=2, padx=(pad, 0))

        caption(5, "② 與 trypan blue 混合（取預稀釋後的液體；一般 10＋10）")
        entry(dil, 6, "取細胞液", "tb_cell", "µL")
        entry(dil, 7, "加 trypan blue", "tb_dye", "µL")
        self.tb_factor = ttk.Label(dil, text="", font=self.headline_font)
        self.tb_factor.grid(row=6, column=3, rowspan=2, padx=(pad, 0))

        self.df_label = ttk.Label(dil, text="", font=self.section_font)
        self.df_label.grid(row=8, column=0, columnspan=4, sticky="w", pady=(8, 0))

        cnt = box("計數（血球計數盤）")
        ttk.Label(cnt, text="每槽大格數").grid(row=0, column=0, sticky="w")
        grids = ttk.Frame(cnt)
        grids.grid(row=0, column=1, columnspan=2, sticky="w", padx=pad)
        for value in ("4", "9"):
            ttk.Radiobutton(grids, text=f"{value} 格", value=value,
                            variable=self.vars["grids"]).pack(side="left", padx=(0, pad))
        ttk.Label(cnt, text="活").grid(row=1, column=1)
        ttk.Label(cnt, text="死").grid(row=1, column=2)
        ttk.Label(cnt, text="上槽").grid(row=2, column=0, sticky="w")
        self.a_live_entry = ttk.Entry(cnt, textvariable=self.vars["a_live"], width=9, justify="right")
        self.a_live_entry.grid(row=2, column=1, padx=2, pady=2)
        ttk.Entry(cnt, textvariable=self.vars["a_dead"], width=9, justify="right").grid(row=2, column=2, padx=2)
        ttk.Label(cnt, text="下槽").grid(row=3, column=0, sticky="w")
        ttk.Entry(cnt, textvariable=self.vars["b_live"], width=9, justify="right").grid(row=3, column=1, padx=2, pady=2)
        ttk.Entry(cnt, textvariable=self.vars["b_dead"], width=9, justify="right").grid(row=3, column=2, padx=2)
        ttk.Label(cnt, text="填各槽總顆數；下槽可留空（只用上槽）",
                  foreground=self.COLORS["muted"]).grid(row=4, column=0, columnspan=3, sticky="w", pady=(4, 0))

        exp = box("實驗設定")
        ttk.Label(exp, text="快速帶入").grid(row=0, column=0, sticky="w")
        ttk.Combobox(exp, textvariable=self.vars["preset"], values=list(PRESETS) + [CUSTOM],
                     state="readonly", width=8).grid(row=0, column=1, padx=pad, pady=2)
        self.setup_entries = [
            entry(exp, 1, "每孔細胞數", "cells_per_well", "cells"),
            entry(exp, 2, "每孔體積", "well_ul", "µL"),
            entry(exp, 3, "孔數", "wells", "孔"),
        ]
        entry(exp, 4, "多抓的餘量", "margin", "%")
        ttk.Label(exp, text="三格都可直接改，改了會自動變成 custom",
                  foreground=self.COLORS["muted"], wraplength=380, justify="left"
                  ).grid(row=5, column=0, columnspan=3, sticky="w", pady=(4, 0))

        buttons = ttk.Frame(left)
        buttons.pack(fill="x")
        ttk.Button(buttons, text="存入紀錄", command=self.save_log).pack(fill="x", pady=2)
        ttk.Button(buttons, text="複製結果", command=self.copy_report).pack(fill="x", pady=2)
        ttk.Button(buttons, text="清除計數", command=self.clear_counts).pack(fill="x", pady=2)

    def _build_output(self) -> None:
        right = ttk.Frame(self.root, padding=self.PAD)
        right.grid(row=0, column=1, sticky="nsew")
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)

        self.headline = ttk.Label(right, font=self.headline_font, wraplength=520, justify="left")
        self.headline.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, self.PAD))
        self.report = tk.Text(right, width=64, height=26, wrap="word", relief="flat",
                              borderwidth=0, padx=6, pady=4, state="disabled")
        scroll = ttk.Scrollbar(right, command=self.report.yview)
        self.report.configure(yscrollcommand=scroll.set)
        self.report.grid(row=1, column=0, sticky="nsew")
        scroll.grid(row=1, column=1, sticky="ns")
        for level, color in self.COLORS.items():
            self.report.tag_configure(level, foreground=color)
        self.report.tag_configure("head", font=self.section_font)
        self.report.tag_configure("step", lmargin1=10, lmargin2=28)
        self.status = ttk.Label(right, text="", foreground=self.COLORS["muted"])
        self.status.grid(row=2, column=0, columnspan=2, sticky="w", pady=(self.PAD, 0))

    # --- 事件 ---
    def _on_preset(self, *_) -> None:
        name = self.vars["preset"].get()
        if name in PRESETS:
            base = PRESETS[name]
            self._filling = True
            try:
                self.vars["cells_per_well"].set(fmt_num(base.cells_per_well))
                self.vars["well_ul"].set(fmt_num(base.well_ul))
                self.vars["wells"].set(fmt_num(base.wells))
            finally:
                self._filling = False

    def _on_setup_edit(self, *_) -> None:
        name = self.vars["preset"].get()
        if self._filling or name not in PRESETS:
            return
        base = PRESETS[name]
        try:
            same = (float(self.vars["cells_per_well"].get()) == base.cells_per_well
                    and float(self.vars["well_ul"].get()) == base.well_ul
                    and float(self.vars["wells"].get()) == base.wells)
        except ValueError:
            same = False
        if not same:
            self.vars["preset"].set(CUSTOM)

    def _update_factors(self) -> None:
        """即時顯示兩個稀釋倍數與總倍數（輸入未完成時顯示 ?）。"""
        def factor(cell_key: str, diluent_key: str) -> Optional[float]:
            try:
                return dilution_factor(float(self.vars[cell_key].get()), float(self.vars[diluent_key].get()))
            except ValueError:
                return None

        pre, tb = factor("pre_cell", "pre_dmem"), factor("tb_cell", "tb_dye")
        show = lambda f: "?" if f is None else fmt_num(round(f, 3))
        self.pre_factor.configure(text=f"＝ ×{show(pre)}")
        self.tb_factor.configure(text=f"＝ ×{show(tb)}")
        total = "?" if pre is None or tb is None else fmt_num(round(pre * tb, 3))
        self.df_label.configure(text=f"總稀釋倍數 DF ＝ {show(pre)} × {show(tb)} ＝ {total}")

    def refresh(self, *_) -> None:
        self._update_factors()
        try:
            self.analysis = analyze(inputs_from_form({k: v.get() for k, v in self.vars.items()}))
        except ValueError as exc:
            self.analysis = None
            self.headline.configure(text=str(exc), foreground=self.COLORS["muted"])
            self._set_report("")
            return
        level, text = summary_line(self.analysis)
        self.headline.configure(text=text, foreground=self.COLORS[level])
        self._set_report(render_report(self.analysis))

    def _set_report(self, text: str) -> None:
        self.report.configure(state="normal")
        self.report.delete("1.0", "end")
        for line in text.splitlines():
            tag = None
            if line.startswith("=== "):
                tag = "head"
            elif line.startswith("[OK]"):
                tag = OK
            elif line.startswith("[!]"):
                tag = WARN
            elif line.startswith("[X]"):
                tag = ERR
            elif line.startswith("  ") and line.strip()[:1].isdigit():
                tag = "step"
            self.report.insert("end", line + "\n", tag or ())
        self.report.configure(state="disabled")

    def _say(self, text: str, level: str = "muted") -> None:
        self.status.configure(text=text, foreground=self.COLORS[level])

    def save_log(self) -> bool:
        if self.analysis is None:
            self._say("目前的輸入有誤，無法存入紀錄", ERR)
            return False
        try:
            append_log(self.log_path, build_log_row(self.analysis, dt.datetime.now()))
        except OSError as exc:
            self._say(f"寫入失敗：{exc}", ERR)
            return False
        self._say(f"{dt.datetime.now():%H:%M:%S} 已存入 {self.log_path}", OK)
        return True

    def copy_report(self) -> None:
        if self.analysis is None:
            self._say("目前的輸入有誤，沒有結果可複製", ERR)
            return
        self.root.clipboard_clear()
        self.root.clipboard_append(render_report(self.analysis))
        self._say("已複製結果到剪貼簿", OK)

    def clear_counts(self) -> None:
        for key, value in (("a_live", ""), ("a_dead", "0"), ("b_live", ""), ("b_dead", "")):
            self.vars[key].set(value)
        self.a_live_entry.focus_set()


def run_gui(args: argparse.Namespace) -> int:
    if tk is None:
        print("這個 Python 沒有 tkinter，無法開啟圖形介面。請改用命令列模式。", file=sys.stderr)
        return 1
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        print(f"無法開啟視窗：{exc}", file=sys.stderr)
        return 1
    SeedingApp(root, form_defaults(args), args.log_file or default_log_path())
    root.minsize(820, 560)
    root.mainloop()
    return 0


def main(argv: Optional[List[str]] = None, input_fn: InputFn = input) -> int:
    for stream in (sys.stdout, sys.stderr):
        # 導向到檔案/管線時 Windows 預設是 cp950，µ、×、≤ 會變問號；改用 UTF-8。
        # 在終端機直接顯示則維持原樣（Python 走 Unicode console API，本來就正常）。
        try:
            if stream.isatty():
                stream.reconfigure(errors="replace")
            else:
                stream.reconfigure(encoding="utf-8")
        except (AttributeError, ValueError):
            pass
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.gui:
        return run_gui(args)
    try:
        if args.interactive or args.a_live is None:
            inputs = collect_interactive(args, input_fn)
        else:
            inputs = collect_batch(args)

        if args.mode == "harvest":
            plate = args.harvest_plate or "96"
            vol = args.harvest_vol_ul
            if vol is None:
                vol = HARVEST_PRESETS[plate].default_vol_ul
            result = analyze_harvest(
                inputs,
                plate_type=plate,
                harvest_vol_ul=vol,
                harvested_wells=args.harvest_wells,
                custom_area=args.custom_area,
                custom_wells=args.custom_wells,
                seeded_per_well=args.seeded_cells,
            )
            report_text = render_harvest_report(result)
        else:
            result = analyze(inputs)
            report_text = render_report(result)
    except (EOFError, KeyboardInterrupt):
        print("\n已取消。")
        return 130
    except ValueError as exc:
        parser.error(str(exc))
        return 2  # parser.error 已經 exit，這行只是讓型別檢查安心

    print()
    print(report_text)


    if not args.no_log:
        path = args.log_file or default_log_path()
        try:
            append_log(path, build_log_row(result, dt.datetime.now()))
            print(f"\n（已附加到 {path}）")
        except OSError as exc:
            print(f"\n[!] 無法寫入紀錄檔 {path}：{exc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
