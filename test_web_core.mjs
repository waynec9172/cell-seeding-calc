// 測試 web/index.html 內的計算核心：node --test test_web_core.mjs
// 只抽出 <script> 中 "core:start" 到 "core:end" 之間的純函式來測試（不需要瀏覽器）。
import { readFileSync } from "node:fs";
import { test } from "node:test";
import assert from "node:assert/strict";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const html = readFileSync(join(dirname(fileURLToPath(import.meta.url)), "web", "index.html"), "utf8");
const m = html.match(/\/\/ ==== core:start ====[^\n]*\n([\s\S]*?)\/\/ ==== core:end ====/);
assert.ok(m, "找不到 core:start / core:end 標記");
const C = new Function(`"use strict";\n${m[1]}\nreturn CORE;`)();

const base = (o = {}) => ({
  tubeMl: 2, preCell: 10, preDil: 0, tbCell: 10, tbDye: 10, grids: 4,
  aLive: 200, aDead: 8, bLive: 208, bDead: 10,
  preset: "96d1", cellsPerWell: 10000, wellUl: 100, wells: 39, marginPct: 25, ...o,
});
const near = (a, b, d = 1e-6) => assert.ok(Math.abs(a - b) <= d, `${a} ≉ ${b} (±${d})`);

test("案例 1：96d1 一段稀釋", () => {
  const a = C.analyze(base());
  near(a.df, 2);
  near(a.stock, 1.02e6, 1);
  near(a.diffPct, 3.92, 0.005);
  near(a.viability, 95.77, 0.005);
  assert.equal(a.plan.totalUl, 5000);
  near(a.plan.v1, 490, 1);
  near(a.plan.dmem, 4510, 1);
  assert.equal(a.plan.twoStep, null);
  assert.deepEqual(a.findings.filter((f) => f.level !== "OK"), []);
});

test("案例 2：24tw", () => {
  const a = C.analyze(base({ preset: "24tw", cellsPerWell: 5000, wellUl: 500, wells: 12 }));
  near(a.plan.targetConc, 1e4);
  assert.equal(a.plan.totalUl, 7500);
  near(a.plan.v1, 73.5, 0.1);
  assert.equal(a.plan.twoStep, null);
});

test("案例 3：原液 1.2×10^7 → 兩段稀釋", () => {
  const p = C.planSeeding(1.2e7, 10000, 100, 39, 25);
  near(p.v1, 41.7, 0.05);
  const t = p.twoStep;
  assert.equal(t.factor, 5);
  assert.equal(t.intermediate, 1000);
  assert.equal(t.stock, 200);
  assert.equal(t.diluent1, 800);
  near(t.transfer, 208, 0.5);
  near(t.diluent2, 4792, 0.5);
  near((t.transfer * (1.2e7 / t.factor)) / 1000, p.cellsNeeded, 1);
});

test("數字工具", () => {
  assert.deepEqual([72.5, 73.5, 0.4, 490.2].map(C.roundHalfUp), [73, 74, 0, 490]);
  assert.equal(C.ceilToMultiple(4875, 500), 5000);
  assert.equal(C.ceilToMultiple(5000.0000000001, 500), 5000);
  assert.equal(C.prepareVolume(40, 100, 25), 5000);
  assert.equal(C.prepareVolume(41, 100, 25), 5500);
  assert.equal(C.sci(1.02e6), "1.02×10^6");
  assert.equal(C.sci(1e5), "1×10^5");
  assert.equal(C.sci(999999), "1×10^6");
  assert.equal(C.sci(1.2e7), "1.2×10^7");
  assert.equal(C.sci(0), "0");
  near(C.quickK(2, 4), 2500);
});

test("兩段稀釋的倍數選擇與門檻", () => {
  for (const [v1, f] of Object.entries({ 45: 5, 30: 5, 20: 5, 19: 10, 9: 20, 4: 50, 49: 5 }))
    assert.equal(C.twoStepDilution(Number(v1), 5000).factor, f, `V1=${v1}`);
  assert.equal(C.planSeeding(1e7, 10000, 100, 39, 25).twoStep, null);            // V1 = 50 不改兩段
  assert.ok(C.planSeeding(1.0001e7, 10000, 100, 39, 25).twoStep);
  const t = C.twoStepDilution(1, 5000);
  assert.equal(t.factor, 50);
  assert.equal(t.transferOk, false);
});

test("驗證：負數、格數、稀釋倍數", () => {
  for (const k of ["aLive", "aDead", "bLive", "bDead", "tubeMl", "preDil", "tbDye", "cellsPerWell", "wellUl", "wells", "marginPct"])
    assert.throws(() => C.analyze(base({ [k]: -1 })), /負數|大於 0/, k);
  for (const g of [0, 3, 5, 16]) assert.throws(() => C.analyze(base({ grids: g })), /4 或 9/);
  assert.throws(() => C.chamberConcentration(100, 4, 0.5), /不能 < 1/);
  assert.throws(() => C.dilutionFactor(0, 10), /大於 0/);
  assert.throws(() => C.dilutionFactor(10, -1), /負數/);
  assert.throws(() => C.analyze(base({ bLive: 200, bDead: null })), /一起輸入/);
  assert.throws(() => C.analyze(base({ aLive: 0, aDead: 0, bLive: null, bDead: null })), /沒有輸入/);
});

test("判斷：太濃、太稀、差異、存活率、細胞不夠、V1", () => {
  const msgs = (o, level) => C.analyze(base(o)).findings.filter((f) => f.level === level).map((f) => f.message).join(" ");
  assert.match(msgs({ aLive: 480, bLive: 470 }, "!"), /太濃/);
  assert.match(msgs({ aLive: 60, aDead: 2, bLive: 62, bDead: 2 }, "!"), /太稀/);
  assert.doesNotMatch(msgs({ grids: 9, aLive: 100, aDead: 3, bLive: 100, bDead: 3 }, "!"), /改數 9 格/);
  assert.match(msgs({ aLive: 150, bLive: 250 }, "X"), /沒打散/);
  assert.doesNotMatch(msgs({ aLive: 180, aDead: 0, bLive: 180, bDead: 0 }, "X"), /存活率/);
  assert.match(msgs({ aLive: 180, aDead: 20, bLive: 180, bDead: 20 }, "X"), /存活率/);           // 剛好 90% → 錯誤
  assert.match(msgs({ aLive: 20, bLive: 20, aDead: 0, bDead: 0 }, "X"), /細胞不夠/);
  const big = C.analyze(base({ aLive: 5, bLive: 5, aDead: 0, bDead: 0, tubeMl: 100 }));
  assert.ok(big.findings.some((f) => f.level === "X" && /V1/.test(f.message) && /還稀/.test(f.message)));
  assert.equal(big.plan.feasible, false);
  const thin = C.analyze(base({ aLive: 24, bLive: 24, aDead: 0, bDead: 0, tubeMl: 100 }));
  assert.ok(thin.findings.some((f) => f.level === "!" && /偏稀/.test(f.message)));
  assert.equal(thin.plan.feasible, true);
  assert.match(msgs({ aDead: 0, bDead: 0 }, "!"), /死細胞總數為 0/);
});

test("零活細胞不會當機", () => {
  const a = C.analyze(base({ aLive: 0, bLive: 0 }));
  assert.equal(a.plan, null);
  assert.equal(a.hasError, true);
  assert.match(C.renderReport(a), /無法配製/);
});

test("只用上槽", () => {
  const a = C.analyze(base({ bLive: null, bDead: null }));
  near(a.stock, 1.0e6, 1e-6);
  assert.equal(a.diffPct, null);
  near(a.k, 5000);
});

test("報告與配方文字", () => {
  const a = C.analyze(base());
  const text = C.renderReport(a);
  for (const s of ["1.02×10^6", "3.92 %", "95.77 %", "取原液 490 µL", "加 DMEM 4510 µL", "[OK]", "(200 + 208) × 2500"])
    assert.ok(text.includes(s), s);
  const steps = C.recipeSteps(a.inputs, a.plan);
  assert.equal(steps.length, 3);
  assert.ok(C.stepText(steps.at(-1)).startsWith("每孔加 100 µL，每加 3 孔回管子混勻一次"));

  const two = C.analyze(base({ aLive: 5000, bLive: 5000, aDead: 0, bDead: 0, grids: 9 }));
  const s2 = C.recipeSteps(two.inputs, two.plan);
  assert.equal(s2.length, 5);
  const tr = Number(C.stepText(s2[2]).match(/取中間稀釋液 (\d+) µL/)[1]);
  const dm = Number(C.stepText(s2[3]).match(/加 DMEM (\d+) µL/)[1]);
  assert.equal(tr + dm, 5000);

  const sm = C.summaryLine(a);
  assert.equal(sm.level, "OK");
  assert.match(sm.text, /取原液 490 µL ＋ DMEM 4510 µL/);
  assert.equal(C.summaryLine(C.analyze(base({ aLive: 150, bLive: 250 }))).level, "X");
  assert.equal(C.summaryLine(two).level, "!");
});

test("表單字串 → 輸入", () => {
  const f = (o = {}) => ({
    dmem_ml: "2", pre_cell: "10", pre_dmem: "0", tb_cell: "10", tb_dye: "10", grids: "4",
    a_live: "200", a_dead: "8", b_live: "208", b_dead: "10",
    preset: "96d1", cells_per_well: "10000", well_ul: "100", wells: "39", margin: "25", ...o,
  });
  const a = C.analyze(C.inputsFromForm(f()));
  near(a.stock, 1.02e6, 1);
  assert.equal(a.inputs.preset, "96d1");
  assert.throws(() => C.inputsFromForm(f({ a_live: "" })), /請輸入上槽活細胞數/);
  assert.equal(C.inputsFromForm(f({ b_live: "", b_dead: "" })).bLive, null);
  assert.equal(C.inputsFromForm(f({ b_dead: "" })).bDead, 0);
  assert.throws(() => C.inputsFromForm(f({ b_live: "", b_dead: "5" })), /下槽/);
  for (const [k, v, re] of [["a_live", "-3", /不能是負數/], ["grids", "5", /4 或 9/], ["a_live", "abc", /要輸入數字/],
                            ["wells", "2.5", /整數/], ["tb_cell", "0", /大於 0/]])
    assert.throws(() => C.inputsFromForm(f({ [k]: v })), re, `${k}=${v}`);
  // 自己填每孔細胞數 → 記為 custom
  const own = C.inputsFromForm(f({ cells_per_well: "8000" }));
  assert.equal(own.preset, "custom");
  assert.equal(own.cellsPerWell, 8000);
  assert.equal(C.inputsFromForm(f({ preset: "24tw", cells_per_well: "5000", well_ul: "500", wells: "12" })).preset, "24tw");
});

test("紀錄列與 CSV", () => {
  const a = C.analyze(base());
  const row = C.buildLogRow(a, new Date(2026, 9, 11, 9, 30, 0));
  assert.equal(row["日期時間"], "2026-10-11 09:30:00");
  assert.equal(row["上槽活"], 200);
  assert.equal(row["存活率_%"], 95.77);
  assert.equal(Number(row["原液濃度C_per_mL"]), 1.02e6);
  assert.equal(row["稀釋方式"], "一段");
  const only = C.buildLogRow(C.analyze(base({ bLive: null, bDead: null })), new Date());
  assert.equal(only["下槽活"], "");
  const csv = C.toCsv([row, only]);
  const lines = csv.trim().split("\r\n");
  assert.equal(lines.length, 3);
  assert.equal(lines[0].split(",").length, C.LOG_COLUMNS.length);
  assert.ok(lines[0].startsWith("日期時間,"));
  // 含逗號、引號的欄位要加引號並跳脫
  assert.ok(C.toCsv([{ "判斷": 'a,"b"' }]).includes('"a,""b"""'));
  const warn = C.buildLogRow(C.analyze(base({ aLive: 150, bLive: 250 })), new Date());
  assert.match(warn["判斷"], /^\[X\] 兩槽差異/);
});
