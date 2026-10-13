# -*- coding: utf-8 -*-
"""cell_seeding_calc 的單元測試：python -m unittest test_cell_seeding_calc -v"""

import contextlib
import csv
import datetime as dt
import io
import os
import tempfile
import unittest

import cell_seeding_calc as c


def make_inputs(**overrides):
    base = dict(
        tube_ml=2.0, pre_cell_ul=10, pre_diluent_ul=0, tb_cell_ul=10, tb_dye_ul=10,
        grids=4, a_live=200, a_dead=8, b_live=208, b_dead=10,
        preset="96d1", cells_per_well=10000, well_ul=100, wells=39, margin_pct=25,
    )
    base.update(overrides)
    return c.Inputs(**base)


def preset_kwargs(name):
    p = c.PRESETS[name]
    return dict(preset=name, cells_per_well=p.cells_per_well, well_ul=p.well_ul, wells=p.wells)


class SpecCases(unittest.TestCase):
    """使用者指定的三組測試案例。"""

    def test_case1_96d1_single_step(self):
        a = c.analyze(make_inputs())
        self.assertAlmostEqual(a.df, 2.0)
        self.assertAlmostEqual(a.stock_conc, 1.02e6, delta=1)
        self.assertAlmostEqual(a.diff_pct, 3.92, places=2)
        self.assertAlmostEqual(a.viability_pct, 95.77, places=2)
        self.assertEqual(a.plan.total_ul, 5000)
        self.assertAlmostEqual(a.plan.v1_ul, 490, delta=1)
        self.assertAlmostEqual(a.plan.dmem_ul, 4510, delta=1)
        self.assertIsNone(a.plan.two_step)
        self.assertFalse([f for f in a.findings if f.level != c.OK])

    def test_case2_24tw(self):
        a = c.analyze(make_inputs(**preset_kwargs("24tw")))
        self.assertAlmostEqual(a.plan.target_conc, 1e4)
        self.assertEqual(a.plan.total_ul, 7500)
        self.assertAlmostEqual(a.plan.v1_ul, 73.5, delta=0.1)
        self.assertIsNone(a.plan.two_step)

    def test_case3_two_step(self):
        plan = c.plan_seeding(1.2e7, 10000, 100, 39, 25)
        self.assertAlmostEqual(plan.v1_ul, 41.7, delta=0.05)
        t = plan.two_step
        self.assertIsNotNone(t)
        self.assertEqual(t.factor, 5)
        self.assertEqual(t.intermediate_ul, 1000)
        self.assertEqual(t.stock_ul, 200)
        self.assertEqual(t.diluent1_ul, 800)
        self.assertAlmostEqual(t.transfer_ul, 208, delta=0.5)
        self.assertAlmostEqual(t.diluent2_ul, 4792, delta=0.5)
        # 物料守恆：第二段吸的中間液帶了多少細胞 = 需要的細胞數
        self.assertAlmostEqual(t.transfer_ul * (1.2e7 / t.factor) / 1000, plan.cells_needed, delta=1)


class Formulas(unittest.TestCase):
    def test_dilution_factors(self):
        self.assertEqual(c.dilution_factor(10, 0), 1)
        self.assertEqual(c.dilution_factor(10, 10), 2)
        self.assertEqual(c.total_dilution(10, 30, 10, 10), 8)

    def test_chamber_concentration(self):
        self.assertAlmostEqual(c.chamber_concentration(200, 4, 2), 1.0e6)
        self.assertAlmostEqual(c.chamber_concentration(90, 9, 2), 2.0e5)

    def test_quick_k_matches_mean_concentration(self):
        a = c.analyze(make_inputs())
        self.assertAlmostEqual(a.k, 2500)
        self.assertAlmostEqual((200 + 208) * a.k, a.stock_conc, delta=1e-6)

    def test_prepare_volume_rounds_up_to_500(self):
        self.assertEqual(c.prepare_volume(39, 100, 25), 5000)   # 4875
        self.assertEqual(c.prepare_volume(40, 100, 25), 5000)   # 剛好 5000 不能多進位
        self.assertEqual(c.prepare_volume(41, 100, 25), 5500)   # 5125
        self.assertEqual(c.prepare_volume(1, 10, 0), 500)

    def test_ceil_ignores_float_noise(self):
        self.assertEqual(c.ceil_to_multiple(5000.0000000001, 500), 5000)

    def test_round_half_up(self):
        self.assertEqual([c.round_half_up(x) for x in (72.5, 73.5, 0.4, 490.2)], [73, 74, 0, 490])

    def test_scientific_notation(self):
        self.assertEqual(c.sci(1.02e6), "1.02×10^6")
        self.assertEqual(c.sci(1e5), "1×10^5")
        self.assertEqual(c.sci(999999), "1×10^6")
        self.assertEqual(c.sci(1.2e7), "1.2×10^7")
        self.assertEqual(c.sci(0), "0")

    def test_single_chamber_uses_upper_only(self):
        a = c.analyze(make_inputs(b_live=None, b_dead=None))
        self.assertAlmostEqual(a.stock_conc, 1.0e6)
        self.assertIsNone(a.diff_pct)
        self.assertAlmostEqual(a.viability_pct, 200 / 208 * 100)
        self.assertAlmostEqual(a.k, 5000)


class TwoStepRules(unittest.TestCase):
    def test_threshold_is_50(self):
        # V1 = Ct*Vt/C；C 取 5e5 → 1e5*5000/5e5 = 1000；取 1e7 → 50（不改兩段）
        self.assertIsNone(c.plan_seeding(1e7, 10000, 100, 39, 25).two_step)
        self.assertIsNotNone(c.plan_seeding(1.0001e7, 10000, 100, 39, 25).two_step)

    def test_factor_selection(self):
        cases = {45: 5, 30: 5, 20: 5, 19: 10, 9: 20, 4: 50}   # V1 → f
        for v1, expected in cases.items():
            self.assertEqual(c.two_step_dilution(v1, 5000).factor, expected, v1)
        self.assertEqual(c.two_step_dilution(49, 5000).factor, 5)   # 49*2=98 < 100

    def test_factor_exhausted_flags_transfer(self):
        t = c.two_step_dilution(1.0, 5000)
        self.assertEqual(t.factor, 50)
        self.assertFalse(t.transfer_ok)

    def test_intermediate_volume_covers_transfer(self):
        for v1 in (1.5, 10, 25, 49.9):
            t = c.two_step_dilution(v1, 5000)
            self.assertGreaterEqual(t.intermediate_ul, t.transfer_ul + 100 - 1e-9)
            self.assertAlmostEqual(t.stock_ul + t.diluent1_ul, t.intermediate_ul)
            self.assertAlmostEqual(t.transfer_ul + t.diluent2_ul, 5000)


class Validation(unittest.TestCase):
    def test_negative_rejected(self):
        for field in ("a_live", "a_dead", "b_live", "b_dead", "tube_ml", "pre_diluent_ul",
                      "tb_dye_ul", "cells_per_well", "well_ul", "wells", "margin_pct"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                make_inputs(**{field: -1})

    def test_grids_only_4_or_9(self):
        for bad in (0, 1, 3, 5, 16):
            with self.subTest(grids=bad), self.assertRaises(ValueError):
                make_inputs(grids=bad)
        for good in (4, 9):
            make_inputs(grids=good)

    def test_dilution_factor_below_one_rejected(self):
        with self.assertRaises(ValueError):
            c.chamber_concentration(100, 4, 0.5)
        with self.assertRaises(ValueError):
            c.dilution_factor(0, 10)       # 細胞液 0 → 分母為 0
        with self.assertRaises(ValueError):
            c.dilution_factor(10, -1)      # 會讓倍數 < 1

    def test_half_filled_lower_chamber_rejected(self):
        with self.assertRaises(ValueError):
            make_inputs(b_live=200, b_dead=None)

    def test_no_cells_at_all_rejected(self):
        with self.assertRaises(ValueError):
            make_inputs(a_live=0, a_dead=0, b_live=None, b_dead=None)

    def test_parsers(self):
        self.assertEqual(c.parse_grids("9"), 9)
        self.assertEqual(c.parse_wells("39"), 39)
        self.assertEqual(c.parse_preset(" 96D1 "), "96d1")
        for parser, bad in ((c.parse_grids, "5"), (c.parse_wells, "2.5"), (c.parse_preset, "48w"),
                            (c.number_parser("x"), "abc"), (c.number_parser("x"), "-3"),
                            (c.number_parser("x"), "nan"), (c.number_parser("x", True), "0")):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                parser(bad)


class Findings(unittest.TestCase):
    def levels(self, **overrides):
        return [f.level for f in c.analyze(make_inputs(**overrides)).findings]

    def messages(self, level, **overrides):
        return " ".join(f.message for f in c.analyze(make_inputs(**overrides)).findings if f.level == level)

    def test_too_dense(self):
        self.assertIn("太濃", self.messages(c.WARN, a_live=480, b_live=470))

    def test_too_sparse(self):
        self.assertIn("太稀", self.messages(c.WARN, a_live=60, a_dead=2, b_live=62, b_dead=2))

    def test_sparse_hint_with_9_grids_does_not_suggest_9(self):
        msg = self.messages(c.WARN, grids=9, a_live=100, a_dead=3, b_live=100, b_dead=3)
        self.assertIn("太稀", msg)
        self.assertNotIn("改數 9 格", msg)

    def test_difference_over_20_is_error(self):
        self.assertIn("沒打散", self.messages(c.ERR, a_live=150, b_live=250))

    def test_viability_boundary(self):
        self.assertNotIn("存活率", self.messages(c.ERR, a_live=180, a_dead=0, b_live=180, b_dead=0))
        # 剛好 90% → 錯誤（≤ 90）
        self.assertIn("存活率", self.messages(c.ERR, a_live=180, a_dead=20, b_live=180, b_dead=20))

    def test_not_enough_cells(self):
        # 原液 ≈ 1e5/mL → 整管只有 2e5，需要 5e5
        a = c.analyze(make_inputs(a_live=20, b_live=20, a_dead=0, b_dead=0, tube_ml=2))
        self.assertTrue(any("細胞不夠" in f.message for f in a.findings if f.level == c.ERR))

    def test_v1_greater_than_vt_is_error(self):
        a = c.analyze(make_inputs(a_live=5, b_live=5, a_dead=0, b_dead=0, tube_ml=100))
        self.assertTrue(any("V1" in f.message and "還稀" in f.message for f in a.findings if f.level == c.ERR))
        self.assertFalse(a.plan.feasible)

    def test_v1_over_60_percent_is_warning_not_error(self):
        # C = 1.2e5/mL → V1 = 1e5*5000/1.2e5 ≈ 4167 (> 3000 但 < 5000)
        a = c.analyze(make_inputs(a_live=24, b_live=24, a_dead=0, b_dead=0, tube_ml=100))
        self.assertTrue(any("偏稀" in f.message for f in a.findings if f.level == c.WARN))
        self.assertFalse(any("V1" in f.message for f in a.findings if f.level == c.ERR))
        self.assertTrue(a.plan.feasible)

    def test_two_step_is_reported(self):
        a = c.analyze(make_inputs(a_live=5000, b_live=5000, a_dead=0, b_dead=0, grids=9, tube_ml=2))
        self.assertIsNotNone(a.plan.two_step)
        self.assertTrue(any("兩段稀釋" in f.message for f in a.findings))

    def test_zero_dead_is_flagged(self):
        self.assertIn("死細胞總數為 0", self.messages(c.WARN, a_dead=0, b_dead=0))

    def test_zero_live_gives_error_not_crash(self):
        a = c.analyze(make_inputs(a_live=0, b_live=0))
        self.assertIsNone(a.plan)
        self.assertTrue(a.has_error)
        self.assertIn("無法配製", c.render_report(a))


class Report(unittest.TestCase):
    def test_single_step_report(self):
        text = c.render_report(c.analyze(make_inputs()))
        self.assertIn("1.02×10^6", text)
        self.assertIn("3.92 %", text)
        self.assertIn("95.77 %", text)
        self.assertIn("取原液 490 µL", text)
        self.assertIn("加 DMEM 4510 µL", text)
        self.assertIn("[OK]", text)
        self.assertIn("(200 + 208) × 2500", text)
        last_step = [l for l in text.splitlines() if l.strip().startswith(("1.", "2.", "3.", "4.", "5."))][-1]
        self.assertIn("每孔加 100 µL，每加 3 孔回管子混勻一次", last_step)

    def test_two_step_report_steps(self):
        a = c.analyze(make_inputs(a_live=5000, b_live=5000, a_dead=0, b_dead=0, grids=9))
        steps = c.recipe_steps(a.inputs, a.plan)
        self.assertEqual(len(steps), 5)
        self.assertIn("每孔加 100 µL，每加 3 孔回管子混勻一次", steps[-1])
        # 第二段「中間液 + DMEM」加起來等於總體積
        transfer = int(steps[2].split("取中間稀釋液 ")[1].split(" µL")[0])
        dmem = int(steps[3].split("加 DMEM ")[1].split(" µL")[0])
        self.assertEqual(transfer + dmem, 5000)

    def test_symbols(self):
        a = c.analyze(make_inputs(a_live=150, b_live=250))
        text = c.render_report(a)
        self.assertIn("[X]", text)
        self.assertIn("沒打散", text)


class Logging(unittest.TestCase):
    def test_creates_file_with_header_then_appends(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "seeding_log.csv")
            a = c.analyze(make_inputs())
            when = dt.datetime(2026, 10, 11, 9, 30, 0)
            c.append_log(path, c.build_log_row(a, when))
            c.append_log(path, c.build_log_row(c.analyze(make_inputs(b_live=None, b_dead=None)), when))
            with open(path, encoding="utf-8-sig", newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 2)
            self.assertEqual(list(rows[0].keys()), c.LOG_COLUMNS)
            self.assertEqual(rows[0]["日期時間"], "2026-10-11 09:30:00")
            self.assertEqual(rows[0]["上槽活"], "200")
            self.assertEqual(rows[0]["存活率_%"], "95.77")
            self.assertEqual(float(rows[0]["原液濃度C_per_mL"]), 1.02e6)
            self.assertEqual(rows[0]["稀釋方式"], "一段")
            self.assertEqual(rows[1]["下槽活"], "")
            with open(path, "rb") as fh:       # BOM 只在檔案開頭出現一次
                self.assertEqual(fh.read().count(b"\xef\xbb\xbf"), 1)

    def test_empty_existing_file_gets_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "seeding_log.csv")
            open(path, "w").close()
            c.append_log(path, c.build_log_row(c.analyze(make_inputs()), dt.datetime(2026, 1, 1)))
            with open(path, encoding="utf-8-sig") as fh:
                self.assertEqual(fh.readline().strip().split(",")[0], "日期時間")


class Cli(unittest.TestCase):
    def run_main(self, argv, inputs=()):
        feed = iter(inputs)
        out, err = io.StringIO(), io.StringIO()

        def fake_input(prompt):
            out.write(prompt)
            try:
                return next(feed)
            except StopIteration:
                raise EOFError from None

        code = None
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                code = c.main(argv, input_fn=fake_input)
            except SystemExit as exc:
                code = exc.code
        return code, out.getvalue(), err.getvalue()

    def test_one_liner_from_spec(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = os.path.join(tmp, "log.csv")
            code, out, _ = self.run_main(
                ["--a-live", "200", "--a-dead", "8", "--b-live", "208", "--b-dead", "10",
                 "--preset", "96d1", "--log-file", log])
            self.assertEqual(code, 0)
            self.assertIn("1.02×10^6", out)
            self.assertTrue(os.path.exists(log))

    def test_spec_example_without_dead_counts_works_and_warns(self):
        code, out, _ = self.run_main(["--a-live", "200", "--b-live", "208", "--preset", "96d1", "--no-log"])
        self.assertEqual(code, 0)
        self.assertIn("死細胞總數為 0", out)

    def test_invalid_argument_exits_2(self):
        code, _, err = self.run_main(["--a-live", "-5", "--no-log"])
        self.assertEqual(code, 2)
        self.assertIn("不能是負數", err)
        code, _, err = self.run_main(["--a-live", "10", "--grids", "5", "--no-log"])
        self.assertEqual(code, 2)
        self.assertIn("4 或 9", err)

    def test_custom_requires_all_three(self):
        code, _, err = self.run_main(["--a-live", "200", "--preset", "custom", "--wells", "10", "--no-log"])
        self.assertEqual(code, 2)
        self.assertIn("custom", err)

    def test_custom_by_overrides_alone(self):
        code, out, _ = self.run_main(["--a-live", "200", "--b-live", "208", "--cells-per-well", "5000",
                                      "--well-ul", "200", "--wells", "20", "--no-log"])
        self.assertEqual(code, 0)
        self.assertIn("custom", out)
        self.assertIn("Ct = 2.5×10^4", out)

    def test_interactive_all_defaults_except_counts(self):
        answers = ["",              # DMEM mL
                   "", "",          # 預稀釋
                   "", "",          # trypan
                   "",              # 格數
                   "200", "8",      # 上槽
                   "208", "10",     # 下槽
                   "",              # preset → 96d1
                   ""]              # 餘量
        with tempfile.TemporaryDirectory() as tmp:
            code, out, _ = self.run_main(["--log-file", os.path.join(tmp, "l.csv")], answers)
        self.assertEqual(code, 0)
        self.assertIn("[預設 2]", out)
        self.assertIn("[預設 96d1]", out)
        self.assertIn("原液濃度 C：1.02×10^6", out)

    def test_interactive_reprompts_on_invalid_input(self):
        answers = ["", "", "", "", "", "7",   # 格數 7 → 被拒絕
                   "4",
                   "-1", "200", "8",           # 上槽活 -1 → 被拒絕
                   "", "",                     # 下槽留空
                   "24tw", ""]
        code, out, _ = self.run_main(["--no-log"], answers)
        self.assertEqual(code, 0)
        self.assertIn("格數只能是 4 或 9", out)
        self.assertIn("不能是負數", out)
        self.assertIn("只用上槽", out)

    def test_interactive_custom_preset(self):
        answers = ["", "", "", "", "", "", "200", "8", "", "custom", "3000", "150", "24", ""]
        code, out, _ = self.run_main(["--no-log"], answers)
        self.assertEqual(code, 0)
        self.assertIn("每孔 3000 cells / 150 µL", out)

    def test_eof_cancels_cleanly(self):
        code, out, _ = self.run_main(["--no-log"], [])
        self.assertEqual(code, 130)
        self.assertIn("已取消", out)

    def test_all_presets_match_spec(self):
        self.assertEqual({k: tuple(v) for k, v in c.PRESETS.items()}, {
            "96d1": (10000, 100, 39), "96d3": (4000, 100, 30), "96d5": (1000, 100, 30),
            "24tw": (5000, 500, 12), "24d1": (20000, 500, 12)})


class GuiForm(unittest.TestCase):
    """GUI 背後的純函式：表單字串 → Inputs，與一行摘要。"""

    def form(self, **overrides):
        values = c.form_defaults()
        values.update(a_live="200", a_dead="8", b_live="208", b_dead="10")
        values.update(overrides)
        return values

    def test_form_matches_spec_case1(self):
        a = c.analyze(c.inputs_from_form(self.form()))
        self.assertAlmostEqual(a.stock_conc, 1.02e6, delta=1)
        self.assertEqual(a.plan.total_ul, 5000)

    def test_defaults_ask_for_upper_live(self):
        with self.assertRaisesRegex(ValueError, "上槽活細胞數"):
            c.inputs_from_form(c.form_defaults())

    def test_blank_lower_chamber_means_upper_only(self):
        inp = c.inputs_from_form(self.form(b_live="", b_dead=""))
        self.assertIsNone(inp.b_live)

    def test_lower_dead_blank_is_zero_but_dead_without_live_is_rejected(self):
        self.assertEqual(c.inputs_from_form(self.form(b_dead="")).b_dead, 0.0)
        with self.assertRaises(ValueError):
            c.inputs_from_form(self.form(b_live="", b_dead="5"))

    def test_bad_values_give_readable_errors(self):
        for field, value, message in (("a_live", "-3", "不能是負數"), ("grids", "5", "4 或 9"),
                                      ("a_live", "abc", "要輸入數字"), ("tb_cell", "0", "大於 0")):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                c.inputs_from_form(self.form(**{field: value}))
        with self.assertRaisesRegex(ValueError, "整數"):          # 孔數只在 custom 時才會被讀取
            c.inputs_from_form(self.form(preset="custom", wells="2.5"))

    def test_custom_preset_uses_typed_values(self):
        inp = c.inputs_from_form(self.form(preset="custom", cells_per_well="3000", well_ul="150", wells="24"))
        self.assertEqual((inp.cells_per_well, inp.well_ul, inp.wells), (3000, 150, 24))

    def test_typed_cells_per_well_wins_and_marks_custom(self):
        inp = c.inputs_from_form(self.form(preset="96d1", cells_per_well="8000"))
        self.assertEqual((inp.cells_per_well, inp.well_ul, inp.wells), (8000, 100, 39))
        self.assertEqual(inp.preset, "custom")
        same = c.inputs_from_form(self.form(preset="24tw", cells_per_well="5000", well_ul="500", wells="12"))
        self.assertEqual(same.preset, "24tw")

    def test_defaults_from_cli_args_override_to_custom(self):
        args = c.build_parser().parse_args(["--preset", "96d1", "--wells", "20", "--a-live", "150"])
        values = c.form_defaults(args)
        self.assertEqual((values["preset"], values["wells"], values["a_live"]), ("custom", "20", "150"))

    def test_summary_line(self):
        level, text = c.summary_line(c.analyze(make_inputs()))
        self.assertEqual(level, c.OK)
        self.assertIn("取原液 490 µL", text)
        self.assertIn("DMEM 4510 µL", text)
        level, text = c.summary_line(c.analyze(make_inputs(a_live=150, b_live=250)))
        self.assertEqual(level, c.ERR)
        level, text = c.summary_line(c.analyze(make_inputs(a_live=5000, b_live=5000, grids=9)))
        self.assertEqual(level, c.WARN)
        self.assertIn("兩段稀釋", text)


@unittest.skipIf(c.tk is None, "沒有 tkinter")
class GuiWindow(unittest.TestCase):
    def setUp(self):
        try:
            self.root = c.tk.Tk()
        except c.tk.TclError as exc:
            self.skipTest(f"沒有可用的顯示環境：{exc}")
        self.root.withdraw()
        self.tmp = tempfile.TemporaryDirectory()
        self.log = os.path.join(self.tmp.name, "log.csv")
        self.app = c.SeedingApp(self.root, c.form_defaults(), self.log)

    def tearDown(self):
        self.root.destroy()
        self.tmp.cleanup()

    def set_counts(self, a_live, a_dead, b_live, b_dead):
        for key, value in (("a_live", a_live), ("a_dead", a_dead), ("b_live", b_live), ("b_dead", b_dead)):
            self.app.vars[key].set(value)

    def test_starts_asking_for_input(self):
        self.assertIsNone(self.app.analysis)
        self.assertIn("上槽活細胞數", self.app.headline.cget("text"))

    def test_live_recalculation(self):
        self.set_counts("200", "8", "208", "10")
        self.assertAlmostEqual(self.app.analysis.stock_conc, 1.02e6, delta=1)
        self.assertIn("取原液 490 µL", self.app.headline.cget("text"))
        self.assertIn("1.02×10^6", self.app.report.get("1.0", "end"))
        self.app.vars["preset"].set("24tw")           # 換 preset：欄位同步、結果重算
        self.assertEqual(self.app.vars["wells"].get(), "12")
        self.assertEqual(self.app.analysis.plan.total_ul, 7500)

    def test_dilution_factor_labels_follow_inputs(self):
        self.assertIn("1 × 2 ＝ 2", self.app.df_label.cget("text"))
        self.app.vars["pre_dmem"].set("10")
        self.assertIn("＝ ×2", self.app.pre_factor.cget("text"))
        self.assertIn("2 × 2 ＝ 4", self.app.df_label.cget("text"))
        self.app.vars["tb_cell"].set("")
        self.assertIn("?", self.app.tb_factor.cget("text"))
        self.assertIn("?", self.app.df_label.cget("text"))

    def test_cells_per_well_always_editable_and_flips_to_custom(self):
        self.set_counts("200", "8", "208", "10")
        entry = self.app.setup_entries[0]
        self.assertEqual(str(entry.cget("state")), "normal")
        self.app.vars["cells_per_well"].set("20000")
        self.assertEqual(self.app.vars["preset"].get(), "custom")
        self.assertEqual(self.app.analysis.plan.target_conc, 2e5)
        self.assertEqual(self.app.analysis.inputs.preset, "custom")
        self.app.vars["preset"].set("96d1")           # 重新選 preset → 數值被帶回
        self.assertEqual(self.app.vars["cells_per_well"].get(), "10000")
        self.assertEqual(self.app.vars["preset"].get(), "96d1")

    def test_save_and_clear(self):
        self.assertFalse(self.app.save_log())         # 輸入不完整不能存
        self.set_counts("200", "8", "208", "10")
        self.assertTrue(self.app.save_log())
        with open(self.log, encoding="utf-8-sig", newline="") as fh:
            self.assertEqual(len(list(csv.DictReader(fh))), 1)
        self.app.clear_counts()
        self.assertEqual(self.app.vars["a_live"].get(), "")
        self.assertIsNone(self.app.analysis)


class HarvestCases(unittest.TestCase):
    """切孔／收孔細胞量計算測試。"""

    def test_harvest_96_well(self):
        h = c.plan_harvest(4.05e5, 95.0, 100.0, 1, "96", seeded_per_well=10000.0)
        self.assertAlmostEqual(h.total_harvest_live, 40500.0, delta=1.0)
        self.assertAlmostEqual(h.cells_per_well, 40500.0, delta=1.0)
        self.assertAlmostEqual(h.density, 40500.0 / 0.32, delta=1.0)
        self.assertAlmostEqual(h.plate_total, 40500.0 * 96, delta=1.0)
        self.assertAlmostEqual(h.fold, 4.05, delta=0.01)
        self.assertAlmostEqual(h.doublings, 2.018, delta=0.01)

    def test_harvest_24_well(self):
        h = c.plan_harvest(5.0e5, 96.0, 500.0, 1, "24", seeded_per_well=50000.0)
        self.assertAlmostEqual(h.total_harvest_live, 250000.0, delta=1.0)
        self.assertAlmostEqual(h.cells_per_well, 250000.0, delta=1.0)
        self.assertAlmostEqual(h.density, 250000.0 / 1.9, delta=1.0)
        self.assertAlmostEqual(h.plate_total, 250000.0 * 24, delta=1.0)
        self.assertAlmostEqual(h.fold, 5.0, delta=0.01)

    def test_harvest_pooled_wells(self):
        h = c.plan_harvest(4.0e5, 95.0, 200.0, 2, "96")
        self.assertAlmostEqual(h.total_harvest_live, 80000.0, delta=1.0)
        self.assertAlmostEqual(h.cells_per_well, 40000.0, delta=1.0)

    def test_harvest_cli_mode(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            ret = c.main([
                "--mode", "harvest",
                "--a-live", "80", "--a-dead", "4",
                "--harvest-plate", "96", "--harvest-vol-ul", "100",
                "--seeded-cells", "10000",
                "--no-log"
            ])
        self.assertEqual(ret, 0)
        text = out.getvalue()
        self.assertIn("=== 收孔細胞量分析 ===", text)
        self.assertIn("96 孔盤", text)
        self.assertIn("單孔活細胞數", text)

    def test_harvest_log_row(self):
        inp = make_inputs(a_live=80, a_dead=4, b_live=82, b_dead=4)
        a = c.analyze_harvest(inp, plate_type="96", harvest_vol_ul=100.0, harvested_wells=1, seeded_per_well=10000.0)
        row = c.build_log_row(a, dt.datetime(2026, 10, 14, 12, 0, 0))
        self.assertEqual(row["實驗設定"], "收孔:96 孔盤")
        self.assertEqual(row["每孔細胞數"], 40500)
        self.assertEqual(row["稀釋方式"], "收孔分析")


if __name__ == "__main__":
    unittest.main()

