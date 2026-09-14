"""
make_derivation_doc.py
生成 IMM-KF → IMM-PF 数学推导 Word 文档
运行：  python make_derivation_doc.py
"""

from docx import Document
from docx.shared import Pt, Inches, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_ALIGN_VERTICAL
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import os

OUT_PATH = os.path.join(os.path.dirname(__file__), "IMM_KF_to_PF_Derivation.docx")

# ─────────────────────────── 样式辅助函数 ──────────────────────────────────

def set_cell_bg(cell, hex_color: str):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_color)
    tcPr.append(shd)

def set_cell_border(cell, sides=("top","bottom","left","right"), color="AAAAAA", sz="4"):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    tcBorders = OxmlElement("w:tcBorders")
    for side in sides:
        el = OxmlElement(f"w:{side}")
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), sz)
        el.set(qn("w:color"), color)
        tcBorders.append(el)
    tcPr.append(tcBorders)

def add_heading(doc, text, level=1):
    p = doc.add_heading(text, level=level)
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    return p

def add_para(doc, text="", bold=False, italic=False, size=11, indent=0, color=None, spacing_after=6):
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(spacing_after)
    p.paragraph_format.space_before = Pt(2)
    if indent:
        p.paragraph_format.left_indent = Inches(indent)
    if text:
        run = p.add_run(text)
        run.bold = bold
        run.italic = italic
        run.font.size = Pt(size)
        if color:
            run.font.color.rgb = RGBColor(*bytes.fromhex(color))
    return p

def add_equation(doc, eq_text, comment=""):
    """灰色背景公式块，仿 LaTeX 排版效果"""
    table = doc.add_table(rows=1, cols=1)
    table.style = "Table Grid"
    cell = table.rows[0].cells[0]
    set_cell_bg(cell, "F2F2F2")
    set_cell_border(cell, color="CCCCCC")
    para = cell.paragraphs[0]
    para.paragraph_format.space_before = Pt(4)
    para.paragraph_format.space_after  = Pt(4)
    para.paragraph_format.left_indent  = Inches(0.15)
    run = para.add_run(eq_text)
    run.font.name = "Cambria Math"
    run.font.size = Pt(11)
    if comment:
        run2 = para.add_run(f"        ← {comment}")
        run2.font.name = "Calibri"
        run2.font.size = Pt(9)
        run2.font.color.rgb = RGBColor(0x66, 0x66, 0x66)
    doc.add_paragraph().paragraph_format.space_after = Pt(4)
    return table

def add_step_box(doc, step_no, step_title, content_lines):
    """蓝色标题步骤框"""
    # 标题行
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(8)
    p.paragraph_format.space_after  = Pt(2)
    r = p.add_run(f"  Step {step_no}  {step_title}")
    r.bold = True
    r.font.size = Pt(11)
    r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    # 蓝色背景段落通过 XML 实现
    pPr = p._p.get_or_add_pPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear"); shd.set(qn("w:color"), "auto"); shd.set(qn("w:fill"), "2E75B6")
    pPr.append(shd)
    # 内容
    for line in content_lines:
        ip = doc.add_paragraph(line)
        ip.paragraph_format.left_indent = Inches(0.3)
        ip.paragraph_format.space_before = Pt(1)
        ip.paragraph_format.space_after  = Pt(1)
        ip.runs[0].font.size = Pt(10.5)

def make_sym_table(doc):
    headers = ["符号", "KF 含义", "PF 对应"]
    rows = [
        ["x = [pₓ, pᵧ, vₓ, vᵧ]ᵀ", "4 维状态向量", "单个粒子状态（同维）"],
        ["x̂ⱼ", "模型 j 后验均值", "粒子加权均值  x̂ⱼ = Σᵢ wⱼ⁽ⁱ⁾ xⱼ⁽ⁱ⁾"],
        ["Pⱼ", "模型 j 后验协方差", "粒子加权散度（隐式）"],
        ["{xⱼ⁽ⁱ⁾, wⱼ⁽ⁱ⁾}ᵢ₌₁ᴺ", "—（KF 无此量）", "N 粒子集，Σᵢ wⱼ⁽ⁱ⁾ = 1"],
        ["Λⱼ", "KF 归一化常数（高斯）", "粒子似然均值  Λⱼ = mean{ℓⱼ(i)}"],
        ["μⱼ", "IMM 后验模型概率", "完全不变（标量）"],
        ["c̄ⱼ", "IMM 预测模型概率", "完全不变（标量）"],
        ["πᵢⱼ", "Markov 转移概率", "完全不变"],
        ["β_{t,d}", "JPDA 关联概率", "完全不变（JPDA 输出）"],
        ["β_{t,0}", "漏检概率", "完全不变"],
        ["F, Q", "状态转移 / 过程噪声", "完全不变（用于粒子传播）"],
        ["H, R", "量测矩阵 / 噪声", "仅雷达使用；声学改为概率场"],
    ]
    t = doc.add_table(rows=1+len(rows), cols=3)
    t.style = "Table Grid"
    # 表头
    for i, h in enumerate(headers):
        c = t.rows[0].cells[i]
        set_cell_bg(c, "2E75B6")
        pr = c.paragraphs[0]
        pr.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = pr.add_run(h)
        r.bold = True; r.font.color.rgb = RGBColor(0xFF,0xFF,0xFF); r.font.size = Pt(10.5)
    # 内容
    for ri, row in enumerate(rows):
        for ci, val in enumerate(row):
            c = t.rows[ri+1].cells[ci]
            set_cell_bg(c, "EBF3FB" if ri % 2 == 0 else "FFFFFF")
            p = c.paragraphs[0]
            p.paragraph_format.left_indent = Inches(0.05)
            run = p.add_run(val)
            if ci == 0:
                run.font.name = "Cambria Math"; run.bold = True
            run.font.size = Pt(10)

def make_kf_to_pf_table(doc):
    headers = ["现有函数/字段", "KF 版本行为", "PF 替换/保留", "是否新建"]
    rows = [
        ["imm_mix.m",        "高斯参数混合(均值+协方差)", "imm_pf_mix.m — 粒子重混合",    "新建"],
        ["imm_predict.m",    "KF 线性预测",               "imm_pf_predict.m — 粒子传播", "新建"],
        ["ekf_pda_update.m", "KF PDA 更新",               "pf_jpda_update.m — 粒子权重更新", "新建"],
        ["imm_fuse.m",       "加权均值/协方差融合",        "完全不变（接收粒子均值）",     "保留"],
        ["imm_update_mu.m",  "模型概率贝叶斯更新",         "完全不变（接收 Λⱼ）",          "保留"],
        ["imm_init_models.m","初始化 x,P,F,Q",             "新增 .particles,.weights,.N", "修改"],
        ["track_create.m",   "创建轨迹结构",               "当 cfg.use_pf=true 时初始化粒子","修改"],
        ["track_imm_pda_update.m","调用 ekf_pda_update",  "当 use_pf=true 时调用 pf_jpda_update","修改"],
        ["main_imm_jpda_ekf.m","主循环 KF 流程",           "PF 分支切换+声学场输入",      "修改"],
        ["JPDA 全部",        "β_{t,d} 计算",               "完全不变",                    "保留"],
        ["track_manage.m",   "M/N 管理",                   "完全不变",                    "保留"],
    ]
    t = doc.add_table(rows=1+len(rows), cols=4)
    t.style = "Table Grid"
    col_colors = ["2E75B6","2E75B6","2E75B6","2E75B6"]
    for i, h in enumerate(headers):
        c = t.rows[0].cells[i]
        set_cell_bg(c, "2E75B6")
        pr = c.paragraphs[0]; pr.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = pr.add_run(h); r.bold = True
        r.font.color.rgb = RGBColor(0xFF,0xFF,0xFF); r.font.size = Pt(10)
    status_colors = {"新建":"D9F0D3","保留":"EBF3FB","修改":"FFF2CC","完全不变":"EBF3FB"}
    for ri, row in enumerate(rows):
        for ci, val in enumerate(row):
            c = t.rows[ri+1].cells[ci]
            if ci == 3:
                bg = status_colors.get(val, "FFFFFF")
            else:
                bg = "F9F9F9" if ri % 2 == 0 else "FFFFFF"
            set_cell_bg(c, bg)
            p = c.paragraphs[0]
            r = p.add_run(val)
            r.font.size = Pt(10)
            if ci == 3:
                r.bold = True

# ─────────────────────────── 主文档构建 ────────────────────────────────────

def build_doc():
    doc = Document()

    # 全局字体
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(11)

    # 页边距
    for section in doc.sections:
        section.left_margin   = Cm(2.5)
        section.right_margin  = Cm(2.5)
        section.top_margin    = Cm(2.5)
        section.bottom_margin = Cm(2.0)

    # ── 封面 ────────────────────────────────────────────────────────────────
    doc.add_paragraph()
    doc.add_paragraph()
    t = doc.add_paragraph("IMM-KF → IMM-PF")
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    t.runs[0].font.size = Pt(28); t.runs[0].bold = True
    t.runs[0].font.color.rgb = RGBColor(0x1F, 0x49, 0x7D)

    t2 = doc.add_paragraph("完整数学推导与工程实现指南")
    t2.alignment = WD_ALIGN_PARAGRAPH.CENTER
    t2.runs[0].font.size = Pt(16)
    t2.runs[0].font.color.rgb = RGBColor(0x40, 0x40, 0x40)

    doc.add_paragraph()
    sub = doc.add_paragraph("低空无人机声学-雷达协同多目标跟踪系统\nIMM-JPDA 框架下的粒子滤波升级方案")
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sub.runs[0].font.size = Pt(12)
    sub.runs[0].font.color.rgb = RGBColor(0x55, 0x55, 0x55)

    doc.add_paragraph()
    doc.add_paragraph()
    date_p = doc.add_paragraph("2026 年 6 月")
    date_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    date_p.runs[0].font.size = Pt(11)

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第一章  符号体系
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第一章  符号体系", 1)
    add_para(doc, "下表给出全文统一符号定义。粒子滤波引入粒子集合记号；所有 IMM 层符号（μ、c̄、π、Λ）保持与 KF 版完全一致。", size=11)
    doc.add_paragraph()
    make_sym_table(doc)
    doc.add_paragraph()

    add_para(doc, "状态向量", bold=True, size=11)
    add_equation(doc, "x = [ pₓ,  pᵧ,  vₓ,  vᵧ ]ᵀ  ∈ ℝ⁴", "二维位置 + 二维速度，本系统保持不变")

    add_para(doc, "运动模型（双 DWNA-IMM）", bold=True, size=11)
    add_equation(doc,
        "F = | 1  0  Δt  0  |      Q_j = q_j · | Δt³/3   0      Δt²/2   0    |\n"
        "    | 0  1  0   Δt |                   | 0       Δt³/3  0       Δt²/2|\n"
        "    | 0  0  1   0  |                   | Δt²/2   0      Δt      0    |\n"
        "    | 0  0  0   1  |                   | 0       Δt²/2  0       Δt   |\n"
        "\n"
        "Model 1（低机动）: q₁ = 0.5 m²/s³\n"
        "Model 2（高机动）: q₂ = 8.0 m²/s³")

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第二章  IMM-KF 标准流程回顾
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第二章  IMM-KF 标准流程回顾", 1)
    add_para(doc, "以下为本系统已实现的 IMM-KF 基线流程（Bar-Shalom 2001, §11.6.6）。"
                  "升级为 PF 时，Step 1–2 替换为粒子版本，Step 3 以 JPDA→PF 接口替换 KF 更新，"
                  "Step 4–5 完全保留不变。", size=11)
    doc.add_paragraph()

    steps_kf = [
        ("0", "初始条件",
         ["x̂ⱼ^(k-1|k-1) ∈ ℝ⁴   Pⱼ^(k-1|k-1) ∈ ℝ⁴ˣ⁴   μⱼ^(k-1) ∈ [0,1]",
          "j = 1,…,M   （本系统 M=2）"]),
        ("1", "IMM 混合（Interaction）",
         ["c̄ⱼ    = Σᵢ πᵢⱼ · μᵢ^(k-1)               ← 预测模型概率",
          "μᵢ|ⱼ  = πᵢⱼ · μᵢ^(k-1) / c̄ⱼ             ← 混合系数",
          "x̂₀ⱼ  = Σᵢ μᵢ|ⱼ · x̂ᵢ^(k-1|k-1)           ← 混合均值",
          "P₀ⱼ   = Σᵢ μᵢ|ⱼ · [ Pᵢ + Δxᵢ Δxᵢᵀ ]      ← 混合协方差",
          "        Δxᵢ = x̂ᵢ^(k-1|k-1) − x̂₀ⱼ"]),
        ("2", "KF 预测",
         ["x̂ⱼ⁻ = Fⱼ · x̂₀ⱼ",
          "Pⱼ⁻  = Fⱼ · P₀ⱼ · Fⱼᵀ + Qⱼ"]),
        ("3", "KF-PDA 更新（JPDA 输出 βt,d 输入）",
         ["v̄   = Σd βd · vd,   vd = zd − H x̂ⱼ⁻       ← PDA 合成新息",
          "S̄   = H Pⱼ⁻ Hᵀ + R̄                         ← 新息协方差",
          "K   = Pⱼ⁻ Hᵀ S̄⁻¹                            ← 卡尔曼增益",
          "x̂ⱼ⁺ = x̂ⱼ⁻ + K v̄                           ← 状态更新",
          "Pⱼ⁺ = β₀ Pⱼ⁻ + (1−β₀)(Pⱼ⁻−K S̄ Kᵀ) + K [Σd βd vd vdᵀ − v̄ v̄ᵀ] Kᵀ"]),
        ("4", "IMM 模型概率更新（imm_update_mu，保留不变）",
         ["μⱼ^(k) = Λⱼ · c̄ⱼ / Σₗ ( Λₗ · c̄ₗ )",
          "Λⱼ: KF 中为归一化常数（高斯密度值）"]),
        ("5", "IMM 融合输出（imm_fuse，保留不变）",
         ["x̂^(k) = Σⱼ μⱼ^(k) · x̂ⱼ⁺",
          "P^(k)  = Σⱼ μⱼ^(k) · [ Pⱼ⁺ + (x̂ⱼ⁺−x̂^(k))(x̂ⱼ⁺−x̂^(k))ᵀ ]"]),
    ]
    for sno, stitle, scontent in steps_kf:
        add_step_box(doc, sno, stitle, scontent)
        doc.add_paragraph().paragraph_format.space_after = Pt(4)

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第三章  IMM-PF 完整推导
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第三章  IMM-PF 完整推导", 1)

    # 3.1 核心思路
    add_heading(doc, "3.1  核心思路", 2)
    add_para(doc,
        "将每个 IMM 子模型的【高斯假设】（KF）替换为【非参数粒子集】（PF）。"
        "IMM 框架的三个核心操作完全保留：", size=11)
    items = [
        "模型概率更新（imm_update_mu）：接口不变，仅将 Λⱼ 替换为粒子似然均值",
        "融合输出（imm_fuse）：接口不变，由粒子集提取均值/协方差后传入",
        "Markov 转移矩阵 π：不变",
    ]
    for it in items:
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(it).font.size = Pt(11)
        p.paragraph_format.left_indent = Inches(0.3)

    # 3.2 粒子表示
    add_heading(doc, "3.2  粒子表示", 2)
    add_para(doc,
        "模型 j 在帧 k 的状态用 N 个带权粒子表示：", size=11)
    add_equation(doc,
        "𝒳ⱼ^(k) = { ( xⱼ^(k,i),  wⱼ^(k,i) ) }ᵢ₌₁ᴺ\n"
        "\n"
        "约束:  Σᵢ wⱼ^(k,i) = 1,   wⱼ^(k,i) ≥ 0")

    add_para(doc, "粒子均值（与 KF model.x 对齐）：", size=11)
    add_equation(doc, "x̂ⱼ = Σᵢ wⱼ^(i) · xⱼ^(i)")
    add_para(doc, "粒子协方差（与 KF model.P 对齐，供 imm_fuse 使用）：", size=11)
    add_equation(doc,
        "Pⱼ = Σᵢ wⱼ^(i) · ( xⱼ^(i) − x̂ⱼ ) ( xⱼ^(i) − x̂ⱼ )ᵀ")

    # 3.3 Step 1: IMM 粒子混合
    add_heading(doc, "3.3  Step 1：IMM 粒子混合（imm_pf_mix）", 2)
    add_para(doc,
        "目标：为模型 j 构造混合后粒子集，使其边际分布近似为：", size=11)
    add_equation(doc,
        "p₀ⱼ(x) = Σᵢ μᵢ|ⱼ · pᵢ^(k-1|k-1)(x)\n"
        "\n"
        "μᵢ|ⱼ = πᵢⱼ · μᵢ^(k-1) / c̄ⱼ,    c̄ⱼ = Σᵢ πᵢⱼ · μᵢ^(k-1)")

    add_para(doc, "算法实现（对每个粒子 m = 1,…,N 独立执行）：", bold=True, size=11)
    algo_lines = [
        "① 以概率 { μ₁|ⱼ, μ₂|ⱼ, …, μₘ|ⱼ } 选择源模型 i*",
        "② 从模型 i* 的粒子集，按归一化权重 { wᵢ*^(l) } 重采样一个粒子 x^(lₘ)",
        "③ 令 x₀ⱼ^(m) = xᵢ*^(lₘ)",
        "④ 令 w₀ⱼ^(m) = 1/N   （混合后均匀初始化）",
    ]
    for line in algo_lines:
        p = doc.add_paragraph(line)
        p.paragraph_format.left_indent = Inches(0.4)
        p.runs[0].font.size = Pt(10.5)
        p.paragraph_format.space_after = Pt(2)

    add_para(doc,
        "精度说明：N 足够大时，此操作等效于 KF 混合步的均值/协方差计算，"
        "且天然处理非高斯分布，无需线性化。", size=10,
        color="444444", indent=0.2)

    # 3.4 Step 2: 粒子预测
    add_heading(doc, "3.4  Step 2：粒子预测（imm_pf_predict）", 2)
    add_para(doc, "对模型 j 的每个混合粒子施加动力学传播：", size=11)
    add_equation(doc,
        "xⱼ^(k,−,m) = Fⱼ · x₀ⱼ^(m) + wⱼ^(m)\n"
        "\n"
        "wⱼ^(m) ~ 𝒩( 0, Qⱼ )   独立同分布过程噪声\n"
        "\n"
        "预测后权重保持均匀:  wⱼ^(k,−,m) = 1/N")
    add_para(doc,
        "向量化实现：一次 Cholesky 分解 Qⱼ = LLᵀ，"
        "生成全部 N 列噪声  noise = L · randn(4, N)，"
        "避免逐粒子循环，计算效率提升约 N 倍。",
        size=10, color="444444", indent=0.2)

    # 3.5 JPDA → PF 接口
    add_heading(doc, "3.5  Step 3：JPDA → PF 权重更新（pf_jpda_update）", 2)
    add_heading(doc, "3.5.1  JPDA 输出的使用", 3)
    add_para(doc,
        "JPDA 模块（jpda_run.m）完全不修改，输出边缘关联概率矩阵：", size=11)
    add_equation(doc,
        "β_{t,d} = P( 量测 d 来自轨迹 t )\n"
        "β_{t,0} = P( 轨迹 t 漏检 )\n"
        "\n"
        "满足:  β_{t,0} + Σd β_{t,d} = 1")
    add_para(doc,
        "β 由 JPDA 在轨迹层计算，与 IMM 子模型无关。"
        "每个子模型 PF 使用相同的 β（轨迹级关联概率），"
        "不区分模型索引 j。", size=11)

    add_heading(doc, "3.5.2  雷达量测似然（每粒子）", 3)
    add_para(doc, "对模型 j 的粒子 m，雷达 PDA 加权似然：", size=11)
    add_equation(doc,
        "ℓⱼ^radar(m) = β₀ · λc + Σd β_{t,d} · 𝒩( zd ;  H xⱼ^(m), Rd )\n"
        "\n"
        "𝒩( z ; μ, Σ ) = (2π)^(−m/2) |Σ|^(−1/2) exp( −½ (z−μ)ᵀ Σ⁻¹ (z−μ) )\n"
        "\n"
        "λc = 杂波空间密度（来自 cfg.lambda_c，不变）")
    add_para(doc,
        "推导：将 JPDA 的 β 视为数据关联事件后验概率，"
        "粒子 m 的联合似然 = Σθ P(θ|Z)·p(Z|θ, x^(m))，"
        "在单轨迹关联假设下近似为上式（与 KF-PDA 推导路径一致）。",
        size=10, color="444444", indent=0.2)

    add_heading(doc, "3.5.3  声学概率场似然（每粒子）", 3)
    add_para(doc,
        "声学模块输出二维概率场（build_acoustic_prob_field 的 pack.post），"
        "定义在网格 (xg, zg) 上的概率质量函数：", size=11)
    add_equation(doc,
        "𝐩_ac = { p_ac(u,v) }_{u,v},     Σ_{u,v} p_ac(u,v) = 1\n"
        "\n"
        "网格来自 acoustic_config:  xc (1×nₓ mm),  zc (1×nz mm)")
    add_para(doc, "粒子 m 的声学似然（双线性插值）：", size=11)
    add_equation(doc,
        "ℓⱼ^ac(m) = interp2( 𝐩_ac, pₓ^(m) · κ,  pᵧ^(m) · κ )\n"
        "\n"
        "κ = cfg.pf_m2ac_scale = 1000   （m → mm 单位转换）\n"
        "\n"
        "越界粒子:  ℓⱼ^ac(m) ← cfg.pf_ac_bg_likelihood（背景似然，≈ 0 但非零）")

    add_heading(doc, "3.5.4  融合似然与声学权重 γ", 3)
    add_para(doc, "动态融合权重由声学探测置信度决定：", size=11)
    add_equation(doc,
        "γ = pack.p_detect_joint    （pf_ac_gamma_mode = 'snr'，推荐）\n"
        "γ = cfg.pf_ac_gamma_fixed  （固定权重，对照实验）\n"
        "γ ∈ [0, 1]")
    add_para(doc, "融合似然（乘法融合）：", size=11)
    add_equation(doc,
        "ℓⱼ(m) = ℓⱼ^radar(m) · [ γ · ℓⱼ^ac(m)  +  (1−γ) ]\n"
        "\n"
        "当 γ → 0（声学未检测）: ℓⱼ(m) ≈ ℓⱼ^radar(m)  退化为纯雷达\n"
        "当 γ → 1（声学强信号）: ℓⱼ(m) = ℓⱼ^radar(m) · ℓⱼ^ac(m)  完全融合",
        "各传感器相互独立假设下的贝叶斯乘法融合")

    add_heading(doc, "3.5.5  权重更新", 3)
    add_equation(doc,
        "w̃ⱼ^(k,m) = wⱼ^(k,−,m) · ℓⱼ(m)        ← 未归一化更新\n"
        "\n"
        "wⱼ^(k,m) = w̃ⱼ^(k,m) / Σₗ w̃ⱼ^(k,l)  ← 归一化")

    # 3.6 模型似然
    add_heading(doc, "3.6  Step 4：PF 模型似然 Λⱼ（供 imm_update_mu）", 2)
    add_para(doc,
        "KF 中 Λⱼ 是高斯归一化常数。PF 中用粒子似然均值替代：", size=11)
    add_equation(doc,
        "Λⱼ = mean{ ℓⱼ(m) }ₘ₌₁ᴺ\n"
        "\n"
        "等效：Λⱼ = N · mean{ w̃ⱼ^(k,m) }   （w̃ 是未归一化权重，预测权重 = 1/N）\n"
        "\n"
        "实用下限保护:  Λⱼ ← max( Λⱼ, cfg.pf_lambda_floor )",
        "仅需相对大小，与 KF Λ 量纲无关")
    add_para(doc,
        "将此 Λⱼ 代入 imm_update_mu（接口完全不变）：", size=11)
    add_equation(doc,
        "μⱼ^(k) = Λⱼ · c̄ⱼ / Σₗ ( Λₗ · c̄ₗ )     ← imm_update_mu.m 原有公式，不改动")

    # 3.7 ESS 与重采样
    add_heading(doc, "3.7  Step 5：ESS 计算与系统重采样", 2)
    add_para(doc, "有效样本量（Effective Sample Size）：", size=11)
    add_equation(doc,
        "ESSⱼ = 1 / Σᵢ ( wⱼ^(i) )²\n"
        "\n"
        "ESSⱼ ∈ [1, N]：\n"
        "  → ESSⱼ ≈ N  权重均匀（理想状态）\n"
        "  → ESSⱼ ≈ 1  权重退化（单粒子支配）")
    add_para(doc, "重采样触发条件：", size=11)
    add_equation(doc,
        "触发:  ESSⱼ < α · N,    α = cfg.pf_resample_thresh = 0.5（推荐）")
    add_para(doc, "系统重采样算法（O(N)，方差最小）：", bold=True, size=11)
    add_equation(doc,
        "① 计算累积权重:  Cₘ = Σₗ₌₁ᵐ w^(l),   C_N = 1\n"
        "② 采样起点:      U₁ ~ Uniform[0, 1/N)\n"
        "③ 等间距格点:    Uₘ = U₁ + (m−1)/N,   m = 1,…,N\n"
        "④ 查找索引:      i*(m) = min{ j : Cⱼ ≥ Uₘ }\n"
        "⑤ 更新:          x^(m) ← x^(i*(m)),   w^(m) ← 1/N",
        "等价于分层重采样，比 multinomial 方差小约 N 倍")

    # 3.8 融合输出
    add_heading(doc, "3.8  Step 6：IMM 融合输出（imm_fuse，保留不变）", 2)
    add_para(doc,
        "由粒子集提取均值与协方差，写入 model.x 和 model.P，"
        "再传入原有的 imm_fuse（完全不修改）：", size=11)
    add_equation(doc,
        "x̂ⱼ = Σᵢ wⱼ^(i) · xⱼ^(i)                   → model(j).x\n"
        "Pⱼ  = Σᵢ wⱼ^(i) (xⱼ^(i)−x̂ⱼ)(xⱼ^(i)−x̂ⱼ)ᵀ  → model(j).P\n"
        "\n"
        "imm_fuse( models, μ ) 原有公式（不变）：\n"
        "  x̂   = Σⱼ μⱼ · x̂ⱼ\n"
        "  P    = Σⱼ μⱼ · [ Pⱼ + (x̂ⱼ−x̂)(x̂ⱼ−x̂)ᵀ ]")

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第四章  声学概率场推导
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第四章  声学概率场建模", 1)

    add_heading(doc, "4.1  从点量测到概率场", 2)
    add_para(doc,
        "原始系统将声学建模为高斯点量测 z_ac ~ N(Hx, R_ac)，"
        "存在以下问题：\n"
        "① 声学 DOA 估计误差在高 SNR 时才近似高斯，低 SNR 时严重非高斯；\n"
        "② 单点量测损失了空间分布信息；\n"
        "③ 与 PF 框架的粒子权重更新不自然对齐。\n\n"
        "升级方案：声学输出二维概率场 P_ac(x,y)，"
        "由麦克风节点阵列的声压级差特征构建。", size=11)

    add_heading(doc, "4.2  单节点 SPL 物理模型", 2)
    add_equation(doc,
        "Lp_m(r) = A_m − n_m · log₁₀(r) − α · r    [dB]\n"
        "\n"
        "r       = 声源到节点 m 的距离 (m)\n"
        "A_m     = 标定参考声压级 (dB)\n"
        "n_m     = 球面扩散斜率（理想自由场 ≈ 20，含地面反射 8~9）\n"
        "α       = 大气吸收系数 (dB/m)，可选",
        "来自 acoustic_config.m 标定参数")
    add_para(doc, "地面反射（镜像声源）：", size=11)
    add_equation(doc,
        "Lp_m^total = 10 log₁₀ ( 10^(Lp_direct/10)  +  R² · 10^(Lp_image/10) )\n"
        "\n"
        "R = cfg.R_ground ∈ [0,1]，反射系数；R=0 退化为自由场")

    add_heading(doc, "4.3  多节点协同得分场", 2)
    add_equation(doc,
        "对候选源位置 (u,v) ∈ 网格:\n"
        "\n"
        "SNR_m(u,v) = Lp_m(u,v) − NF_m       （各节点 SNR）\n"
        "\n"
        "w_m = softmax( SNR_m / τ_w )         （可靠性权重，τ_w = cfg.weight_tau）\n"
        "\n"
        "S(u,v) = Σₘ w_m · ΔLp_m(u,v)        （加权对数似然得分场）\n"
        "\n"
        "ΔLp_m(u,v) = −½ ( Lp_m^meas − Lp_m(u,v) )² / σ_m²",
        "build_score_grid.m")

    add_heading(doc, "4.4  概率场归一化（三种方案）", 2)
    schemes = [
        ("方案一  Softmax（玻尔兹曼分布）",
         "P_softmax(u,v) = exp( S(u,v) / τ ) / Σ_{u',v'} exp( S(u',v') / τ )\n"
         "τ = τ₀ · exp( c · max(0, SNR_ref − SNR̄) )  （低 SNR 时温度升高，分布变平）"),
        ("方案二  能量域强度归一化",
         "P_intensity(u,v) ∝ max( 10^(S/10) − I_noise, 0 )^α\n"
         "归一化使 Σ P = 1"),
        ("方案三  NMS-GMM 峰值混合",
         "① NMS 提取 K 个局部极大值峰  { (uₖ, vₖ, sₖ) }\n"
         "② 各峰权重:  ωₖ = softmax(sₖ / τ)\n"
         "③ GMM:  P_nms(u,v) = Σₖ ωₖ · 𝒩( (u,v); (uₖ,vₖ), σ_k²I )"),
    ]
    for name, eq in schemes:
        add_para(doc, name, bold=True, size=11)
        add_equation(doc, eq)

    add_heading(doc, "4.5  递归贝叶斯更新（帧间平滑）", 2)
    add_equation(doc,
        "预测（Chapman-Kolmogorov 扩散）:\n"
        "  p̄^(k)(u,v) = ∫ p^(k-1)(u',v') · 𝒩( (u,v); (u',v'), σ_diff² I ) du'dv'\n"
        "  实现：高斯卷积，σ_diff = cfg.bayes_diffusion_mm\n"
        "\n"
        "更新（贝叶斯）:\n"
        "  p^(k)(u,v) ∝ p̄^(k)(u,v) · L_sel(u,v)   （选定似然场）\n"
        "\n"
        "熵门控（低信息量时跳过更新，防止错误累积）:\n"
        "  H_norm = −Σ L log L / log(N_grid)   > cfg.bayes_entropy_gate → 跳过",
        "recursive_bayes_field.m")

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第五章  文件修改清单与接口对齐表
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第五章  文件修改清单与接口对齐", 1)

    add_heading(doc, "5.1  KF → PF 接口映射表", 2)
    make_kf_to_pf_table(doc)
    doc.add_paragraph()

    add_heading(doc, "5.2  新建文件清单", 2)
    new_files = [
        ("PF/pf_init.m",               "从高斯初始分布采样 N 粒子，初始化单模型 PF struct"),
        ("PF/pf_extract_state.m",       "从粒子集计算均值/协方差，同步 pf.x, pf.P"),
        ("PF/pf_systematic_resample.m", "系统重采样，O(N)，权重均匀化"),
        ("PF/pf_jpda_update.m",         "粒子权重更新（雷达 PDA + 声学场插值融合）"),
        ("PF/pf_compute_lambda.m",      "计算 Λⱼ = mean(ℓⱼ) 供 imm_update_mu 使用"),
        ("IMM/imm_pf_mix.m",            "粒子重混合（替代 imm_mix 的高斯参数版本）"),
        ("IMM/imm_pf_predict.m",        "粒子传播（替代 imm_predict 的 KF 版本）"),
        ("Scenario/acoustic_field_to_tracker.m", "将 simulate_node_spl + build_acoustic_prob_field 封装为跟踪器接口"),
        ("Evaluation/plot_results_pf.m","8 幅完整可视化（轨迹/模型概率/粒子云/声学场/ESS/JPDA/RMSE/OSPA）"),
    ]
    t = doc.add_table(rows=1+len(new_files), cols=2)
    t.style = "Table Grid"
    for i, h in enumerate(["文件路径", "功能说明"]):
        c = t.rows[0].cells[i]; set_cell_bg(c, "2E75B6")
        r = c.paragraphs[0].add_run(h); r.bold = True
        r.font.color.rgb = RGBColor(0xFF,0xFF,0xFF); r.font.size = Pt(10.5)
    for ri, (fp, desc) in enumerate(new_files):
        bg = "D9F0D3" if ri % 2 == 0 else "F0FAF0"
        c0 = t.rows[ri+1].cells[0]; c1 = t.rows[ri+1].cells[1]
        set_cell_bg(c0, bg); set_cell_bg(c1, bg)
        r0 = c0.paragraphs[0].add_run(fp); r0.font.name = "Courier New"; r0.font.size = Pt(9.5)
        r1 = c1.paragraphs[0].add_run(desc); r1.font.size = Pt(10)

    add_heading(doc, "5.3  修改文件清单", 2)
    mod_files = [
        ("config_tracker.m",            "追加 PF 参数：use_pf, pf_N, pf_resample_thresh 等"),
        ("IMM/imm_init_models.m",       "当 cfg.use_pf=true 时，为每个 model 追加 .particles/.weights/.N 字段"),
        ("TrackManagement/track_create.m","调用 imm_init_models 后，若 use_pf 则用 pf_init 初始化粒子"),
        ("TrackManagement/track_imm_pda_update.m","若 use_pf 则调用 pf_jpda_update 替代 ekf_pda_update"),
        ("main_imm_jpda_ekf.m",         "① 调用 acoustic_field_to_tracker 生成 ac_pack\n② 混合/预测步调用 PF 版本\n③ 传 ac_pack 至 pf_jpda_update"),
    ]
    t2 = doc.add_table(rows=1+len(mod_files), cols=2)
    t2.style = "Table Grid"
    for i, h in enumerate(["文件路径", "修改内容"]):
        c = t2.rows[0].cells[i]; set_cell_bg(c, "C55A11")
        r = c.paragraphs[0].add_run(h); r.bold = True
        r.font.color.rgb = RGBColor(0xFF,0xFF,0xFF); r.font.size = Pt(10.5)
    for ri, (fp, desc) in enumerate(mod_files):
        bg = "FFF2CC" if ri % 2 == 0 else "FFFBEE"
        c0 = t2.rows[ri+1].cells[0]; c1 = t2.rows[ri+1].cells[1]
        set_cell_bg(c0, bg); set_cell_bg(c1, bg)
        r0 = c0.paragraphs[0].add_run(fp); r0.font.name = "Courier New"; r0.font.size = Pt(9.5)
        r1 = c1.paragraphs[0].add_run(desc); r1.font.size = Pt(10)

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第六章  IMM-PF 完整流程图
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第六章  IMM-PF 完整帧处理流程", 1)

    flow_text = (
        "帧 k-1 后验:  { 粒子集ⱼ, μⱼ }  (j=1..M)\n"
        "         │\n"
        "         ▼\n"
        "┌─────────────────────────────────────────────────────────────┐\n"
        "│  Step 1: IMM 粒子混合  [imm_pf_mix.m]                       │\n"
        "│  c̄ⱼ = Σ πᵢⱼ μᵢ                                            │\n"
        "│  μᵢ|ⱼ = πᵢⱼμᵢ / c̄ⱼ                                       │\n"
        "│  x₀ⱼ^(m) ← 按 μᵢ|ⱼ 选源模型 i*，按 wᵢ* 重采样粒子           │\n"
        "│  w₀ⱼ^(m) = 1/N                                              │\n"
        "└─────────────────────────────────────────────────────────────┘\n"
        "         │\n"
        "         ▼\n"
        "┌─────────────────────────────────────────────────────────────┐\n"
        "│  Step 2: 粒子预测  [imm_pf_predict.m]                        │\n"
        "│  xⱼ^(k,−,m) = Fⱼ x₀ⱼ^(m) + 𝒩(0, Qⱼ)                      │\n"
        "│  w = 1/N  (预测后均匀)                                        │\n"
        "└─────────────────────────────────────────────────────────────┘\n"
        "         │\n"
        "         ▼  同时\n"
        "┌────────────────────────────────┐    ┌───────────────────────┐\n"
        "│  Step 2b: 融合预测（供门控）    │    │  声学场生成            │\n"
        "│  x̂⁻ = Σⱼ c̄ⱼ x̂ⱼ^(−)         │    │  simulate_node_spl    │\n"
        "│  P⁻  = Σⱼ c̄ⱼ [Pⱼ+Δxⱼ Δxⱼᵀ] │    │  → build_acoustic_    │\n"
        "│  ↓ imm_fuse(models, c_bar)     │    │    prob_field         │\n"
        "│  送 JPDA 门控（track.x, P）     │    │  → ac_pack.post       │\n"
        "└────────────────────────────────┘    └───────────────────────┘\n"
        "         │                                        │\n"
        "         ▼                                        │\n"
        "┌─────────────────────────────────────────────────┘\n"
        "│  Step 3: JPDA（完全不变，jpda_run.m）\n"
        "│  输出 β_{t,d}, β_{t,0}\n"
        "└─────────────────────────────────────────────────┐\n"
        "         │                                        │\n"
        "         ▼                                        ▼\n"
        "┌─────────────────────────────────────────────────────────────┐\n"
        "│  Step 4: PF 权重更新  [pf_jpda_update.m]                     │\n"
        "│  ℓ^radar(m) = β₀λc + Σd βd · 𝒩(zd; Hx^(m), Rd)            │\n"
        "│  ℓ^ac(m)    = interp2( pack.post, pₓ·κ, pᵧ·κ )             │\n"
        "│  ℓ(m)       = ℓ^radar(m) · [γ·ℓ^ac(m) + (1−γ)]            │\n"
        "│  w̃(m) = w(m)·ℓ(m),  归一化 → w(m)                          │\n"
        "│  Λⱼ = mean(ℓ(m))                                            │\n"
        "└─────────────────────────────────────────────────────────────┘\n"
        "         │\n"
        "         ▼\n"
        "┌─────────────────────────────────────────────────────────────┐\n"
        "│  Step 5: ESS 检验 + 系统重采样  [pf_systematic_resample.m]   │\n"
        "│  ESS = 1/Σ(w²),  若 ESS < 0.5N → 重采样                    │\n"
        "└─────────────────────────────────────────────────────────────┘\n"
        "         │\n"
        "         ▼\n"
        "┌─────────────────────────────────────────────────────────────┐\n"
        "│  Step 6: IMM 模型概率更新  [imm_update_mu.m，不变]           │\n"
        "│  μⱼ^(k) = Λⱼ c̄ⱼ / Σₗ Λₗ c̄ₗ                             │\n"
        "└─────────────────────────────────────────────────────────────┘\n"
        "         │\n"
        "         ▼\n"
        "┌─────────────────────────────────────────────────────────────┐\n"
        "│  Step 7: 粒子均值/协方差提取 → imm_fuse  [不变]              │\n"
        "│  x̂ = Σⱼ μⱼ x̂ⱼ,   P = Σⱼ μⱼ[Pⱼ + (x̂ⱼ−x̂)(x̂ⱼ−x̂)ᵀ]     │\n"
        "└─────────────────────────────────────────────────────────────┘\n"
        "         │\n"
        "         ▼\n"
        "  轨迹管理 track_manage.m（不变）→ 帧 k 后验"
    )
    p = doc.add_paragraph()
    r = p.add_run(flow_text)
    r.font.name = "Courier New"; r.font.size = Pt(8.5)
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after  = Pt(8)

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第七章  可视化方案
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第七章  可视化方案（8 幅图）", 1)
    fig_items = [
        ("图 1", "真实轨迹 vs 估计轨迹",
         "plot_results_pf.m → subplot(2,4,1)\n"
         "• 真值：colored solid line（每目标一色）\n"
         "• PF 估计：dashed line（同色）\n"
         "• 诞生/消亡帧标记"),
        ("图 2", "IMM 模型概率曲线",
         "subplot(2,4,2)\n"
         "• y轴：μ₁(k)（蓝）, μ₂(k)（红）\n"
         "• 标注机动区间（μ₂ > 0.5 时填色）"),
        ("图 3", "粒子云演化（选定帧）",
         "subplot(2,4,3)\n"
         "• scatter：粒子位置，颜色编码权重\n"
         "• 叠加真值与均值估计\n"
         "• 选 k = 50, 100, 150 三帧快照"),
        ("图 4", "声学概率场热力图",
         "subplot(2,4,4)\n"
         "• imagesc(xc, zc, pack.post')\n"
         "• colorbar，叠加粒子位置\n"
         "• 麦克风节点位置标注"),
        ("图 5", "ESS 曲线",
         "subplot(2,4,5)\n"
         "• ESSⱼ(k)（每模型一条线）\n"
         "• 红色虚线：阈值 0.5N\n"
         "• 重采样帧以竖线标记"),
        ("图 6", "JPDA 关联概率矩阵",
         "subplot(2,4,6)\n"
         "• imagesc(β)\n"
         "• 行=轨迹，列=量测（含 β₀）\n"
         "• 选取量测数最多的帧显示"),
        ("图 7", "RMSE 曲线",
         "subplot(2,4,7)\n"
         "• RMSE_pos(k) = mean_{j} ||p̂ⱼ−pⱼ||\n"
         "• RMSE_vel(k) 同步\n"
         "• 对比 KF 基线（虚线）"),
        ("图 8", "OSPA 指标",
         "subplot(2,4,8)\n"
         "• OSPA(k)（p=2, c=100 m）\n"
         "• 分解：距离项 + 势项（面积填色）"),
    ]
    for fn, ftitle, fdesc in fig_items:
        add_para(doc, f"{fn}：{ftitle}", bold=True, size=11, spacing_after=2)
        add_para(doc, fdesc, size=10, indent=0.3, color="333333", spacing_after=6)

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第八章  麦克风布设接口说明
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第八章  麦克风节点布设接口对齐", 1)

    add_heading(doc, "8.1  当前声学系统（matlab/ 文件夹）", 2)
    add_para(doc,
        "matlab/ 文件夹已实现完整的声学概率场计算链路（已验证），"
        "只需对齐以下两处接口：\n"
        "① 麦克风坐标单位：acoustic_config 使用 mm，跟踪器使用 m → 需在接口层统一；\n"
        "② 概率场网格范围：需覆盖跟踪场景的空间范围（与 run_scenario 场景对齐）。",
        size=11)

    add_heading(doc, "8.2  acoustic_config 修改说明", 2)
    add_para(doc,
        "仅修改以下两处，其余参数不变：", size=11)
    add_equation(doc,
        "【修改 1】麦克风坐标（与场景空间对齐，单位 mm）\n"
        "cfg.mic_coords = [   % 根据 run_scenario 实际 acousticPos(m) × 1000\n"
        "    x₁·1000,  y₁·1000,  z₁·1000;  % mic 1\n"
        "    x₂·1000,  y₂·1000,  z₂·1000;  % mic 2\n"
        "    ...                             % 按 acousticPos 行数扩充\n"
        "];\n"
        "\n"
        "【修改 2】概率场网格范围与高度（mm）\n"
        "cfg.src_height_mm = mean_altitude × 1000;  % 场景飞行高度均值\n"
        "cfg.placement_area_mm = [xmin, xmax, zmin, zmax] × 1000;  % 场景范围")

    add_heading(doc, "8.3  接口封装函数说明", 2)
    add_para(doc,
        "acoustic_field_to_tracker.m 封装全部调用，主循环只需一行：", size=11)
    add_equation(doc,
        "% 主循环中（每帧调用）：\n"
        "ac_pack = acoustic_field_to_tracker( truth_positions_k, ac_cfg, prev_post );\n"
        "\n"
        "% 函数内部流程：\n"
        "% 1. src_xyz (m) → src_xyz_mm = src_xyz × 1000\n"
        "% 2. measured_lp = simulate_node_spl(src_xyz_mm, ac_cfg)  （仿真帧）\n"
        "% 3. pack = build_acoustic_prob_field(measured_lp, ac_cfg, 'softmax', prev_post)\n"
        "% 4. 返回 pack（含 post, detected, p_detect_joint 等）")

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第九章  稳定性与数值保护
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第九章  数值稳定性设计", 1)

    stab_items = [
        ("权重退化保护",
         "若 Σ w̃ < 1e-300（所有粒子似然接近零）\n→ 均匀分布 w = 1/N，记录警告"),
        ("模型似然下限",
         "Λⱼ = max( mean(ℓ), cfg.pf_lambda_floor = 1e-300 )\n→ 防止 imm_update_mu 分母为零"),
        ("粒子越界",
         "声学场插值越界粒子赋背景似然 pf_ac_bg_likelihood = 1e-6\n→ 不完全消除但权重极低，保留先验信息"),
        ("协方差对称性",
         "pf_extract_state: P = (P + Pᵀ)/2\nimm_pf_mix: P_mix = (P + Pᵀ)/2\n→ 防止浮点误差累积导致非对称"),
        ("粒子多样性",
         "重采样后添加轻微扰动（可选）：x^(m) += ε · randn，ε 较小\n→ 防止粒子贫化，尤其在长时间无量测帧"),
        ("Markov 概率保护",
         "c̄ⱼ = max(c̄ⱼ, 1e-300)（来自 imm_mix，已有实现）\n→ 不修改此保护逻辑"),
    ]
    t3 = doc.add_table(rows=1+len(stab_items), cols=2)
    t3.style = "Table Grid"
    for i, h in enumerate(["问题", "处理方案"]):
        c = t3.rows[0].cells[i]; set_cell_bg(c, "375623")
        r = c.paragraphs[0].add_run(h); r.bold = True
        r.font.color.rgb = RGBColor(0xFF,0xFF,0xFF); r.font.size = Pt(10.5)
    for ri, (prob, sol) in enumerate(stab_items):
        bg = "E2EFDA" if ri % 2 == 0 else "F4F9F1"
        c0 = t3.rows[ri+1].cells[0]; c1 = t3.rows[ri+1].cells[1]
        set_cell_bg(c0, bg); set_cell_bg(c1, bg)
        r0 = c0.paragraphs[0].add_run(prob); r0.bold = True; r0.font.size = Pt(10)
        r1 = c1.paragraphs[0].add_run(sol); r1.font.name = "Courier New"; r1.font.size = Pt(9)

    doc.add_page_break()

    # ════════════════════════════════════════════════════════════════
    # 第十章  测试方案
    # ════════════════════════════════════════════════════════════════
    add_heading(doc, "第十章  完整测试方案", 1)

    add_heading(doc, "10.1  单元测试", 2)
    unit_tests = [
        "test_pf_init：验证粒子均值接近 x0，协方差接近 P0×scale",
        "test_pf_mix：验证混合后粒子均值 = 手算均值，误差 < 5%",
        "test_pf_predict：无量测帧，验证均值传播 = F·x，协方差单调增大",
        "test_pf_weight：给定已知位置粒子 + 精确量测，验证权重集中于真实粒子",
        "test_pf_resample：验证 ESS 触发后粒子均匀分布，权重 = 1/N",
        "test_acoustic_field：验证 build_acoustic_prob_field 输出 sum(post)≈1，峰值在真实位置附近",
        "test_imm_update_mu：单模型主导时，对应 Λ 高，该模型概率应提升",
    ]
    for it in unit_tests:
        p = doc.add_paragraph(it, style="List Bullet")
        p.runs[0].font.size = Pt(10.5)
        p.paragraph_format.left_indent = Inches(0.3)

    add_heading(doc, "10.2  集成测试", 2)
    add_para(doc, "运行 main_imm_jpda_ekf.m（cfg.use_pf = true），检验：", size=11)
    integ_tests = [
        "无报错完成全部 N 帧循环",
        "轨迹数量与 KF 基线相当（≥3 条确认轨迹）",
        "RMSE < 2×KF_RMSE（PF 初期略差于 KF 属正常，N=500 时应接近）",
        "ESS 曲线不全为 1（说明重采样正常触发）",
        "模型概率曲线有切换（说明 Λⱼ 正常反映量测支持度）",
        "plot_results_pf 生成 8 幅图无报错",
    ]
    for it in integ_tests:
        p = doc.add_paragraph(it, style="List Bullet")
        p.runs[0].font.size = Pt(10.5)
        p.paragraph_format.left_indent = Inches(0.3)

    add_heading(doc, "10.3  消融实验", 2)
    ablation = [
        ("cfg.use_pf = false", "KF 基线，对比 RMSE/OSPA 参考值"),
        ("cfg.pf_N = 100 / 500 / 2000", "粒子数敏感性，观察 RMSE 与计算时间权衡"),
        ("cfg.pf_ac_gamma_mode = 'fixed', gamma = 0", "纯雷达 PF，验证声学场贡献"),
        ("cfg.pf_ac_gamma_mode = 'fixed', gamma = 1", "纯声学 PF，验证声学场单独性能"),
        ("cfg.pf_resample_thresh = 0.3 / 0.5 / 0.8", "重采样频率对稳定性的影响"),
    ]
    t4 = doc.add_table(rows=1+len(ablation), cols=2)
    t4.style = "Table Grid"
    for i, h in enumerate(["实验配置", "观测指标"]):
        c = t4.rows[0].cells[i]; set_cell_bg(c, "7030A0")
        r = c.paragraphs[0].add_run(h); r.bold = True
        r.font.color.rgb = RGBColor(0xFF,0xFF,0xFF); r.font.size = Pt(10.5)
    for ri, (cfg_str, obs) in enumerate(ablation):
        bg = "F4E6FF" if ri % 2 == 0 else "FAF4FF"
        c0 = t4.rows[ri+1].cells[0]; c1 = t4.rows[ri+1].cells[1]
        set_cell_bg(c0, bg); set_cell_bg(c1, bg)
        r0 = c0.paragraphs[0].add_run(cfg_str); r0.font.name = "Courier New"; r0.font.size = Pt(9.5)
        r1 = c1.paragraphs[0].add_run(obs); r1.font.size = Pt(10)

    # 保存
    doc.save(OUT_PATH)
    print(f"[OK] 文档已保存: {OUT_PATH}")

if __name__ == "__main__":
    build_doc()
