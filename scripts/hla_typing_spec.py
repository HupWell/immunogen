# -*- coding: utf-8 -*-
"""
hla_typing.json 键名与解析辅助（契约见 docs/hla_typing.md）。

供 validate_input、后续 NetMHCIIpan 等模块复用，避免键名写散在多处。
"""

# I 类：BioDriver 必须提供
REQUIRED_HLA_I_KEYS = ("HLA-A", "HLA-B", "HLA-C")

# II 类：可选；用于 NetMHCIIpan 与 MHC-II 表位，OptiType 不产出 II，常由 HLA-HD 等补充
OPTIONAL_HLA_II_KEYS = (
    "HLA-DRB1",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DPA1",
    "HLA-DPB1",
)

# 已知且参与校验的键；其余键名打印警告后忽略，便于未来扩展
KNOWN_HLA_KEYS = REQUIRED_HLA_I_KEYS + OPTIONAL_HLA_II_KEYS


def flatten_hla_class_ii(hla_json: dict) -> list:
    """
    从契约 JSON 中收集 MHC-II 等位基因列表（去空白，失败键跳过）。
    顺序：DRB1 → DQA1 → DQB1 → DPA1 → DPB1。
    """
    out = []
    for k in OPTIONAL_HLA_II_KEYS:
        vals = hla_json.get(k)
        if not isinstance(vals, list):
            continue
        for x in vals:
            s = str(x).strip()
            if s:
                out.append(s)
    return out
