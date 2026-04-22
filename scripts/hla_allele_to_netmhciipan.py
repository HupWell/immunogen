# -*- coding: utf-8 -*-
"""
将 BioDriver / hla_typing.json 中的等位基因串，转为 NetMHCIIpan 更常见的写法。

- 先查 data/hla_allele_map_netmhciipan.json 中 manual_overrides 手工表；
- 再按简规则补全 HLA- 前缀（如 DRB1*15:01 → HLA-DRB1*15:01）。

无生物学先验时只需：在 JSON 里加一对 from/to，或保持上游已带 HLA- 的写法即可。
"""
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional


def _default_json_path() -> str:
    return os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "data", "hla_allele_map_netmhciipan.json")
    )


def load_mapping(path: Optional[str] = None) -> Dict[str, Any]:
    """读入映射 JSON；path 缺省为仓库内 data/hla_allele_map_netmhciipan.json。"""
    p = path or _default_json_path()
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def _overrides_list(data: Dict[str, Any]) -> List[Dict[str, str]]:
    raw = data.get("manual_overrides")
    if not isinstance(raw, list):
        return []
    return [x for x in raw if isinstance(x, dict)]


def to_netmhciipan(allele: str, mapping_data: Optional[Dict[str, Any]] = None) -> str:
    """
    将一条 MHC-II 等位基因从「契约里可能看到的写法」转为 NetMHC 常用形式。

    - 非 II 类（如 HLA-A*…）不强行改写，原样返回，供调用方与 MHC-I 工具区分使用。
    - 已带 HLA- 的 II 类：仅统一基因名大小写，不动数字部分。
    """
    s = (allele or "").strip()
    if not s:
        return s

    data = mapping_data if mapping_data is not None else load_mapping()
    for item in _overrides_list(data):
        if item.get("from_biodriver") == s:
            return (item.get("to_netmhciipan") or s).strip()

    # HLA-DRB1*… 等：归一化基因名大写
    m_hla = re.match(
        r"^HLA-((?:DRB1|DQA1|DQB1|DPA1|DPB1))(\*.*)$",
        s,
        re.I,
    )
    if m_hla:
        g = m_hla.group(1).upper()
        return f"HLA-{g}{m_hla.group(2).strip()}"

    # 简写 DRB1*15:01
    m_short = re.match(
        r"^((?:DRB1|DQA1|DQB1|DPA1|DPB1))(\*.*)$",
        s,
        re.I,
    )
    if m_short:
        g = m_short.group(1).upper()
        return f"HLA-{g}{m_short.group(2).strip()}"

    return s


def batch_to_netmhciipan(alleles: List[str], mapping_data: Optional[Dict[str, Any]] = None) -> List[str]:
    data = mapping_data if mapping_data is not None else load_mapping()
    return [to_netmhciipan(a, data) for a in alleles]


def main():
    if len(sys.argv) < 2:
        print("用法: python hla_allele_to_netmhciipan.py <一条等位基因> [再多条…]")
        print("示例: python hla_allele_to_netmhciipan.py DRB1*15:01 HLA-DRB1*15:01")
        sys.exit(0)
    data = load_mapping()
    for a in sys.argv[1:]:
        print(f"{a} -> {to_netmhciipan(a, data)}")


if __name__ == "__main__":
    main()
