from collections import defaultdict
import re

# =========================================
# ここに元のマップを丸ごと貼る
# 例：
map_text = """
{ "pi", 3.141592653589793238462643383279502884},
{ "Pi", 3.141592653589793238462643383279502884},
{ "tau", 6.283185307179586476925286766559005768},
{ "Tau", 6.283185307179586476925286766559005768},
{ "c", 299792458.0},
{ "speed_of_light", 299792458.0},
{ "g0", 9.80665},
{ "g", 9.80665},
{ "kgf", 9.80665},
{ "stefan_boltzmann", 5.670374419e-8},
{ "sigma_sb", 5.670374419e-8},
{ "sigma", 5.670374419e-8},
{ "stefan_boltzmann_constant", 5.670374419e-8},
{ "wien_b", 2.897771955e-3},
{ "b_wien", 2.897771955e-3},
{ "wien", 2.897771955e-3},
{ "wien_displacement", 2.897771955e-3},
"""
# =========================================

# 1. マップから (key, value) ペアを抽出
pattern = re.compile(r'{\s*"([^"]+)"\s*,\s*([0-9.eE+\-]+)\s*}')
map_list = [(m.group(1), float(m.group(2))) for m in pattern.finditer(map_text)]

# 2. 数値ごとのキーリストを作る
val_to_keys = defaultdict(list)
for k, v in map_list:
    val_to_keys[v].append(k)

# 3. 重複している値だけ cnst_XXX を作る
cnst_lines = []
replacements = {}
used_names = set()

def make_safe_name(name):
    base = f"cnst_{name}"
    if base not in used_names:
        used_names.add(base)
        return base
    i = 2
    while f"{base}_{i}" in used_names:
        i += 1
    base = f"{base}_{i}"
    used_names.add(base)
    return base

for v, keys in val_to_keys.items():
    if len(keys) > 1:
        cnst_name = make_safe_name(keys[0])
        cnst_lines.append(f"const double {cnst_name} = {v};")
        for k in keys:
            replacements[k] = cnst_name

# 4. 元のマップを書き換える
new_map_lines = []
for k, v in map_list:
    if k in replacements:
        new_map_lines.append(f'{{ "{k}", {replacements[k]} }},')
    else:
        new_map_lines.append(f'{{ "{k}", {v} }},')  # 重複しないものはそのまま

# 5. 出力
print("// ==== 定数定義 ====")
print("\n".join(cnst_lines))
print("\n// ==== マップ ====")
print("\n".join(new_map_lines))