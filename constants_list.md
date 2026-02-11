# mmCalculator – Mechanical & Manufacturing Calculator

© 2021–2026 mmKreuzef (aka Daiki.NIIMI) Project mme
Licensed under the BSD 3-Clause License

---

戻る / Back <[English](README.md) | [日本語](README.ja.md)>

---

## cnst関数の定数一覧 / List of Constants for the cnst Function

- 値は基本的にSI(MKS)で統一(m`, `kg`, `s`, `K`, `A`, `mol`, `cd)，角度についてradで表記しているが，本電卓はdeg既定に注意
- Values are fundamentally unified in SI (MKS) units (m`, `kg`, `s`, `K`, `A`, `mol`, `cd)`, `Angles are displayed in radians`, `but note that this calculator defaults to degrees.

---

### 基本数学 / Basic Math

| 定数名 / Name   | 値 / Value                             | 説明 / Explanation                |
| --------------- | -------------------------------------- | --------------------------------- |
| `pi`, `Pi`      | 3.141592653589793238462643383279502884 | 円周率 / π                        |
| `tau`, `Tau`    | 6.28318530717958647692528676655900577  | 2π                                |
| `e`, `E`        | 2.718281828459045235360287471352662498 | ネイピア数 /Napiere               |
| `phi`, `Phi`    | 1.618033988749894848204586834365638118 | 黄金比 / Golden Ratio             |
| `sqrt2`         | 1.414213562373095048801688724209698079 | √2                                |
| `sqrt3`         | 1.732050807568877293527446341505872367 | √3                                |
| `ln2`           | 0.693147180559945309417232121458176568 | 自然対数2 / Natural logarithm 2   |
| `ln10`          | 2.302585092994045684017991454684364208 | 自然対数10 / Natural logarithm 10 |
| `log2e`         | 1.442695040888963407359924681001892137 | log2(e)                           |
| `log10e`        | 0.434294481903251827651128918916605082 | log10(e)                          |
| `inv_pi`        | 0.318309886183790671537767526745028724 | 1/π                               |
| `inv_tau`       | 0.159154943091895335768883763372514362 | 1/τ                               |
| `pi_over_2`     | 1.570796326794896619231321691639751442 | π/2                               |
| `pi_over_3`     | 1.047197551196597746154214461093167628 | π/3                               |
| `pi_over_4`     | 0.785398163397448309615660845819875721 | π/4                               |
| `pi_over_6`     | 0.523598775598298873077107230546583814 | π/6                               |
| `two_pi_over_3` | 2.094395102393195492308428922186335257 | 2π/3                              |
| `two_pi_over_5` | 1.256637061435917295385057353311801154 | 2π/5                              |

---

### SI接頭語(無次元) / SI Prefix (dimensionless)

| 定数名 / Name | 値 / Value | 説明 / Explanation |
| ------------- | ---------- | ------------------ |
| `kilo`        | 1e3        | キロ               |
| `mega`        | 1e6        | メガ               |
| `giga`        | 1e9        | ギガ               |
| `tera`        | 1e12       | テラ               |
| `peta`        | 1e15       | ペタ               |
| `exa`         | 1e18       | エクサ             |
| `zetta`       | 1e21       | ゼタ               |
| `yotta`       | 1e24       | ヨタ               |
| `milli`       | 1e-3       | ミリ               |
| `micro`       | 1e-6       | マイクロ           |
| `nano`        | 1e-9       | ナノ               |
| `pico`        | 1e-12      | ピコ               |
| `femto`       | 1e-15      | フェムト           |
| `atto`        | 1e-18      | アット             |
| `zepto`       | 1e-21      | ゼプト             |
| `yocto`       | 1e-24      | ヨクト             |

---

### 物理定数 / Physics Constants

| 定数名 / Name                                               | 値 / Value       | 説明 / Explanation                    |
| ----------------------------------------------------------- | ---------------- | ------------------------------------- |
| `c`, `speed_of_light`, `celeritas`                          | 299792458.0      | 光速 [m/s]                            |
| `g0`, `standard_gravity`, `g`                               | 9.80665          | 標準重力 [m/s²]                       |
| `G`, `grav_const`, `gravitational_constant`                 | 6.67430e-11      | 万有引力定数 [m³/(kg·s²)]             |
| `h`, `planck`, `planck_constant`                            | 6.62607015e-34   | プランク定数 [J·s]                    |
| `hbar`, `reduced_planck`, `planck_reduced`                  | 1.054571817e-34  | ディラック定数 [J·s]                  |
| `e_charge`, `qe`, `elementary_charge`                       | 1.602176634e-19  | 電気素量 [C]                          |
| `na`, `avogadro`, `avogadro_constant`                       | 6.02214076e23    | アボガドロ定数 [1/mol]                |
| `kB`, `boltzmann`, `boltzmann_constant`                     | 1.380649e-23     | ボルツマン定数 [J/K]                  |
| `R`, `gas_constant`, `molar_gas_constant`                   | 8.31446261815324 | 気体定数 [J/(mol·K)]                  |
| `mu0`, `vacuum_mu`, `vacuum_permeability`                   | 1.25663706212e-6 | 真空透磁率 [N/A²]                     |
| `eps0`, `vacuum_eps`, `vacuum_permittivity`                 | 8.8541878128e-12 | 真空誘電率 [F/m]                      |
| `Z0`, `vacuum_impedance`, `free_space_impedance`            | 376.730313668    | 真空インピーダンス [Ω]                |
| `sigma_sb`, `stefan_boltzmann`, `stefan_boltzmann_constant` | 5.670374419e-8   | ステファン=ボルツマン定数 [W/(m²·K⁴)] |
| `b_wien`, `wien`, `wien_displacement`                       | 2.897771955e-3   | ウィーン変位定数 [m·K]                |

---

### 熱・気体 / Thermodynamics & Gas

| 定数名 / Name                       | 値 / Value         | 説明 / Explanation     |
| ----------------------------------- | ------------------ | ---------------------- |
| `atm`, `standard_atm`, `atmosphere` | 101325.0           | 標準大気圧 [Pa]        |
| `bar`                               | 100000.0           | bar [Pa]               |
| `torr`                              | 133.32236842105263 | torr [Pa]              |
| `mmhg`                              | 133.322387415      | mmHg [Pa]              |
| `cal`                               | 4.184              | カロリー [J]           |
| `kcal`                              | 4184.0             | kcal [J]               |
| `btu`                               | 1055.05585262      | BTU [J]                |
| `wien_b`                            | 2.897771955e-3     | ウィーン定数 [m·K]     |
| `c_p_water`, `cp_water`             | 4181.3             | 水の比熱 [J/(kg·K)]    |
| `c_p_air`, `cp_air`                 | 1005.0             | 空気の比熱 [J/(kg·K)]  |
| `gamma_air`                         | 1.4                | 比熱比 (air)           |
| `pr_air`                            | 0.71               | プラントル数 (air)     |
| `k_water`                           | 0.598              | 熱伝導率水 [W/(m·K)]   |
| `k_air`                             | 0.0257             | 熱伝導率空気 [W/(m·K)] |
| `alpha_water`                       | 1.43e-7            | 熱拡散率水 [m²/s]      |
| `alpha_air`                         | 2.2e-5             | 熱拡散率空気 [m²/s]    |

---

### 放射・光学 / Radiation & Optics

| 定数名 / Name | 値 / Value      | 説明 / Explanation         |
| ------------- | --------------- | -------------------------- |
| `c1_planck`   | 3.741771852e-16 | プランク第1放射定数 [W·m²] |
| `c2_planck`   | 1.438776877e-2  | プランク第2放射定数 [m·K]  |

---

### 電気・磁気 / Electromagnetism

| 定数名 / Name                          | 値 / Value       | 説明 / Explanation     |
| -------------------------------------- | ---------------- | ---------------------- |
| `faraday`, `F`, `faraday_constant`     | 96485.33212      | ファラデー定数 [C/mol] |
| `k_e`, `coulomb_k`, `coulomb_constant` | 8.9875517923e9   | クーロン定数 [N·m²/C²] |
| `epsilon0`, `epsilon_0`                | 8.8541878128e-12 | 真空誘電率 [F/m]       |
| `mu_0`, `mu0_exact`                    | 1.25663706212e-6 | 真空透磁率 [N/A²]      |
| `z0`                                   | 376.730313668    | 真空インピーダンス [Ω] |
| `mu_b`, `bohr_magneton`                | 9.2740100783e-24 | ボーア磁子 [J/T]       |
| `mu_n`, `nuclear_magneton`             | 5.0507837461e-27 | 核磁子 [J/T]           |

---

### 天文 / Astronomy

| 定数名 / Name                     | 値 / Value            | 説明 / Explanation |
| --------------------------------- | --------------------- | ------------------ |
| `au`, `astronomical_unit`, `ua`   | 149597870700.0        | 天文単位 [m]       |
| `ly`, `light_year`, `lightyear`   | 9.4607304725808e15    | 光年 [m]           |
| `pc`, `parsec`, `parallax_second` | 3.0856775814913673e16 | パーセク [m]       |
| `earth_mass`, `m_earth`, `Me`     | 5.9722e24             | 地球質量 [kg]      |
| `earth_radius`, `r_earth`, `Re`   | 6.371e6               | 地球半径 [m]       |
| `sun_mass`, `m_sun`, `Ms`         | 1.98847e30            | 太陽質量 [kg]      |
| `sun_radius`, `r_sun`, `Rs`       | 6.957e8               | 太陽半径 [m]       |
| `moon_mass`, `m_moon`, `Mm`       | 7.342e22              | 月質量 [kg]        |
| `moon_radius`, `r_moon`, `Rm`     | 1.7374e6              | 月半径 [m]         |
| `solar_constant`                  | 1361.0                | 太陽定数 [W/m²]    |
| `earth_gm`                        | 3.986004418e14        | 地球GM [m³/s²]     |
| `sun_gm`                          | 1.32712440018e20      | 太陽GM [m³/s²]     |
| `moon_gm`                         | 4.9048695e12          | 月GM [m³/s²]       |

---

### 材料・機械 / Materials & Mechanical

| 定数名 / Name                                          | 値 / Value | 説明 / Explanation        |
| ------------------------------------------------------ | ---------- | ------------------------- |
| `rho_water`, `water_density`                           | 1000.0     | 水密度 [kg/m³]            |
| `rho_air`, `air_density`                               | 1.225      | 空気密度 [kg/m³]          |
| `nu_steel`, `poisson_steel`, `steel_poisson`           | 0.30       | ポアソン比(鋼)            |
| `E_steel`, `young_steel`, `steel_young`                | 2.05e11    | ヤング率(鋼) [Pa]         |
| `E_al`, `young_al`, `al_young`                         | 6.9e10     | ヤング率(アルミ) [Pa]     |
| `E_cu`, `young_cu`, `cu_young`                         | 1.10e11    | ヤング率(銅) [Pa]         |
| `E_ti`, `young_ti`, `ti_young`                         | 1.05e11    | ヤング率(チタン) [Pa]     |
| `rho_steel`, `steel_density`                           | 7850.0     | 鋼密度 [kg/m³]            |
| `rho_al`, `al_density`                                 | 2700.0     | アルミ密度 [kg/m³]        |
| `rho_cu`, `cu_density`                                 | 8960.0     | 銅密度 [kg/m³]            |
| `rho_ti`, `ti_density`                                 | 4500.0     | チタン密度 [kg/m³]        |
| `yield_steel_mild`, `sigma_y_mild`, `mild_steel_yield` | 2.35e8     | 降伏点 [Pa]               |
| `uts_steel_mild`, `sigma_u_mild`, `mild_steel_uts`     | 4.0e8      | 引張強さ [Pa]             |
| `g_steel`                                              | 7.9e10     | 剪断弾性係数(鉄) [Pa]     |
| `g_al`                                                 | 2.6e10     | 剪断弾性係数(アルミ) [Pa] |
| `g_cu`                                                 | 4.4e10     | 剪断弾性係数(銅) [Pa]     |
| `g_ti`                                                 | 4.1e10     | 剪断弾性係数(チタン) [Pa] |
| `k_steel`                                              | 50.0       | 熱伝導率 [W/(m·K)]        |
| `k_al`                                                 | 205.0      | 熱伝導率 [W/(m·K)]        |
| `k_cu`                                                 | 385.0      | 熱伝導率 [W/(m·K)]        |
| `k_ti`                                                 | 21.9       | 熱伝導率 [W/(m·K)]        |
| `alpha_steel`                                          | 12e-6      | 線膨張係数 [1/K]          |
| `alpha_al`                                             | 23e-6      | 線膨張係数 [1/K]          |
| `alpha_cu`                                             | 17e-6      | 線膨張係数 [1/K]          |
| `alpha_ti`                                             | 8.6e-6     | 線膨張係数 [1/K]          |
| `nu_al`                                                | 0.33       | ポアソン比(アルミ)        |
| `nu_cu`                                                | 0.34       | ポアソン比(銅)            |
| `nu_ti`                                                | 0.34       | ポアソン比(チタン)        |

---

### 流体 / Fluid

| 定数名 / Name                         | 値 / Value | 説明 / Explanation |
| ------------------------------------- | ---------- | ------------------ |
| `mu_water_20c`, `viscosity_water_20c` | 1.002e-3   | 水動粘度 [Pa·s]    |
| `nu_water_20c`, `kinematic_water_20c` | 1.004e-6   | 水動粘度 [m²/s]    |
| `mu_air_20c`, `viscosity_air_20c`     | 1.81e-5    | 空気動粘度 [Pa·s]  |
| `nu_air_20c`, `kinematic_air_20c`     | 1.48e-5    | 空気動粘度 [m²/s]  |

---

### 化学 / Chemistry

| 定数名 / Name                  | 値 / Value        | 説明 / Explanation     |
| ------------------------------ | ----------------- | ---------------------- |
| `f`, `faraday_const`           | 96485.33212       | ファラデー定数 [C/mol] |
| `amu`, `u`, `atomic_mass_unit` | 1.66053906660e-27 | 原子質量単位 [kg]      |

---

### 換算 / conversion

| 定数名 / Name           | 値 / Value                             | 説明 / Explanation     |
| ----------------------- | -------------------------------------- | ---------------------- |
| `inch`, `in`            | 0.0254                                 | インチ -> m            |
| `foot`, `ft`            | 0.3048                                 | フィート -> m          |
| `yard`, `yd`            | 0.9144                                 | ヤード -> m            |
| `mile`, `mi`            | 1609.344                               | マイル -> m            |
| `nautical_mile`, `nmi`  | 1852.0                                 | 海里 -> m              |
| `angstrom`, `aa`        | 1e-10                                  | Å -> m                 |
| `micron`, `um`          | 1e-6                                   | μm -> m                |
| `lb`, `pound`           | 0.45359237                             | ポンド -> kg           |
| `oz`, `ounce`           | 0.028349523125                         | オンス -> kg           |
| `ton_us`, `short_ton`   | 907.18474                              | ショートトン -> kg     |
| `ton_uk`, `long_ton`    | 1016.0469088                           | ロングトン -> kg       |
| `lb`, `pound`           | 0.45359237                             | ポンド -> kg           |
| `oz`, `ounce`           | 0.028349523125                         | オンス -> kg           |
| `ton_us`, `short_ton`   | 907.18474                              | ショートトン -> kg     |
| `ton_uk`, `long_ton`    | 1016.0469088                           | ロングトン -> kg       |
| `lbf`, `pound_force`    | 4.4482216152605                        | ポンド力 -> N          |
| `kgf`, `kilogram_force` | 9.80665                                | キログラム力 -> N      |
| `dyn`, `dyne`           | 1e-5                                   | ダイン -> N            |
| `psi`                   | 6894.757293168                         | psi -> Pa              |
| `ksi`                   | 6.894757293168e6                       | ksi -> Pa              |
| `mpa`                   | 1e6                                    | MPa -> Pa              |
| `gpa`                   | 1e9                                    | GPa -> Pa              |
| `kpa`                   | 1e3                                    | kPa -> Pa              |
| `mmh2o`                 | 9.80665                                | mmH2O -> Pa            |
| `inh2o`                 | 249.08891                              | inH2O -> Pa            |
| `minute`                | 60.0                                   | 分 -> s                |
| `hour`                  | 3600.0                                 | 時間 -> s              |
| `day`                   | 86400.0                                | 日 -> s                |
| `week`                  | 604800.0                               | 週 -> s                |
| `year`                  | 31557600.0                             | 年 (Julian) -> s       |
| `minute`                | 60.0                                   | 分 -> s                |
| `hour`                  | 3600.0                                 | 時間 -> s              |
| `day`                   | 86400.0                                | 日 -> s                |
| `week`                  | 604800.0                               | 週 -> s                |
| `year`                  | 31557600.0                             | 年 (Julian) -> s       |
| `arcmin`                | 2.908882086657216e-4                   | 角分 -> rad            |
| `arcsec`                | 4.84813681109536e-6                    | 角秒 -> rad            |
| `deg_to_rad`            | 0.017453292519943295769236907684886127 | 度 -> ラジアン (π/180) |
| `rad_to_deg`            | 57.29577951308232087679815481410517033 | ラジアン -> 度 (180/π) |
