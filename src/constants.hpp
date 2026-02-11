#pragma once
// 定数辞書(MKS: m, kg, s)
//   - 値は基本的に SI(MKS) で統一(m, kg, s, K, A, mol, cd)
//   - unordered_map は array から構築する(初期化順の問題回避)
// ============================================================

#include <string_view>
#include <unordered_map>

#include "core.hpp"

namespace mm::cal {

 inline constexpr double cnst_pi = 3.141592653589793238462643383279502884;
 inline constexpr double cnst_e = 2.718281828459045235360287471352662498;
 inline constexpr double cnst_phi = 1.618033988749894848204586834365638118;

 inline constexpr double cnst_c = 299792458.0;
 inline constexpr double cnst_g0 = 9.80665;
 inline constexpr double cnst_G = 6.67430e-11;
 inline constexpr double cnst_h = 6.62607015e-34;
 inline constexpr double cnst_hbar = 1.054571817e-34;
 inline constexpr double cnst_qe = 1.602176634e-19;
 inline constexpr double cnst_NA = 6.02214076e23;
 inline constexpr double cnst_kB = 1.380649e-23;
 inline constexpr double cnst_R = 8.31446261815324;

 inline constexpr double cnst_mu0 = 1.25663706212e-6;
 inline constexpr double cnst_eps0 = 8.8541878128e-12;
 inline constexpr double cnst_Z0 = 376.730313668;

 inline constexpr double cnst_sigma_sb = 5.670374419e-8;
 inline constexpr double cnst_b_wien = 2.897771955e-3;
 inline constexpr double cnst_F = 96485.33212;
 inline constexpr double cnst_k_e = 8.9875517923e9;
 inline constexpr double cnst_au = 149597870700.0;

 // ============================================================
 // internal: raw entries
 // ============================================================

 struct ConstEntry {
   std::string_view key;
   double value;
 };

 // - 値は SI(MKS) の数値
 // - 角度は rad
 // - 温度は K
 // - 圧力は Pa
 // - エネルギーは J
 // - 電気は A, C, V, Ω, F, H
 //
 // ============================================================

 inline constexpr ConstEntry kConstTable[] = {
     // ============================================================
     // 基本数学 / Basic Math
     // ============================================================
     {                       "pi",                                cnst_pi},
     {                       "Pi",                                cnst_pi},
     {                      "tau",                            cnst_pi * 2}, // 2*pi
     {                      "Tau",                            cnst_pi * 2},
     {                        "e",                                 cnst_e},
     {                        "E",                                 cnst_e},
     {                      "phi",                               cnst_phi}, // golden ratio
     {                      "Phi",                               cnst_phi},
     {                    "sqrt2", 1.414213562373095048801688724209698079},
     {                    "sqrt3", 1.732050807568877293527446341505872367},
     {                      "ln2", 0.693147180559945309417232121458176568},
     {                     "ln10", 2.302585092994045684017991454684364208},
     {                    "log2e", 1.442695040888963407359924681001892137},
     {                   "log10e", 0.434294481903251827651128918916605082},
     {                   "inv_pi", 0.318309886183790671537767526745028724}, // 1/pi
     {                  "inv_tau", 0.159154943091895335768883763372514362}, // 1/tau
     {                "pi_over_2", 1.570796326794896619231321691639751442},
     {                "pi_over_3", 1.047197551196597746154214461093167628},
     {                "pi_over_4", 0.785398163397448309615660845819875721},
     {                "pi_over_6", 0.523598775598298873077107230546583814},
     {            "two_pi_over_3", 2.094395102393195492308428922186335257},
     {            "two_pi_over_5", 1.256637061435917295385057353311801154},

     // ============================================================
     // SI接頭語(無次元) / SI Prefix (dimensionless)
     // ============================================================
     {                     "kilo",                                    1e3},
     {                     "mega",                                    1e6},
     {                     "giga",                                    1e9},
     {                     "tera",                                   1e12},
     {                     "peta",                                   1e15},
     {                      "exa",                                   1e18},
     {                    "zetta",                                   1e21},
     {                    "yotta",                                   1e24},
     {                    "milli",                                   1e-3},
     {                    "micro",                                   1e-6},
     {                     "nano",                                   1e-9},
     {                     "pico",                                  1e-12},
     {                    "femto",                                  1e-15},
     {                     "atto",                                  1e-18},
     {                    "zepto",                                  1e-21},
     {                    "yocto",                                  1e-24},

     // ============================================================
     // 物理定数 / Physics Constants (CODATA style)
     // ============================================================
     {                        "c",                                 cnst_c}, // 光速, speed of light [m/s]
     {           "speed_of_light",                                 cnst_c},
     {                "celeritas",                                 cnst_c},
     {                       "g0",                                cnst_g0}, // 標準重力, standard gravity [m/s^2]
     {         "standard_gravity",                                cnst_g0},
     {                        "g",                                cnst_g0},
     {                        "G",                                 cnst_G}, // 万有引力定数, gravitational constant [m^3/(kg s^2)]
     {               "grav_const",                                 cnst_G},
     {   "gravitational_constant",                                 cnst_G},
     {                        "h",                                 cnst_h}, // プランク定数, Planck constant [J s]
     {                   "planck",                                 cnst_h},
     {          "planck_constant",                                 cnst_h},
     {                     "hbar",                              cnst_hbar}, // ディラック定数, reduced Planck constant [J s]
     {           "reduced_planck",                              cnst_hbar},
     {           "planck_reduced",                              cnst_hbar},
     {                 "e_charge",                                cnst_qe}, // 電気素量, elementary charge [C]
     {                       "qe",                                cnst_qe},
     {        "elementary_charge",                                cnst_qe},
     {                       "na",                                cnst_NA}, // アボガドロ定数, Avogadro constant [1/mol]
     {                 "avogadro",                                cnst_NA},
     {        "avogadro_constant",                                cnst_NA},
     {                       "kB",                                cnst_kB}, // ボルツマン定数, Boltzmann constant [J/K]
     {                "boltzmann",                                cnst_kB},
     {       "boltzmann_constant",                                cnst_kB},
     {                        "R",                                 cnst_R}, // 気体定数, gas constant [J/(mol K)]
     {             "gas_constant",                                 cnst_R},
     {       "molar_gas_constant",                                 cnst_R},
     {                      "mu0",                               cnst_mu0}, // 真空の透磁率, vacuum permeability [N/A^2]
     {                "vacuum_mu",                               cnst_mu0},
     {      "vacuum_permeability",                               cnst_mu0},
     {                     "eps0",                              cnst_eps0}, // 真空の誘電率, vacuum permittivity [F/m]
     {               "vacuum_eps",                              cnst_eps0},
     {      "vacuum_permittivity",                              cnst_eps0},
     {                       "Z0",                                cnst_Z0}, // 真空インピーダンス, impedance of free space [ohm]
     {         "vacuum_impedance",                                cnst_Z0},
     {     "free_space_impedance",                                cnst_Z0},
     {                 "sigma_sb",                          cnst_sigma_sb}, // ステファン=ボルツマン, Stefan-Boltzmann [W/(m^2 K^4)]
     {         "stefan_boltzmann",                          cnst_sigma_sb},
     {"stefan_boltzmann_constant",                          cnst_sigma_sb},
     {                   "b_wien",                            cnst_b_wien}, // ウィーン変位則, Wien displacement [m K]
     {                     "wien",                            cnst_b_wien},
     {        "wien_displacement",                            cnst_b_wien},

     // ============================================================
     // 熱・気体 / Thermodynamics & Gas
     // ============================================================
     {                      "atm",                               101325.0}, // 標準大気圧, standard atmosphere [Pa]
     {             "standard_atm",                               101325.0},
     {               "atmosphere",                               101325.0},
     {                      "bar",                               100000.0}, // bar [Pa]
     {                     "torr",                     133.32236842105263}, // torr [Pa]
     {                     "mmhg",                          133.322387415}, // mmHg [Pa] (approx)
     {                      "cal",                                  4.184}, // カロリー, calorie [J]
     {                     "kcal",                                 4184.0}, // kcal [J]
     {                      "btu",                          1055.05585262}, // BTU (IT) [J]
     {                   "wien_b",                            cnst_b_wien}, // ウィーン定数, Wien displacement [m K]
     {                "c_p_water",                                 4181.3}, // 水の比熱, specific heat of water [J/(kg K)]
     {                 "cp_water",                                 4181.3},
     {                  "c_p_air",                                 1005.0}, // 空気の比熱, specific heat of air [J/(kg K)]
     {                   "cp_air",                                 1005.0},
     {                "gamma_air",                                    1.4}, // 比熱比, heat capacity ratio of air
     {                   "pr_air",                                   0.71}, // プラントル数, Prandtl number (air)
     {                  "k_water",                                  0.598}, // 熱伝導率, thermal conductivity water [W/(m K)]
     {                    "k_air",                                 0.0257}, // thermal conductivity air [W/(m K)]
     {              "alpha_water",                                1.43e-7}, // 熱拡散率, thermal diffusivity water [m^2/s]
     {                "alpha_air",                                 2.2e-5}, // thermal diffusivity air [m^2/s]

     // ============================================================
     // 放射・光学 / Radiation & Optics
     // ============================================================
     {                "c1_planck",                        3.741771852e-16}, // プランク第1放射定数, 1st radiation constant [W m^2]
     {                "c2_planck",                         1.438776877e-2}, // プランク第2放射定数, 2nd radiation constant [m K]

     // ============================================================
     // 電気・磁気 / Electromagnetism
     // ============================================================
     {                  "faraday",                                 cnst_F}, // ファラデー定数, Faraday constant [C/mol]
     {                        "F",                                 cnst_F},
     {         "faraday_constant",                                 cnst_F},
     {                      "k_e",                               cnst_k_e}, // クーロン定数, Coulomb constant [N m^2 / C^2]
     {                "coulomb_k",                               cnst_k_e},
     {         "coulomb_constant",                               cnst_k_e},
     {                 "epsilon0",                              cnst_eps0}, // 真空誘電率, vacuum permittivity [F/m]
     {                     "mu_0",                               cnst_mu0}, // vacuum permeability [N/A^2]
     {                "epsilon_0",                              cnst_eps0},
     {                "mu0_exact",                               cnst_mu0}, // 近似(2019以降は厳密ではない)
     {                       "z0",                                cnst_Z0}, // 真空インピーダンス, impedance of free space [ohm]
     {                     "mu_b",                       9.2740100783e-24}, // ボーア磁子, Bohr magneton [J/T]
     {            "bohr_magneton",                       9.2740100783e-24},
     {                     "mu_n",                       5.0507837461e-27}, // 核磁子, nuclear magneton [J/T]
     {         "nuclear_magneton",                       5.0507837461e-27},

     // ============================================================
     // 天文 / Astronomy
     // ============================================================
     {                       "au",                                cnst_au}, // 天文単位, astronomical unit [m]
     {        "astronomical_unit",                                cnst_au},
     {                       "ua",                                cnst_au},
     {                       "ly",                     9.4607304725808e15}, // 光年, light-year [m]
     {               "light_year",                     9.4607304725808e15},
     {                "lightyear",                     9.4607304725808e15},

     {                       "pc",                  3.0856775814913673e16}, // パーセク, parsec [m]
     {                   "parsec",                  3.0856775814913673e16},
     {          "parallax_second",                  3.0856775814913673e16},

     {               "earth_mass",                              5.9722e24}, // 地球質量, Earth mass [kg]
     {                  "m_earth",                              5.9722e24},
     {                       "Me",                              5.9722e24},

     {             "earth_radius",                                6.371e6}, // 地球半径, Earth radius [m]
     {                  "r_earth",                                6.371e6},
     {                       "Re",                                6.371e6},

     {                 "sun_mass",                             1.98847e30}, // 太陽質量, solar mass [kg]
     {                    "m_sun",                             1.98847e30},
     {                       "Ms",                             1.98847e30},

     {               "sun_radius",                                6.957e8}, // 太陽半径, solar radius [m]
     {                    "r_sun",                                6.957e8},
     {                       "Rs",                                6.957e8},

     {                "moon_mass",                               7.342e22}, // 月質量, Moon mass [kg]
     {                   "m_moon",                               7.342e22},
     {                       "Mm",                               7.342e22},

     {              "moon_radius",                               1.7374e6}, // 月半径, Moon radius [m]
     {                   "r_moon",                               1.7374e6},
     {                       "Rm",                               1.7374e6},

     {           "solar_constant",                                 1361.0}, // 太陽定数, solar constant [W/m^2]
     {                 "earth_gm",                         3.986004418e14}, // 地球標準重力パラメータ, GM Earth [m^3/s^2]
     {                   "sun_gm",                       1.32712440018e20}, // GM Sun [m^3/s^2]
     {                  "moon_gm",                           4.9048695e12}, // GM Moon [m^3/s^2]

     // ============================================================
     // 材料・機械(代表値) / Materials & Mechanical (typical)
     // ============================================================
     {                "rho_water",                                 1000.0}, // 水密度, water density [kg/m^3]
     {            "water_density",                                 1000.0},
     {                  "rho_air",                                  1.225}, // 空気密度, air density [kg/m^3] at sea level
     {              "air_density",                                  1.225},

     {                 "nu_steel",                                   0.30}, // ポアソン比(鋼代表), Poisson ratio (steel typical)
     {            "poisson_steel",                                   0.30},
     {            "steel_poisson",                                   0.30},

     {                  "E_steel",                                2.05e11}, // ヤング率, Young's modulus steel [Pa]
     {              "young_steel",                                2.05e11},
     {              "steel_young",                                2.05e11},

     {                     "E_al",                                 6.9e10}, // ヤング率, Young's modulus aluminum [Pa]
     {                 "young_al",                                 6.9e10},
     {                 "al_young",                                 6.9e10},

     {                     "E_cu",                                1.10e11}, // ヤング率, Young's modulus copper [Pa]
     {                 "young_cu",                                1.10e11},
     {                 "cu_young",                                1.10e11},

     {                     "E_ti",                                1.05e11}, // ヤング率, Young's modulus titanium [Pa]
     {                 "young_ti",                                1.05e11},
     {                 "ti_young",                                1.05e11},

     {                "rho_steel",                                 7850.0}, // 密度, steel density [kg/m^3]
     {            "steel_density",                                 7850.0},
     {                   "rho_al",                                 2700.0}, // density aluminum
     {               "al_density",                                 2700.0},
     {                   "rho_cu",                                 8960.0}, // density copper
     {               "cu_density",                                 8960.0},
     {                   "rho_ti",                                 4500.0}, // density titanium
     {               "ti_density",                                 4500.0},

     {         "yield_steel_mild",                                 2.35e8}, // 降伏点, mild steel yield [Pa]
     {             "sigma_y_mild",                                 2.35e8},
     {         "mild_steel_yield",                                 2.35e8},

     {           "uts_steel_mild",                                  4.0e8}, // 引張強さ, mild steel UTS [Pa]
     {             "sigma_u_mild",                                  4.0e8},
     {           "mild_steel_uts",                                  4.0e8},

     // ============================================================
     // 流体 / Fluid
     // ============================================================
     {             "mu_water_20c",                               1.002e-3}, // 粘度, water dynamic viscosity [Pa s]
     {      "viscosity_water_20c",                               1.002e-3},
     {             "nu_water_20c",                               1.004e-6}, // 動粘度, kinematic viscosity [m^2/s]
     {      "kinematic_water_20c",                               1.004e-6},

     {               "mu_air_20c",                                1.81e-5}, // air dynamic viscosity [Pa s]
     {        "viscosity_air_20c",                                1.81e-5},
     {               "nu_air_20c",                                1.48e-5}, // air kinematic viscosity [m^2/s]
     {        "kinematic_air_20c",                                1.48e-5},

     {                    "nu_al",                                   0.33}, // ポアソン比, Poisson ratio aluminum
     {                    "nu_cu",                                   0.34}, // Poisson ratio copper
     {                    "nu_ti",                                   0.34}, // Poisson ratio titanium

     {                  "g_steel",                                 7.9e10}, // せん断弾性係数, shear modulus steel [Pa]
     {                     "g_al",                                 2.6e10}, // shear modulus aluminum [Pa]
     {                     "g_cu",                                 4.4e10}, // shear modulus copper [Pa]
     {                     "g_ti",                                 4.1e10}, // shear modulus titanium [Pa]

     {                  "k_steel",                                   50.0}, // 熱伝導率, thermal conductivity steel [W/(m K)]
     {                     "k_al",                                  205.0}, // aluminum
     {                     "k_cu",                                  385.0}, // copper
     {                     "k_ti",                                   21.9}, // titanium

     {              "alpha_steel",                                  12e-6}, // 線膨張係数, thermal expansion steel [1/K]
     {                 "alpha_al",                                  23e-6}, // aluminum
     {                 "alpha_cu",                                  17e-6}, // copper
     {                 "alpha_ti",                                 8.6e-6}, // titanium

     // ============================================================
     // 化学 / Chemistry
     // ============================================================
     {                        "f",                                 cnst_F}, // ファラデー定数, Faraday constant [C/mol]
     {            "faraday_const",                                 cnst_F},

     {                      "amu",                      1.66053906660e-27}, // 原子質量単位, atomic mass unit [kg]
     {                        "u",                      1.66053906660e-27},
     {         "atomic_mass_unit",                      1.66053906660e-27},

     // ============================================================
     // 換算(長さ) / Length conversion (to meters)
     // ============================================================
     {                     "inch",                                 0.0254}, // inch -> m
     {                       "in",                                 0.0254},
     {                     "foot",                                 0.3048}, // ft -> m
     {                       "ft",                                 0.3048},
     {                     "yard",                                 0.9144}, // yd -> m
     {                       "yd",                                 0.9144},
     {                     "mile",                               1609.344}, // mi -> m
     {                       "mi",                               1609.344},
     {            "nautical_mile",                                 1852.0}, // nmi -> m
     {                      "nmi",                                 1852.0},
     {                 "angstrom",                                  1e-10}, // Å -> m
     {                       "aa",                                  1e-10},
     {                   "micron",                                   1e-6}, // μm -> m
     {                       "um",                                   1e-6},

     // ============================================================
     // 換算(質量) / Mass conversion (to kg)
     // ============================================================
     {                       "lb",                             0.45359237}, // pound -> kg
     {                    "pound",                             0.45359237},
     {                       "oz",                         0.028349523125}, // ounce -> kg
     {                    "ounce",                         0.028349523125},
     {                   "ton_us",                              907.18474}, // short ton -> kg
     {                "short_ton",                              907.18474},
     {                   "ton_uk",                           1016.0469088}, // long ton -> kg
     {                 "long_ton",                           1016.0469088},

     // ============================================================
     // 換算(力) / Force conversion (to N)
     // ============================================================
     {                      "lbf",                        4.4482216152605}, // pound-force -> N
     {              "pound_force",                        4.4482216152605},
     {                      "kgf",                                cnst_g0}, // kilogram-force -> N
     {           "kilogram_force",                                cnst_g0},
     {                      "dyn",                                   1e-5}, // dyne -> N
     {                     "dyne",                                   1e-5},

     // ============================================================
     // 換算(圧力) / Pressure conversion (to Pa)
     // ============================================================
     {                      "psi",                         6894.757293168}, // psi -> Pa
     {                      "ksi",                       6.894757293168e6}, // ksi -> Pa
     {                      "mpa",                                    1e6},
     {                      "gpa",                                    1e9},
     {                      "kpa",                                    1e3},
     {                    "mmh2o",                                cnst_g0}, // mmH2O -> Pa (approx)
     {                    "inh2o",                              249.08891}, // inH2O -> Pa (approx)

     // ============================================================
     // 換算(時間) / Pressure conversion (to time)
     // ============================================================
     {                   "minute",                                   60.0}, // 分, minute [s]
     {                     "hour",                                 3600.0}, // 時, hour [s]
     {                      "day",                                86400.0}, // 日, day [s]
     {                     "week",                               604800.0}, // 週, week [s]
     {                     "year",                             31557600.0}, // 年, Julian year [s]

     // ============================================================
     // 換算(角度) / Pressure conversion (to angle)
     // ============================================================
     {                   "arcmin",                   2.908882086657216e-4}, // 角分, arcminute [rad]
     {                   "arcsec",                    4.84813681109536e-6}, // 角秒, arcsecond [rad]
     {               "deg_to_rad", 0.017453292519943295769236907684886127}, // pi/180
     {               "rad_to_deg", 57.29577951308232087679815481410517033}, // 180/pi
 };

 // ============================================================
 // Build unordered_map from array
 // ============================================================

 inline constexpr size_t kConstTableSize = sizeof(kConstTable) / sizeof(kConstTable[0]);

 inline std::unordered_map<std::string, Value> buildConstantsDic() {
  std::unordered_map<std::string, Value> m;
  m.reserve(kConstTableSize * 2);

  for (auto &e : kConstTable)
   m.emplace(e.key, Value(e.value));

  return m;
 }

 // ============================================================
 // exported dictionary
 // ============================================================

 std::unordered_map<std::string, Value> constants_dic = buildConstantsDic();

} // namespace mm::cal
