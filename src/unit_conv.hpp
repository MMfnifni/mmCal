#pragma once
#include "core.hpp"

namespace mm::cal {
 inline Value deg2rad_v(const Value &v, size_t pos) {
  constexpr double k = PI / 180.0;
  return mulReal(v, k, pos);
 }

 inline Value rad2deg_v(const Value &v, size_t pos) {
  constexpr double k = 180.0 / PI;
  return mulReal(v, k, pos);
 }

 inline Value deg2grad_v(const Value &v, size_t pos) { return mulReal(v, 400.0 / 360.0, pos); }

 inline Value grad2deg_v(const Value &v, size_t pos) { return mulReal(v, 360.0 / 400.0, pos); }

 inline Value rad2grad_v(const Value &v, size_t pos) {
  constexpr double k = 200.0 / PI;
  return mulReal(v, k, pos);
 }

 inline Value grad2rad_v(const Value &v, size_t pos) {
  constexpr double k = PI / 200.0;
  return mulReal(v, k, pos);
 }

 const std::unordered_map<std::string_view, std::unordered_map<std::string_view, double>> unitTable = {

     {               "length", {{"m", 1.0}, {"mm", 0.001}, {"cm", 0.01}, {"km", 1000.0}, {"inch", 0.0254}, {"ft", 0.3048}, {"mi", 1609.344}, {"micron", 0.000001}, {"\u00B5m", 0.000001}, {"microm", 0.000001}, {"nm", 1e-9}, {"pm", 1e-12}}},
     {                 "mass",                                 {{"kg", 1.0}, {"g", 0.001}, {"mg", 0.000001}, {"ton", 1000.0}, {"lb", 0.453592}, {"oz", 0.0283495}, {"gr", 0.0000647989}, {"\u00B5g", 1e-9}, {"mg", 1e-6}, {"tonne", 1000.0}}},
     {                 "time",                                                     {{"s", 1.0}, {"ms", 0.001}, {"\u00B5s", 1e-6}, {"ns", 1e-9}, {"min", 60.0}, {"hour", 3600.0}, {"day", 86400.0}, {"week", 604800.0}, {"year", 31536000.0}}},
     {          "temperature",                                                                                                                                                                          {{"K", 1.0}, {"C", 1.0}, {"F", 1.0}}},
     {     "electric_current",                                                                                                                                                  {{"A", 1.0}, {"mA", 0.001}, {"\u00B5A", 1e-6}, {"nA", 1e-9}}},
     {  "amount_of_substance",                                                                                                                                                          {{"mol", 1.0}, {"mmol", 0.001}, {"\u00B5mol", 1e-6}}},
     {   "luminous_intensity",                                                                                                                                                                                 {{"cd", 1.0}, {"mcd", 0.001}}},
     {                 "area",                                                                                                    {{"m2", 1.0}, {"cm2", 0.0001}, {"mm2", 1e-6}, {"km2", 1000000.0}, {"acre", 4046.86}, {"hectare", 10000.0}}},
     {               "volume",                                          {{"m3", 1.0}, {"cm3", 0.000001}, {"mm3", 1e-9}, {"L", 0.001}, {"mL", 0.000001}, {"gal", 0.00378541}, {"qt", 0.000946353}, {"pt", 0.000473176}, {"cup", 0.000236588}}},
     {             "velocity",                                                                                                                    {{"m/s", 1.0}, {"km/h", 0.277778}, {"mph", 0.44704}, {"ft/s", 0.3048}, {"knot", 0.514444}}},
     {         "acceleration",                                                                                                                                                                               {{"m/s2", 1.0}, {"g", 9.80665}}},
     {                "force",                                                                                                                            {{"N", 1.0}, {"kN", 1000.0}, {"dyn", 0.00001}, {"lbf", 4.44822}, {"kgf", 9.80665}}},
     {             "pressure",                                                              {{"Pa", 1.0}, {"kPa", 1000.0}, {"MPa", 1000000.0}, {"bar", 100000.0}, {"psi", 6894.76}, {"atm", 101325.0}, {"mmHg", 133.322}, {"torr", 133.322}}},
     {               "energy",                                                        {{"J", 1.0}, {"kJ", 1000.0}, {"MJ", 1000000.0}, {"cal", 4.184}, {"kcal", 4184.0}, {"Wh", 3600.0}, {"kWh", 3600000.0}, {"BTU", 1055.06}, {"erg", 1e-7}}},
     {                "power",                                                                                                                               {{"W", 1.0}, {"kW", 1000.0}, {"MW", 1000000.0}, {"hp", 745.7}, {"cv", 735.499}}},
     {      "electric_charge",                                                                                                                                                  {{"C", 1.0}, {"mC", 0.001}, {"\u00B5C", 1e-6}, {"nC", 1e-9}}},
     {   "electric_potential",                                                                                                                                                                   {{"V", 1.0}, {"kV", 1000.0}, {"mV", 0.001}}},
     { "electric_capacitance",                                                                                                                                                  {{"F", 1.0}, {"\u00B5F", 1e-6}, {"nF", 1e-9}, {"pF", 1e-12}}},
     {  "electric_resistance",                                                                                                                                                               {{"Ω", 1.0}, {"kΩ", 1000.0}, {"MΩ", 1000000.0}}},
     { "electric_conductance",                                                                                                                                                                {{"S", 1.0}, {"mS", 0.001}, {"\u00B5S", 1e-6}}},
     {        "magnetic_flux",                                                                                                                                                             {{"Wb", 1.0}, {"mWb", 0.001}, {"\u00B5Wb", 1e-6}}},
     {"magnetic_flux_density",                                                                                                                                                                {{"T", 1.0}, {"mT", 0.001}, {"\u00B5T", 1e-6}}},
     {            "frequency",                                                                                                                                     {{"Hz", 1.0}, {"kHz", 1000.0}, {"MHz", 1000000.0}, {"GHz", 1000000000.0}}},
     {        "heat_capacity",                                                                                                                                                                              {{"J/K", 1.0}, {"kJ/K", 1000.0}}},
     {        "specific_heat",                                                                                                                                                                    {{"J/(kg·K)", 1.0}, {"kJ/(kg·K)", 1000.0}}},
     { "thermal_conductivity",                                                                                                                                                                       {{"W/(m·K)", 1.0}, {"mW/(m·K)", 0.001}}},
     {            "viscosity",                                                                                                                                                                                {{"Pa·s", 1.0}, {"cP", 0.001}}},
     {        "ampere_second",                                                                                                                                                                                                  {{"C", 1.0}}},
     {  "molar_concentration",                                                                                                                                                                           {{"mol/L", 1.0}, {"mol/m3", 0.001}}}
 };
} // namespace mm::cal
