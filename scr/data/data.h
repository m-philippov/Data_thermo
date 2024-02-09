#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cmath>

#define PROPERTIES_OF_AIR(func) OpenPropertiesOFAir(func, #func)
//начальные параметры
const float g_free = 9.8;
const float T0 = 273.15;
const float sigma = 5.67 * static_cast<float>(pow(10.0, -8));


//считывает файл и отдает данные в вектор пар
std::map<float, float> AirCriterium(std::ifstream&);

//Ћинеализаци€
std::map<float, float> Linealize(const std::map<float, float>&);


void OpenPropertiesOFAir(std::map<float, float>&, const std::string&);

//поставл€ет данные "“еплопроводность воздуха в зависимости от температуры Ч таблица"
//AirLambdaPatametr


//поставл€ет данные "ѕлотность воздуха в зависимости от температуры Ч таблица"
// AirDensity

//поставл€ет данные "ƒинамическа€  в€зкость воздуха в зависимости от температуры Ч таблица"
// AirDynamicViscosity

//поставл€ет данные " инематическа€  в€зкость воздуха в зависимости от температуры Ч таблица"
//AirKinematicViscosity

//поставл€ет данные "”дельна€ теплоемкость воздуха при различных температурах Ч таблица"
//AirSpecificHeat

//поставл€ет данные "ѕлотность воздуха в зависимости от температуры Ч таблица"
//AirThermalConductivity
// 
//поставл€ет данные "ƒинамическа€  в€зкость воздуха в зависимости от температуры Ч таблица"
//AirThermalDiffusivity

//поставл€ет данные " инематическа€  в€зкость воздуха в зависимости от температуры Ч таблица"
//AirVolumetric—oefficientofThermalExpansion

//поставл€ет данные "”дельна€ теплоемкость воздуха при различных температурах Ч таблица"
//AirPrandtl

/**/
