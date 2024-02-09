#pragma once
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "../data/data.h"
#include "../interpolate/interpolate.h"
#include "../read_begin_conditional.h"

class Nu { 
public:
	//Ra ~ Nu^4
	Nu(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionalsNu4>& bC, const std::string& x_name, const std::string& y_name, int i, float lenght_plate, float weight_plate, float epsilon_kraska, float epsilon_material, float diameter_jet);
	//Ra ~ Nu^3
	//Nu(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionalsNu4>& bC, const std::string& x_name, const std::string& y_name, int i, float lenght_plate, float weight_plate, float epsilon_kraska, float epsilon_material, float diameter_jet);
	void AreaNu();
	void AreaIntegralNu();
	void AlongXbNu();
	void AlongXNu(float x_0);
	void AlongYNu(float y_0);
	void AlongXbIntegralNu(float dx);
	void AlongXIntegralNu(float dx, float x_0);
	void AlongYIntegralNu(float dx, float y_0);
	void RoundIntegralNu(float x_0, float y_0, float dx);
	~Nu() = default;
private:
	int width_pix, height_pix, margin_from_edges_pix; 
	float linear_ax, linear_b;
	float x_begin_coord, y_begin_coord, non_dimensionalize_coeff;
	std::vector<beginConditionalsNu4> bC; 
	std::string x_name, y_name;
	int i;
	std::vector<float> X;
	std::vector<float> Y;
	std::vector<float> T_;
	std::vector<float> Tin_;
	std::vector<float> Nu_;
	std::vector<float> Qizl_;
	std::vector<float> Qel_;
	std::vector<float> Qck1_;
	std::vector<float> Qck2_;
	std::vector<float> Qw_;
	std::vector<float> Ra_;
	std::vector<float> Ra2_;
	std::vector<float> Alpha_;
	
	std::map<float, float> AirKinematicViscosity;
	std::map<float, float> AirDensity;
	std::map<float, float> AirSpecificHeat;
	std::map<float, float> AirThermalDiffusivity;
	std::map<float, float> AirPrandtl;
	std::map<float, float> AirLambdaPatametr;
	std::map<float, float> AirVolumetric—oefficientofThermalExpansion;
};

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
