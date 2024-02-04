#pragma once
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "../data/data.h"
#include "../interpolate/interpolate.h"
/*
class Nu_1Builder {
public:
	//operator Nu_1() const;

	Nu_1Builder& SetTair(float Tair);
	Nu_1Builder& SetEpsilonIn(float epsilon_in);
	Nu_1Builder& SetEpsilonOut(float epsilon_out);
	Nu_1Builder& SetDiameter(float diameter);
	Nu_1Builder& SetTbegin(float Tbegin);
	Nu_1Builder& SetN(float N);
	Nu_1Builder& SetCharacteristicSize(float characteristic_size);
	Nu_1Builder& SetXName(float X_Name);
	Nu_1Builder& SetYName(float Y_Name);
	Nu_1Builder& SetXbegin(float X_begin);
	Nu_1Builder& SetHigh(float high);
	Nu_1Builder& SetWeight(float weight);
	Nu_1Builder& SetLenght(float lenght);
	Nu_1Builder& SetVelocity(float velocity);
	Nu_1Builder& SetFilename(float filename);
	Nu_1Builder& SetNumberOFWall(float NumberofWall);
	Nu_1Builder& SetT(std::vector<float>& T);
	Nu_1Builder& SetX(std::vector<float>& X);
	Nu_1Builder& SetY(std::vector<float>& Y);
private:
	float _Tair;
	float _epsilon_in;
	float _epsilon_out;
	float _diameter;
	float _Tbegin;
	float _N;
	float _characteristic_size;
	std::string _X_Name;
	std::string _Y_Name;
	float _X_begin;
	float _high;
	float _weight;
	float _lenght; 
	float _velocity;
	std::string _filename;
	int _NumberofWall;

	std::vector<float> _T;
	std::vector<float> _X; // 
	std::vector<float> _Y; //
	
};

class Nu_1 { //Ra ~ T3
public:
	friend class Nu_1Builder;

	
	void PrintNu_1(const std::string& delimetr);
protected:
	Nu_1(float Tair, float epsilon_in, float epsilon_out, float diameter, float Tbegin, float N, 
		const std::vector<float>& T, const std::vector<float>& X, const std::vector<float>& Y, float characteristic_size,
		std::string X_Name, std::string Y_Name, float X_begin, float high, float weight, float lenght, 
		float velocity, std::string filename, int NumberofWall);
	std::vector<float> _T;
	std::vector<float> _Tin;
	std::vector<float> _Nu;
	std::vector<float> _Alpha;
	std::vector<float> _Qw; // общий тепловой поток
	std::vector<float> _Qel; // электрический тепловай поток 
	std::vector<float> _Qizl; // тепловой поток излучение
	std::vector<float> _Qck_in; // конвективный тепловой поток внутри объекта
	std::vector<float> _Qck_out; // конвективный тепловой поток снаружи объекта
	std::vector<float> _Ra_in; // 
	std::vector<float> _Ra_out; // 
	std::vector<float> _X; // 
	std::vector<float> _Y; //
	float _characteristic_size;
	std::string _X_Name;
	std::string _Y_Name;
	std::map<float, float> AirLambdaPatametr;
	std::map<float, float> AirKinematicViscosity;
	std::map<float, float> AirDensity;
	std::map<float, float> AirSpecificHeat;
	std::map<float, float> AirThermalDiffusivity;
	std::map<float, float> AirPrandtl;
	std::map<float, float> AirLambdaPatametr;
	std::string _filename;
	int _NumberofWall;
};

*/
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
