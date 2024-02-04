#pragma once
// теплофая эффективность
#include <vector>
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>

// поверхность
#define weight_plate = 0.2;
#define lenght_plate = 0.292;

//канал
#define high0 = 0.02;
#define lenght0 = 0.36;
#define weight0 = 0.15;

enum class InitialConditions : int {
	NUM_OF_EXPERIMENT = 0,
	T_IN_JET_EXIT,
	T_AIR,
	T_LOWER_THERMOCOUPLE,
	T_UPPER_TERMOCOUPLE,
	VELOCITY_IN_CHANNEL,
	JET_FLOW,
	FLOW_RATIO,
	X0,
	YO,
	X_LOWER_INDICATOR,
	Y_LOWER_INDICATOR,
	X_UPPER_INDICATOR,
	Y_UPPER_INDICATOR,
};
float LinearInterpolation(std::map<float, float> a, float temperature);